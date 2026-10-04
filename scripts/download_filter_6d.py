"""Download every GaiaSource_*.csv.gz bulk chunk IN PARALLEL, stream-filter
each to only the rows with a measured radial_velocity (the full 6D
phase-space subset this project needs), and append those rows to a single
growing output CSV. Resumable: tracks completed chunks in progress.txt.

Gaia's bulk CSVs are actually ECSV (astropy): ~1000 lines of "#"-prefixed
YAML column metadata precede the real CSV header line, and null values are
the literal string "null", not an empty field -- both handled below.
"""
import csv
import gzip
import io
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

BASE = "https://gaia.eu-1.cdn77-storage.com/"
LIST_FILE = "gaia_file_list.txt"
OUT_FILE = "gaia_dr3_6d.csv"
PROGRESS_FILE = "progress.txt"
RV_COL = "radial_velocity"
N_WORKERS = 16

write_lock = threading.Lock()
progress_lock = threading.Lock()
header_written = threading.Event()


def load_file_list():
    files = []
    with open(LIST_FILE) as f:
        for line in f:
            key, size = line.strip().split("\t")
            if key.endswith(".csv.gz"):
                files.append((key, int(size)))
    return files


def load_progress():
    if not os.path.exists(PROGRESS_FILE):
        return set()
    with open(PROGRESS_FILE) as f:
        return set(line.strip() for line in f if line.strip())


def mark_done(key):
    with progress_lock:
        with open(PROGRESS_FILE, "a") as f:
            f.write(key + "\n")


def download_bytes(key):
    url = BASE + key
    r = requests.get(url, timeout=180)
    r.raise_for_status()
    return r.content


def filter_chunk(raw_gz_bytes):
    """Returns (header_row_or_None, [filtered_rows], total_rows)."""
    buf = io.BytesIO(raw_gz_bytes)
    with gzip.GzipFile(fileobj=buf) as gz:
        text_stream = io.TextIOWrapper(gz, encoding="utf-8", newline="")
        line = text_stream.readline()
        while line.startswith("#"):
            line = text_stream.readline()
        reader = csv.reader([line] + list(text_stream))
        header = next(reader)
        rv_idx = header.index(RV_COL)
        kept_rows = []
        total = 0
        for row in reader:
            total += 1
            if row[rv_idx] != "null" and row[rv_idx] != "":
                kept_rows.append(row)
        return header, kept_rows, total


def process_one(key):
    attempt = 0
    while True:
        try:
            raw = download_bytes(key)
            header, rows, total = filter_chunk(raw)
            return key, header, rows, total
        except Exception as e:
            attempt += 1
            if attempt >= 5:
                print(f"  GIVING UP on {key} after 5 attempts: {e}", flush=True)
                return key, None, [], 0
            time.sleep(5)


def main():
    files = load_file_list()
    done = load_progress()
    remaining = [(k, s) for k, s in files if k not in done]
    print(f"{len(files)} total chunks, {len(done)} already done, {len(remaining)} remaining, "
          f"{N_WORKERS} parallel workers", flush=True)

    if os.path.exists(OUT_FILE) and os.path.getsize(OUT_FILE) > 0:
        header_written.set()
    out_mode = "a" if header_written.is_set() else "w"
    outf = open(OUT_FILE, out_mode, newline="")
    writer = csv.writer(outf)

    t0 = time.time()
    total_kept = 0
    completed = 0

    with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(process_one, key): key for key, size in remaining}
        for future in as_completed(futures):
            key, header, rows, total = future.result()
            with write_lock:
                if header is not None and not header_written.is_set():
                    writer.writerow(header)
                    header_written.set()
                for row in rows:
                    writer.writerow(row)
                outf.flush()
            mark_done(key)
            completed += 1
            total_kept += len(rows)
            elapsed = time.time() - t0
            rate = completed / elapsed if elapsed > 0 else 0
            eta_min = (len(remaining) - completed) / rate / 60 if rate > 0 else float('inf')
            print(f"[{completed}/{len(remaining)}] {key}: {total} rows, {len(rows)} with RV "
                  f"(running kept total: {total_kept}) -- {rate*3600:.1f} chunks/hr, ETA {eta_min:.0f} min",
                  flush=True)

    outf.close()
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
