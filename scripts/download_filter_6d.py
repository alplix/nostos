"""Download every GaiaSource_*.csv.gz bulk chunk, stream-filter to only the
rows with a measured radial_velocity (the full 6D phase-space subset this
project needs), append those rows to a single growing output CSV, then
delete the downloaded chunk before moving to the next. Resumable: tracks
completed chunks in progress.txt so an interrupted run can pick back up
without re-downloading or re-counting anything already done.

~753GB passes through temporarily (one ~230MB chunk at a time), but nothing
but the filtered output (expected far smaller -- only ~33M of ~1.8B stars
have a measured radial velocity) is kept.
"""
import csv
import gzip
import io
import os
import sys
import time
import requests

BASE = "https://gaia.eu-1.cdn77-storage.com/"
LIST_FILE = "gaia_file_list.txt"
OUT_FILE = "gaia_dr3_6d.csv"
PROGRESS_FILE = "progress.txt"
TMP_FILE = "_chunk_tmp.csv.gz"

RV_COL = "radial_velocity"

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
    with open(PROGRESS_FILE, "a") as f:
        f.write(key + "\n")

def download(key, dest):
    url = BASE + key
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)

def process_one(key, writer, header_written):
    download(key, TMP_FILE)
    kept = 0
    total = 0
    with gzip.open(TMP_FILE, "rt", newline="") as gz:
        # Gaia's bulk CSVs are actually ECSV (astropy): ~1000 lines of "#"-prefixed
        # YAML column metadata precede the real CSV header line. Skip past those.
        line = gz.readline()
        while line.startswith("#"):
            line = gz.readline()
        reader = csv.reader([line] + list(gz))
        header = next(reader)
        rv_idx = header.index(RV_COL)
        if not header_written[0]:
            writer.writerow(header)
            header_written[0] = True
        for row in reader:
            total += 1
            # Nulls are the literal string "null", not an empty field.
            if row[rv_idx] != "null" and row[rv_idx] != "":
                writer.writerow(row)
                kept += 1
    os.remove(TMP_FILE)
    return total, kept

def main():
    files = load_file_list()
    done = load_progress()
    remaining = [(k, s) for k, s in files if k not in done]
    print(f"{len(files)} total chunks, {len(done)} already done, {len(remaining)} remaining", flush=True)

    header_written = [os.path.exists(OUT_FILE) and os.path.getsize(OUT_FILE) > 0]
    out_mode = "a" if header_written[0] else "w"

    with open(OUT_FILE, out_mode, newline="") as outf:
        writer = csv.writer(outf)
        t0 = time.time()
        total_kept = 0
        for i, (key, size) in enumerate(remaining):
            attempt = 0
            while True:
                try:
                    total, kept = process_one(key, writer, header_written)
                    outf.flush()
                    break
                except Exception as e:
                    attempt += 1
                    print(f"  retry {attempt} for {key}: {e}", flush=True)
                    if os.path.exists(TMP_FILE):
                        os.remove(TMP_FILE)
                    if attempt >= 5:
                        print(f"  GIVING UP on {key} after 5 attempts", flush=True)
                        total, kept = 0, 0
                        break
                    time.sleep(5)
            total_kept += kept
            mark_done(key)
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta_min = (len(remaining) - i - 1) / rate / 60 if rate > 0 else float('inf')
            print(f"[{i+1}/{len(remaining)}] {key}: {total} rows, {kept} with RV "
                  f"(running kept total: {total_kept}) -- {rate*3600:.1f} chunks/hr, ETA {eta_min:.0f} min",
                  flush=True)

    print("ALL DONE", flush=True)

if __name__ == "__main__":
    main()
