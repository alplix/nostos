"""List every GaiaSource_*.csv.gz chunk in the bulk download bucket (paginated),
and report total count + total size."""
import requests
import xml.etree.ElementTree as ET

BASE = "https://gaia.eu-1.cdn77-storage.com/"
PREFIX = "Gaia/gdr3/gaia_source/"
NS = {"s3": "http://s3.amazonaws.com/doc/2006-03-01/"}

files = []
marker = None
while True:
    params = {"prefix": PREFIX, "delimiter": "/"}
    if marker:
        params["marker"] = marker
    r = requests.get(BASE, params=params, timeout=60)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    for c in root.findall("s3:Contents", NS):
        key = c.find("s3:Key", NS).text
        size = int(c.find("s3:Size", NS).text)
        files.append((key, size))
    truncated = root.find("s3:IsTruncated", NS).text == "true"
    if not truncated:
        break
    marker = root.find("s3:NextMarker", NS).text

total_size = sum(s for _, s in files)
print(f"Total files: {len(files)}")
print(f"Total size: {total_size / 1e9:.1f} GB")
print(f"First file: {files[0]}")
print(f"Last file: {files[-1]}")

with open("gaia_file_list.txt", "w") as f:
    for key, size in files:
        f.write(f"{key}\t{size}\n")
print("Wrote gaia_file_list.txt")
