// Byte-range row iteration for splitting one streaming pass over the 43 GB Gaia CSV across several
// processes (src/gaia_scan.c --byte-range=START:END, added 2026-10-08). A row belongs to the range
// that contains its FIRST byte, so adjacent ranges [a,b) [b,c) cover every row exactly once,
// whatever the split offsets are (tests/test_csv_range.c checks every split point of a small file).
// The FILE must be opened in binary mode ("rb"): on Windows a text-mode ftell is not a byte offset.
#ifndef NOSTOS_CSV_RANGE_H
#define NOSTOS_CSV_RANGE_H

#include <stdio.h>
#include <stdbool.h>

// Positions f at the first row starting at offset >= start. start <= 0 is a no-op (the stream stays
// where it is, e.g. just after a header the caller has already read; a `start` inside the header
// line also lands on the first data row). Uses buf (size n) as scratch. Returns false on a seek
// error or if nothing follows offset start-1.
bool csv_range_begin(FILE *f, long long start, char *buf, int n);

// Reads the next row into buf if it starts before `end` (end < 0: no limit). Returns false at the
// end of the range or the file.
bool csv_range_next(FILE *f, long long end, char *buf, int n);

#endif
