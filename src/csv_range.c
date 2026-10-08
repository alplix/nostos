#if !defined(_WIN32)
#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200112L
#endif
#include "csv_range.h"

#if defined(_WIN32)
#define NOSTOS_FSEEK _fseeki64
#define NOSTOS_FTELL _ftelli64
typedef long long nostos_off_t;
#else
#include <sys/types.h>
#define NOSTOS_FSEEK fseeko
#define NOSTOS_FTELL ftello
typedef off_t nostos_off_t;
#endif

bool csv_range_begin(FILE *f, long long start, char *buf, int n) {
    if (start <= 0) return true; // keep the current position (e.g. just after a header)
    // Seek one byte before `start` and discard the rest of that line: if byte start-1 is the
    // newline ending the previous row, only that newline is discarded and the row beginning exactly
    // at `start` is kept; otherwise the partial row (owned by the previous range) is skipped.
    if (NOSTOS_FSEEK(f, (nostos_off_t)(start - 1), SEEK_SET) != 0) return false;
    return fgets(buf, n, f) != NULL;
}

bool csv_range_next(FILE *f, long long end, char *buf, int n) {
    if (end >= 0 && (long long)NOSTOS_FTELL(f) >= end) return false;
    return fgets(buf, n, f) != NULL;
}
