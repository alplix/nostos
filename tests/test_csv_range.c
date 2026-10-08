// Checks src/csv_range.c (gaia_scan.c --byte-range) exhaustively on a small file: for EVERY split
// offset s, reading [header_end, s) and then [s, EOF) must return every data row exactly once and
// in order -- including splits inside the header, exactly at a row start, on a newline, and at EOF
// -- for both LF and CRLF line endings (the real catalogue CSV is CRLF). Also a three-way split.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/csv_range.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else printf("  ok:   %s\n", what);
}

static const char *ROWS[] = { "1,a,0.5", "22,bb,1.25", "333,c,", "4,dddddddd,7", "5,e,8e-3", "66,f,9" };
#define N_ROWS 6

// Reads one range into out (concatenated row texts with line endings stripped, '|'-separated).
static void read_range(FILE *f, long long hdr_end, long long start, long long end, char *out) {
    char buf[256];
    out[0] = '\0';
    fseek(f, 0, SEEK_SET);
    if (!fgets(buf, sizeof buf, f)) return; // header
    (void)hdr_end;
    if (!csv_range_begin(f, start, buf, sizeof buf)) return;
    while (csv_range_next(f, end, buf, sizeof buf)) {
        size_t L = strlen(buf);
        while (L && (buf[L - 1] == '\n' || buf[L - 1] == '\r')) buf[--L] = '\0';
        strcat(out, buf);
        strcat(out, "|");
    }
}

static void run(const char *eol) {
    const char *path = "build/test_csv_range_tmp.csv";
    FILE *w = fopen(path, "wb");
    fprintf(w, "id,name,value%s", eol);
    long long hdr_end = (long long)ftell(w);
    char expect[512] = "";
    for (int i = 0; i < N_ROWS; i++) {
        fprintf(w, "%s%s", ROWS[i], eol);
        strcat(expect, ROWS[i]);
        strcat(expect, "|");
    }
    long long size = (long long)ftell(w);
    fclose(w);

    FILE *f = fopen(path, "rb");
    int bad = 0;
    char a[512], b[512], c[512], joined[1536];
    for (long long s = 0; s <= size; s++) {
        read_range(f, hdr_end, 0, s, a);
        read_range(f, hdr_end, s, -1, b);
        snprintf(joined, sizeof joined, "%s%s", a, b);
        if (strcmp(joined, expect) != 0) {
            if (bad < 3) printf("  split at %lld: got '%s' expected '%s'\n", s, joined, expect);
            bad++;
        }
    }
    char msg[128];
    snprintf(msg, sizeof msg, "%s: every one of %lld two-way splits returns each row exactly once",
             eol[0] == '\r' ? "CRLF" : "LF", size + 1);
    check(bad == 0, msg);

    int bad3 = 0;
    for (long long s1 = 0; s1 <= size; s1 += 3)
        for (long long s2 = s1; s2 <= size; s2 += 5) {
            if (s2 == s1) continue;
            read_range(f, hdr_end, 0, s1, a);
            read_range(f, hdr_end, s1, s2, b);
            read_range(f, hdr_end, s2, -1, c);
            snprintf(joined, sizeof joined, "%s%s%s", a, b, c);
            if (strcmp(joined, expect) != 0) bad3++;
        }
    snprintf(msg, sizeof msg, "%s: three-way splits too", eol[0] == '\r' ? "CRLF" : "LF");
    check(bad3 == 0, msg);
    fclose(f);
    remove(path);
}

int main(void) {
    printf("test_csv_range:\n");
    run("\n");
    run("\r\n");
    if (g_failures) { printf("test_csv_range: %d FAILURE(S)\n", g_failures); return 1; }
    printf("test_csv_range: all passed\n");
    return 0;
}
