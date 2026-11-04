#!/usr/bin/env python3
"""Summary of a run_bench.sh results.tsv: per label the median, minimum and maximum wall time, CPU time, memory, GPU
use, speed relative to the first label, the verdict of the acceptance rule below, and whether the label's results are
the same as the first label's (or as an earlier label's: backends and compilers differ at rounding level, so a CUDA
build is compared with the CUDA build before it).

Acceptance rule (as in docs/BENCHMARK.md): a candidate is faster when its median is at least 5 % below the first
label's and its slowest run is faster than the first label's fastest, or when its median is at least 10 % below. With
fewer than 3 runs of either label the verdict is "too few runs".

The last column of results.tsv is the smallest min_pc of the run's `out` (the closest encounter of the work unit); it is
only used to say how far apart two different outputs are.

usage: summarize_bench.py <results.tsv>
"""
import statistics, sys


def main(path):
    runs, order = {}, []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        label = f[0]
        if label not in runs:
            runs[label] = []
            order.append(label)
        runs[label].append(dict(wall=float(f[2] or "nan"), user=float(f[3] or "nan"), sys=float(f[4] or "nan"),
                                rss=f[5], code=f[6], hash=f[7], gpu=f[8], best=f[9]))
    if not order:
        print("no runs")
        return
    ref = order[0]
    ref_hashes = {r["hash"] for r in runs[ref]}
    ref_best = [float(r["best"]) for r in runs[ref] if r["best"] not in ("-", "")]
    rmed = statistics.median(r["wall"] for r in runs[ref])
    rmin = min(r["wall"] for r in runs[ref])
    print("%-12s %4s %9s %8s %8s %8s %7s %5s  %7s  %s" % ("label", "runs", "median s", "min s", "max s", "cpu s",
                                                         "RSS MB", "GPU%", "vs 1st", "verdict / output"))
    for lb in order:
        rs = runs[lb]
        walls = [r["wall"] for r in rs]
        med, lo, hi = statistics.median(walls), min(walls), max(walls)
        cpu = statistics.median(r["user"] + r["sys"] for r in rs)
        rss = statistics.median(float(r["rss"] or 0) for r in rs)
        gpus = [float(r["gpu"]) for r in rs if r["gpu"] not in ("-", "")]
        gpu = "%.0f" % statistics.median(gpus) if gpus else "-"
        bad = [r["code"] for r in rs if r["code"] != "0"]
        # run_bench.sh writes "none" for a run without out and "incomplete" for an out without its "# done N" trailer
        unfinished = [r["hash"] for r in rs if r["hash"] in ("none", "incomplete")]
        hashes = {r["hash"] for r in rs}
        if bad:
            out = "FAILED (exit %s)" % ",".join(sorted(set(bad)))
        elif unfinished:
            out = "FAILED (out %s)" % ",".join(sorted(set(unfinished)))
        elif len(hashes) > 1:
            out = "results differ between its own runs"
        elif hashes == ref_hashes:
            out = "same output"
        elif any({r["hash"] for r in runs[o]} == hashes for o in order[:order.index(lb)]):
            # e.g. a CUDA build after the CPU builds: compare like with like
            out = "same output as " + next(o for o in order if {r["hash"] for r in runs[o]} == hashes)
        else:
            c = [float(r["best"]) for r in rs if r["best"] not in ("-", "")]
            d = (statistics.median(c) - statistics.median(ref_best)) if c and ref_best else float("nan")
            out = "OTHER output (closest min_pc %+.6f)" % d
        if lb == ref:
            verdict = "baseline"
        else:
            gain = 1 - med / rmed
            if min(len(rs), len(runs[ref])) < 3:
                # one or two runs say nothing about the spread, so the rule cannot be applied
                verdict = "too few runs (%+.1f %%)" % (100 * gain)
            elif (gain >= 0.05 and hi < rmin) or gain >= 0.10:
                verdict = "ACCEPT (%+.1f %%)" % (100 * gain)
            elif gain <= -0.05:
                verdict = "slower (%+.1f %%)" % (100 * gain)
            else:
                verdict = "no clear gain (%+.1f %%)" % (100 * gain)
        print("%-12s %4d %9.2f %8.2f %8.2f %8.1f %7.0f %5s  %6.3fx  %s; %s" % (lb, len(rs), med, lo, hi, cpu, rss, gpu,
                                                                                rmed / med, verdict, out))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    main(sys.argv[1])
