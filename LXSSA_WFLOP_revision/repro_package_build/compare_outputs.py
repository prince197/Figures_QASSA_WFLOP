#!/usr/bin/env python3
"""Compare regenerated files with the stored copies.

Usage:  python3 compare_outputs.py STORED_DIR REGEN_DIR FILE [FILE ...]   (FILE may be a glob, relative to both dirs)

Each file is reported as
  IDENTICAL            byte-identical,
  IDENTICAL-NORMALIZED identical after removing only volatile content: generation time stamps
                       (YYYY-MM-DD HH:MM:SS), the PDF /CreationDate, the gzip header time, absolute paths of the
                       machine, run-time statements ("in 12 s", "12.3 s", JSON "seconds": 27.2) and the
                       environment version strings recorded in JSON outputs ("python": "3.11.15", "numpy": ...,
                       "pandas": ...); the normalized lines are counted,
  DIFFERENT            anything else (the first differing lines are printed),
  MISSING              absent in one of the two directories.
Exit code 1 if any file is DIFFERENT or MISSING.
"""
import difflib, glob, gzip, os, re, sys

TS = re.compile(rb"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}")
PDFDATE = re.compile(rb"/CreationDate \(D:[^)]*\)")
PATH = re.compile(rb"(/[A-Za-z0-9_.\-]+)+/(analysis|figures_mpce|rev2_out|out)\b")
SECS = re.compile(rb"\b(in |done in |elapsed )?\d+(\.\d+)? ?s\b")
JSECS = re.compile(rb'"(seconds|elapsed_s|runtime_s)": [0-9.eE+-]+')
JVERS = re.compile(rb'"(python|numpy|scipy|pandas|matplotlib)": "[0-9][^"]*"')


def norm(b):
    b = PDFDATE.sub(b"/CreationDate ()", b)
    b = TS.sub(b"<TIME>", b)
    b = PATH.sub(b"<PATH>", b)
    b = SECS.sub(b"<SECONDS>", b)
    b = JSECS.sub(rb'"\1": <SECONDS>', b)
    b = JVERS.sub(rb'"\1": "<VERSION>"', b)
    return b


def main(argv):
    if len(argv) < 3:
        print(__doc__); return 2
    sd, rd, pats = argv[0], argv[1], argv[2:]
    files = []
    for p in pats:
        hits = sorted(set(os.path.relpath(f, sd) for f in glob.glob(os.path.join(sd, p))) |
                      set(os.path.relpath(f, rd) for f in glob.glob(os.path.join(rd, p))))
        files += hits or [p]
    bad = 0
    for f in files:
        a, b = os.path.join(sd, f), os.path.join(rd, f)
        if not (os.path.exists(a) and os.path.exists(b)):
            print(f"MISSING              {f} (stored: {os.path.exists(a)}, regenerated: {os.path.exists(b)})"); bad += 1; continue
        x, y = open(a, "rb").read(), open(b, "rb").read()
        if x == y:
            print(f"IDENTICAL            {f}"); continue
        if f.endswith(".gz"):  # the gzip header holds the write time: compare the decompressed content
            x, y = gzip.decompress(x), gzip.decompress(y)
            if x == y:
                print(f"IDENTICAL-NORMALIZED {f} (decompressed content byte-identical; gzip header time differs)"); continue
        nx, ny = norm(x), norm(y)
        if nx == ny:
            n = sum(1 for u, v in zip(x.split(b"\n"), y.split(b"\n")) if u != v)
            print(f"IDENTICAL-NORMALIZED {f} ({n} line(s) differ only in time stamps / paths / run times / version strings)"); continue
        bad += 1
        lx, ly = nx.decode("utf-8", "replace").splitlines(), ny.decode("utf-8", "replace").splitlines()
        d = [l for l in difflib.unified_diff(lx, ly, "stored", "regenerated", n=0, lineterm="") if not l.startswith("@@")]
        print(f"DIFFERENT            {f} ({sum(1 for l in d if l[:1] in '+-') - 2} changed lines)")
        for l in d[2:12]:
            print("      " + l[:200])
    print(f"{len(files) - bad} of {len(files)} files identical (byte-identical or after normalization)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
