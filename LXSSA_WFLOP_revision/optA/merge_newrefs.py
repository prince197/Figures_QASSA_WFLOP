#!/usr/bin/env python3
"""Merge new \\bibitem entries into optA/11_back.tex and put the bibliography in first-citation order.

Usage (from LXSSA_WFLOP_revision/):
    python3 optA/merge_newrefs.py            # dry run: report only, nothing is written
    python3 optA/merge_newrefs.py --apply    # rewrite optA/11_back.tex
    python3 optA/merge_newrefs.py --out F    # write the merged file to F instead (for testing)

Steps
1. Read the \\bibitem entries of the thebibliography environment in optA/11_back.tex.
2. Add every \\bibitem found in optA/newrefs/*.tex (files in name order); a key that already exists is
   skipped (the existing entry wins; a differing text is reported).
3. Collect \\cite / \\nocite keys in reading order: MPCE_PSO_VNS.tex, following its \\input / \\include commands
   recursively (optA/01..11 in \\input order, and the analysis/*.tex tables they input); comments are ignored,
   and the thebibliography environment itself is not scanned.
4. Write all cited entries in first-citation order and drop entries that are no longer cited (printed).
   Keys cited but without an entry are printed as errors (the file is then not written unless --force).
Each entry keeps its lines verbatim (including % VERIFY-REF comments after the \\bibitem line).
"""
import argparse
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # LXSSA_WFLOP_revision
MAIN = os.path.join(ROOT, "MPCE_PSO_VNS.tex")
BACK = os.path.join(HERE, "11_back.tex")
NEWREFS = os.path.join(HERE, "newrefs")

BIBITEM = re.compile(r"\\bibitem\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}")
CITE = re.compile(r"\\(?:no)?cite[a-zA-Z]*\*?\s*(?:\[[^\]]*\]\s*)*\{([^}]*)\}")
INPUT = re.compile(r"\\(?:input|include)\s*\{([^}]+)\}")
BEGIN_BIB = re.compile(r"\\begin\{thebibliography\}\{[^}]*\}")
END_BIB = "\\end{thebibliography}"


def strip_comments(text):
    """Remove TeX comments (an unescaped % to the end of the line)."""
    return "\n".join(re.sub(r"(?<!\\)%.*", "", line) for line in text.split("\n"))


def split_entries(body):
    """Split a bibliography body into (key, entry_text) in order; text before the first \\bibitem is ignored."""
    starts = [m for m in BIBITEM.finditer(body) if not _in_comment(body, m.start())]
    entries = []
    for i, m in enumerate(starts):
        end = starts[i + 1].start() if i + 1 < len(starts) else len(body)
        entries.append((m.group(1).strip(), body[m.start():end].strip()))
    return entries


def _in_comment(text, pos):
    line_start = text.rfind("\n", 0, pos) + 1
    return re.search(r"(?<!\\)%", text[line_start:pos]) is not None


def resolve(name):
    for cand in (name, name + ".tex"):
        p = os.path.join(ROOT, cand)
        if os.path.isfile(p):
            return p
    return None


def cited_keys(path, seen_files, out, missing_inputs):
    """Append cited keys of path (recursively through \\input) to out, in reading order."""
    if path in seen_files:
        return
    seen_files.add(path)
    text = strip_comments(open(path, encoding="utf-8").read())
    # do not scan the bibliography itself
    m = BEGIN_BIB.search(text)
    if m:
        e = text.find(END_BIB, m.end())
        text = text[:m.start()] + (text[e + len(END_BIB):] if e >= 0 else "")
    tokens = sorted([(m.start(), "cite", m) for m in CITE.finditer(text)] +
                    [(m.start(), "input", m) for m in INPUT.finditer(text)], key=lambda t: t[0])
    for _, kind, m in tokens:
        if kind == "cite":
            for k in m.group(1).split(","):
                k = k.strip()
                if k:
                    out.append(k)
        else:
            p = resolve(m.group(1).strip())
            if p is None:
                missing_inputs.append(m.group(1))
            else:
                cited_keys(p, seen_files, out, missing_inputs)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--apply", action="store_true", help="rewrite optA/11_back.tex")
    ap.add_argument("--out", help="write the result to this file instead of optA/11_back.tex")
    ap.add_argument("--force", action="store_true", help="write even if cited keys have no entry")
    args = ap.parse_args()

    back = open(BACK, encoding="utf-8").read()
    mb = BEGIN_BIB.search(back)
    me = back.find(END_BIB, mb.end()) if mb else -1
    if not mb or me < 0:
        sys.exit("thebibliography environment not found in optA/11_back.tex")
    head, body, tail = back[:mb.start()], back[mb.end():me], back[me + len(END_BIB):]

    entries = {}
    old_order = []
    for k, t in split_entries(body):
        if k in entries:
            print(f"WARNING: duplicate key {k} in 11_back.tex; first entry kept")
            continue
        entries[k] = t
        old_order.append(k)

    added, skipped = [], []
    for f in sorted(glob.glob(os.path.join(NEWREFS, "*.tex"))):
        for k, t in split_entries(open(f, encoding="utf-8").read()):
            src = os.path.relpath(f, ROOT)
            if k in entries:
                same = " ".join(entries[k].split()) == " ".join(t.split())
                skipped.append((k, src, same))
                continue
            entries[k] = t
            added.append((k, src))

    cites, missing_inputs = [], []
    cited_keys(MAIN, set(), cites, missing_inputs)
    order = list(dict.fromkeys(cites))
    undefined = [k for k in order if k not in entries]
    kept = [k for k in order if k in entries]
    dropped = [k for k in list(entries) if k not in set(order)]

    print(f"entries in 11_back.tex: {len(old_order)}; added from newrefs: {len(added)}; "
          f"cited keys: {len(order)}; output entries: {len(kept)}")
    for k, src in added:
        print(f"  added   {k:24s} from {src}")
    for k, src, same in skipped:
        print(f"  skipped {k:24s} from {src} (key exists{'' if same else '; TEXT DIFFERS, existing entry kept'})")
    for k in dropped:
        print(f"  DROPPED {k} (no longer cited)")
    for k in undefined:
        print(f"  ERROR   {k} is cited but has no \\bibitem")
    for n in missing_inputs:
        print(f"  note: \\input{{{n}}} not found, skipped")
    old_kept = [k for k in old_order if k in kept]
    new_old = [k for k in kept if k in old_order]
    print("order of existing entries:", "unchanged" if old_kept == new_old else "changed")
    if old_kept != new_old:
        for i, k in enumerate(kept, 1):
            was = old_order.index(k) + 1 if k in old_order else None
            if was != i:
                print(f"  [{i}] {k}" + (f" (was [{was}])" if was else " (new)"))

    if not (args.apply or args.out):
        print("dry run: nothing written (use --apply to rewrite optA/11_back.tex)")
        return
    if undefined and not args.force:
        sys.exit("not written: cited keys without \\bibitem (use --force to write anyway)")
    width = "99" if len(kept) < 100 else "999"
    new = (head + "\\begin{thebibliography}{" + width + "}\n\n" +
           "\n\n".join(entries[k] for k in kept) + "\n\n" + END_BIB + tail)
    dest = args.out or BACK
    with open(dest, "w", encoding="utf-8") as fh:
        fh.write(new)
    print("written:", os.path.relpath(dest, ROOT) if dest.startswith(ROOT) else dest)


if __name__ == "__main__":
    main()
