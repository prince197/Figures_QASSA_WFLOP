"""Strip LaTeX comments from every .tex file under a folder (in place). A full-line comment is removed; a trailing
comment keeps its '%' (so end-of-line whitespace suppression is unchanged). Escaped \\% is respected."""
import re, sys, pathlib
root = pathlib.Path(sys.argv[1])
def first_comment(line):
    i = 0
    while True:
        j = line.find('%', i)
        if j < 0: return -1
        k = j - 1; bs = 0
        while k >= 0 and line[k] == '\\': bs += 1; k -= 1
        if bs % 2 == 0: return j
        i = j + 1
nfile = nline = 0
for f in root.rglob('*.tex'):
    out = []; changed = False
    for line in f.read_text(encoding='utf-8').splitlines(keepends=True):
        j = first_comment(line)
        if j < 0: out.append(line); continue
        changed = True; nline += 1
        if line[:j].strip() == '':
            continue                                   # full-line comment: drop the line
        out.append(line[:j + 1] + '\n')                # keep code and the bare %
    if changed:
        nfile += 1; f.write_text(''.join(out), encoding='utf-8')
print(f'stripped comments: {nline} lines in {nfile} files')
