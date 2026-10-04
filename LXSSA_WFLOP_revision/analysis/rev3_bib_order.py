"""Order a manual thebibliography by first citation (revision 3, decision D13).

Usage:
    python3 analysis/rev3_bib_order.py <manuscript .aux> <bib .tex file> [--check] [--drop-uncited]

<manuscript .aux>  the .aux of a complete compilation of the document (e.g. SWEVO_manuscript.aux); the order of its
                   \\citation{...} lines (following \\@input{...} sub-aux files in place) is the order in which the
                   citations appear in the source.
<bib .tex file>    the file holding the \\begin{thebibliography}...\\end{thebibliography} (e.g. optA/swevo_back.tex).

The entries are rewritten in order of first citation; everything outside the list (declarations, wrapper such as
{\\footnotesize ...}, text before the first \\bibitem) and the text of every entry, including its % comment lines, is
kept unchanged. Comment lines placed directly above a \\bibitem (no blank line between) move with that entry.

Reported:
  * cited but missing: keys cited in the .aux that have no \\bibitem (the script cannot create them; exit status 2);
  * uncited: entries never cited. They are kept at the end of the list in their present order (Elsevier requires
    every listed reference to be cited, so remove them or cite them), or dropped with --drop-uncited.
--check only reports (exit status 1 if the file is not in citation order, 2 if keys are missing) and writes nothing.
The file is written only when its content changes. Rerunning on an ordered file changes nothing.
"""
import os
import re
import sys


def citation_order(aux_path, seen=None):
    """Keys in order of first \\citation in the .aux, following \\@input sub-aux files."""
    order = [] if seen is None else seen
    base = os.path.dirname(os.path.abspath(aux_path))
    with open(aux_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            m = re.match(r'\\@input\{([^}]+)\}', line)
            if m:
                sub = os.path.join(base, m.group(1))
                if os.path.exists(sub):
                    citation_order(sub, order)
                continue
            for c in re.findall(r'\\citation\{([^}]*)\}', line):
                for k in c.split(','):
                    k = k.strip()
                    if k and k != '*' and k not in order:
                        order.append(k)
    return order


def split_entries(body):
    """Split the body of a thebibliography into (preamble, [(key, block)]).

    A block starts at the comment lines directly above its \\bibitem (if any) and ends before the next block."""
    lines = body.splitlines(keepends=True)
    starts = []  # (start line index, key)
    for i, ln in enumerate(lines):
        m = re.match(r'\s*\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}', ln)
        if m:
            j = i
            while j > 0 and lines[j - 1].lstrip().startswith('%'):
                j -= 1
            starts.append((j, m.group(1).strip()))
    if not starts:
        return body, []
    pre = ''.join(lines[:starts[0][0]])
    ents = []
    for n, (j, key) in enumerate(starts):
        end = starts[n + 1][0] if n + 1 < len(starts) else len(lines)
        ents.append((key, ''.join(lines[j:end])))
    return pre, ents


def main(argv):
    args = [a for a in argv if not a.startswith('--')]
    flags = {a for a in argv if a.startswith('--')}
    unknown = flags - {'--check', '--drop-uncited'}
    if len(args) != 2 or unknown:
        sys.exit(__doc__)
    aux, bib = args
    src = open(bib, encoding='utf-8').read()
    m = re.search(r'(\\begin\{thebibliography\}\{[^}]*\}[^\n]*\n)(.*?)(\\end\{thebibliography\})', src, re.S)
    if not m:
        sys.exit(f'no thebibliography environment in {bib}')
    pre, ents = split_entries(m.group(2))
    keys = [k for k, _ in ents]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    if dup:
        sys.exit('duplicate \\bibitem keys: ' + ', '.join(dup))
    blocks = dict(ents)
    order = citation_order(aux)
    missing = [k for k in order if k not in blocks]
    cited = [k for k in order if k in blocks]
    uncited = [k for k in keys if k not in order]

    print(f'{bib}: {len(keys)} entries; {len(order)} distinct keys cited in {aux}')
    print('cited but missing (no \\bibitem): ' + (', '.join(missing) if missing else 'none'))
    print('uncited entries: ' + (', '.join(uncited) if uncited else 'none')
          + (' (dropped)' if uncited and '--drop-uncited' in flags else ' (kept at the end)' if uncited else ''))
    new_keys = cited + ([] if '--drop-uncited' in flags else uncited)
    moved = sum(1 for a, b in zip(keys, new_keys) if a != b)
    print(f'entries not at their citation-order position: {moved}')

    def norm(block):  # one blank line after every entry
        return block.rstrip('\n') + '\n\n'
    body = pre + ''.join(norm(blocks[k]) for k in new_keys)
    body = body.rstrip('\n') + '\n\n'
    new = src[:m.start(2)] + body + src[m.end(2):]
    if '--check' in flags:
        sys.exit(2 if missing else (1 if new_keys != keys else 0))
    if new != src:
        open(bib, 'w', encoding='utf-8').write(new)
        print(f'written: {bib}')
    else:
        print('unchanged')
    for i, k in enumerate(new_keys, 1):
        print(f'  [{i}] {k}')
    sys.exit(2 if missing else 0)


if __name__ == '__main__':
    main(sys.argv[1:])
