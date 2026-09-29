"""Merge new references into optA/swevo_back.tex and order the bibliography by first citation in the SWEVO main text.

Usage: python3 analysis/swevo_bib.py <SWEVO_manuscript.aux> [extra .tex files with \\bibitem entries ...]
Entries not cited in the main text are reported and dropped (Elsevier: every reference must be cited).
"""
import re, sys

BACK = 'optA/swevo_back.tex'

def entries(text):
    """Split a thebibliography body into {key: block}; a block runs from its \\bibitem to the next one."""
    parts = re.split(r'(?=^\\bibitem\{)', text, flags=re.M)
    out = {}
    for p in parts:
        m = re.match(r'\\bibitem\{([^}]+)\}', p)
        if m and m.group(1) not in out:
            out[m.group(1)] = p.rstrip('\n') + '\n'
    return out

src = open(BACK).read()
m = re.search(r'(\\begin\{thebibliography\}\{\d+\}\n)(.*?)(\\end\{thebibliography\})', src, re.S)
head, body, tail = src[:m.start()], m.group(2), src[m.end():]
# a previous run wrapped the list in {\footnotesize ...}: remove that wrapper before adding it again
head = re.sub(r'\{\\footnotesize[^\n]*\n$', '', head)
tail = re.sub(r'^\n\}', '', tail)
pre = body[:body.find('\\bibitem')]
ents = entries(body)
for f in sys.argv[2:]:
    for k, v in entries(open(f).read()).items():
        ents.setdefault(k, v)
order = []
for c in re.findall(r'\\citation\{([^}]+)\}', open(sys.argv[1]).read()):
    for k in c.split(','):
        k = k.strip()
        if k not in order:
            order.append(k)
missing = [k for k in order if k not in ents]
uncited = [k for k in ents if k not in order]
if missing:
    sys.exit('cited but no entry: ' + ', '.join(missing))
print('entries', len(order), '; dropped (uncited):', ', '.join(uncited) or 'none')
blocks = [ents[k] if ents[k].endswith('\n\n') else ents[k] + '\n' for k in order]
new = head + '{\\footnotesize\\setlength{\\bibsep}{1pt plus 0.3ex}\n\\begin{thebibliography}{99}\n' + pre + ''.join(blocks).rstrip('\n') + '\n\n\\end{thebibliography}\n}' + tail
open(BACK, 'w').write(new)
