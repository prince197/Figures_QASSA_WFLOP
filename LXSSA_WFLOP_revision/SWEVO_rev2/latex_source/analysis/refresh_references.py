#!/usr/bin/env python3
"""Refresh portable external labels after compilation, without parsing scientific data."""
from pathlib import Path
import re
ROOT = Path(__file__).resolve().parents[1]
for aux, prefix, target in [('SWEVO_manuscript.aux','M-','manuscript_labels.tex'),('SWEVO_supplement.aux','S-','supplement_labels.tex')]:
    source = ROOT/aux
    if not source.exists():
        continue
    labels=[]
    for line in source.read_text().splitlines():
        m=re.match(r'\\newlabel\{([^}]+)\}(.*)',line)
        if m and not m.group(1).startswith(('M-','S-')):
            labels.append('\\newlabel{'+prefix+m.group(1)+'}'+m.group(2))
    (ROOT/'xref'/target).write_text('\\makeatletter\n'+'\n'.join(labels)+'\n\\makeatother\n')
    print(f"Refreshed {len(labels)} labels in {target}")
