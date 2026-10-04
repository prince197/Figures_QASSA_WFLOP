#!/usr/bin/env python3
"""Generate the full manuscript from canonical modular sources. No scientific edits."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]

def expand(path, stack=()):
    path = path.resolve()
    if path in stack:
        raise ValueError(f"Cyclic input: {path}")
    body = path.read_text()
    def include(match):
        rel = match.group(1)
        child = ROOT / rel
        if not child.suffix:
            child = child.with_suffix('.tex')
        if not child.is_file():
            raise FileNotFoundError(child)
        return f"% begin {rel}\n" + expand(child, stack+(path,)) + f"\n% end {rel}\n"
    # Expand file-level includes only, leaving parameterized macro bodies intact.
    return re.sub(r'^\\input\{([^}#]+)\}', include, body, flags=re.M)

if __name__ == '__main__':
    out = ROOT/'SWEVO_manuscript_full.tex'
    out.write_text('% Generated from SWEVO_manuscript.tex. Edit modular sources, then regenerate.\n'+expand(ROOT/'SWEVO_manuscript.tex'))
    print(f"Generated {out.name}")
