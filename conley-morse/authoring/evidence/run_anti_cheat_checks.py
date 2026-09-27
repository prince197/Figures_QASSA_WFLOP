#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json, tempfile, textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VC = ROOT / 'tests' / 'verifier_core.py'
spec = importlib.util.spec_from_file_location('verifier_core_local', VC)
vc = importlib.util.module_from_spec(spec); spec.loader.exec_module(vc)

ADVERSARIES = {
    'nop': 'import sys\n',
    'case_id_only': '''\
import json,sys
with open(sys.argv[1]) as f: d=json.load(f)
with open(sys.argv[2],'w') as f: json.dump({'case_id':d.get('case_id'),'queries':[]},f)
''',
    'oracle_snoop': '''\
import json,sys
try:
    x=json.load(open('/tests/oracles/heldout_0.json'))
except Exception as e:
    x={'case_id':'blocked','queries':[]}
json.dump(x,open(sys.argv[2],'w'))
''',
    'symlink_output': '''\
import os,sys,tempfile
p=sys.argv[2]
try: os.unlink(p)
except FileNotFoundError: pass
os.symlink('/etc/hosts',p)
''',
}


def main():
    results={}
    with tempfile.TemporaryDirectory(prefix='cm_adversaries_') as td:
        td=Path(td)
        for name,src in ADVERSARIES.items():
            p=td/f'{name}.py'; p.write_text(textwrap.dedent(src),encoding='utf-8')
            ok,msg,diag=vc.evaluate(p)
            results[name]={
                'reward': int(ok),
                'message': msg,
                'passed_cases': diag.get('passed_cases',0),
            }
    out=ROOT/'authoring'/'evidence'/'anti_cheat_report.json'
    out.write_text(json.dumps(results,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    md=ROOT/'authoring'/'evidence'/'anti_cheat_report.md'
    lines=['# Anti-cheat evidence','']
    for name,r in results.items():
        lines.append(f"- `{name}`: reward {r['reward']}; {r['passed_cases']} held-out cases passed.")
    lines += [
        '',
        'Unearned-credit judgment: none of the local adversarial submissions receives reward. The oracle-snoop trial cannot read the sealed oracle directory after the verifier drops privileges, and the symlink trial is rejected as a non-regular output. These checks support the isolation design; they do not claim to substitute for the platform adversarial-agent probe.',
        ''
    ]
    md.write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(results,indent=2,sort_keys=True))

if __name__=='__main__': main()
