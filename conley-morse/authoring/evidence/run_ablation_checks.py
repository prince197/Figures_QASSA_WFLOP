#!/usr/bin/env python3
from __future__ import annotations
import copy, importlib.util, json, subprocess, sys, tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('vc',ROOT/'tests'/'verifier_core.py')
vc=importlib.util.module_from_spec(spec); spec.loader.exec_module(vc)
ORACLE_DIR=ROOT/'tests'/'oracles'


def greedy_single_branch(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['selected_morse_ids_by_frame']=[s[:1] for s in q['selected_morse_ids_by_frame']]
    return x


def endpoint_dimension_as_rank(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        dims=q['conley_node_dimensions_F2']
        for r in q['conley_generalized_rank_queries']:
            a,b=r['start_node'],r['end_node']
            r['rank_F2']=[min(dims[d][a],dims[d][b]) for d in range(3)]
        gdims=q['morse_graph_node_dimensions_F2']
        for r in q['morse_graph_generalized_rank_queries']:
            a,b=r['start_node'],r['end_node']
            r['rank_F2']=[min(gdims[d][a],gdims[d][b]) for d in range(2)]
    return x


def identity_loop_transport(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        for rec in q['conley_loop_signatures']:
            occ=2*rec['start_occurrence']
            dims=q['conley_node_dimensions_F2']
            for d,deg in enumerate(rec['degrees']):
                n=dims[d][occ]
                deg.update(relation_dimension_F2=n,left_projection_rank_F2=n,right_projection_rank_F2=n,
                           regular_rank_F2=n,fixed_subspace_dimension_F2=n,power_regular_ranks_F2=[n,n,n,n],
                           kronecker_invariants_F2=_identity_kronecker(n))
        for rec in q['morse_graph_loop_signatures']:
            occ=2*rec['start_occurrence']
            dims=q['morse_graph_node_dimensions_F2']
            for d,deg in enumerate(rec['degrees']):
                n=dims[d][occ]
                deg.update(relation_dimension_F2=n,left_projection_rank_F2=n,right_projection_rank_F2=n,
                           regular_rank_F2=n,fixed_subspace_dimension_F2=n,
                           kronecker_invariants_F2=_identity_kronecker(n))
    return x


def zero_graph_queries(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['morse_graph_node_dimensions_F2']=[[],[]]
        q['morse_graph_generalized_rank_queries']=[]
        q['morse_graph_loop_signatures']=[]
        q['morse_graph_holonomy_word_signatures']=[]
    return x


def zero_holonomy(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['conley_holonomy_word_signatures']=[]
        q['morse_graph_holonomy_word_signatures']=[]
    return x


def _identity_kronecker(n):
    return {'right_minimal_indices':[],'left_minimal_indices':[],'infinite_elementary_divisors':[],
            'finite_elementary_divisors':[[3,1]]*n}


def _all_relation_records(q):
    for rec in q['conley_loop_signatures']+q['morse_graph_loop_signatures']:
        for deg in rec['degrees']:
            yield deg
    for rec in q['conley_holonomy_word_signatures']+q['morse_graph_holonomy_word_signatures']:
        for deg in rec['degrees']:
            for w in deg['word_signatures']:
                yield w


def graph_assumption_kronecker(ans):
    """Treat every continuation relation as the graph of an automorphism on its
    regular part: report only finite blocks x+1, one per regular-rank dimension."""
    x=copy.deepcopy(ans)
    for q in x['queries']:
        for w in _all_relation_records(q):
            w['kronecker_invariants_F2']={'right_minimal_indices':[],'left_minimal_indices':[],
                'infinite_elementary_divisors':[],'finite_elementary_divisors':[[3,1]]*w['regular_rank_F2']}
    return x


def singular_blocks_as_zero(ans):
    """Keep the regular part but collapse every singular block into zero-relation
    blocks (left index 0), i.e. ignore minimal indices."""
    x=copy.deepcopy(ans)
    for q in x['queries']:
        for w in _all_relation_records(q):
            k=w['kronecker_invariants_F2']
            used=sum(e for e in k['right_minimal_indices'])+sum(h+1 for h in k['left_minimal_indices'])
            k['right_minimal_indices']=[]; k['left_minimal_indices']=[0]*used
    return x


def drop_converse_words(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        for rec in q['conley_holonomy_word_signatures']+q['morse_graph_holonomy_word_signatures']:
            for deg in rec['degrees']:
                deg['word_signatures']=[w for w in deg['word_signatures'] if w['word'].isupper()]
    return x


def no_barcodes(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['conley_zigzag_barcodes_F2']=[[],[],[]]
    return x


def barcode_from_node_dimensions(ans):
    """Every node class is its own length-zero bar (ignores all maps)."""
    x=copy.deepcopy(ans)
    for q in x['queries']:
        dims=q['conley_node_dimensions_F2']
        q['conley_zigzag_barcodes_F2']=[[[i,i,v] for i,v in enumerate(dims[d]) if v] for d in range(3)]
    return x


_ZERO_SHORTCUT_CACHE={}

def absorbing_zero_shortcut(ans):
    """Run the retired oracle logic (zero relation treated as absorbing)."""
    case=ans['case_id']
    if case not in _ZERO_SHORTCUT_CACHE:
        inp=ROOT/'tests'/'heldout'/f'{case}.json'
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out.json'
            subprocess.run([sys.executable,str(ROOT/'authoring'/'evidence'/'zero_shortcut_variant.py'),str(inp),str(out)],check=True)
            _ZERO_SHORTCUT_CACHE[case]=json.loads(out.read_text())
    return _ZERO_SHORTCUT_CACHE[case]


def no_graph_barcodes(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['morse_graph_zigzag_barcodes_F2']=[[],[]]
    return x


def graph_barcode_capped_sweep(ans, cap=1024):
    """What a per-start sweep truncated at `cap` nodes (to fit the time limit)
    reports: every interval longer than the cap is cut at start+cap-1."""
    x=copy.deepcopy(ans)
    for q in x['queries']:
        out=[]
        for bars in q['morse_graph_zigzag_barcodes_F2']:
            acc={}
            for s,e,m in bars:
                e2=min(e,s+cap-1); acc[(s,e2)]=acc.get((s,e2),0)+m
            out.append([[s,e,m] for (s,e),m in sorted(acc.items())])
        q['morse_graph_zigzag_barcodes_F2']=out
    return x


def graph_barcode_from_node_dimensions(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        dims=q['morse_graph_node_dimensions_F2']
        q['morse_graph_zigzag_barcodes_F2']=[[[i,i,v] for i,v in enumerate(dims[d]) if v] for d in range(2)]
    return x


def no_index_maps(ans):
    x=copy.deepcopy(ans)
    for q in x['queries']:
        q['conley_index_maps']=[]
    return x


def identity_index_maps(ans):
    """Report every index map as the identity (homology of the pair only)."""
    x=copy.deepcopy(ans)
    for q in x['queries']:
        for r in q['conley_index_maps']:
            for d in r['degrees']:
                n=d['space_dimension_F2']; d['rank_F2']=n; d['elementary_divisors_F2']=[[3,1]]*n
    return x


_EXITSET_CACHE={}

def naive_exit_set(ans):
    """Index pair (P, P minus Morse boxes) instead of the non-reaching exit set."""
    case=ans['case_id']
    if case not in _EXITSET_CACHE:
        inp=ROOT/'tests'/'heldout'/f'{case}.json'
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'out.json'
            subprocess.run([sys.executable,str(ROOT/'authoring'/'evidence'/'index_exitset_variant.py'),str(inp),str(out)],check=True)
            _EXITSET_CACHE[case]=json.loads(out.read_text())
    return _EXITSET_CACHE[case]

ROUTES={
 'full_reference': lambda x: copy.deepcopy(x),
 'greedy_single_branch': greedy_single_branch,
 'endpoint_dimension_as_rank': endpoint_dimension_as_rank,
 'identity_loop_transport': identity_loop_transport,
 'zero_graph_queries': zero_graph_queries,
 'zero_holonomy': zero_holonomy,
 'absorbing_zero_shortcut': absorbing_zero_shortcut,
 'graph_assumption_kronecker': graph_assumption_kronecker,
 'singular_blocks_as_zero': singular_blocks_as_zero,
 'drop_converse_words': drop_converse_words,
 'no_barcodes': no_barcodes,
 'barcode_from_node_dimensions': barcode_from_node_dimensions,
 'no_graph_barcodes': no_graph_barcodes,
 'graph_barcode_capped_sweep': graph_barcode_capped_sweep,
 'graph_barcode_from_node_dimensions': graph_barcode_from_node_dimensions,
 'no_index_maps': no_index_maps,
 'identity_index_maps': identity_index_maps,
 'naive_exit_set': naive_exit_set,
}


def main():
    cases=sorted(ORACLE_DIR.glob('heldout_*.json'))
    report={'heldout_cases':len(cases),'required_cases':4,'routes':[]}
    for name,fn in ROUTES.items():
        case_pass={}
        for p in cases:
            expected=json.loads(p.read_text())
            got=fn(expected)
            ok,_=vc.compare_case(got,expected)
            case_pass[p.name]=int(ok)
        passed=sum(case_pass.values())
        report['routes'].append({'route':name,'passed_cases':passed,'case_pass':case_pass,'passes_verifier':passed>=4})
    (ROOT/'authoring'/'evidence'/'ablation_report.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=='__main__': main()
