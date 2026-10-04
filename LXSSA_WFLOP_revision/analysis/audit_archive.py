"""Independent archive audit and label-based reliability summaries.

Does not import absent optimizer/analysis modules and does not rerun optimizers.
Run from any directory: python experiments/audit_archive.py
Use --out DIR to choose the derived-output directory.
Existing raw CSVs are read-only. Legacy feasibility labels are not relabelled.
"""
from pathlib import Path
import argparse, ast, csv, hashlib, json, math, platform, sys
import numpy as np
import pandas as pd
from record_io import decode_coordinates
import rev2_site_model as lg

ROOT = Path(__file__).resolve().parent.parent
FAMILIES = dict(ga=2040, gahr=60, gaiea=120, grad=240, laplace=2160,
                spacing=6480, lg16=270, lg16b=270)
LABELS = dict(PSOBV="PSO-VNS", PSOC="PSO", SSABV="SSA-VNS", SSA="SSA",
              LXSSA="LX-SSA", DE="DE", BVNS="VNS", SLSQP="MS-SLSQP", RSDVNS="RSD-VNS")


def wilson(successes, n, z=1.959963984540054):
    p=successes/n; den=1+z*z/n
    mid=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return mid-half,mid+half


def sign_p(wins, losses):
    n=wins+losses
    return min(1.0,2*sum(math.comb(n,k) for k in range(min(wins,losses)+1))/2**n) if n else 1.0


def holm(values):
    order=np.argsort(values); out=np.ones(len(values)); previous=0.0
    for rank,i in enumerate(order):
        previous=max(previous,min(1.0,(len(values)-rank)*values[i]));out[i]=previous
    return out


def require(condition, message):
    if not condition: raise ValueError(message)


def load_family(family, expected):
    files=sorted((ROOT/'experiments').glob(f'rev2_{family}_s*of*.csv'))
    require(bool(files),f'{family}: no shards')
    shard_ids=[]; shard_totals=[]
    import re
    for p in files:
        m=re.search(r'_s(\d+)of(\d+)\.csv$',p.name)
        shard_ids.append(int(m[1]));shard_totals.append(int(m[2]))
    require(len(set(shard_totals))==1 and sorted(shard_ids)==list(range(shard_totals[0])),f'{family}: incomplete shards')
    df=pd.concat([pd.read_csv(p) for p in files],ignore_index=True)
    require(len(df)==expected,f'{family}: expected {expected} records, found {len(df)}')
    require(df.Feasible.isin([True,False]).all(),f'{family}: invalid feasibility labels')
    key=['Algorithm','Dataset','Radius','Turbines','Seed','Budget','Init']+(['Spacing'] if family=='spacing' else [])
    require(not df.duplicated(key).any(),f'{family}: duplicate run keys')
    require((df.Calls==df.Budget).all(),f'{family}: budget mismatch')
    require(np.isfinite(df[['Objective','Ideal','WakeLoss','MinSpacing','Seconds']]).all().all(),f'{family}: nonfinite scalar')
    error=float(abs(df.Ideal-df.Objective-df.WakeLoss).max())
    require(error<1e-7,f'{family}: objective identity mismatch')
    if family=='grad': require((df.FunCalls+df.CG*df.GradCalls+df.Unused==df.Budget).all(),'gradient charge mismatch')
    info=dict(records=len(df),shards=len(files),loss_identity_max_abs_error=error,
              feasible_labels=int(df.Feasible.sum()),rounded_spacing_failures=0,
              rounded_boundary_failures=0,rounded_any_failures=0,boundary_check='not available for Horns Rev 1' if family=='gahr' else 'checked')
    replay=[]
    for row in df.itertuples():
        xy=np.array(decode_coordinates(row.Coordinates));require(len(xy)==row.Turbines,f'{family}: coordinate count mismatch')
        dist=np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2);np.fill_diagonal(dist,np.inf)
        smin=row.SMin if family=='spacing' else (260 if family in ('gaiea','grad') else (279 if family.startswith('lg16') else (320 if family=='gahr' else 308)))
        spacing_bad=dist.min()<smin-1e-6
        if family.startswith('lg16'): boundary_bad=lg.outside_distance(xy,lg.site()[1]).max()>1e-6
        elif family=='gahr': boundary_bad=False
        else: boundary_bad=np.linalg.norm(xy,axis=1).max()>row.Radius+1e-6
        if row.Feasible:
            info['rounded_spacing_failures']+=int(spacing_bad)
            info['rounded_boundary_failures']+=int(boundary_bad)
            info['rounded_any_failures']+=int(spacing_bad or boundary_bad)
            if family.startswith('lg16'): replay.append((abs(lg.aep_gwh(xy)-row.Objective),row.Algorithm,int(row.Seed)))
    if replay:
        largest=max(replay);info['rounded_aep_max_abs_GWh']=largest[0]
        info['largest_replay_difference_run']=dict(algorithm=largest[1],seed=largest[2])
        info['rounded_aep_errors_above_0_001_GWh']=int(sum(v[0]>.001 for v in replay))
    return df,info


def reliability(frames):
    rows=[];pairs=[]
    for family in ('lg16','lg16b'):
        df=frames[family];budget=int(df.Budget.iloc[0])
        for alg,x in df.groupby('Algorithm'):
            k=int(x.Feasible.sum());lo,hi=wilson(k,len(x))
            rows.append(dict(algorithm=alg,label=LABELS[alg],budget=budget,runs=len(x),feasible=k,
                             success=k/len(x),wilson95_lower=lo,wilson95_upper=hi,
                             conditional_mean_GWh=float(x.loc[x.Feasible,'Objective'].mean()) if k else None))
        focus=df[df.Algorithm=='PSOBV'].set_index('Seed')
        for alg in sorted(set(df.Algorithm)-{'PSOBV'}):
            other=df[df.Algorithm==alg].set_index('Seed')
            require(set(focus.index)==set(other.index),'paired seeds do not match')
            scores=[]
            for seed in sorted(focus.index):
                a,b=focus.loc[seed],other.loc[seed]
                if a.Feasible!=b.Feasible: score=float(a.Feasible)
                elif not a.Feasible: score=0.5
                else: score=0.5 if abs(a.Objective-b.Objective)<=1e-9 else float(a.Objective>b.Objective)
                scores.append(score)
            scores=np.array(scores);wins=int((scores==1).sum());losses=int((scores==0).sum());ties=int((scores==.5).sum())
            rng=np.random.default_rng(20261004+budget)
            boot=scores[rng.integers(0,len(scores),size=(20000,len(scores)))].mean(axis=1)
            lo,hi=np.quantile(boot,[.025,.975])
            pairs.append(dict(budget=budget,first='PSOBV',second=alg,label=LABELS[alg],wins=wins,ties=ties,losses=losses,
                              superiority=float(scores.mean()),bootstrap95_lower=float(lo),bootstrap95_upper=float(hi),
                              exact_sign_p=sign_p(wins,losses)))
    adjusted=holm([p['exact_sign_p'] for p in pairs])
    for p,v in zip(pairs,adjusted):p['holm16_p']=float(v)
    return rows,pairs


def latex_reliability(rows,pairs):
    by={(r['algorithm'],r['budget']):r for r in rows}
    text=[r'\begin{table}[!htbp]',r'\centering',r'\caption{Lillgrund reliability and conditional AEP from the supplied run records.}',r'\label{tab:archive-reliability}',r'\small\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{lcccccc}',r'\toprule',r'& \multicolumn{3}{c}{6,030 evaluations} & \multicolumn{3}{c}{30,030 evaluations} \\',r'Method & Feas. & 95\% CI & Mean AEP & Feas. & 95\% CI & Mean AEP \\',r'\midrule']
    for alg in ['PSOBV','SLSQP','PSOC','BVNS','SSABV','RSDVNS','SSA','LXSSA','DE']:
        cells=[LABELS[alg]]
        for b in (6030,30030):
            r=by[alg,b];cells.extend([f"{r['feasible']}/{r['runs']}",f"[{r['wilson95_lower']:.2f}, {r['wilson95_upper']:.2f}]",f"{r['conditional_mean_GWh']:.3f}" if r['conditional_mean_GWh'] is not None else '--'])
        text.append(' & '.join(cells)+r' \\')
    text.extend([r'\bottomrule',r'\end{tabular}',r'\par\smallskip\parbox{\textwidth}{\footnotesize Feas.: archived feasible labels; CI: marginal Wilson interval for success probability; mean AEP: conditional on those labels, GWh/yr. These intervals concern independent seeds at this block, not variation over sites. Rounded archival coordinates do not independently reproduce all labels; Section~\ref{sec:archive-precision}.}',r'\end{table}'])
    text.extend([r'\begin{table}[!htbp]',r'\centering',r'\caption{All-run paired outcomes for PSO-VNS on Lillgrund, using the archived feasibility labels.}',r'\label{tab:archive-paired}',r'\small\setlength{\tabcolsep}{4pt}',r'\begin{tabular}{rlcccc}',r'\toprule',r'Budget & Compared with & W/T/L & Superiority & 95\% CI & Holm $p$ \\',r'\midrule'])
    for p in pairs:
        text.append(f"{p['budget']:,} & {p['label']} & {p['wins']}/{p['ties']}/{p['losses']} & {p['superiority']:.3f} & [{p['bootstrap95_lower']:.3f}, {p['bootstrap95_upper']:.3f}] & {p['holm16_p']:.3g}"+r' \\')
    text.extend([r'\bottomrule',r'\end{tabular}',r'\par\smallskip\parbox{\textwidth}{\footnotesize A feasible run beats an infeasible run; if both are feasible, larger AEP wins (ties within $10^{-9}$ GWh); two infeasible runs tie. Superiority is $(W+T/2)/30$. CI: 20,000 paired-seed percentile bootstrap samples, seed $20261004+B$. Two-sided exact sign tests exclude ties; Holm correction covers all 16 comparisons in this table. This additional endpoint does not replace the original signed-rank tests and gives no ordering between failed layouts.}',r'\end{table}'])
    return '\n'.join(text)+'\n'


def local_missing_modules():
    missing=set()
    core={'mpce_results','mpce_inference_extra','mpce_experiments','iea37_experiments','iea37_model','init_hook','authors_optimizers','authors_objective','wflop_model','feasible_init','extra_baselines','original_vns','hornsrev_model','hybrid_lxssa_bvns','rs_vns'}
    for name in core:
        if not any((base/(name+'.py')).exists() for base in (ROOT/'experiments',ROOT/'latex_source'/'analysis')): missing.add(name)
    return sorted(missing)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=ROOT/'validation')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    frames={};audit={}
    for family,n in FAMILIES.items():frames[family],audit[family]=load_family(family,n)
    rows,pairs=reliability(frames)
    for name,data in [('reliability.csv',rows),('paired_outcomes.csv',pairs)]:
        with (args.out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=data[0].keys());writer.writeheader();writer.writerows(data)
    summary=dict(record_audit=audit,total_records=sum(len(d) for d in frames.values()),reported_total_optimization_runs=47830,
                 missing_original_optimization_records=36190,missing_modules=local_missing_modules(),
                 python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,
                 limitations=['Stored-record verification only; optimizers not rerun.','Original analysis modules and primary records absent.',
                              'Legacy rounded coordinates cannot recover original full precision.','Horns Rev boundary geometry cannot be checked from supplied model files.'])
    (args.out/'archive_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    (args.out/'reliability_tables.tex').write_text(latex_reliability(rows,pairs))
    manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((ROOT/'experiments').glob('*.csv'))]
    (args.out/'raw_data_checksums.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f"PASS: {summary['total_records']:,} stored records; budgets, unique keys, complete shards and objective identities.")
    print(f"OPEN: {len(summary['missing_modules'])} project modules absent; 36,190 original optimization records absent.")
    print('OPEN: rounded-coordinate replay limitations are reported, not repaired by relabelling.')
    print('Derived outputs:',args.out)


if __name__=='__main__':main()
