#!/usr/bin/env python3
"""Focused checks of new serialization, reliability arithmetic and PSO correction."""
import math, random, struct, hashlib, json
from pathlib import Path
import numpy as np
from audit_archive import wilson, sign_p, holm
from record_io import encode_coordinates, decode_coordinates

def main():
    rng=random.Random(20261004)
    coordinates=[(0.0,-0.0),(1e-320,-1e-320),(1e308,-1e308)]
    coordinates += [(rng.uniform(-1e6,1e6),rng.uniform(-1e6,1e6)) for _ in range(1000)]
    restored=decode_coordinates(encode_coordinates(coordinates))
    assert len(restored)==len(coordinates)
    assert all(struct.pack('!d',a)==struct.pack('!d',b) for row1,row2 in zip(coordinates,restored) for a,b in zip(row1,row2))
    for bad in [[(math.nan,1)],[(math.inf,1)],[(1,)],[]]:
        try: encode_coordinates(bad)
        except ValueError: pass
        else: raise AssertionError('Invalid coordinates accepted')
    assert abs(sign_p(17,0)-2**-16)<1e-15
    assert sign_p(0,0)==1 and sign_p(10,20)==sign_p(20,10)
    assert abs(sign_p(28,2)-2*(1+30+435)/2**30)<1e-15
    assert np.allclose(holm([0.01,0.04,0.03]),[0.03,0.06,0.06])
    lo,hi=wilson(0,30);assert abs(lo)<1e-14 and 0.11<hi<0.12
    lo2,hi2=wilson(30,30);assert abs(hi2-1)<1e-14 and abs(lo2-(1-hi))<1e-14
    settings=[]
    for w,c1,c2 in [(0.7,2,2),(0.7298,1.49618,1.49618),(0,0,1)]:
        mu=(c1+c2)/2;sigma2=(c1*c1+c2*c2)/12;m=1+w-mu
        M=np.array([[m*m+sigma2,-2*w*m,w*w],[m,-w,0],[1,0,0]])
        radius=float(max(abs(np.linalg.eigvals(M))))
        settings.append(dict(w=w,c1=c1,c2=c2,spectral_radius=radius))
    assert settings[0]['spectral_radius']>1 and settings[1]['spectral_radius']<1
    assert abs(settings[2]['spectral_radius']-1/3)<1e-14
    # With c1=0, c2=1, p!=g and w=0, q=0 despite p!=g.
    assert 0**2*1**2*(4-7)**2/(6*(0+1)**2)==0
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'validation/raw_data_checksums.json').read_text())
    for entry in manifest:
        assert hashlib.sha256((root/entry['path']).read_bytes()).hexdigest()==entry['sha256']
    result=dict(status='PASS',float_roundtrips=len(coordinates)*2,raw_csv_files_verified=len(manifest),
                checks=['nonfinite/malformed coordinates rejected','exact sign-test reference values','Holm reference family',
                        'Wilson boundary symmetry','PSO stability settings and degenerate-coefficient counterexample'],
                pso_settings=settings,limitations=['No optimizer reruns: required original modules absent.'])
    (root/'validation/focused_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"PASS: {len(coordinates)*2} binary float round-trips, statistical reference values, PSO correction and {len(manifest)} raw CSV hashes.")

if __name__=='__main__': main()
