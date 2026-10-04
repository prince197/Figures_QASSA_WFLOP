# C4 precision — local results and text templates (from the precision agent's final local report)

Local reruns (class iii): 1,649 records: 1,160 bit-identical, 412 float-noise, 11 near, 66 diverged; 1,572 labels
confirmed (1,550 feasible, 22 infeasible), 0 changed; 77 undecided, all SLSQP-based (SciPy SLSQP not bit-reproducible
across CPU/BLAS builds). One diverged MS-SLSQP rerun (mpce_iea36 row 21) ended infeasible — says nothing about the stored
layout; counted as undecided. Cloud: 4,065 records in 35 shards (1,362 SLSQP-based). Local + cloud = 5,714.

After the cloud shards are in analysis/: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 rev3_precision_rerun.py collect`
then `python3 rev3_precision_rerun.py table`; check errors 0, class-(iii) sum 5,714, Open = near + diverged.

TBD text template (B_c, F_c, C_c, X_c, U_c = cloud counts):
"All 5,714 class-(iii) records were rerun through the unchanged drivers with a 17-digit coordinate writer (1,649 locally,
4,065 on cloud workers): 1,160+B_c reproduced the stored record bit for bit and 412+F_c to float noise; for all of them
the strict label from the full-precision layout confirmed the stored label (1,572+C_c confirmed, 0+X_c changed). The
remaining 77+U_c records are MS-SLSQP-based runs (SciPy SLSQP is not bit-reproducible across CPU/BLAS builds), whose
stored labels remain undecidable at 1e-6 m."

Why 56,540 audited > 52,930 runs in Table 3: 56,540 = 36,190 + 11,640 + 5,500 + 3,210 (5,500 further rows in used files:
2,040 modified-VNS and 2,040 earlier MS-SLSQP runs in fresh_grid, 720 LX-SSA-VNS split runs, 700 Horns Rev 1 rows with the
earlier direction binning; 3,210 runs of superseded studies: fresh_hgrid 2,040, fresh_hsplit 720, fresh_hr 390,
mpce_hr16new 60). The 5,100 revision-3 runs store 17-digit coordinates (no audit needed).
