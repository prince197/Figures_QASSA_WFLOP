# conley-morse: changes after the failed step-9 judge review

The judge's diagnosis was correct. The sealed oracles were wrong, not the agents.
The reference solver, and the "independent" solver that shared its code, stopped
composing linear relations as soon as the running relation became zero
(`if not R: break` in `_relation_between_nodes` and `holonomy_word_signature`).
The zero relation is not absorbing, because {(0,0)} composed with S equals
{(0,w) : (0,w) in S}. Agents that composed exactly were marked wrong on 2 probes
in heldout_0, 1 probe in heldout_1, and 6 loops plus 9 probes in heldout_3. That
single defect explains the failures on *verifier correctness* (8/8),
*solution discoverability* (6/8) and *meaningful difficulty* (5/8).

## 1. Correctness fixes

* Removed both invalid early exits from the reference. The one remaining
  `if not R: break` is in the sparse-rank loop, where it is valid (only the
  regular rank is read, and it stays 0). It is now commented as such.
* Regenerated every sealed oracle and the public expected output.
* The retired behaviour is now an explicit ablation
  (`authoring/evidence/zero_shortcut_variant.py`, route `absorbing_zero_shortcut`).
  It fails **all four** hidden families and the public example.
* heldout_2's old parameter family had no windows that pass through the zero
  relation, so the defect could not be detected there. I replaced it with a
  34x34 / 16x14 family that does have such windows (r0=2.885, s0=2.920, c=0.027).
* Rewrote the independent solver so that it shares **no** relation, Kronecker or
  barcode code with the reference:
  * integer box maps and Kosaraju SCCs;
  * window relations obtained as projections of section spaces (inverse limits),
    glued chunk by chunk with pullbacks;
  * Kronecker invariants from kernel dimensions of the pencil beta - s*alpha over
    F2[s]/(q) and from polynomial-kernel ranks;
  * barcodes from per-start rank sweeps.

  It reproduces the public expected output and all four sealed oracles byte for byte.

## 2. Added difficulty (same task and format, with extra graded fields)

1. **Kronecker invariants** (`kronecker_invariants_F2`) for every Conley and
   Morse-graph loop relation and every holonomy word: right and left minimal
   indices, infinite elementary divisors, and finite elementary divisors over F2.
   Together these are a complete conjugacy invariant of a linear relation. On the
   Morse-graph layer, relations of dimension up to 333 split into mixtures of
   full, zero, multivalued, kernel (nilpotent) and identity blocks.
2. **Converse holonomy words**: `Ab`, `aB`, `ABab`, `abAB` were added to the
   original ten words (lower-case letters are converse relations).
3. **Full Conley zigzag barcodes** (`conley_zigzag_barcodes_F2`) for H0, H1 and H2.
4. **Structure-aware case generation**: loops and probes are now chosen to cover
   rare relation types and windows that pass through the zero relation. Every
   hidden family and the public example contain right, left, infinite, nilpotent
   and identity blocks, plus zero-crossing windows.
5. There are now 10 Conley holonomy probes per hidden query (previously 8).

The spec for all of this is written out precisely in `environment/data/README.txt`,
with block definitions, code conventions, reference points and the barcode
formula. The public example exercises every block kind and the zero-crossing trap.

## 3. Evidence (regenerated)

* Reference in the verifier sandbox: 4/4, 78-99 s per case (limit 300 s).
* Independent solver: exact match on public + 4/4 hidden.
* Ablations: every wrong route (11 of them) passes 0/4.
* Anti-cheat: nop, case-id-only, oracle snoop and symlink output all get reward 0.
* Kronecker routes: both fuzz-tested on 1500+ randomly conjugated block sums.
  Barcode sweep: fuzz-tested against brute-force limit/colimit ranks.

The bundle is in `conley-morse/`. The same bundle, zipped for upload, is
`conley-morse.zip`.

## 4. Round 2: the easiness probe passed 3/3, so the task was hardened again

* **New graded field `morse_graph_zigzag_barcodes_F2`**: the full interval
  decomposition of the Morse-graph zigzag in H0 and H1.
* **Scale:** the H1 zigzags have up to 9599 nodes on spaces of dimension up to 333.
  Each query has 6309 to 13797 distinct intervals (total multiplicity up to
  114565), and some span the whole trajectory.
* **Why naive fails:** sweeping a composed relation from every start node needs
  15.7 to 46.1 million relation compositions per query, which cannot fit in the
  300 s limit. A real zigzag-persistence algorithm is required. The spec
  already defines the barcode exactly, and it says the barcode covers every
  interval, not just the listed windows.
* **Reference:** a single left-to-right flag sweep takes about 12 s for the
  largest query. The whole case runs in 93-132 s in the verifier sandbox.
* **Independent check:** a right-to-left sweep over domain and kernel flags
  (different subspaces, opposite direction) matches all four oracles and the
  public example exactly. On 2000 random zigzags the reference sweep, the
  reverse sweep, a per-start sweep and brute-force limit/colimit ranks agree.
* **New ablations:** no graph barcode, graph barcode read from node dimensions,
  and a per-start sweep capped at 1024 nodes. All score 0/4. 14 ablations now
  score 0/4 in total.
* **Unchanged:** all previously graded records. Only the new field was added.
