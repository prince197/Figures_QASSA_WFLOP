# Conley Morse continuation trajectories

## Difficulty

This task targets a difficult part of parameterized Conley Morse computation: transporting exact topological information through long bifurcation trajectories when continuation is genuinely relation valued. A solver that locally matches recurrent components, keeps only Betti numbers, or replaces a relation by a convenient map can remain plausible for many steps and still be wrong after repeated split and merge events.

The hidden families use phase grids from 34 by 34 to 36 by 36 and parameter lattices of 15 by 15, 14 by 16, 16 by 14 and 15 by 14. Each family has three low-backtracking trajectories ranging from 2800 to 4800 occurrences, so the largest zigzag has 9599 nodes. Across the hidden set, the trajectories visit 98.3 to 100 percent of lattice edges, immediate reversals stay between 0.47 and 1.45 percent, and the propagated Morse selection changes on 94.6 to 97.7 percent of steps. Selections reach 25 Morse sets and the Morse reachability-graph H1 dimension reaches 333. Every hidden query asks for 384 Conley generalized ranks, the complete Conley zigzag barcode in all three degrees, the complete Morse-graph zigzag barcode in degrees 0 and 1 (351 to 773 distinct intervals in each of H0 and H1; H2 vanishes on these families), 128 Morse-graph generalized ranks, 20 Conley loops, 10 graph loops, 10 Conley holonomy probes and 6 graph holonomy probes.

Five layers carry the difficulty.

Relation composition must be exact. Continuation along a window is the composite of the linear relations of all zigzag arrows. The zero relation is not absorbing under relation composition: {(0,0)} composed with an edge relation S is {(0,w) : (0,w) in S}, which is nonzero as soon as a later backward map has a kernel. Hidden loop windows and holonomy words are deliberately chosen so that some compositions pass through the zero relation and become nonzero again; a solver that stops composing at a zero relation reports wrong dimensions, projection ranks and Kronecker types. The public example contains the same phenomenon.

Loop relations must be classified up to conjugacy, not only by ranks. For every Conley and Morse-graph loop relation and for every holonomy word, the grader asks for the Kronecker invariants of the endorelation: right and left minimal indices, infinite elementary divisors, and finite elementary divisors over F2. These form a complete conjugacy invariant of a linear relation. On the Morse-graph layer the relations live on spaces of dimension up to 333 and decompose into mixtures of full-relation blocks, zero-relation blocks, kernel and multivalued blocks, nilpotent blocks and identity blocks, so a classifier that only recognises graphs of maps, or only the zero and identity relations, fails.

Holonomy words include converse relations. Each probe yields two closed endorelations A and B on the same space. The grader evaluates A, B, AB, BA, ABA, BAB, AABB, BBAA, ABAB, BABA and the mixed words Ab, aB, ABab and abAB, where lower-case letters are converse relations. These words depend on how the transported subspaces sit inside the common space and are not determined by individual loop ranks or dimensions.

The Conley zigzag barcode is requested in full. The interval decomposition of each 2L-1 node Conley zigzag must be reported exactly, which requires the generalized rank function on all intervals rather than on a sparse sample.

The Morse-graph H1 barcode is requested in full, and this is the algorithmic core of the task. The Morse-graph H1 zigzag of each hidden query has up to 9599 nodes carrying spaces of dimension up to 333, and its interval decomposition has 6309 to 13797 distinct intervals (total multiplicity up to 114565), including intervals that span the whole trajectory. Because intervals are long, the straightforward route of sweeping a composed relation from every start node needs between 15.7 and 46.1 million compositions of relations on spaces of dimension up to 333 per query, which cannot finish within the 300 second limit. A solver must use a genuine zigzag persistence algorithm, for example a single sweep that maintains, for all start nodes at once, the range and multivalued-part flags of the composed relations (as the reference does), or a right-filtration or right-to-left domain and kernel sweep (as the independent implementation does). Correctness can be checked against brute force on small zigzags and against the public example; efficiency is what the hidden trajectories test.

The data is synthetic but realistic: each frame is an exact closed outer approximation of a coupled nonlinear map, and the generated families exhibit recurrent SCCs, split and merge continuation, nontrivial relative homology, and nontrivial continuation holonomy. Researchers in rigorous numerics, computational dynamics, computational topology, and topological data analysis perform closely related work.

## Reference solution

The reference solver reconstructs every box map with exact rational interval bounds. It computes recurrent strongly connected components, their reachability order, and the relative cubical homology groups associated with each propagated Morse selection. For an adjacent parameter pair it builds the union box map and derives the continuation relation from union graph strongly connected components. Selected Morse ids are propagated through that relation in the direction of travel.

For the Conley trajectory, each occurrence becomes a relative pair and each adjacent pair contributes a bridge pair. Inclusion induced maps are computed over F2 with bitset Gaussian elimination that retains homology coordinates. Repeated relative complexes are cached. Sparse generalized ranks compose endpoint relations incrementally. The full barcode is computed by a single left-to-right sweep that maintains, for every start node b, the range flag ran(R_{b,k}) and the multivalued-part flag mul(R_{b,k}) inside the current node space; the generalized rank is their dimension difference, and interval multiplicities follow by inclusion-exclusion.

Closed loop windows keep the complete relation between the initial and final homology spaces, composed arrow by arrow without any early exit. From that subspace the solver computes both coordinate projection ranks, the regular rank, the diagonal fixed subspace, the regular ranks of the first four powers for Conley loops, and the Kronecker invariants. Kronecker invariants are read from four functorial chains of subspaces, R^j(0), R^-j(0), R^k(U) and R^-k(U), whose dimensions are additive over the Kronecker decomposition; the finite divisors other than x come from the automorphism induced on the regular core, via characteristic-polynomial factorisation over F2 and kernel dimensions of p(T)^j. Holonomy probes retain two closed endorelations and their converses and compose the prescribed words.

The Morse graph layer is deliberately different. At each frame it forms the undirected comparability graph of the recurrent Morse order and computes H0 and H1. Endpoint vertices map through the union frame SCC quotient; graph edges either map to the corresponding union edge or collapse to zero. The resulting graph homology maps are used for their own generalized rank, loop relation and holonomy queries. Other mathematically correct strategies are allowed; the verifier does not require this implementation.

## Verification

The only graded artifact is `/app/trajectory_solver.py`. The verifier compiles its source in memory, copies the file into a fresh temporary sandbox and runs that copy as the unprivileged `nobody` user. Only the current hidden input is copied into that sandbox. The held-out directory and sealed oracle directory are owned by root and are not readable by the submitted process. Network access is disabled in the verifier image. The verifier runtime matches the agent runtime: Python 3.11 with numpy 2.1.3, scipy 1.14.1 and networkx 3.4.2. Each hidden invocation has the 300 second wall-clock limit stated in the instruction.

The parser accepts harmless JSON presentation differences, including whitespace, key order and integral spellings such as `100.0` for `100`. Duplicate keys, non-finite numbers, non-integral numbers, symbolic links and non-regular output files are rejected. Redundant echoes such as parameter paths are not used as scientific gates. The scientific comparisons are made on transported selections, Conley dimensions, generalized-rank records, Conley and Morse-graph barcodes, loop records, holonomy-word records (including Kronecker invariants), and the corresponding Morse-graph quantities.

The sealed oracles were regenerated after removing a defect in the earlier oracle generator, which stopped relation composition as soon as the running relation became zero. That shortcut is invalid for linear relations, and it is now an explicit ablation (`absorbing_zero_shortcut`) that fails all four hidden families. The oracle generator composes every arrow of every window.

Gate: transported selections and node dimensions. Threshold: exact on every queried trajectory. Reference worst case: exact on all four held-out families. Independent implementation worst case: exact on all four. Nearest wrong route: greedy one-to-one component matching, which diverges after split-merge continuation and fails the held-out gate.

Gate: generalized-rank, barcode, loop, and holonomy records. Threshold: exact equality for every record in every query. Reference worst case: 100 percent. Independent implementation worst case: 100 percent. Endpoint-dimension rank inference, identity-loop transport, the absorbing-zero shortcut, graph-only Kronecker classification, ignoring minimal indices, omitting converse words, and barcodes read from node dimensions each fail every held-out family.

Gate: held-out coverage. Threshold: all four sealed families must pass every scientific gate. Reference result: four of four. Independent result: four of four. The four held-out families are the complete graded calibration set; their grids, parameter values, path lengths, loop locations, holonomy probes, and query windows differ from the public case.

Gate: runtime. Threshold: 300 seconds for each hidden invocation. Reference worst case in local calibration is reported in `authoring/evidence/validation_report.json`; the independent implementation is reported in `authoring/evidence/independent_report.json`. Both are well below the limit.

The independent implementation (`authoring/evidence/independent_solver.py`) shares no relation, Kronecker or barcode code with the reference. It builds box maps with integer arithmetic, uses a Kosaraju SCC route, obtains every window relation as the projection of the section space (inverse limit) of the restricted zigzag, computed chunk-wise and glued by pullbacks, computes Kronecker invariants from kernel dimensions of the pencil beta - s alpha over polynomial quotient rings F2[s]/(q) and from polynomial-kernel ranks, and computes barcodes from per-start generalized-rank sweeps. Both Kronecker routes were additionally fuzz-tested on thousands of randomly conjugated direct sums of known blocks.

## Ablation ladder

Full reference: exact outer approximation, relation valued continuation, chain level induced maps, exact long range relation composition, Conley barcodes, loop subspace invariants, Kronecker invariants, converse holonomy words, and Morse graph homology. Held out result: 4 of 4 pass.

Replace continuation by a greedy one to one Morse matching: propagated selections diverge after split merge events; loop and generalized rank records disagree. Held out result: fail.

Keep only Betti dimensions and infer interval ranks from endpoint dimensions: long generalized ranks and closed loop regular ranks disagree. Held out result: fail.

Treat every loop as identity transport: loop, power and Kronecker records disagree. Held out result: fail.

Stop composing when the running relation is zero (the retired oracle defect): loop and holonomy records disagree wherever a composition passes through zero. Held out result: fail.

Classify relations as if they were graphs of automorphisms, or collapse all singular Kronecker blocks to zero-relation blocks: Kronecker records disagree. Held out result: fail.

Omit the converse words, or omit the barcode, or read the barcode from node dimensions: fail.

Omit the Morse-graph barcode, read it from node dimensions, or compute it with a per-start sweep capped at 1024 nodes so that it fits the time limit: fail on every hidden family.

Omit the Morse reachability graph layer: node-dimension and graph-transport gates fail even if all Conley outputs are correct. Held out result: 0 of 4 pass.

Return the visible example or hard code public trajectory values: sealed case identifiers, grids, parameters, seeds, paths, and windows differ. Held out result: 0 of 4 pass.
