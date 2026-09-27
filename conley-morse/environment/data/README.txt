COUPLED LOGISTIC CONTINUATION TRAJECTORY DATA

The verifier invokes the submitted program on unseen files with the same schema as example.json. All arithmetic-defining values in the JSON are integers. Work over F2 for homology.

MODEL AND INDEXING

phase_grid gives nx and ny. Phase box (i,j) has id i + nx*j and is the closed rectangle [i/nx,(i+1)/nx] x [j/ny,(j+1)/ny]. parameter_grid gives na, nb, arrays r_num and s_num, and model.parameter_scale. Parameter (i,j) has id i + na*j. Use r=r_num[i]/parameter_scale, s=s_num[j]/parameter_scale, and c=coupling_num/coupling_den.

The map is the two-component coupled logistic map written in model.formula. Its box map is the standard closed outer approximation: bound each scalar logistic factor exactly on its source interval, including the critical point 1/2 whenever present, combine the component bounds exactly, and include every closed phase box whose rectangle intersects the closed image rectangle. Boundary contact counts.

parameter_edges lists neighboring parameter ids. An edge may be traversed in either direction by a query.

MORSE AND CONTINUATION CONVENTIONS

At one parameter, recurrent Morse sets are the nontrivial strongly connected components of the box map together with singleton components having a self loop. Order them by minimum phase-box id and number them from zero.

For a recurrent Morse set M, P(M) is the set of every phase box reachable from M and E(M)=P(M) minus M. The Conley group in degree d is H_d(K(P(M)),K(E(M));F2), where K(S) is the closed cubical complex carried by box set S.

For a parameter edge p,q, form the union box map by taking the union of the two endpoint target sets box by box. A source Morse id A and target Morse id B are continuation-related exactly when A and B lie in one strongly connected component of this union map. Reverse the ordered pairs when a query traverses the edge backward.

A query starts with the recurrent Morse set containing seed_box_id. Propagate the selected set of Morse ids through each oriented continuation relation. Input cases guarantee a nonempty propagated selection.

CONLEY TRAJECTORY ZIGZAG

At occurrence i, union P(M) over all selected Morse ids, union their Morse boxes B, and use E=P minus B. Between consecutive occurrences use the union of the two endpoint P sets and the union of the two endpoint E sets. Thus a query of L occurrences has 2*L-1 zigzag nodes, alternating endpoint and bridge relative pairs.

The maps are the inclusion-induced homology maps from each endpoint pair into its neighboring bridge pair. A generalized rank on node interval [a,b] is the rank of the canonical map from the inverse limit of that restricted zigzag to its direct limit.

MORSE REACHABILITY GRAPH ZIGZAG

For each parameter, make an undirected graph whose vertices are recurrent Morse ids. Put an edge {u,v} exactly when u and v are strictly comparable in the recurrent reachability order. Build the corresponding graph for every union frame. The endpoint-to-union SCC quotient induces the graph chain map; an endpoint graph edge maps to the union edge between the quotient vertices, or to zero if both vertices collapse. Use H0 and H1 over F2. The node order is again endpoint, bridge, endpoint.

LOOP RELATIONS

loop_windows and graph_loop_windows use occurrence indices, not zigzag-node indices. Each listed window [a,b] is closed in parameter space. Conley loop windows also return to the same propagated Morse selection.

Each zigzag arrow is a linear relation between its two nodes: a forward map f from node k to node k+1 gives {(x, f(x))}, and a backward map g from node k+1 to node k gives {(g(y), y)}. Relations compose exactly as relations: (x,z) lies in XY if and only if there is some y with (x,y) in X and (y,z) in Y. Compose the linear relations along zigzag nodes 2*a through 2*b. For the resulting relation R subset U direct-sum U report dim R, the ranks of its two coordinate projections, regular_rank = left_projection_rank + right_projection_rank - dim R, fixed_subspace_dimension = dim(R intersect {(x,x)}), and the Kronecker invariants of R defined below. Conley loop records also report the regular ranks of R, R^2, R^4, and R^8 in that order.

CONTINUATION HOLONOMY WORDS

holonomy_probes and graph_holonomy_probes contain triples [a,b,c] of occurrence indices. For a Conley holonomy probe, occurrences a, b, and c have the same parameter id and the same propagated Morse selection. For a graph holonomy probe they have the same parameter id. Let A be the closed-loop linear relation from a to b and B the closed-loop relation from b to c on the common endpoint homology space.

Lower-case letters denote converse relations: a = A^T and b = B^T, where X^T = {(y,x) : (x,y) in X}. For each degree, evaluate the fourteen relation words

A, B, AB, BA, ABA, BAB, AABB, BBAA, ABAB, BABA, Ab, aB, ABab, abAB

by composition from left to right, so AB is A followed by B and Ab = A B^T. For every word report relation_dimension_F2, left_projection_rank_F2, right_projection_rank_F2, regular_rank_F2, fixed_subspace_dimension_F2, and kronecker_invariants_F2. Preserve the listed word order exactly.

KRONECKER INVARIANTS OF AN ENDORELATION

Let R be a linear subspace of U direct-sum U over F2 with n = dim U. Two relations R and R' on U are conjugate when R' = {(gx, gy) : (x,y) in R} for some invertible linear g on U. Viewing R as a representation of the Kronecker quiver through its two coordinate projections R -> U, the Kronecker-Weierstrass theorem gives a direct-sum decomposition U = U_1 + ... + U_m such that R is the direct sum of the pieces R_i = R intersect (U_i direct-sum U_i), each R_i is conjugate (inside U_i) to one of the following indecomposable relations, and the multiset of piece types is uniquely determined by R. In each description u_i is a basis of the piece.

right block of index e >= 1, on a piece with basis u_1..u_e: R is spanned by (u_1,0), (u_2,u_1), (u_3,u_2), ..., (u_e,u_(e-1)), and (0,u_e). Its dimension is e+1. The index-1 block is the full relation on a line.

left block of index h >= 0, on a piece with basis u_0..u_h: R is spanned by (u_0,u_1), (u_1,u_2), ..., (u_(h-1),u_h). Its dimension is h. The index-0 block is the zero relation {(0,0)} on a line.

finite block (p, t), with p a monic irreducible polynomial over F2 and t >= 1, on a piece of dimension t*deg(p): R is the graph {(x, Cx)} of a linear map C with a single elementary divisor p^t, for example the companion matrix of p^t.

infinite block of size t >= 1, on a piece of dimension t: R is the converse {(Nx, x)} of the graph of a nilpotent map N with N^t = 0 and N^(t-1) != 0 (a single nilpotent Jordan block).

Report kronecker_invariants_F2 as an object with four lists, each sorted ascending and listing repeated pieces repeatedly:

right_minimal_indices: the indices e of all right blocks.
left_minimal_indices: the indices h of all left blocks.
infinite_elementary_divisors: the sizes t of all infinite blocks.
finite_elementary_divisors: pairs [code(p), t] of all finite blocks, sorted lexicographically, where code(p) is the integer whose binary digit i is the coefficient of x^i in p. Thus x has code 2, x+1 has code 3, x^2+x+1 has code 7, and x^3+x+1 has code 11.

The four lists always satisfy n = sum(e) + sum(h+1) + sum(t*deg p) + sum(t) and dim R = sum(e+1) + sum(h) + sum(t*deg p) + sum(t). Reference points: the identity relation on F2^n has n finite blocks [3,1]; the graph {(x,0)} of the zero map has n finite blocks [2,1]; its converse {(0,y)} has n infinite blocks of size 1; the zero relation {(0,0)} has n left blocks of index 0; the full relation U direct-sum U has n right blocks of index 1; the zero space U = 0 gives four empty lists.

CONLEY ZIGZAG BARCODES

For every query and every degree d in 0, 1, 2, the Conley trajectory zigzag of degree d (all 2*L-1 nodes) decomposes as a direct sum of interval representations, uniquely up to order. An interval [s,e] with 0 <= s <= e <= 2*L-2 is F2 on nodes s..e, zero elsewhere, with identity maps between consecutive nodes inside the interval. Report every interval as [start_node, end_node, multiplicity], sorted by start_node and then end_node, with multiplicity >= 1 and each (start_node, end_node) listed once. Intervals of length zero are included. The multiplicity of [s,e] equals rho(s,e) - rho(s-1,e) - rho(s,e+1) + rho(s-1,e+1), where rho is the generalized rank defined above and rho is zero outside 0 <= s <= e <= 2*L-2.

INPUT QUERY FIELDS

queries[].rank_windows contains zigzag-node intervals for Conley generalized-rank queries.
queries[].graph_rank_windows contains zigzag-node intervals for Morse-graph generalized-rank queries.
queries[].loop_windows contains Conley closed loops in occurrence coordinates.
queries[].graph_loop_windows contains Morse-graph closed loops in occurrence coordinates.
queries[].holonomy_probes contains Conley triples [a,b,c] in occurrence coordinates.
queries[].graph_holonomy_probes contains Morse-graph triples [a,b,c] in occurrence coordinates.
Preserve all listed orders.

OUTPUT

The program is called as

python3 /app/trajectory_solver.py INPUT OUTPUT

It must write one regular UTF-8 JSON file at OUTPUT with this shape:

{
  "case_id": string,
  "queries": [
    {
      "query_id": int,
      "parameter_ids": int[],
      "seed_box_id": int,
      "selected_morse_ids_by_frame": int[][],
      "conley_node_dimensions_F2": [[int,...],[int,...],[int,...]],
      "conley_generalized_rank_queries": [
        {"window_id": int, "start_node": int, "end_node": int, "rank_F2": [int,int,int]}
      ],
      "conley_zigzag_barcodes_F2": [ [[int,int,int],...], [[int,int,int],...], [[int,int,int],...] ],
      "conley_loop_signatures": [
        {"loop_id": int, "start_occurrence": int, "end_occurrence": int,
         "degrees": [
           {"dimension": int, "relation_dimension_F2": int,
            "left_projection_rank_F2": int, "right_projection_rank_F2": int,
            "regular_rank_F2": int, "fixed_subspace_dimension_F2": int,
            "power_regular_ranks_F2": [int,int,int,int],
            "kronecker_invariants_F2": {"right_minimal_indices": int[], "left_minimal_indices": int[],
                                       "infinite_elementary_divisors": int[], "finite_elementary_divisors": [[int,int],...]}}
         ]}
      ],
      "conley_holonomy_word_signatures": [
        {"probe_id": int, "occurrences": [int,int,int],
         "degrees": [
           {"dimension": int, "space_dimension_F2": int,
            "word_signatures": [
              {"word": string, "relation_dimension_F2": int,
               "left_projection_rank_F2": int, "right_projection_rank_F2": int,
               "regular_rank_F2": int, "fixed_subspace_dimension_F2": int,
               "kronecker_invariants_F2": {"right_minimal_indices": int[], "left_minimal_indices": int[],
                                       "infinite_elementary_divisors": int[], "finite_elementary_divisors": [[int,int],...]}}
            ]}
         ]}
      ],
      "morse_graph_node_dimensions_F2": [[int,...],[int,...]],
      "morse_graph_generalized_rank_queries": [
        {"window_id": int, "start_node": int, "end_node": int, "rank_F2": [int,int]}
      ],
      "morse_graph_loop_signatures": [
        {"loop_id": int, "start_occurrence": int, "end_occurrence": int,
         "degrees": [
           {"dimension": int, "relation_dimension_F2": int,
            "left_projection_rank_F2": int, "right_projection_rank_F2": int,
            "regular_rank_F2": int, "fixed_subspace_dimension_F2": int,
            "kronecker_invariants_F2": {"right_minimal_indices": int[], "left_minimal_indices": int[],
                                       "infinite_elementary_divisors": int[], "finite_elementary_divisors": [[int,int],...]}}
         ]}
      ],
      "morse_graph_holonomy_word_signatures": [
        {"probe_id": int, "occurrences": [int,int,int],
         "degrees": [
           {"dimension": int, "space_dimension_F2": int,
            "word_signatures": [
              {"word": string, "relation_dimension_F2": int,
               "left_projection_rank_F2": int, "right_projection_rank_F2": int,
               "regular_rank_F2": int, "fixed_subspace_dimension_F2": int,
               "kronecker_invariants_F2": {"right_minimal_indices": int[], "left_minimal_indices": int[],
                                       "infinite_elementary_divisors": int[], "finite_elementary_divisors": [[int,int],...]}}
            ]}
         ]}
      ]
    }
  ]
}


Array orientation is part of the output contract. `selected_morse_ids_by_frame` has one outer entry per trajectory occurrence; each inner list is the sorted Morse-id selection at that occurrence. `conley_node_dimensions_F2` has exactly three outer lists in degree order H0, H1, H2; each of those lists has one integer per endpoint-bridge zigzag node. `conley_zigzag_barcodes_F2` has exactly three outer lists in degree order H0, H1, H2; each is the sorted interval list of that degree. `morse_graph_node_dimensions_F2` has exactly two outer lists in degree order H0, H1, again with one integer per zigzag node. Do not transpose these arrays.

The visible file `example_expected.json` is the expected output for `example.json`; it may be used to check both orientation and serialization before submission.

All dimensions, ranks, ids, and coordinates are nonnegative integers. Preserve query, window, loop, and degree order from the input. JSON key order is irrelevant. Do not emit NaN or infinities.
