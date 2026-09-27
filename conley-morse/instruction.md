Use `/app/data/example.json` as the visible practice input and `/app/data/example_expected.json` as its expected output. The complete scientific conventions, field definitions, array orientation, indexing, and JSON schema are in `/app/data/README.txt`. Hidden cases use the same contract with different grids, parameter families, trajectories, seeds, and query windows.

Create `/app/trajectory_solver.py`. The grader runs `python3 /app/trajectory_solver.py INPUT OUTPUT`. Read only `INPUT` and write one regular UTF-8 JSON file to `OUTPUT`. Python 3.11, `numpy==2.1.3`, `scipy==1.14.1`, and `networkx==3.4.2` are available. Each hidden invocation has a 1500 second wall-clock limit. Network access is unavailable.

The grader checks transported Morse selections, relative Conley homology dimensions, requested generalized ranks and complete zigzag barcodes on the Conley and Morse-graph endpoint-bridge zigzags, closed-loop continuation invariants, continuation-holonomy word signatures (including Kronecker invariants of every loop and word relation), the corresponding Morse reachability-graph homology quantities, and Conley index maps of Morse intervals.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
