# Hub-labeling

This C++ code computes a hub labeling (aka 2-hop labeling) of a weighted directed graph. A hub-labeling is a compact representation of the distance matrix of the graph. The algorithm is based on [Robust Distance Queries on Massive Networks](https://link.springer.com/chapter/10.1007/978-3-662-44777-2_27) which itself is built upon [Fast Exact Shortest-Path Distance Queries on Large Networks by Pruned Landmark Labeling](https://arxiv.org/abs/1304.4661) for which an unweighted version is available at [pruned landmark labeling](https://github.com/iwiwi/pruned-landmark-labeling).

Author : Laurent Viennot, Inria 2020

## Compile:
```
make
```

## Try it (with download of data):
```
make test
```

## Example: all distances within a subset of nodes

Graph example (a directed path with 6 nodes, all edges have length 1):
```
cat graph.txt
1 2 1
2 3 1
3 4 1
4 5 1
5 6 1
```

Selection example (nodes 1, 3, and 6):
```
cat selection.txt
1
3
6
```

All distances between selected nodes (INT64_MAX represents infinity):
```
_build/hl_trans distances graph.txt selection.txt 2> /dev/null
1 1 0
1 3 2
1 6 5
3 1 9223372036854775807
3 3 0
3 6 3
6 1 9223372036854775807
6 3 9223372036854775807
6 6 0
```




## Usage

```
Usage: _build/hl_trans [command] [graph] [OPT [subset]]

Computes a hub labeling of the weighted directed graph G in file [graph] and
 performs some operations based on a command (see list below and selected nodes
 listed in the optional file [subset] (all nodes if no file is specified).

With command 'hubs', it outputs the list of arcs from selected nodes to their
 out-hubs and from in-hubs to selected nodes.
With command 'distances', it outputs the distances between selected nodes as
 triples 'u v dist_uv' for all pairs 'u, v' of selected nodes.
The promise is that gathering transitive arcs obtained by following one arc of
 the 'out-hubs' list and one arc of the 'in-hubs' list provides the subgrap G*
 induced by selected nodes in the transitive closure of G.
With command 'hubs-next-hop', it outputs in addition the next hop to reach each
 hub.
With command 'test', it computes a hub labeling and checks with 100 random nodes
 that distances from these nodes are correct.
With command 'benchmark', it computes a hub labeling and computes 1000 x 1000
 distances between random nodes.
With command 'rank', it computes a hub labeling and outputs the rank ordering
 used (most important hubs first).
With command 'stats_rank_threshold', it computes a hub labeling and outputs the
 size of labels if cut at rank [thr] for all [thr].
Command 'closure' computes G* and outputs its arcs.

Inputs:
A '-' for [graph] or [subset] stands for standard input.
 Graph format: one arc [src_id] [dst_id] [arc_length] per line
 Subset format: one node [id] per line.
 Command format: 'hubs', 'closure',... (see the full list above)
 Output format: arcs [type] [node/hub id] [OPT next/hop id] [hub/node id]
 [length] where [type] is either 'i' or 'o' or 'c' for in-hub to node arc, or
 node to out-hub arc, or transitive closure arc respectively.
```