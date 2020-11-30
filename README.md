# Hub-labeling

This C++ code computes a hub labeling (aka 2-hop labeling) of a weighted directed graph. A hub-labeling is a compact representation of the distance matrix of the graph. The algorithm is based on [Robust Distance Queries on Massive Networks](https://link.springer.com/chapter/10.1007/978-3-662-44777-2_27) which itself is built upon [Fast Exact Shortest-Path Distance Queries on Large Networks by Pruned Landmark Labeling](https://arxiv.org/abs/1304.4661) for which an unweighted version is available at [pruned landmark labeling](https://github.com/iwiwi/pruned-landmark-labeling).

## Compile:
```
make
```

## Try it (with download of data):
```
make test
```

Author : Laurent Viennot, Inria 2020