#include <iostream>

#include "mgraph.hh"

typedef mgraph<int, int> graph;
    
void mgraph_test(int n, int deg) {
    std::vector<graph::edge> edg;
    for (int u = 0; u < n; ++u) {
        for (int d = 0; d < deg; ++d) {
            edg.push_back(graph::edge(u, rand() % n, 1));
        }
    }
    std::cerr << "mgraph_test: scan graph\n";
    graph g(n, edg);
    for (int u : g) {
        for (int v : g[u]) std::cerr << u <<","<< v <<" ";
    }
    std::cerr << "\n";
    for (graph::edge e : g.edges()) {
        std::cerr << e.src <<","<< e.dst <<" ";
    }
    std::cerr << "\n";
    for (int u : g.nodes_rev()) {
        for (int v : g.neighbors_rev(u)) std::cerr << u <<","<< v <<" ";
    }
    std::cerr << "\n";

    // subgraph:
    auto sub = g.subgraph([](int u) { return u % 2 == 0; });
    graph h = sub.first;
    std::vector<int> vtx = sub.second;
    std::cerr << "mgraph_test: subgraph\n";
    for (int u : h) {
        for (int v : h[u]) {
            std::cerr << vtx[u] <<","<< vtx[v] <<" ";
        }
    }
    std::cerr << "\n";
    g = g.reverse().reverse();
    int hm = 0;
    for (int u : g) {
        for (int v : g[u])
            if (u % 2 == 0 && v % 2 == 0) ++hm;
    }
    assert(h.m() == hm);
    for (int u : h) {
        for (int v : h[u])
            assert(g.has_edge(vtx[u], vtx[v]));
    }

    // subgraph of subgraph:
    sub = h.subgraph([&vtx](int u) { return vtx[u] % 4 == 0; }, vtx);
    h = sub.first;
    vtx = sub.second;
    std::cerr << "mgraph_test: ";
    for (int u : h) {
        int iv = 0;
        for (int v : h[u]) {
            std::cerr << vtx[u] <<","<< vtx[v] <<" ";
        }
    }
    std::cerr << "\n";
    hm = 0;
    for (int u : g) {
        for (int v : g[u])
            if (u % 4 == 0 && v % 4 == 0) ++hm;
    }
    assert(h.m() == hm);
    for (int u : h) {
        for (int v : h[u])
            assert(g.has_edge(vtx[u], vtx[v]));
    }
}

    
