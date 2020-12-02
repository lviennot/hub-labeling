#include <iostream>

#include "mgraph.hh"
#include "unit.hh"

typedef mgraph<int, int> graph;
    
void mgraph_test(int n, int deg) {
    std::vector<graph::edge> edg;
    for (int u = 0; u < n; ++u) {
        for (int d = 0; d < deg; ++d) {
            edg.push_back(graph::edge(u, rand() % n, 1));
        }
    }
    std::cout << "mgraph_test: scan graph\n";
    graph g(n, edg);
    for (int u : g) {
        for (int v : g[u]) std::cout << u <<","<< v <<" ";
    }
    std::cout << "\n";
    for (graph::edge e : g.edges()) {
        std::cout << e.src <<","<< e.dst <<" ";
    }
    std::cout << "\n";
    for (int u : g.nodes_rev()) {
        for (int v : g.neighbors_rev(u)) std::cout << u <<","<< v <<" ";
    }
    std::cout << "\n";

    // subgraph:
    auto sub = g.subgraph([](int u) { return u % 2 == 0; });
    graph h = sub.first;
    std::vector<int> vtx = sub.second;
    std::cout << "mgraph_test: subgraph\n";
    for (int u : h) {
        for (int v : h[u]) {
            std::cout << vtx[u] <<","<< vtx[v] <<" ";
        }
    }
    std::cout << "\n";
    g = g.reverse().reverse();
    unsigned int hm = 0;
    for (int u : g) {
        for (int v : g[u])
            if (u % 2 == 0 && v % 2 == 0) ++hm;
    }
    CHECK(h.m() == hm);
    for (int u : h) {
        for (int v : h[u])
            CHECK(g.has_edge(vtx[u], vtx[v]));
    }

    // subgraph of subgraph:
    sub = h.subgraph([&vtx](int u) { return vtx[u] % 4 == 0; }, vtx);
    h = sub.first;
    vtx = sub.second;
    std::cout << "mgraph_test: ";
    for (int u : h) {
        for (int v : h[u]) {
            std::cout << vtx[u] <<","<< vtx[v] <<" ";
        }
    }
    std::cout << "\n";
    hm = 0;
    for (int u : g) {
        for (int v : g[u])
            if (u % 4 == 0 && v % 4 == 0) ++hm;
    }
    CHECK(h.m() == hm);
    for (int u : h) {
        for (int v : h[u])
            CHECK(g.has_edge(vtx[u], vtx[v]));
    }
}

    
