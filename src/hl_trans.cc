#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>

#include "mgraph.hh"
#include "traversal.hh"
#include "pruned_landmark_labeling.hh"
#include "logging.hh"

typedef mgraph<int, int64_t> graph;
typedef pruned_landmark_labeling<graph> pl_lab;
typedef mgraph<int, pl_lab::hubinfo> graphL;

void usage_exit (char **argv) {
    auto paragraph = [](std::string s, int width=80) -> std::string {
        std::string acc;
        while (s.size() > 0) {
            int pos = s.size();
            if (pos > width) pos = s.rfind(' ', width);
            std::string line = s.substr(0, pos);
            acc += line + "\n";
            s = s.substr(pos);
        }
        return acc;
    };
    
    std::cerr <<"Usage: "<< argv[0] <<" [command] [graph] [OPT [subset]]\n"
              << paragraph (
        "Computes a hub labeling of the graph G in file [graph] and "
        "outputs labels of selected nodes listed in the optional "
        "file [subset] (all nodes if no file is specified)." )
              << paragraph (
        "With command 'hubs', it outputs the list "
        "of arcs from selected nodes to their out-hubs and from "
        "in-hubs to selected nodes." )
              << paragraph (
        "The promise is that gathering transitive arcs obtained by "
        "following one arc of the 'out-hubs' list and one arc of the "
        "'in-hubs' list provides the subgrap G* induced by selected nodes "
        "in the transitive closure of G." )
              << paragraph (
        "With command 'hubs-next-hop', it outputs in addition the next "
        "hop to reach each hub." )
              << paragraph (
        "With command 'test', it computes a hub labeling and checks with "
        "100 random nodes that distances from these nodes are correct." )
              <<
        "Command 'closure' computes G* and outputs its arcs.\n"
        "A '-' for [graph] or [subset] stands for standard input.\n"
        " Graph format: one arc [src_id] [dst_id] [length] per line\n"
        " Subset format: one node [id] per line.\n"
        " Command format: 'hubs' or 'closure'\n"
              <<paragraph (
         " Output format: arcs "
         "[type] [node/hub id] [OPT next/hop id] [hub/node id] [length] "
         "where [type] is either 'i' or 'o' or 'c' for in-hub to node arc, "
         "or node to out-hub arc, or transitive closure arc respectively." );
        exit(1);
}


int main (int argc, char **argv) {
    logging main_log("--");

    // ------------------------ usage -------------------------
    std::string cmd(argc >= 2 ? argv[1] : "");
    if (argc < 3 || (cmd != "hubs" && cmd != "hubs-next-hop"
                     && cmd != "test" && cmd != "closure")) {
        usage_exit(argv);
    }

    // ------------------------ time -------------------------
    main_log.cerr() << "start\n";
    double t = main_log.lap();

    // ------------------------- load graph ----------------------
    std::unordered_map<std::string,int> vi; // vertex index
    std::vector<std::string> lab;
    graph g;
    size_t n = 0;
    {
        std::vector<graph::edge> edg;
        FILE *in = (std::string("-") != argv[2]) ? fopen(argv[2], "r") : stdin;
        char u[1024], v[1024];
        long long int w;
        for ( ; fscanf(in, " %s %s %lld \n", u, v, &w) >= 3 ; ) {
            if (vi[u] == 0) { lab.push_back(u); vi[u] = 1+n++; }
            if (vi[v] == 0) { lab.push_back(v); vi[v] = 1+n++; }
            edg.push_back(graph::edge(vi[u]-1, vi[v]-1, w));
            //if(symmetrize) edg.push_back(graph::edge(vi[v]-1, vi[u]-1, w));
            if (main_log.progress()) {
                main_log.cerr(t) << "read "<< edg.size() << " edges\n";
            }
        }
        // graph
        g.set_edges(edg);
    }
    g = g.simple();
    assert(n == g.n());
    size_t m = g.m();
    main_log.cerr(t)
        << "loaded graph with n=" << n << " nodes, m=" << m <<" edges\n";
    t = main_log.lap();
    
    // ------------------------- load subset -----------------------
    std::vector<int> sel;
    if (argc > 3) {
        FILE *in = (std::string("-") != argv[3]) ? fopen(argv[3], "r") : stdin;
        char u[1024];
        for ( ; fscanf(in, " %s \n", u) >= 1 ; ) {
            if (vi[u] != 0) {
                int v = vi[u] - 1;
                sel.push_back(v);
            } else {
                std::cerr <<"  stop "<< u <<" not in the graph\n";
            }
        }
    } else {
        for (int u = 0; u < n; ++u) {
            sel.push_back(u);
        }
    } 
    std::vector<bool> is_sel(n, false);
    for (int v : sel) { is_sel[v] = true; }
    int n_sel = sel.size();
    main_log.cerr(t) << "loaded subset of "<< n_sel <<" nodes\n";
    t = main_log.lap();
    
    // ------------------------- hub labeling -----------------------
    pl_lab hl(g);
    hl.print_stats(std::cerr, is_sel, is_sel);
    main_log.cerr(t) << "hub lab\n";
    t = main_log.lap();

    // ----------------------------- output ------------------------
    if (cmd == "test") {
        // hub graphs
        std::vector<pl_lab::edgeL> edg_out = hl.out_hub_edges(is_sel, is_sel);
        std::vector<pl_lab::edgeL> edg_in = hl.in_hub_edges(is_sel, is_sel);
        graphL g_out(edg_out), g_in(edg_in);
        g_in = g_in.reverse();
        assert(g_out.is_ID_sorted());
        assert(g_in.is_ID_sorted());
        g = g.reverse().reverse(); // sort neighbors by ID
        traversal<graph> trav(g.n());
        std::vector<int> src = {};
        for (int i = 0; i < 100; ++i) src.push_back(rand() % g.n());
        for (int s : src) {
            trav.clear();
            trav.dijkstra(g, s);
            for (int u : g) {
                assert(trav.dist(u) == hl.distance(s, u));
                if (trav.dist(u) == INT64_MAX) continue;
                std::pair<int, int> hub = hl.common_hub(s, u);
                graphL::edge_head sx = g_out.neighbor(s, hub.first);
                graphL::edge_head xu = g_in.neighbor(u, hub.second);
                assert(sx.dst == xu.dst); // hub x
                int x = hl.hub_ID(sx.dst);
                assert(trav.dist(u) == sx.wgt.dist + xu.wgt.dist);
                int64_t len = 0L;
                for (int s2 = s; s2 != x; ) {
                    pl_lab::hubinfo hi = g_out.edge_weight(s2, sx.dst);
                    len += g.edge_weight(s2, hi.next_hop);
                    s2 = hi.next_hop;
                }
                assert(len == sx.wgt.dist);
                for (int u2 = u; u2 != x; ) {
                    pl_lab::hubinfo hi = g_in.edge_weight(u2, xu.dst);
                    len += g.edge_weight(hi.next_hop, u2);
                    u2 = hi.next_hop;
                }
                assert(len == trav.dist(u));
            }
        }
    } else if (cmd == "hubs") {
        std::vector<pl_lab::edgeL> edg; 
        edg = hl.in_hub_edges(is_sel, is_sel);
        for (const pl_lab::edgeL &e : edg) {
            std::cout <<"i "<< lab[e.wgt.hub]
                      <<" "<< lab[e.dst] <<" "<< e.wgt.dist <<"\n";
        }
        edg = hl.out_hub_edges(is_sel, is_sel);
        for (const pl_lab::edgeL &e : edg) {
            std::cout <<"o "<< lab[e.src]
                      <<" "<< lab[e.wgt.hub] <<" "<< e.wgt.dist <<"\n";
        }
    } else if (cmd == "hubs-next-hop") {
        assert(sel.size() == n); // Makes sense if next hops are also selected.
        std::vector<pl_lab::edgeL> edg; 
        edg = hl.in_hub_edges(is_sel, is_sel);
        for (const pl_lab::edgeL &e : edg) {
            std::cout <<"i "<< lab[e.wgt.hub] <<" "<< lab[e.wgt.next_hop]
                      <<" "<< lab[e.dst] <<" "<< e.wgt.dist <<"\n";
        }
        edg = hl.out_hub_edges(is_sel, is_sel);
        for (const pl_lab::edgeL &e : edg) {
            std::cout <<"o "<< lab[e.src] <<" "<< lab[e.wgt.next_hop]
                      <<" "<< lab[e.wgt.hub] <<" "<< e.wgt.dist <<"\n";
        }
    } else if (cmd == "closure") {
        // hub graphs
        std::vector<pl_lab::edgeL> edg_out = hl.out_hub_edges(is_sel, is_sel);
        std::vector<pl_lab::edgeL> edg_in = hl.in_hub_edges(is_sel, is_sel);
        graphL g_out(edg_out), g_in(edg_in);
        main_log.cerr() << "hub graphs\n";
        // selection
        std::vector<int> sel(n_sel), sel_inv(n);
        for (int i = 0, u = 0; u < n; ++u) {
            if (is_sel[u]) {
                sel[i] = u;
                sel_inv[u] = i;
                ++i;
            } else {
                sel_inv[u] = graph::not_vertex;
            }
        }
        main_log.cerr() << "selection\n";
        // trans arcs
        std::vector<graph::edge> edg;
        int u = 0;
        for ( ; u < n; ++u) {
            if (is_sel[u]) {
                int i = sel_inv[u];
                for (auto e : g_out[u]) {
                    if (e.wgt.dist < INT64_MAX) {
                        for (auto f : g_in[e.dst]) {
                            if (f.wgt.dist < INT64_MAX) {
                                int j = sel_inv[f.dst];
                                int64_t d_ij = e.wgt.dist + f.wgt.dist;
                                edg.push_back(graph::edge(i, j, d_ij));
                            }
                        }
                    }
                }
                if (edg.size() > 2 * n_sel * n_sel) break;
            }
        }
        if (u >= n) {
            main_log.cerr() << "closure\n";
            for (const graph::edge &e : graph(edg).simple().edges()) {
                std::cout <<"c "<< lab[e.src]<<" "<< lab[e.dst]
                          <<" "<< e.wgt <<"\n";
            }
        } else {
            // adj matrix is smaller
            main_log.cerr() << "adj matrix\n";
            std::vector<std::vector<int64_t> > mat(n_sel);
            for (int i = 0; i < n_sel; ++i) {
                mat[i].reserve(n_sel);
                for (int j = 0; j < n_sel; ++j) {
                    mat[i].push_back(INT64_MAX);
                }
            }
            // reuse found arcs
            for (const graph::edge &e : edg) {
                if (e.wgt < mat[e.src][e.dst]) mat[e.src][e.dst] = e.wgt;
            }
            edg = {}; // release memory
            // end closure
            for ( ; u < n; ++u) {
                if (is_sel[u]) {
                    int i = sel_inv[u];
                    for (auto e : g_out[u]) {
                        if (e.wgt.dist < INT64_MAX) {
                            for (auto f : g_in[e.dst]) {
                                if (f.wgt.dist < INT64_MAX) {
                                    int j = sel_inv[f.dst];
                                    int64_t d_ij = e.wgt.dist + f.wgt.dist;
                                    if (d_ij < mat[i][j]) mat[i][j] = d_ij; 
                                }
                            }
                        }
                    }
                }
            }
            main_log.cerr() << "closure\n";
            for (int i = 0; i < n_sel; ++i) {
                for (int j = 0; j < n_sel; ++j) {
                    std::cout <<"c "<< lab[sel[i]] <<" "<< lab[sel[j]]
                              <<" "<< mat[i][j] <<"\n";
                }
            }
        }
    } else {
        usage_exit(argv);
    }
    
    main_log.cerr() << "end\n";
    exit(0);
}

