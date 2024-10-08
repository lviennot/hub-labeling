#include <cstdlib>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>

#include "mgraph.hh"
#include "traversal.hh"
#include "pruned_landmark_labeling.hh"
#include "logging.hh"
#include "unit.hh"

typedef unsigned int uint;
typedef mgraph<uint, int64_t> graph;
typedef pruned_landmark_labeling<graph> pl_lab;
typedef mgraph<uint, pl_lab::hubinfo> graphL;

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
    
    std::cerr <<"\nUsage: "<< argv[0] <<" [command] [graph] [OPT [subset]]\n\n"
              << paragraph (
        "Computes a hub labeling of the weighted directed graph G in file "
        "[graph] and performs some operations based on a command (see list "
        "below and selected nodes listed in the optional file [subset] "
        "(all nodes if no file is specified).\n" )
              << paragraph (
        "With command 'hubs', it outputs the list "
        "of arcs from selected nodes to their out-hubs and from "
        "in-hubs to selected nodes." )
              << paragraph (
        "With command 'distances', it outputs the distances between selected "
        "nodes as triples 'u v dist_uv' for all pairs 'u, v' of selected nodes." )
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
              << paragraph (
        "With command 'benchmark', it computes a hub labeling and "
        "computes 1000 x 1000 distances between random nodes." )
              << paragraph (
        "With command 'rank', it computes a hub labeling and outputs the "
        "rank ordering used (most important hubs first)." )
              << paragraph (
        "With command 'stats_rank_threshold', it computes a hub labeling "
        "and outputs the size of labels if cut at rank [thr] for all [thr]." )
              << paragraph (
        "Command 'closure' computes G* and outputs its arcs." )
              << "\nInputs:\n"
        "A '-' for [graph] or [subset] stands for standard input.\n"
        " Graph format: one arc [src_id] [dst_id] [arc_length] per line\n"
        " Subset format: one node [id] per line.\n"
        " Command format: 'hubs', 'closure',... (see the full list above)\n"
              <<paragraph (
         " Output format: arcs "
         "[type] [node/hub id] [OPT next/hop id] [hub/node id] [length] "
         "where [type] is either 'i' or 'o' or 'c' for in-hub to node arc, "
         "or node to out-hub arc, or transitive closure arc respectively." );
        exit(1);
}


int main (int argc, char **argv) {
    logging main_log("--");

    // -------------------- options ----------------------
    auto i_arg = [&argc,&argv](std::string a) {
        //for (int i = 1; i < argc; ++i) std::cerr << argv[i] <<" "; std::cerr <<" for "<< a <<"\n";
        for (int i = 1; i < argc; ++i)
            if (a == argv[i])
                return i;
        return -1;
    };
    auto del_arg = [&argc,&argv,i_arg](std::string a) {
        int i = i_arg(a);
        if (i >= 0) {
            for (int j = i+1; j < argc; ++j)
                argv[j-1] = argv[j];
            --argc;
            return true;
        }
        return false;
    };

    bool symmetrize = del_arg("--symmetrize");

    // ------------------------ usage -------------------------
    std::string cmd(argc >= 2 ? argv[1] : "");
    if (argc < 3 || (cmd != "hubs" && cmd != "hubs-next-hop" && cmd != "distances"
                     && cmd != "test" && cmd != "benchmark" && cmd != "rank"
                     && cmd != "stats_rank_threshold" && cmd != "closure")) {
        usage_exit(argv);
    }

    // ------------------------ time -------------------------
    main_log.cerr() << "start\n";
    double t = main_log.lap();

    // ------------------------- load graph ----------------------
    std::unordered_map<std::string,int> vi; // vertex index
    std::vector<std::string> lab;
    graph g;
    bool weighted=false;
    size_t n = 0;
    {
        std::vector<graph::edge> edg;
        FILE *in = (std::string("-") != argv[2]) ? fopen(argv[2], "r") : stdin;
        char u[1024], v[1024];
        long long int w;
        for ( ; fscanf(in, " %s %s %lld \n", u, v, &w) >= 3 ; ) {
            if ( ! vi.contains(u)) { lab.push_back(u); vi[u] = n++; }
            if ( ! vi.contains(v)) { lab.push_back(v); vi[v] = n++; }
            edg.push_back(graph::edge(vi[u], vi[v], w));
            if(symmetrize) edg.push_back(graph::edge(vi[v], vi[u], w));
            if (w != 1) weighted = true;
            if (main_log.progress()) {
                main_log.cerr(t) << "read "<< edg.size() << " edges\n";
            }
        }
        // graph
        g.set_edges(edg);
    }
    g = g.simple();
    CHECK(n == g.n());
    size_t m = g.m();
    main_log.cerr(t)
        << "loaded graph with n=" << n << " nodes, m=" << m <<" edges " 
        << "weighted=" << weighted << "\n";
    t = main_log.lap();
    
    // ------------------------- load subset -----------------------
    std::vector<uint> sel;
    if (argc > 3) {
        FILE *in = (std::string("-") != argv[3]) ? fopen(argv[3], "r") : stdin;
        char u[1024];
        for ( ; fscanf(in, " %s \n", u) >= 1 ; ) {
            if (vi.contains(u)) {
                uint v = vi[u];
                sel.push_back(v);
            } else {
                std::cerr <<"  stop "<< u <<" not in the graph\n";
                exit(2);
            }
        }
    } else {
        for (uint u = 0; u < n; ++u) {
            sel.push_back(u);
        }
    } 
    std::vector<bool> is_sel(n, false);
    for (uint v : sel) { is_sel[v] = true; }
    uint n_sel = sel.size();
    main_log.cerr(t) << "loaded subset of "<< n_sel <<" nodes\n";
    t = main_log.lap();
    
    // ------------------------- hub labeling -----------------------
    pl_lab hl(g, {}, {}, 1, weighted);
    hl.print_stats(std::cerr, is_sel, is_sel);
    main_log.cerr(t) << "hub lab 1\n";
    t = main_log.lap();
    //pl_lab hl(g, 123, hl1.rank_order()); // does not help much
    //hl.print_stats(std::cerr, is_sel, is_sel);
    //main_log.cerr(t) << "hub lab 2\n";
    //t = main_log.lap();

    // ----------------------------- output ------------------------
    if (cmd == "distances") {
        for (uint u : sel) {
            for (uint v : sel) {
                int64_t dist_uv = hl.distance(u, v);
                std::cout << lab[u] <<" "<< lab[v] <<" "<< dist_uv <<"\n";
            }
        }
    } else if (cmd == "test") {
        // hub graphs
        std::vector<pl_lab::edgeL> edg_out = hl.out_hub_edges({}, {});
        std::vector<pl_lab::edgeL> edg_in = hl.in_hub_edges({}, {});
        graphL g_out(edg_out), g_in(edg_in);
        g_in = g_in.reverse();
        CHECK(g_out.is_ID_sorted());
        CHECK(g_in.is_ID_sorted());
        g = g.reverse().reverse(); // sort neighbors by ID
        traversal<graph> trav(g.n());
        std::vector<uint> src = {}, dst = {};
        for (uint i = 0; i < 100; ++i) src.push_back(rand() % g.n());
        for (uint i = 0; i < 100; ++i) dst.push_back(rand() % g.n());
        for (uint s : src) {
            trav.clear();
            trav.dijkstra(g, s);
            for (uint u : dst) {
                CHECK(trav.dist(u) == hl.distance(s, u));
                if (trav.dist(u) == INT64_MAX) continue;
                std::pair<uint, uint> hub = hl.common_hub(s, u);
                graphL::edge_head sx = g_out.neighbor(s, hub.first);
                graphL::edge_head xu = g_in.neighbor(u, hub.second);
                CHECK(sx.dst == xu.dst); // hub x
                uint x = hl.hub_ID(sx.dst);
                CHECK(trav.dist(u) == sx.wgt.dist + xu.wgt.dist);
                int64_t len = 0L;
                for (uint s2 = s; s2 != x; ) {
                    pl_lab::hubinfo hi = g_out.edge_weight(s2, sx.dst);
                    len += g.edge_weight(s2, hi.next_hop);
                    s2 = hi.next_hop;
                }
                CHECK(len == sx.wgt.dist);
                for (uint u2 = u; u2 != x; ) {
                    pl_lab::hubinfo hi = g_in.edge_weight(u2, xu.dst);
                    len += g.edge_weight(hi.next_hop, u2);
                    u2 = hi.next_hop;
                }
                CHECK(len == trav.dist(u));
            }
        }
    } else if (cmd == "benchmark") {
        std::vector<uint> src = {}, dst = {};
        for (uint i = 0; i < 1000; ++i) src.push_back(rand() % g.n());
        for (uint i = 0; i < 1000; ++i) dst.push_back(rand() % g.n());
        main_log.cerr(t) << "benchmark beg: "
                         << src.size() <<" x "<< dst.size() <<"\n";
        t = main_log.lap();
        uint dum_sum=0;
        for (uint s : src) {
            for (uint t : dst) {
                dum_sum += hl.distance(s, t);
            }
        }
        main_log.cerr(t) << "benchmark end: "<< dum_sum <<" "
                         << src.size() <<" x "<< dst.size() <<"\n";
        t = main_log.lap();
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
        CHECK(sel.size() == n); // Makes sense if next hops are also selected.
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
    } else if (cmd == "rank") {
        for (auto v : hl.rank_order()) {
            std::cout << lab[v] <<"\n";
        }
    } else if (cmd == "stats_rank_threshold") {
        hl.print_stats_rank_threshold(std::cout, is_sel, is_sel);
    } else if (cmd == "closure") {
        // hub graphs
        std::vector<pl_lab::edgeL> edg_out = hl.out_hub_edges(is_sel, is_sel);
        std::vector<pl_lab::edgeL> edg_in = hl.in_hub_edges(is_sel, is_sel);
        graphL g_out(edg_out), g_in(edg_in);
        main_log.cerr() << "hub graphs\n";
        // selection
        std::vector<uint> sel(n_sel), sel_inv(n);
        for (uint i = 0, u = 0; u < n; ++u) {
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
        uint u = 0;
        for ( ; u < n; ++u) {
            if (is_sel[u]) {
                uint i = sel_inv[u];
                for (auto e : g_out[u]) {
                    if (e.wgt.dist < INT64_MAX) {
                        for (auto f : g_in[e.dst]) {
                            if (f.wgt.dist < INT64_MAX) {
                                uint j = sel_inv[f.dst];
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
            for (graph::edge e : graph(edg).simple().edges()) {
                std::cout <<"c "<< lab[e.src]<<" "<< lab[e.dst]
                          <<" "<< e.wgt <<"\n";
            }
        } else {
            // adj matrix is smaller
            main_log.cerr() << "adj matrix\n";
            std::vector<std::vector<int64_t> > mat(n_sel);
            for (uint i = 0; i < n_sel; ++i) {
                mat[i].reserve(n_sel);
                for (uint j = 0; j < n_sel; ++j) {
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
                    uint i = sel_inv[u];
                    for (auto e : g_out[u]) {
                        if (e.wgt.dist < INT64_MAX) {
                            for (auto f : g_in[e.dst]) {
                                if (f.wgt.dist < INT64_MAX) {
                                    uint j = sel_inv[f.dst];
                                    int64_t d_ij = e.wgt.dist + f.wgt.dist;
                                    if (d_ij < mat[i][j]) mat[i][j] = d_ij; 
                                }
                            }
                        }
                    }
                }
            }
            main_log.cerr() << "closure\n";
            for (uint i = 0; i < n_sel; ++i) {
                for (uint j = 0; j < n_sel; ++j) {
                    std::cout <<"c "<< lab[sel[i]] <<" "<< lab[sel[j]]
                              <<" "<< mat[i][j] <<"\n";
                }
            }
        }
    } else {
        std::cerr << "Unkwnown command: "<< cmd <<"\n\n";
        usage_exit(argv);
    }
    
    main_log.cerr() << "end\n";
    exit(0);
}

