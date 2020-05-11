#ifndef MGRAPH_HH
#define MGRAPH_HH

#include <assert.h>
#include <stdio.h>
#include <limits>
#include <vector>
#include <utility>
#include <functional>
#include <algorithm>

#include "edge.hh"
#include "int_util.hh"

/**
 * Minimalist graph implementation when vertices are ints from 0 to n-1.
 *
 * Example:
 *   typedef mgraph<int> graph;
 *   std:vector<graph::edge> edges;
 *   edges.push_back(graph::edge(1,2,100)); // edge 1 --> 2 with weight 100
 *   edges.push_back(graph::edge(2,3,200));
 *   graph g(4, edges); // 0 is also a vertex
 *   for (int u : g)
 *      for (auto e : g[u]) 
 *         std::cout << u << " " << e.dst << " " << e.wgt << std::endl;
 *
 *   // ignoring weights :
 *   for (int u : g)
 *      for (int v : g[u]) 
 *         std::cout << u << " " << v << std::endl;
 */


template<typename V, // vertex-number type
         typename W, // weight type
         V nb_not_vertex = std::numeric_limits<V>::max()>

class mgraph {
    // Minimalist graph as a flat array of dst,wgt pairs.

public:
    static const V not_vertex = nb_not_vertex;
    typedef V vertex;
    typedef W weight;
    typedef edge::dst_wgt<V,W> edge_head;  // A dst,wgt pair.
    
private:
    V n_;            // number of vertices
    std::vector<size_t> sdeg;    // prefix sum of degrees
    std::vector<edge_head> adj;  // neighbors of u are at index sdeg[u]
    
public:

    V n() const { return n_; }
    size_t m() const { return adj.size(); }

    typedef edge::src_dst_wgt<V,W> edge;

    // n must be greater than any vertex number (src or dst) in edg.
    mgraph(V n, const std::vector<edge> &edg) {
        init_from_edges(n, edg);
    }
    
    mgraph(const std::vector<edge> &edg) {
        set_edges(edg);
    }

    void set_edges(const std::vector<edge> &edg, V n = 0) {
        if (n == 0) {
            for (auto const e : edg)
                n = std::max(n, std::max(e.src, e.dst) + 1);
        }
        init_from_edges(n, edg);
    }

    mgraph() {
        static std::vector<edge> edg;
        init_from_edges(0, edg);
    }

    V degree(V u) const { return (int) sdeg[u+1] - sdeg[u]; }

    V max_degree() const {
        if (n_ <= 0) return 0;
        V d = degree(0);
        for (V u = 1; u < n_; ++u) {
            V du = degree(u);
            if (du > d) d = du;
        }
        return d;
    }
    
    edge_head neighbor(V u, size_t i) {
        size_t e = sdeg[u] + i;
        assert(e < sdeg[u+1]);
        return adj[e];
    }
    
    size_t degree_sum(V u) const { return sdeg[u]; }
    
    // asserts sorted adjacency lists (use g.reverse() or g.reverse().reverse())
    bool has_edge(V u, V v) const {
        size_t e1 = sdeg[u], e2 = sdeg[u+1];
        // is v in adj[e1 .. e2-1] ?
        while (e1 < e2) {
            size_t m = (e1 + e2) / 2;
            V w = adj[m].dst;
            if (w == v) return true;
            if (v < w) e2 = m; // in adj[e1 .. m-1]
            else e1 = m + 1; // in adj[m+1 .. e2-1]
        }
        return false;
    }

    // asserts sorted adjacency lists (use g.reverse() or g.reverse().reverse())
    W edge_weight(V u, V v) const {
        size_t e1 = sdeg[u], e2 = sdeg[u+1];
        // is v in adj[e1 .. e2-1] ?
        while (e1 < e2) {
            size_t m = (e1 + e2) / 2;
            V w = adj[m].dst;
            if (w == v) return adj[m].wgt;
            if (v < w) e2 = m; // in adj[e1 .. m-1]
            else e1 = m + 1; // in adj[m+1 .. e2-1]
        }
        throw std::invalid_argument("mgraph.edge_weight(): edge not found");
    }

    mgraph reverse() const {
        std::vector<edge> edg;
        edg.reserve(m());
        for (V u = 0; u < n_; ++u) {
            for (size_t e = sdeg[u]; e < sdeg[u+1]; ++e) {
                edg.push_back(edge(adj[e].dst, u, adj[e].wgt));
            }
        }
        return mgraph(n_, edg);
    }

    bool is_symmetric(bool same_weight = true) const {
        return asymmetry(same_weight) == 0;
    }

    size_t asymmetry(bool same_weight) const {
        mgraph r = reverse();
        size_t n_asym = 0;
        for (V u = 0; u < n_; ++u) {
            for (size_t e = sdeg[u]; e < sdeg[u+1]; ++e) {
                W w = adj[e].wgt;
                V v = adj[e].dst;
                if ( (! r.has_edge(u, v))
                     || (same_weight && r.edge_weight(u, v) != w) ) {
                    /*std::cerr <<"... "<< v <<","<< u <<"\n";
                    if (r.has_edge(u, v))
                    std::cerr << w <<" "<< r.edge_weight(u, v) <<"\n"; */
                    ++n_asym;
                } 
            }
        }
        return n_asym;
    }

    
    static W aggregate_min (W w, W x) { return std::min(w, x); }
    static W aggregate_sum (W w, W x) { return w + x; }
    
    mgraph simple(std::function<W(W,W)> aggr = aggregate_min) const {
        mgraph g = reverse().reverse(); // sort adjacencies
        std::vector<edge> edg;
        edg.reserve(m());
        for (V u = 0; u < g.n_; ++u) {
            for (size_t e = g.sdeg[u]; e < g.sdeg[u+1]; ++e) {
                W w = g.adj[e].wgt;
                V v = g.adj[e].dst;
                while (e+1 < g.sdeg[u+1] && g.adj[e+1].dst == v) {
                    ++e;
                    w = aggr(w, g.adj[e].wgt);
                }
                edg.push_back(edge(u, v, w));
            }
        }
        return mgraph(n_, edg);
    }

    std::vector<edge> edges_vector() const {
        std::vector<edge> edg;
        edg.reserve(m());
        size_t i = 0;
        for (V u = 0; u < n_; ++u) {
            for (size_t e = sdeg[u]; e < sdeg[u+1]; ++e) {
                edg.emplace_back(u, adj[e].dst, adj[e].wgt);
            }
        }
        return edg;
    }

    void sort_neighbors_by_weight() {
        for (V u = 0; u < n_; ++u) {
            size_t e = sdeg[u], f = sdeg[u+1];
            std::sort(adj.begin() + e, adj.begin() + f,
                      [](const edge_head &e, const edge_head &f) {
                          return e.wgt < f.wgt;
                      });
        }
    }

    bool is_ID_sorted() { // asymptotically slower than reverse() twice
        for (V u = 0; u < n_; ++u) {
            size_t e = sdeg[u], f = sdeg[u+1];
            for (size_t e = sdeg[u] + 1; e < sdeg[u+1]; ++e) {
                if (adj[e-1].dst > adj[e].dst) return false;
            }
        }
        return true;
    }

    std::pair<mgraph<V, W, nb_not_vertex>, std::vector<V> >
    subgraph(std::function<bool(V)> is_in, const std::vector<V> &vtx = {})
        const {
        bool use_vtx = vtx.size() > 0;
        if (use_vtx) assert(vtx.size() == n_);

        std::vector<V> vtx_sub, vtx_inv(n_, nb_not_vertex);
        V n_sub = 0;
        for (V u = 0; u < n_; ++u) { if (is_in(u)) ++n_sub; }
        vtx_sub.reserve(n_sub);
        for (V u = 0; u < n_; ++u) {
            if (is_in(u)) {
                vtx_inv[u] = vtx_sub.size();
                vtx_sub.push_back(use_vtx ? vtx[u] : u);
            }
        }
        
        std::vector<edge> edg;
        for (V u = 0; u < n_; ++u) {
            if (is_in(u)) {
                for (size_t e = sdeg[u]; e < sdeg[u+1]; ++e) {
                    V v = adj[e].dst;
                    if (is_in(v)) {
                        edg.push_back(edge(vtx_inv[u], vtx_inv[v], adj[e].wgt));
                    }
                }
            }
        }

        return std::make_pair(mgraph(n_sub, edg), vtx_sub);
    }

    
    //  --------------------- iterators : -----------------------
    //
    // for (V u : g)
    //     for (auto e : g[u])
    //         f(u, e.dst, e.wgt);
    //

    int_iterator<V> begin() const { return int_iterator<V>(0); }
    int_iterator<V> end() const { return int_iterator<V>(n_); }

    const mgraph& nodes() const { return *this; }
    range_rev<V> nodes_rev() const { return range_rev<V>(n_, 0); }

    
    class neighborhood {
        const mgraph &g;
        const V u;
     public:
        neighborhood(const mgraph &g, V u) : g(g), u(u) {}
        typedef typename std::vector<edge_head>::const_iterator eh_iterator;
        eh_iterator begin() const { return g.adj.cbegin() + g.sdeg[u]; }
        eh_iterator end() const { return g.adj.cbegin() + g.sdeg[u+1]; }
        const edge_head *b, *e;
    };

    neighborhood operator[](V u) const {
        if (u < 0 || u >= n_)
            throw std::invalid_argument("mgraph: not a vertex index: "
                                        + std::to_string(u));
        return neighborhood(*this, u);
    }
    neighborhood neighbors(V u) const { return (*this)[u]; }

    class neighborhood_rev {
    public:
        class rev_eh_iterator {
            const edge_head *eh;
        public:
            rev_eh_iterator(const  edge_head *eh) : eh(eh) {}
            const edge_head &operator*() const { return *eh; }
            rev_eh_iterator &operator++() { --eh; return *this; }
            bool operator!=(const rev_eh_iterator& o) { return eh != o.eh; }
        };
    private:
        rev_eh_iterator b, e;
    public:
        neighborhood_rev(const mgraph &g, V u)
            : b(g.adj.data() + g.sdeg[u+1] - 1),
              e(g.adj.data() + g.sdeg[u] - 1) { assert(u >= 0 && u < g.n_); }
        rev_eh_iterator begin() const { return b; }
        rev_eh_iterator end() const { return e; }
    };

    neighborhood_rev neighbors_rev(V u) const {
        return neighborhood_rev(*this, u);
    }

    class edg_iterator {
        const mgraph &g;
        V u;
        size_t e;
    public:
        edg_iterator(const mgraph &g, V u, size_t e) : g(g), u(u), e(e) {}
        edge operator*() const { return edge(u, g.adj[e].dst, g.adj[e].wgt); }
        edg_iterator &operator++() {
            ++e;
            if (e >= g.sdeg[u+1]) ++u;
            return *this;
        }
        bool operator!=(const edg_iterator& o) { return e != o.e; }
    };

    class edge_set {
        const mgraph &g;
    public:
        edge_set(const mgraph &g) : g(g) {}
        edg_iterator begin() const { return edg_iterator(g, 0, 0); }
        edg_iterator end() const { return edg_iterator(g, g.n_, g.m()); }
    };

    edge_set edges() const { return edge_set(*this); }
    
    
private:

    void init_from_edges(V n, const std::vector<edge> &edg) {
        n_ = n;
        size_t m_ = edg.size();

        // degrees:
        sdeg.resize(n_+1);
        for (V u = 0; u <= n_; ++u) sdeg[u] = 0;
        for (size_t i = 0; i < m_; ++i) {
            assert(0 <= edg[i].src && edg[i].src < n_);
            assert(0 <= edg[i].dst && edg[i].dst < n_);
            sdeg[edg[i].src] += 1;
        }
        for (V u = 1; u <= n_; ++u) {
            sdeg[u] += sdeg[u-1];
        }
        
        // adjacencies:
        adj.resize(m_);
        for (size_t i = m_; i > 0; ) {
            --i;
            size_t u = edg[i].src;
            size_t e = sdeg[u] - 1;
            adj[e].wgt = edg[i].wgt;
            adj[e].dst = edg[i].dst;
            sdeg[u] = e;
        }
    }

    
}; // mgraph




#endif // MGRAPH_HH
