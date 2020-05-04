#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <assert.h>
#include <stdint.h>
#include <climits>
#include <vector>
#include <queue>          // std::priority_queue
#include <algorithm>      // std::min
#include <functional>
#include <unordered_map>

#include "heap.hh"
#include "mgraph.hh"

/**
 * Classical graph traversals: BFS, Dijkstra, ...
 * Example :
 *    mgraph<> g;
 *    traversal<graph> trav(n);
 *    int u=0; v=1;
 *    trav.dijkstra(g, u);
 *    int dist_u_v = trav.dist(v);
 *    trav.bfs(g, u);
 *    int hops_u_v = trav.dist(v);
 */

enum dfs_option { NONE = 0, SCC = 1, CHECK_DAG = 2 };

template<typename G,  // Graph type with G::vertex = int
         typename WL = int64_t, // long weight for summing weights, 
                     // with appropriate casting for the following constants:
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class traversal {
public:
    typedef typename G::vertex V;
    typedef typename G::weight W;
    typedef WL long_weight;
private:
    
    std::vector<V> q_vec_, scc_stack_;
    std::vector<V> visit_;
    std::vector<bool> visited_, in_scc_stack_;
    std::vector<int> visited_at_;
    std::vector<int> lowlink_, visit_end_; // for DFS
    std::vector<V> parent_; 
    std::vector<WL> dist_, tree_ecc_; // for Dijkstra and BFS
    heap queue_;
    std::vector<int> size_; // for tree centroid, and component sizes
    
    int n_, nvis_;
    
public:

    const V not_vertex = G::not_vertex;

    // Handle graphs with vertices in [0..n-1].
    traversal (int n)
        : dist_(n, max_weight),
          queue_([this](V u, V v){ return dist_[u] < dist_[v]; }, n), 
          visit_(n), visited_(n, false), visited_at_(n, n),
          size_(n, 0), parent_(n, n), n_(n), nvis_(0),
          in_scc_stack_(n, false), lowlink_(n, n), visit_end_(n, n),
          tree_ecc_(n, zero_weight)
    {
        q_vec_.reserve(n_);
    }

    int n() const { return n_; }
    int nvis() const { return nvis_; }
    V visit(int i) const { assert(i < nvis_); return visit_[i]; }
    bool visited(V u) { return visited_[u]; }
    int visited_at(V u) { return visited_at_[u]; }
    V parent(V u) const { return parent_[u]; }
    WL dist(V u) const { return dist_[u]; }
    const std::vector<WL>& distances() const { return dist_; }
    int size(V u) const { return size_[u]; }
    WL tree_ecc(V u) const { return tree_ecc_[u]; }
    V first_visited (int i_start = 0) const { return visit_[i_start]; }
    V last_visited () const { return visit_[nvis_ - 1]; }

    typedef typename G::edge edge;
    
    G graph() const {
        int n = nvis_;
        std::vector<edge> edg;
        edg.reserve(2*n);
        for (int i = 0; i < n; ++i) {
            V u = visit_[i];
            int p = visited_at_[parent_[u]];
            edg.push_back(edge(i, p, dist_[u] - dist_[p]));
            edg.push_back(edge(p, i, dist_[u] - dist_[p]));
        }
        return G(n, edg);
    }

    G digraph_to_sons(std::function<W(WL)> dist_to_weight
                      = [](WL w) { return (W)w; }) const { // FIXME !!!
        int n = nvis_;
        std::vector<edge> edg;
        edg.reserve(n);
        for (int i = 1; i < n; ++i) {
            V u = visit_[i];
            int p = visited_at_[parent_[u]];
            edg.push_back(edge(p, i, dist_to_weight(dist_[u] - dist_[p])));
        }
        return G(n, edg);
    }


    struct tree_node { // for compact representation of a tree
        int vtx; // associated vertex in the graph (of the tree of the node)
        int parent; // visit number of the parent
        //WL dist; // distance from root to vtx
        //WL dist_max; // distance from root of furthest descendant of vtx
        int size; // number of nodes in subtree
        tree_node()
            : vtx(-1), parent(-1), //dist(0), dist_max(zero_weight),
              size(0) {}
    };

    typedef tree_node node;

    std::vector<node> digraph_nodes() const { // node i cooresponds to visit_[i]
        std::vector<node> nds(nvis_);
        for (int i = nvis_ - 1; i != -1; --i) {
            V u = visit_[i];
            nds[i].vtx = u;
            int p = visited_at_[parent_[u]];
            nds[i].parent = p;
            //nds[i].dist = dist_[u];
            //nds[i].dist_max = std::max(nds[i].dist_max, dist_[u]);
            nds[i].size += 1;
            if (p != i) {
                //nds[p].dist_max = std::max(nds[p].dist_max, nds[i].dist_max);
                nds[p].size += nds[i].size;
            }
        }
        return nds;
    }

    
    typedef std::unordered_map<V,int> index;
    
    index digraph_index() const {
        index idx;
        for (int i = 0; i < nvis_; ++i) {
            V u = visit_[i];
            idx[u] = i;
        }
        return idx;
    }


    void clear(WL dft_wgt = max_weight, int n = 0) {
        int up_to = std::max(n, nvis_);
        for (int i = 0; i < nvis_; ++i)  {
            V u = visit_[i];
            visit_[i] = not_vertex;
            visited_[u] = false;
            visited_at_[u] = n_;
            in_scc_stack_[u] = false;
            lowlink_[u] = n_;
            visit_end_[u] = n_;
            dist_[u] = dft_wgt;
            parent_[u] = not_vertex;
            size_[u] = 0;
            tree_ecc_[u] = zero_weight;
        }
        scc_stack_.clear();
        nvis_ = 0;
    }

    static bool visit_all(V v, WL d, V p, WL dp) { return true; }
    
    // returns number of visited nodes
    int bfs(const G &g, V s,
            std::function<bool(V, WL, V, WL)> filtr
               = [](V v, WL d, V p, WL dp) { return true; },
            std::vector<V> more_sources = {}, WL step_weight = 1) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        int head = 0, tail = 0;
        more_sources.push_back(s);
        for (int s : more_sources) {
            dist_[s] = zero_weight;
            parent_[s] = s;
            q_vec_[tail++] = s;
        }

        while (head < tail) {
            V u = q_vec_[head++];
            WL du = dist_[u];
            if ( ! visited_[u]) {
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                WL dv = du + step_weight;
                for (V v : g[u]) {
                    if (dist_[v] == max_weight && filtr(v, dv, u, du)) {
                        dist_[v] = dv;
                        parent_[v] = u;
                        q_vec_[tail++] = v;
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    static bool filter_all(V v, WL d, V p, WL dp) { return true; }
    
    // returns number of visited nodes
    int dijkstra (const G &g, V s,
                  std::function<bool(V, WL, V, WL)> filtr
                  = [](V v, WL d, V p, WL dp) { return true; },
                  std::vector<V> more_sources = {},
                  std::vector<WL> more_sources_dist = {},
                  std::function<WL(WL,V,V,W)> add_arc_weight
                  = [](WL du, V u, V v, W w) { return du + (WL)w; }) { // FIXME !!!
        assert(g.n() <= n_);
        queue_.clear();
        queue_.set_compare([this](V u, V v){ return dist_[u] < dist_[v]; });
        int nvis0 = nvis_;
        more_sources.push_back(s);
        for (int i = 0; i < more_sources.size(); ++i) {
            V s = more_sources[i];
            if (s != not_vertex) {
                dist_[s] = i < more_sources_dist.size() ? more_sources_dist[i]
                               : zero_weight;
                parent_[s] = s;
                queue_.push(s);
            }
        }
            
        while ( ! queue_.empty()) {
            V u = queue_.pop();
            WL du = dist_[u];
            if (! visited_[u]) {
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                for (auto e : g[u]) {
                    V v = e.dst;
                    WL dv = add_arc_weight(du, u, v, e.wgt);
                    if((! visited_[v]) && dv < dist_[v] && filtr(v, dv, u, du)){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(v);
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // returns number of visited nodes
    int dijkstra_weight (const G &g, V s, WL ds = zero_weight,
                  std::function<WL(WL,V,V,W)> add_arc_weight
                  = [](WL du, V u, V v, W w) { assert(false); },
                  std::function<bool(V, WL, V, WL)> filtr
                  = [](V v, WL d, V p, WL dp) { return true; },
                  std::vector<V> more_sources = {},
                  std::vector<WL> more_sources_dist = {}) {
        assert(g.n() <= n_);
        queue_.clear();
        queue_.set_compare([this](V u, V v){ return dist_[u] < dist_[v]; });
        int nvis0 = nvis_;
        more_sources.push_back(s);
        for (int i = 0; i < more_sources.size(); ++i) {
            V s = more_sources[i];
            if (s != not_vertex) {
                dist_[s] = i < more_sources_dist.size() ? more_sources_dist[i]
                               : ds;
                parent_[s] = s;
                queue_.push(s);
            }
        }
            
        while ( ! queue_.empty()) {
            V u = queue_.pop();
            WL du = dist_[u];
            if (! visited_[u]) {
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                for (auto e : g[u]) {
                    V v = e.dst;
                    WL dv = add_arc_weight(du, u, v, e.wgt);
                    if((! visited_[v]) && dv < dist_[v] && filtr(v, dv, u, du)){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(v);
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // returns number of visited nodes
    int a_star (const G &g, V s, V t,
                std::function<WL(V)> dist_to_t_lb,
                std::function<bool(V, WL, V, WL)> filtr
                = [](V v, WL d, V p, WL dp) { return true; }) {
        assert(g.n() <= n_);
        queue_.clear();
        queue_.set_compare([this, dist_to_t_lb](V u, V v){
                return dist_[u] + dist_to_t_lb(u) < dist_[v] + dist_to_t_lb(v);
            });
        int nvis0 = nvis_;
        dist_[s] = zero_weight;
        parent_[s] = s;
        queue_.push(s);

        while ( ! queue_.empty()) {
            V u = queue_.pop();
            if (! visited_[u]) {
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                if (u == t) { break; }
                for (auto e : g[u]) {
                    WL dv = dist_[u] + e.wgt;
                    V v = e.dst;
                    if((! visited_[v]) && dv < dist_[v]
                       && filtr(v, dv, u, dist_[u])){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(v);
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    void clear_a_star(const G &g, WL dft_wgt = max_weight, int n = 0) {
        int up_to = std::max(n, nvis_);
        for (int i = 0; i < nvis_; ++i)  {
            V u = visit_[i];
            visit_[i] = not_vertex;
            visited_[u] = false;
            visited_at_[u] = n_;
            lowlink_[u] = n_;
            visit_end_[u] = n_;
            dist_[u] = dft_wgt;
            parent_[u] = not_vertex;
            size_[u] = 0;
            tree_ecc_[u] = zero_weight;
            for (V v : g[u]) {
                dist_[v] = dft_wgt;
                parent_[v] = not_vertex;
            }
        }
        nvis_ = 0;
    }


    // returns number of visited nodes
    int dfs (const G &g, V s, dfs_option opt=NONE) {
        assert(g.n() <= n_);
        int nvis0 = nvis_, vis_end = nvis_;
        parent_[s] = s;
        int tail = 0;
        q_vec_[tail++] = s;

        while (tail > 0) {
            V u = q_vec_[--tail];
            /* lowlink_[u] (see Tarjan alg. for strongly connected components)
             * is the smallest visit number of a node accessible from u through
             * a sequence of forward edges plus eventually one backward edge. 
             * Here we consider more paths, possibly with several backward
             * edges; the important point is that all paths considered by
             * Tarjan algorithm are considered here, and the first node visited
             * in a strongly connected component (scc) will thus be the only 
             * node of the component that has its own visit  number as lowlink 
             * here also. (Not storing the visit number of nodes allows to
             * avoid allocating one more array).
             */
            if ( ! visited_[u]) { // begin visit of [u]
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                lowlink_[u] = i_vis;
                if (opt == SCC) {
                    scc_stack_.push_back(u);
                    in_scc_stack_[u] = true;
                }
                
                // add u to stack for detecting end of visit when popped again
                if (tail >= q_vec_.size()) q_vec_.push_back(not_vertex);
                q_vec_[tail++] = u;
                
                for (V v : g[u]) {
                    if ( ! visited_[v]) {
                        // possible forward edge
                        parent_[v] = u; // possibly overwritten
                        if (tail >= q_vec_.size()) q_vec_.push_back(not_vertex);
                        q_vec_[tail++] = v;
                    } else if (visited_[v]) {
                        // backward edge
                        if (opt == SCC && in_scc_stack_[v]) {
                            //lowlink_[u] = std::min(lowlink_[u], visit_beg[v]);
                            lowlink_[u] = std::min(lowlink_[u], lowlink_[v]);
                        } else if (opt == CHECK_DAG && visit_end_[v] == n_) {
                            throw "cycle";
                        }
                    }
                }
            } else if (visit_end_[u] == n_) { // end visit of [u]
                visit_end_[u] = vis_end++;
                if (opt == SCC) {
                    V p = parent_[u];
                    int lowlink_u = lowlink_[u];
                    lowlink_[p] = std::min(lowlink_[p], lowlink_u);
                    assert(0 <= lowlink_u && lowlink_u < nvis_);
                    if (visit_[lowlink_u] == u) {
                        // u is the first visited node in its scc,
                        // pop scc and
                        // set all lowlinks in the scc to u's visit nb
                        while ( ! scc_stack_.empty() ) {
                            V v = scc_stack_.back();
                            scc_stack_.pop_back();
                            lowlink_[v] = lowlink_u;
                            in_scc_stack_[v] = false;
                            if (v == u) break; // all scc scanned
                        }
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // Compute strongly connected components of the input graph and returns
    // the number of components found.
    int strongly_connected_components(const G &g) {
        clear();
        for (V s = 0; s < g.n(); ++s) {
            if ( ! visited_[s]) {
                dfs(g, s, SCC);
            }
        }
        int nb = 0;
        for (int i = 0; i < nvis_; ++i) {
            V u = visit_[i];
            int ll = lowlink_[u];
            assert(0 <= ll && ll < nvis_);
            ++(size_[ll]);
            if (ll == i) ++nb;
        }
        return nb;
    }

    // Returns the scc number of the strongly connected component of [v],
    // asserts that [strongly_connected_components(g)] has been called
    // and that [v <= g.n()].
    // The return number is in the range 0 .. g.n() - 1.
    int scc_number(V v) {
        assert(v < nvis_);
        return lowlink_[v];
    }
    int scc_largest() {
        assert(nvis_ > 0);
        int m = 0;
        for (int i = 1; i < nvis_; ++i) {
            if (size_[i] > size_[m]) m = i;
        }
        return m;
    }
    int scc_size(int i) { return size_[i]; }
    std::vector<int> scc_vector() {
        return lowlink_;
    }
    // Returns a node of scc number i.
    V scc_node(int i) {
        assert(i < nvis_);
        return visit_[i];
    }
    
    // Returns a topological ordering [ord] of [g]
    // such that for all edge [uv] of [g], [u] appears before [v] in [ord].
    std::vector<int> topological_ordering (G &g) {
        clear();
        for (V u : g) {
            if ( ! visited_[u]) {
                dfs(g, u, true);
            }
        }
        int n = g.n();
        std::vector<int> ord;
        ord.reserve(n);
        for (V u : g) ord[n - 1 - visit_end_[u]] = u;
        return ord;
    }


    
    /** The following functions assume that WL is an int type. */
    
    // returns number of visited nodes, assumes 
    int max_card_search (const G &g, V s,
                 std::function<bool(V, V)> filtr
                 = [](V v, V p) { return true; }) {
        assert(g.n() <= n_);
        queue_.clear();
        // number of visited neighbors of u is stored in dist_[u]
        queue_.set_compare([this](V u, V v){ return dist_[u] > dist_[v]; });
        int nvis0 = nvis_;
        dist_[s] = 0;
        parent_[s] = s;
        queue_.push(s);

        while ( ! queue_.empty()) {
            V u = queue_.pop();
            if (! visited_[u]) {
                int i_vis = nvis_++;
                visit_[i_vis] = u;
                visited_[u] = true;
                visited_at_[u] = i_vis;
                for (V v : g[u]) {
                    if ( (! visited_[v]) && filtr(v, u)) {
                        dist_[v] += 1;
                        queue_.push(v);
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }


    // Asserts that the graph is a tree (with symmetric links).
    int tree_size(const G &g, int root,
                  std::function<int(V)> node_size = [](V v) { return 1; }) {
        return tree_size(g, root, root);
    }
    int tree_size(const G &g, int u, int par,
                  std::function<int(V)> node_size = [](V v) { return 1; }) {
        int s = 1; // for u
        size_[u] = node_size(u);
        assert (size_[u] > 0);
        for (V v : g[u]) {
            if (size_[v] <= 0) {
                s += tree_size(g, v, u);
            } else assert(v == par || v == u); // Check that we are in a tree.
        }
        size_[u] = s;
        return s ;
    }

    void tree_size_clear(const G &g, int r) { tree_size_clear(g, r, r); }
    void tree_size_clear(const G &g, int u, int par) {
        if (size_[u] > 0) {
            size_[u] = 0;
            for (V v : g[u]) {
                if (size_[v] > 0) {
                    tree_size_clear(g, v, u);
                } else assert(v == par || v == u); // Check tree.
            }
        }
    }

    // Asserts that tree_size(g, root) has been called.
    V centroid(const G &g, int r, int n_threshold) {
        return centroid(g, r, r, n_threshold);
    }
    V centroid(const G &g, int r, int par, int n_threshold) {
        for (V v : g[r]) {
            if (v != par && v != r && size_[v] > n_threshold)
                return centroid (g, v, r,  n_threshold);
        }
        return r;
    }

    

}; // traversal

#endif // TRAVERSAL_HH
