#ifndef HUB_LABELING_HH
#define HUB_LABELING_HH

#include <assert.h>
#include <vector>
#include <queue>          // std::priority_queue
#include <set>

#include "heap.hh"
#include "traversal.hh"
#include "pruned_landmark_labeling.hh"
#include "logging.hh"

/** Compute hub labels (alias a two-hop labeling) using random sampling and
 *  pruned landmark labeling. */

template<typename G,
         typename WL = int64_t, // long weigth (not necessarily T::long_weigth)
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class hub_labeling {
private:
    typedef typename G::vertex V;
    typedef typename G::weight W;
    typedef traversal<G, WL, max_weight, zero_weight> trav_G;
    typedef pruned_landmark_labeling<G, WL, max_weight, zero_weight> pll_G;

    std::vector<G> trees;
    std::vector<std::vector<typename trav_G::node> > nodes;
    std::vector<std::vector<int> > in_trees, in_trees_index;

    int n_, desired_samples;
    int64_t n_all_trees, i_sample, avg_tree_size, sum_all_trees;
    std::vector<V> sample_fwd, sample_bwd;
    std::vector<int> n_labs; // how many labels of random sample would have i
    std::vector<int64_t> n_pairs; // how many pairs from/to random sample would be covered by i
    trav_G trav;
    pll_G pll;
    std::vector<bool> is_hub;
    heap queue;
    inline bool covers_more (V u, V v) {
        //return n_pairs[u] > n_pairs[v];
        assert(n_labs[u] != 0 && n_labs[v] != 0);
        return n_pairs[u] * n_labs[v] > n_pairs[v] * n_labs[u];
    }

    static const int const_log_n = 4; // 2->10 ==> max lab -= 20%, 2->1 => +20%
    
public:

    hub_labeling(const G &g, bool weighted = true)
        : n_(g.n()), i_sample(1), avg_tree_size(g.n()), sum_all_trees(0),
          sample_fwd(g.n()), sample_bwd(g.n()), is_hub(g.n(), false),
          queue([this](int u, int v){ return covers_more(u, v); }, g.n()),
          in_trees(g.n()), in_trees_index(g.n()),
          n_labs(g.n()), n_pairs(g.n()), trav(g.n()), pll() {

        // sample size:
        int n = n_, log_n = 0;
        while (n >>= 1) ++log_n;
        desired_samples = const_log_n * log_n;
        G g_rev = g.reverse();

        // random perms:
        for (V v = 0; v < n_; ++v) {
            sample_fwd[v] = v;
            sample_bwd[v] = v;
        }
        for (V v = n_ - 1; v > 0; --v) {
            std::swap(sample_fwd[v], sample_fwd[rand() % (v+1)]);
            std::swap(sample_bwd[v], sample_bwd[rand() % (v+1)]);
        }

        // go:
        logging log("--hl");
        pll.init_index(n_);
        for (int i_hub = 0; i_hub < n_; ++i_hub) {
            //std::cerr << i_hub <<"    \n";
            update_samples(g, g_rev, weighted);
            V u = G::not_vertex;;
            while ( (! queue.empty()) && is_hub[u = queue.pop()] ); 
            if (u == G::not_vertex || is_hub[u]) {
                log.cerr() << i_hub <<" hubs\n";
                break;
            }
            is_hub[u] = true;
            // u as hub :
            pll.forward(g, u, trav, weighted, true);
            pll.backward(g_rev, u, trav, weighted, true);
            pll.incr_i_hub();
            if (log.progress()) {
                log.cerr() << i_hub <<" q="<< queue.size() * 100 / n_ <<"%"
                           <<" p/l="<< (n_pairs[u] / std::max(n_labs[u],1))
                           <<" l="<< n_labs[u]
                           <<" ts="<< avg_tree_size<<"\n  "
                           <<"smpl="<< (n_ - sample_fwd.size()) * 2
                           <<" all="<< sum_all_trees * 36 / 1000000 <<"m"
                           <<"\n  ";
                pll.progress(i_hub, true);
                std::cerr <<"\n";
            }
            update_trees(u);
        }
    }

    void print_stats(std::ostream &cout) { pll.print_stats(cout); }
    inline WL distance(int u, int v) { return pll.distance(u, v); }


    
private:
    
    void sample(const G &g, bool fwd, bool weighted) {
        assert((fwd ? sample_fwd.size() : sample_bwd.size()) > 0);
        V u = fwd ? sample_fwd.back() : sample_bwd.back();
        if (fwd) sample_fwd.pop_back(); else sample_bwd.pop_back();
        
        if (fwd) pll.forward(g, u, trav, weighted);
        else pll.backward(g, u, trav, weighted);

        n_all_trees += trav.nvis();
        int t = trees.size();
        trees.push_back(trav.digraph_to_sons());
        nodes.push_back(trav.digraph_nodes());
        //indexes.push_back(trav.digraph_index());
        auto const &nds = nodes.back();
        avg_tree_size = (4 * avg_tree_size + nds.size()) / 5;
        sum_all_trees += nds.size();
        for (int i = nds.size() - 1; i != -1; --i) {
            V u = nds[i].vtx;
            in_trees[u].push_back(t);
            in_trees_index[u].push_back(i);
            // skel: if (nds[i].dist >= alpha * (nds[i].dist_max - nds[i].dist))
            n_labs[u] += 1;
            n_pairs[u] += nds[i].size;
            assert(n_labs[u] > 0);
            queue.push(u);
        }
    }

    void update_samples(const G &g, const G &g_rev, bool weighted) {
        while (sample_fwd.size() > 0
               && n_all_trees < desired_samples * sample_fwd.size()) {
            assert(sample_fwd.size() == sample_bwd.size());
            sample(g, true, weighted);
            sample(g_rev, false, weighted);
            
            desired_samples = 1;
            int n = avg_tree_size;
            while (n >>= 1) desired_samples += const_log_n;
            //if (sample_fwd.size() < n_/2) break; // 1 at each time at most
        }
    }

    int dfs_del(int i, G &tre,
                 std::vector<typename trav_G::node> &nds) {
        int count = 0;
        if (nds[i].size > 0) {
            count = 1;
            V v = nds[i].vtx;
            //idx.erase(v);
            n_labs[v] -= 1;
            n_pairs[v] -= nds[i].size;
            assert(n_labs[v] >= 0);
            if (n_labs[v] == 0) queue.remove(v);
            else queue.push(v);
            nds[i].size = 0;
            //nds[i].dist_max = zero_weight;
            for (int j : tre[i]) { count += dfs_del(j, tre, nds); }
        }
        return count;
    }

    void update_trees(V x) { // x was taken as hub
        int i = 0;
        for (int t : in_trees[x]) {
            int ix = in_trees_index[x][i++];
            auto &nds = nodes[t];
            assert(nds[ix].vtx == x);
            //auto search = indexes[t].find(x);
            if (nds[ix].size > 0) {//search != indexes[t].end()) {
                //assert(ix == search->second);
                auto &tre = trees[t];
                int x_subtree_size = nds[ix].size;
                n_all_trees -= x_subtree_size;
                //WL x_dmax = nds[ix].dist_max;

                // Descendants are not possible hubs for the root
                int count = dfs_del(ix, tre, nds);
                assert(count == x_subtree_size);

                // Update treesize up to the root
                int ip = nds[ix].parent;
                while(ip != ix) {
                    nds[ip].size -= x_subtree_size;
                    V v = nds[ip].vtx;
                    n_pairs[v] -= x_subtree_size;
                    /* skeleton
                    if (nds[ip].dist_max == x_dmax) {
                        WL dm = zero_weight;
                        for (int j : tre[ip])
                            dm = std::max(dm, nds[j].dist_max);
                        nds[ip].dist_max = dm;
                    }
                    */
                    assert(n_labs[v] > 0);
                    queue.push(v);
                    ix = ip;
                    ip = nds[ix].parent;
                }
                }
        }
        in_trees[x] = {};
        in_trees_index[x] = {};
    }

};

#endif // HUB_LABELING_HH
