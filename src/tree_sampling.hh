#ifndef TREE_SAMPLING_HH
#define TREE_SAMPLING_HH

#include <assert.h>
#include <vector>
#include <queue>          // std::priority_queue
#include <set>

#include "heap.hh"
#include "traversal.hh"
#include "logging.hh"

/** Compute hub labels (alias a two-hop labeling) using random sampling and
 *  pruned landmark labeling. */

template<typename G,
         typename WL = int64_t, // long weigth (not necessarily T::long_weigth)
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class tree_sampling {
private:
    typedef typename G::vertex V;
    typedef typename G::weight W;
    typedef traversal<G, WL, max_weight, zero_weight> trav_G;
    typedef std::pair<int,int> tree_index;
    
    std::vector<G> trees;
    std::vector<std::vector<V> > trees_vtx;
    std::vector<const std::vector<bool>* > trees_sel;
    std::vector<std::vector<tree_index> > in_trees;
    std::vector<int> n_subtree;
    std::vector<int> n_labs; // how many labels of random sample would have i
    std::vector<int64_t> n_pairs; // how many pairs from/to random sample would
    int64_t n_all_trees, n_deleted, n_all_pairs;
    int boost_selection;

    heap queue;
    inline bool covers_more (V u, V v) {
        //return n_pairs[u] > n_pairs[v];
        assert(n_labs[u] != 0 && n_labs[v] != 0);
        return n_pairs[u] * n_labs[v] > n_pairs[v] * n_labs[u];
    }
    
public:
    tree_sampling(int n, int boost_sel = 1) : trees(), trees_vtx(), trees_sel(),
                           in_trees(n), n_subtree(n), n_labs(n), n_pairs(n),
                           n_all_trees(0), n_deleted(0), n_all_pairs(0),
                           boost_selection(boost_sel),
             queue([this](int u, int v){ return covers_more(u, v); }, n) {
    }

    void clear() {
        trees.clear();
        trees_vtx.clear();
        trees_sel.clear();
        int n = in_trees.size();
        for (int v = 0; v < n; ++v) {
            in_trees[v] = {};
            n_labs[v] = 0;
            n_pairs[v] = 0;
        }
        n_all_trees = 0;
        n_deleted = 0;
        n_all_pairs = 0;
        queue.clear();
    }

    int64_t size_all() { return n_all_trees; }
    int64_t size_deleted() { return n_deleted; }
    int64_t size_pairs() { return n_all_pairs; }
    bool queue_empty() { return queue.empty(); }
    V queue_pop_best_cover() { return queue.pop(); }

    void add_tree(const trav_G &trav, const std::vector<bool> *sel_ptr) {
        const std::vector<bool> &sel = *sel_ptr;
        int t = trees.size();
        std::vector<V> vtx(trav.nvis());
        int tn = trav.nvis();
        //bool reach_sel = false;
        for (int i = 0; i < tn; ++i) {
            V v = trav.visit(i);
            vtx[i] = v ;
            n_subtree[v] = 0;
            //if (sel[v]) reach_sel = true;
        }
        //if ( ! reach_sel) { return; } // no selected node to cover!
        trees.push_back(trav.graph());
        trees_vtx.push_back(vtx);
        assert(trees.size() == trees_vtx.size());
        trees_sel.push_back(sel_ptr);
        assert(trees.size() == trees_sel.size());
        n_all_trees += tn;
        for (int i = tn - 1; i != -1; --i) {
            V u = trav.visit(i);
            in_trees[u].push_back(tree_index(t, i));
            // skel: if (nds[i].dist >= alpha * (nds[i].dist_max - nds[i].dist))
            if (sel[u]) n_subtree[u] += boost_selection; else n_subtree[u] += 1;
            V p = trav.parent(u);
            if (p != u) n_subtree[p] += n_subtree[u];
            if (n_subtree[u] > 0) {
                n_labs[u] += 1;
                n_pairs[u] += n_subtree[u];
                n_all_pairs += n_subtree[u];
                //queue_.push(hub_stat(u, n_pairs[u], n_labs[u]));
                queue.push(u);
            }
        }
    }

    void add_tree_prune(const trav_G &trav, const std::vector<bool> *sel_ptr) {
        const std::vector<bool> &sel = *sel_ptr;
        int t = trees.size();
        std::vector<V> vtx(trav.nvis());
        int tn = trav.nvis();
        bool reach_sel = false;
        for (int i = 0; i < tn; ++i) {
            V v = trav.visit(i);
            vtx[i] = v ;
            n_subtree[v] = 0;
            if (sel[v]) reach_sel = true;
        }
        if ( ! reach_sel) { return; } // no selected node to cover!

        for (int i = tn - 1; i != -1; --i) {
            V u = trav.visit(i);
            // skel: if (nds[i].dist >= alpha * (nds[i].dist_max - nds[i].dist))
            if (sel[u]) n_subtree[u] += boost_selection;
            V p = trav.parent(u);
            if (p != u) n_subtree[p] += n_subtree[u];
        }

        G tr = trav.graph();
        auto sub = tr.subgraph([this, &vtx](int i){
                return n_subtree[vtx[i]] > 0;
            }, vtx);
        tr = sub.first;
        vtx = sub.second;
        
        trees.push_back(tr);
        trees_vtx.push_back(vtx);
        assert(trees.size() == trees_vtx.size());
        trees_sel.push_back(sel_ptr);
        assert(trees.size() == trees_sel.size());
        tn = tr.n();
        n_all_trees += tn;
        for (int i = tn - 1; i != -1; --i) {
            V u = sub.second[i];
            in_trees[u].push_back(tree_index(t, i));
            // skel: if (nds[i].dist >= alpha * (nds[i].dist_max - nds[i].dist))
            assert(n_subtree[u] > 0);
            n_labs[u] += 1;
            n_pairs[u] += n_subtree[u];
            n_all_pairs += n_subtree[u];
            //queue_.push(hub_stat(u, n_pairs[u], n_labs[u]));
            queue.push(u);
        }
    }

    int dfs_del(G &tree, std::vector<V> &vtx, const std::vector<bool> &sel,
                int i) {
        int count = 0;
        if (vtx[i] != G::not_vertex) {
            V v = vtx[i];
            vtx[i] = G::not_vertex; // mark as pruned
            if(sel[v]) count = boost_selection; else count = 1;
            for (int j : tree[i]) {
                if (j > i) // a son, not parent
                    count += dfs_del(tree, vtx, sel, j);
            }
            n_labs[v] -= 1;
            n_pairs[v] -= count;
            n_all_pairs -= count;
            assert(n_labs[v] >= 0 && n_pairs[v] >= 0);
            if (n_labs[v] == 0 || n_pairs[v] == 0) queue.remove(v);
            else queue.push(v);
        }
        return count;
    }
            
    int dfs_count(G &tree, std::vector<V> &vtx, int i) {
        int count = 0;
        if (vtx[i] != G::not_vertex) {
            //V v = vtx[i];
            count = 1;
            for (int j : tree[i]) {
                if (j > i) // a son, not parent
                    count += dfs_count(tree, vtx, j);
            }
        }
        return count;
    }
    void remove_subtrees(V x) {
        //int nd = 0;
        for (auto ti : in_trees[x]) {
            G &tre = trees[ti.first];
            std::vector<V> &vtx =  trees_vtx[ti.first];
            const std::vector<bool> &sel = *(trees_sel[ti.first]);
            int ix = ti.second;
            if (vtx[ix] == G::not_vertex) continue; // already removed from tre
            //++nd;
            int x_subtree_size = dfs_count(tre, vtx, ix);
            int x_subtree_score = dfs_del(tre, vtx, sel, ix);
            n_all_trees -= x_subtree_size;
            n_deleted += x_subtree_size;
            // Update treesize up to the root
            int ip = ix;
            while (true) {
                for (int j : tre[ix]) {
                    if (j < ix) { ip = j; break; }
                }
                if (ip >= ix) break; // at root
                V v = vtx[ip];
                assert(v != G::not_vertex);
                n_pairs[v] -= x_subtree_score;
                n_all_pairs -= x_subtree_score;
                assert(n_labs[v] > 0 && n_pairs[v] >= 0);
                if (n_labs[v] == 0 || n_pairs[v] == 0) queue.remove(v);
                else queue.push(v);
                ix = ip;
            }
        }
        assert(n_labs[x] == 0 && n_pairs[x] == 0);
        in_trees[x] = {};
    }
};

#endif // TREE_SAMPLING_HH
