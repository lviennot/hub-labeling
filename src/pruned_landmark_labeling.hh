#ifndef PRUNED_LANDMARK_LABELING_HH
#define PRUNED_LANDMARK_LABELING_HH

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <xmmintrin.h> // prefetch

#include "edge.hh"
#include "traversal.hh"
#include "tree_sampling.hh"

template<typename G,
         typename WL = int64_t, // long weights
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class pruned_landmark_labeling {

public:
    typedef typename G::vertex V;

    struct hubinfo {
        V hub; // original ID
        WL dist;
        V next_hop;
        //int next_hub_index;
        hubinfo(V x, WL d, V nh) : hub(x), dist(d), next_hop(nh) {}
        hubinfo() : hub(-1), dist(INT64_MAX), next_hop(-1) {}
    };

    typedef typename edge::src_dst_wgt<V, hubinfo> edgeL;

private:
    typedef typename G::weight W;
    typedef typename edge::src_dst_wgt<V, W> edge;
    
    struct label_t {
        std::vector<V> out_v; // out hubs
        std::vector<WL> out_d;   // dist to them
        char pad1[16]; // 24+24+16=64
        std::vector<V> in_v;  // in hubs
        std::vector<WL> in_d;    // dist from them
        char pad2[16];
        std::vector<V> out_nh; // out next hop
        std::vector<V> in_nh;  // in next hop
        char pad3[16];
    }; // __attribute__((aligned(64)));

    static inline void * memalign (size_t bytes) {
        void *ptr;
        assert(0 == posix_memalign(&ptr, 64, bytes));
        return ptr;
    }
    
    int n_; // number of nodes
    std::vector<V> ranked_hubs;
    std::vector<label_t> index_;

    int i_hub; // next hub number
    int64_t sum_nvis, last_nvis, last_r; // stats for progress


public:
    pruned_landmark_labeling(const G &g,
                             std::vector<V> sel_from = {},
                             std::vector<V> sel_to = {},
                             int boost_sel = 1,
                             bool weighted = true) {
        logging hl_log("- hl -");
        
        init_index(g.n());
        tree_sampling<G,WL,max_weight,zero_weight> tsampl(n_, boost_sel);

        // create labels for all pairs u,v with u in sel_from, v in sel_to
        if (sel_from.size() == 0) { // put all nodes
            sel_from.reserve(n_);
            for (int v = 0; v < n_; ++v) sel_from.push_back(v);
        }
        if (sel_to.size() == 0) { // put all nodes
            sel_to.reserve(n_);
            for (int v = 0; v < n_; ++v) sel_to.push_back(v);
        }
        // random perms:
        for (V i = sel_from.size() - 1; i > 0; --i) {
            std::swap(sel_from[i], sel_from[rand() % (i+1)]);
        }
        for (V i = sel_to.size() - 1; i > 0; --i) {
            std::swap(sel_to[i], sel_to[rand() % (i+1)]);
        }
        // bool vects:
        std::vector<bool> is_sel_from(n_, false), is_sel_to(n_, false);
        for (V v : sel_from) { is_sel_from[v] = true; }
        for (V v : sel_to) { is_sel_to[v] = true; }
        // Complete:
        int n_from = sel_from.size(), n_to = sel_to.size();
        sel_from.reserve(n_);
        for (int v = 0; v < n_; ++v) {
            if (! is_sel_from[v]) sel_from.push_back(v);
        }
        sel_to.reserve(n_);
        for (int v = 0; v < n_; ++v) {
            if (! is_sel_to[v]) sel_to.push_back(v);
        }
        // random perm:
        for (V i = sel_from.size() - 1; i > n_from; --i) {
            std::swap(sel_from[i], sel_from[n_from + (rand() % (i+1-n_from))]);
        }
        for (V i = sel_to.size() - 1; i > n_to; --i) {
            std::swap(sel_to[i], sel_to[n_to + (rand() % (i+1-n_to))]);
        }
        
        // sample size:
        const int64_t const_log_n = 4; // 2->10 ==> max lab -= 20%, 2->1 => +20%
        int64_t n = n_, log_n = 0;
        while (n >>= 1) ++log_n;
        int64_t desired_samples = log_n * const_log_n * n_;
        int64_t desired_pairs = log_n * const_log_n * boost_sel;

        G g_rev = g.reverse();
        trav_G trav(n_);
        std::vector<bool> is_hub(n_, false);
 
        for (V i_from = 0, i_to = 0;
             i_from < sel_from.size() || i_to < sel_to.size()
                 || ! tsampl.queue_empty()  ; ) {
            // release memory:
            if (tsampl.size_deleted()
                  > std::max(tsampl.size_all(), desired_samples) ) {
                hl_log.cerr() << "clear tsampl "<< tsampl.size_deleted()
                              << " "<< tsampl.size_all() <<" "<< desired_samples
                              << "\n";
                tsampl.clear();
                i_from = 0;
                i_to = 0;
            }
            // more random sampling:
            while ((tsampl.size_all() < desired_samples
                    || tsampl.size_pairs() < desired_pairs)
                   && (i_from < sel_from.size() || i_to < sel_to.size())) {
                while (i_from < sel_from.size() && is_hub[sel_from[i_from]])
                    ++i_from;
                if (i_from < sel_from.size()) {
                    forward(g, sel_from[i_from++], trav, weighted, false);
                    for (int i = trav.nvis() - 1 ; i > 0 ; --i) {
                        assert ( ! is_hub[trav.visit(i)] );
                    }
                    tsampl.add_tree(trav, &is_sel_to);
                }
                while (i_to < sel_to.size() && is_hub[sel_to[i_to]])
                    ++i_to;
                if (i_to < sel_to.size()) {
                    backward(g_rev, sel_to[i_to++], trav, weighted, false);
                    for (int i = trav.nvis() - 1 ; i > 0 ; --i) {
                        assert ( ! is_hub[trav.visit(i)] );
                    }
                    tsampl.add_tree(trav, &is_sel_from);
                }
            }
            // best covering node:
            if (tsampl.queue_empty()) {
                if (i_from < sel_from.size() || i_to < sel_to.size()) {
                    // need to resample
                    hl_log.cerr() << "clear tsampl, no sel pair covered\n";
                    tsampl.clear();
                    i_from = 0;
                    i_to = 0;
                }
                continue;
            }
            V x = tsampl.queue_pop_best_cover();
            assert( ! is_hub[x]);
            is_hub[x] = true;
            forward(g, x, trav, weighted, true);
            backward(g_rev, x, trav, weighted, true);
            // update sampled trees:
            tsampl.remove_subtrees(x);
            assert(ranked_hubs.size() == i_hub);
            ranked_hubs.push_back(x);
            ++i_hub;
            //progress(i_hub);
            if (hl_log.progress()) {
                hl_log.cerr() <<"hub "<< i_hub 
                              <<" : avg_nvis="<< (sum_nvis / (2*(i_hub+1)))
                              <<" lst_nvis="<< (last_nvis / 2 /(i_hub - last_r))
                    //<<" n=" << n_ <<" "
                              <<" avg_hs=" << (sum_nvis / (2*n_))
                              <<" (" << (sum_nvis * 12 / 1000000) <<"m)"
                              <<"\n";
                last_r = i_hub - 1;
                last_nvis = 0;
            }
        }
        //std::cerr << "\n";
        print_stats(std::cerr, is_sel_from, is_sel_to);
    }

    V hub_ID(int i) { return ranked_hubs[i]; }

    inline WL distance(int u, int v) {
        const label_t &lab_u = index_[u];
        const label_t &lab_v = index_[v];

        WL d = max_weight;

        //__mm_prefetch(&(lab_u.out_v[0], _MM_HINT_T0));
        //__mm_prefetch(&(lab_v.in_v[0], _MM_HINT_T0));

        for (int iu = 0, iv = 0 ; ; ) {
            V xu = lab_u.out_v[iu], xv = lab_v.in_v[iv];
            if (xu == xv) { // common hub
                if (xu == n_) break; // sentinel
                WL d_uxv = lab_u.out_d[iu] + lab_v.in_d[iv];
                if (d_uxv < d) d = d_uxv;
                ++iu;
                ++iv;
            } else {
                if (xu < xv) ++iu;
                else ++iv;
            }
        }

        return d;
    }

    inline std::pair<int,int> common_hub(int u, int v) {
        const label_t &lab_u = index_[u];
        const label_t &lab_v = index_[v];

        WL d = max_weight;

        //__mm_prefetch(&(lab_u.out_v[0], _MM_HINT_T0));
        //__mm_prefetch(&(lab_v.in_v[0], _MM_HINT_T0));

        int iu_opt = -1, iv_opt = -1;
        
        for (int iu = 0, iv = 0 ; ; ) {
            V xu = lab_u.out_v[iu], xv = lab_v.in_v[iv];
            if (xu == xv) { // common hub
                if (xu == n_) { break; } // sentinel
                WL d_uxv = lab_u.out_d[iu] + lab_v.in_d[iv];
                if (d_uxv < d) { d = d_uxv; iu_opt = iu; iv_opt = iv; }
                ++iu;
                ++iv;
            } else {
                if (xu < xv) ++iu;
                else ++iv;
            }
        }

        return std::pair<int,int>(iu_opt, iv_opt);
    }

    std::vector<edgeL> out_hub_edges(const std::vector<bool> &is_sel_from = {},
                                    const std::vector<bool> &is_sel_to = {}) {
        std::vector<edgeL> edg;
        int64_t fwd = 0;
        bool fwd_sel = is_sel_from.size() == 0;
        for (int u = 0; u < n_; ++u) {
            if (fwd_sel || is_sel_from[u]) {
                fwd += index_[u].out_v.size();
            }
        }
        edg.reserve(fwd);
        for (int u = 0; u < n_; ++u) {
            if (fwd_sel || is_sel_from[u]) {
                const label_t &lab_u = index_[u];
                for (int i = 0; i < lab_u.out_v.size(); ++i) {
                    if (lab_u.out_v[i] == n_) break; // sentinel
                    V x = ranked_hubs[lab_u.out_v[i]];
                    hubinfo hi(x, lab_u.out_d[i], lab_u.out_nh[i]);
                    edg.push_back(edgeL(u, lab_u.out_v[i], hi));
                }
            }
        }
        return edg;
    }

    std::vector<edgeL> in_hub_edges(const std::vector<bool> &is_sel_from = {},
                                   const std::vector<bool> &is_sel_to = {}) {
        std::vector<edgeL> edg;
        int64_t bwd = 0;
        bool bwd_sel = is_sel_to.size() == 0;
        for (int u = 0; u < n_; ++u) {
            if (bwd_sel || is_sel_to[u]) {
                bwd += index_[u].in_v.size();
            }
        }
        edg.reserve(bwd);
        for (int u = 0; u < n_; ++u) {
            if (bwd_sel || is_sel_to[u]) {
                const label_t &lab_u = index_[u];
                for (int i = 0; i < lab_u.in_v.size(); ++i) {
                    if (lab_u.in_v[i] == n_) break; // sentinel
                    V x = ranked_hubs[lab_u.in_v[i]];
                    hubinfo hi(x, lab_u.in_d[i], lab_u.in_nh[i]);
                    edg.push_back(edgeL(lab_u.in_v[i], u, hi));
                }
            }
        }
        return edg;
    }

    void print_stats(std::ostream &cout,
                     const std::vector<bool> &is_sel_from = {},
                     const std::vector<bool> &is_sel_to = {}) {
        int64_t fwd = 0, bwd = 0, fmax = 0, bmax = 0;
        int n_fwd = 0, n_bwd = 0;
        bool fwd_sel = is_sel_from.size() == 0, bwd_sel = is_sel_to.size() == 0;
        for (int u = 0; u < n_; ++u) {
            if (fwd_sel || is_sel_from[u]) {
                if (index_[u].out_v.size() > fmax)
                    fmax = index_[u].out_v.size();
                fwd += index_[u].out_v.size();
                ++n_fwd;
            }
        }
        for (int u = 0; u < n_; ++u) {
            if (bwd_sel || is_sel_to[u]) {
                if (index_[u].in_v.size() > bmax)
                    bmax = index_[u].in_v.size();
                bwd += index_[u].in_v.size();
                ++n_bwd;
            }
        }
        cout << n_fwd <<" fwd "<< n_bwd  <<" bwd / "<< n_ << " labels\n";
        if (n_ > 0) cout << "  avg fwd = " << (n_fwd == 0 ? 0 : fwd / n_fwd)
                         << "  avg bwd = " << (n_bwd == 0 ? 0 : bwd / n_bwd)
                         << "  max fwd = " << fmax
                         << "  max bwd = " << bmax
                         << "\n";
    }

    pruned_landmark_labeling(const G &g, int  dummy,
                             const std::vector<V> &rank_order,
                             bool wgted = true) {
        construct_index(g, rank_order, wgted);
    }

    // restrict labels to cover only pair from [sel_from] to [sel_to].
    pruned_landmark_labeling(pruned_landmark_labeling
                             <G, WL, max_weight, zero_weight> pll,
                             const std::vector<V> &sel_from,
                             const std::vector<V> &sel_to) {
        n_ = pll.n_;
        index_ = std::vector<label_t>(n_);
        
        std::vector<bool> keep(n_, false);
        for (V u : sel_from) {
            int lab_len = pll.index_[u].out_v.size(), new_len = 1;
            for (V v : sel_to) {
                std::pair<int,int> h = pll.common_hub(u,v);
                if (h.first >= 0) { keep[h.first] = true; ++new_len; }
            }
            label_t &lab_u = index_[u], &lab_pll = pll.index_[u];
            lab_u.out_v.reserve(new_len);
            lab_u.out_d.reserve(new_len);
            for (int iu = 0; iu < lab_len; ++iu) {
                if (keep[iu]) {
                    keep[iu] = false;
                    lab_u.out_v.push_back(lab_pll.out_v[iu]);
                    lab_u.out_d.push_back(lab_pll.out_d[iu]);
                }
            }
            lab_u.out_v.push_back(n_);
            lab_u.out_d.push_back(max_weight);
        }
        
        for (V u : sel_to) {
            int lab_len = pll.index_[u].in_v.size(), new_len = 1;
            for (V v : sel_from) {
                std::pair<int,int> h = pll.common_hub(v, u);
                if (h.second >= 0) { keep[h.second] = true; ++new_len; }
            }
            label_t &lab_u = index_[u], &lab_pll = pll.index_[u];
            lab_u.in_v.reserve(new_len);
            lab_u.in_d.reserve(new_len);
            for (int iu = 0; iu < lab_len; ++iu) {
                if (keep[iu]) {
                    keep[iu] = false;
                    lab_u.in_v.push_back(lab_pll.in_v[iu]);
                    lab_u.in_d.push_back(lab_pll.in_d[iu]);
                }
            }
            lab_u.in_v.push_back(n_);
            lab_u.in_d.push_back(max_weight);
        }
    }

    
    void restrict(const std::vector<V> &sel_from,
                  const std::vector<V> &sel_to) {
        std::vector<bool> keep(n_, false);

        for (V u : sel_from) {
            int lab_len = index_[u].out_v.size(), new_len = 1;
            for (V v : sel_to) {
                std::pair<int,int> h = common_hub(u,v);
                if (h.first >= 0) { keep[h.first] = true; ++new_len; }
            }
            label_t &lab_u = index_[u];
            std::vector<V> out_v; out_v.reserve(new_len);
            std::vector<WL> out_d; out_d.reserve(new_len);
            for (int iu = 0; iu < lab_len; ++iu) {
                if (keep[iu]) {
                    keep[iu] = false;
                    out_v.push_back(lab_u.out_v[iu]);
                    out_d.push_back(lab_u.out_d[iu]);
                }
            }
            out_v.push_back(n_);
            out_d.push_back(max_weight);
            std::swap(lab_u.out_v, out_v);
            std::swap(lab_u.out_d, out_d);
        }
        
        for (V u : sel_to) {
            int lab_len = index_[u].in_v.size(), new_len = 1;
            for (V v : sel_from) {
                std::pair<int,int> h = common_hub(v, u);
                if (h.second >= 0) { keep[h.second] = true; ++new_len; }
            }
            label_t &lab_u = index_[u];
            std::vector<V> in_v; in_v.reserve(new_len);
            std::vector<WL> in_d; in_d.reserve(new_len);
            for (int iu = 0; iu < lab_len; ++iu) {
                if (keep[iu]) {
                    keep[iu] = false;
                    in_v.push_back(lab_u.in_v[iu]);
                    in_d.push_back(lab_u.in_d[iu]);
                }
            }
            in_v.push_back(n_);
            in_d.push_back(max_weight);
            std::swap(lab_u.in_v, in_v);
            std::swap(lab_u.in_d, in_d);
        }
    }    


    void init_index(int n) {
        // Number of vertices
        n_ = n;
        i_hub = 0;
        sum_nvis = 0; last_nvis = 0; last_r = -1;

        // Allocate index
        index_ = std::vector<label_t>(n_);
        for (int u = 0; u < n_; ++u) {
            label_t &lab_u = index_[u];
            lab_u.in_v.push_back(n_);
            lab_u.in_nh.push_back(n_);
            lab_u.in_d.push_back(max_weight);
            lab_u.out_v.push_back(n_);
            lab_u.out_nh.push_back(n_);
            lab_u.out_d.push_back(max_weight);
        }
    }

    typedef traversal<G, WL, max_weight, zero_weight> trav_G;

    void forward(const G &g, V u, trav_G &trav,
                 bool wgted, bool add_hub=false) {
        auto filter = [this, u](V v, WL dv, V prt, WL dprt){
            return dv < distance(u, v);
        };
        trav.clear();
        if (wgted) trav.dijkstra(g, u, filter);
        else trav.bfs(g, u, filter);
        if (add_hub) {
            sum_nvis += trav.nvis(); last_nvis += trav.nvis();
            for (int i = trav.nvis() - 1; i != -1; --i) {
                int v = trav.visit(i);
                label_t &lab_v = index_[v];
                lab_v.in_v.back() = i_hub;
                lab_v.in_nh.back() = trav.parent(v);
                lab_v.in_d.back() = trav.dist(v);
                lab_v.in_v.push_back(n_);
                lab_v.in_nh.push_back(n_);
                lab_v.in_d.push_back(max_weight);
            }
        }
    }
    
    void backward(const G &g_rev, V u, trav_G &trav,
                  bool wgted, bool add_hub=false) {
        auto filter = [this, u](V v, WL dv, V prt, WL dprt){
            return dv < distance(v, u);
        };
        trav.clear();
        if (wgted) trav.dijkstra(g_rev, u, filter);
        else trav.bfs(g_rev, u, filter);
        if (add_hub) {
            sum_nvis += trav.nvis(); last_nvis += trav.nvis();
            for (int i = trav.nvis() - 1; i != -1; --i) {
                int v = trav.visit(i);
                label_t &lab_v = index_[v];
                lab_v.out_v.back() = i_hub;
                lab_v.out_nh.back() = trav.parent(v);
                lab_v.out_d.back() = trav.dist(v);
                lab_v.out_v.push_back(n_);
                lab_v.out_nh.push_back(n_);
                lab_v.out_d.push_back(max_weight);
            }
        }
    }

    void progress(int r, bool force = false) {
        if ((force && r > last_r)
            || r <= 10 || (r <= 100 && r % 10 == 0) || r % 100 == 0) {
            std::cerr << r  << " avg_nvis="<< (sum_nvis / (2*(r+1)))
                      <<" lst_nvis="<< (last_nvis / 2 / (r - last_r))
                      << " n=" << n_ <<" "
                      << "avg_hs=" << (sum_nvis / (2*n_))
                      << " (" << (sum_nvis * 12 / 1000000) <<"m)"
                      <<"         \r";
            std::cerr.flush();
            last_r = r;
            last_nvis = 0;
        }
    }

    void incr_i_hub() { ++i_hub; }
    
    void construct_index(const G &g, std::vector<V> rank_order, bool wgted) {
        init_index(g.n());
        assert(n_ == rank_order.size());
        G g_rev = g.reverse();
        trav_G trav(n_);
        for ( ; i_hub < n_; ++i_hub) {
            V u = rank_order[i_hub];
            forward(g, u, trav, wgted, true);
            backward(g_rev, u, trav, wgted, true);
            assert(ranked_hubs.size() == i_hub);
            ranked_hubs.push_back(u);
            progress(i_hub);
        }
        std::cerr << "\n";
    }
    
};

#endif //PRUNED_LANDMARK_LABELING_HH
