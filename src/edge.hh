#ifndef EDGE_HH
#define EDGE_HH

#include <assert.h>

namespace edge { // Various edge information.

    template<typename V> // vertex number type
    struct dst_only {
        V dst;
        dst_only(V v) : dst(v) {}
        dst_only() {}
        operator V() const { return dst; }
    };

    template<typename V> // vertex number type
    struct src_dst : public dst_only<V> {
        V src;
        src_dst(V u, V v) : dst_only<V>::dst_only(v), src(u) {}
        src_dst() {}
    };


    template<typename V, // vertex number type
             typename W> // weight type
    struct dst_wgt : public dst_only<V> {
        W wgt;
        dst_wgt(V v, W w) : dst_only<V>::dst_only(v), wgt(w) {}
        dst_wgt() {}
    };

    template<typename V, // vertex number type
             typename W> // weight type
    struct src_dst_wgt : public src_dst<V> {
        W wgt;
        src_dst_wgt(V u, V v, W w) : src_dst<V>::src_dst(u,v), wgt(w) {}
        src_dst_wgt(const src_dst<V> &e) : src_dst<V>::src_dst(e), wgt(1) {}
        src_dst_wgt() {}
    };

    
} // edge

#endif // EDGE_HH
