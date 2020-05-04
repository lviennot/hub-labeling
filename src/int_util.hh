#ifndef INT_UTIL_HH
#define INT_UTIL_HH


template<typename I>
class int_iterator {
    I i;
    const bool incr;
public:
    int_iterator(I i, bool incr=true) : i(i), incr(incr) {}
    I operator*() const { return i; }
    int_iterator &operator++() { if(incr) ++i; else --i; return *this; }
    int_iterator &operator--() { assert( ! incr); --i; return *this; }
    bool operator!=(const int_iterator& o) { return i != o.i; }
};


template<typename I>
class range {
    const I from, to;
public:
    range(I from, I to) : from(from), to(to) {} 
    int_iterator<I> begin() const { return int_iterator<I>(from); }
    int_iterator<I> end() const { return int_iterator<I>(to); }
};

template<typename I>
class range_rev {
    const I from, to;
public:
    range_rev(I from, I to) : from(from), to(to) {}
    int_iterator<I> begin() const { return int_iterator<I>(from-1, false); }
    int_iterator<I> end() const { return int_iterator<I>(to-1, false); }
};

typedef range<int> irange;
typedef range_rev<int> irange_rev;

    
#endif // INT_UTIL_HH
