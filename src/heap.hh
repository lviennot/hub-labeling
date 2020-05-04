#ifndef HEAP_HH
#define HEAP_HH

#include <assert.h>
#include <functional>
#include <vector>
#include <iostream>

/** Priority queue for ints in 0..n-1 through a heap. 
 *  Smallest element according to [cmp_less] is popped first. */


static const int not_pos = -1;

class heap {
private:
    std::function<bool(const int &, const int &)> cmp_less;
    std::vector<int> elts, pos;
    
public:
    heap(std::function<bool(const int &, const int &)> cmp, int n = 0)
        : cmp_less(cmp), elts(), pos(n, not_pos) {
       elts.reserve(n);
   }

    inline bool empty() { return elts.size() == 0; }
    inline int size() { return elts.size(); }
    inline int top() { assert(elts.size() > 0); return elts[0]; }
    void clear() {
        elts.clear();
        for (int e = 0; e < pos.size(); ++e) pos[e] = not_pos;
    }
    void set_compare(std::function<bool(const int &, const int &)> cmp) {
        cmp_less = cmp;
    }
    
    /** Push or update an element: if the [e]th element was not in the queue,
     *  add it, otherwise consider its order position according to [cmp_less]
     *  has been updated. */ 
    void push(int e) {
        if (pos[e] == not_pos) {
            pos[e] = elts.size();
            elts.push_back(e);
        }
        if ( ! move_up(pos[e])) {
            move_down(pos[e]);
        }
    }

    /** Pop the minimum element according to cmp_less. */
    int pop() {
        assert(elts.size() > 0);
        int e = elts[0];
        assert(pos[e] == 0);
        remove(e);
        return e;
    }

    void remove(int e) {
        int i = pos[e];
        if (i != not_pos) {
            pos[e] = not_pos;
            int j = elts.size() - 1; // index of elts.back()
            if (j != i) {
                pos[elts[j]] = i;
                std::swap(elts[i], elts[j]);
                elts.pop_back();
                move_down(i);
            } else {
                elts.pop_back();
            }
        }
    }

    void print(std::ostream &cout) {
        for(int i : elts) cout << i <<" ";
        cout << "\n";
    }

private:
    bool move_up(int i) {
        assert(i >= 0 && i < elts.size());
        int p = (i+1)/2 - 1; // parent
        int ui, up;
        bool goes_up = false;
        while (p >= 0 && cmp_less(elts[i], elts[p])) {
            goes_up = true;
            std::swap(elts[i], elts[p]);
            std::swap(pos[elts[i]], pos[elts[p]]);
            i = p;
            p = (i+1)/2 - 1;
        }
        return goes_up;
    }

    void move_down(int i) {
        int size = elts.size();
        assert(i >= 0 && i < size);
        while (true) {
            int c1 = 2*i + 1, c2 = 2*i + 2; // children
            if (c1 >= size) { break; }
            int cmin = c1;
            if (c2 < size) assert(pos[elts[c2]] >= 0);
            if (c2 < size && cmp_less(elts[c2], elts[c1])) { cmin = c2; }
            if (cmp_less(elts[cmin], elts[i])) {
                std::swap(elts[i], elts[cmin]);
                std::swap(pos[elts[i]], pos[elts[cmin]]);
                i = cmin;
            } else { break; }
        }
    }

};

namespace unit {
    
    void heap_test(int n, int rnd = 1000) {
        std::cerr <<"heap_test:";
        std::vector<int> key(n);
        for (int i = 0; i < n; ++i) { key[i] = rand() % rnd; }
        auto cmp = [&key](int i, int j) -> bool { return key[i] < key[j]; };
        auto get = [&key](int i) -> int { return key[i]; };
        heap h(cmp, n);
        for (int i = 0; i < n; ++i) {
            if (i % 5 == 3) {
                key[i] = rnd; // greater than others
                h.push(i); // at back of elts
                h.remove(i); // test of remove back
            } else {
                h.push(i);
            }
            for (int k = 0; k < 3 && i > 0; ++k) {
                int j = rand() % i;
                key[j] = rand() % rnd;
                h.push(j); // update
                h.push(j); // update
            }
            if (i % 5 == 4) h.remove(rand() % i);
            //h.print(std::cerr, get);
        }
        int last = INT_MIN;
        while ( ! h.empty()) {
            int i = h.pop();
            std::cerr /* << i */ <<" "<< key[i];
            assert(last <= key[i]);
            last = key[i];
        }
        std::cerr <<"\n";
    }

}


#endif // HEAP_HH
