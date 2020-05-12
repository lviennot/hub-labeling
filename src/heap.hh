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
    //-std=c++17: inline static const int not_pos = -1; 
    
    heap(std::function<bool(const int &, const int &)> cmp, int n = 0)
        : cmp_less(cmp) {
       elts.reserve(n);
       pos.assign(n, not_pos);
   }

    inline bool empty() { return elts.size() == 0; }
    inline size_t size() { return elts.size(); }
    inline int top() { assert(elts.size() > 0); return elts[0]; }
    void clear() {
        elts.clear();
        for (size_t e = 0; e < pos.size(); ++e) pos[e] = not_pos;
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
    bool move_up(size_t i) {
        if (i == 0) return false;
        assert(i > 0 && i < elts.size());
        size_t p = (i+1)/2 - 1; // parent
        bool goes_up = false;
        while (cmp_less(elts[i], elts[p])) {
            goes_up = true;
            std::swap(elts[i], elts[p]);
            std::swap(pos[elts[i]], pos[elts[p]]);
            if (p == 0) break;
            i = p;
            p = (i+1)/2 - 1;
        }
        return goes_up;
    }

    void move_down(size_t i) {
        size_t size = elts.size();
        assert(i >= 0 && i < size);
        while (true) {
            size_t c1 = 2*i + 1, c2 = 2*i + 2; // children
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



#endif // HEAP_HH
