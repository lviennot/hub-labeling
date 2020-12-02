#include <iostream>

#include "heap.hh"
#include "unit.hh"

void heap_test(int n, int rnd = 1000) {
    std::cout <<"heap_test:";
    std::vector<int> key(n);
    for (int i = 0; i < n; ++i) { key[i] = rand() % rnd; }
    auto cmp = [&key](int i, int j) -> bool { return key[i] < key[j]; };
    //auto get = [&key](int i) -> int { return key[i]; };
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
        //h.print(std::cout, get);
    }
    int last = INT_MIN;
    while ( ! h.empty()) {
        int i = h.pop();
        std::cout /* << i */ <<" "<< key[i];
        CHECK(last <= key[i]);
        last = key[i];
    }
    std::cout <<"\n";
}

