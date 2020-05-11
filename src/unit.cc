#include "heap.hh"
void heap_test(int, int);
#include "mgraph.hh"
void mgraph_test(int, int);

int main (int argc, char **argv) {
    for (int n : {0, 1, 2, 3, 5, 20, 60})
        for (int r : {2, 5, 1000}) heap_test(n, r);
    for (int n : {1, 5, 10, 20, 60})
        for (int d : {1, 2, 3, 5}) mgraph_test(n, d);
}
