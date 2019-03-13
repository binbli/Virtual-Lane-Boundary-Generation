#ifndef UNION_FIND_H_
#define UNION_FIND_H_
#include<vector>
//taken from Competitive Programming 3 by Steven Halim
namespace LD {

    class UnionFind {
        // OOP style
    private:
        std::vector<long long int> p, rank;

    public:
        UnionFind(long long int N) {
            rank.assign(N, 0);
            p.assign(N, 0);
            for (long long int i = 0; i < N; i++) p[i] = i;
        }

        long long int findSet(long long int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }

        bool isSameSet(long long int i, long long int j) { return findSet(i) == findSet(j); }

        void unionSet(long long int i, long long int j) {
            if (!isSameSet(i, j)) {
                // if from different set
                long long int x = findSet(i), y = findSet(j);
                if (rank[x] > rank[y]) p[y] = x;
                    // rank keeps the tree short
                else {
                    p[x] = y;
                    if (rank[x] == rank[y]) rank[y]++;
                }
            }
        }
    };

}

#endif
