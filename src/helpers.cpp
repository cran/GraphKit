#include "helpers.h"

int edgeNum(int i, int j) {
    return i - j + (j * (j + 1))/2;
    
}

UnionFind::UnionFind (int n) {
    count = n;
    parent.reserve(n);
    size.reserve(n);
    
    for(int i = 0; i < n; i++) {
        parent[i] = i;
        size[i] = 1;
    }
}

int UnionFind::getCount() {
    return count;
}

int UnionFind::find(int p) {
    int root = p;
    
    while (root != parent[root])
        root = parent[root];
    while (p != root) {
        int newp = parent[p];
        parent[p] = root;
        p = newp;
    }
    
    return root;
}

bool UnionFind::connected(int p, int q) {
    return find(p) == find(q);
}

void UnionFind::connect(int p, int q) {
    int rootP =  find(p);
    int rootQ = find(q);
    if (rootP == rootQ) return;
    
    if (size[rootP] < size[rootQ]) {
        parent[rootP] = rootQ;
        size[rootQ] += size[rootP];
    }
    
    else {
        parent[rootQ] = rootP;
        size[rootP] += size[rootQ];
    }
    
    count--;
}

int main() {
    UnionFind uf = UnionFind(10);
    uf.connect(3,5);
    uf.connect(5,6);
}
