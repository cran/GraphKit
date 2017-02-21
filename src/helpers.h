#ifndef HELPERS_H
#define HELPERS_H

#include <vector>
#include <iostream>
using namespace std;

int edgeNum(int i, int j);

class UnionFind {
    
vector<int> parent;
vector<int> size;
int count;

public:
    int find(int p);
    int getCount();
    bool connected(int p, int q);
    void connect(int p, int q);
    UnionFind(int n);
};

#endif
