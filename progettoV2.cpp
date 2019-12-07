#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <vector>

using namespace std;

// Vector to store the euler walk
vector<int> euler;

// Vector to store the depths of 
// vertexes in the euler walk
vector<int> depths;

// This matrix will compute the access in O(1) to 
// find operation of the LCA (lowest common ancestor)
// of a pair of vertexes
vector<vector<int>> sparseMatrix;

// Hash map to map the vertex to the position in euler walk (index)
unordered_map<int, int> nodesMapping;

// Matrix of cliques where each row is aclique and 
// its columns are the nodes of the clique
vector<vector<int>> cliques;

// Value used to identify the cliques: 0 for node which 
// don't belong to any clique, 1 for the nodes in the 
// first clique, 2 for the nodes in the second clique and so on
vector<int> isInCliqueNumber;

int counterCliques = 1;

void printSet(set<int> s) {
    if (s.size() > 2) {
        vector<int> clique;
        for (auto x : s) {
            clique.push_back(x);
            isInCliqueNumber[x] = counterCliques;
        }
        cliques.push_back(clique);
        counterCliques++;
    }
}

/**
 * Utility function to compute the union of two sets
 */
set<int> setUnion(set<int> a, set<int> b) {
    set<int> c;
    set_union(a.begin(), a.end(), b.begin(), b.end(), inserter(c, c.end()));
    return c;
}

/**
 * Utility function to compute the difference of two sets
 */
set<int> setDifference(set<int> a, set<int> b) {
    set<int> c;
    set_difference(a.begin(), a.end(), b.begin(), b.end(),
                   inserter(c, c.end()));
    return c;
}

/**
 * Utility function to compute the intersection of two sets
 */
set<int> setIntersection(set<int> a, set<int> b) {
    set<int> c;
    set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                     inserter(c, c.end()));
    return c;
}

/**
 * Bron-Kerbosh algorithm to find cliques in an undirected graph
 */
void bronKerbosh(set<int> R, set<int> P, set<int> X,
                 vector<vector<int>> &graph) {
    if (P.empty() && X.empty()) {
        printSet(R);
    }
    set<int>::iterator v = P.begin();
    while (!P.empty() && v != P.end()) {
        set<int> singleton = {(*v)};
        set<int> neigh;
        for (auto node : graph[(*v)]) {
            neigh.insert(node);
        }
        bronKerbosh(setUnion(R, singleton), setIntersection(P, neigh),
                    setIntersection(X, neigh), graph);
        P = setDifference(P, singleton);
        X = setUnion(X, singleton);
        if (!P.empty()) {
            v = P.begin();
        }
    }
}

/**
 * Recursive part of dfs with pre-order visit; this will
 * populate euler walk and depths for nodes of the graph
 */
void dfsRec(vector<vector<int>> &graph, int &node, vector<bool> &visited,
            int depth) {
    visited[node] = true;
    euler.push_back(node);
    depths.push_back(depth);
    nodesMapping[node] = euler.size() - 1;
    for (int a : graph[node]) {
        if (!visited[a]) {
            dfsRec(graph, a, visited, depth + 1);
            euler.push_back(node);
            depths.push_back(depth);
            nodesMapping[node] = euler.size() - 1;
        }
    }
}

/**
 * This procedure just creates a vector to keep track of
 * which node gets visited and calls the recursion
 */
void dfs(vector<vector<int>> &graph, int node) {
    vector<bool> visited;
    visited.resize(graph.size(), false);
    dfsRec(graph, node, visited, 0);
}

/**
 * This procedure fills the matrix[log(euler.size()) + 1][euler.size()]
 * in complexity O(n*log(n)); this will allow to compute the
 * LCA (lowest common ancestor) of a pair of vertexes
 * in O(1). It takes O(n*log(n)) to store the matrix
 */
void sparseMatrixComputation() {
    int n = euler.size();
    int h = floor(log2(n));
    sparseMatrix.resize(log2(n) + 1);
    for (int i = 0; i < sparseMatrix.size(); ++i) {
        sparseMatrix[i].resize(n);
    }
    for (int i = 0; i < n; ++i) {
        sparseMatrix[0][i] = depths[i];
    }
    for (int i = 1; i <= h; i++) {
        for (int j = 0; j + (1 << i) <= n; j++) {
            sparseMatrix[i][j] = min(sparseMatrix[i - 1][j],
                                     sparseMatrix[i - 1][j + (1 << (i - 1))]);
        }
    }

    /* cout << "matrice" << endl; */
    /* cout << "-----------------------------------" << endl; */
    /* for (vector<int> vec : sparseMatrix) { */
    /*     for (int v : vec) { */
    /*         cout << v << " "; */
    /*     } */
    /*     cout << endl; */
    /* } */

    /* cout << "-----------------------------------" << endl; */
}

/**
 * The __builtin_clz function returns how many zeros occur
 * before the first occurence of 1 in the binary representation 
 * of a given number. Since integers in c++ are 32 bits values, 
 * we subtract from 31 from the number of zeros to get p, defined
 * as the position of the most significant bit set to 1. 
 * The function uses p to access the matrix and 
 * find the depth of the Lowest Common Ancestor
 */
int query(int start, int end) {
    int p = 31 - __builtin_clz(end - start);
    return min(sparseMatrix[p][start], sparseMatrix[p][end - (1 << p)]);
}

/**
 * This takes as input a pair containing the start and end nodes of a request.
 * It computes in O(1) the position in the euler walk and passes it to
 * query function, which returns the depths of the LCA and 
 * allow us to have the request done
 */
int checkRequest(pair<int, int> request) {
    int start, end;
    start = nodesMapping[request.first];
    end = nodesMapping[request.second];
    if (start > end) {
        swap(start, end);
    }
    /* cout << "start: " << start << endl; */
    /* cout << "end: " << end << endl; */
    /* cout << "query: " << query(start, end) << endl; */

    int returnValue = depths[start] + depths[end] - 2 * query(start, end + 1);

    return returnValue;
}

int main() {
    ifstream in("input.txt");
    vector<vector<int>> graph;
    set<int> R, P, X;
    int N, M, Q;
    vector<pair<int, int>> request;
    in >> N >> M >> Q;
    graph.resize(N);
    isInCliqueNumber.resize(N);
    int node1, node2;

    for (int i = 0; i < M; ++i) {
        in >> node1 >> node2;
        graph[node1].push_back(node2);
        graph[node2].push_back(node1);
    }

    int start, end;
    for (int i = 0; i < Q; ++i) {
        in >> start >> end;
        request.push_back(pair<int, int>(start, end));
    }

    for (int i = 0; i < graph.size(); i++) {
        P.insert(i);
    }

    dfs(graph, 0);
    sparseMatrixComputation();
    bronKerbosh(R, P, X, graph);
    /* cout << "---------------------------------------------------" << endl; */
    /* cout << "cricche: " << endl; */
    /* for (auto x : cliques) { */
    /*     for (auto y : x) { */
    /*         cout << y << " "; */
    /*     } */
    /*     cout << endl; */
    /* } */
    /* cout << "---------------------------------------------------" << endl; */
    /* cout << "vettori in cricca: " << endl; */
    /* for (auto x : isInCliqueNumber) { */

    /*     cout << x << " "; */
    /* } */
    /* cout << endl; */
    /* cout << "---------------------------------------------------" << endl; */
    /* cout << "eulero e profondita" << endl; */
    /* int ca = 0; */
    /* for (int v : euler) { */
    /*     cout << ca++ << " "; */
    /* } */
    /* cout << endl; */
    /* for (int v : euler) { */
    /*     cout << v << " "; */
    /* } */
    /* cout << endl; */
    /* for (int d : depths) { */
    /*     cout << d << " "; */
    /* } */
    /* cout << endl; */
    /* cout << "---------------------------------------------------" << endl; */
    /* cout << "mapping" << endl; */
    /* for (auto x : nodesMapping) { */
    /*     cout << x.first << " " << x.second << endl; */
    /* } */
    /* cout << "---------------------------------------------------" << endl; */
    /* for (int i : power2) { */
    /*     cout << i << " "; */
    /* } */
    /* cout << endl; */
    /* for (int i : logN) { */
    /*     cout << i << " "; */
    /* } */
    /* cout << endl; */
    /* cout << "---------------------------------------------------" << endl; */
    in.close();
    ofstream out("output.txt");
    for (pair<int, int> r : request) {
        out << checkRequest(r) << endl;
    }
    out.close();
    return 0;
}
