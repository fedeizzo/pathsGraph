#include <algorithm>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <math.h>
#include <queue>
#include <stack>
#include <unordered_map>
#include <vector>

using namespace std;

void bfs(vector<vector<int>> &graph, int &node, vector<int> &erdos, int &end) {
    queue<int> Q;
    Q.push(node);
    erdos.resize(0);
    for (int i = 0; i < graph.size(); i++) {
        erdos.push_back(-1);
    }
    erdos[node] = 0;

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        for (int value : graph[u]) {
            if (erdos[value] == -1) {
                erdos[value] = erdos[u] + 1;
                Q.push(value);
            }
            if (value == end) {
                return;
            }
        }
    }
}

vector<int> euler;
vector<int> depths;
vector<vector<int>> sparseMatrix;
vector<int> power2;
vector<int> logN;
unordered_map<int, int> nodesMapping;

/* void dfs(vector<vector<int>> graph, int node) { */
/*     stack<int> S; */
/*     S.push(node); */
/*     vector<bool> visited; */
/*     for (int i = 0; i < graph.size(); i++) { */
/*         visited.push_back(false); */
/*     } */
/*     while (!S.empty()) { */
/*         int u = S.top(); */
/*         S.pop(); */
/*         if (!visited[u]) { */
/*             /1* cout << "nodo: " << u << endl; *1/ */
/*             euler.push_back(u); */
/*             visited[u] = true; */
/*             for (int a : graph[u]) { */
/*                 euler.push_back(a); */
/*                 /1* cout << "   arco (" << u << ", " << a << ")" << endl; *1/
 */
/*                 S.push(a); */
/*             } */
/*         } */
/*     } */
/* } */

void tarjan(vector<vector<int>> &graph) {}

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

void dfs(vector<vector<int>> &graph, int node) {
    vector<bool> visited;
    visited.resize(graph.size(), false);
    dfsRec(graph, node, visited, 0);
}

void preprocess(int N) {
    power2.resize(17);
    power2[0] = 1;
    for (int i = 1; i < 18; ++i) {
        power2[i] = power2[i - 1] * 2;
    }

    int val = 1, ptr = 0;
    logN.resize(N);
    for (int i = 1; i < N; ++i) {
        logN[i] = ptr - 1;
        if (val == i) {
            val *= 2;
            logN[i] = ptr;
            ptr++;
        }
    }
}

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

int query(int start, int end) {
    int p = 31 - __builtin_clz(end - start);
    return min(sparseMatrix[p][start], sparseMatrix[p][end - (1 << p)]);
}

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
    int N, M, Q;
    vector<int> erdos;
    vector<int> parent;
    vector<pair<int, int>> request;
    in >> N >> M >> Q;
    graph.resize(N);
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

    dfs(graph, 0);
    sparseMatrixComputation();
    preprocess(graph.size());
    cout << "eulero e profondita" << endl;
    int ca = 0;
    for (int v : euler) {
        cout << ca++ << " ";
    }
    cout << endl;
    for (int v : euler) {
        cout << v << " ";
    }
    cout << endl;
    for (int d : depths) {
        cout << d << " ";
    }
    cout << endl;
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
