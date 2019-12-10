#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

// ------------------- SOURCES ------------------
/**
 * (1) BFS --> Montresor's slides
 * (2) DFS --> Montresor's slides
 * (3) hash_pair -->
 * https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
 * (4) sparse matrix --> https://brilliant.org/wiki/sparse-table/
 */

/**
 * Structure for allow the use of a pair inside an unordered_map (used to store
 * siblings) (source 3)
 */
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

/**
 * Bfs to travers all node and find a coverTree. During the execution this
 * function also populate other support structures (source 1)
 *
 * @params
 * vector<vector<int>> graph --> graph
 * vector<vector<int>> coverTree --> cover tree of the graph
 * vector<int> parent --> mapping of nodes to the parent
 * unordered_map<pair<int, int>, bool, hash_pair> &siblings --> list of siblings
 * vector<bool> nodeIsCricca --> vector to check if a node belong to a clique
 *
 * COMPLEXITY: O(n + m)
 */
void bfs(vector<vector<int>> &graph, vector<vector<int>> &coverTree,
         vector<int> &parent,
         unordered_map<pair<int, int>, bool, hash_pair> &siblings,
         vector<bool> &nodeIsCricca) {
    queue<int> Q;
    vector<bool> visited;
    visited.resize(graph.size(), false);
    parent[0] = 0;
    Q.push(0);

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        for (int value : graph[u]) {
            if (!visited[value]) {
                visited[value] = true;
                coverTree[u].push_back(value);
                parent[value] = u;
                Q.push(value);
            } else if (value != parent[u]) {
                if (parent[value] == parent[u]) {
                    siblings[pair<int, int>(u, value)] = true;
                } else {
                    siblings[pair<int, int>(value, u)] = true;
                }
                nodeIsCricca[value] = true;
                nodeIsCricca[u] = true;
                nodeIsCricca[parent[u]] = true;
            }
        }
    }
}

/**
 * Recursive part of dfs with pre-order visit; this will
 * populate euler walk and other support structures for nodes of the graph
 * (source 2)
 *
 * @params
 * vector<vector<int>> graph --> graph
 * int node --> node used during a step
 * vector<bool> --> vector to map if a node is already visited during the dfs
 * int depth --> depth used during a step
 * vector<int> euler --> euler tour
 * vector<pair<int,int>> depths--> vector of the depths for the euler tour
 * vector<int> nodesMapping --> mapping of the nodes to the euler tour
 * vector<int> depthGraph --> list of depths of the coverTree
 *
 * COMPLEXITY: O(n + m)
 */
void dfsRec(vector<vector<int>> &graph, int &node, vector<bool> &visited,
            int depth, vector<int> &euler, vector<pair<int, int>> &depths,
            vector<int> &nodesMapping, vector<int> &depthGraph) {
    euler.push_back(node);
    depths.push_back(make_pair(depth, euler.size() - 1));
    if (!visited[node]) {
        nodesMapping[node] = euler.size() - 1;
    }
    visited[node] = true;

    for (int a : graph[node]) {
        if (!visited[a]) {
            depthGraph[a] = depthGraph[node] + 1;
            dfsRec(graph, a, visited, depth + 1, euler, depths, nodesMapping,
                   depthGraph);
            euler.push_back(node);
            depths.push_back(make_pair(depth, euler.size() - 1));
        }
    }
}

/**
 * This procedure just creates a vector to keep track of
 * which node gets visited and calls the recursion (source 2)
 *
 * @params
 * vector<vector<int>> graph --> graph
 * int node --> node used during a step
 * vector<int> euler --> euler tour
 * vector<pair<int,int>> depths--> vector of the depths for the euler tour
 * vector<int> nodesMapping --> mapping of the nodes to the euler tour
 * vector<int> depthGraph --> list of depths of the coverTree
 */
void dfs(vector<vector<int>> &graph, int node, vector<int> &euler,
         vector<pair<int, int>> &depths, vector<int> &nodesMapping,
         vector<int> &depthGraph) {
    nodesMapping.resize(graph.size());
    vector<bool> visited;
    int indexMapping = 0;
    depthGraph[node] = 0;
    visited.resize(graph.size(), false);
    dfsRec(graph, node, visited, indexMapping, euler, depths, nodesMapping,
           depthGraph);
}

/**
 * This procedure fills the matrix[log(euler.size()) + 1][euler.size()]
 * in complexity O(n*log(n)); this will allow to compute the
 * LCA (lowest common ancestor) of a pair of vertexes
 * in O(1). It takes O(n*log(n)) to store the matrix (source 4)
 *
 * @params
 * vector<vector<pair<int,int>>> sparseMatrix --> matrix to compute LCA
 * vector<pair<int,int>> depths--> vector of the depths for the euler tour
 */
void sparseMatrixComputation(vector<vector<pair<int, int>>> &sparseMatrix,
                             vector<pair<int, int>> &depths) {
    int n = depths.size();
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
}

/**
 * The __builtin_clz function returns how many zeros occur
 * before the first occurence of 1 in the binary representation
 * of a given number. Since integers in c++ are 32 bits values,
 * we subtract from 31 from the number of zeros to get p, defined
 * as the position of the most significant bit set to 1.
 * The function uses p to access the matrix and
 * find the Lowest Common Ancestor. (source 4)
 *
 * @params
 * vector<vector<pair<int,int>>> sparseMatrix --> matrix to compute LCA
 * vector<int> euler --> euler tour
 * int start --> start node
 * int end --> end node
 *
 * @return
 * int LCA --> LCA of the start, end nodes
 * */
int query(vector<vector<pair<int, int>>> &sparseMatrix, vector<int> &euler,
          int start, int end) {
    int p = 31 - __builtin_clz(end - start);
    pair<int, int> firstValue = sparseMatrix[p][start];
    pair<int, int> secondValue = sparseMatrix[p][end - (1 << p)];
    if (firstValue.first > secondValue.first) {
        return euler[secondValue.second];
    } else {
        return euler[firstValue.second];
    }
}

/**
 * This function given a request calculate the position inside of the euler tour
 * of the elements of the request. Then it uses them to calculate LCA. If LCA is
 * a clique there are some possibilities:
 *     1. the start node of the request is the LCA => the cost is the amount of;
 *        step to arrive to the LCA from the end node of the request
 *     2. same as 1 but with start and end swapped;
 *     3. nor start and end node of the request is the LCA => the cost is the
 *        amount of step to arrive to the LCA from both point of the request.
 *        Then is necessary to check if the last node before the LCA are
 *        siblings, if so the real cost is less than 1 of the path through the
 *        LCA.
 * If LCA is not a clique the cost is evaluated from the depths of the start,
 * end, LCA nodes.
 *
 * @params
 * int start --> start node of the request
 * int end --> end node of the request
 * vector<int> nodesMapping --> mapping of the nodes to the euler tour
 * vector<vector<pair<int,int>>> sparseMatrix --> matrix to compute LCA
 * vector<int> euler --> euler tour
 * vector<pair<int,int>> depths--> vector of the depths for the euler tour
 * vector<int> parent --> mapping of nodes to the parent
 * unordered_map<pair<int, int>, bool, hash_pair> &siblings --> list of siblings
 * vector<bool> nodeIsCricca --> vector to check if a node belong to a clique
 * vector<int> depthCoverTree --> list of depths of the coverTree
 *
 * @return
 * int cost --> the cost to go from start to end
 */
int handleRequest(int &start, int &end, vector<int> &nodesMapping,
                  vector<vector<pair<int, int>>> &sparseMatrix,
                  vector<int> &eueler, vector<pair<int, int>> &depths,
                  vector<int> &parent,
                  unordered_map<pair<int, int>, bool, hash_pair> &siblings,
                  vector<bool> &nodeIsCricca, vector<int> &depthCoverTree) {
    int cost = 0;
    int startMapped = nodesMapping[start];
    int endMapped = nodesMapping[end];
    int lowestCommonAncestor;
    if (startMapped > endMapped) {
        swap(startMapped, endMapped);
    }
    lowestCommonAncestor =
        query(sparseMatrix, eueler, startMapped, endMapped + 1);
    if (nodeIsCricca[lowestCommonAncestor]) {
        if (start == lowestCommonAncestor) {
            cost = depthCoverTree[end] - depthCoverTree[lowestCommonAncestor];
            return cost;

        } else if (end == lowestCommonAncestor) {
            cost = depthCoverTree[start] - depthCoverTree[lowestCommonAncestor];
            return cost;
        } else {
            while (parent[start] != lowestCommonAncestor) {
                start = parent[start];
                cost++;
            }
            while (parent[end] != lowestCommonAncestor) {
                end = parent[end];
                cost++;
            }
            if (start > end) {
                if (siblings.find(pair<int, int>(end, start)) !=
                    siblings.end()) {
                    cost++;
                } else {
                    cost = cost + 2;
                }
            } else {
                if (siblings.find(pair<int, int>(start, end)) !=
                    siblings.end()) {
                    cost++;
                } else {
                    cost = cost + 2;
                }
            }
        }
    } else {
        if (start == lowestCommonAncestor) {
            cost = (depthCoverTree[end] - depthCoverTree[start]);
        } else if (end == lowestCommonAncestor) {
            cost = (depthCoverTree[start] - depthCoverTree[end]);
        } else {
            cost = (depthCoverTree[start] + depthCoverTree[end] -
                    2 * (depthCoverTree[lowestCommonAncestor]));
        }
        return cost;
    }
    return cost;
}

int main() {
    ifstream in("input.txt");
    // Adjacency list for save graph from input
    vector<vector<int>> graph;
    // Adjacency list for the a cover tree of the graph
    vector<vector<int>> coverTree;
    // Vector to save the depth of the nodes inside the cover tree
    vector<int> depthCoverTree;
    // Vector to save parents of nodes (relation node<->parent)
    vector<int> parent;
    // Vector that map to true if a node belong to a clique
    vector<bool> nodeIsCricca;
    // Hash_map to save if there is a couple of nodes that are siblings
    unordered_map<pair<int, int>, bool, hash_pair> siblings;
    // Vector of requests read from input
    vector<pair<int, int>> request;
    // Vector to store the euler walk
    vector<int> euler;
    // Vector to store the depths of vertexes in the euler walk
    vector<pair<int, int>> depths;
    // This matrix will compute the access in O(1) to find operation of the LCA
    // (lowest common ancestor) of a pair of vertexes
    vector<vector<pair<int, int>>> sparseMatrix;
    // Vector to map nodes to positions iniside euler tour
    vector<int> nodesMapping;

    int N, M, Q;
    int node1, node2;

    in >> N >> M >> Q;

    graph.resize(N);
    coverTree.resize(N);
    depthCoverTree.resize(N);
    parent.resize(N, -1);
    nodeIsCricca.resize(N, false);

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

    in.close();

    // bfs to compute coverTree
    bfs(graph, coverTree, parent, siblings, nodeIsCricca);
    // dfs to compute euler tour and to populate other support structure
    dfs(coverTree, 0, euler, depths, nodesMapping, depthCoverTree);
    // Computation of the sparse matrix
    sparseMatrixComputation(sparseMatrix, depths);

    ofstream out("output.txt");
    for (pair<int, int> r : request) {
        out << handleRequest(r.first, r.second, nodesMapping, sparseMatrix,
                             euler, depths, parent, siblings, nodeIsCricca,
                             depthCoverTree)
            << endl;
    }
    out.close();
    return 0;
}
