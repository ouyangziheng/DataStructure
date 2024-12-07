#include <assert.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <iostream>
#include <queue>
#include <stack>
#include <vector>
using namespace std;

// 顶点
class Node {
public:
    Node(char data = 0) {
        this->data = data;
        is_visited = false;
    }
    Node(const Node &node) {
        if (this == &node)
            return;
        *this = node;
    }

    Node &operator=(const Node &node) {
        if (this == &node)
            return *this;
        this->data = node.data;
        this->is_visited = node.is_visited;
        return *this;
    }

public:
    char data;       // 数据
    bool is_visited; // 是否访问，用于遍历
};

// 边
class Edge {
public:
    Edge(int nodeIndexA = 0, int nodeIndexB = 0, int weightValue = 0)
        : node_index_a(nodeIndexA), node_index_b(nodeIndexB),
          weight_value(weightValue), is_selected(false) {}
    Edge(const Edge &edge) {
        if (this == &edge)
            return;
        *this = edge;
    }

    Edge &operator=(const Edge &edge) {
        if (this == &edge)
            return *this;
        this->node_index_a = edge.node_index_a;
        this->node_index_b = edge.node_index_b;
        this->weight_value = edge.weight_value;
        this->is_selected = edge.is_selected;
        return *this;
    }

public:
    int node_index_a; // 头顶点
    int node_index_b; // 尾顶点
    int weight_value; // 权重
    bool is_selected; // 是否被选中
};

// 图
class Mymap {
private:
    int capacity;          // 顶点总数
    int node_count;        // 当前顶点数量
    Node *node_array;      // 顶点集合
    int *adjacency_matrix; // 邻接距阵
    Edge *edge_array;      // 最小生成树边集合
public:
    Mymap(int iCapacity) {
        capacity = iCapacity;
        node_count = 0;
        node_array = new Node[capacity];
        adjacency_matrix = new int[capacity * capacity];
        memset(adjacency_matrix, 0, capacity * capacity * sizeof(int));
        edge_array = new Edge[capacity - 1];
    }
    ~Mymap(void) {
        delete[] node_array;
        delete[] adjacency_matrix;
        delete[] edge_array;
    }

private:
    // 广度遍历具体实现
    void breadthFirstTraverseImpl(vector<int> pre_vec) {
        int val = 0;
        vector<int> cur_vec;
        for (int i = 0; i < pre_vec.size(); i++) {
            for (int j = 0; j < capacity; j++) {
                getValueFromMatrix(pre_vec[i], j, val);
                if (val != 0) {
                    if (node_array[j].is_visited)
                        continue;
                    cout << node_array[j].data << " ";
                    node_array[j].is_visited = true;
                    cur_vec.push_back(j);
                } else
                    continue;
            }
        }
        if (cur_vec.empty())
            return;
        else
            breadthFirstTraverseImpl(cur_vec);
    }

    // 取最小边
    int getMinEdge(const vector<Edge> &edgeVec) {
        int min = 0, minEdge = 0;

        for (int i = 0; i < (int)edgeVec.size(); i++) {
            if (edgeVec[i].is_selected)
                continue;
            min = edgeVec[i].weight_value;
            minEdge = i;
        }

        for (int i = 0; i < (int)edgeVec.size(); i++) {
            if (edgeVec[i].is_selected)
                continue;
            if (min > edgeVec[i].weight_value) {
                min = edgeVec[i].weight_value;
                minEdge = i;
            }
        }

        if (min == 0)
            return -1;

        return minEdge;
    }

    bool isInSet(const vector<int> &nodeSet, int target) {
        for (int i = 0; i < (int)nodeSet.size(); i++) {
            if (nodeSet[i] == target)
                return true;
        }

        return false;
    }

    void mergeNodeSet(vector<int> &nodeSetA, const vector<int> &nodeSetB) {
        for (size_t i = 0; i < (int)nodeSetB.size(); i++) {
            nodeSetA.push_back(nodeSetB[i]);
        }
    }

public:
    // 添加顶点
    void addNode(Node *node) {
        node_array[node_count].data = node->data;
        node_count++;
    }
    // 将顶点访问设置默认
    void resetNode() {
        for (int i = 0; i < node_count; i++)
            node_array[i].is_visited = false;
    }
    // 设置权重-有向图
    bool setValueToMatrixForDirectedGraph(int row, int col, int val = 1) {
        if (row < 0 || row >= capacity)
            return false;
        if (col < 0 || col >= capacity)
            return false;
        adjacency_matrix[row * capacity + col] = val;
        return true;
    }

    // 设置权重-无向图
    bool setValueToMatrixForUndirectedGraph(int row, int col, int val = 1) {
        if (row < 0 || row >= capacity)
            return false;
        if (col < 0 || col >= capacity)
            return false;
        adjacency_matrix[row * capacity + col] = val;
        adjacency_matrix[col * capacity + row] = val;
        return true;
    }
    // 获取权重，进行基础判断
    bool getValueFromMatrix(int row, int col, int &val) {
        if (row < 0 || row >= capacity)
            return false;
        if (col < 0 || col >= capacity)
            return false;
        val = adjacency_matrix[row * capacity + col];
        return true;
    }
    // 打印矩阵
    void printMatrix() {
        for (int i = 0; i < capacity; i++) {
            for (int j = 0; j < capacity; j++)
                cout << adjacency_matrix[i * capacity + j] << " ";
            cout << endl;
        }
    }

    // 深度遍历
    void depthFirstTraverse(int index) {
        int val = 0;
        cout << node_array[index].data << " ";
        node_array[index].is_visited = true;

        for (int i = 0; i < capacity; i++) {
            getValueFromMatrix(index, i, val);
            if (val != 0) {
                if (node_array[i].is_visited)
                    continue;
                depthFirstTraverse(i);
            } else
                continue;
        }
    }
    void topologicalSortWithStack() {
        stack<int> s;
        vector<bool> visited(capacity, false);

        for (int i = 0; i < capacity; i++) {
            if (!visited[i]) {
                depthFirstSearchWithStack(i, s, visited);
            }
        }

        while (!s.empty()) {
            cout << node_array[s.top()].data << " ";
            s.pop();
        }
    }

    void depthFirstSearchWithStack(int node, stack<int> &s,
                                   vector<bool> &visited) {
        visited[node] = true;

        for (int i = 0; i < capacity; i++) {
            int val = 0;
            getValueFromMatrix(node, i, val);

            if (val != 0 && !visited[i]) {
                depthFirstSearchWithStack(i, s, visited);
            }
        }
        s.push(node);
    }

    void dijkstra(int start) {
        vector<int> dist(capacity, INT_MAX);
        vector<bool> visited(capacity, false);
        dist[start] = 0;

        priority_queue<pair<int, int>, vector<pair<int, int>>,
                       greater<pair<int, int>>>
            pq;
        pq.push({0, start});

        while (!pq.empty()) {
            int currentNode = pq.top().second;
            int currentDist = pq.top().first;
            pq.pop();

            if (visited[currentNode])
                continue;

            visited[currentNode] = true;

            for (int i = 0; i < capacity; i++) {
                int val = 0;
                getValueFromMatrix(currentNode, i, val);

                if (val != 0 && !visited[i] && currentDist + val < dist[i]) {
                    dist[i] = currentDist + val;
                    pq.push({dist[i], i});
                }
            }
        }

        for (int i = 0; i < capacity; i++) {
            if (dist[i] == INT_MAX) {
                cout << "从节点 " << start << " 到节点 " << i << " 不可达"
                     << endl;
            } else {
                cout << "从节点 " << start << " 到节点 " << i
                     << " 的最短路径为 " << dist[i] << endl;
            }
        }
    }

    bool isDirectedGraph() {
        for (int i = 0; i < capacity; i++) {
            for (int j = 0; j < capacity; j++) {
                if (adjacency_matrix[i * capacity + j] !=
                    adjacency_matrix[j * capacity + i]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool dfsForCycle(int node, vector<bool> &visited, vector<bool> &recStack) {
        visited[node] = true;
        recStack[node] = true;

        for (int i = 0; i < capacity; i++) {
            if (adjacency_matrix[node * capacity + i] != 0) {
                if (!visited[i] && dfsForCycle(i, visited, recStack)) {
                    return true;
                } else if (recStack[i]) {
                    return true;
                }
            }
        }

        recStack[node] = false;
        return false;
    }

    bool isDAG() {
        vector<bool> visited(capacity, false);
        vector<bool> recStack(capacity, false);

        for (int i = 0; i < capacity; i++) {
            if (!visited[i]) {
                if (dfsForCycle(i, visited, recStack)) {
                    return false;
                }
            }
        }
        return true;
    }

    void dfsForConnectivity(int node, vector<bool> &visited) {
        visited[node] = true;
        for (int i = 0; i < capacity; i++) {
            if (adjacency_matrix[node * capacity + i] != 0 && !visited[i]) {
                dfsForConnectivity(i, visited);
            }
        }
    }

    bool isConnectedGraph() {
        vector<bool> visited(capacity, false);
        dfsForConnectivity(0, visited); // 从节点0开始DFS遍历

        // 检查是否所有节点都被访问过
        for (int i = 0; i < capacity; i++) {
            if (!visited[i]) {
                return false; // 有节点未被访问，说明不连通
            }
        }
        return true; // 所有节点都被访问过，图是连通的
    }

    // 广度遍历
    void breadthFirstTraverse(int index) {
        cout << node_array[index].data << " ";
        node_array[index].is_visited = true;

        vector<int> curVec;
        curVec.push_back(index);

        breadthFirstTraverseImpl(curVec);
    }

    // 求最小生成树-普里斯算法
    void primTree(int index) {
        int val = 0;
        int iEdgeCount = 0;
        vector<Edge> edgeVec;

        vector<int> nodeIndexVec;
        nodeIndexVec.push_back(index);

        while (iEdgeCount < capacity - 1) {
            int row = nodeIndexVec.back();
            cout << node_array[row].data << endl;
            node_array[row].is_visited = true;

            for (int i = 0; i < capacity; i++) {
                getValueFromMatrix(row, i, val);
                if (val == 0)
                    continue;
                if (node_array[i].is_visited)
                    continue;
                Edge edge(row, i, val);
                edgeVec.push_back(edge);
            }

            // 取出最小边
            int retIndex = getMinEdge(edgeVec);
            if (retIndex != -1) {
                edgeVec[retIndex].is_selected = true;
                edge_array[iEdgeCount] = edgeVec[retIndex];
                cout << node_array[edge_array[iEdgeCount].node_index_a].data
                     << " - ";
                cout << node_array[edge_array[iEdgeCount].node_index_b].data
                     << " (";
                cout << edge_array[iEdgeCount].weight_value << ") " << endl;
                iEdgeCount++;

                int iNodeIndex = edgeVec[retIndex].node_index_b;
                node_array[iNodeIndex].is_visited = true;
                nodeIndexVec.push_back(iNodeIndex);
            }
        }
    }
    void bellmanFord(int start) {
        vector<int> dist(capacity, INT_MAX);
        dist[start] = 0;

        for (int i = 1; i < capacity; i++) {
            for (int u = 0; u < capacity; u++) {
                for (int v = 0; v < capacity; v++) {
                    int val = 0;
                    getValueFromMatrix(u, v, val);
                    if (val != 0 && dist[u] != INT_MAX &&
                        dist[u] + val < dist[v]) {
                        dist[v] = dist[u] + val;
                    }
                }
            }
        }

        for (int u = 0; u < capacity; u++) {
            for (int v = 0; v < capacity; v++) {
                int val = 0;
                getValueFromMatrix(u, v, val);
                if (val != 0 && dist[u] != INT_MAX && dist[u] + val < dist[v]) {
                    cout << "图中存在负权环！" << endl;
                    return;
                }
            }
        }

        for (int i = 0; i < capacity; i++) {
            if (dist[i] == INT_MAX) {
                cout << "从节点 " << start << " 到节点 " << i << " 不可达"
                     << endl;
            } else {
                cout << "从节点 " << start << " 到节点 " << i
                     << " 的最短路径为 " << dist[i] << endl;
            }
        }
    }

    // 最小生成树-克鲁斯卡尔算法
    void kruskalTree() {
        int val = 0;
        int edgeCount = 0;

        vector<vector<int>> nodeSets;

        vector<Edge> edgeVec;
        for (int i = 0; i < capacity; i++) {
            for (int j = i + 1; j < capacity; j++) {
                getValueFromMatrix(i, j, val);
                if (0 == val)
                    continue;
                if (node_array[i].is_visited)
                    continue;
                Edge edge(i, j, val);
                edgeVec.push_back(edge);
            }
        }

        while (edgeCount < capacity - 1) {
            int retIndex = getMinEdge(edgeVec);
            if (-1 != retIndex) {
                edgeVec[retIndex].is_selected = true;

                int nodeAIndex = edgeVec[retIndex].node_index_a;
                int nodeBIndex = edgeVec[retIndex].node_index_b;

                bool nodeAInSet = false;
                bool nodeBInSet = false;
                int nodeAInSetLabel = -1;
                int nodeBInSetLabel = -1;

                for (int i = 0; i < (int)nodeSets.size(); i++) {
                    nodeAInSet = isInSet(nodeSets[i], nodeAIndex);
                    if (nodeAInSet)
                        nodeAInSetLabel = i;
                }

                for (int i = 0; i < (int)nodeSets.size(); i++) {
                    nodeBInSet = isInSet(nodeSets[i], nodeBIndex);
                    if (nodeBInSet)
                        nodeBInSetLabel = i;
                }

                if (nodeAInSetLabel == -1 && nodeBInSetLabel == -1) {
                    vector<int> vec;
                    vec.push_back(nodeAIndex);
                    vec.push_back(nodeBIndex);
                    nodeSets.push_back(vec);
                } else if (nodeAInSetLabel == -1 && nodeBInSetLabel != -1) {
                    nodeSets[nodeBInSetLabel].push_back(nodeAIndex);
                } else if (nodeAInSetLabel != -1 && nodeBInSetLabel == -1) {
                    nodeSets[nodeAInSetLabel].push_back(nodeBIndex);
                } else if (-1 != nodeAInSetLabel && -1 != nodeBInSetLabel &&
                           nodeAInSetLabel != nodeBInSetLabel) {

                    nodeSets[nodeAInSetLabel].insert(
                        nodeSets[nodeAInSetLabel].end(),
                        nodeSets[nodeBInSetLabel].begin(),
                        nodeSets[nodeBInSetLabel].end());
                    for (int k = nodeBInSetLabel; k < (int)nodeSets.size() - 1;
                         k++) {
                        nodeSets[k] = nodeSets[k + 1];
                    }
                } else if (nodeAInSetLabel != -1 && nodeBInSetLabel != -1 &&
                           nodeAInSetLabel == nodeBInSetLabel) {
                    continue;
                }

                edge_array[edgeCount] = edgeVec[retIndex];
                edgeCount++;

                cout << node_array[edgeVec[retIndex].node_index_a].data
                     << " - ";
                cout << node_array[edgeVec[retIndex].node_index_b].data << " (";
                cout << edgeVec[retIndex].weight_value << ") " << endl;
            }
        }
    }
};

void implementAlgorithmsSet(Mymap *myMap) {
    cout << "打印矩阵: " << endl;
    myMap->printMatrix();
    cout << endl;

    myMap->resetNode();

    cout << "是否为有向图？" << endl;
    cout << endl;
    if (myMap->isDirectedGraph()) {
        cout << "是有向图" << endl;
        cout << endl;

        cout << "普里姆算法: " << endl;
        myMap->primTree(0);
        cout << endl;

        myMap->resetNode();

        cout << "克鲁斯卡尔算法: " << endl;
        myMap->kruskalTree();
        cout << endl;

        myMap->resetNode();

    } else {
        cout << "是无向图" << endl;
        cout << endl;
    }

    cout << "是否为有向无环图(DAG)？" << endl;
    cout << endl;

    if (myMap->isDAG()) {
        cout << "是DAG" << endl;
        cout << endl;

        cout << "对DAG进行拓扑排序:" << endl;

        myMap->topologicalSortWithStack();

        cout << endl;

    } else {
        cout << "不是DAG" << endl;
    }

    cout << "是否为连通图？" << endl;
    if (myMap->isConnectedGraph()) {
        cout << "是连通图" << endl;
    } else {
        cout << "不是连通图" << endl;
    }

    cout << "深度优先遍历: " << endl;
    myMap->depthFirstTraverse(0);
    cout << endl;

    myMap->resetNode();

    myMap->resetNode();

    cout << "广度优先遍历: " << endl;
    myMap->breadthFirstTraverse(0);
    cout << endl;

    myMap->resetNode();
    cout << endl;

    cout << "Dijkstra 算法最短路径：" << endl;
    myMap->dijkstra(0);
    cout << endl;

    cout << "BellmanFord 算法最短路径：" << endl;
    myMap->bellmanFord(0);
    cout << endl;
    cout << endl;
}

void test01() {
    Mymap *myMap = new Mymap(6);

    Node *pNodeA = new Node('A');
    Node *pNodeB = new Node('B');
    Node *pNodeC = new Node('C');
    Node *pNodeD = new Node('D');
    Node *pNodeE = new Node('E');
    Node *pNodeF = new Node('F');

    myMap->addNode(pNodeA);
    myMap->addNode(pNodeB);
    myMap->addNode(pNodeC);
    myMap->addNode(pNodeD);
    myMap->addNode(pNodeE);
    myMap->addNode(pNodeF);

    myMap->setValueToMatrixForDirectedGraph(0, 1, 7);  // 0 -> 1
    myMap->setValueToMatrixForDirectedGraph(0, 2, 1);  // 0 -> 2
    myMap->setValueToMatrixForDirectedGraph(0, 3, 9);  // 0 -> 3
    myMap->setValueToMatrixForDirectedGraph(1, 2, 2);  // 1 -> 2
    myMap->setValueToMatrixForDirectedGraph(1, 4, 3);  // 1 -> 4
    myMap->setValueToMatrixForDirectedGraph(2, 3, 11); // 2 -> 3
    myMap->setValueToMatrixForDirectedGraph(2, 4, 8);  // 2 -> 4
    myMap->setValueToMatrixForDirectedGraph(2, 5, 4);  // 2 -> 5
    myMap->setValueToMatrixForDirectedGraph(3, 5, 5);  // 3 -> 5
    myMap->setValueToMatrixForDirectedGraph(4, 5, 15); // 4 -> 5

    implementAlgorithmsSet(myMap);
}

void test02() {
    Mymap *myMap = new Mymap(6);

    Node *pNodeA = new Node('A');
    Node *pNodeB = new Node('B');
    Node *pNodeC = new Node('C');
    Node *pNodeD = new Node('D');
    Node *pNodeE = new Node('E');
    Node *pNodeF = new Node('F');

    myMap->addNode(pNodeA);
    myMap->addNode(pNodeB);
    myMap->addNode(pNodeC);
    myMap->addNode(pNodeD);
    myMap->addNode(pNodeE);
    myMap->addNode(pNodeF);

    // 添加有向边
    myMap->setValueToMatrixForDirectedGraph(0, 1, 5); // A -> B
    myMap->setValueToMatrixForDirectedGraph(1, 2, 7); // B -> C
    myMap->setValueToMatrixForDirectedGraph(2, 3, 3); // C -> D
    myMap->setValueToMatrixForDirectedGraph(3, 4, 4); // D -> E
    myMap->setValueToMatrixForDirectedGraph(4, 5, 6); // E -> F
    implementAlgorithmsSet(myMap);
}

void test03() {
    Mymap *myMap = new Mymap(6);

    Node *pNodeA = new Node('A');
    Node *pNodeB = new Node('B');
    Node *pNodeC = new Node('C');
    Node *pNodeD = new Node('D');
    Node *pNodeE = new Node('E');
    Node *pNodeF = new Node('F');

    myMap->addNode(pNodeA);
    myMap->addNode(pNodeB);
    myMap->addNode(pNodeC);
    myMap->addNode(pNodeD);
    myMap->addNode(pNodeE);
    myMap->addNode(pNodeF);

    // 添加无向边（双向）
    myMap->setValueToMatrixForUndirectedGraph(0, 1, 5); // A -- B
    myMap->setValueToMatrixForUndirectedGraph(1, 2, 7); // B -- C
    myMap->setValueToMatrixForUndirectedGraph(2, 3, 3); // C -- D
    myMap->setValueToMatrixForUndirectedGraph(3, 4, 4); // D -- E
    myMap->setValueToMatrixForUndirectedGraph(4, 5, 6); // E -- F
    implementAlgorithmsSet(myMap);
}

void test04() {
    Mymap *myMap = new Mymap(4);

    Node *pNodeA = new Node('A');
    Node *pNodeB = new Node('B');
    Node *pNodeC = new Node('C');
    Node *pNodeD = new Node('D');

    myMap->addNode(pNodeA);
    myMap->addNode(pNodeB);
    myMap->addNode(pNodeC);
    myMap->addNode(pNodeD);

    // 添加有向边，所有边的权重为5
    myMap->setValueToMatrixForDirectedGraph(0, 1, 5); // A -> B
    myMap->setValueToMatrixForDirectedGraph(1, 2, 5); // B -> C
    myMap->setValueToMatrixForDirectedGraph(2, 3, 5); // C -> D
    myMap->setValueToMatrixForDirectedGraph(3, 0, 5); // D -> A (形成环)
    implementAlgorithmsSet(myMap);
}

int main() {
    test01();
    test02();
    test03();
    test04();

    system("pause");
    return 0;
}