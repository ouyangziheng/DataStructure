#include <algorithm>
#include <assert.h>
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
        data = data;
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
class CMap {
private:
    int capacity;          // 顶点总数
    int node_count;        // 当前顶点数量
    Node *node_array;      // 顶点集合
    int *adjacency_matrix; // 邻接距阵
    Edge *edge_array;      // 最小生成树边集合
public:
    CMap(int iCapacity) {
        capacity = iCapacity;
        node_count = 0;
        node_array = new Node[capacity];
        adjacency_matrix = new int[capacity * capacity];
        memset(adjacency_matrix, 0, capacity * capacity * sizeof(int));
        edge_array = new Edge[capacity - 1];
    }
    ~CMap(void) {
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
        assert(node);
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
        stack<int> s; // 创建一个栈来存储拓扑排序的结果
        vector<bool> visited(capacity,
                             false); // 访问标记数组，记录节点是否已访问

        // 遍历所有节点，如果节点未被访问，则进行深度优先遍历
        for (int i = 0; i < capacity; i++) {
            if (!visited[i]) {
                depthFirstSearchWithStack(i, s, visited); // 从未访问节点开始DFS
            }
        }

        // 打印拓扑排序的结果（栈中存储的是逆拓扑排序）
        while (!s.empty()) {
            cout << node_array[s.top()].data << " "; // 输出栈中的节点数据
            s.pop();                                 // 弹出栈顶元素
        }
    }

    // 深度优先遍历辅助函数，用于拓扑排序
    void depthFirstSearchWithStack(int node, stack<int> &s,
                                   vector<bool> &visited) {
        visited[node] = true; // 标记当前节点为已访问

        // 遍历当前节点的所有邻接节点
        for (int i = 0; i < capacity; i++) {
            int val = 0;
            getValueFromMatrix(node, i, val); // 获取从node到i的边权重

            // 如果存在边且i未被访问，递归访问该邻接节点
            if (val != 0 && !visited[i]) {
                depthFirstSearchWithStack(i, s,
                                          visited); // 深度优先遍历邻接节点
            }
        }
        // 将当前节点压入栈中，表示当前节点及其所有邻接节点都已处理完
        s.push(node);
    }
    void dijkstra(int start) {
        vector<int> dist(capacity, INT_MAX); // 存储每个节点的最短距离
        vector<bool> visited(capacity, false); // 标记每个节点是否已访问
        dist[start] = 0; // 起始节点到自己的距离为0

        // 使用优先队列（最小堆），存储<距离, 节点>对
        priority_queue<pair<int, int>, vector<pair<int, int>>,
                       greater<pair<int, int>>>
            pq;
        pq.push({0, start}); // 起始节点入队

        while (!pq.empty()) {
            int currentNode = pq.top().second;
            int currentDist = pq.top().first;
            pq.pop();

            if (visited[currentNode])
                continue; // 如果该节点已访问，则跳过

            visited[currentNode] = true;

            // 遍历所有邻接节点，进行松弛操作
            for (int i = 0; i < capacity; i++) {
                int val = 0;
                getValueFromMatrix(currentNode, i,
                                   val); // 获取当前节点到i的边权重

                if (val != 0 && !visited[i] && currentDist + val < dist[i]) {
                    dist[i] = currentDist + val; // 松弛操作
                    pq.push({dist[i], i}); // 将更新后的节点推入优先队列
                }
            }
        }

        // 输出从起点到各个节点的最短距离
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
                // 如果发现有邻接关系的对称点不一致，则为有向图
                if (adjacency_matrix[i * capacity + j] !=
                    adjacency_matrix[j * capacity + i]) {
                    return true; // 是有向图
                }
            }
        }
        return false; // 是无向图
    }

    bool dfsForCycle(int node, vector<bool> &visited, vector<bool> &recStack) {
        visited[node] = true;
        recStack[node] = true;

        for (int i = 0; i < capacity; i++) {
            if (adjacency_matrix[node * capacity + i] != 0) {
                if (!visited[i] && dfsForCycle(i, visited, recStack)) {
                    return true; // 检测到环
                } else if (recStack[i]) {
                    return true; // 检测到环
                }
            }
        }

        recStack[node] = false;
        return false;
    }

    bool isDAG() {
        vector<bool> visited(capacity, false);
        vector<bool> recStack(capacity, false); // 记录递归栈，判断环

        for (int i = 0; i < capacity; i++) {
            if (!visited[i]) {
                if (dfsForCycle(i, visited, recStack)) {
                    return false; // 如果有环，说明不是DAG
                }
            }
        }
        return true; // 没有环，说明是DAG
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
        vector<Edge> edgeVec; // 待选边集合

        // 从传入点开始找
        vector<int> nodeIndexVec;
        nodeIndexVec.push_back(index);

        // 结束条件：边数=顶点数-1
        while (iEdgeCount < capacity - 1) {
            // 查找传入点的符合要求（权重不为0且目的点没有被访问）边
            int row = nodeIndexVec.back();
            cout << node_array[row].data << endl;
            node_array[row].is_visited = true;

            for (int i = 0; i < capacity; i++) {
                getValueFromMatrix(row, i, val);
                if (0 == val)
                    continue;
                if (node_array[i].is_visited)
                    continue;
                Edge edge(row, i, val);
                edgeVec.push_back(edge);
            }

            // 取出最小边
            int retIndex = getMinEdge(edgeVec);
            if (-1 != retIndex) {
                // 存储选中边
                edgeVec[retIndex].is_selected = true;
                edge_array[iEdgeCount] = edgeVec[retIndex];
                cout << node_array[edge_array[iEdgeCount].node_index_a].data
                     << " - ";
                cout << node_array[edge_array[iEdgeCount].node_index_b].data
                     << " (";
                cout << edge_array[iEdgeCount].weight_value << ") " << endl;
                iEdgeCount++;

                int iNodeIndex = edgeVec[retIndex].node_index_b;
                // 设置点被访问
                node_array[iNodeIndex].is_visited = true;
                // 存入目的点递归查找
                nodeIndexVec.push_back(iNodeIndex);
            }
        }
    }
    void bellmanFord(int start) {
        vector<int> dist(capacity, INT_MAX); // 存储每个节点的最短距离
        dist[start] = 0; // 起始节点到自己的距离为0

        // 放松操作，最多进行 V-1 次
        for (int i = 1; i < capacity; i++) {
            for (int u = 0; u < capacity; u++) {
                for (int v = 0; v < capacity; v++) {
                    int val = 0;
                    getValueFromMatrix(u, v, val);
                    if (val != 0 && dist[u] != INT_MAX &&
                        dist[u] + val < dist[v]) {
                        dist[v] = dist[u] + val; // 松弛操作
                    }
                }
            }
        }

        // 检测是否存在负权环
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

        // 输出从起点到各个节点的最短路径
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

        // 定义存放节点集合数组
        vector<vector<int>> nodeSets;

        // 第一步、取出所有边
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

        // 第二步、从所有边中取出组成最小生成树的边
        // 1、算法结束条件：边数=顶点数-1
        while (edgeCount < capacity - 1) {
            // 2、从边集合中找出最小边
            int retIndex = getMinEdge(edgeVec);
            if (-1 != retIndex) {
                edgeVec[retIndex].is_selected = true;

                // 3、找出最小边连接点
                int nodeAIndex = edgeVec[retIndex].node_index_a;
                int nodeBIndex = edgeVec[retIndex].node_index_b;

                // 4、找出点所在集合
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

                // 5、根据点集合的不同做不同处理
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
                    // mergeNodeSet(nodeSets[nodeAInSetLabel],
                    // nodeSets[nodeBInSetLabel]);
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

int main(int argc, char **argv) {

    CMap *pMap = new CMap(6);

    Node *pNodeA = new Node('A');
    Node *pNodeB = new Node('B');
    Node *pNodeC = new Node('C');
    Node *pNodeD = new Node('D');
    Node *pNodeE = new Node('E');
    Node *pNodeF = new Node('F');

    pMap->addNode(pNodeA);
    pMap->addNode(pNodeB);
    pMap->addNode(pNodeC);
    pMap->addNode(pNodeD);
    pMap->addNode(pNodeE);
    pMap->addNode(pNodeF);

    pMap->setValueToMatrixForDirectedGraph(0, 1, 7);  // 0 -> 1
    pMap->setValueToMatrixForDirectedGraph(0, 2, 1);  // 0 -> 2
    pMap->setValueToMatrixForDirectedGraph(0, 3, 9);  // 0 -> 3
    pMap->setValueToMatrixForDirectedGraph(1, 2, 2);  // 1 -> 2
    pMap->setValueToMatrixForDirectedGraph(1, 4, 3);  // 1 -> 4
    pMap->setValueToMatrixForDirectedGraph(2, 3, 11); // 2 -> 3
    pMap->setValueToMatrixForDirectedGraph(2, 4, 8);  // 2 -> 4
    pMap->setValueToMatrixForDirectedGraph(2, 5, 4);  // 2 -> 5
    pMap->setValueToMatrixForDirectedGraph(3, 5, 5);  // 3 -> 5
    pMap->setValueToMatrixForDirectedGraph(4, 5, 15); // 4 -> 5

    cout << "打印矩阵: " << endl;
    pMap->printMatrix();
    cout << endl;

    pMap->resetNode();

    cout << "是否为有向图？" << endl;
    if (pMap->isDirectedGraph()) {
        cout << "是有向图" << endl;
    } else {
        cout << "是无向图" << endl;
    }

    // 继续判断是否是DAG
    cout << "是否为有向无环图(DAG)？" << endl;
    if (pMap->isDAG()) {
        cout << "是DAG" << endl;
    } else {
        cout << "不是DAG" << endl;
    }

    // 判断是否是连通图
    cout << "是否为连通图？" << endl;
    if (pMap->isConnectedGraph()) {
        cout << "是连通图" << endl;
    } else {
        cout << "不是连通图" << endl;
    }

    cout << "深度优先遍历: " << endl;
    pMap->depthFirstTraverse(0);
    cout << endl;

    pMap->resetNode();

    cout << "拓扑排序: " << endl;
    pMap->topologicalSortWithStack();
    cout << endl;

    pMap->resetNode();

    cout << "广度优先遍历: " << endl;
    pMap->breadthFirstTraverse(0);
    cout << endl;

    pMap->resetNode();

    cout << "普里姆算法: " << endl;
    pMap->primTree(0);
    cout << endl;

    pMap->resetNode();

    cout << "克鲁斯卡尔算法: " << endl;
    pMap->kruskalTree();
    cout << endl;

    pMap->resetNode();

    cout << "Dijkstra 算法最短路径：" << endl;
    pMap->dijkstra(0);

    cout << "BellmanFord 算法最短路径：" << endl;
    pMap->bellmanFord(0);

    system("pause");
    return 0;
}