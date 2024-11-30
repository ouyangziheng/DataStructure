// 错误提示
// 1. 使用==进行比较
// 2. 在交换的时候和运算的时候 对参数本身和参数的引用混淆
// 3. 注意一次循环的结束，注意单if对后面的影响，思考出循环的条件，防止出现死循环

#include <algorithm>
#include <climits> 
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <functional>

using namespace std;

// 计数器
int compareCount = 0;
int swapCount = 0;

// 重置计数器
void resetCounters() {
    compareCount = 0;
    swapCount = 0;
}

// 随机生成数据并保存到文件
void generateRandomData(const string &filename, int n) {
    srand(time(0));
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "无法打开文件!" << endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        outFile << rand() % 6666 << " "; 
    }
    outFile.close();
    cout << "数据已成功生成并保存在 " << filename << " 文件中。" << endl;
}

// 从文件中读取数据
void readDataFromFile(const string &filename, vector<int> &vec) {
    ifstream inFile(filename);
    int num;
    while (inFile >> num) {
        vec.push_back(num);
    }
    inFile.close();
}

void saveToFile(const vector<int> &vec, const string &filename) {
    if (remove(filename.c_str()) != 0) {
        cout << "无法删除原文件 " << filename << endl;
    } else {
        cout << "原文件已删除。" << endl;
    }

    ofstream outFile(filename);
    if (outFile.is_open()) {
        for (int num : vec) {
            outFile << num << " ";
        }
        outFile.close();
        cout << "排序结果已保存到文件：" << filename << endl;
    } else {
        cout << "无法打开文件 " << filename << " 进行写入" << endl;
    }
}

// 检查文件是否被排序
bool isFileSorted(const std::string &filename) {
    std::ifstream infile(filename);

    if (!infile.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return false;
    }

    int prev, current;
    if (!(infile >> prev)) {
        infile.close();
        return true;
    }

    while (infile >> current) {
        if (current < prev) {
            cout << current << " " << prev;
            infile.close();
            return false; 
        }
        prev = current;
    }

    infile.close();
    return true; 
}

// 冒泡排序
void bubbleSort(vector<int> &vec) {
    int n = vec.size();

    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            compareCount++;
            if (vec[j] > vec[j + 1]) {
                swapCount++;
                swap(vec[j], vec[j + 1]);
            }
        }
    }
}

// 优化选择排序
void optimizedBubbleSort(vector<int> &vec) {
    int n = vec.size();
    for (int i = 0; i < n - 1; i++) {
        bool ifSwapped = 0;
        for (int j = 0; j < n - 1 - i; j++) {
            compareCount++;
            if (vec[j] > vec[j + 1]) {
                swap(vec[j], vec[j + 1]);
                ifSwapped = 1;
                swapCount++;
            }
        }
        if (ifSwapped == 0) {
            break;
        }
    }
}

// 双向冒泡排序
void cocktailShakerSort(vector<int> &vec) {
    int n = vec.size();
    bool ifSwapped = true;
    int begin = 0;
    int end = n - 1;
    while (ifSwapped == true) {
        ifSwapped = false;
        for (int i = begin; i < end; i++) {
            compareCount++;
            if (vec[i] > vec[i + 1]) {
                swap(vec[i + 1], vec[i]);
                swapCount++;
                ifSwapped = true;
            }
        }

        end--;

        for (int i = end; i > begin; i--) {
            compareCount++;
            if (vec[i - 1] > vec[i]) {
                swap(vec[i], vec[i - 1]);
                swapCount++;
                ifSwapped = true;
            }
        }
        begin++;
    }
}

// 选择排序
void selectionSort(vector<int> &vec) {
    int n = vec.size();
    for (int i = 0; i < n - 1; ++i) {
        int minIndex = i;
        for (int j = i + 1; j < n; ++j) {
            compareCount++; 
            if (vec[j] < vec[minIndex]) {
                minIndex = j;
            }
        }
        if (minIndex != i) {
            swapCount++; 
            swap(vec[i], vec[minIndex]);
        }
    }
}

// 二分查找，返回插入位置
int binarySearch(const vector<int> &vec, int left, int right, int key) {
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == key) {
            compareCount++;
            return mid;
        } else if (vec[mid] < key) {
            compareCount++;
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left; 
}

// 使用二分查找优化的插入排序
void insertionSortWithBinarySearch(vector<int> &vec) {
    int n = vec.size();
    for (int i = 1; i < n; i++) {
        int key = vec[i]; 
        int insertPos =
            binarySearch(vec, 0, i - 1, key); 

        for (int j = i; j > insertPos; j--) {
            vec[j] = vec[j - 1];
            swapCount++;
        }
        vec[insertPos] = key;
    }
}

// 插入排序
void insertionSort(vector<int> &vec) {
    int n = vec.size();
    int index = 1;
    for (int i = 0; i < n - 1; i++) {
        for (int j = index; j > 0; j--) {
            compareCount++;
            if (vec[j] < vec[j - 1]) {
                swap(vec[j], vec[j - 1]);
                swapCount++;
            }
        }
        index++;
    }
}

// 快排辅助，获取中间元素索引
int partition(vector<int> &vec, int left, int right) {
    int index = left - 1;
    int pivot = vec[right];
    for (int i = left; i < right; i++) {
        compareCount++;
        if (vec[i] < pivot) {
            swap(vec[i], vec[++index]);
            swapCount++;
        }
    }
    swap(vec[right], vec[index + 1]);
    return index + 1;
}

// 快排
void quickSort(vector<int> &vec, int left, int right) {
    if (right <= left) {
        compareCount++;
        return;
    }
    int mid = partition(vec, left, right);
    quickSort(vec, left, mid - 1);
    quickSort(vec, mid + 1, right);
    return;
}

// 治：合并操作
void merge(vector<int> &vec, int left, int mid, int right) {
    if (left == right)
        return;
    int n1 = mid - left + 1;
    int n2 = right - mid;
    vector<int> leftVec, rightVec;

    for (int i = 0; i < n1; i++) {
        leftVec.push_back(vec[left + i]);
    }
    for (int i = 0; i < n2; i++) {
        rightVec.push_back(vec[i + mid + 1]);
    }

    int p = 0;
    int q = 0;

    for (int i = left; i < right + 1; i++) {
        compareCount++;
        if (p < n1 && q < n2) {
            if (leftVec[p] <= rightVec[q]) {
                vec[i] = leftVec[p++];
            } else {
                vec[i] = rightVec[q++];
            }
        } else if (p == n1)
            vec[i] = rightVec[q++];
        else
            vec[i] = leftVec[p++];
    }
}

// 分：分开数据，默认左边数组携带 mid 元素
void mergeSort(vector<int> &vec, int left, int right) {
    if (left >= right) {
        return;
    }
    int mid = left + (right - left) / 2;
    mergeSort(vec, left, mid);
    mergeSort(vec, mid + 1, right);
    merge(vec, left, mid, right);
}

// 最小堆实现
class MinHeap {
  public:
    vector<int> heapElements;

    MinHeap(const vector<int> &inputVector) { buildHeap(inputVector); }

    int getLeftChildIndex(int currentIndex) const {
        return 2 * currentIndex + 1;
    }

    int getRightChildIndex(int currentIndex) const {
        return 2 * currentIndex + 2;
    }

    bool isLeafNode(int currentIndex) const {
        int heapSize = heapElements.size();
        return getLeftChildIndex(currentIndex) >= heapSize;
    }

    int getParentIndex(int currentIndex) const {
        if (currentIndex == 0) {
            return -1;
        }
        return (currentIndex - 1) / 2;
    }

    void siftUp(int currentIndex) {
        while (currentIndex != 0) {
            int parentIndex = getParentIndex(currentIndex);
            compareCount++;
            if (heapElements[currentIndex] < heapElements[parentIndex]) {
                swapCount++;
                swap(heapElements[currentIndex], heapElements[parentIndex]);
                currentIndex = parentIndex;
            } else {
                break;
            }
        }
    }

    void siftDown(int currentIndex) {
        int heapSize = heapElements.size();
        while (true) {
            int leftIndex = getLeftChildIndex(currentIndex);
            int rightIndex = getRightChildIndex(currentIndex);
            int smallestIndex = currentIndex;

            // 比较左
            if (leftIndex < heapSize &&
                heapElements[leftIndex] < heapElements[smallestIndex]) {
                smallestIndex = leftIndex;
                compareCount++;
            }

            // 比较右
            if (rightIndex < heapSize &&
                heapElements[rightIndex] < heapElements[smallestIndex]) {
                smallestIndex = rightIndex;
                compareCount++;
            }

            // 下沉
            if (smallestIndex != currentIndex) {
                swapCount++;
                swap(heapElements[currentIndex], heapElements[smallestIndex]);
                currentIndex = smallestIndex;
            } else {
                break;
            }
        }
    }

    // 插入
    void insert(int newElement) {
        heapElements.push_back(newElement);
        siftUp(heapElements.size() - 1);
    }

    int top() const { return heapElements[0]; }

    void pop() {
        int heapSize = heapElements.size();
        if (heapSize == 1) {
            heapElements.pop_back();
            return;
        }
        swap(heapElements[0], heapElements[heapSize - 1]);
        heapElements.pop_back();
        siftDown(0);
    }

    void buildHeap(const vector<int> &inputVector) {
        heapElements = inputVector;
        int heapSize = heapElements.size();
        // 从最后一个非叶子节点开始向上调用 siftDown
        for (int i = (heapSize / 2) - 1; i >= 0; i--) {
            siftDown(i);
        }
    }
};

// 堆排序
void heapSort(vector<int> &vec) {
    MinHeap minHeap(vec);

    int index = 0;
    while (!minHeap.heapElements.empty()) {
        int minElement = minHeap.top();
        vec[index++] = minElement;
        minHeap.pop();
    }
}

// 希尔排序
void shellSortDonald(vector<int> &vec) {
    int n = vec.size();
    vector<int> recordGaps;
    int gap = n / 2;

    // gap
    while (gap > 0) {
        recordGaps.push_back(gap);
        gap = gap / 2;
    }

    for (int i = 0; i < recordGaps.size(); i++) {
        gap = recordGaps[i];
        for (int j = gap; j < n; j++) {
            int temp = j;
            while (temp >= gap) {
                compareCount++;
                if (vec[temp - gap] > vec[temp]) {
                    swapCount++;
                    swap(vec[temp], vec[temp - gap]);
                    temp = temp - gap;
                } else {
                    break;
                }
            }
        }
    };
}

void shellSortHibbard(vector<int> &vec) {
    int n = vec.size();
    vector<int> recordGaps;
    int gap = 1;
    
    // gap
    while (gap <= n) {
        recordGaps.push_back(gap);
        gap = 2 * gap + 1;
    }
    reverse(recordGaps.begin(), recordGaps.end());
    for (int i = 0; i < recordGaps.size(); i++) {
        gap = recordGaps[i];
        for (int j = gap; j < n; j++) {
            int temp = j;
            while (temp >= gap) {
                compareCount++;
                if (vec[temp - gap] > vec[temp]) {
                    swapCount++;
                    swap(vec[temp], vec[temp - gap]);
                    temp = temp - gap;
                } else {
                    break;
                }
            }
        }
    };
}

void shellSortSedgewick(vector<int> &vec) {
    int n = vec.size();
    vector<int> recordGaps;
    recordGaps.push_back(1);
    int i = 0;

    // gap
    while (true) {
        int gap = 4 * pow(2, i) - 3 * pow(2, i / 2) + 1;
        if (gap >= n)
            break;
        recordGaps.push_back(gap);
        i++;
    }
    reverse(recordGaps.begin(), recordGaps.end());
    for (int i = 0; i < recordGaps.size(); i++) {
        int gap = recordGaps[i];
        for (int j = gap; j < n; j++) {
            int temp = j;
            while (temp >= gap) {
                compareCount++;
                if (vec[temp - gap] > vec[temp]) {
                    swapCount++;
                    swap(vec[temp], vec[temp - gap]);
                    temp = temp - gap;
                } else {
                    break;
                }
            }
        }
    };
}

// 桶排序
void bucketSort(vector<int> &vec, int numBuckets = 10) {
    if (vec.empty())
        return;

    int minValue = *min_element(vec.begin(), vec.end());
    int maxValue = *max_element(vec.begin(), vec.end());

    vector<vector<int>> buckets(numBuckets);

    for (auto num : vec) {
        int bucketIndex =
            int(((num - minValue) / (maxValue - minValue) * numBuckets));
        if (bucketIndex == numBuckets)
            bucketIndex--;
        buckets[bucketIndex].push_back(num);
    }

    for (auto &bucket : buckets) {
        mergeSort(bucket, 0, bucket.size() - 1);
    }

    int index = 0;
    for (auto &bucket : buckets) {
        for (auto num : bucket) {
            vec[index++] = num;
        }
    }
}

// 计数排序
void countingSortWithoutCounting(vector<int> &vec) {
    if (vec.empty())
        return;

    int minValue = *min_element(vec.begin(), vec.end());
    int maxValue = *max_element(vec.begin(), vec.end());

    // 频率数组
    vector<int> frequency(maxValue - minValue + 1, 0); 

    for (int num : vec) {
        frequency[num - minValue]++;
    }

    int insertPosition = 0;

    for (int countIndex = 0; countIndex < frequency.size(); countIndex++) {
        while (frequency[countIndex] > 0) {
            vec[insertPosition++] = countIndex + minValue;
            frequency[countIndex]--;
        }
    }
}

// 规定计数排序的 hashFuction
int hashFuctionInRadix(int num, int exp) { return (num / exp) % 10; }
int hashFuctionInCounting(int num, int exp) { return num; }

// 计数排序
void countingSortWithCounting(vector<int> &vec, int exp) {
    int n = vec.size();

    if (n <= 1)
        return;

    int maxVal = *max_element(vec.begin(), vec.end());
    int minVal = *min_element(vec.begin(), vec.end());

    vector<int> count(maxVal - minVal + 1, 0);

    for (int i = 0; i < n; i++) {
        int temp = hashFuctionInCounting(vec[i] - minVal, exp);
        count[temp]++;
    }

    for (int i = 1; i < count.size(); i++) {
        count[i] += count[i - 1];
    }

    // 直接按照元素的顺序获取对应的位置，
    // 从后往前进行索引 从而保证元素的相对位置不发生改变，保证了稳定性
    vector<int> output(n);
    for (int i = n - 1; i >= 0; i--) {
        // 通过hash获取元素位置
        int temp = hashFuctionInCounting(vec[i] - minVal, exp);
        output[count[temp] - 1] = vec[i];
        count[vec[i] - minVal]--;
    }

    for (int i = 0; i < n; i++) {
        vec[i] = output[i];
    }
}

// 基数排序
void radixSort(vector<int> &vec) {
    int n = vec.size();

    if (n <= 1)
        return;

    int maxVal = *max_element(vec.begin(), vec.end());

    // 逐位排序
    for (int exp = 1; maxVal / exp > 0; exp *= 10) {
        countingSortWithCounting(vec, exp);
    }
}

bool isSorted(const vector<int> &vec) {
    for (size_t i = 1; i < vec.size(); ++i) {
        if (vec[i] < vec[i - 1]) {
            return false;
        }
    }
    return true;
}


void testSortingAlgorithm(function<void(vector<int>&)> sortFunc, const string& algoName, vector<int>& vec) {
    vector<int> vecCopy = vec;  
    resetCounters();
    sortFunc(vecCopy);  

    cout << algoName << ": 比较次数 = " << compareCount
         << ", 交换次数 = " << swapCount;
    if (isSorted(vecCopy))
        cout << "，排序成功。" << endl;
    else
        cout << "，排序失败。" << endl;
}

void testAllAlgorithms(vector<int>& vec, const string& filename) {
    testSortingAlgorithm([](vector<int>& vecCopy) { bubbleSort(vecCopy); }, "冒泡排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { optimizedBubbleSort(vecCopy); }, "优化冒泡排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { cocktailShakerSort(vecCopy); }, "双向冒泡排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { selectionSort(vecCopy); }, "选择排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { insertionSort(vecCopy); }, "插入排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { insertionSortWithBinarySearch(vecCopy); }, "二叉搜索插入排序", vec);

    testSortingAlgorithm([](vector<int>& vecCopy) { quickSort(vecCopy, 0, vecCopy.size() - 1); }, "快速排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { mergeSort(vecCopy, 0, vecCopy.size() - 1); }, "归并排序", vec);

    testSortingAlgorithm([](vector<int>& vecCopy) { heapSort(vecCopy); }, "堆排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { shellSortDonald(vecCopy); }, "希尔排序 (Donald)", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { shellSortHibbard(vecCopy); }, "希尔排序 (Hibbard)", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { shellSortSedgewick(vecCopy); }, "希尔排序 (Sedgewick)", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { bucketSort(vecCopy); }, "桶排序", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { countingSortWithoutCounting(vecCopy); }, "计数排序 (不使用计数器)", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { countingSortWithCounting(vecCopy, 1); }, "计数排序 (使用计数器)", vec);
    testSortingAlgorithm([](vector<int>& vecCopy) { radixSort(vecCopy); }, "基数排序", vec);
}
// 主函数
int main() {
    string inputFilename = "random_data.txt"; 
    string outFilename = "sorted_data.txt";

    vector<int> vec;

    generateRandomData(inputFilename, 10000);
    readDataFromFile(inputFilename, vec);
    testAllAlgorithms(vec, outFilename);

    return 0;
}
