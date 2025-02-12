# 排序算法

**感悟：先实现内部的具体逻辑，在从外部进行串线，先实现核心逻辑，在进行相关的检测**

1. 冒泡排序
    1. 优化1：增加一个ifSwapped的优化，在没有更新的循环中，直接退出循环
    2. 优化2：对数据进行前后排序，交替地从数组的两端逐步将最大元素和最小元素“冒泡”
2. 选择排序
    1. 不稳定性——找到最小元素与当前元素交换
        1. [{3, 'A'}, {2, 'B'}, {3, 'C'}, {1, 'D'}]
        2. 第一次选择排序——[{1, 'D'}, {2, 'B'}, {3, 'C'}, {3, 'A'}]
        3. AC位置互换，故不稳定
3. 插入排序
    1. 从未排序部分中取出元素，将取出的元素与已排序部分的元素从后向前比较，找到插入位置，排序部分中，比取出元素大的元素逐一向后移动，为插入的元素腾出空间
    2. 改进：将增加一个位置的二分查找，减少比较的次数（无法减少交换的次数）
4. 快速排序
    1. 快速排序分两个函数，partition和quickSort
    2. 在partition中实现数据元素的交换，注意维护一个index，注意循环中的参数为left，right，而不是0和size
    3. 在quickSort中实现递归进行partition
5. 归并排序
    1. 归并排序同样分两个函数merge和mergeSort
    2. 在mergeSort中实现对原先数组进行拆分，将left，right， mid传入到merge之中
    3. 在merge之中，分进行数组的拆分和组合，实现排序
    4. 无需进行数据的交换，数据重排序
6. 堆排序
    1. 堆内部用一个数组来存贮树的结构，在其中，可以存储数据，取得其中的最大值和最小值
    2. 在实现上面，主要有两个重要点——上浮和下沉
        1. 上浮过程只需要和parent节点进行对比，如果大，进行更换
        2. 下沉通过递归实现，对比左右子节点，如果根节点大于子节点，进行更换。直到无法找到左右节点进行更换的时候或者进入叶节点的时候，退出递归
    3. 通过堆，每次获取最上面的节点即可，在删除最上面的节点时，将最上面的节点和最后一个节点进行交换，然后删除最后一个节点，对最上面的节点进行下沉操作即可
7. 希尔排序
    1. 实际上希尔排序是对插入排序的优化，在通过利用变量gap，来将原先的数组分成若干组，然后对这些组进行插入排序的操作
    2. 优化—使用不同的gap
        1. $counts/2$
        2. $a_n = 2a_{n-1} + 1$
        3. $a_n = 9\cdot 4^n - 9\cdot 2^n +1$
8. 计数排序
    1. 类似于图的方法，将数据根据最大值和最小值建立映射
    2. 在处理上面，对于仅元素排序，可以不用counting数组，直接将对应遍历元素的位置进行输出即可
    3. 对于更广泛的情况，增加counting数组，用于获取每个元素的位置，然后从后往前进行索引，获取位置进行输出，在其过程中，可以使用一个hash函数对数据进行映射的排序
    4. 增加counting后，排序具有稳定性
9. 基数排序
    1. 对每一位进行排序，最后获取全部的排序结果
    2. 实现可以基于计数排序，将hash函数设置为`(num / exp) % 10` 即可
10.  桶排序
    1. 将数据分到不同的桶之中，然后在桶的内部直接排序即可
    2. 分桶的过程需要首先获取最大值和最小值，然后设定一个和桶个数的超参数，最后确定每个通的范围即可

## API解释

1. 随机数据生成
- `generateRandomData(filename, n)`：生成 n 个随机整数并保存到指定文件中
2. 数据读取
- `readDataFromFile(filename, vec)`：从指定文件中读取数据
3. 排序算法实现
实现了以下排序算法：
- 冒泡排序 (`bubbleSort`)
- 优化冒泡排序 (`optimizedBubbleSort`)：加入提前退出机制的冒泡排序，减少不必要的遍历
- 双向冒泡排序 (`cocktailShakerSort`)：从两端向中间遍历，优化排序过程
- 选择排序 (`selectionSort`)：找到未排序部分的最小值并交换
- 插入排序 (`insertionSort`)：将元素插入到已排序的部分
- 二叉搜索插入排序 (`insertionSortWithBinarySearch`)：通过二叉搜索优化插入排序
- 快速排序 (`quickSort`)：找出基准元素，分区然后递归快排
- 归并排序 (`mergeSort`)：分割数组，然后合并
- 堆排序 (`MinHeap`)：实现了最小堆的类，用于排序
- 希尔排序：
  - `shellSortDonald`：使用 Donald 的间隔序列
  - `shellSortHibbard`：使用 Hibbard 的间隔序列
  - `shellSortSedgewick`：使用 Sedgewick 的间隔序列
- 桶排序 (`bucketSort`)：将元素分配到不同的桶中，然后对每个桶进行排序，最后将桶中的元素合并
- 计数排序：
  - `countingSortWithoutCounting`：不使用额外计数器的计数排序
  - `countingSortWithCounting`：基于计数排序的变体，采用哈希函数
- 基数排序 (`radixSort`)：通过逐位进行排序来实现排序，从最低位到最高位排序
4. 性能计数和结果分析
- `compareCount`：算法在执行时进行的元素比较操作的次数
- `swapCount`：算法在执行时交换元素的次数
- `testAllAlgorithms` ：检测算法


## 排序算法结果记录

注：使用10000个随机数，记录比较次数和交换次数

---

| 排序算法                | 比较次数    | 交换次数   | 平均时间复杂度          | 空间复杂度        |
| ----------------------- | ----------- | ---------- | ----------------------- | ------------------ |
| **冒泡排序**            | 49995000    | 25105229   | O(n^2)                 | O(1)               |
| **优化冒泡排序**        | 49984989    | 25105229   | O(n^2)                 | O(1)               |
| **双向冒泡排序**        | 37597290    | 25105229   | O(n^2)                 | O(1)               |
| **选择排序**            | 49995000    | 9987       | O(n^2)                 | O(1)               |
| **插入排序**            | 49995000    | 25105229   | O(n^2)                 | O(1)               |
| **二叉搜索插入排序**    | 59974       | 25111429   | O(n^2)                 | O(n)               |
| **快速排序**            | 172770      | 81213      | O(n log n)             | O(log n)           |
| **归并排序**            | 133616      | 0          | O(n log n)             | O(n)               |
| **堆排序**              | 166455      | 114244     | O(n log n)             | O(1)               |
| **希尔排序 (Donald)**   | 259926      | 145039     | O(n^(3/2))             | O(1)               |
| **希尔排序 (Hibbard)**  | 243869      | 135458     | O(n^(3/2))             | O(1)               |
| **希尔排序 (Sedgewick)**| 204181      | 95252      | O(n^(3/2))             | O(1)               |
| **桶排序**              | 133601      | 0          | O(n + k)               | O(n + k)           |
| **计数排序 (不使用计数器)** | 0         | 0          | O(n + k)               | O(k)               |
| **计数排序 (使用计数器)** | 0         | 0          | O(n + k)               | O(k)               |
| **基数排序**            | 0           | 0          | O(nk)                  | O(n + k)           |

---