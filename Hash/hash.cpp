#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

const int g = 131;        // 基数
const int p = 1000000007; // 素数

// 哈希表类：链式哈希
class HashTableChaining {
public:
    HashTableChaining(int size) : size(size) {
        table.resize(size); // 初始化哈希表，大小为 size
    }

    int hash(int key) { return key % size; } // 哈希函数

    void insert(int key) {
        int index = hash(key);
        table[index].push_back(key); // 插入元素
    }

    int search(int key) {
        int index = hash(key);
        int steps = 0;
        for (int item : table[index]) {
            steps++;
            if (item == key)
                return steps;
        }
        return 0;
    }

    double average_search_length() {
        int total_search_length = 0, total_elements = 0;
        for (const auto &bucket : table) {
            for (int item : bucket) {
                total_elements++;
                total_search_length += search(item);
            }
        }
        return total_elements == 0
                   ? 0
                   : (double)total_search_length / total_elements;
    }

private:
    int size;
    vector<list<int>> table;
};

// 哈希表类：开放定址哈希
class HashTableOpenAddressing {
public:
    HashTableOpenAddressing(int size) : size(size) {
        table.resize(size, -1); // 空槽用 -1 表示
    }

    int hash(int key) { return key % size; }

    void insert(int key) {
        int index = hash(key);
        int original_index = index;

        while (table[index] != -1) {
            index = (index + 1) % size; // 线性探查
            if (index == original_index) {
                cout << "哈希表已满" << endl;
                return;
            }
        }
        table[index] = key; // 找到空槽，插入
    }

    int search(int key) {
        int index = hash(key);
        int original_index = index;
        int steps = 0;

        while (table[index] != -1) {
            steps++;
            if (table[index] == key)
                return steps;
            index = (index + 1) % size;
            if (index == original_index)
                break;
        }
        return 0;
    }

    double average_search_length() {
        int total_search_length = 0, total_elements = 0;
        for (int key : table) {
            if (key != -1) {
                total_elements++;
                total_search_length += search(key);
            }
        }
        return total_elements == 0
                   ? 0
                   : (double)total_search_length / total_elements;
    }

private:
    int size;
    vector<int> table;
};

long long poly_hash(const string &s) {
    long long hash_value = 0;
    long long power = 1;
    for (int i = s.size() - 1; i >= 0; --i) {
        hash_value = (hash_value + (s[i] * power) % p) % p;
        power = (power * g) % p;
    }
    return hash_value;
}

bool is_vowel(char c) {
    c = tolower(c);
    return (c == 'a' || c == 'e' || c == 'i' || c == 'o' || c == 'u');
}

vector<string> split_into_syllables_by_rules(const string &word) {
    vector<string> syllables;
    string current_syllable;

    for (int i = 0; i < word.length(); ++i) {
        current_syllable += word[i];

        if (is_vowel(word[i])) {
            if (i + 1 < word.length() && !is_vowel(word[i + 1]) ||
                i == word.length() - 1) {
                syllables.push_back(current_syllable);
                current_syllable = "";
            }
        }
    }

    return syllables;
}

long long get_combined_hash_with_add(const string &word) {
    vector<string> syllables = split_into_syllables_by_rules(word); 
    long long combined_hash = 0;

    for (const string &syllable : syllables) {
        long long syllable_hash = poly_hash(syllable); 
        combined_hash =
            (combined_hash + syllable_hash) % p; 
    }

    return combined_hash; 
}

long long get_combined_hash_with_poly(const string &word) {
    vector<string> syllables = split_into_syllables_by_rules(word); 
    long long combined_hash = 0;

    int power = 1;

    for (int i = 0; i < syllables.size(); i++) {
        string syllable = syllables[i];
        long long syllable_hash = poly_hash(syllable);

        combined_hash = (combined_hash + (syllable_hash * power) % p) % p;

        power = (power * g) % p;
    }

    return combined_hash;
}

// 对于链式哈希的构建
void build_hash_table(HashTableChaining &hash_table, const string &article,
                      int key) {
    stringstream ss(article);
    string word;
    int temp;
    while (ss >> word) {
        if (key == 0) {
            temp = get_combined_hash_with_add(word);
        } else if (key == 1) {
            temp = get_combined_hash_with_poly(word);
        } else {
            temp = poly_hash(word);
        }
        hash_table.insert(temp);
    }
}

// 对于开放定址哈希的构建
void build_hash_table(HashTableOpenAddressing &hash_table,
                      const string &article, int key) {
    stringstream ss(article);
    string word;
    int temp;
    while (ss >> word) {
        if (key == 0) {
            temp = get_combined_hash_with_add(word);
        } else if (key == 1) {
            temp = get_combined_hash_with_poly(word);
        } else {
            temp = poly_hash(word);
        }
        hash_table.insert(temp);
    }
}

// 计算链式哈希的平均查找长度
double calculate_average_search_length(HashTableChaining &hash_table,
                                       const string &article, int key) {
    stringstream ss(article);
    string word;
    double total_search_length = 0;
    int total_elements = 0;
    while (ss >> word) {
        if (key == 0) {
            key = get_combined_hash_with_add(word);
        } else if (key == 1) {
            key = get_combined_hash_with_poly(word);
        } else {
            key = poly_hash(word);
        }
        hash_table.insert(key);
        total_search_length += hash_table.search(key);
        total_elements++;
    }
    return total_elements == 0 ? 0 : total_search_length / total_elements;
}

// 计算开放定址哈希的平均查找长度
double calculate_average_search_length(HashTableOpenAddressing &hash_table,
                                       const string &article,
                                       bool use_poly_hash) {
    stringstream ss(article);
    string word;
    double total_search_length = 0;
    int total_elements = 0;
    while (ss >> word) {
        int key = use_poly_hash ? poly_hash(word) : stoi(word);
        total_search_length += hash_table.search(key);
        total_elements++;
    }
    return total_elements == 0 ? 0 : total_search_length / total_elements;
}

string read_file(const string &filename) {
    ifstream file(filename); 
    if (!file.is_open()) {  
        cerr << "无法打开文件: " << filename << endl;
        return "";
    }

    string content;
    string line;
    while (getline(file, line)) { 
        content += line + "\n";  
    }

    file.close(); 
    return content;
}


int main() {
    string article1 = read_file("1.txt");
    string article2 = read_file("2.txt");

  // 1. 链式哈希 + 音节哈希
    HashTableChaining hash_table_chaining_syllable_add(100);
    build_hash_table(hash_table_chaining_syllable_add, article1, 0); 
    double avg_search_length = hash_table_chaining_syllable_add.average_search_length();
    cout << "Average Search Length (Chaining + Syllable add Hash): " << avg_search_length << endl;

    // 2. 链式哈希 + 音节多项式哈希
    HashTableChaining hash_table_chaining_syllable_poly(100);
    build_hash_table(hash_table_chaining_syllable_poly, article1, 1); 
    double avg_search_length_ = hash_table_chaining_syllable_poly.average_search_length();
    cout << "Average Search Length (Chaining + Syllable poly Hash): " << avg_search_length_ << endl;

    // 3. 链式哈希 + 多项式哈希
    HashTableChaining hash_table_chaining_poly(100);
    build_hash_table(hash_table_chaining_poly, article1, 2); 
    double avg_search_length_poly = hash_table_chaining_poly.average_search_length();
    cout << "Average Search Length (Chaining + Poly Hash): " << avg_search_length_poly << endl;

    // 4. 开放定址哈希 + 音节哈希
    HashTableOpenAddressing hash_table_open_syllable_add(200);
    build_hash_table(hash_table_open_syllable_add, article1, 0); 
    double avg_search_length_open = hash_table_open_syllable_add.average_search_length();
    cout << "Average Search Length (Open Addressing + Syllable add Hash): " << avg_search_length_open << endl;

    // 5. 开放定址哈希 + 音节多项式哈希
    HashTableOpenAddressing hash_table_open_syllable_poly(200);
    build_hash_table(hash_table_open_syllable_poly, article1, 1); 
    double avg_search_length_open_ = hash_table_open_syllable_poly.average_search_length();
    cout << "Average Search Length (Open Addressing + Syllable poly Hash): " << avg_search_length_open_ << endl;

    // 6. 开放定址哈希 + 多项式哈希
    HashTableOpenAddressing hash_table_open_poly(200);
    build_hash_table(hash_table_open_poly, article1, 2); 
    double avg_search_length_open_poly = hash_table_open_poly.average_search_length();
    cout << "Average Search Length (Open Addressing + Poly Hash): " << avg_search_length_open_poly << endl;

    return 0;
}
