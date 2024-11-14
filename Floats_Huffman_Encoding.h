#ifndef FLOATS_HUFFMAN_ENCODING_H
#define FLOATS_HUFFMAN_ENCODING_H

#include <stdio.h>
#include <queue>
#include <vector>

// Removed <stdbool.h> as C++ has a built-in bool type

#define MAXLENGTH (128 * 1024 * 1024) // Added parentheses for safety

/* Struct that represents a node in the Huffman tree */
struct Node {
    int index;              // The bucket or value this node represents
    int freq;               // Frequency of the bucket or value
    Node* left;             // Left child in the Huffman tree
    Node* right;            // Right child in the Huffman tree
    unsigned int encode;    // Huffman encoding
    int encode_length;      // Length of the Huffman encoding

    Node(int idx, int frequency) 
        : index(idx), freq(frequency), left(nullptr), right(nullptr), encode(0), encode_length(0) {}
};

/* Comparator for the priority queue (min-heap based on frequency) */
struct LessThanByCnt {
    bool operator()(const Node* lhs, const Node* rhs) const {
        return lhs->freq > rhs->freq; // Min-heap based on frequency
    }
};

// Function Declarations
int quantize(float* src, size_t size, float error, int* quantized_src, float* min);
void record_frequencies(const int* quantized_src, size_t fsize, size_t num_buckets, int* frqs);
void make_queue(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& phtree, int* frqs, size_t num_buckets);
Node* build_tree(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& tree);
void assign_encode(Node* root);
void store_encodings(Node* root, int* encodings, size_t num_buckets);
void store_encodings_helper(Node* root, int* encodings, size_t num_buckets);
void compress_quantization_levels(const int* src, char* dest, size_t num_elements, int* destsize, int* encodings);
double calc_speed(long original_size, double compression_time);
void free_tree(Node* root);
int compress(float* src, char* dest, size_t num_elements, int* destsize, float error);
void writefile(const char* fname, const char* buff, size_t size);
int decompress(const char* src, float* dest, size_t src_size, int* destsize);
//void huffman_decompress(const char* src, float* dest, Node* root, int num_elements, float min, float bucket_size);


#endif // FLOATS_HUFFMAN_ENCODING_H
