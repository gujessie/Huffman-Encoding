#if !defined(FLOATS_HUFFMAN_ENCODING_H)
#define FLOATS_HUFFMAN_ENCODING_H

#include <stdio.h>
#include <queue> 
#include <vector>    
#include <stdbool.h>

#define MAXLENGTH 128 * 1024 * 1024


/* struct that represents a letter */
struct Node {
    int index; // The bucket or value this node represents
    int freq;  // Frequency of the bucket or value
    Node* left;  // Left child in the Huffman tree
    Node* right; // Right child in the Huffman tree
    unsigned int encode;
    int encode_length;

    Node(int idx, int frequency) : index(idx), freq(frequency), left(NULL), right(NULL) {}
};


struct LessThanByCnt {
    bool operator()(const Node* lhs, const Node* rhs) const {
        return lhs->freq > rhs->freq;  // Min-heap based on frequency
    }
};

int quantize(float* src, size_t size, float error, int* quantized_src, float* min);
void record_frequencies(const int* quantized_src, int fsize, int num_buckets, int* frqs);
void make_queue(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& phtree, int* frqs, int num_buckets);
Node* build_tree(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& tree);
void assign_encode(Node *root);
void store_encodings(Node* root, int encodings[][2], int num_buckets);
void store_encodings_helper(Node* root, int encodings[][2], int num_buckets);
void compress(const int* src, char* dest, int fsize, int* destsize, int encodings[][2]);
void reverse_quantize(const int* quantized_data, int size, float min, float bucket_size, float* data);
void decompress(const char* src, float* dest, Node* root, int num_elements, float min, float bucket_size);
double calc_speed(long original_size, double compression_time);
void free_tree(Node* root);
#endif
