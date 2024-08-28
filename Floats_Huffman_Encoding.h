#if !defined(FLOATS_HUFFMAN_ENCODING_H)
#define FLOATS_HUFFMAN_ENCODING_H

#include <stdio.h>
#include <queue>     
#include <stdbool.h>

#define MAXLENGTH 128 * 1024 * 1024


/* struct that represents a letter */
struct node{
   int index;
   int frq;
   struct node *left;
   struct node *right;
   unsigned int encode;
   int encode_length;
 

   bool operator<(const struct node& comp_node) const {
        if (frq > comp_node.frq) {
            return true;
        }
        return false;
    }
};
typedef struct node Node;

int quantize(float* src, size_t size, float error, int* quantized_src, float* min);
void record_frequencies(const int* quantized_src, int fsize, int num_buckets, int* frqs);
void make_queue(std::priority_queue<Node> &tree, int *frqs, int num_buckets);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode(Node *root);
void store_encodings_helper(Node* root, int encodings[][2], int num_buckets);
void compress(const int* src, char* dest, int fsize, int* destsize, int encodings[][2]);
void reverse_quantize(const int* quantized_data, int size, float min, float bucket_size, float* data);
void decompress(const char* src, float* dest, Node* root, int num_elements, float min, float bucket_size);
double calc_speed(long original_size, double compression_time);
void free_tree(Node* root);
#endif
