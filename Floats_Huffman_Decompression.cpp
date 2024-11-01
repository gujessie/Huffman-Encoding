#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include "Floats_Huffman_Encoding.h"

//seperate compression and decompression into seperate files

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <vector>
#include "Floats_Huffman_Encoding.h"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "error in arguments");
        return 1;
    }

    // Opening input file
    FILE *input = fopen(argv[1], "rb");
    if (NULL == input) { 
        perror("Error opening input file");
        return 1;
    }

    // figure out the size of compressed file.
   // Calculate compressed data size, read to end of file
    fseek(input, 0, SEEK_END); // Move to the end of the file
    long file_size = ftell(input); // Get current position (file size)
    printf("file size = %lu\n", file_size);    

    char * compressed_buffer = (char *) malloc(file_size);
    if (NULL == compressed_buffer) {
        perror("Error allocating buffer");
        return 1;
    }

    dest = (float *) malloc (file_size);
    if (NULL == dest) {
        perror("Error allocating buffer");
        return 1;
    }
    decompress(compressed_buffer, dest, file_size, int * destsize){

#if 0
    // Read: min, bucket_size, fsize, and frequencies
    float min, bucket_size;
    size_t num_elements;
    fread(&num_elements, sizeof(size_t), 1, input);
    fread(&min, sizeof(float), 1, input);
    fread(&bucket_size, sizeof(float), 1, input);
    int frq_count; //elements in frq array
    fread(&frq_count, sizeof(int), 1, input); //num elements in frqs, aka number of buckets

    printf("num_elements: %d\n", num_elements);
    printf ("min = %f, bucket_size = %f\n", min, bucket_size);

    // Read frequency array from the file
    int *frqs = (int *)malloc(frq_count * sizeof(int));
    fread(frqs, sizeof(int), frq_count, input);

    printf("frq_count = %d\n", frq_count);

    printf ("before make queue\n");
    // Priority queue to rebuild the Huffman tree
    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
    make_queue(tree, frqs, frq_count);

    printf ("before build tree\n");

    Node *root = build_tree(tree);

    printf ("after build tree\n");

    // Calculate compressed data size, read to end of file
    fseek(input, 0, SEEK_END); // Move to the end of the file
    long file_size = ftell(input); // Get current position (file size)
    printf("file size = %lu\n", file_size);
    //Move to where compression starts

    fseek(input, frq_count * sizeof(int) 
                + sizeof(num_elements) 
                + sizeof(min) 
                + sizeof(bucket_size) 
                + sizeof(frq_count), 
                SEEK_SET); 

    
    //sizeof(size_t) + sizeof(float) + sizeof(float) + sizeof(int) + frq_count * sizeof(int)

    printf ("Header size = %d\n", frq_count * sizeof(int) 
                + sizeof(num_elements) 
                + sizeof(min) 
                + sizeof(bucket_size) 
                + sizeof(frq_count));

    long compressed_size = file_size - ftell(input);  // compressed_size is remaining size

printf("compressed_size = %ld\n", compressed_size);
    // Read compressed data from file
    char *compressed_data = (char *)malloc(compressed_size);  // Adjust size as needed
    fread(compressed_data, 1, compressed_size, input);
    fclose(input);

    // Decompress the data
    float *dest = (float *)malloc(num_elements * sizeof(float));

    if (NULL == dest) {
        perror("Error in malloc for decomp data\n");
        free(frqs);
        free(compressed_data);
        free_tree(root);
        return 1;
    }

    printf("Jessie\n");
    
    decompress(compressed_data, dest, root, num_elements, min, bucket_size);

    printf("Jessie 2\n");

    // Write decompressed data to output file
    FILE *output = fopen(argv[2], "wb");
    if (output == NULL) {
        perror("Error opening output file");
        free(frqs);
        free(compressed_data);
        free(dest);
        free_tree(root);
        return 1;
    }
    fwrite(dest, sizeof(float), num_elements, output);
    fclose(output);

    // Free memory
    free(frqs);
    free(compressed_data);
    free(dest);
    free_tree(root);
#endif
    return 0;
}


int quantize(float* src, size_t size, float error, int* quantized_src, float* min) {
    // Find min and max values in data
    *min = src[0];
    float max = src[0];
    for (int i = 1; i < size; i++) {
        if (src[i] < *min) *min = src[i];
        if (src[i] > max) max = src[i];
    }

    // calc the interval and bucket width
    float interval = max - *min;
    int num_buckets = (int)(interval/(error * 2) + 1);
    float bucket_size = interval / num_buckets;

    // Assign each data point to a bucket
    for (int i = 0; i < size; i++) {
        int bucket = (int)((src[i] - *min) / bucket_size);
        if (bucket == num_buckets) {
            bucket--;  // 
        }

        // Store the bucket index
        quantized_src[i] = bucket;
    }

    return num_buckets;
}

void record_frequencies(const int* quantized_src, int fsize, int num_buckets, int* frqs) {

    // Count frqs
    for (int i = 0; i < fsize; i++) {
        int bucket = quantized_src[i];
        if (bucket >= 0 && bucket < num_buckets) {
            frqs[bucket]++; // Increment the frequency for the bucket
        }
    }
}

//inserts nodes that have a frequency(exists in the data) into the priority queue tree
void make_queue(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& phtree, int* frqs, int num_buckets) {
    // Iterate through all buckets
    for (int i = 0; i < num_buckets; ++i) {
        if (frqs[i] != 0) {
            // Create a new node with the bucket index and its frequency
            Node* new_node = new Node(i, frqs[i]);
            // Push the pointer to the node into the priority queue
            phtree.push(new_node);
        }
    }
}

//creates huffman tree
Node* build_tree(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& tree) {
    // Continue until there is only one node left in the priority queue
    while (tree.size() > 1) {
        // Extract the two nodes with the smallest frequencies
        Node* left = tree.top();
        tree.pop();
        Node* right = tree.top();
        tree.pop();

        // Create a new parent node with a combined frequency
        Node* parent = new Node(-1, left->freq + right->freq);
        parent->left = left;
        parent->right = right;

        // Push the new parent node back into the priority queue
        tree.push(parent);
    }

    // The last remaining node is the root of the Huffman tree
    return tree.top();
}

//assigns the binary code to each node
void assign_encode_helper(Node *node, unsigned int encode, int length) {
    if (NULL == node){
        return;
    }

    
    // If the node is a leaf node (no left or right children)
    if (node->left == NULL && node->right == NULL) {
        // Assign the encoding and its length to the leaf node
        node->encode = encode;
        node->encode_length = length;
    } else {
        // Traverse the left subtree; append 0 to the encode
        assign_encode_helper(node->left, (encode << 1), length + 1);

        // Traverse the right subtree; append 1 to the encode
        assign_encode_helper(node->right, (encode << 1) | 1, length + 1);
    }
}

void assign_encode(Node *root) {
     if (root != NULL) {
        assign_encode_helper(root, 0, 0);
    }
}

void store_encodings_helper(Node* root, int encodings[][2], int num_buckets) {
    if (!root) return;

    if (root->index >= 0) {
        encodings[root->index][0] = root->encode;
        encodings[root->index][1] = root->encode_length;
    }

    if (root->left) {
        store_encodings_helper(root->left, encodings, num_buckets);
    }

    if (root->right) {
        store_encodings_helper(root->right, encodings, num_buckets);
    }
}

void store_encodings(Node* root, int encodings[][2], int num_buckets) {
    for (int i = 0; i < num_buckets; ++i) {
        encodings[i][0] = 0;   // Initialize encode
        encodings[i][1] = 0;   // Initialize encode length
    }
    store_encodings_helper(root, encodings, num_buckets);
}

//STORE ELEMENTS OF DECOMPRESS IN THE BEGINNING OF FILE
// ZSTD: size_t ZSTD_decompress( void* dst, size_t dstCapacity,
//                  const void* src, size_t compressedSize);

// src here points to a buffer that has metadata (min, frqs...) and compressed data.
// void decompress(const char * fname, ....)
void decompress(const char* src, float* dest, size_t src_size, int * destsize){
    float min;
    // inside here, the code figures out min, frqs, building a tree.

    // Not like this:

 
}


void decompress_1(const char* src, float* dest, Node* root, int num_elements, float min, float bucket_size) {
    Node* current = root;
    unsigned int buffer = 0;
    int bits_in_buffer = 0;
    int src_index = 0;
    int dest_index = 0;

    while (dest_index < num_elements) {
 
        // Load bits from the source buffer into the buffer
        unsigned char byte = src[src_index++];
        buffer = (buffer << 8) | byte;
        bits_in_buffer += 8;

        // Decode bits from the buffer into original values
        while (bits_in_buffer > 0) {
            int bit = (buffer >> (bits_in_buffer - 1)) & 1;
            bits_in_buffer--;

            if (bit == 0) {
                current = current->left;
            } else {
                current = current->right;
            }

            // If we've reached a leaf node, decode the quantized value
            if (current->left == NULL && current->right == NULL) {
                int bucket_index = current->index;
                // Convert bucket index back to the quantized float value
                float quantized_value = min + bucket_index * bucket_size + bucket_size / 2.0;
                dest[dest_index++] = quantized_value; // Store the reconstructed float value
                current = root; // Reset to the root for the next decoding
            }
        }
    }
}

double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / compression_time;
}

void free_tree(Node* root) { 
    if (!root) return;

    // Recursively free left and right subtrees
    if (root->left){
        free_tree(root->left);
    }
    
    if (root->right){
        free_tree(root->right);
    }

    // Free the current node
    free(root);
}



void reverse_quantize(const int* quantized_data, int size, float min, float bucket_size, float* data) {
    for (int i = 0; i < size; i++) {
        // Calculate the approximate original value by finding the midpoint of the bucket
        data[i] = min + quantized_data[i] * bucket_size + (bucket_size / 2.0);
    }
}
