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

using namespace std;

// Implement writefile function
void writefile(const char* fname, const char* buff, size_t size) {
    FILE* output = fopen(fname, "wb");
    if (NULL == output) {
        perror("Error opening output file");
        return;
    }

    size_t written = fwrite(buff, 1, size, output);
    if (written != size) {
        printf("Error writing to output file\n");
    }

    fclose(output);
}

// Implement main function
int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    // Opening input file
    FILE* input = fopen(argv[1], "rb"); // Changed to "rb" for binary read
    if (NULL == input) {
        perror("Error opening input file");
        return 1;
    }

    // Get input file size
    fseek(input, 0, SEEK_END);
    size_t num_elements = ftell(input) / sizeof(float); // Number of elements in file
    fseek(input, 0, SEEK_SET);

    printf("num_elements: %ld\n", num_elements);

    // Allocate memory for source data
    float* src = (float*)malloc(num_elements * sizeof(float));
    if (NULL == src) {
        printf("Error allocating memory for src\n");
        fclose(input);
        return 1;
    }

    // Read data from input file
    size_t read_elements = fread(src, sizeof(float), num_elements, input);
    if (read_elements != num_elements) {
        printf("Error reading input file\n");
        fclose(input);
        free(src);
        return 1;
    }
    fclose(input);

    // Allocate memory for compressed data (estimate)
    // Assuming maximum compressed size is roughly 2 bytes per float
    char* dest = (char*)malloc(num_elements * 2);
    if (NULL == dest) {
        perror("Error allocating memory for dest");
        free(src);
        return 1;
    }

    // Define quantization error
    const float error = 0.25f;

    // Variables to hold compression metadata
    float min = 0.0f;
    float bucket_size = 0.0f;
    int frq_count = 0;
    int* frqs = NULL;
    int** encodings = NULL;
    int destsize = 0;

    // Perform compression
    int compress_status = compress(src, dest, num_elements, &destsize, error,
                                   &min, &bucket_size, &frq_count, &frqs, &encodings);
    if (compress_status != 0) {
        printf("Compression failed\n");
        free(src);
        free(dest);
        return 1;
    }


    size_t header_size = sizeof(size_t) + sizeof(float) + sizeof(float) + sizeof(int) +
                         frq_count * sizeof(int) + frq_count * 2 * sizeof(int);

    // Calculate total size: header + compressed data
    size_t total_size = header_size + ((destsize + 7) / 8);

    // Allocate buffer for header and data
    char* file_buffer = (char*)malloc(total_size);
    if (NULL == file_buffer) {
        perror("Error allocating memory for file_buffer");
        free(src);
        free(dest);
        return 1;
    }


    // Serialize header into file_buffer
    size_t offset_ser = 0;
    memcpy(file_buffer + offset_ser, &num_elements, sizeof(size_t));
    offset_ser += sizeof(size_t);
    memcpy(file_buffer + offset_ser, &min, sizeof(float));
    offset_ser += sizeof(float);
    memcpy(file_buffer + offset_ser, &bucket_size, sizeof(float));
    offset_ser += sizeof(float);
    memcpy(file_buffer + offset_ser, &frq_count, sizeof(int));
    offset_ser += sizeof(int);
    memcpy(file_buffer + offset_ser, frqs, frq_count * sizeof(int));
    offset_ser += frq_count * sizeof(int);

    // Serialize encodings (int[frq_count][2])
    // Since encodings are int** pointing to contiguous encodings_block
    // Iterate and copy each pair
    for (int i = 0; i < frq_count; ++i) {
        memcpy(file_buffer + offset_ser, encodings[i], 2 * sizeof(int));
        offset_ser += 2 * sizeof(int);
    }

    // Serialize compressed data
    memcpy(file_buffer + offset_ser, dest, (destsize + 7) / 8);
    offset_ser += (destsize + 7) / 8;

    // Write to output file
    writefile(argv[2], file_buffer, total_size);

    // Free allocated memory
    free(src);
    free(dest);
    free(file_buffer);
    // Free frqs and encodings
    free(frqs);
    free(encodings); // encodings_block was allocated as a single block and pointed by encodings[i]

    return 0;
}

// Implement quantize function
int quantize(float* src, size_t size, float error, int* quantized_src, float* min) {
    // Find min and max values in data
    *min = src[0];
    float max = src[0];
    for (size_t i = 1; i < size; i++) {
        if (src[i] < *min)
            *min = src[i];
        if (src[i] > max)
            max = src[i];
    }

    // Calculate the interval and bucket width
    float interval = max - *min;
    if (interval == 0) {
        // All values are the same
        for (size_t i = 0; i < size; i++) {
            quantized_src[i] = 0;
        }
        return 1; // Only one bucket
    }
    int frq_count = static_cast<int>(interval / (error * 2) + 1);
    float bucket_size = interval / frq_count;

    // Assign each data point to a bucket
    for (size_t i = 0; i < size; i++) {
        int bucket = static_cast<int>((src[i] - *min) / bucket_size);
        if (bucket >= frq_count) {
            bucket = frq_count - 1; // Handle edge case
        }

        // Store the bucket index
        quantized_src[i] = bucket;
    }

    return frq_count;
}


void record_frequencies(const int* quantized_src, size_t fsize, size_t num_buckets, int* frqs) {
    // Count frequencies
    for (size_t i = 0; i < fsize; i++) {
        int bucket = quantized_src[i];
        if (bucket >= 0 && static_cast<size_t>(bucket) < num_buckets) {
            frqs[bucket]++; // Increment the frequency for the bucket
        }
    }
}


void make_queue(std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt>& phtree, int* frqs, size_t num_buckets) {
    // Iterate through all buckets
    for (size_t i = 0; i < num_buckets; ++i) {
        if (frqs[i] != 0) {
            // Create a new node with the bucket index and its frequency
            Node* new_node = new Node(static_cast<int>(i), frqs[i]);
            // Push the pointer to the node into the priority queue
            phtree.push(new_node);
        }
    }
}


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


void assign_encode_helper(Node* node, unsigned int encode, int length) {
    if (NULL == node) {
        return;
    }

    // If the node is a leaf node (no left or right children)
    if (NULL == node->left && NULL == node->right) {
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

// Implement assign_encode function
void assign_encode(Node* root) {
    if (NULL != root) {
        assign_encode_helper(root, 0, 0);
    }
}


void store_encodings_helper(Node* root, int encodings[][2], size_t frq_count) {
    if (!root)
        return;

    if (root->index >= 0 && static_cast<size_t>(root->index) < frq_count) {
        encodings[root->index][0] = root->encode;
        encodings[root->index][1] = root->encode_length;
    }

    if (root->left) {
        store_encodings_helper(root->left, encodings, frq_count);
    }

    if (root->right) {
        store_encodings_helper(root->right, encodings, frq_count);
    }
}


void store_encodings(Node* root, int encodings[][2], size_t frq_count) {
    // Initialize encodings
    for (size_t i = 0; i < frq_count; ++i) {
        encodings[i][0] = 0; // Initialize encode
        encodings[i][1] = 0; // Initialize encode length
    }
    store_encodings_helper(root, encodings, frq_count);
}


void compress_quantization_levels(const int* src, char* dest, size_t num_elements, int* destsize, int encodings[][2]) {
    int dest_index = 0;           // Tracks current position in the destination buffer
    unsigned int buffer = 0;      // Stores bits before writing to the destination buffer
    size_t src_index = 0;         // Tracks current position in the source buffer
    int bits_in_buffer = 0;       // Counts bits currently in the buffer
    *destsize = 0;                 // Size of the destination in bits

    // Iterate over each element in the source array
    while (src_index < num_elements) {
        int bucket_index = src[src_index];
        int encode = encodings[bucket_index][0];    // Get Huffman encoding for the value
        int length = encodings[bucket_index][1];    // Get the length of the encoding

        // Add the current value's encoding to the bit buffer
        buffer = (buffer << length) | encode;
        bits_in_buffer += length;

        // Write to the destination buffer byte by byte until there are less than 8 bits left
        while (bits_in_buffer >= 8) {
            unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
            dest[dest_index++] = byte;    // Add the byte to the destination buffer
            bits_in_buffer -= 8;          // Update the number of bits in the buffer
            *destsize += 8;
        }
        src_index++;
    }

    // Handle any remaining bits in the buffer that don't fit into a full byte
    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        dest[dest_index++] = byte;
        *destsize += bits_in_buffer;
    }
    // No need to null-terminate binary data
}


void reverse_quantize(const int* quantized_data, size_t size, float min, float bucket_size, float* data) {
    for (size_t i = 0; i < size; i++) {
        // Calculate the approximate original value by finding the midpoint of the bucket
        data[i] = min + quantized_data[i] * bucket_size + (bucket_size / 2.0f);
    }
}



double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / (compression_time * 1024 * 1024);
}

void free_tree(Node* root) {
    if (!root)
        return;

    // Recursively free left and right subtrees
    if (root->left) {
        free_tree(root->left);
    }

    if (root->right) {
        free_tree(root->right);
    }

    // Free the current node
    delete root; // Use delete instead of free for C++ objects
}


int compress(float* src, char* dest, size_t num_elements, int* destsize, float error,
            float* min, float* bucket_size, int* frq_count, int** frqs, int*** encodings) {
    // Quantization
    int* quantized_src = (int*)malloc(num_elements * sizeof(int));
    if (NULL == quantized_src) {
        printf("Error allocating memory for quantized_src\n");
        return 1;
    }

    *frq_count = quantize(src, num_elements, error, quantized_src, min);
    *bucket_size = error * 2; // As per your original code

    // Allocate and initialize frequency array
    *frqs = (int*)malloc(*frq_count * sizeof(int));
    if (NULL == *frqs) {
        printf("Error allocating memory for frqs\n");
        free(quantized_src);
        return 1;
    }
    memset(*frqs, 0, (*frq_count) * sizeof(int));

    // Record frequencies
    record_frequencies(quantized_src, num_elements, *frq_count, *frqs);

    // Build Huffman Tree
    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
    make_queue(tree, *frqs, *frq_count);
    Node* root = build_tree(tree);
    assign_encode(root);

    // Allocate encodings array as a contiguous block
    int* encodings_block = (int*)malloc(*frq_count * 2 * sizeof(int));
    if (NULL == encodings_block) {
        printf("Error allocating memory for encodings_block\n");
        free(quantized_src);
        free(*frqs);
        free_tree(root);
        return 1;
    }

    // Allocate encodings as array of pointers to each pair
    *encodings = (int**)malloc(*frq_count * sizeof(int*));
    if (NULL == *encodings) {
        printf("Error allocating memory for encodings pointers\n");
        free(encodings_block);
        free(quantized_src);
        free(*frqs);
        free_tree(root);
        return 1;
    }
    for (int i = 0; i < *frq_count; ++i) {
        (*encodings)[i] = encodings_block + (i * 2);
    }

    // Store encodings
    store_encodings(root, encodings_block, *frq_count);

    // Perform compression
    compress_quantization_levels(quantized_src, dest, num_elements, destsize, encodings_block);

    // Free allocated memory
    free(quantized_src);
    free(*frqs);
    free(encodings_block);
    free(*encodings);
    free_tree(root);

    return 0;
}
