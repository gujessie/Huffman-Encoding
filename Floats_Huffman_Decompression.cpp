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

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "error in arguments", argv[0]);
        return 1;
    }

    // Opening compressed input file
    FILE *input = fopen(argv[1], "rb");
    if (input == NULL) { 
        perror("Error opening input file");
        return 1;
    }

    // Determine the size of the compressed file
    fseek(input, 0, SEEK_END); // Move to the end of the file
    long file_size = ftell(input); // Get current position (file size)
    fseek(input, 0, SEEK_SET); // Reset to beginning of the file

    printf("file size = %lu\n", file_size);

    // Allocate memory for the entire compressed buffer
    char *compressed_buffer = (char *) malloc(file_size);
    if (compressed_buffer == NULL) {
        perror("Error allocating buffer for compressed data");
        fclose(input);
        return 1;
    }

    // Read the entire compressed file into the buffer
    size_t read_bytes = fread(compressed_buffer, 1, file_size, input);
    if (read_bytes != file_size) {
        fprintf(stderr, "Error reading compressed file\n");
        free(compressed_buffer);
        fclose(input);
        return 1;
    }
    fclose(input);

    size_t num_elements;
    memcpy(&num_elements, compressed_buffer, sizeof(size_t));

    float* dest = (float *) malloc (num_elements * sizeof(float)); // Over-allocate
    if (dest == NULL) {
        perror("Error allocating buffer for decompressed data");
        free(compressed_buffer);
        return 1;
    }

    int destsize = 0; // To hold the size of decompressed data

    // Perform decompression
    int decompress_status = decompress(compressed_buffer, dest, file_size, &destsize);
    if (decompress_status != 0) {
        fprintf(stderr, "Decompression failed\n");
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    // Replace writefile with direct fwrite
    FILE* output = fopen(argv[2], "wb");
    if (output == NULL) {
        perror("Error opening output file");
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    size_t written = fwrite(dest, sizeof(float), destsize, output);
    if (written != (size_t)destsize) {
        fprintf(stderr, "Error writing to output file\n");
        fclose(output);
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    fclose(output);

    // Free allocated memory
    free(compressed_buffer);
    free(dest);

    return 0;
}

// Decompress helper function
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
        while (bits_in_buffer > 0 && dest_index < num_elements) {
            int bit = (buffer >> (bits_in_buffer - 1)) & 1;
            bits_in_buffer--;

            if (bit == 0) {
                current = current->left;
            } else {
                current = current->right;
            }

            //for leaf node, decode the quantized value
            if (current->left == NULL && current->right == NULL) {
                int bucket_index = current->index;
                // Convert bucket index back to the quantized float value
                float quantized_value = min + bucket_index * bucket_size + (bucket_size / 2.0f);
                dest[dest_index++] = quantized_value; // Store the reconstructed float value
                current = root; // Reset to the root for the next decoding
            }
        }
    }
}

// Decompress function implementation
#if 0
file size = 129623407
num_elements: 134217728
min = 0.000000, bucket_size = 0.500000
frq_count = 201
Compressed data size = 129620975 bytes
offset size = 2432 bytes
#endif
int decompress(const char* src, float* dest, size_t src_size, int* destsize){
    // Initialize pointer to traverse the src buffer
    size_t offset = 0;

    // Read metadata from the src buffer
    size_t num_elements;
    memcpy(&num_elements, src + offset, sizeof(size_t));
    offset += sizeof(size_t);

    float min;
    memcpy(&min, src + offset, sizeof(float));
    offset += sizeof(float);

    float bucket_size;
    memcpy(&bucket_size, src + offset, sizeof(float));
    offset += sizeof(float);

    int frq_count;
    memcpy(&frq_count, src + offset, sizeof(int));
    offset += sizeof(int);

    printf("num_elements: %zu\n", num_elements);
    printf("min = %f, bucket_size = %f\n", min, bucket_size);
    printf("frq_count = %d\n", frq_count);

    // Read frequency array
    int *frqs = (int *)malloc(frq_count * sizeof(int));
    if (frqs == NULL) {
        perror("Error allocating memory for frequencies");
        return 1;
    }
    memcpy(frqs, src + offset, frq_count * sizeof(int));
    offset += frq_count * sizeof(int);

    // Read encodings array
    int* encodings = (int*)malloc(frq_count * 2 * sizeof(int));
    if (encodings == NULL) {
        perror("Error allocating memory for encodings");
        free(frqs);
        return 1;
    }
    memcpy(encodings, src + offset, frq_count * 2 * sizeof(int));
    offset += frq_count * 2 * sizeof(int);

    // Rebuild Huffman Tree
    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
    make_queue(tree, frqs, frq_count);
    Node* root = build_tree(tree);
    assign_encode(root);
    store_encodings(root, encodings, frq_count);

    // Calculate the size of compressed data
    size_t compressed_size = src_size - offset;
    printf("Compressed data size = %zu bytes\n", compressed_size);

    printf("offset size = %u bytes\n", offset);

    // Read compressed data
    const char* compressed_data = src + offset;

    // Perform decompression
    decompress_1(compressed_data, dest, root, num_elements, min, bucket_size);

    // Update destsize 
    *destsize = num_elements;

    // Free allocated memory
    free(frqs);
    free(encodings);
    free_tree(root);

    return 0;
}
