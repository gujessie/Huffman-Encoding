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

int main (int argc, char *argv[])
{
    if ( argc != 3 ) {
        fprintf(stderr, "error in arguments");
        return 1;
    }
   
   //opening input file
    FILE *input = fopen(argv[1], "r");
    if (NULL == input) { 
        perror("Error opening input file");
        return 1;
    }
    
    //get input file size 
    fseek(input, 0, SEEK_END);
    size_t fsize = ftell(input) / sizeof(float); //num elmts in file
    fseek(input, 0, SEEK_SET);

    float *src = (float *)malloc(fsize * sizeof(float));
    if (NULL == src) {
        printf("Error allocating memory\n");
        fclose(input);
        return 1;
    }
    fread(src, sizeof(float), fsize, input);
    fclose(input);

    /* Quantization */
    //calculate num buckets
    const float error = 0.25; // can change this number

    float min, bucket_size = error * 2;
    int *quantized_src = (int *)malloc(fsize * sizeof(int));
    if (NULL == quantized_src) {
        printf("Error allocating memory\n");
        free(src);
        return 1;
    }
    int num_buckets = quantize(src, fsize, error, quantized_src, &min);

    int *frqs = (int *)malloc(num_buckets * sizeof(int));
    if (NULL == frqs) {
        printf("Error allocating memory\n");
        free(src);
        free(quantized_src);
        return 1;
    }


    record_frequencies(quantized_src, fsize, num_buckets, frqs);

    // Priority queue to build the Huffman tree

    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
    make_queue(tree, frqs, num_buckets);
    Node* root = build_tree(tree);
    assign_encode(root);
    int encodings[num_buckets][2];
    store_encodings(root, encodings, num_buckets);

    /* Compression */

    //open output file 
    FILE *output = fopen(argv[2], "wb");
    if (NULL == output) {
        perror("Error opening output file");
        free(src);
        free(quantized_src);
        free(frqs);
        return 1;
    }

    //intialize dest buffer
    char *dest = (char *)malloc(fsize * 2); //create buffer 
    if (NULL == dest) {
       perror("Error allocating memory for dest \n");
        fclose(output);
        free(src);
        free(quantized_src);
        free(frqs);
        return 1;
    }

using namespace std;
 
int destsize; //size of dest in bits
//timer for compression + compress function

//   using namespace std::chrono;
//   high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
//   compress(quantized_src, dest, fsize, &destsize, encodings);

//   std::cout << std::endl;
//   high_resolution_clock::time_point t2 = high_resolution_clock::now();
//   duration<double> time = duration_cast<duration<double>>(t2 - t1);
//   std::cout << "It took me " << time.count() << " seconds.";
//   std::cout << std::endl;g


    // int destsize; //size of dest in bits
    //timer for compression + compress function


    clock_t start = clock();
    compress(quantized_src, dest, fsize, &destsize, encodings);    
    clock_t end = clock();

    //write data to output file

    fwrite(&min, sizeof(float), 1, output);
    fwrite(&bucket_size, sizeof(float), 1, output);
    // Store the number of elements in the output file
    fwrite(&fsize, sizeof(int), 1, output);
    fwrite(frqs, sizeof(frqs), 1, output);
    //write compressed data to ouput file
    fwrite(dest, 1, (destsize + 7) / 8, output);

    fclose(output);

    free(src);
    free(quantized_src);
    free(frqs);
    free(dest);
    free_tree(root);

    //calculate compressions speed
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    double comp_speed = calc_speed(fsize, time);
    printf("Throughput = %f\n", comp_speed);

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

void compress(const int* src, char* dest, int fsize, int* destsize, int encodings[][2]) {
    int dest_index = 0; // Tracks current position in the destination buffer
    unsigned int buffer = 0; // Stores bits before writing to the destination buffer
    int src_index = 0; // Tracks current position in the source buffer
    int bits_in_buffer = 0; // Counts bits currently in the buffer
    *destsize = 0; // Size of the destination in bits

    // Iterate over each element in the source array
    while (src_index < fsize) {
        int bucket_index = src[src_index]; //
        int encode = encodings[bucket_index][0]; // Get Huffman encoding for the value
        int length = encodings[bucket_index][1]; // Get the length of the encoding

        // Add the current value's encoding to the bit buffer
        buffer = (buffer << length) | encode;
        bits_in_buffer += length;

        // Write to the destination buffer byte by byte until there are less than 8 bits left
        while (bits_in_buffer >= 8) {
            unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
            dest[dest_index++] = byte; // Add the byte to the destination buffer
            bits_in_buffer -= 8; // Update the number of bits in the buffer
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
    dest[dest_index] = '\0';
}

double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / (compression_time * 1024 * 1024);
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

// float interpolate(float x0, float y0, float x1, float y1, float x2) {
//     float slope = (y1 - y0) / (x1 - x0); //find slope
//     float y2 = y0 + slope * (x2 - x0); //find y2
//     return y2;
// }

// void compress_interpolation_error(float actual, float predicted, float error_bound, float *compressed_val) {
//     float difference = actual - predicted;
//     if (difference < 0) {
//         difference = -difference;
//     }
//     if (difference <= error_bound) {
//         *compressed_val = predicted;
//     } else {
//         *compressed_val = actual;
//     }
// }
