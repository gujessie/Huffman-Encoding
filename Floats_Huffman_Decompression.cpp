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

    // Read: min, bucket_size, fsize, and frequencies
    float min, bucket_size;
    int num_elements;
    fread(&num_elements, sizeof(int), 1, input);
    fread(&min, sizeof(float), 1, input);
    fread(&bucket_size, sizeof(float), 1, input);
    int frq_count;
    fread(&frq_count, sizeof(int), 1, input);

    printf("num_elements: %d\n", num_elements);

    // Read frequency array from the file
    int *frqs = (int *)malloc(frq_count * sizeof(int));
    fread(frqs, sizeof(int), frq_count, input);

    // Priority queue to rebuild the Huffman tree
    std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
    make_queue(tree, frqs, frq_count);
    Node *root = build_tree(tree);



    // Calculate compressed data size, read to end of file
    fseek(input, 0, SEEK_END); // Move to the end of the file
    long file_size = ftell(input); // Get current position (file size)
    //Move to where compression starts
    fseek(input, ftell(input) - frq_count * sizeof(int) - sizeof(num_elements) - sizeof(min) - sizeof(bucket_size) - sizeof(frq_count), SEEK_SET); 

    long compressed_size = file_size - ftell(input);  // compressed_size is remaining size

    // Read compressed data from file
    char *compressed_data = (char *)malloc(compressed_size);  // Adjust size as needed
    fread(compressed_data, 1, compressed_size * 2, input);
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
    decompress(compressed_data, dest, root, num_elements, min, bucket_size);


    // Reverse quantization
    //float *decompressed_data = (float *)malloc(fsize * sizeof(float));
    // reverse_quantize(dest, fsize, min, bucket_size, decompressed_data);

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

    return 0;
}


// int main (int argc, char *argv[])
// {
//     if ( argc != 3 ) {
//         fprintf(stderr, "error in arguments");
//         return 1;
//     }
   
//    //opening input file
//     FILE *input = fopen(argv[1], "rb");
//     if (NULL == input) { 
//         perror("Error opening input file");
//         return 1;
//     }
    
//     //get input file size 
//     fseek(input, 0, SEEK_END);
//     size_t fsize = ftell(input) / sizeof(float); //num elmts in file
//     fseek(input, 0, SEEK_SET);

//     float *src = (float *)malloc(fsize * sizeof(float));
//     if (NULL == src) {
//         printf("Error allocating memory\n");
//         fclose(input);
//         return 1;
//     }
//     fread(src, sizeof(float), fsize, input);
//     fclose(input);

//     /* Quantization */
//     //calculate num buckets
//     const float error = 0.25; // can change this number

//     float min, bucket_size = error * 2;
//     int *quantized_src = (int *)malloc(fsize * sizeof(int));
//     if (NULL == quantized_src) {
//         printf("Error allocating memory\n");
//         free(src);
//         return 1;
//     }
//     int num_buckets = quantize(src, fsize, error, quantized_src, &min);

//     int *frqs = (int *)malloc(num_buckets * sizeof(int));
//     if (NULL == frqs) {
//         printf("Error allocating memory\n");
//         free(src);
//         free(quantized_src);
//         return 1;
//     }

//     record_frequencies(quantized_src, fsize, num_buckets, frqs);

//     // Priority queue to build the Huffman tree

//     std::priority_queue<Node*, std::vector<Node*>, LessThanByCnt> tree;
//     make_queue(tree, frqs, num_buckets);
//     Node* root = build_tree(tree);


//     assign_encode(root);
//     int encodings[num_buckets][2];
//     store_encodings(root, encodings, num_buckets);

 
//     /* decompression */

//     // Read compressed data from file into a string
//     FILE *comp_input = fopen(argv[2], "rb");
//     if (comp_input == NULL) {
//         printf("Error opening file\n");
//         return 1;
//     }

//     fseek(comp_input, 0, SEEK_END);
//     fsize = ftell(comp_input);
//     fseek(comp_input, 0, SEEK_SET);

//     // Read min value and bucket size
//     int num_elements;
//     fread(&min, sizeof(float), 1, comp_input);
//     fread(&bucket_size, sizeof(float), 1, comp_input);
//     // Read quantized data size
//     fread(&num_elements, sizeof(int), 1, comp_input);
//     fread(frqs, sizeof(frqs), 1, comp_input);

//     make_queue(tree, frqs, num_buckets);
//     root = build_tree(tree);

//     char *compressed_data = (char *)malloc(fsize);
//     if (NULL == compressed_data) {
//         printf("Error allocating memory\n");
//         fclose(input);
//         return 1;
//     }

//     //read compressed data from file 
//     fread(compressed_data, fsize, 1, comp_input);
//     fclose(comp_input);


//     // Initialize decompressed data buffer
//     float *decompressed_data = (float *)malloc(fsize + 1);
//     if (NULL == decompressed_data) {
//         printf("Error allocating memory\n");
//         free(compressed_data);
//         return 1;
//     }


//     // Decompress
//     decompress(compressed_data, decompressed_data, root, num_elements, min, bucket_size);


//     // Write decompressed data to output file
//     FILE *decomp_output = fopen(argv[3], "w");
//     if (decomp_output == NULL) {
//         printf("Error opening file\n");
//         free(compressed_data);
//         free(decompressed_data);
//         return 1;
//     }
//     fwrite(decompressed_data, sizeof(float), num_elements, decomp_output);
//     fclose(decomp_output);

//     free(compressed_data);
//     free(decompressed_data);
//     free_tree(root);

//     return 0;
// }


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
void decompress(const char* src, float* dest, Node* root, int num_elements, float min, float bucket_size) {
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
