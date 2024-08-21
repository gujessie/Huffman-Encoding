#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <cstdlib>
#include <iostream>

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

int quantize(float* src, int size, float error, int* quantized_src, float* min, float* bucket_size);
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
    long fsize = ftell(input) / sizeof(float);
    fseek(input, 0, SEEK_SET);

    float *src = (float *)malloc(fsize);
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

    float min, bucket_size;
    int *quantized_src = (int *)malloc(fsize);
    if (NULL == quantized_src) {
        printf("Error allocating memory\n");
        free(src);
        return 1;
    }
    int num_buckets = quantize(src, fsize, error, quantized_src, &min, &bucket_size);

    int *frqs = (int *)malloc(num_buckets * sizeof(int *));
    if (NULL == frqs) {
        printf("Error allocating memory\n");
        free(src);
        free(quantized_src);
        return 1;
    }

    record_frequencies(quantized_src, fsize, num_buckets, frqs);

    std::priority_queue<Node> tree;
    make_queue(tree, frqs, num_buckets);
    Node* root = build_tree(tree); 
    assign_encode(root);
    int encodings[num_buckets][2];
    store_encodings(root, encodings, num_buckets);
    


    /* compression */


    //open output file 
    FILE *output = fopen(argv[2], "wb");
    if (output == NULL) {
        perror("Error opening output file");
        fclose(input);
        return 1;
    }

    //intialize dest buffer
    char *dest = (char *)malloc(fsize); //create buffer 
    if (NULL == dest) {
       perror("Error allocating memory for dest \n");
        fclose(output);
        free(src);
        return 1;
    }

    int destsize; //size of dest in bits
    //timer for compression + compress function
    clock_t start = clock();
    compress(quantized_src, dest, fsize, &destsize, encodings);
    clock_t end = clock();

    //write data to output file
    fwrite(&min, sizeof(float), 1, output);
    fwrite(&bucket_size, sizeof(float), 1, output);

    // Store the number of elements in the output file
    fwrite(&fsize, sizeof(int), 1, output);

    //write compressed data to ouput file
    fwrite(dest, 1, (destsize + 7) / 8, output);

    fclose(output);
    free(src);
    free(dest);

    //calculate compressions speed
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    double comp_speed = calc_speed(fsize, time);


    /* decompression */

    // Read compressed data from file into a string
    FILE *comp_input = fopen(argv[2], "rb");
    if (comp_input == NULL) {
        printf("Error opening file\n");
        return 1;
    }

    fseek(comp_input, 0, SEEK_END);
    long fsize = ftell(comp_input);
    fseek(comp_input, 0, SEEK_SET);

    // Read min value and bucket size
    float min, bucket_size;
    fread(&min, sizeof(float), 1, input);
    fread(&bucket_size, sizeof(float), 1, input);

    // Read quantized data size
    int num_elements;
    fread(&num_elements, sizeof(int), 1, input);


    char *compressed_data = (char *)malloc(fsize);
    if (NULL == compressed_data) {
        printf("Error allocating memory\n");
        fclose(input);
        return 1;
    }

    //read compressed data from file 
    fread(compressed_data, fsize, 1, comp_input);
    fclose(comp_input);


    // Initialize decompressed data buffer
    float *decompressed_data = (float *)malloc(fsize + 1);
    if (NULL == decompressed_data) {
        printf("Error allocating memory\n");
        free(compressed_data);
        return 1;
    }

    // Decompress
    decompress(compressed_data, decompressed_data, root, num_elements, min, bucket_size);

    // Write decompressed data to output file
    FILE *decomp_output = fopen(argv[3], "w");
    if (decomp_output == NULL) {
        printf("Error opening file\n");
        free(compressed_data);
        free(decompressed_data);
        return 1;
    }
    fwrite(decompressed_data, sizeof(float), num_elements, output);
    fclose(decomp_output);

    free(compressed_data);
    free(decompressed_data);
    free_tree(root);

    return 0;


    //HOMEWORK
    //test larger inputs, 100 MB!!!!!!!!!!!!!
    //floating point data: data generator, generate random floating point data and write to file
    //  input data: floating point data, and will be converted to a series of bucket numbers, then use huffman to encode bucket numbers

 
}


int quantize(float* src, int size, float error, int* quantized_src, float* min, float* bucket_size) {
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
    *bucket_size = interval / num_buckets;

    // Assign each data point to a bucket
    for (int i = 0; i < size; i++) {
        int bucket = (int)((src[i] - *min) / *bucket_size);
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
void make_queue(std::priority_queue<Node> &tree, int *frqs, int num_buckets) {
    for (int i = 0; i < num_buckets; ++i) {
        if (frqs[i] > 0) { // Check if the bucket has a frequency greater than 0
            Node node;
            node.frq = frqs[i];  // Frequency
            node.index = i; // Quantized bucket value
            node.left = NULL;
            node.right = NULL;
            node.encode = 0;
            node.encode_length = 0;
            tree.push(node);  // Push the node into the priority queue
        }
    }
}

//creates huffman tree
Node* build_tree(std::priority_queue<Node> &tree){

    while(tree.size() > 1){
        //allocates space for left_node
       Node* left_node = (Node*)malloc(sizeof(Node));
        if (!left_node) {
            exit(EXIT_FAILURE);
        }
        *left_node = tree.top();
        tree.pop();

        //allocates space for right_node
        Node* right_node = (Node*)malloc(sizeof(Node));
        if (!right_node) {
            free(left_node);
            exit(EXIT_FAILURE);
        }
        *right_node = tree.top();
        tree.pop();

        //allocates space for new_node
        Node* new_node = (Node*)malloc(sizeof(Node));
        if (NULL == new_node) {
            free(left_node);
            free(right_node);
            exit(EXIT_FAILURE);
        }
        //creates new node
        new_node->frq = left_node->frq + right_node->frq;
        new_node->index = -1;
        new_node->left = left_node;
        new_node->right = right_node;
        new_node->encode = 0;
        new_node->encode_length = 0;
        tree.push(*new_node);
    }

    //returns root node
    Node* root = (Node*)malloc(sizeof(Node));
    if (root == NULL) {
        exit(EXIT_FAILURE);
    }
    *root = tree.top();
    tree.pop();
    return root;

    return NULL;
}

//assigns the binary code to each node
void assign_encode_helper(Node *root, unsigned int encode, int length) {
    if (NULL == root) return;

    if (root->index >= 0) {
        root->encode = encode;
        root->encode_length = length;
    }

    if (root->left) {
        assign_encode_helper(root->left, (encode << 1), length + 1);  
    }

    if (root->right) {
        assign_encode_helper(root->right, ((encode << 1) | 1), length + 1);  
    }
}

void assign_encode(Node *root) {
    assign_encode_helper(root, 0, 0);
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
    for (int i = 0; i <= num_buckets; ++i) {
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
        int value = src[src_index];
        int encode = encodings[value][0]; // Get Huffman encoding for the value
        int length = encodings[value][1]; // Get the length of the encoding

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
