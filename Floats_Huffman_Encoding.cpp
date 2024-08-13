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
   float val;
   float frq;
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

void quantize(float* data, int size, int num_buckets, float* quantized_data);
void record_frequencies(const float* quantized_data, int fsize, int num_buckets, float** frqs);
void make_queue(std::priority_queue<Node> &tree, float **frqs, int num_buckets);
Node* build_tree(std::priority_queue<Node> &tree);


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

    float *data = (float *)malloc(fsize);
    if (NULL == data) {
        printf("Error allocating memory\n");
        fclose(input);
        return 1;
    }
    fread(data, fsize, 1, input);
    fclose(input);

     // Quantization
    const int num_buckets = 10; // can change this number
    float *quantized_data = (float *)malloc(fsize);
    if (NULL == quantized_data) {
        printf("Error allocating memory\n");
        free(data);
        return 1;
    }
    quantize(data, fsize / sizeof(float), num_buckets, quantized_data);

    float **frqs = (float **)malloc(num_buckets * sizeof(float *));
    for (int i = 0; i < num_buckets; ++i) {
        frqs[i] = (float *)malloc(2 * sizeof(float)); // bucket_val, frq
    }
    record_frequencies(quantized_data, fsize / sizeof(float), num_buckets, frqs);

    std::priority_queue<Node> tree;
    make_queue(tree, frqs, num_buckets);
    Node* root = build_tree(tree); 


    //HOMEWORK
    //test larger inputs, 100 MB!!!!!!!!!!!!!
    //floating point data: data generator, generate random floating point data and write to file
    //  quantization: try to write a quantization routine
    //  divide min an max into buckets between the interval (btwn min and max) and then do quantization
    //  in quantization step, convert all values in the bucket into the bucket number
    //  for example, convert all numbers btwn 0 and 1 into 0 , 1 and 2 into 1
    //  dont want bucket numbers to be chars so change huffman encoding to work with integers and not characters
    //  change to work w integers 
    //  also change file generator 
    //  input data: floating point data, and will be converted to a series of bucket numbers, then use huffman to encode bucket numbers

    //makefile go back to -o0, makes it easier to debug
 
}


void quantize(float* data, int size, int num_buckets, float* quantized_data) {
    // Find min and max values in data
    float min = data[0];
    float max = data[0];
    for (int i = 1; i < size; i++) {
        if (data[i] < min) min = data[i];
        if (data[i] > max) max = data[i];
    }

    // calc the interval and bucket width
    float interval = max - min;
    float bucket_width = interval / num_buckets;

    // Assign each data point to a bucket
    for (int i = 0; i < size; i++) {
        float bucket_index = (data[i] - min) / bucket_width;

        // Store the bucket index as a float
        quantized_data[i] = bucket_index;
    }
}

void record_frequencies(const float* quantized_data, int fsize, int num_buckets, float** frqs) {
    // make frq array 
    for (int i = 0; i < num_buckets; i++) {
        frqs[i][0] = i; // Bucket val
        frqs[i][1] = 0; // frq
    }

    // Count frqs
    for (int i = 0; i < fsize; i++) {
        int bucket = (int)(quantized_data[i]);
        if (bucket >= 0 && bucket < num_buckets) {
            frqs[bucket][1]++; // Increment the frequency for the bucket
        }
    }
}

//inserts nodes that have a frequency(exists in the data) into the priority queue tree
void make_queue(std::priority_queue<Node> &tree, float **frqs, int num_buckets) {
    for (int i = 0; i < num_buckets; ++i) {
        if (frqs[i][1] > 0) { // Check if the bucket has a frequency greater than 0
            Node node;
            node.frq = frqs[i][1];        // Frequency
            node.val = frqs[i][0]; // Quantized bucket value
            node.left = nullptr;
            node.right = nullptr;
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






