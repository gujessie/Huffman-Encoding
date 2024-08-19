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

int quantize(float* src, int size, float error, int* quantized_src);
void record_frequencies(const int* quantized_src, int fsize, int num_buckets, int* frqs);
void make_queue(std::priority_queue<Node> &tree, int *frqs, int num_buckets);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode(Node *root);
void store_encodings_helper(Node* root, int encodings[][2], int num_buckets);

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

    int *quantized_src = (int *)malloc(fsize);
    if (NULL == quantized_src) {
        printf("Error allocating memory\n");
        free(src);
        return 1;
    }
    int num_buckets = quantize(src, fsize, error, quantized_src);

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
    
}


int quantize(float* src, int size, float error, int* quantized_src) {
    // Find min and max values in data
    float min = src[0];
    float max = src[0];
    for (int i = 1; i < size; i++) {
        if (src[i] < min) min = src[i];
        if (src[i] > max) max = src[i];
    }

    // calc the interval and bucket width
    float interval = max - min;
    int num_buckets = (int)(interval/(error * 2) + 1);
    float bucket_size = interval / num_buckets;

    // Assign each data point to a bucket
    for (int i = 0; i < size; i++) {
        int bucket = (int)((src[i] - min) / bucket_size);
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
