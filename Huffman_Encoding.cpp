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
   struct node *left;
   struct node *right;
   int frq;
   char letter_name;
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

void record_letters(char *string, int ascii_letters[]);
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode_helper(Node *root, unsigned int encode, int length);
void assign_encode(Node *root);
void compress(const char* src, char* dest, int fsize, int* destsize, int encodings[128][2]);
void decompress(const char* src, char* dest, Node* root, int fsize);
void free_tree(Node* root);
Node* find_node(Node* root, char letter);
void store_encodings(Node* root, int encodings[128][2]);
void assign_encode_helper(Node *root, unsigned int encode, int length);
double calc_speed(long input_size, double time);

int main (int argc, char *argv[])
{
    int ascii_letters[128] = { 0 };
    Node letters[128];
    long fsize;
    std::priority_queue<Node> tree;

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
    fsize = ftell(input);
    fseek(input, 0, SEEK_SET);

    //record frequency of letters
   
    //intitialize src buffer
    char *src = (char *)malloc(fsize);
    if (NULL == src) {
        printf("Error allocating memory for src\n");
        fclose(input);
        return 1;
    }
    fread(src, 1, fsize, input);
    fclose(input);

   //record letter frqs
    record_letters(src, ascii_letters);
    fseek(input, 0, SEEK_SET);

    //make tree and encodings array
    make_queue(tree, ascii_letters);
    Node* root = build_tree(tree);
    assign_encode(root);
    int encodings[128][2];
    store_encodings(root, encodings);

    /* compression */


    //open output file 
    FILE *output = fopen(argv[2], "wb");
    if (output == NULL) {
        perror("Error opening output file");
        fclose(input);
        return 1;
    }

    //intitialize src buffer
    char *src = (char *)malloc(fsize + 1);
    if (NULL == src) {
        printf("Error allocating memory for src\n");
        fclose(input);
        fclose(output);
        return 1;
    }
    fread(src, 1, fsize, input);
    fclose(input);
    src[fsize] = '\0';

    //intialize dest buffer
    char *dest = (char *)malloc(fsize + 1); //create buffer 
    if (NULL == dest) {
       perror("Error allocating memory for dest \n");
        fclose(output);
        free(src);
        return 1;
    }

    //store the ascii_letters array in output file
    unsigned char array_size = 128; 
    fwrite(&array_size, 1, 1, output);  // write size as first byte
    fwrite(ascii_letters, sizeof(int), 128, output);  // write ascii_letters array next

    int destsize; //size of dest in bits
    //timer for compression + compress function
    clock_t start = clock();
    compress(src, dest, fsize, &destsize, encodings);
    clock_t end = clock();

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
    long compressed_size = ftell(comp_input);
    fseek(comp_input, 0, SEEK_SET);

    char *compressed_data = (char *)malloc(compressed_size + 1);
    if (NULL == compressed_data) {
        printf("Error allocating memory\n");
        fclose(input);
        return 1;
    }


    // read ascii_letters array size and skip to encode part of file
    fread(&array_size, 1, 1, input);
    fseek(input, array_size * sizeof(int), SEEK_CUR);

    //read compressed data from file 
    fread(compressed_data, compressed_size, 1, comp_input);
    fclose(comp_input);
    compressed_data[compressed_size] = '\0';

    // Initialize decompressed data buffer
    char *decompressed_data = (char *)malloc(fsize + 1);
    if (NULL == decompressed_data) {
        printf("Error allocating memory\n");
        free(compressed_data);
        return 1;
    }

    // Decompress
    decompress(compressed_data, decompressed_data, root, fsize);

    // Write decompressed data to output file
    FILE *decomp_output = fopen(argv[3], "w");
    if (decomp_output == NULL) {
        printf("Error opening file\n");
        free(compressed_data);
        free(decompressed_data);
        return 1;
    }
    fwrite(decompressed_data, 1, strlen(decompressed_data), output);
    fclose(decomp_output);

    free(compressed_data);
    free(decompressed_data);
    free_tree(root);

    return 0;
}




/* iterates through a string and records the frequncy of the characters by 
incrementing the value corresponding index (ascii value) */
void record_letters(char *string, int ascii_letters[]){
    int index = 0, ascii_value;
    char letter;

    while (string[index] != '\0'){
        letter = string[index];
        ascii_value = (int)letter;
        ascii_letters[ascii_value] = ascii_letters[ascii_value] + 1;
        index++;    
    }
}

//inserts nodes that have a frequency(exists in the data) into the priority queue tree
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters) { //change ascii_letters parameter to pointer
    // Insert nodes 


    //change in the future
    int i;
    for (i = 0; i < 128; ++i) {
        if (ascii_letters[i] > 0) { 
            Node node;
            node.frq = ascii_letters[i];
            node.letter_name = (char)(i);
            node.left = NULL;
            node.right = NULL;
            node.encode = 0;
            node.encode_length = 0;
            tree.push(node);   //use pointers
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
        new_node->letter_name = '\0';
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

//assigns the binary code to each letter
void assign_encode_helper(Node *root, unsigned int encode, int length) {
    if (NULL == root) return;

    if (root->letter_name != '\0') {
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

Node* find_node(Node* root, char letter) {
    if (!root) {
        return NULL;
    }

    if (root->letter_name == letter) {
        return root;
    }

    Node* node = find_node(root->left, letter);
    if (!node) {
        node = find_node(root->right, letter);
    }

    return node;
}

/* store_encode and its helper creates 2D array where index is a char's 
ascii value, first col is encode, second is encode length*/
void store_encodings_helper(Node* root, int encodings[128][2]) {
    if (!root) return;

    if (root->letter_name != '\0') {
        encodings[(int)root->letter_name][0] = root->encode;
        encodings[(int)root->letter_name][1] = root->encode_length;
    }

    if (root->left) {
        store_encodings_helper(root->left, encodings);
    }

    if (root->right) {
        store_encodings_helper(root->right, encodings);
    }
}

void store_encodings(Node* root, int encodings[128][2]) {
    for (int i = 0; i < 128; ++i) {
        encodings[i][0] = 0;   // Initialize encode
        encodings[i][1] = 0;   // Initialize encode length
    }
    store_encodings_helper(root, encodings);
}

/* note: return 1 if error and check in main*/
void compress(const char* src, char* dest, int fsize, int* destsize, int encodings[128][2]) {
    int dest_index = 0; //tracks current position in buffer
    unsigned int buffer = 0; //stores bits before writing to buffer
    int src_index = 0; //tracks current position in src 
    int bits_in_buffer = 0; //counts bits in buffer
    char letter;
    *destsize = 0; //size of dest in bits

    //read each char in line 
    while (src_index < fsize) {
        letter = src[src_index];
        int encode = encodings[(int)letter][0]; // get encode
        int length = encodings[(int)letter][1]; //get encode_length

        //add line[i]'s encoding to bit_buffer
        buffer = (buffer << length) | encode;
        bits_in_buffer += length; //bits_in_buffer updated w/ length of bit_buffer

        //write to file byte by byte until bits_in_buffer is less than 8
        while (bits_in_buffer >= 8) {
            /* byte stores 8 bits from bit_buffer, bit_ buffer shifted accordingly;
            bit_buffer keeps any bits that do not fit into the 8 bytes for 
            the next iteration */
            unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
            dest[dest_index++] = byte; //adds a byte to buffer
            bits_in_buffer -= 8; //update bits_in_buffer
            *destsize += 8;
        }
        src_index++;
    }

    //handles last byte; no extra 0's
    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        dest[dest_index++] = byte;
        *destsize += bits_in_buffer;
    }
    dest[dest_index] = '\0';
}

/* note: want to avoid code duplication for node traversal, write helper */
void decompress(const char* src, char* dest, Node* root, int fsize) {
    Node* current = root;
    unsigned int buffer = 0;
    int bits_in_buffer = 0;
    int src_index = 0;
    int dest_index = 0;

    while (dest_index < fsize) {
        // Load bits from src into buffer
        unsigned char byte = src[src_index++];
        buffer = (buffer << 8) | byte;
        bits_in_buffer += 8;

        // Decode bits into characters
        while (bits_in_buffer > 0) {
            int bit = (buffer >> (bits_in_buffer - 1)) & 1;
            bits_in_buffer--;

            if (bit == 0) {
                current = current->left;
            } else {
                current = current->right;
            }

            if (current->left == NULL && current->right == NULL) {
                dest[dest_index++] = current->letter_name;
                current = root;
            }
        }
    }

    // Null-terminate the output string
    dest[dest_index] = '\0';
}

double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / compression_time;
}


