#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <cstdlib>
#include <iostream>

#define MAXLENGTH 1000000


/* struct that represts a letter */
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
void compress(FILE *input, Node* root, const char* output_filename);
void decompress(const char* input_filename, Node* root, const char* output_filename);
void free_tree(Node* root);
Node* find_node(Node* root, char letter);
void store_encodings(Node* root, int encodings[128][2]);
void assign_encode_helper(Node *root, unsigned int encode, int length);

int main (int argc, char *argv[])
{
    int ascii_letters[128] = { 0 };
    Node letters[128];
    std::priority_queue<Node> tree;

    
    if ( argc != 2 ) 
    {
        return 1;
    }
   
    FILE *input = fopen( argv[1], "r" );
    if (NULL == input) { 
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    
    char line[MAXLENGTH];
    while (fgets(line, sizeof(line), input)) {
        record_letters(line, ascii_letters);
    }
    fseek(input, 0, SEEK_SET);

    make_queue(tree, ascii_letters);
    Node* root = build_tree(tree);
    assign_encode(root);
    compress(input, root, argv[2]);
    fclose(input);

    decompress(argv[2], root, argv[3]);
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
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters) { //change function name, change ascii_letters parameter to pointer
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

void compress(FILE *input, Node* root, const char* output_filename) {
    FILE *output = fopen(output_filename, "wb");
    if (output == NULL) {
        printf("Error opening output file\n");
        return;
    }

    char line[MAXLENGTH];
    unsigned int buffer = 0;
    int bits_in_buffer = 0;

    while (fgets(line, sizeof(line), input)) {
        int index = 0;
        char letter;
        Node* node;

        while (line[index] != '\0') {
            letter = line[index];
            node = find_node(root, letter);
            if (node) {
                buffer = (buffer << node->encode_length) | node->encode;
                bits_in_buffer += node->encode_length;

                while (bits_in_buffer >= 8) {
                    unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
                    fwrite(&byte, 1, 1, output);
                    bits_in_buffer -= 8;
                }
            }
            index++;
        }
    }

    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        fwrite(&byte, 1, 1, output);
    }

    fclose(output);
}

/* store_encode and its helper creates 2D array where index is a char's 
ascii value, first col is encode, second is encode length*/
void store_encodings(Node* root, int encodings[128][2]) {
    for (int i = 0; i < 128; ++i) {
        encodings[i][0] = 0;   // Initialize encode
        encodings[i][1] = 0;   // Initialize encode length
    }
    store_encodings_helper(root, encodings);
}

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

void compress(FILE *input, int encodings[128][2], const char* output_filename) {
    FILE *output = fopen(output_filename, "wb");
    if (output == NULL) {
        printf("Error opening output file\n");
        return;
    }

    char line[MAXLENGTH];
    unsigned int buffer = 0; //buffer holds up to 32 bits
    int bits_in_buffer = 0; //counts bits in buffer

    //read input liine by line
    while (fgets(line, sizeof(line), input)) {
        int index = 0;
        char letter;

        //read each char in line 
        while (line[index] != '\0') {
            letter = line[index];
            int encode = encodings[(int)letter][0]; // get encode
            int length = encodings[(int)letter][1]; //get encode_length

            //add line[i]'s encoding to buffer
            buffer = (buffer << length) | encode;
            bits_in_buffer += length; //bits_in_buffer updated w/ length

            //write to file byte by byte until bits_in_buffer is less than 8
            while (bits_in_buffer >= 8) {
                /* byte stores 8 bits from buffer and buffer is shifted accordingly;
                buffer keeps any bits that do not fit into the 8 bytes for 
                the next iteration */
                unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
                fwrite(&byte, 1, 1, output);
                bits_in_buffer -= 8; //update bits_in_buffer
            }
            index++;
        }
    }

    //handles last byte; no extra 0's
    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        fwrite(&byte, 1, 1, output);
    }

    fclose(output);
}



void decompress(const char* input_filename, Node* root, const char* output_filename) {
    FILE *input = fopen(input_filename, "rb");
    if (input == NULL) {
        printf("Error opening input file\n");
        return;
    }

    FILE *output = fopen(output_filename, "w");
    if (output == NULL) {
        printf("Error opening output file\n");
        fclose(input);
        return;
    }

    unsigned int buffer = 0; //buffer 
    int bits_in_buffer = 0; //counts bits in buffer
    unsigned char byte;
    Node* current = root; 

    //read byte by byte from input file
    while (fread(&byte, 1, 1, input) == 1) {
        buffer = (buffer << 8) | byte; //add a byte to buffer
        bits_in_buffer += 8;

        //read but by bit in buffer
        while (bits_in_buffer > 0) {
            int bit = (buffer >> (bits_in_buffer - 1)) & 1;
            bits_in_buffer--;

            //traverse tree according to bits
            if (bit == 0) {
                current = current->left;
            } else {
                current = current->right;
            }

            //if at a leaf node, write character to file
            if (current->left == NULL && current->right == NULL) {
                fputc(current->letter_name, output);
                current = root; //reset tree
            }
        }
    }

    fclose(input);
    fclose(output);
}
