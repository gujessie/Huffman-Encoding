#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ncurses.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <cstdlib>
#include <iostream>

#define MAXLENGTH 1000000


//write makefile for this code, autotools, cmake, gnumake


/* struct that represts a letter */
struct node{
   struct node *left;
   struct node *right;
   int frq;
   char letter_name;
   int encode;
   int encode_length;

   bool operator<(const node& comp_node) const {
        if (frq > comp_node.frq) {
            return true;
        }
        return false;
    }
};
typedef struct node Node;


void record_letters(char *string, int ascii_letters[]);
int compare(const void* a, const void* b);
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode(std::priority_queue<Node> &tree);
void assign_encode_helper(Node *root, unsigned int encode, int length);
char* compress(char *stringtemp, Node letters[]);
void free_tree(Node* root);
Node* find_node(Node* root, char letter);


int main() {
    // Test string
    char input_string[] = "aaabbcadddighseoirhgioaerhgiohaeroighoerhaghoiaehrogherhgihaeriohgioaherioghoiaerhoighaeoirhgihaeroighoiaerhgiohaeroighoiaehrogiheoirhgiohaeroighoaiehrgioaheroihgioaherioghaoierhgiohaeroighaoiehrgoihearioghaioerhgioaerhgiohaeroihgaiohreogihaeiorhgoaiheroighaeoirhgoehrgohaeiorhgioheroihgoiheroighaoeirhgoiherdikfhgkdhfksldjbo;jeropjgiohaerkhdgjkndksbjoajerighkaihgbajolerjh";

    // Array to store ASCII frequencies
    int ascii_letters[128];

    // Record the frequencies of each letter in the string
    record_letters(input_string, ascii_letters);
        int i = 0;
        for (i = 0; i < 128; i++){
            if(ascii_letters[i] != 0){
                //printf("%c: %d\n", (char)i, ascii_letters[i]);
            }
        }
       

    // Array to hold nodes for Huffman tree construction
    std::priority_queue<Node> tree;
    make_queue(tree, ascii_letters);

    // Number of nodes with non-zero frequency
    int num_nodes = 0;
    for (i = 0; i < 128; ++i) {
        if (ascii_letters[i] != 0) {
            num_nodes++;
        }
    }


    // Build the Huffman tree from the nodes
    Node* root = build_tree(tree);

    assign_encode_helper(root, 0, 0);

    // Print the Huffman encodings for each character in the input string
    // printf("Huffman Encodings:\n");
    // for (char& c : input_string) {
    //     Node* node = find_node(root, c);
    //     if (node != nullptr) {
    //         printf("%c: %d, %d\n", c, node->encode, node->encode_length);
    //     }
    // }

    // free
    free(root);

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

//fix, dont push struct
//inserts nodes that have a frequency(exists in the data) into the priority queue tree
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters) { //change function name, change ascii_letters parameter to pointer
    // Insert nodes 
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
            tree.push(node); //pass pointer
        }
    }
    // Node* node = tree.top();
    // printf("%c, %d" ,node->frq, node->letter_name);

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
    if (!root){
        return;
    }
    //printf("root freq: %d", root->frq);
    if (root->letter_name != '\0') {
        //printf("%c;\n", root->letter_name);
        root->encode = encode;
        root->encode_length = length;
        //printf("Assigned encode %x with length %d to letter %c\n", root->encode, root->encode_length, root->letter_name);

        printf ("encode = %c, len = %d\n", root->letter_name, length);
    }

    if (root->left) {
        if (root->left->letter_name != '\0') {
            //printf("%c, %d, %d\n", root->left->letter_name, root->left->encode, root->left->encode_length);
        }
        assign_encode_helper(root->left, (encode << 1), length + 1); 
    }

    if (root->right) {
        if (root->right->letter_name != '\0') {
            //printf("%c, %d, %d\n", root->right->letter_name, root->right->encode, root->right->encode_length);
        }
        assign_encode_helper(root->right, ((encode << 1) | 1), length + 1);  
    }
}

// void assign_encode(std::priority_queue<Node> &tree) {
//     Node root = tree.top();
//     assign_encode_helper(&root, 0, 0);
// }

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
