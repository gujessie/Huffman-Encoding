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


/* struct that represts a letter */
struct node{
   struct node *left;
   struct node *right;
   int frq;
   char letter_name;
   int encode;
   int encode_length;

   bool operator<(const Node& comp_node) const {
        if (frq < comp_node.frq) {
            return true;
        }
        return false;
    }
};
typedef struct node Node;


void record_letters(char *string, int ascii_letters[]);
int compare(const void* a, const void* b);
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters);
void build_tree(std::priority_queue<Node> &tree);
void assign_encode(std::priority_queue<Node> &tree);
void assign_encode_helper(Node *root, unsigned int encode, int length);
char* compress(char *stringtemp, Node letters[]);

int main (int argc, char *argv[])
{
    int ascii_letters[128];
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
    
    fseek(input, 0, SEEK_END);
    long fsize = ftell(input);
    fseek(input, 0, SEEK_SET);  

    char *stringtemp = (char*)std::malloc(fsize + 1);
    //check return type, make sure stringtemp is not null, is so malloc has failed
    if (NULL == stringtemp){
        return 1;
    }

    fread(stringtemp, fsize, 1, input);
    fclose(input);
    stringtemp[fsize] = '\0';
 
    record_letters(stringtemp, ascii_letters);

    //sort_letters(letters, ascii_letters);
    make_queue(tree, ascii_letters);
    build_tree(tree);

    assign_encode(tree);
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
    for (int i = 0; i < 128; ++i) {
        if (ascii_letters[i] > 0) { 
            Node node;
            node.frq = ascii_letters[i];
            node.letter_name = (char)(i);
            node.left = NULL;
            node.right = NULL;
            node.encode = 0;
            node.encode_length = 0;
            tree.push(node); 
        }
    }
}

//creates huffman tree
void build_tree(std::priority_queue<Node> &tree){
    while(tree.size() > 1){
        Node left_node = tree.top();
        tree.pop();
        Node right_node = tree.top();
        tree.pop();

        Node new_node;
        new_node.letter_name = '\0'; 
        new_node.frq = left_node.frq + right_node.frq;
        new_node.left = &left_node;
        new_node.right = &right_node;
        new_node.encode = 0;  
        new_node.encode_length = 0;
        tree.push(new_node);
    }
}

void assign_encode(std::priority_queue<Node> &tree) {
    Node root = tree.top();
    assign_encode_aux(&root, 0, 0);
}

//assigns the binary code to each letter
void assign_encode_aux(Node *root, unsigned int encode, int length) {
    if (NULL == root){
       return;
    }

    if (root->letter_name != '\0') {
        root->encode = encode;
        root->encode_length = length;
    }

    if (NULL != root->left) {
        assign_encode_helper(root->left, (encode << 1), length + 1);  
    }

    if (NULL != root->right) {
        assign_encode_helper(root->right, ((encode << 1) | 1), length + 1);  
    }
}

// //write to file
// char* compress(char *stringtemp, Node letters[]){
//     char* compdata;
//     int index = 0, i;
//     char letter;

//     while (stringtemp[index] != '\0'){
//         letter = stringtemp[index];
//         for(i = 0; i < 128; i++){
//             if(letters[i].letter_name == letter){
//                 break;
//             }
//             i++;
//         }
//     //...
//     }

//     return compdata;
// }
