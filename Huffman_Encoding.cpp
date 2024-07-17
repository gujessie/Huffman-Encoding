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
        if (frq > comp_node.frq) {
            return true;
        }
        return false;
    }
};
typedef struct node Node;

void record_letters(const char *string, int ascii_letters[]);
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode_helper(Node *root, unsigned int encode, int length);
void assign_encode(Node *root);
void compress(FILE *input, Node* root, const char* output_filename);
void decompress(const char* input_filename, Node* root, const char* output_filename);
void free_tree(Node* root);
Node* find_node(Node* root, char letter);

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
            tree.push(node); 
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

void assign_encode(std::priority_queue<Node> &tree) {
    Node root = tree.top();
    assign_encode_helper(&root, 0, 0);
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
        printf("Error opening file\n");
        return;
    }

    char line[MAXLENGTH];
    unsigned char buffer = 0;
    int buffer_index = 0;

    while (fgets(line, sizeof(line), input)) {
        int index = 0;
        char letter;
        Node* node;

        while (line[index] != '\0') {
            letter = line[index];
            node = find_node(root, letter);
            if (node) {
                int i;
                for (i = 0; i < node->encode_length; i++) {
                    buffer <<= 1;
                    buffer |= (node->encode >> (node->encode_length - i - 1)) & 1;
                    buffer_index++;
                    if (buffer_index == 8) {
                        fwrite(&buffer, 1, 1, output);
                        buffer = 0;
                        buffer_index = 0;
                    }
                }
            }
            index++;
        }
    }

    if (buffer_index > 0) {
        buffer <<= (8 - buffer_index);
        fwrite(&buffer, 1, 1, output);
    }

    fclose(output);
}


void decompress(const char* input_filename, Node* root, const char* output_filename) {
    FILE *input = fopen(input_filename, "rb");
    FILE *output = fopen(output_filename, "w");
    if (input == NULL || output == NULL) {
        printf("Error opening file\n");
        return;
    }

    unsigned char buffer;
    int buffer_index = 0;
    Node* current = root;

    while (fread(&buffer, 1, 1, input) == 1) {
        int i;
        for (i = 7; i >= 0; i--) {
            if (buffer & (1 << i)) {
                current = current->right;
            } else {
                current = current->left;
            }

            if (current->left == NULL && current->right == NULL) {
                fputc(current->letter_name, output);
                current = root;
            }
        }
    }

    fclose(input);
    fclose(output);
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
//     
//     }

//     return compdata;
// }



// //sorts letters from greatest frequency to least, and records the frequency and name
// int compare(const void* a, const void* b) {
//     Node* nodeA = (Node*)a;
//     Node* nodeB = (Node*)b;
//     return nodeB->frequency - nodeA->frequency; // For descending order
// }

// void sort_letters(Node letters[128], int ascii_letters[128]) {
//     int count = 0;

//     for (int i = 0; i < 128; i++) {
//         if (ascii_letters[i] > 0) {
//             letters[count].letter_name = (char)i;
//             letters[count].frequency = ascii_letters[i];
//             count++;
//         }
//     }

//     qsort(letters, count, sizeof(Node), compare);
// }


// //build the tree
// void build_tree(Node letters[128], Node nodes[192]) {
//     int count = 0, i;

//     //make nodes array with character nodes
//     for (i = 0; i < 128; i++) {
//         if (letters[i].frq > 0){
//             nodes[count].frq = letters[i].frq;
//             nodes[count].letter_name = letters[i].letter_name;
//             nodes[count].left = NULL;
//             nodes[count].right = NULL;
//             count++;
//         }
//     }

//     // build tree
//     while (count > 1) {
//         Node* left = &nodes[count - 1]; //call pop on priority queue to get two least frequency nodes
//         Node* right = &nodes[count - 2]; 

//         //malloc new node and then set left and right tree
//         // make sum node
//         Node* new_node = &nodes[count];
//         new_node->letter_name = NULL; 
//         new_node->frq = left->frq + right->frq;
//         new_node->left = left;
//         new_node->right = right;

//         count--;

//         //insert
//         int c = count -1;  
//         while (c >= 0 && nodes[c].frq > new_node->frq) {
//             nodes[c + 1] = nodes[c];
//             c--;
//         }
//         nodes[c + 1] = *new_node;
//         count++;
//     }
// }
  

// while (phtree->size() > 1) {

//     htree_node *top_node1 = phtree->top();

//     phtree->pop();

//     htree_node *top_node2 = phtree->top();

//     phtree->pop();

 

//     htree_node *new_node = new_htree_node(-1, top_node1->cnt + top_node2->cnt);

//     new_node->left = top_node1;

//     new_node->right = top_node2;

//     phtree->push(new_node);
