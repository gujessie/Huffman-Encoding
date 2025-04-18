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
void compress(FILE *input, int encodings[128][2], const char* output_filename, int ascii_letters[128]);
void decompress(const char* input_filename, Node* root, const char* output_filename);
void free_tree(Node* root);
Node* find_node(Node* root, char letter);
void store_encodings(Node* root, int encodings[128][2]);
void assign_encode_helper(Node *root, unsigned int encode, int length);
double calc_speed(long input_size, double time);

int main (int argc, char *argv[])
{
    int ascii_letters[128] = { 0 };
    Node letters[128];
    long input_size;
    std::priority_queue<Node> tree;

    if ( argc != 2 ) {
        return 1;
    }
   
   //opening input file
    FILE *input = fopen(argv[1], "r");
    if (NULL == input) { 
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    
    //get input file size
    fseek(input, 0, SEEK_END);
    input_size = ftell(input);
    fseek(input, 0, SEEK_SET);

    //record frequency of letters
    char line[MAXLENGTH];
    while (fgets(line, sizeof(line), input)) {
        record_letters(line, ascii_letters);
    }
    fseek(input, 0, SEEK_SET);

    //make tree and encodings array
    make_queue(tree, ascii_letters);
    Node* root = build_tree(tree);
    assign_encode(root);
    int encodings[128][2];
    store_encodings(root, encodings);

    //timer for compression + compress function
    clock_t start = clock();
    compress(input, encodings, argv[2], ascii_letters);
    clock_t end = clock();

    //calculate compressions speed
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    double comp_speed = calc_speed(input_size, time);
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

void compress(FILE *input, int encodings[128][2], const char* output_filename, int ascii_letters[128]) {
    FILE *output = fopen(output_filename, "wb");
    if (output == NULL) {
        printf("Error opening output file\n");
        return;
    }

    //use input file size as upper bound for data malloc
    //what if file is very large and allocation fails: allocate as much as you can and then 
    //only measure compression
    //take buffer, give compressed buffer
    //break up function to only do compression

    //store the ascii_letters array
    unsigned char array_size = 128; 
    fwrite(&array_size, 1, 1, output);  // write size as first byte
    fwrite(ascii_letters, sizeof(int), 128, output);  // write ascii_letters array next

    //allocate buffer memory
    char *data = (char *)malloc(MAXLENGTH); //buffer 
    if (data == NULL) {
        printf("Error allocating memory\n");
        fclose(output);
        return;
    }
    int data_index = 0; //tracks current position in buffer
    unsigned int buffer = 0; //stores bits before writing to buffer
    char line[MAXLENGTH]; //buffer that holds each lineof text in input 
    int bits_in_buffer = 0; //counts bits in buffer

    //put this in other routine, anything operating on file should be moved
    //read data in main
    //return buffer and then call compress on buffer

    //read input line by line
    while (fgets(line, sizeof(line), input)) { //call fread only once
        int index = 0;
        char letter;

        //read each char in line 
        while (line[index] != '\0') {
            letter = line[index];
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
                data[data_index++] = byte; //adds a byte to buffer
                bits_in_buffer -= 8; //update bits_in_buffer

                //write entire buffer to file
                if (data_index >= MAXLENGTH) {
                    fwrite(data, 1, data_index, output);
                    data_index = 0;
                }
            }
            index++;
        }
    }

    //handles last byte; no extra 0's
    if (bits_in_buffer > 0) {
        unsigned char byte = buffer << (8 - bits_in_buffer);
        data[data_index++] = byte;
        data[data_index++] = bits_in_buffer; //stores num valid bits in last byte
    }
    else {
        data[data_index++] = 8;
    }
     
    //writes the remaining buffer to data
    if (data_index > 0) {
        fwrite(data, 1, data_index, output);
    }

    free(data);
    fclose(output);
}

/* note: want to avoid code duplication for node traversal, write helper */
void decompress(const char* input_filename, Node* root, const char* output_filename) { //change signature, give input buffer, decompress, then return decompressed buffer
    FILE *input = fopen(input_filename, "rb");
    if (input == NULL) {
        printf("Error opening input file\n");
        return;
    }
 
    //move file operation outside


    //counter for file size, stop when number of chars reached, remove valid_bits


    /* get last file byte, which stores number of valid bytes in last encode byte */
    fseek(input, -1, SEEK_END); //go to last byte of file
    long file_size = ftell(input); // get file size
    unsigned char valid_bits;
    fread(&valid_bits, 1, 1, input); // read num in last byte
    fseek(input, 0, SEEK_SET); // back to the beginning of file

    FILE *output = fopen(output_filename, "w");
    if (output == NULL) {
        printf("Error opening output file\n");
        fclose(input);
        return;
    }

    // read ascii_letters array size and skip to encode part of file
    //dont need to write size as first byte, can calculate
    unsigned char array_size;
    fread(&array_size, 1, 1, input);
    fseek(input, array_size * sizeof(int), SEEK_CUR);

    unsigned int buffer = 0; //buffer 
    int bits_in_buffer = 0; //counts bits in buffer
    unsigned char byte;
    Node* current = root; 

    //read byte by byte from input file until last 2 bytes
    while (ftell(input) < file_size - 1) {
        fread(&byte, 1, 1, input);
        buffer = (buffer << 8) | byte; //add a byte to buffer
        bits_in_buffer += 8;

        //read bit by bit in buffer
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

            //move outside routine
            if (current->left == NULL && current->right == NULL) {
                fputc(current->letter_name, output);
                current = root; //reset tree
            }
        }
    }

/* reads rest of input, making sure to stop at the right number of valid bites in 
last byte, as found above in valid_bits */
fread(&byte, 1, 1, input);
buffer = (buffer << 8) | byte;
bits_in_buffer += valid_bits;
while (bits_in_buffer > 0) {
    int bit = (buffer >> (bits_in_buffer - 1)) & 1;
    bits_in_buffer--;

    if (bit == 0) {
        current = current->left;
    } else {
        current = current->right;
    }

    if (current->left == NULL && current->right == NULL) {
        fputc(current->letter_name, output);
        current = root;
    }
}


    fclose(input);
    fclose(output);
}


double calc_speed(long original_size, double compression_time) {
    if (compression_time == 0) {
        return 0;
    }
    return original_size / compression_time;
}


// void compress(FILE *input, Node* root, const char* output_filename) {
//     FILE *output = fopen(output_filename, "wb");
//     if (output == NULL) {
//         printf("Error opening file\n");
//         return;
//     }

//     char line[MAXLENGTH];
//     unsigned int buffer = 0; 
//     int buffer_index = 0;

//     while (fgets(line, sizeof(line), input)) {
//         int index = 0;
//         char letter;
//         Node* node;

//         while (line[index] != '\0') {
//             letter = line[index];
//             node = find_node(root, letter); //too slow, use another method like 2D array, store encode and length
//             //write seperate function to store

//             //traversing tree can be inefficient
//             //unisgned int a['a'] =...

//             if (node) {
//                 int i;
//                 for (i = 0; i < node->encode_length; i++) {
//                     buffer <<= 1;
//                     buffer |= (node->encode >> (node->encode_length - i - 1)) & 1;
//                     buffer_index++;
//                     if (buffer_index == node->encode_length - 1) {
//                         fwrite(&buffer, 1, 1, output);
//                         buffer = 0;
//                         buffer_index = 0;
//                     }
//                 }
//             }
//             index++;
//         }
//     }

//     if (buffer_index > 0) {
//         buffer <<= (8 - buffer_index);
//         fwrite(&buffer, 1, 1, output);
//     }

//     fclose(output);
// }


// void decompress(const char* input_filename, Node* root, const char* output_filename) {
//     FILE *input = fopen(input_filename, "rb");
//     FILE *output = fopen(output_filename, "w");
//     if (input == NULL || output == NULL) {
//         printf("Error opening file\n");
//         return;
//     }

//     unsigned char buffer;
//     int buffer_index = 0;
//     Node* current = root;

//     while (fread(&buffer, 1, 1, input) == 1) {
//         int i;
//         for (i = 7; i >= 0; i--) {
//             if (buffer & (1 << i)) {
//                 current = current->right;
//             } else {
//                 current = current->left;
//             }

//             if (current->left == NULL && current->right == NULL) {
//                 fputc(current->letter_name, output);
//                 current = root;
//             }
//         }
//     }

//     fclose(input);
//     fclose(output);
// }




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
