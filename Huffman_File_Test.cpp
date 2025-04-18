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

struct node {
    struct node *left;
    struct node *right;
    int frq;
    char letter_name;
    unsigned int encode;
    int encode_length;

    bool operator<(const struct node& comp_node) const {
        return frq > comp_node.frq;
    }
};
typedef struct node Node;

void record_letters(const char *string, int ascii_letters[]);
void make_queue(std::priority_queue<Node> &tree, int *ascii_letters);
Node* build_tree(std::priority_queue<Node> &tree);
void assign_encode_helper(Node *root, unsigned int encode, int length);
void assign_encode(Node *root);
void compress(const char *input, int encodings[128][2], const char* output_filename);
void decompress(const char* input_filename, Node* root, const char* output_filename);
void free_tree(Node* root);
void store_encodings(Node* root, int encodings[128][2]);
void store_encodings_helper(Node* root, int encodings[128][2]);
void print_encodings(int encodings[128][2]);

int main(void) {
    const char *test_string = "aaaabbbbaaaccddddee";
    int ascii_letters[128] = {0};
    std::priority_queue<Node> tree;

    record_letters(test_string, ascii_letters);
    make_queue(tree, ascii_letters);
    Node* root = build_tree(tree);
    assign_encode(root);

    int encodings[128][2] = {0};
    store_encodings(root, encodings);
    
    print_encodings(encodings);

    // Compress
    compress(test_string, encodings, "compressed.bin");

    // Print compressed content in base two
    printf("compressed: ");
    FILE *comp_file = fopen("compressed.bin", "rb");
    if (comp_file != NULL) {
        unsigned char byte;
        while (fread(&byte, 1, 1, comp_file) == 1) {
            for (int i = 7; i >= 0; i--) {
                printf("%d", (byte >> i) & 1);
            }
            printf(" ");
        }
        fclose(comp_file);
        printf("\n");
    } else {
        printf("Error opening compressed file\n");
    }

    // Decompress
    decompress("compressed.bin", root, "decompressed.txt");

    // Print decompressed content
    FILE *output = fopen("decompressed.txt", "r");
    if (output != NULL) {
        char result[MAXLENGTH];
        fgets(result, sizeof(result), output);
        printf("Decompressed: %s\n", result);
        fclose(output);
    } else {
        printf("Error opening decompressed file\n");
    }

    free_tree(root);

    return 0;
}

void record_letters(const char *string, int ascii_letters[]) {
    int index = 0;
    while (string[index] != '\0') {
        char letter = string[index];
        int ascii_value = (int)letter;
        ascii_letters[ascii_value]++;
        index++;
    }
}

void make_queue(std::priority_queue<Node> &tree, int *ascii_letters) {
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

Node* build_tree(std::priority_queue<Node> &tree) {
    while (tree.size() > 1) {
        Node* left_node = (Node*)malloc(sizeof(Node));
        if (!left_node) exit(EXIT_FAILURE);
        *left_node = tree.top();
        tree.pop();

        Node* right_node = (Node*)malloc(sizeof(Node));
        if (!right_node) {
            free(left_node);
            exit(EXIT_FAILURE);
        }
        *right_node = tree.top();
        tree.pop();

        Node* new_node = (Node*)malloc(sizeof(Node));
        if (NULL == new_node) {
            free(left_node);
            free(right_node);
            exit(EXIT_FAILURE);
        }
        new_node->letter_name = '\0';
        new_node->frq = left_node->frq + right_node->frq;
        new_node->left = left_node;
        new_node->right = right_node;
        new_node->encode = 0;
        new_node->encode_length = 0;
        tree.push(*new_node);
    }

    Node* root = (Node*)malloc(sizeof(Node));
    if (root == NULL) exit(EXIT_FAILURE);
    *root = tree.top();
    tree.pop();
    return root;
}

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
    if (root->left) free_tree(root->left);
    if (root->right) free_tree(root->right);
    free(root);
}

void store_encodings(Node* root, int encodings[128][2]) {
    for (int i = 0; i < 128; ++i) {
        encodings[i][0] = 0;
        encodings[i][1] = 0;
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

void print_encodings(int encodings[128][2]) {
    printf("Encodings:\n");
    for (int i = 0; i < 128; ++i) {
        if (encodings[i][1] > 0) {
            printf("%c: ", (char)i);
            for (int j = encodings[i][1] - 1; j >= 0; --j) {
                printf("%d", (encodings[i][0] >> j) & 1);
            }
            printf("\n");
        }
    }
}

void compress(const char *input, int encodings[128][2], const char* output_filename) {
    FILE *output = fopen(output_filename, "wb");
    if (output == NULL) {
        printf("Error opening output file\n");
        return;
    }

    unsigned int buffer = 0;
    int bits_in_buffer = 0;
    int index = 0;

    while (input[index] != '\0') {
        char letter = input[index];
        int encode = encodings[(int)letter][0];
        int length = encodings[(int)letter][1];

        buffer = (buffer << length) | encode;
        bits_in_buffer += length;

        while (bits_in_buffer >= 8) {
            unsigned char byte = (buffer >> (bits_in_buffer - 8)) & 0xFF;
            fwrite(&byte, 1, 1, output);
            bits_in_buffer -= 8;
        }
        index++;
    }

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

    unsigned int buffer = 0;
    int bits_in_buffer = 0;
    unsigned char byte;
    Node* current = root;

    while (fread(&byte, 1, 1, input) == 1) {
        buffer = (buffer << 8) | byte;
        bits_in_buffer += 8;

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
    }

    fclose(input);
    fclose(output);
}
