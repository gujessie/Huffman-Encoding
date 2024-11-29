#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <functional>
#include <queue>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <chrono>
#include "Floats_Huffman_Encoding.h"

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


    //HOMEWORK
    //test larger inputs, 100 MB!!!!!!!!!!!!!
    //floating point data: data generator, generate random floating point data and write to file
    //  quantization: try to write a quantization routine
    //  divide min an max into buckets between the interval (btwn min and max) and then do quantization
    //  in quantization step, convert all values in the bucket into the bucket number
    //  for example, convert all numbers btwn 0 and 1 into 0 , 1 and 2 into 1
    //  dont want bucket numbers to be chars s change huffman encoding to work with integers and not characters
    //  change to work w integers 
    //  also change file generator 
    //  input data: floating point data, and will be converted to a series of bucket numbers, then use huffman to encode bucket numbers

    //makefile go back to -o0, makes it easier to debug
    //record frequency of letters
    //read all of data, no fgets, 


    //intitialize src buffer

    char *src = (char *)malloc(fsize);
    if (NULL == src) {
        printf("Error allocating memory for src\n");
        fclose(input);
        return 1;
    }
    fread(src, 1, fsize, input);
    fclose(input);

    record_letters(src, ascii_letters);

    fseek(input, 0, SEEK_SET);

    //.....

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

    //intialize dest buffer
    char *dest = (char *)malloc(fsize); //create buffer 
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

    char *compressed_data = (char *)malloc(compressed_size);
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
