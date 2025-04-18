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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    // Opening input file
    FILE* input = fopen(argv[1], "rb"); // Changed to "rb" for binary read
    if (NULL == input) {
        perror("Error opening input file");
        return 1;
    }

    // Get input file size
    fseek(input, 0, SEEK_END);
    size_t num_elements = ftell(input) / sizeof(float); // Number of elements in file
    fseek(input, 0, SEEK_SET);

    printf("num_elements: %ld\n", num_elements);

    // Allocate memory for source data
    float* src = (float*)malloc(num_elements * sizeof(float));
    if (NULL == src) {
        printf("Error allocating memory for src\n");
        fclose(input);
        return 1;
    }

    // Read data from input file
    size_t read_elements = fread(src, sizeof(float), num_elements, input);
    if (read_elements != num_elements) {
        printf("Error reading input file\n");
        fclose(input);
        free(src);
        return 1;
    }
    fclose(input);

    char* dest = (char*)malloc(MAXLENGTH);
    if (NULL == dest) {
        perror("Error allocating memory for dest");
        free(src);
        return 1;
    }

    // Define quantization error
    const float error = 0.25f;

    // Variable to hold compression metadata
    int destsize = 0;

    // Perform compression
    int compress_status = compress(src, dest, num_elements, &destsize, error);
    if (compress_status != 0) {
        printf("Compression failed\n");
        free(src);
        free(dest);
        return 1;
    }

    // Write to output file
    writefile(argv[2], dest, destsize);

    // Free allocated memory
    free(src);
    free(dest);

    return 0;
}
