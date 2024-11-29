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

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "error in arguments", argv[0]);
        return 1;
    }

    // Opening compressed input file
    FILE *input = fopen(argv[1], "rb");
    if (input == NULL) { 
        perror("Error opening input file");
        return 1;
    }

    // Determine the size of the compressed file
    fseek(input, 0, SEEK_END); // Move to the end of the file
    long file_size = ftell(input); // Get current position (file size)
    fseek(input, 0, SEEK_SET); // Reset to beginning of the file

    printf("file size = %lu\n", file_size);

    // Allocate memory for the entire compressed buffer
    char *compressed_buffer = (char *) malloc(file_size);
    if (compressed_buffer == NULL) {
        perror("Error allocating buffer for compressed data");
        fclose(input);
        return 1;
    }

    // Read the entire compressed file into the buffer
    size_t read_bytes = fread(compressed_buffer, 1, file_size, input);
    if (read_bytes != file_size) {
        fprintf(stderr, "Error reading compressed file\n");
        free(compressed_buffer);
        fclose(input);
        return 1;
    }
    fclose(input);

    size_t num_elements;
    memcpy(&num_elements, compressed_buffer, sizeof(size_t));

    float* dest = (float *) malloc (num_elements * sizeof(float)); // Over-allocate
    if (dest == NULL) {
        perror("Error allocating buffer for decompressed data");
        free(compressed_buffer);
        return 1;
    }

    int destsize = 0; // To hold the size of decompressed data

    // Perform decompression
    int decompress_status = decompress(compressed_buffer, dest, file_size, &destsize);
    if (decompress_status != 0) {
        fprintf(stderr, "Decompression failed\n");
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    // Replace writefile with direct fwrite
    FILE* output = fopen(argv[2], "wb");
    if (output == NULL) {
        perror("Error opening output file");
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    size_t written = fwrite(dest, sizeof(float), destsize, output);
    if (written != (size_t)destsize) {
        fprintf(stderr, "Error writing to output file\n");
        fclose(output);
        free(compressed_buffer);
        free(dest);
        return 1;
    }

    fclose(output);

    // Free allocated memory
    free(compressed_buffer);
    free(dest);

    return 0;
}
