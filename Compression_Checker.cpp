#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Floats_Huffman_Encoding.h"
#define ERROR 0.005  //error

bool compare_arrays(float* arr1, float* arr2, int length) {
    for (int i = 0; i < length; i++) {
        if (fabs(arr1[i] - arr2[i]) > ERROR) {
            return false;
        }
    }
    return true;
}

int main() {
    // Sample input array
    float input[] = {1.23, 4.56, 7.89, 3.14, 2.71}; // Example input for compression
    int num_elements = sizeof(input) / sizeof(input[0]);

    // Variables for compression
    const char* compressed_file = "compressed.bin";

    // Compress the input (finish compression)
    // compress
    // void compress_quantization_levels(const int* src, char* dest, int fsize, int* destsize, in
    compress_quantization_levels(input, compressed_file, num_elements);

#if 0
    // Decompressed array to hold results
    float decompressed[num_elements];
    float min = 1.23;  // Replace with actual minimum used in your quantization
    float bucket_size = 0.01;  // Replace with the actual bucket size used

    // Decompress the compressed file
    decompress(compressed_file, decompressed, num_elements, min, bucket_size);

    // Compare original and decompressed data
    if (compare_arrays(input, decompressed, num_elements)) {
        printf("Decompression is correct.\n");
    } else {
        printf("Decompression is incorrect.\n");
    }
#endif
    return 0;
}

