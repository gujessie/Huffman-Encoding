#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_file(const char* filename, long size);

int main() {
    const char* filename = "test_input.bin"; 
    long mb = 512;
    long size = mb * 1024 * 1024; 

    generate_file(filename, size);
    printf("Generated: %s\n", filename);
    return 0;
}

void generate_file(const char* filename, long size) {
    FILE *file = fopen(filename, "wb"); 
    if (file == NULL) {
        printf("Error opening file for writing\n");
        exit(EXIT_FAILURE);
    }

    srand((unsigned int)time(NULL));
    long bytes_written = 0;


    long num_floats = size / sizeof(float);
    float value;

    while (bytes_written < size) {

        value = (float)rand() / RAND_MAX * 100.0f; 


        if (fwrite(&value, sizeof(float), 1, file) != 1) {
            printf("Error writing to file\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        bytes_written += sizeof(float);
    }

    fclose(file);
}
