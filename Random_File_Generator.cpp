#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_file(const char* filename, long size);

int main() {
    const char* filename = "test_input.txt";
    long mb = 512; 
    long size = mb * 1024 * 1024;

    generate_file(filename, size);
    printf("generated: %s\n", filename);
    return 0;
}

void generate_file(const char* filename, long size) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing\n");
        exit(EXIT_FAILURE);
    }

    srand((unsigned int)time(NULL));
    long bytes_written = 0;

    while (bytes_written < size) {
        int ascii_value = rand() % 128;  // Generate random ASCII value
        int repeat_count = rand() % 100 + 1;  // Generate random repeat count between 1 and 100

        for (int i = 0; i < repeat_count; i++) {
            if (bytes_written < size) {
                fputc(ascii_value, file);
                bytes_written++;
            } else {
                break;
            }
        }
    }

    fclose(file);
}
