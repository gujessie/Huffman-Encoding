#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ncurses.h>
#include <ctype.h>
#include <time.h>
//define a constant
#define MAXLENGTH 100000

int calculate_length(char* data);
char* compress(char* data);
char* decompress(char* data);
double calculate_compratio(char* data, char* compdata);
double calculate_throughput(char* data, double time);

int main(void)
{
    char c_d;
    double compratio = 0;
    char y_n;
    clock_t start, end;
    double time;

    //take file input
    FILE *input;
    char input_filename[MAXLENGTH + 1], output_filename[MAXLENGTH + 1];

    //no prompts in the future, use argc, argv and pass to main funtion
    //check for return types for every function
    printf("Enter input file name: ");
    scanf("%s", input_filename);
    input = fopen(input_filename, "r");
    if (NULL == input) { //make sure to write null before value compared, makes easier to debug if only one = is written
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
        return 1;
    }
    //read input into a string
    fseek(input, 0, SEEK_END);
    long fsize = ftell(input);
    fseek(input, 0, SEEK_SET);  


    char *stringtemp = malloc(fsize + 1);
    //check return type, make sure stringtemp is not null, is so malloc has failed
    fread(stringtemp, fsize, 1, input);
    fclose(input);
    stringtemp[fsize] = '\0';

    //create file ouput
    FILE * output;
    printf("Provide output file name: ");
    scanf("%s", output_filename);
    output = fopen(output_filename, "w");
    if (NULL == output) {
        printf("Error opening file\n");
        exit(EXIT_FAILURE);
        return 1;
    }

    printf("\nWould you like to compress or decompress data: Type c or d\n");
    scanf(" %c", &c_d);
    if(c_d == 'c' || c_d =='C'){
    //    char data[MAXLENGTH + 1];
    //    printf("Enter data to compress: ");
    //    scanf("%s", data);
       start = clock();
       char* compdata = compress(stringtemp);
       end = clock();
       time = ((double) (end-start)) / CLOCKS_PER_SEC;

       fprintf(output, "%s \n", compdata); 
       fclose(output);

        // look into #ifdef, compiler pragma
       printf("Would you like to calculate the compression ratio? Type y or n\n");
       scanf(" %c", &y_n);
       if(y_n == 'y' || y_n =='Y'){
          compratio = calculate_compratio(stringtemp, compdata);
          printf("%f", compratio);
       }

       double throughput = 0;
       printf("\nWould you like to calculate the throughput? Type y or n\n");
       scanf(" %c", &y_n);
       if(y_n == 'y' || y_n == 'Y'){
          throughput = calculate_throughput(stringtemp, time);
          printf("The throughput is: %f MB per second", throughput);
       }
    }
    if(c_d == 'd' || c_d =='D'){
    //    char data[MAXLENGTH + 1];
    //    printf("Enter data to decompress: ");
    //    scanf("%s", data);
       const char* decompdata = decompress(stringtemp);
       fprintf(output, "%s", decompdata);
    }
    free(stringtemp);
    return 0;

}

char* compress(char* data){
    int index = 0, rindex = 0;
    char *result;
    int datalength = strlen(data);

    //allocates size
    result = malloc(sizeof(char)* (2 * datalength + 1));
    if(NULL == data){
        printf("Error processing compression");
        return NULL;
    }
    int length = calculate_length(data);
    while (index<length){
        int count = 0;
        char letter = data[index];

        while(data[index] == letter){
            count++;
            index++;
        }
        strcat(result, &letter);
        rindex++;

        //make portable
        char countstr[12]; 
        snprintf(countstr, sizeof(countstr), "%d", count);
        for (int i = 0; countstr[i] != '\0'; i++) {
            result[rindex++] = countstr[i];
        }
        // char buf[2];
        // snprintf(buf, 2, "%d", count);
        // strcat(result, buf);
    }
    result[rindex] = '\0';
    //check return type
    result = realloc(result, sizeof(char) * (rindex + 1));
    return result;
}

char* decompress(char* data){
    int index = 0, rindex = 0, countindex = 0, count;
    char* result;
    
    result = malloc(sizeof(char)*(MAXLENGTH + 1));

    if(data == NULL){
        printf("Error processing decompression");
        return NULL;
    }
    // int length = calculate_length(data);
    // while(index < length){
    //     //memset to get rid of loop
    //     char letter[] = {data[index], '\0'};
    //     int count = data[index + 1] - '0';
    //     int i = 0;
    //     for(i=0; i<count; i++){
    //         strcat(result, letter);
    //     }
    //     index = index + 2;
    // }
    // int resultlength = calculate_length(result);
    // result[resultlength] = '\0';
    // return result; 
    int length = calculate_length(data);
    while (index < length) {
        char letter = data[index];
        index++;
        
        //gets number of times to repeat letter
        char countstr[12]; 

        while (isdigit(data[index])) {
            countstr[countindex] = data[index];
            countindex++;
            index++;
        }
        countstr[countindex] = '\0';
        count = atoi(countstr);

        // memset to replace loop
        memset(result + rindex, letter, count);
        rindex += count;
    }
    result[rindex] = '\0';
    return result;
}

double calculate_compratio(char* data, char* compdata){
    int datalength = calculate_length(data);
    int complength = calculate_length(compdata);
    //checks parameters
    if (datalength == 0 || complength == 0){
        printf("Error: invalid parameter provided");
        return 0;
    }
    double ratio = (double)datalength/(double)complength;
    return ratio;
}

double calculate_throughput(char* data, double time){
    int datalength = calculate_length(data);
    //checks parameters
    if (datalength == 0 || time == 0){
        printf("Error: invalid parameter provided");
        return 0;
    }
    double throughput = (double)datalength/time;
    throughput = (throughput/1024)/1024;
    return throughput;
}

int calculate_length(char* data){
    int index = 0;
    int count = 0;
    while(data[index] != 0){
        count++;
        index++;
    }
    return count;
}
