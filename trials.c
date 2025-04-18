#include<stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void){
int i = 6;
int* p;
p = &i;
printf("%d %d\n", *p, *(&i));
printf("%d", *p);
}
