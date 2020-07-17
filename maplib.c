// Compile as
// gcc -shared -o maplib.so.1 maplib.c
#include <stdio.h>
void makemaps(double* mapout, int nx, int ny){
    printf("hello\n");
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            printf("%d%d\n", i, j);
            int idx = i*ny + j;
            mapout[idx] += i + 0.5*j;
        }
    }
}