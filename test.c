# include <stdio.h>

void test(double *in_array, int n, int m){
    int i, j;
    for(i = 0; i<n; i++) {
            for(j = 0; j<m; j++) {
                printf("%e \t", in_array[i][j]);
        }
        printf("\n");
    }    
}
/*
int test(float *x, int n){
    int i;
    int count = 0;
    for (i = 0; i < n; i++){
        count += x[i];
        printf("%g  ", x[i]);
    }
    return count;
}
*/