// Compile as
// gcc -shared maputilslib.c -o maputilslib.so.1
#include <stdio.h>
#include <time.h>
#include <math.h>

void coadd4D(float* map1, int* nhit1, float* rms1, 
             float* map2, int* nhit2, float* rms2,
             float* map,  int* nhit,  float* rms,
             int n0,      int n1,     int n2, 
             int n3                
             ){
    int prod = n2 * n3;
    int prod1 = prod * n1;
    int idx;
    float inv_var1, inv_var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    idx = prod1 * i + prod * j + n3 * k + l;
                    if (nhit1[idx] * nhit2[idx] > 0){
                        inv_var1  = 1.0 / (rms1[idx] * rms1[idx]);
                        inv_var2  = 1.0 / (rms2[idx] * rms2[idx]);
                        map[idx]  = map1[idx] * inv_var1 + map2[idx] * inv_var2;
                        map[idx] /=  inv_var1 + inv_var2;

                        nhit[idx] = nhit1[idx] + nhit2[idx];
                        rms[idx] = 1 / sqrt(inv_var1 + inv_var2);
                    }
                    else {
                        map[idx] = 0;
                        nhit[idx] = 0;
                        rms[idx] = 0;
                    }       
                }
            }
        }
    }
}


void coadd5D(float* map1, int* nhit1, float* rms1, 
             float* map2, int* nhit2, float* rms2,
             float* map,  int* nhit,  float* rms,
             int n0,      int n1,     int n2, 
             int n3,      int n4                
             ){
    int prod = n3 * n4;
    int prod1 = prod * n2;
    int prod2 =  prod1 * n1;
    int idx;
    float inv_var1, inv_var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    for(int m=0; m < n4; m++){
                        idx = prod2 * i + prod1 * j + prod * k + n4 * l + m;
                        if (nhit1[idx] * nhit2[idx] > 0){
                            inv_var1  = 1.0 / (rms1[idx] * rms1[idx]);
                            inv_var2  = 1.0 / (rms2[idx] * rms2[idx]);
                            map[idx]  = map1[idx] * inv_var1 + map2[idx] * inv_var2;
                            map[idx] /=  inv_var1 + inv_var2;

                            nhit[idx] = nhit1[idx] + nhit2[idx];
                            rms[idx] = 1 / sqrt(inv_var1 + inv_var2);
                        }
                        else {
                            map[idx] = 0;
                            nhit[idx] = 0;
                            rms[idx] = 0;
                        }       
                    }
                }
            }
        }
    }
}

void coadd6D(float* map1, int* nhit1, float* rms1, 
             float* map2, int* nhit2, float* rms2,
             float* map,  int* nhit,  float* rms,
             int n0,      int n1,     int n2, 
             int n3,      int n4,     int n5){
    int prod = n4 * n5;
    int prod1 = prod * n3;
    int prod2 =  prod1 * n2;
    int prod3 =  prod2 * n1;
    int idx;
    float inv_var1, inv_var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    for(int m=0; m < n4; m++){
                        for(int n=0; n < n5; n++){
                            idx = prod3 * i + prod2 * j + prod1 * k + prod * l + n5 * m + n;
                            if (nhit1[idx] * nhit2[idx] > 0){
                                inv_var1  = 1.0 / (rms1[idx] * rms1[idx]);
                                inv_var2  = 1.0 / (rms2[idx] * rms2[idx]);
                                map[idx]  = map1[idx] * inv_var1 + map2[idx] * inv_var2;
                                map[idx] /=  inv_var1 + inv_var2;

                                nhit[idx] = nhit1[idx] + nhit2[idx];
                                rms[idx] = 1 / sqrt(inv_var1 + inv_var2);
                            }
                            else {
                                map[idx] = 0;
                                nhit[idx] = 0;
                                rms[idx] = 0;
                            }
                        }      
                    }
                }
            }
        }
    }
 }


/*
void add5D(float* mapout, float* mapin1, float* mapin2, int n0, int n1,
              int n2, int n3, int n4){
    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    for(int m=0; m < n4; m++){
                        int idx = (n4 * (n3 * (n2 * (i*n1 + j) + k) + l) + m);
                        mapout[idx] = mapin1[idx] + mapin2[idx];
                    }
                }
            }
        }
    }
    
    clock_t t = clock();
    t = clock() - t;
    
    double time_taken = ((double) t)/ CLOCKS_PER_SEC;
    printf("Run time in C: %g sec\n", time_taken);
}
*/
    /*
    Python:
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        self.maputilslib.add5D.argtypes = [float32_array5, float32_array5, float32_array5,
                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                    ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4 = map1.shape
        map = np.zeros((n0, n1, n2, n3, n4), dtype = ctypes.c_float)
        start = time.time()
        self.maputilslib.add5D(map, map1, map2, n0, n1, n2, n3, n4)
        print("Run time outside C: ", time.time() - start, " sec")

        start = time.time()
        sum = map1 + map2
        print("Run time python/numpy: ", time.time() - start, " sec")
        print(np.allclose(map, sum))
        
    */