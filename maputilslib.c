// Compile as
// gcc -shared maputilslib.c -o maputilslib.so.1
#include <stdio.h>
#include <time.h>
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
}