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

void subtract4D(float* map1, int* nhit1, float* rms1, 
                float* map2, int* nhit2, float* rms2,
                float* map,  int* nhit,  float* rms,
                int n0,      int n1,     int n2, 
                int n3                
                ){
    int prod = n2 * n3;
    int prod1 = prod * n1;
    int idx;
    float var1, var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    idx = prod1 * i + prod * j + n3 * k + l;
                    if (nhit1[idx] * nhit2[idx] > 0){
                        var1  = rms1[idx] * rms1[idx];
                        var2  = rms2[idx] * rms2[idx];
                        map[idx]  = map1[idx] - map2[idx];

                        nhit[idx] = nhit1[idx] + nhit2[idx];
                        rms[idx]  =  sqrt(var1 + var2);
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

void subtract5D(float* map1, int* nhit1, float* rms1, 
             float* map2, int* nhit2, float* rms2,
             float* map,  int* nhit,  float* rms,
             int n0,      int n1,     int n2, 
             int n3,      int n4                
             ){
    int prod = n3 * n4;
    int prod1 = prod * n2;
    int prod2 =  prod1 * n1;
    int idx;
    float var1, var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    for(int m=0; m < n4; m++){
                        idx = prod2 * i + prod1 * j + prod * k + n4 * l + m;
                        if (nhit1[idx] * nhit2[idx] > 0){
                            var1        = rms1[idx] * rms1[idx];
                            var2        = rms2[idx] * rms2[idx];
                            map[idx]    = map1[idx] - map2[idx];

                            nhit[idx]   = nhit1[idx] + nhit2[idx];
                            rms[idx]    =  sqrt(var1 + var2);
                        }
                        else {
                            map[idx]    = 0;
                            nhit[idx]   = 0;
                            rms[idx]    = 0;
                        }       
                    }
                }
            }
        }
    }
}

void subtract6D(float* map1, int* nhit1, float* rms1, 
             float* map2, int* nhit2, float* rms2,
             float* map,  int* nhit,  float* rms,
             int n0,      int n1,     int n2, 
             int n3,      int n4,     int n5){
    int prod = n4 * n5;
    int prod1 = prod * n3;
    int prod2 =  prod1 * n2;
    int prod3 =  prod2 * n1;
    int idx;
    float var1, var2;

    for(int i=0; i < n0; i++){
        for(int j=0; j < n1; j++){
            for(int k=0; k < n2; k++){
                for(int l=0; l < n3; l++){
                    for(int m=0; m < n4; m++){
                        for(int n=0; n < n5; n++){
                            idx = prod3 * i + prod2 * j + prod1 * k + prod * l + n5 * m + n;
                            if (nhit1[idx] * nhit2[idx] > 0){
                                var1  = rms1[idx] * rms1[idx];
                                var2  = rms2[idx] * rms2[idx];
                                map[idx]  = map1[idx] - map2[idx];

                                nhit[idx] = nhit1[idx] + nhit2[idx];
                                rms[idx]  =  sqrt(var1 + var2);
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

void dgrade4D(float* map_h, int* nhit_h, float* rms_h, 
              float* map_l, int* nhit_l, float* rms_l,
              int n0,       int n1,      int n2, 
              int n3,       int N2,      int N3,       
              int num){

    int p_h = n2 * n3;
    int p_h1 = p_h * n1;
    int idx_h;
    
    int p_l = N2 * N3;
    int p_l1 = p_l * n1;
    int idx_l;

    float inv_var1, inv_var1_sum;
    int i, j, k ,l, a, b;

    for(i=0; i < n0; i++){
        for(j=0; j < n1; j++){
            for(a = 0; a < N2; a++){
                for(b = 0; b < N3; b++){
                    idx_l = p_l1 * i + p_l * j + N3 * a + b;;
                    inv_var1_sum = 0;
                    
                    for(k = a * num; k < (a + 1) * num; k++){
                        for(l = b * num; l < (b + 1) * num; l++){
                            idx_h = p_h1 * i + p_h * j + n3 * k + l;
                            if (nhit_h[idx_h] > 0){
                                inv_var1        = 1.0 / (rms_h[idx_h] * rms_h[idx_h]);
                                map_l[idx_l]     += map_h[idx_h] * inv_var1;
                                nhit_l[idx_l]    += nhit_h[idx_h];
                            }
                            else{
                                inv_var1    = 0;
                                map_l[idx_l]  = 0;
                            }
                            inv_var1_sum   += inv_var1;
                        }
                    }
                    if (nhit_l[idx_l] > 0){
                        map_l[idx_l] /= inv_var1_sum;
                        rms_l[idx_l] = 1 / sqrt(inv_var1_sum);      
                    }
                    else{
                        map_l[idx_l] = 0;
                        rms_l[idx_l] = 0;
                    }
                }      
            }
        }
    }
}

void dgrade5D(float* map_h, int* nhit_h, float* rms_h, 
              float* map_l, int* nhit_l, float* rms_l,
              int n0, int n1, int n2,
              int n3, int n4, int N3,
              int N4, int num){
    int p_h = n3 * n4;
    int p_h1 = p_h * n2;
    int p_h2 =  p_h1 * n1;
    int idx_h;
    
    int p_l = N3 * N4;
    int p_l1 = p_l * n2;
    int p_l2 =  p_l1 * n1;
    int idx_l;

    float inv_var1, inv_var1_sum;
    int i, j, k ,l, m, a, b;

    for(i=0; i < n0; i++){
        for(j = 0; j < n1; j++){
            for(k=0; k < n2; k++){
                for(a = 0; a < N3; a++){
                    for(b = 0; b < N4; b++){
                            
                        idx_l = p_l2 * i + p_l1 * j + p_l * k + N4 * a + b;;
                        inv_var1_sum = 0;
                        
                        for(l = a * num; l < (a + 1) * num; l++){
                            for(m = b * num; m < (b + 1) * num; m++){
                                idx_h = p_h2 * i + p_h1 * j + p_h * k + n4 * l + m;
                                
                                if (nhit_h[idx_h] > 0){
                                    inv_var1        = 1.0 / (rms_h[idx_h] * rms_h[idx_h]);
                                    map_l[idx_l]     += map_h[idx_h] * inv_var1;
                                    nhit_l[idx_l]    += nhit_h[idx_h];
                                }
                                else{
                                    inv_var1    = 0;
                                    map_l[idx_l]  = 0;
                                }
                                inv_var1_sum   += inv_var1;
                            }
                        }
                        if (nhit_l[idx_l] > 0){
                            map_l[idx_l] /= inv_var1_sum;
                            rms_l[idx_l] = 1 / sqrt(inv_var1_sum);      
                        }
                        else{
                            map_l[idx_l] = 0;
                            rms_l[idx_l] = 0;
                        }
                    }      
                }
            }
        }
    }
}

void dgrade6D(float* map_h, int* nhit_h, float* rms_h, 
              float* map_l,  int* nhit_l,  float* rms_l,
              int n0,      int n1,     int n2, 
              int n3,      int n4,     int n5,
              int N4,      int N5,     int num){
    int p_h = n4 * n5;
    int p_h1 = p_h * n3;
    int p_h2 =  p_h1 * n2;
    int p_h3 =  p_h2 * n1;
    int idx_h;

    int p_l = N4 * N5;
    int p_l1 = p_l * n3;
    int p_l2 = p_l1 * n2;
    int p_l3 = p_l2 * n1;
    int idx_l;
    float inv_var1, inv_var1_sum;
    int i, j, k ,l, m, n, a, b;
    for(i=0; i < n0; i++){
        for(j=0; j < n1; j++){
            for(k=0; k < n2; k++){
                for(l=0; l < n3; l++){
                    for(a = 0; a < N4; a++){
                        for(b = 0; b < N5; b++){
                            
                            idx_l = p_l3 * i + p_l2 * j + p_l1 * k + p_l * l + N5 * a + b;

                            for(m = a * num; m < (a + 1) * num; m++){
                                for(n = b * num; n < (b + 1) * num; n++){
                                    idx_h = p_h3 * i + p_h2 * j + p_h1 * k + p_h * l + n5 * m + n;
                                    if (nhit_h[idx_h] > 0){
                                        inv_var1        = 1.0 / (rms_h[idx_h] * rms_h[idx_h]);
                                        map_l[idx_l]     += map_h[idx_h] * inv_var1;
                                        nhit_l[idx_l]    += nhit_h[idx_h];
                                    }
                                    else{
                                        inv_var1    = 0;
                                        map_l[idx_l]  = 0;
                                    }
                                    inv_var1_sum   += inv_var1;
                                }
                            }
                            if (nhit_l[idx_l] > 0){
                                map_l[idx_l] /= inv_var1_sum;
                                rms_l[idx_l] = 1 / sqrt(inv_var1_sum);      
                            }
                            else{
                                map_l[idx_l] = 0;
                                rms_l[idx_l] = 0;
                            }
                        }      
                    }
                }
            }
        }
    }
}
