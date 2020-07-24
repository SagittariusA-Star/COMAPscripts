import numpy as np
import h5py as h5
import ctypes

"""
a = np.arange(6 * 6 * 6 * 6 * 6)
a = a.reshape(6, 6, 6, 6, 6)

num = 2
n0, n1, n2, n3, n4 = a.shape 
m3, m4 = int(n3 / num), int(n4 / num)
b = np.zeros((n0, n1, n2, m3, m4))

for q in range(n0):
    for h in range(n1):
        for p in range(n2):
            for i in range(m3):
                for j in range(m4):
                    for k in range(i * num, (i + 1)* num):
                        for l in range(j * num, (j + 1) * num):
                            b[q, h, p, i, j] += a[q, h, p, k, l]
                                
print("first")
print(a[0, 0, 0, :, :])
print(b[0, 0, 0, :, :]) 

m_h = np.arange(6 * 6 * 6 * 6 * 6)
m_h = m_h.reshape(6, 6, 6, 6, 6).astype(ctypes.c_float)

n_h = np.arange(6 * 6 * 6 * 6 * 6)
n_h = m_h.reshape(6, 6, 6, 6, 6).astype(ctypes.c_int)

r_h = np.arange(6 * 6 * 6 * 6 * 6)
r_h = m_h.reshape(6, 6, 6, 6, 6).astype(ctypes.c_float)



num = int(3)
n0, n1, n2, n3, n4 = m_h.shape
N3, N4 = int(n3 / num), int(n4 / num)

m_l = np.zeros((n0, n1, n2, N3, N4), dtype = ctypes.c_float)
n_l = np.zeros((n0, n1, n2, N3, N4), dtype = ctypes.c_int)
r_l = np.zeros((n0, n1, n2, N3, N4), dtype = ctypes.c_float)

maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library

float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
maputilslib.dgrade5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                 float32_array5, int32_array5, float32_array5,
                                 ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                 ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                 ctypes.c_int, ctypes.c_int]
maputilslib.dgrade5D(m_h, n_h, r_h,
                     m_l, n_l, r_l,
                     n0, n1, n2,
                     n3, n4, N3, N4, num)
print("second")
print(n_h[1, 2, 3, :, :])
print(n_l[1, 2, 3, :, :])
#print(np.allclose(n_h, a))
#print(np.allclose(n_l, b))
"""
"""
a = np.arange(6 * 6 * 6 * 6 * 6 * 6)
a = a.reshape(6, 6, 6, 6, 6, 6)

num = 2
n0, n1, n2, n3, n4, n5 = a.shape 
N3 = int(n3 / num)
b = np.zeros((n0, n1, n2, N3, n4, n5))

for q in range(n0):
    for p in range(n1):
        for k in range(n2):
            for l in range(N3):
                for h in range(l * num, (l + 1)* num):
                    for i in range(n4):
                        for j in range(n5):
                            b[q, p, k, l, i, j] += a[q, p, k, h, i, j]

print(a.shape)
print(b.shape)
print(a[0, 0, 0, 0, :, :], "\n")
print(a[0, 0, 0, 1, :, :], "\n")
print(b[0, 0, 0, 0, :, :])



m_h = np.arange(6 * 6 * 6 * 6 * 6 * 6)
m_h = m_h.reshape(6, 6, 6, 6, 6, 6).astype(ctypes.c_float)

n_h = np.arange(6 * 6 * 6 * 6 * 6 * 6)
n_h = m_h.reshape(6, 6, 6, 6, 6, 6).astype(ctypes.c_int)

r_h = np.arange(6 * 6 * 6 * 6 * 6 * 6)
r_h = m_h.reshape(6, 6, 6, 6, 6, 6).astype(ctypes.c_float)

num = int(2)
n0, n1, n2, n3, n4, n5 = m_h.shape
N3 = int(n3 / num)

m_l = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)
n_l = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_int)
r_l = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)


maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library

float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
maputilslib.dgradeZ6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                  float32_array6, int32_array6, float32_array6,
                                  ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                  ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                  ctypes.c_int, ctypes.c_int]

maputilslib.dgradeZ6D(m_h, n_h, r_h,
                      m_l, n_l, r_l,
                      n0,  n1, n2,
                      n3,  n4, n5,
                      N3, num)

print(np.allclose(a, n_h))
print(np.allclose(b, n_l))
"""
"""
A = np.arange(6 * 6 * 6 * 6)
A = A.reshape(6, 6, 6, 6)

numXY = 2
numZ  = 2
n0, n1, n2, n3 = A.shape 
N1, N2, N3 = int(n1 / numZ), int(n2 / numXY), int(n3 / numXY)
B = np.zeros((n0, N1, N2, N3))

for q in range(n0):
    for l in range(N1):
        for a in range(N2):
            for b in range(N3):

                for h in range(l * numZ, (l + 1)* numZ):
                    for i in range(a * numXY, (a + 1) * numXY):
                        for j in range(b * numXY, (b + 1) * numXY):
                            B[q, l, a, b] += A[q, h, i, j]
print(A[0, 0, :, :], "\n")
print(A[0, 1, :, :], "\n")
print(B[0, 0, :, :])

m_h = np.arange(6 * 6 * 6 * 6)
m_h = m_h.reshape(6, 6, 6, 6).astype(ctypes.c_float)

n_h = np.arange(6 * 6 * 6 * 6)
n_h = m_h.reshape(6, 6, 6, 6).astype(ctypes.c_int)

r_h = np.arange(6 * 6 * 6 * 6)
r_h = m_h.reshape(6, 6, 6, 6).astype(ctypes.c_float)

numZ = int(2)
numXY = int(2)
n0, n1, n2, n3 = m_h.shape
N1, N2, N3 = int(n1 / numZ), int(n2 / numXY), int(n3 / numXY)

m_l = np.zeros((n0, N1, N2, N3), dtype = ctypes.c_float)
n_l = np.zeros((n0, N1, N2, N3), dtype = ctypes.c_int)
r_l = np.zeros((n0, N1, N2, N3), dtype = ctypes.c_float)


maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library

float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
maputilslib.dgradeXYZ4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                    float32_array4, int32_array4, float32_array4,
                                    ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                    ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                    ctypes.c_int,   ctypes.c_int, ctypes.c_int]

maputilslib.dgradeXYZ4D(m_h, n_h,  r_h,
                        m_l, n_l,  r_l,
                        n0,  n1,   n2,
                        n3,  N1,   N2,  
                        N3,  numZ, numXY)

print(np.allclose(A, n_h))
print(np.allclose(B, n_l))
"""
"""
A = np.arange(3 * 3 * 3 * 3 * 3 * 3)
A = A.reshape(3, 3, 3, 3, 3, 3)

numXY = 2
n0, n1, n2, n3, n4, n5 = A.shape 
N4, N5 = n4 * numXY, n5 * numXY
B = np.zeros((n0, n1, n2, n3, N4, N5))

for q in range(n0):
    for l in range(n1):
        for k in range(n2):
            for p in range(n3):
                for a in range(n4):
                    for b in range(n5):
                        for i in range(a * numXY, (a + 1) * numXY):
                            for j in range(b * numXY, (b + 1) * numXY):
                                B[q, l, k, p, i, j] = A[q, l, k, p, a, b]
print(A[0, 0, 0, 0, :, :], "\n")
print(A[0, 0, 0, 1, :, :], "\n")
print(B[0, 0, 0, 0, :, :])


m_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
m_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_float)

n_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
n_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_int)

r_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
r_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_float)

numZ = int(2)
numXY = int(2)
n0, n1, n2, n3, n4, n5 = m_l.shape
N4, N5 = n4 * numXY, n5 * numXY

m_h = np.zeros((n0, n1, n2, n3, N4, N5), dtype = ctypes.c_float)
n_h = np.zeros((n0, n1, n2, n3, N4, N5), dtype = ctypes.c_int)
r_h = np.zeros((n0, n1, n2, n3, N4, N5), dtype = ctypes.c_float)


maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library

float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
maputilslib.ugradeXY6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                   float32_array6, int32_array6, float32_array6,
                                   ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                   ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                   ctypes.c_int,   ctypes.c_int, ctypes.c_int]
maputilslib.ugradeXY6D(m_h, n_h,  r_h,
                       m_l, n_l,  r_l,
                       n0,  n1,   n2,
                       n3,  n4,   n5, 
                       N4,  N5,   numXY)
print(m_l.shape, m_h.shape)

print(np.allclose(A, n_l))
print(np.allclose(B, n_h))
"""

A = np.arange(3 * 3 * 3 * 3 * 3 * 3)
A = A.reshape(3, 3, 3, 3, 3, 3)

numZ = 4
n0, n1, n2, n3, n4, n5 = A.shape 
N3 = n3 * numZ
B = np.zeros((n0, n1, n2, N3, n4, n5))

for q in range(n0):
    for h in range(n1):
        for p in range(n2):
            for k in range(n3):
                for i in range(n4):
                    for j in range(n5):
                        for z in range(k * numZ, (k + 1) * numZ):
                            B[q, h, p, z, i, j] = A[q, h, p, k, i, j]

print(A[0, 0, 0, 0, :, :], "\n")
print(A[0, 0, 0, 1, :, :], "\n")
print(A[0, 0, 0, 2, :, :], "\n")
print(B[0, 0, 0, 0, :, :], "\n")
print(B[0, 0, 0, 1, :, :], "\n")
print(B[0, 0, 0, 2, :, :], "\n")
print(B[0, 0, 0, 3, :, :], "\n")
print(B[0, 0, 0, 4, :, :], "\n")
print(B[0, 0, 0, 5, :, :])


m_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
m_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_float)

n_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
n_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_int)

r_l = np.arange(3 * 3 * 3 * 3 * 3 * 3)
r_l = m_l.reshape(3, 3, 3, 3, 3, 3).astype(ctypes.c_float)

numZ = int(4)
numXY = int(2)
n0, n1, n2, n3, n4, n5 = m_l.shape
N3 = n3 * numZ

m_h = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)
n_h = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_int)
r_h = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)


maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library

float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
maputilslib.ugradeZ6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                   float32_array6, int32_array6, float32_array6,
                                   ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                   ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                   ctypes.c_int,   ctypes.c_int]
print(m_l.shape, m_h.shape)

maputilslib.ugradeZ6D(m_h, n_h,  r_h,
                        m_l, n_l,  r_l,
                        n0,  n1,   n2,
                        n3,  n4,   n4,
                        N3,  numZ)

print(np.allclose(A, n_l))
print(np.allclose(B, n_h))
print("second")
print(m_l[0, 0, 0, 0, :, :], "\n")
print(m_l[0, 0, 0, 1, :, :], "\n")
print(m_l[0, 0, 0, 2, :, :], "\n")
print(m_h[0, 0, 0, 0, :, :], "\n")
print(m_h[0, 0, 0, 1, :, :], "\n")
print(m_h[0, 0, 0, 2, :, :], "\n")
print(m_h[0, 0, 0, 3, :, :], "\n")
print(m_h[0, 0, 0, 4, :, :], "\n")
print(m_h[0, 0, 0, 5, :, :], "\n")
print(m_h[0, 0, 0, 6, :, :], "\n")
print(m_h[0, 0, 0, 7, :, :], "\n")
print(m_h[0, 0, 0, 8, :, :], "\n")
