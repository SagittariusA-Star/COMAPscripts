# import ctypes
import ctypes
import numpy as np

maplib = ctypes.cdll.LoadLibrary("maplib.so.1")  # Load shared library

float64_array2 = np.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, flags="contiguous")

maplib.makemaps.argtypes = [float64_array2, ctypes.c_int, ctypes.c_int]

nx = 5
ny = 7
map = np.zeros((nx, ny))
print(map)
maplib.makemaps(map, nx, ny)
print(map)
