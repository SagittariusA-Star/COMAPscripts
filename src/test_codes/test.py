import numpy as np
import numpy.ctypeslib as npct
import ctypes

liblr = ctypes.CDLL('test.so')
"""
array_2d_double = npct.ndpointer(dtype=np.uint64,ndim=2, flags='CONTIGUOUS')
liblr = npct.load_library('test.so', '.')

liblr.test.restype = None
liblr.test.argtypes = [array_2d_double, ctypes.c_int]

x = np.arange(100).reshape((10,10)).astype(np.float64)
liblr.test(x, 10)

"""
c_float_p = ctypes.POINTER(ctypes.c_float)
a = np.random.normal(loc = 0, scale = 1, size = (200, 100))
a = a.astype(np.float32)
a_p = a.ctypes.data_as(c_float_p)
print(a_p)
liblr.test(a_p, 200, 100)
