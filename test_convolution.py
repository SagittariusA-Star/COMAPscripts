#import numpy.fft as fft
import numpy as np
import h5py as h5
import ctypes
import scipy
from scipy import signal
import matplotlib.pyplot as plt

np.random.seed(123345)

a = np.random.randint(0, 9, (10, 10, 10, 201, 201))
a[:, :, :, 50, :] = np.linspace(0, 50, 201)
a[:, :, :, :, 50] = np.linspace(0, -50, 201)
a[:, :, :, :, 75] = 50 * np.sin(np.linspace(0, -50, 201))

a = a.reshape(10, 10 * 10, 201, 201)

def gaussian_kernel(sigma_x, sigma_y, n_sigma=5.0):
    size_y = int(n_sigma * sigma_y)
    size_x = int(n_sigma * sigma_x)
    x, y = np.mgrid[-size_x:size_x + 1, -size_y:size_y + 1]
    g = np.exp(-(x**2 / (2. * sigma_x**2) + y**2 / (2. * sigma_y**2)))
    return g / g.sum()


def gaussian_smooth(mymap, sigma_x, sigma_y, n_sigma=5.0):
    kernel = gaussian_kernel(sigma_x, sigma_y, n_sigma=n_sigma)
    smoothed_map = signal.fftconvolve(mymap, kernel[np.newaxis, np.newaxis, :, :], mode='same', axes = [len(mymap.shape) - 2, len(mymap.shape) - 1])
    #[len(mymap.shape)-2, len(mymap.shape) - 1])
    return smoothed_map

def gaussian_kernelZ(sigma_z, n_sigma):
        size_z = int(n_sigma * sigma_z)
        z = np.arange(-size_z, size_z + 1)
        g = np.exp(-(z**2 / (2. * sigma_z ** 2)))
        return g / g.sum()

def gaussian_smoothZ(mymap, sigma_z, n_sigma=5.0):
    kernel = gaussian_kernelZ(sigma_z, n_sigma=n_sigma)
    smoothed_map = signal.fftconvolve(mymap, kernel[np.newaxis, :, np.newaxis, np.newaxis], mode='same', axes = len(mymap.shape) - 3)
    #[len(mymap.shape)-2, len(mymap.shape) - 1])
    return smoothed_map

def gaussian_kernelXYZ(sigma_x, sigma_y, sigma_z, n_sigma=5.0):
    size_y = int(n_sigma * sigma_y)
    size_x = int(n_sigma * sigma_x)
    size_z = int(n_sigma * sigma_z)

    z, x, y = np.mgrid[-size_z:size_z + 1, -size_x:size_x + 1, -size_y:size_y + 1]
    g = np.exp(-(x**2 / (2. * sigma_x**2) + y**2 / (2. * sigma_y**2) + z**2 / (2. * sigma_z**2)))
    print(size_x, size_y, size_z)
    print(z[:, 0, 0])
    print(x[0, :, 0])
    print(y[0, 0, :])
    return g / g.sum()


def gaussian_smoothXYZ(mymap, sigma_x, sigma_y, sigma_z, n_sigma = 5.0):
    kernel = gaussian_kernelXYZ(sigma_x, sigma_y, sigma_z, n_sigma = n_sigma)
    smoothed_map = signal.fftconvolve(mymap, kernel[np.newaxis, :, :, :], mode='same', axes = [len(mymap.shape) - 3, len(mymap.shape) - 2, len(mymap.shape) - 1])
    #[len(mymap.shape)-2, len(mymap.shape) - 1])
    return smoothed_map



b = gaussian_smoothXYZ(a, 5, 5, 1)

fig, ax = plt.subplots(1, 2)
ax[0].imshow(a[5, :, :, 0].T)
ax[0].set_title("Test map")
ax[1].imshow(b[5, :, :, 0].T)
ax[1].set_title("Smoothed")
plt.savefig("test_smoothed_map.png")
plt.show()