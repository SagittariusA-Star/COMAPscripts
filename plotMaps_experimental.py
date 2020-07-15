import sys
import getopt
import numpy as np
import h5py as h5
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib import animation, rc
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import textwrap
#import FileDialog
import time

def readMap(filename):
    dfile   = h5.File(filename,'r')
    x       = dfile['x'][:]#;         x  = np.array(x[:]).astype(float)
    y       = dfile['y'][:]#;         y  = np.array(y[:]).astype(float)
    freq    = dfile['freq'][:]#;    freq = np.array(freq[...]).astype(float)
   
    maps    = dfile['map'][:]#;       maps  = np.array(maps[...]).astype(float)
    hit     = dfile['nhit'][:]#;      hit  = np.array(hit[...]).astype(float)
    rms     = dfile['rms'][:]#;       rms  = np.array(rms[...]).astype(float)
    return x, y, maps, hit, rms, freq



def writeMap(outfile, x, y, maps, hit, rms, freq):
    h5_file = h5.File(outfile, "w")
    h5_file.create_dataset("x", data = x)
    h5_file.create_dataset("y", data = y)
    h5_file.create_dataset("freq", data = freq)
    h5_file.create_dataset("map", data = maps)
    h5_file.create_dataset("hit", data = hit)
    h5_file.create_dataset("rms", data = rms)
    b = h5_file["map"][:]
    h5_file.close()


start = time.time()
x, y, maps, hit, rms, freq = readMap("co2_bigmap.h5")
print("Run time: ", time.time() - start, " sec")
writeMap("out_map.hdf5", x, y, maps, hit, rms, freq)