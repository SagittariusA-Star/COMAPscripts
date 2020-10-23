import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(ncols=2)

plt.axes(ax1)
hp.mollview(np.random.random(hp.nside2npix(32)), hold=True)

plt.axes(ax2)
hp.mollview(np.arange(hp.nside2npix(32)), hold=True)
plt.show()