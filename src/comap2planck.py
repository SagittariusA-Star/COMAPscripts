import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import healpy as hp
import h5py as h5
import sys
import getopt
import textwrap
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy_healpix import HEALPix
from astropy.io import fits
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable

def get_ring(center, radius, param):
    phi = np.linspace(0, 2 * np.pi, int(1e4))
    x = radius * np.cos(phi) + center[0]
    y = param * radius * np.sin(phi) + center[1]

    return x, y

fontsize = 12
fonts = {"axes.labelsize": fontsize,
    "font.size": fontsize,
    "legend.fontsize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize
    }
plt.rcParams.update(fonts)


cmb_name = "/home/sagittarius/Documents/Summerjob/COMAPscripts/LFI_SkyMap_30GHz.fits"

healmap = hp.fitsfunc.read_map(cmb_name)

healmap *= 1e6

fig, ax0 = plt.subplots(figsize = (16, 8))


#cmap = copy.copy(plt.cm.get_cmap("RdYlBu"))
cmap = copy.copy(plt.cm.get_cmap("magma"))
#cmap = copy.copy(plt.cm.get_cmap("tab10"))
#cmap = copy.copy(plt.cm.get_cmap("viridis"))
#cmap = copy.copy(plt.cm.get_cmap("rainbow"))
#cmap = copy.copy(plt.cm.get_cmap("Spectral"))
#cmap = copy.copy(plt.cm.get_cmap("CMRmap"))
#cmap = cmap.reversed()

lon_co2, lat_co2 = get_ring((149.0, -60.3), 2, 1)
lon_co6, lat_co6 = get_ring((91.35, 53.22), 2, 1)
lon_co7, lat_co7 = get_ring((150.64, 59.53), 2, 1)
#lon, lat = get_ring((0, 0), 30)


plt.axes(ax0)
skymap = hp.mollview(
    healmap,
    cmap = cmap,
    nest=False,
    title=None,
    xsize = 5000,
    cbar = True,
    return_projected_map = True,
    hold = True,
    unit = r"$\mu K_\mathrm{CMB}$",
    min = 0, 
    max = 1200,
    margins = (0, 1, 0, 0.5)
    )
hp.projplot(lon_co2, lat_co2, c = "w", lw = 1, lonlat = True)
hp.projplot(lon_co6, lat_co6, c = "w", lw = 1, lonlat = True)
hp.projplot(lon_co7, lat_co7, c = "w", lw = 1, lonlat = True)

hp.projtext(149.0 - 5, -60.3, "CO2", c = "w", lonlat = True)
hp.projtext(91.35 - 5, 53.22, "CO6", c = "w", lonlat = True)
hp.projtext(150.64 - 8, 59.53, "CO7", c = "w", lonlat = True)

f = plt.gcf().get_children()
CbAx = f[2]

unit_text_obj = CbAx.get_children()[1]
unit_text_obj.set_fontsize(fontsize)

#hp.graticule(local = False, color = "gray")

plt.savefig("/home/sagittarius/Documents/COMAP_general/COMAP_general/src/sim/fields_in_healpix_magma.pdf", dpi = 150)
#plt.show()

#fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (16, 8))
dummyfig, dummyax = plt.subplots()

fig1 = plt.figure(figsize = (16, 8))

gs = mpl.gridspec.GridSpec(nrows=1, ncols=3, wspace = 0.01)
ax1 = fig1.add_subplot(gs[0])
ax2 = fig1.add_subplot(gs[1])
ax3 = fig1.add_subplot(gs[2])


plt.axes(ax1)

hp.delgraticules()
co2 = hp.cartview(
    healmap,
    cmap = cmap,
    rot = (149.0, -60.3),
    lonra = [- 5, 5],
    latra = [- 5, 5],
    nest=False,
    title=None,
    xsize = 5000,
    min = 0, 
    max = 500,
    cbar = False,
    return_projected_map = True,
    hold = True,
    unit = r"$\mu K_\mathrm{CMB}$",
    notext = True,
    format = "%f",
    aspect = "equal"
)
hp.projplot(lon_co2, lat_co2, c = "w", lw = 1, lonlat = True)
hp.projtext(149.0 - 3, -60.3, "CO2", c = "w", lonlat = True)

hp.graticule(color = "gray", dmer = 0.01, 
    dpar = 0.01)
    
plt.axes(ax2)
hp.delgraticules()

co6 = hp.cartview(
    healmap,
    cmap = cmap,
    rot = (91.35, 53.22),
    lonra = [- 5, 5],
    latra = [- 5, 5],
    nest=False,
    title=None,
    xsize = 5000,
    min = 0, 
    max = 500,
    cbar = False,
    return_projected_map = True,
    hold = True,
    unit = r"$\mu K_\mathrm{CMB}$", 
    notext = True,
    format = "%f",
    aspect = "equal"
)
hp.projplot(lon_co6, lat_co6, c = "w", lw = 1, lonlat = True)
hp.projtext(91.35 - 3, 53.22, "CO6", c = "w", lonlat = True)

hp.graticule(color = "gray", dmer = 0.01, 
    dpar = 0.01)

plt.axes(ax3)
hp.delgraticules()

co7 = hp.cartview(
    healmap,
    cmap = cmap,
    rot = (150.64, 59.53),
    lonra = [- 5, 5],
    latra = [- 5, 5],
    nest=False,
    title=None,
    xsize = 5000,
    min = 0, 
    max = 500,
    cbar = False,
    return_projected_map = True,
    hold = True,
    unit = r"$\mu K_\mathrm{CMB}$", 
    notext = True,
    format = "%.f",
    aspect = "equal"
)

hp.projplot(lon_co7, lat_co7, c = "w", lw = 1, lonlat = True)
hp.projtext(150.64 - 3, 59.53, "CO7", c = "w", lonlat = True)

hp.graticule(color = "gray", dmer = 0.01, 
    dpar = 0.01)



dummyim = dummyax.imshow(co7, vmin = 0, vmax = 500, cmap = cmap)

cb_ax = fig1.add_axes([0.91, 0.24, 0.02, 0.512])

cbar = fig1.colorbar(dummyim, cax = cb_ax, ax = ax3, ticks = [0, 125, 250, 375, 500])
cbar.set_label(r'$\mu K_\mathrm{CMB}$')

#cbar.ax.set_xticklabels(["0", "1"], rotation = 90)
print(cbar)

plt.savefig("/home/sagittarius/Documents/COMAP_general/COMAP_general/src/sim/cutout_in_magma.pdf", dpi = 150)




#im0 = ax0.imshow(skymap[::-1, :], vmin = 0, vmax  = 1200, cmap = cmap)
#im1 = ax1.imshow(co2[::-1, :], vmin = 0, vmax  = 1200, cmap = cmap)
#im2 = ax2.imshow(co6[::-1, :], vmin = 0, vmax  = 1200, cmap = cmap)
#im3 = ax3.imshow(co7[::-1, :], vmin = 0, vmax  = 1200, cmap = cmap)

#divider = make_axes_locatable(ax0)
#cax = divider.append_axes("bottom", size = "5%", pad = 0.05)

#cbar = fig.colorbar(im0, ax = ax0, cax = cax, orientation = "horizontal")
#cbar.set_label(r'$\mu K_\mathrm{CMB}$')
#print(skymap.shape)

"""
fig = plt.gcf()
ax = plt.gca()
image = ax.get_images()[0]
cbar = fig.colorbar(image, ax=ax, orientation = "horizontal")
cbar.ax.tick_params(labelsize=14)
cbar.ax.set_ylabel("Hei")
"""


#cbar.ax.set_ticks([-1, 0, 1]) 

#cbar.ax.set_ticklabels(['negative', 'zero', 'positive']) 

"""
hp.projplot(lon_co2, lat_co2, c = "w", lw = 1, lonlat = True)
hp.projplot(lon_co6, lat_co6, c = "w", lw = 1, lonlat = True)
hp.projplot(lon_co7, lat_co7, c = "w", lw = 1, lonlat = True)

hp.projtext(149.0 - 5, -60.3, "CO2", c = "w", lonlat = True)
hp.projtext(91.35 - 5, 53.22, "CO6", c = "w", lonlat = True)
hp.projtext(150.64 - 8, 59.53, "CO7", c = "w", lonlat = True)
"""

#plt.axes(ax2)



"""
f = plt.gcf().get_children()
HpxAx = f[1]
CbAx = f[2]
print(len(f))
print(f[0])
print(f[1])
print(f[2])

coord_text_obj = HpxAx.get_children()[0]
coord_text_obj.set_fontsize(fontsize)
print(coord_text_obj)

unit_text_obj = CbAx.get_children()[1]
unit_text_obj.set_fontsize(fontsize)
print(CbAx.get_children())

plt.rcParams.update({'font.size':fontsize})
"""

plt.show()

"""

import healpy as hp
import numpy as np
import matplotlib

fontsize = 20

d = np.arange(12*16**2)
hp.mollview(d, title='Hello', unit=r'T', notext=False, coord=['G','C'])

matplotlib.rcParams.update({'font.size':fontsize})
matplotlib.pyplot.show()

f = matplotlib.pyplot.gcf().get_children()
HpxAx = f[1]
CbAx = f[2]

coord_text_obj = HpxAx.get_children()[0]
coord_text_obj.set_fontsize(fontsize)

unit_text_obj = CbAx.get_children()[1]
unit_text_obj.set_fontsize(fontsize)

matplotlib.pyplot.show()
"""