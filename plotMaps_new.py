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
import FileDialog


def usage():
    prefix = ""
    preferredWidth = 150
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=preferredWidth,subsequent_indent=' '*len(prefix))
    m1 = "(Side band as a list, ie. [1,2,4]. Default all)"
    m2 = "(Frequency, given as list, ie. [1,6,26]. Default all) "
    m3 = "(Type of plot. Choices are map, rms, map/rms, sim, rms_sim, hit, feed, and var. Default map) "
    m4 = "(Filename)"
    m5 = "(color_limits of colorbar) as nested list. First plot --> first list in main list. Default none)"
    m6 = "(Which detector as a list, ie. [4,11,18]. Default all)"
    m7 = "(Outfile name, default 'outfile')"
    m8 = "(x range [xmin,xmax])"
    m9 = "(y range [ymin,ymax])"
    m10 = "(Include to split the map as the focal plane)"
    m11 = "(If the field non-stationary)"
    m12 = "(Which simulation, default first element (0))"
    m13 = "(Scale the data by a factor w, default 1)"
    m14 = "(rms limit. Include a mask that removes pixels with the given rms limit.)"
    m15 = "(x vs z for a given y index)"
    m16 = "(y vs z for a given x index)"

    print("\nThis is the usage function\n")
    print("Flags:")
    print("-f ----> optional --filename " + wrapper.fill(m4))
    print("-o ----> optional --out "      + wrapper.fill(m7))
    print("-p ----> optional --plots "    + wrapper.fill(m3))
    print("-d ----> optional --det "      + wrapper.fill(m6))
    print("-s ----> optional --sb "       + wrapper.fill(m1))
    print("-n ----> optional --nu "       + wrapper.fill(m2))
    print("-c ----> optional --colorlim " + wrapper.fill(m5))
    print("-w ----> optional --scale "    + wrapper.fill(m13))
    print("-m ----> optional --mask "     + wrapper.fill(m14))
    print("-x ----> optional --xlim "     + wrapper.fill(m8))
    print("-y ----> optional --ylim "     + wrapper.fill(m9))
    print("-j ----> optional --jupiter "  + wrapper.fill(m11))
    print("-z ----> optional --sim "      + wrapper.fill(m12))
    print("-r ----> optional --deepx"     + wrapper.fill(m15))
    print("-l ----> optional --deepy"     + wrapper.fill(m16))
    sys.exit()

plot_choices = ["map", "rms", "map/rms", "hit", "sim", "rms_sim", "feed", "var"]
plot         = "map"
freq         = "all"
set_xlim     = False
set_ylim     = False
x_lim        = [0,360]
y_lim        = [-90,90]
det_list     = range(1,20)
sb_list      = range(1,5)
freq_list    = range(1,65)
plots        = ["map"]
color_lim    = [None, None]
set_colorlim = False
outfile      = "outfile"
jupiter      = False
sim_numb     = 0
scale        = 1
beam         = False
patch        = ''
rms_lim      = 200000.
deepx         = False
deepy         = False

if len(sys.argv) == 1:
    usage()

try:
    opts, args = getopt.getopt(sys.argv[1:],"s:n:p:f:h:c:d:o:x:y:j:z:w:m:r:l:", ["sb=", "nu=", "plots=", "filename=", "help=", "colorlim", "det", "out", "xlim", "ylim","jupiter", "sim","scale","mask","deepx", "deepy"])
except getopt.GetoptError:
    usage()

for opt, arg in opts:
    if opt in ("-x", "--xlim"):
        set_xlim = True
        x_lim = eval(arg)
        if type(x_lim) != list:
            print("xlim needs to be a list with [min, max]")
            sys.exit()
    elif opt in ("-y", "--ylim"):
        set_ylim = True
        y_lim = eval(arg)
        if type(y_lim) != list:
            print("ylim needs to be a list with [min, max]")
            sys.exit()
    elif opt in ("-o", "--out"):
        outfile = arg
    elif opt in ("-d", "--det"):
        det_list = eval(arg)
        if type(det_list) != list:
            print("Detectors my be inserted as a list, ie. -d [1,2,5,7]")
            sys.exit()
        else:
            if 0 in det_list:
                print("Use 1-base, not 0-base please")
                sys.exit()
    elif opt in ("-s", "--sb"):
        sb_list = eval(arg)
        if type(sb_list) != list:
            print("Side bands my be inserted as a list, ie. -d [1,2,4]")
            sys.exit()
        else:
            if 0 in sb_list:
                print("Use 1-base, not 0-base please")
                sys.exit()
    elif opt in ("-n", "--nu"):
        freq_list = eval(arg)
        if type(freq_list)!= list:
            print("Frequencies my be inserted as a list, ie. -n [1,34,50]")
        else:
            if 0 in freq_list:
                print("Use 1-base, not 0-base please")
                sys.exit()
            for i in freq_list:
                if i > 64:
                    print("There are only 64 frequencies pr. side band")
                    sys.exit()
    elif opt in ("-p", "--plots"):
        plot = arg
        if plot not in plot_choices:
            print("Make sure you have chosen the correct plot choices")                                                                                                   
            sys.exit() 
        if plot == "map/rms":
            plot = "plot_rms"
        #if plot in plot_choices:
        #    plot_list = eval(arg)
        #else:
        #    option = ""
        #    plot_list = []
        #    for i in arg:
        #        if i == "[":
        #            continue
        #        elif i == "]":
        #            plot_list.append(option)
        #        elif i == ",":
        #            plot_list.append(option)
        #            option = ""
        #        else:
        #            option += i
        #for i, op in enumerate(plot_list):
        #    if op in plot_choices:
        #        if op == "map/rms":
        #            plot_list[i] = "map_rms"
        #    else:
        #        print "Make sure you have chosen the correct plot choices"
        #        sys.exit()
    elif opt in ("-f", "--filename"):
        filename = arg
        temp = filename.split('/')[-1]
        patch = temp.split('_')[0]
    elif opt in ("-h", "--help"):
        usage()
    elif opt in ("-c", "--colorlim"):
        set_colorlim = True
        color_lim = eval(arg)
    elif opt in ("-j", "--jupiter"):
        jupiter = True
    elif opt in ("-z", "--sim"):
        sim_numb = int(arg)
    elif opt in ("-w", "--scale"):
        scale = eval(arg)
    elif opt in ("-m", "--mask"):
        rms_lim = eval(arg)
    elif opt in ("-r", "--deepx"):
        deepx = True
        y_index = int(arg)-1
    elif opt in ("-l", "--deepy"):
        deepy = True
        x_index = int(arg)-1
    else:
        usage()

#if set_colorlim == False:
#    color_lim = [None, None]
    #for i in range(len(plot_list)):
    #    color_lim.append([None, None])

def readMap(filename):
    dfile   = h5.File(filename,'r')
    sim_exist = "map_sim" in dfile

    nx      = dfile['n_x'];      nx  = np.array(nx).astype(int)
    ny      = dfile['n_y'];      ny  = np.array(ny).astype(int)
    x       = dfile['x'];         x  = np.array(x[:]).astype(float)
    y       = dfile['y'];         y  = np.array(y[:]).astype(float)
    freq    = dfile['freq'];    freq = np.array(freq[...]).astype(float)
    if len(det_list) == 19 and plot != 'feed':
        maps    = dfile['map_beam'];  maps  = np.array(maps[...]).astype(float)
        hit     = dfile['nhit_beam']; hit  = np.array(hit[...]).astype(float)
        rms     = dfile['rms_beam'];  rms  = np.array(rms[...]).astype(float)
        beam = True
    else:
        maps    = dfile['map'];       maps  = np.array(maps[...]).astype(float)
        hit     = dfile['nhit'];      hit  = np.array(hit[...]).astype(float)
        rms     = dfile['rms'];       rms  = np.array(rms[...]).astype(float)
    if sim_exist:
        maps_sim = dfile['map_sim'];  maps_sim = np.array(maps_sim[...]).astype(float)
        rms_sim  = dfile['rms_sim'];  rms_sim  = np.array(rms_sim[...]).astype(float)
    else:
        maps_sim = []
        rms_sim = []
    
    return x, y, maps, hit, rms, maps_sim, rms_sim, freq

def setup(filename, outfile, x_lim, y_lim):
    print(filename.split("/")[-1])
    x, y, maps, hit, rms, maps_sim, rms_sim, freq = readMap(filename)
    nx = len(x)
    ny = len(y)
    #ndet = maps.shape[0]
    nsb = freq.shape[0]
    nfreq = freq.shape[1]
    freq_long = np.zeros(nsb*nfreq)

    for i in range(nsb):
        for j in range(nfreq):
            freq_long[i*nfreq + j] = freq[i,j]
    
    #if len(plot_list) == 1:                                                                                                             
    fig, ax = plt.subplots(1)
    fig.set_figheight(5)
    fig.set_figwidth(9)
    #axarr = [ax]
    #else:                                                                                                                               
    #    fig, axarr = plt.subplots(len(plot_list))                                                                                       

    #for p, plot in enumerate(plot_list):                                                                                                
    mapData    = np.zeros((ny,nx))
    mapZ       = np.zeros((ny,nx,nsb*nfreq))
    rmsZ       = np.zeros((ny,nx,nsb*nfreq))
    hitZ       = np.zeros((ny,nx,nsb*nfreq))
    map2Data   = np.zeros((ny,nx))
    hitData    = np.zeros((ny,nx))
    rmsData    = np.zeros((ny,nx))
    maprmsData = np.zeros((ny,nx))
    varData    = np.zeros((ny,nx))
    seenbyfeed = np.zeros((ny,nx))
    hit_temp   = np.zeros((19,ny,nx))
    simData    = np.zeros((ny,nx))
    simrmsData = np.zeros((ny,nx))
    for p in range(nx):
        for q in range(ny):
            data_temp = []
            if len(det_list) == 19 and plot != 'feed':
                if deepx or deepy:
                    for j in range(nsb):
                        for k in range(nfreq):
                            if hit[j-1,k-1,q,p] == 0: 
                                continue
                            mapZ[q,p,j*nfreq+k] += maps[j-1,k-1,q,p]/rms[j-1,k-1,q,p]**2
                            rmsZ[q,p,j*nfreq+k] += 1./rms[j-1,k-1,q,p]**2
                            hitZ[q,p,j*nfreq+k] += hit[j-1,k-1,q,p]
                for j in sb_list:
                    for k in freq_list:
                        if hit[j-1,k-1,q,p] == 0: 
                            continue
                        hitData[q,p] += hit[j-1,k-1,q,p]
                        rmsData[q,p] += 1./rms[j-1,k-1,q,p]**2
                        mapData[q,p] += maps[j-1,k-1,q,p]/rms[j-1,k-1,q,p]**2
                        #map2Data[q,p] += maps[j-1,k-1,q,p]
                        data_temp.append(maps[j-1,k-1,q,p] / rms[j-1,k-1,q,p])
                        if plot == 'sim' or plot == 'rms_sim':
                            if len(maps_sim) == 0:
                                print("Could not find any simulated data")
                                sys.exit()
                            simData[q,p] += maps_sim[sim_numb,j-1,k-1,q,p]/rms_sim[sim_numb,j-1,k-1,q,p]**2
                            simrmsData[q,p] += 1./rms_sim[sim_numb,j-1,k-1,q,p]**2
            else:
                if deepx or deepy:
                    for j in range(nsb):
                        for k in range(nfreq):
                            if hit[i-1,j-1,k-1,q,p] == 0: 
                                continue
                            mapZ[q,p,j*nfreq+k] += maps[i-1,j-1,k-1,q,p]/rms[i-1,j-1,k-1,q,p]**2
                            rmsZ[q,p,j*nfreq+k] += 1./rms[i-1,j-1,k-1,q,p]**2
                            hitZ[q,p,j*nfreq+k] += hit[i-1,j-1,k-1,q,p]
                for j in sb_list:
                    for k in freq_list:
                        map_temp = 0.0; rms_temp = 0.0
                        for i in det_list:
                            if hit[i-1,j-1,k-1,q,p] == 0:
                                continue
                            hit_temp[i-1,:,:] += hit[i-1,j-1,k-1,q,p]
                            hitData[q,p] += hit[i-1,j-1,k-1,q,p]
                            rmsData[q,p] += 1./rms[i-1,j-1,k-1,q,p]**2
                            mapData[q,p] += maps[i-1,j-1,k-1,q,p]/rms[i-1,j-1,k-1,q,p]**2
                            map2Data[q,p] += maps[i-1,j-1,k-1,q,p]
                            map_temp += maps[i-1,j-1,k-1,q,p]
                            rms_temp += 1./rms[i-1,j-1,k-1,q,p]**2
                            if plot == 'sim' or plot == 'rms_sim':
                                if len(maps_sim) == 0:
                                    print("Could not find any simulated data")
                                    sys.exit()
                                simData[q,p] += maps_sim[sim_numb,i-1,j-1,k-1,q,p]/rms_sim[sim_numb,i-1,j-1,k-1,q,p]**2
                                simrmsData[q,p] += 1./rms_sim[sim_numb,i-1,j-1,k-1,q,p]**2
                        data_temp.append(map_temp * np.sqrt(rms_temp))
            varData[q,p] = np.var(data_temp) #np.sqrt(np.var(data_temp)) 
    
    if plot == 'feed':
        for i in det_list:
            for m in range(nx):
                for n in range(ny):
                    if hit_temp[i-1,m,n] > 0.01*hitData[m,n]:
                        seenbyfeed[m,n] = seenbyfeed[m,n] + 1
    
    mapData = mapData/rmsData
    #maprmsData = map2Data/rmsData 
    rmsData = np.sqrt(1./rmsData)
    maprmsData = mapData/rmsData
    mapZ = mapZ/rmsZ
    rmsZ = np.sqrt(1./rmsZ)
    
    if plot == 'sim' or plot == 'rms_sim':
        simData = simData/simrmsData
        simrmsData = np.sqrt(1./simrmsData)

    for p in range(nx):
        for q in range(ny):
            if rmsData[q,p] > rms_lim:
                hitData[q,p] = 0

    mapname = patch + ' ' + plot
    if len(det_list) == 1:
        mapname += ' det ' + str(det_list[0])
    if len(sb_list) == 1:
        mapname += ' sb ' + str(sb_list[0])
    if len(freq_list) == 1:
        mapname += ' freq ' + str(freq_list[0])
    if deepx:
        mapname += ' dec ' + str(y[y_index])
    if deepy:
        mapname += ' RA ' + str(x[x_index])


    data = np.zeros((ny,nx))
    
    if plot   == "hit":
        data  = hitData
    elif plot == "rms":
        data  = rmsData
    elif plot == "map_rms":
        data = maprmsData
    elif plot == "var":
        data = varData
    elif plot == "feed":
        data = seenbyfeed
    elif plot == "sim":
        data = simData
    elif plot == "rms_sim":
        data = simrmsdata
    elif deepx:
        data = mapZ[y_index,:,:]; hitData = hitZ[y_index,:,:]
        y = x; x = freq_long
    elif deepy:
        data = mapZ[:,x_index,:]; hitData = hitZ[:,x_index,:]
        y = y; x = freq_long
    else:
        data  = mapData
    
    makeMap(x, y, data, hitData, mapname, ax, x_lim, y_lim)
    
    outfile += '.png'
    plt.savefig(outfile)
    #plt.show()
    plt.close(fig)

    
def makeMap(x, y, data, hits, mapname, ax, x_lim, y_lim):
    
    data = np.ma.masked_where(hits < 1., data)
    if set_xlim == False:
        dx = x[1] - x[0]
        x_lim[0] = x[0] - 0.5*dx; x_lim[1] = x[-1] + 0.5*dx
    if set_ylim == False:
        dy = y[1] - y[0]
        y_lim[0] = y[1] - 0.5*dy; y_lim[1] = y[-1] + 0.5*dy

    if set_colorlim == False:
        if plot == 'map':
            color_lim[1] = 0.1*np.amax(data)*1e6
            color_lim[0] = -color_lim[1]
        if jupiter:
            color_lim[1] = np.amax(data)
            color_lim[0] = 0
    
 
    if jupiter:
        cmap = cm.viridis.set_bad('k',1)
        im = ax.imshow(data, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest',cmap=cm.viridis, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        ax.set_ylabel(r'$\Delta$El [deg]')
        ax.set_xlabel(r'$\Delta$RA [deg]')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label('K')
    else:
        cmap = cm.CMRmap.set_bad('0.8',1) #color of masked elements 
        if plot == 'map' or plot == 'sim':
            im = ax.imshow(data*1e6, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.CMRmap, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        else:
            im = ax.imshow(data*scale, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.CMRmap, origin='lower', vmin=color_lim[0], vmax=color_lim[1])

        if deepx:
            ax.set_ylabel('Right Ascension [deg]')
            ax.set_xlabel('Frequency [GHz]')
        elif deepy:
            ax.set_ylabel('Declination [deg]')
            ax.set_xlabel('Frequency [GHz]')
        else:
            ax.set_ylabel('Declination [deg]')
            ax.set_xlabel('Right Ascension [deg]')
        ax.set_title(mapname)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        if plot == 'hit':
            cbar.set_label('hits')
        elif plot == 'map' or plot == 'sim':
            cbar.set_label('uK')
        elif plot == 'feed':
            cbar.set_label('feeds')
        else:
            if (scale != 1):
                cbar.set_label(str(scale)+'K')
            else:
                cbar.set_label('K')

setup(filename, outfile, x_lim, y_lim)
