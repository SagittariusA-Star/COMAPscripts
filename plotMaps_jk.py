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
    m3 = "(Type of jackknife. Choices are odde, dayn, half, sdbg. Default odde) "
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
    m17 = "(Type of plot. Choices are map, rms, nhit. Default map)"

    print("\nThis is the usage function\n")
    print("Flags:")
    print("-i ----> optional --in "       + wrapper.fill(m4))
    print("-o ----> optional --out "      + wrapper.fill(m7))
    print("-j ----> optional --jk "       + wrapper.fill(m3))
    print("-p ----> optional --plot"      + wrapper.fill(m17))
    print("-d ----> optional --det "      + wrapper.fill(m6))
    print("-s ----> optional --sb "       + wrapper.fill(m1))
    print("-f ----> optional --freq "     + wrapper.fill(m2))
    print("-c ----> optional --colorlim " + wrapper.fill(m5))
    #print("-w ----> optional --scale "    + wrapper.fill(m13))
    print("-m ----> optional --mask "     + wrapper.fill(m14))
    print("-x ----> optional --xlim "     + wrapper.fill(m8))
    print("-y ----> optional --ylim "     + wrapper.fill(m9))
    #print("-r ----> optional --deepx"     + wrapper.fill(m15))
    #print("-l ----> optional --deepy"     + wrapper.fill(m16))
    sys.exit()

jk_choices   = ["odde", "dayn", "half", "sdbg"]
jk           = "odde"
plot_choices = ["map", "rms", "hit"]
plot         = "map"
set_xlim     = False
set_ylim     = False
x_lim        = [0,360]
y_lim        = [-90,90]
det_list     = range(0,19)
sb           = 2
freq         = 30
plots        = ["map"]
color_lim    = [None, None]
set_colorlim = False
outfile      = "outfile"
scale        = 1
beam         = False
patch        = ''
rms_lim      = 200000.
deepx        = False
deepy        = False

if len(sys.argv) == 1:
    usage()

try:
    opts, args = getopt.getopt(sys.argv[1:],"s:f:j:p:i:h:c:d:o:x:y:w:m:r:l:", ["sb=", "freq=", "jk=", "plot=", "in=", "help=", "colorlim", "det", "out", "xlim", "ylim","scale","mask","deepx", "deepy"])
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
        sb = eval(arg)
        if sb == 0:
            print("Use 1-base, not 0-base please")
            sys.exit()
        if freq > 4:
            print("There are only 4 sibands")
    elif opt in ("-f", "--freq"):
        freq = eval(arg)
        if freq == 0:
            print("Use 1-base, not 0-base please")
            sys.exit()
        if freq > 64:
            print("There are only 64 frequencies pr.side band")
            sys.exit()
    elif opt in ("-j", "--jk"):
        jk = arg
        if jk not in jk_choices:
            print("Make sure you have chosen the correct jk choices")                                                                                                   
            sys.exit() 
    elif opt in ("-p", "--plot"):
        plot = arg
        if plot not in plot_choices:
            print("Make sure you have chosen the correct plot choices")
            sys.exit() 

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
    elif opt in ("-i", "--filename"):
        filename = arg
        temp = filename.split('/')[-1]
        patch = temp.split('_')[0]
    elif opt in ("-h", "--help"):
        usage()
    elif opt in ("-c", "--colorlim"):
        set_colorlim = True
        color_lim = eval(arg)
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

def readMap(filename, jk):
    dfile   = h5.File(filename,'r')

    nx      = dfile['n_x'];      nx  = np.array(nx).astype(int)
    ny      = dfile['n_y'];      ny  = np.array(ny).astype(int)
    x       = dfile['x'];         x  = np.array(x[:]).astype(float)
    y       = dfile['y'];         y  = np.array(y[:]).astype(float)
    freqs   = dfile['freq'];   freqs = np.array(freqs[...]).astype(float)
    mapname = 'jackknives/map_'+str(jk)
    hitname = 'jackknives/nhit_'+str(jk)
    rmsname = 'jackknives/rms_'+str(jk)
    maps    = dfile[mapname];  maps  = np.array(maps[...]).astype(float)
    hits    = dfile[hitname];  hits  = np.array(hits[...]).astype(float)
    rms     = dfile[rmsname];  rms   = np.array(rms[...]).astype(float)
    if maps.ndim == 5:
        beam = True
        
    return x, y, maps, hits, rms, freqs

def setup(filename, outfile, x_lim, y_lim, jk):
    print(filename.split("/")[-1])
    x, y, maps, hits, rms, freqs = readMap(filename, jk)
    nx = len(x)
    ny = len(y)
    #ndet = maps.shape[0]
    nsb = freqs.shape[0]
    nfreq = freqs.shape[1]
    freq_long = np.zeros(nsb*nfreq)

    for i in range(nsb):
        for j in range(nfreq):
            freq_long[i*nfreq + j] = freqs[i,j]

    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2,2)#, figsize=(8, 5))
    
    # 0 = map1; 1 = map2; 2 = sum; 3 = diff
    mapData = np.zeros((4,ny,nx))
    hitData = np.zeros((4,ny,nx))
    rmsData = np.zeros((4,ny,nx))
    
    if maps.ndim == 5:
        mapData[0,:,:] = maps[0,sb,freq,:,:]
        mapData[1,:,:] = maps[1,sb,freq,:,:]
        hitData[0,:,:] = hits[0,sb,freq,:,:]
        hitData[1,:,:] = hits[1,sb,freq,:,:]
        rmsData[0,:,:] = rms[0,sb,freq,:,:]
        rmsData[1,:,:] = rms[1,sb,freq,:,:]
        for p in range(nx):
            for q in range(ny):
                mapData[2,q,p] = mapData[0,q,p]/rmsData[0,q,p]**2 + mapData[1,q,p]/rmsData[1,q,p]**2
                rmsData[2,q,p] = 1./rmsData[0,q,p]**2 + 1./rmsData[1,q,p]**2

    else:
        if len(det_list) == 1:
            det = det_list[0]
            mapData[0,:,:] = maps[0,det,sb,freq,:,:]
            mapData[1,:,:] = maps[1,det,sb,freq,:,:]
            hitData[0,:,:] = hits[0,det,sb,freq,:,:]
            hitData[1,:,:] = hits[1,det,sb,freq,:,:]
            rmsData[0,:,:] = rms[0,det,sb,freq,:,:]
            rmsData[1,:,:] = rms[1,det,sb,freq,:,:]
            for p in range(nx):
                for q in range(ny):
                    mapData[2,q,p] = mapData[0,q,p]/rmsData[0,q,p]**2 + mapData[1,q,p]/rmsData[1,q,p]**2
                    rmsData[2,q,p] = 1./rmsData[0,q,p]**2 + 1./rmsData[1,q,p]**2

        else:
            for p in range(nx):
                for q in range(ny):
                    for det in det_list:
                        mapData[0,q,p] += maps[0,det,sb,freq,q,p]/rms[0,det,sb,freq,q,p]**2
                        mapData[1,q,p] += maps[1,det,sb,freq,q,p]/rms[1,det,sb,freq,q,p]**2
                        hitData[0,q,p] += hits[0,det,sb,freq,q,p]
                        hitData[1,q,p] += hits[1,det,sb,freq,q,p]
                        rmsData[0,q,p] += 1./rms[0,det,sb,freq,q,p]**2
                        rmsData[1,q,0] += 1./rms[1,det,sb,freq,q,p]**2
                        mapData[2,q,p] += maps[0,det,sb,freq,q,p]/rms[0,det,sb,freq,q,p]**2 + maps[1,det,sb,freq,q,p]/rms[1,det,sb,freq,q,p]**2
                        rmsData[2,q,p] += 1./rms[0,det,sb,freq,q,p]**2 + 1./rms[1,det,sb,freq,q,p]**2
                        hitData[2,q,p] += hits[0,det,sb,freq,q,p] + hits[1,det,sb,freq,q,p]
                    mapData[0,q,p] = mapData[0,q,p]/rmsData[0,q,p]
                    mapData[1,q,p] = mapData[1,q,p]/rmsData[1,q,p]
                    rmsData[0,q,p] = 1./np.sqrt(rmsData[0,q,p])
                    rmsData[1,q,p] = 1./np.sqrt(rmsData[1,q,p])

    #mapData[0,:,:] = mapData[0,:,:] / rmsData[0,:,:]
    #mapData[1,:,:] = mapData[1,:,:] / rmsData[1,:,:]

    mapData[2,:,:] = mapData[2,:,:] / rmsData[2,:,:]
    rmsData[2,:,:] = 1./np.sqrt(rmsData[2,:,:])
    hitData[2,:,:] = hitData[0,:,:] + hitData[1,:,:]
        
    mapData[3,:,:] = mapData[0,:,:] - mapData[1,:,:] # ???
    rmsData[3,:,:] = rmsData[0,:,:] - rmsData[1,:,:] # ???
    hitData[3,:,:] = hitData[0,:,:] - hitData[1,:,:]

    '''
    for p in range(nx):
        for q in range(ny):
            if rmsData[q,p] > rms_lim:
                hitData[q,p] = 0
    '''
    
    mapname = patch + ' ' + plot
    if len(det_list) == 1:
        mapname += ' det ' + str(det_list[0])
    #if len(sb_list) == 1:
    mapname += ' sb ' + str(sb) #str(sb_list[0])
    #if len(freq_list) == 1:
    mapname += ' freq ' + str(freq) #str(freq_list[0])
    if deepx:
        mapname += ' dec ' + str(y[y_index])
    if deepy:
        mapname += ' RA ' + str(x[x_index])

    '''
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
    '''
    
    if plot == 'map':
        makeMap(x, y, mapData, hitData, mapname, fig, ax0, ax1, ax2, ax3, x_lim, y_lim, jk)
    elif plot == 'hit':
        makeMap(x, y, hitData, hitData, mapname, fig, ax0, ax1, ax2, ax3, x_lim, y_lim, jk)
    elif plot == 'rms':
        makeMap(x, y, rmsData, hitData, mapname, fig, ax0, ax1, ax2, ax3, x_lim, y_lim, jk)

    plt.show(fig)
        
    
    outfile += '.png'
    plt.savefig(outfile)
    plt.close(fig)

    
def makeMap(x, y, plotData, hits, mapname, fig, ax0, ax1, ax2, ax3, x_lim, y_lim, jk):
    plotData = np.ma.masked_where(hits < 1., plotData)
    #plotData[0,:,:] = np.ma.masked_where(hits[0,:,:] < 1, plotData[0,:,:])
    #plotData[1,:,:] = np.ma.masked_where(hits[1,:,:] < 1, plotData[1,:,:])
    #plotData[2,:,:] = np.ma.masked_where(hits[2,:,:] < 1, plotData[2,:,:])
    #plotData[3,:,:] = np.ma.masked_where(hits[2,:,:] < 1, plotData[3,:,:])
    if set_xlim == False:
        dx = x[1] - x[0]
        x_lim[0] = x[0] - 0.5*dx; x_lim[1] = x[-1] + 0.5*dx
    if set_ylim == False:
        dy = y[1] - y[0]
        y_lim[0] = y[1] - 0.5*dy; y_lim[1] = y[-1] + 0.5*dy

    axlist = [ax0, ax1, ax2, ax3]

    
    '''if set_colorlim == False:
        if plot == 'map':
            color_lim[1] = 0.1*np.amax(plotData[2,:,:])*1e6
            color_lim[0] = -color_lim[1]
    '''
 
    cmap = cm.rainbow.set_bad('w',1) #color of masked elements 
    if plot == 'hit':
        im = ax0.imshow(plotData[0,:,:], extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax1.imshow(plotData[1,:,:], extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax2.imshow(plotData[2,:,:], extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax3.imshow(plotData[3,:,:], extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
    else:
        im = ax0.imshow(plotData[0,:,:]*1e6, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax1.imshow(plotData[1,:,:]*1e6, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax2.imshow(plotData[2,:,:]*1e6, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])
        im = ax3.imshow(plotData[3,:,:]*1e6, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest', aspect='equal', cmap=cm.rainbow, origin='lower', vmin=color_lim[0], vmax=color_lim[1])


    #if deepx:
    #    ax.set_ylabel('Right Ascension [deg]')
    #    ax.set_xlabel('Frequency [GHz]')
    #elif deepy:
    #    ax.set_ylabel('Declination [deg]')
    #    ax.set_xlabel('Frequency [GHz]')
    #else:
    ax0.set_ylabel('Declination [deg]')
    ax2.set_ylabel('Declination [deg]')
    ax2.set_xlabel('Right Ascension [deg]')
    ax3.set_xlabel('Right Ascension [deg]')
    ax0.xaxis.set_visible(False)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    ax3.yaxis.set_visible(False)


    fig.suptitle(mapname, fontsize=20)
    plt.subplots_adjust(top=1, wspace=0.1, hspace=-0.5)
    ax2.set_title('Sum'); ax3.set_title('Diff')
    if jk == 'odde':
        ax0.set_title('Even'); ax1.set_title('Odd')
    elif jk == 'dayn':
        ax0.set_title('Night'); ax1.set_title('Day')
    elif jk == 'half':
        ax0.set_title('First half'); ax1.set_title('Second half')
    elif jk == 'sdbg':
        ax0.set_title('Saddlebags 1+2'); ax1.set_title('Saddlebags 3+4')  
    

    divider0 = make_axes_locatable(ax0)
    divider1 = make_axes_locatable(ax1)
    divider2 = make_axes_locatable(ax2)
    divider3 = make_axes_locatable(ax3)
    #divider = make_axes_locatable((ax1,ax3))
    #cax = divider.append_axes('right', size='5%', pad=0.05)
    #'''
    cbar = plt.colorbar(im, ax=axlist, shrink=0.6)
    if plot == 'hit':
        cbar.set_label('hits')
    else:
        cbar.set_label('uK')
    #'''
    #plt.tight_layout()
    
setup(filename, outfile, x_lim, y_lim, jk)
