import numpy as np
import h5py
import sys
import getopt
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from matplotlib import animation, rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import textwrap
import random
from mpl_toolkits.axes_grid1 import make_axes_locatable
import FileDialog
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size


def usage():
    prefix = " "
    preferredWidth = 150
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=preferredWidth,subsequent_indent=' '*len(prefix))
    m1 = "      (Side band as a list, ie. [1,2,4]. Default all)"
    m2 = "      (Frequency, given as list, ie. [1,6,26]. Default all) "
    m3 = "   (Plot type in a list, ie [map]. Choices are map, rms, map/rms, sim, rms_sim and hit. Default map) "
    m4 ="(Filename)"
    m5 ="(color_limits of colorbar) as nested list. First plot --> first list in main list. Default none)"
    m6 = "(Which detector as a list, ie. [4,11,18]. Default all)"
    m7 = "(Outfile name, default 'outfile')"
    m8 = "(x-range)"
    m9 = "(y-range)"
    m10 = "(Include to split the map as the focal plane)"
    m11 = "(If the field non-stationary)"
    m12 = "(Which simulation, default first element (0))"
    m13 = "(Scale the data by a factor w, default 1)"

    print "\nThis is the usage function\n"
    print "Flags:"
    print "-s ----> optional --sb", wrapper.fill(m1)
    print "-n ----> optional --nu", wrapper.fill(m2)
    print "-p ----> optional --plots", wrapper.fill(m3)
    print "-f ----> optional --filename", wrapper.fill(m4)
    print "-c ----> optional --colorlim", wrapper.fill(m5)
    print "-d ----> optional -- det", wrapper.fill(m6)
    print "-o ----> optional -- out", wrapper.fill(m7)
    print "-x ----> optional -- xlim", wrapper.fill(m8)
    print "-y ----> optional -- yim", wrapper.fill(m9)
    print "-t ----> optional --split", wrapper.fill(m10)
    print "-j ----> optional --jupiter", wrapper.fill(m11)
    print "-z ----> optional --simulation", wrapper.fill(m12)
    print "-w ----> optional --scale", wrapper.fill(m13)
    sys.exit()


#plt.switch_backend('agg')

        
plot_choices = ["map", "rms", "map/rms", "hit", "sim", "rms_sim", "feed", "var"]
plot_list = ["map"]
freq = "all"
set_xlim = False
set_ylim  = False
x_lim = []
y_lim = []
sb_list    = range(1,5)
freq_list = range(1,65)
plots = ["map"]
color_lim = [None, None]
set_colorlim = False
det_list = range(1,20)
outfile = "outfile"
split = False
jupiter = False
sim_numb = 0
scale = 1

if len(sys.argv) == 1:
    usage()
    

try:
    opts, args = getopt.getopt(sys.argv[1:],"s:n:p:f:h:c:d:o:x:y:t:j:z:w:", ["sb=", "nu=", "plots=", "filename=", "help=", "colorlim", "det", "out", "xmax", "ymax","split","jupiter", "simulation","scale"])
except getopt.GetoptError:
    usage()

for opt, arg in opts:
    if opt in ("-x", "--xlim"):
        set_xlim = True
        x_lim = eval(arg)
        if type(x_lim) != list:
            print "x-lim needs to be a list with [min, max]"
            sys.exit()
    elif opt in ("-y", "--ylim"):
        set_ylim = True
        y_lim = eval(arg)
        if type(y_lim) != list:
            print "y-lim needs to be a list with [min, max]"
            sys.exit()
    elif opt in ("-o", "--out"):
        outfile = arg
    elif opt in ("-d", "--det"):
        det_list = eval(arg)
        if type(det_list) != list:
            print "Detectors my be inserted as a list, ie. -d [1,2,5,7]"
            sys.exit()
        else:
            if 0 in det_list:
                print "Use 1-base, not 0-base please"
                sys.exit()
    elif opt in ("-s", "--sb"):
        sb_list = eval(arg)
        if type(sb_list) != list:
            print "Side bands my be inserted as a list, ie. -d [1,2,4]"
            sys.exit()
        else:
            if 0 in sb_list:
                print "Use 1-base, not 0-base please"
                sys.exit()
    elif opt in ("-n", "--nu"):
        freq_list = eval(arg)
        if type(freq_list)!= list:
            print "Frequencies my be inserted as a list, ie. -n [1,34,50]"
        else:
            if 0 in freq_list:
                print "Use 1-base, not 0-base please"
                sys.exit()
            for i in freq_list:
                if i > 64:
                    print "There are only 64 frequencies pr. side band"
                    sys.exit()
    elif opt in ("-p", "--plots"):
        plot = arg
        if plot in plot_choices:
            plot_list = eval(arg)
        else:
            option = ""
            plot_list = []
            for i in arg:
                if i == "[":
                    continue
                elif i == "]":
                    plot_list.append(option)
                elif i == ",":
                    plot_list.append(option)
                    option = ""
                else:
                    option += i
        for i, op in enumerate(plot_list):
            if op in plot_choices:
                if op == "map/rms":
                    plot_list[i] = "map_rms"
            else:
                print "Make sure you have chosen the correct plot choices"
                sys.exit()
    elif opt in ("-f", "--filename"):
        filename = arg
    elif opt in ("-h", "--help"):
        usage()
    elif opt in ("-c", "--colorlim"):
        set_colorlim = True
        color_lim = eval(arg)
    elif opt in ("-t", "--split"):
        split = True
    elif opt in ("-j", "--jupiter"):
        jupiter = True
    elif opt in ("-z", "--simulation"):
        sim_numb = int(arg)
    elif opt in ("-w", "--scale"):
        scale = eval(arg)
    else:
        usage()

if set_colorlim == False:
    color_lim = []
    for i in range(len(plot_list)):
        color_lim.append([None, None])
        
def readMap(filename):
        dfile   = h5py.File(filename,'r')
        sim_exist = "map_sim" in dfile 

        nx      = dfile['n_x'];      nx  = np.array(nx).astype(int)
        ny      = dfile['n_y'];      ny  = np.array(ny).astype(int)
        x       = dfile['x'];         x  = np.array(x[:]).astype(float)
        y       = dfile['y'];         y  = np.array(y[:]).astype(float)
        maps    = dfile['map'];    maps  = np.array(maps[...]).astype(float)
        hit     = dfile['nhit'];    hit  = np.array(hit[...]).astype(float)
        rms     = dfile['rms'];     rms  = np.array(rms[...]).astype(float)
        if sim_exist:
            maps_sim = dfile['map_sim'];  maps_sim = np.array(maps_sim[...]).astype(float) 
            rms_sim  = dfile['rms_sim'];  rms_sim  = np.array(rms_sim[...]).astype(float)
        else:
            maps_sim = []
            rms_sim = []

        return x, y, maps, hit, rms, maps_sim, rms_sim

def swap_xy(my_array):
        new_array = np.zeros((len(my_array[0]), len(my_array)))
        for i in range(len(my_array[0])):
            new_array[i, :] = my_array[:, i]
        return new_array

def setup_one(filename, outfile, x_lim, y_lim):

        print filename.split("/")[-1]
        x, y, maps, hit, rms, maps_sim, rms_sim = readMap(filename)
        nx = len(x)
        ny = len(y)
        ndet = maps.shape[0]
        nsb = maps.shape[1]
        nfreq = maps.shape[2]

        if len(plot_list) == 1:
            fig, ax = plt.subplots(1)
            axarr = [ax]
        else:
            fig, axarr = plt.subplots(len(plot_list))

        for p, plot in enumerate(plot_list):
            mapData = np.zeros((nx,ny))
            hitData = np.zeros((nx,ny))
            sb_nr = str()
            if plot   == "hit":
                data  = hit
            elif plot == "rms":
                data  = rms
            elif plot == "map_rms" or plot == "var":
                data  = maps/rms
            elif plot == "sim":
                if len(maps_sim) == 0:
                    print "Could not find any simulated data."
                    sys.exit()
                else:
                    data = maps_sim[sim_numb,:,:,:,:]
            elif plot == "rms_sim":
                if len(rms_sim) == 0:
                    print "Could not find any simulated data."
                    sys.exit()
                else:
                    data = rms_sim[sim_numb,:,:,:,:]
            else:
                #print(np.shape(maps))
                data  = maps
            
            #outfile = outfile
            #makeMap_gif(x, y, maps, hit, outfile)
            #plt.close(fig)
            #ss


            mapData = np.zeros((ny, nx))
            hitData = np.zeros((ny, nx))
            rms_m1_sum    = np.zeros((ny, nx))
            rms_m1        = np.zeros((ny, nx))
            count = 0
            seenbyfeed = np.zeros((nx,ny))
            hit_temp = np.zeros((19,ny,nx))


            for i in det_list:
                #print "det %i" %i
                for j in sb_list:
                    for k in freq_list:
                        map_check = np.sum(data[i-1,j-1,k-1], axis=-1)
                        map_check = np.sum(map_check)
                        if map_check == 0:
                            count += 1
                            continue
                        else:
                            if plot == "map":
                                rms_m1        = np.ma.fix_invalid(1./rms[i-1,j-1,k-1]**2, fill_value=0)
                                rms_m1_sum   += rms_m1
                                mapData     += np.ma.fix_invalid(data[i-1,j-1,k-1]/rms[i-1,j-1,k-1]**2, fill_value=0)
                            else:
                                mapData      += data[i-1,j-1,k-1]
                                rms_m1_sum = 1.
                            hitData     += hit[i-1,j-1,k-1]
                            hit_temp[i-1,:,:]    += hit[i-1,j-1,k-1,:,:]


            #mapData /= ((ndet-2)*nsb*nfreq-count)
        
            mapData = np.ma.fix_invalid(mapData/rms_m1_sum, fill_value=0)
            
            mapData = swap_xy(mapData) 
            hitData = swap_xy(hitData) 

            #print np.max(mapData)
            mapname = 'd'

            if plot == "feed":
                for i in det_list:
                    for m in range(nx):
                        for n in range(ny):
                            if hit_temp[i-1,n,m] > 0.01*hitData[m,n]:
                                seenbyfeed[m,n] = seenbyfeed[m,n] + 1
                mapData = seenbyfeed

            
            makeMap(x, y, np.transpose(mapData), np.transpose(hitData), mapname, axarr, p, x_lim, y_lim)


        outfile += ".png"
        plt.savefig(outfile)#,bbox_inches='tight')                                             
        #print "Plot saved as", outfile.split("/")[-1]
        #plt.show()
        plt.close(fig)
        #print 'Done!'


def setup_split(filename, outfile, x_lim, y_lim):
    print "Focal plane map"
    x, y, maps, hit, rms, maps_sim, rms_sim = readMap(filename)
    nx = len(x)
    ny = len(y)
    ndet = maps.shape[0]
    nsb = maps.shape[1]
    nfreq = maps.shape[2]

    tweek_beam = 1.#5.
    tweek = 1.#2.
    
    # read offsetfile
    offset_file = '/mn/stornext/d16/cmbco/comap/protodir/auxiliary/Receiver_offset.dat'
    det_nr, dxx, dyy, dAz, dEl = np.loadtxt(offset_file, unpack=True, comments='#')
    
    # center pixels and resolution of grid
    cx = int(nx/2); cy = int(ny/2)
    dx = x[1]-x[0]; dy = abs(y[1]-y[0])

    x = (x - x[cx]*np.ones(nx))
    y = (y - y[cy]*np.ones(ny))

    # pixels in the beam
    beam = tweek_beam*5./60.; beam_x = int(beam/dx); beam_y = int(beam/dy)
    #beam_x = 2; beam_y = 4
    rad_x = int(beam_x/2.)+1; rad_y = int(beam_y/2)+1.
    cx1 = max(min(cx - rad_x,nx-beam_x),1); cy1 = max(min(cy - rad_y,ny-beam_y),1)
    cx2 = max(min(cx + rad_x,nx-beam_x),1); cy2 = max(min(cy + rad_y,ny-beam_y),1)
    
    #print cx, cy, beam_x, beam_y

    if len(plot_list) == 1:
        fig, ax = plt.subplots(1)
        axarr = [ax]
    else:
        fig, axarr = plt.subplots(len(plot_list))

    for p, plot in enumerate(plot_list):
        #mapData = np.zeros((nx,ny))
        #hitData = np.zeros((nx,ny))
        #sb_nr = str()
        if plot   == "hit":
            data  = hit
        elif plot == "rms":
            data  = rms
        elif plot == "map_rms":
            data  = maps/rms
        elif plot == "sim":
            if len(maps_sim) == 0:
                print "Could not find any simulated data."
                sys.exit()
            else: 
                data = maps_sim[sim_numb,:,:,:,:,:]
        elif plot == "rms_sim":
            if len(rms_sim) == 0:
                print "Could not find any simulated data."
                sys.exit()
            else:
                data = rms_sim[sim_numb,:,:,:,:,:]
        else:
            data  = maps


        mapData = np.zeros((ny, nx))
        hitData = np.zeros((ny, nx))
        count = 0
        
        for i in det_list:
            #print "det %i" %i
            dRA = int(tweek*dAz[i-1]/60./dx)-rad_x; dDec = int(tweek*dEl[i-1]/60./dy)-rad_y
            #print "offsets", dRA, dDec
            for j in sb_list:
                for k in freq_list:
                    map_check = np.sum(data[i-1,j-1,k-1], axis=-1)
                    map_check = np.sum(map_check)
                    if map_check == 0:
                        count += 1
                        continue
                    else:
                        mapData[cy1-dDec:cy2-dDec,cx1-dRA:cx2-dRA] += data[i-1,j-1,k-1,cy1:cy2,cx1:cx2]
                        hitData[cy1-dDec:cy2-dDec,cx1-dRA:cx2-dRA] += hit[i-1,j-1,k-1,cy1:cy2,cx1:cx2]
        mapData /= ((len(det_list))*len(sb_list)*len(freq_list)-count)
                        
        mapData = np.ma.fix_invalid(mapData, fill_value=0)

        mapData = swap_xy(mapData) 
        hitData = swap_xy(hitData) 
        
        mapname = 'd'
        #makeMap(x, y, np.transpose(mapData), np.transpose(hitData), mapname, axarr, p, x_lim, y_lim)
        makeMap(x, y, mapData, hitData, mapname, axarr, p, x_lim, y_lim)
        outfile += ".png"
        plt.savefig(outfile)#,bbox_inches='tight')                                             
        print "Plot saved as", outfile
        #plt.show()
        plt.close(fig)
        #print 'Done!'





def makeMap(x, y, mapData, hitData, mapname, axarr,p, x_lim, y_lim):
        mapData = np.ma.masked_where(hitData < 1., mapData)
        #mapData = mapData
        if set_xlim == False:
            x_min = x[1]
            x_max = x[-1] + (x[1]-x[0])
            x_lim = [x_min, x_max]
        if set_ylim == False:
            y_min = y[1]
            y_max = y[-1] + (y[1]-y[0]) 
            y_lim = [y_min, y_max]

        x_cp = int(len(x)/2.); y_cp = int(len(y)/2.)#; print len(x), len(y)
        x_mid = x[x_cp]; x_aim = x_mid - 2.
        y_mid = y[y_cp]; y_aim = y_mid - 2.


        if set_colorlim == False:
            if plot_list[p] == 'map':
                color_lim[p][1] = 0.1*np.amax(mapData)
                color_lim[p][0] = -color_lim[p][1]
            else:
                color_lim[p][:] = [None, None]
            #if plots[p] == 'hit':
            #    color_lim[p][1] = 0.9*np.amax(mapData)
            #    color_lim[p][0] = 0.

        if split:
            cmap = cm.gist_gray.set_bad('k',1) #color of masked elements
            im = axarr[p].imshow(mapData, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest',cmap=cm.gist_gray, origin='lower', vmin=color_lim[p][0], vmax=color_lim[p][1])
            axarr[p].set_ylabel(r'$\Delta$deg')
            axarr[p].set_xlabel(r'$\Delta$deg')
            divider = make_axes_locatable(axarr[p])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('K')
        elif jupiter:
            cmap = cm.viridis.set_bad('k',1)
            im = axarr[p].imshow(mapData/0.7266, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest',cmap=cm.viridis, origin='lower', vmin=color_lim[p][0], vmax=color_lim[p][1])
            axarr[p].set_ylabel(r'$\Delta$El [deg]')
            axarr[p].set_xlabel(r'$\Delta$RA [deg]')            
            divider = make_axes_locatable(axarr[p])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('K')
        else:
            cmap = cm.CMRmap.set_bad('0.8',1) #color of masked elements
            im = axarr[p].imshow(mapData*scale, interpolation='nearest', aspect='equal', cmap=cm.CMRmap, origin='lower', vmin=color_lim[p][0], vmax=color_lim[p][1])
            #im = axarr[p].imshow(hitData, extent=(x_lim[0],x_lim[1],y_lim[0],y_lim[1]), interpolation='nearest',cmap=cm.CMRmap, origin='lower', vmin=color_lim[p][0], vmax=color_lim[p][1])
            axarr[p].set_ylabel('Declination [deg]')
            axarr[p].set_xlabel('Right Ascension [deg]')
            #axarr[p].set_ylabel('Elevation [deg]')
            #axarr[p].set_xlabel('Azimuth [deg]')
	    axarr[p].set_title(plot_list[p])
            divider = make_axes_locatable(axarr[p])
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cbar = plt.colorbar(im, cax=cax)
            if plot_list[p] == 'feed':
                cbar.set_label('feeds')
            if plot_list[p] == 'hit':
                cbar.set_label('hits')
            if plot_list[p] == 'map':
                if (scale != 1):
                    cbar.set_label(str(scale)+'K')
                else:
                    cbar.set_label('K')
#            rect = patches.Rectangle((x_aim,y_aim),4.,4.,linewidth=2,edgecolor='w',facecolor='none')
#            axarr[p].add_patch(rect)

        #im = axarr[p].imshow(mapData, extent=(x_min,x_max,y_min,y_max), interpolation='nearest',cmap=cm.gist_gray, origin='lower', vmin=color_lim[p][0], vmax=color_lim[p][1])

        #print 'std:', np.std(mapData[mapData!=0])


def makeMap_gif(x, y, maps, hitData, mapname):

        cm.viridis.set_bad('w',1.) #color of masked elements                                                                                    
        nsb = 4
        nfreq = 12

        if set_xlim == False:
            x_min = x[1]
            x_max = x[-1] + (x[1]-x[0])
            x_lim = [x_min, x_max]
        if set_ylim == False:
            y_min = y[1]
            y_max = y[-1] + (y[1]-y[0]) 
            y_lim = [y_min, y_max]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.ylabel('Declination [deg]')
        plt.xlabel('Right Ascension [deg]')

        frame = []
        ims = []
        for sb in range(nsb):
            for freq in range(nfreq):
                mapData = np.nanmean(maps[:, sb,freq,:,:], axis=0)  
                #mapData = np.ma.masked_where(np.nanmean(hitData[:,sb,freq], axis=0) < 1., mapData)
                                            
                frame.append(mapData)
        cv0 = frame[0]
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        im = ax.imshow(cv0, origin='lower') # Here make an AxesImage rather than contour                                                        
        cb = fig.colorbar(im,cax=cax)
        tx = ax.set_title('Frame 0')
        
        
        def animate(i):
            arr = frame[i]
            im.set_data(arr)
            im.set_clim(0, 1)
            tx.set_text('Frame {0}'.format(i))
            
        ani = animation.FuncAnimation(fig, animate, frames=32)
        #plt.show()
        print "Gif saved as ", mapname

        ani.save(mapname,writer='imagemagick', fps=2)


if split:
    setup_split(filename, outfile, x_lim, y_lim)
else:
    setup_one(filename, outfile, x_lim, y_lim)
