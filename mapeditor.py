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


class Atlas:
    def __init__(self):
        self.map_mode_choices = ["map", "rms", "map/rms", "hit", "var"]
        self.map_mode         = "map"
        self.jk_choices   = ["odde", "dayn", "half", "sdbg"]
        self.jk           = "odde"
        self.tool_choices   = ["add", "coadd", "sub", "divide", "coaddfeed"]
        self.tool           = "add"
        self.freq         = "all"
        self.det_list     = np.arange(1,20)
        self.sb_list      = np.arange(1,5)
        self.freq_list    = np.arange(1,65)
        self.outfile      = None
        self.scale        = 1
        self.beam        = False
        self.patch1       = ''
        self.patch2       = ''
        self.rms_lim      = 200000.
        self.deepx        = False
        self.deepy        = False
        self.infile1      = None
        self.infile2      = None
        self.input()        
        if self.infile1 != None and self.infile2 == None:
            self.data     = self.readMap(self.infile1)
        elif self.infile1 != None and self.infile2 != None:
            self.data1        = self.readMap(self.infile1)
            self.data2        = self.readMap(self.infile2)
        
        if self.outfile != None:
            self.writeMap(self.outfile)
        else:
            print("To save result to outfile, please provide a valid outfile name!")

    def input(self):
        if len(sys.argv) == 1:
            self.usage()

        try:
            opts, args = getopt.getopt(sys.argv[1:],"s:f:m:i:h:d:o:I:r:l:j:bt:", ["sb=", "freq=", "mode=", "infile1=", "help=", "det",
                                                                                 "out","infile2=","deepx", "deepy", "jk", "beam", "tool"])
        except getopt.GetoptError:
            self.usage()

        for opt, arg in opts:
            if opt in ("-j", "--jk"):
                self.jk = arg
                if self.jk not in self.jk_choices:
                    print("Make sure you have chosen the correct jk choices")                                                                                                   
                    sys.exit() 
            elif opt in ("-b", "--beam"):
                self.beam = True
            elif opt in ("-o", "--out"):
                self.outfile = arg
            elif opt in ("-d", "--det"):
                self.det_list = eval(arg)
                if type(self.det_list) != list and type(self.det_list) != np.ndarray:
                    print("Detectors my be inserted as a list or array, ie. -d [1,2,5,7]")
                    sys.exit()
                else:
                    if 0 in self.det_list:
                        print("Use 1-base, not 0-base please")
                        sys.exit()
                self.det_list = np.array(self.det_list, dtype = int)
            elif opt in ("-s", "--sb"):
                self.sb_list = eval(arg)
                if type(self.sb_list) != list and type(self.sb_list) != np.ndarray:
                    print("Side bands my be inserted as a list or array, ie. -d [1,2,4]")
                    sys.exit()
                else:
                    if 0 in self.sb_list:
                        print("Use 1-base, not 0-base please")
                        sys.exit()
                self.sb_list = np.array(self.sb_list, dtype = int)
            elif opt in ("-f", "--freq"):
                self.freq_list = eval(arg)
                if type(self.freq_list)!= list and type(self.freq_list) != np.ndarray:
                    print("Frequencies my be inserted as a list or array, ie. -n [1,34,50]")
                else:
                    if 0 in self.freq_list:
                        print("Use 1-base, not 0-base please")
                        sys.exit()
                    if np.any(np.array(self.freq_list) > 64):
                        print("There are only 64 frequencies pr. side band")
                        sys.exit()
                self.freq_list = np.array(self.freq_list, dtype = int)
            elif opt in ("-m", "--mode"):
                self.map_mode = arg
                if self.map_mode not in self.map_mode_choices:
                    print("Make sure you have chosen the correct map_mode choices")                                                                                                   
                    sys.exit()
            elif opt in ("-i", "--infile1"):
                self.infile1 = arg
                temp = self.infile1.split('/')[-1]
                self.patch1 = temp.split('_')[0]
            elif opt in ("-I", "--infile2"):
                self.infile2 = arg
                temp = self.infile1.split('/')[-1]
                self.patch2 = temp.split('_')[0]
            elif ("-i" in opt or "--infile1" in opt) and ("-I" in opt or "--infile2" in opt):
                if self.patch1 != self.patch2:
                    print("Can only perform operations on two maps if they belong to the same sky-patch.") 
                    sys.exit()
            elif opt in ("-h", "--help"):
                self.usage()
            elif opt in ("-z", "--sim"):
                self.sim_numb = int(arg)
            elif opt in ("-w", "--scale"):
                self.scale = eval(arg)
            elif opt in ("-m", "--mask"):
                self.rms_lim = eval(arg)
            elif opt in ("-r", "--deepx"):
                self.deepx = True
                self.y_index = int(arg)-1
            elif opt in ("-l", "--deepy"):
                self.deepy = True
                self.x_index = int(arg)-1
            else:
                self.usage()

    def usage(self):
        prefix = ""
        preferredWidth = 150
        wrapper = textwrap.TextWrapper(initial_indent=prefix, width=preferredWidth,subsequent_indent=' '*len(prefix))
        m1 = "(Side band as a list, ie. [1,2,4]. Default all)"
        m2 = "(Frequency, given as list, ie. [1,6,26]. Default all) "
        m3 = "(Type of map_mode. Choices are map, rms, map/rms, sim, rms_sim, hit, and var. Default map) "
        m4 = "(Filename)"
        m5 = "(color_limits of colorbar) as nested list. First map_mode --> first list in main list. Default none)"
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
        print("-p ----> optional --map_modes "    + wrapper.fill(m3))
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

    def readMap(self, infile):
        dfile   = h5.File(infile,'r')
        #self.sim_exist = "map_sim" in dfile
        #self.x       = dfile['x'][:]
        #self.y       = dfile['y'][:]
        #self.freq    = dfile['freq'][self.sb_list - 1, self.freq_list - 1]
        freq_start = self.freq_list[0] - 1
        freq_end   = self.freq_list[-1] - 1
        sb_start = self.sb_list[0] - 1
        sb_end   = self.sb_list[-1] - 1
        if self.beam:
            if self.map_mode == "map":
                print(self.sb_list - 1)
                data =  dfile['map_beam'][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                #data =  dfile['map_beam'][0, 10, ...]
            
            elif self.map_mode == "hit":      
                data =  dfile['nhit_beam'][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            
            elif self.map_mode == "rms":
                data =  dfile['rms_beam'][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
        else:
            if self.map_mode == "map":
                data = dfile['map'][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            elif self.map_mode == "hit":
                data  = dfile['nhit'][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            elif self.map_mode == "rms":
                data  = dfile['rms'][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
        dfile.close()
        return data

    def writeMap(self, outfile):
        dfile   = h5.File(outfile,'w')
        dfile.create_dataset(self.map_mode, data = self.data)
        dfile.close()

if __name__ == "__main__":
    map1 = Atlas()
    data = map1.data
    start = time.time()
    data = np.sum(data + data, axis = 0) / np.sum(data, axis = 0)
    print("Run time: ", time.time() - start, " sek ", data.shape)





    
