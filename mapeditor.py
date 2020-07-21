import time
import sys
import getopt
import numpy as np
import h5py as h5
import textwrap
#import FileDialog
import ctypes

class Atlas:
    def __init__(self):
        self.jk_choices   = ["odde", "dayn", "half", "sdlb"]
        self.jk           = ["odde"]
        self.jack         = False
        self.tool_choices   = ["coadd", "subtract"]
        self.tool           = "coadd"
        self.freq         = "all"
        self.det_list     = np.arange(1,20)
        self.sb_list      = np.arange(1,5)
        self.freq_list    = np.arange(1,65)
        self.outfile      = None
        self.scale       = None
        self.beam        = False
        self.full        = False
        self.everything  = False 
        self.patch1       = ''
        self.patch2       = ''
        self.rms_lim      = 200000.
        self.deepx        = False
        self.deepy        = False
        self.maputilslib = ctypes.cdll.LoadLibrary("maputilslib.so.1")  # Load shared library
        self.infile1      = None
        self.infile2      = None
        self.access       = "a"
        self.input()    
        if self.outfile == None:
            print("To save result to outfile, please provide a valid outfile name!")
            sys.exit() 
        if self.jack and ("jackknives" not in self.dfile1 or "jackknives" not in self.dfile2):
            print("One or both of the input files does not contain any jackknife information!")
            sys.exit()
        if not self.full and not self.beam and not self.jack:
            self.everything = True
            if "jackknives" in self.dfile1 and "jackknives" in self.dfile2:
                nhit_lst = [i for i in self.dfile1["jackknives"].keys() if "nhit_" in i]
                self.jk =  [i.split("_")[1] for i in nhit_lst]

        self.operation()
        self.dfile1.close()
        self.dfile2.close()
        self.ofile.close()

    def input(self):
        if len(sys.argv) == 1:
            self.usage()

        try:
            opts, args = getopt.getopt(sys.argv[1:],"s:f:i:h:d:o:I:r:l:j:t:bw:a::F", ["sb=", "freq=", "infile1=", "help=", "det",
                                                                                 "out","infile2=","deepx", "deepy", "jk", "tool", "beam", "scale","access", "full"])
        except getopt.GetoptError:
            self.usage()

        for opt, arg in opts:
            if opt in ("-j", "--jk"):
                self.jack = True
                self.jk = arg.split(",")
                self.jk = list(self.jk)
                for jk in self.jk:
                    if jk not in self.jk_choices:
                        print("Make sure you have chosen the correct jk choices")                                                                                                   
                        sys.exit() 
            elif opt in ("-t", "--tool"):
                self.tool = arg
                if self.tool not in self.tool_choices:
                    print("Make sure you have chosen the correct tool choices")                                                                                                   
                    sys.exit() 
            elif opt in ("-b", "--beam"):
                self.beam = True
            elif opt in ("-F", "--full"):
                self.full = True
            elif opt in ("-a", "--overwrite"):
                self.access = arg
            elif opt in ("-o", "--out"):
                self.outfile = arg
                self.ofile = h5.File(self.outfile, self.access)
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
            elif opt in ("-i", "--infile1"):
                self.infile1 = arg
                self.dfile1  = h5.File(self.infile1,'r')
                temp = self.infile1.split('/')[-1]
                self.patch1 = temp.split('_')[0]
            elif opt in ("-I", "--infile2"):
                self.infile2 = arg
                self.dfile2        = h5.File(self.infile2,'r')
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
        m3 = "(Type of map_mode. Choices are map, rms, map/rms, sim, rms_sim, nhit, and var. Default map) "
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

    def readMap(self, first_file = True, jackmode = None):
        if first_file:
            dfile = self.dfile1
        else:
            dfile = self.dfile2

        freq_start = self.freq_list[0] - 1
        freq_end   = self.freq_list[-1] - 1
        sb_start = self.sb_list[0] - 1
        sb_end   = self.sb_list[-1] - 1

        if jackmode != None:
            if jackmode == "dayn" or jackmode == "half":
                map  =  dfile["jackknives/map_" + jackmode][:, self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                nhit =  dfile["jackknives/nhit_" + jackmode][:, self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                rms  =  dfile["jackknives/rms_" + jackmode][:, self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            else: 
                map  =  dfile["jackknives/map_" + jackmode][:, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                nhit =  dfile["jackknives/nhit_" + jackmode][:, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                rms  =  dfile["jackknives/rms_" + jackmode][:, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
        else:
            if self.beam:
                map =  dfile["map_beam"][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                nhit =  dfile["nhit_beam"][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                rms =  dfile["rms_beam"][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            elif self.full:
                map =  dfile["map"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                nhit =  dfile["nhit"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                rms =  dfile["rms"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
        return map, nhit, rms
        
    def writeMap(self, jackmode = None, write_the_rest = False):
        if jackmode != None:
            map_name    = "jackknives/map_" + jackmode
            nhit_name   = "jackknives/nhit_" + jackmode
            rms_name    = "jackknives/rms_" + jackmode
            if map_name in self.ofile and self.access == "a":
                map_data        = self.ofile[map_name]
                nhit_data       = self.ofile[nhit_name]
                rms_data        = self.ofile[rms_name]
                map_data[...]   = self.map
                nhit_data[...]  = self.nhit
                rms_data[...]   = self.rms
            else:
                self.ofile.create_dataset(map_name, data = self.map)
                self.ofile.create_dataset(nhit_name, data = self.nhit)
                self.ofile.create_dataset(rms_name, data = self.rms)

        elif self.beam or self.full:
            if self.beam:
                map_name    = "map_beam"
                nhit_name   = "nhit_beam"
                rms_name    = "rms_beam"
            elif self.full: 
                map_name    = "map"
                nhit_name   = "nhit"
                rms_name    = "rms"
                
            if map_name in self.ofile and self.access == "a":
                map_data        = self.ofile[map_name]
                nhit_data       = self.ofile[nhit_name]
                rms_data        = self.ofile[rms_name]
                map_data[...]   = self.map
                nhit_data[...]  = self.nhit
                rms_data[...]   = self.rms
            else:
                self.ofile.create_dataset(map_name, data = self.map)
                self.ofile.create_dataset(nhit_name, data = self.nhit)
                self.ofile.create_dataset(rms_name, data = self.rms)
        if write_the_rest:
            data_not_to_copy = ["jackknives", "map", "map_beam", "nhit", 
                                "nhit_beam", "rms", "rms_beam"]
            jk_data_not_to_copy = ["map_dayn",  "map_half",  "map_odde",  "map_sdlb",
                                   "nhit_dayn", "nhit_half", "nhit_odde", "nhit_sdlb",
                                   "rms_dayn",  "rms_half",  "rms_odde",  "rms_sdlb"]
            for name in self.dfile1.keys():
                if name not in self.ofile.keys() and name not in data_not_to_copy:
                    self.ofile.create_dataset(name, data = self.dfile1[name])    
            
            if "jackknives" in self.dfile1 and "jackknives" in self.dfile2 and "jackknives" in self.ofile:
                for name in self.dfile1["jackknives"].keys():
                    if name not in self.ofile["jackknives"].keys() and name not in jk_data_not_to_copy:
                        self.ofile.create_dataset("jackknives/" + name, 
                                                    data = self.dfile1["jackknives/" + name])   

    def operation(self):                
        if self.infile1 != None and self.infile2 != None:
            if self.everything:
                if "jackknives" in self.dfile1 and "jackknives" in self.dfile2:
                    for jack in self.jk:
                        self.map1, self.nhit1, self.rms1 = self.readMap(True, jack)
                        self.map2, self.nhit2, self.rms2 = self.readMap(False, jack)
                        
                        if self.tool == "coadd":
                            if len(self.map1.shape) == 6:
                                self.C_coadd6D(self.map1, self.nhit1, self.rms1,
                                              self.map2, self.nhit2, self.rms2)  
                    
                            elif len(self.map1.shape) == 5: 
                                self.C_coadd5D(self.map1, self.nhit1, self.rms1,
                                              self.map2, self.nhit2, self.rms2)  

                        elif self.tool == "subtract": 
                            self._map   = self.subtract(self.map1, self.map2)
                            self._nhit   = self.subtract(self.nhit1, self.nhit2)
                            self._rms   = self.add_rms(self.rms1, self.rms2)
                            if len(self.map1.shape) == 6:
                                self.C_subtract6D(self.map1, self.nhit1, self.rms1,
                                              self.map2, self.nhit2, self.rms2)  
                            
                            if len(self.map1.shape) == 5:
                                self.C_subtract5D(self.map1, self.nhit1, self.rms1,
                                              self.map2, self.nhit2, self.rms2) 
                            print(np.allclose(self._map, self.map))
                            print(np.allclose(self._nhit, self.nhit))
                            print(np.allclose(self._rms, self.rms))
                        self.writeMap(jack)
                
                self.full = True
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                self.map2, self.nhit2, self.rms2 = self.readMap(False)
                
                if self.tool == "coadd":
                    self.C_coadd5D(self.map1, self.nhit1, self.rms1,
                                   self.map2, self.nhit2, self.rms2)

                elif self.tool == "subtract": 
                    self.map   = self.subtract(self.map1, self.map2)
                    self.nhit   = self.subtract(self.nhit1, self.nhit2)
                    self.rms   = self.add_rms(self.rms1, self.rms2)
                self.writeMap()
                
                self.full = False
                self.beam = True
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                self.map2, self.nhit2, self.rms2 = self.readMap(False)
                
                if self.tool == "coadd":
                    self.C_coadd4D(self.map1, self.nhit1, self.rms1,
                                   self.map2, self.nhit2, self.rms2)

                elif self.tool == "subtract": 
                    self.map   = self.subtract(self.map1, self.map2)
                    self.nhit   = self.subtract(self.nhit1, self.nhit2)
                    self.rms   = self.add_rms(self.rms1, self.rms2)
                self.writeMap()
            
            elif self.jack:
                for jack in self.jk:
                    self.map1, self.nhit1, self.rms1 = self.readMap(True, jack)
                    self.map2, self.nhit2, self.rms2 = self.readMap(False, jack)
                    
                    if self.tool == "coadd":
                        if len(self.map1.shape) == 6:
                            self.C_coadd6D(self.map1, self.nhit1, self.rms1,
                                           self.map2, self.nhit2, self.rms2)  
                
                        elif len(self.map1.shape) == 5: 
                            self.C_coadd5D(self.map1, self.nhit1, self.rms1,
                                           self.map2, self.nhit2, self.rms2)

                    elif self.tool == "subtract": 
                        if len(self.map1.shape) == 6:
                            self.C_subtract6D(self.map1, self.nhit1, self.rms1,
                                        self.map2, self.nhit2, self.rms2)  
                        
                        if len(self.map1.shape) == 5:
                            self.C_subtract5D(self.map1, self.nhit1, self.rms1,
                                        self.map2, self.nhit2, self.rms2) 
                    self.writeMap(jack)

            elif self.full or self.beam:
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                self.map2, self.nhit2, self.rms2 = self.readMap(False)
                
                if self.tool == "coadd":
                    if self.beam:
                        self.C_coadd4D(self.map1, self.nhit1, self.rms1,
                                       self.map2, self.nhit2, self.rms2)
                    else:
                        self.C_coadd5D(self.map1, self.nhit1, self.rms1,
                                       self.map2, self.nhit2, self.rms2)

                elif self.tool == "subtract": 
                    self.map   = self.subtract(self.map1, self.map2)
                    self.nhit   = self.subtract(self.nhit1, self.nhit2)
                    self.rms   = self.subtract_rms(self.rms1, self.rms2)
                self.writeMap()
        self.writeMap(write_the_rest = True)

    def add(self, data1, data2):
        return data1 + data2
    
    def subtract(self, data1, data2):
        return data1 - data2

    def multiply(self, data, factor):
        return factor * data

    def coadd(self, map1, map2, inv_rms1, inv_rms2):
        inv_var1 = np.square(inv_rms1)
        inv_var2 = np.square(inv_rms2)
        var_inv = inv_var1 + inv_var2
        coadded = map1 * inv_var1 + map2 * inv_var2
        coadded /= var_inv
        return coadded

    def C_coadd4D(self, map1, nhit1, rms1,
                        map2, nhit2, rms2):
        float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
        int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
        self.maputilslib.coadd4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                             float32_array4, int32_array4, float32_array4,
                                             float32_array4, int32_array4, float32_array4,
                                             ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                             ctypes.c_int]
        n0, n1, n2, n3  = self.map1.shape
        self.map        = np.zeros_like(map1,   dtype = ctypes.c_float)
        self.nhit       = np.zeros_like(nhit1,  dtype = ctypes.c_int)
        self.rms        = np.zeros_like(rms1,   dtype = ctypes.c_float)
        self.maputilslib.coadd4D(map1, nhit1, rms1,
                                 map2, nhit2, rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0, n1, n2, n3)
    

    def C_coadd5D(self, map1, nhit1, rms1,
                        map2, nhit2, rms2):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.coadd5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                            float32_array5, int32_array5, float32_array5,
                                            float32_array5, int32_array5, float32_array5,
                                            ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                            ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4 = self.map1.shape
        self.map = np.zeros_like(map1, dtype = ctypes.c_float)
        self.nhit = np.zeros_like(nhit1, dtype = ctypes.c_int)
        self.rms = np.zeros_like(rms1, dtype = ctypes.c_float)
        self.maputilslib.coadd5D(map1, nhit1, rms1,
                                 map2, nhit2, rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0, n1, n2, n3, n4)
                            
    def C_coadd6D(self, map1, nhit1, rms1,
                        map2, nhit2, rms2):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.coadd6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                            float32_array6, int32_array6, float32_array6,
                                            float32_array6, int32_array6, float32_array6,
                                            ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                            ctypes.c_int, ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = self.map1.shape
        self.map = np.zeros_like(map1, dtype = ctypes.c_float)
        self.nhit = np.zeros_like(nhit1, dtype = ctypes.c_int)
        self.rms = np.zeros_like(rms1, dtype = ctypes.c_float)
        self.maputilslib.coadd6D(map1,     nhit1,     rms1,
                                 map2,     nhit2,     rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0,       n1,        n2, 
                                 n3,       n4,        n5)

    def C_subtract4D(self, map1, nhit1, rms1,
                           map2, nhit2, rms2):
        float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
        int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
        self.maputilslib.subtract4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                                float32_array4, int32_array4, float32_array4,
                                                float32_array4, int32_array4, float32_array4,
                                                ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                ctypes.c_int]
        n0, n1, n2, n3  = self.map1.shape
        self.map        = np.zeros_like(map1,   dtype = ctypes.c_float)
        self.nhit       = np.zeros_like(nhit1,  dtype = ctypes.c_int)
        self.rms        = np.zeros_like(rms1,   dtype = ctypes.c_float)
        self.maputilslib.subtract4D(map1, nhit1, rms1,
                                 map2, nhit2, rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0, n1, n2, n3)
    

    def C_subtract5D(self, map1, nhit1, rms1,
                           map2, nhit2, rms2):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.subtract5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                                float32_array5, int32_array5, float32_array5,
                                                float32_array5, int32_array5, float32_array5,
                                                ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3, n4 = self.map1.shape
        self.map = np.zeros_like(map1, dtype = ctypes.c_float)
        self.nhit = np.zeros_like(nhit1, dtype = ctypes.c_int)
        self.rms = np.zeros_like(rms1, dtype = ctypes.c_float)
        self.maputilslib.subtract5D(map1, nhit1, rms1,
                                 map2, nhit2, rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0, n1, n2, n3, n4)
                            
    def C_subtract6D(self, map1, nhit1, rms1,
                           map2, nhit2, rms2):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.subtract6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                                float32_array6, int32_array6, float32_array6,
                                                float32_array6, int32_array6, float32_array6,
                                                ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                ctypes.c_int,   ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = self.map1.shape
        self.map = np.zeros_like(map1, dtype = ctypes.c_float)
        self.nhit = np.zeros_like(nhit1, dtype = ctypes.c_int)
        self.rms = np.zeros_like(rms1, dtype = ctypes.c_float)
        self.maputilslib.subtract6D(map1,     nhit1,     rms1,
                                 map2,     nhit2,     rms2, 
                                 self.map, self.nhit, self.rms,
                                 n0,       n1,        n2, 
                                 n3,       n4,        n5)
                            
    def add_rms(self, rms1, rms2):
        var1 = np.square(rms1)
        var2 = np.square(rms2)
        sum = var1 + var2 
        sum = np.sqrt(sum)
        return sum


if __name__ == "__main__":
    t = time.time()
    map = Atlas()
    print("Run time: ", time.time() - t, " sec")




    
