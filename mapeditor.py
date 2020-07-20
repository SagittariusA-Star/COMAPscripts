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
        self.tool_choices   = ["coadd", "subtract", "divide", "coaddfeed, scale"]
        self.tool           = "coadd"
        self.freq         = "all"
        self.det_list     = np.arange(1,20)
        self.sb_list      = np.arange(1,5)
        self.freq_list    = np.arange(1,65)
        self.outfile      = None
        self.scale       = None
        self.beam        = False
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
        self.operation()
        self.dfile1.close()
        self.dfile2.close()
        self.ofile.close()

    def input(self):
        if len(sys.argv) == 1:
            self.usage()

        try:
            opts, args = getopt.getopt(sys.argv[1:],"s:f:i:h:d:o:I:r:l:j:t:bw:a:", ["sb=", "freq=", "infile1=", "help=", "det",
                                                                                 "out","infile2=","deepx", "deepy", "jk", "tool", "beam", "scale","access"])
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

    def readMap(self, first_file = True, jackmode = "odde"):
        t = time.time()
        if first_file:
            dfile = self.dfile1
        else:
            dfile = self.dfile2

        freq_start = self.freq_list[0] - 1
        freq_end   = self.freq_list[-1] - 1
        sb_start = self.sb_list[0] - 1
        sb_end   = self.sb_list[-1] - 1

        if self.jack:
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
            else:
                map =  dfile["map"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                nhit =  dfile["nhit"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                rms =  dfile["rms"][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
        self.time += time.time() - t
        #print("Run time Read: ", time.time() - t)
        return map, nhit, rms

        """
        dfile   = h5.File(infile,'r')
        freq_start = self.freq_list[0] - 1
        freq_end   = self.freq_list[-1] - 1
        sb_start = self.sb_list[0] - 1
        sb_end   = self.sb_list[-1] - 1
        if self.jack:
            dname = "jackknives/" + mode + "_" + self.jk
            print(dname)
            if self.jk == "dayn" or self.jk == "half":                
                data0 = dfile[dname][0, self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                data1 = dfile[dname][1, self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                
            elif self.jk == "odde" or self.jk == "sdlb":
                data0 = dfile[dname][0, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                data1 = dfile[dname][1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
            
            dfile.close()
            return data0, data1
        else:
            dname = mode
            if self.beam:
                dname += "_beam"
                data =  dfile[dname][sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                dfile.close()
                return data
            else:
                data =  dfile[dname][self.det_list - 1, sb_start:sb_end + 1, freq_start:freq_end + 1, ...]
                dfile.close()
                return data
        """
                
    def writeMap(self, jackmode = "odde"):
        t = time.time()
        if self.jack:
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
        else:
            if self.beam:
                map_name    = "map_beam"
                nhit_name   = "nhit_beam"
                rms_name    = "rms_beam"
            else: 
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
        self.time += time.time() - t
        #print("Run time Write: ", time.time() - t, "\n")
        
        
    def operation(self):                
        if self.infile1 != None and self.infile2 != None:
            if self.jack:
                self.time = 0
                for jack in self.jk:
                    self.map1, self.nhit1, self.rms1 = self.readMap(True, jack)
                    self.map2, self.nhit2, self.rms2 = self.readMap(False, jack)
                    t = time.time()
                    
                    self.mask = self.nhit1 * self.nhit2 > 0
                    self.map1 = np.where(self.mask, self.map1, 0)
                    self.map2 = np.where(self.mask, self.map2, 0)
                    self.nhit1 = np.where(self.mask, self.nhit1, 0)
                    self.nhit2 = np.where(self.mask, self.nhit2, 0)
                    self.inv_rms1 = np.where(self.mask, 1 / self.rms1, 0)
                    self.inv_rms2 = np.where(self.mask, 1 / self.rms2, 0)
                    print("Mask time: ", time.time() - t)
                    
                    self.time += time.time() - t
                    
                    if self.tool == "coadd":
                        self.map = np.zeros_like(self.map1)
                        self._map = np.zeros_like(self.map1)
                        self._map[self.mask] = self.coadd(self.map1[self.mask], self.inv_rms1[self.mask],
                                                        self.map2[self.mask], self.inv_rms2[self.mask])
                        self.nhit = self.add(self.nhit1, self.nhit2)
                        self._nhit = self.add(self.nhit1, self.nhit2)
                        self.rms  = self.add_rms(self.inv_rms1, self.inv_rms2) 
                        self._rms  = self.add_rms(self.inv_rms1, self.inv_rms2) 
                        if len(self.map1.shape) == 5:   
                            float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
                            int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
                            self.maputilslib.coadd5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                                                float32_array5, int32_array5, float32_array5,
                                                                float32_array5, int32_array5, float32_array5,
                                                                ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                                                ctypes.c_int, ctypes.c_int]
                            n0, n1, n2, n3, n4 = self.map1.shape
                            self.map = np.zeros_like(self.map1, dtype = ctypes.c_float)
                            self.nhit = np.zeros_like(self.nhit1, dtype = ctypes.c_int)
                            self.rms = np.zeros_like(self.rms1, dtype = ctypes.c_float)
                            start = time.time()
                            self.maputilslib.coadd5D(self.map1, self.nhit1, self.rms1,
                                                     self.map2, self.nhit2, self.rms2, 
                                                     self.map, self.nhit, self.rms,
                                                     n0, n1, n2, n3, n4)
                            print("Run time outside C: ", time.time() - start, " sec")
                            print(np.allclose(np.nan_to_num(self._rms, posinf = 0, neginf=0), self.rms))
                            print(self.rms1[1, 3, 63, 79, 30], self._rms[1, 3, 63, 79, 30] - self.rms[1, 3, 63, 79, 30])
                
                    elif self.tool == "subtract": 
                        self.map   = self.subtract(self.map1, self.map2)
                        self.nhit   = self.subtract(self.nhit1, self.nhit2)
                        self.rms   = self.subtract_rms(self.rms1, self.rms2)
                    self.writeMap(jack)
                print("Run time jack-loop: ", self.time)   
            else:
                self.map1, self.nhit1, rms1 = self.readMap(True)
                self.map2, self.nhit2, rms2 = self.readMap(False)
                self.mask = self.nhit1 * self.nhit2 > 0
                
                self.map1 = np.where(self.mask, self.map1, 0)
                self.map2 = np.where(self.mask, self.map2, 0)
                self.nhit1 = np.where(self.mask, self.nhit1, 0)
                self.nhit2 = np.where(self.mask, self.nhit2, 0)
                self.inv_rms1 = np.where(self.mask, 1 / rms1, 0)
                self.inv_rms2 = np.where(self.mask, 1 / rms2, 0)
                
                if self.tool == "coadd":
                    self.map = np.zeros_like(self.map1)
                    self.map[self.mask] = self.coadd(self.map1[self.mask], self.inv_rms1[self.mask],
                                                      self.map2[self.mask], self.inv_rms2[self.mask])
                    self.nhit = self.add(self.nhit1, self.nhit2)
                    self.rms  = self.add_rms(self.inv_rms1, self.inv_rms2) 

                elif self.tool == "subtract": 
                    self.map   = self.subtract(self.map1, self.map2)
                    self.nhit   = self.subtract(self.nhit1, self.nhit2)
                    self.rms   = self.subtract_rms(self.rms1, self.rms2)
                self.writeMap()
                

        """
            self.data1 = self.readMap(self.infile1, mode = self.map_mode)
            self.data2 = self.readMap(self.infile2, mode = self.map_mode)
            self.hits1 = self.readMap(self.infile1, mode = "nhit")
            self.hits2 = self.readMap(self.infile2, mode = "nhit")
            self.mask = self.hits1 * self.hits2 > 0
            self.data1 = np.where(self.mask, self.data1, 0)
            self.data2 = np.where(self.mask, self.data2, 0)

            if self.map_mode == "map":
                if self.tool == "add":
                    self.data = self.add(self.data1, self.data2)

                elif self.tool == "coadd":
                    rms1 = self.readMap(self.infile1, mode = "rms")
                    rms2 = self.readMap(self.infile2, mode = "rms")
                    self.data = np.zeros_like(rms1)
                    self.data[self.mask] = self.coadd(self.data1[self.mask], 1 / rms1[self.mask],
                                                      self.data2[self.mask], 1 / rms2[self.mask])

                elif self.tool == "subtract": 
                    self.data   = self.subtract(self.data1, self.data2)

            elif self.map_mode == "nhit":
                if self.tool == "add":
                    self.data = self.add(self.hits1, self.hits2)

                elif self.tool == "subtract": 
                    self.data = self.subtract(self.hits1, self.hits2)

            elif self.map_mode == "rms":
                self.data = np.zeros_like(self.data1)
                if self.tool == "coadd":
                    self.data[self.mask] = self.add_rms(1 / self.data1[self.mask],
                                                        1 / self.data2[self.mask])
                elif self.tool == "subtract": 
                    self.data[self.mask] = self.subtract_rms(1 / self.data1[self.mask], 
                                                             1 / self.data2[self.mask])
        """

    def add(self, data1, data2):
        return data1 + data2
    
    def subtract(self, data1, data2):
        return data1 - data2

    def multiply(self, data, factor):
        return factor * data

    def coadd(self, map1, map2, inv_rms1, inv_rms2):
        t = time.time()
        inv_var1 = np.square(inv_rms1)
        inv_var2 = np.square(inv_rms2)
        var_inv = inv_var1 + inv_var2
        coadded = map1 * inv_rms1 + map2 * inv_rms2
        coadded /= var_inv
        self.time += time.time() - t
        #print("Run time Coadd: ", time.time() - t)
        return coadded

    def add_rms(self, inv_rms1, inv_rms2):
        inv_var1 = np.square(inv_rms1)
        inv_var2 = np.square(inv_rms2)
        sum = inv_var1 + inv_var2 
        sum = 1 / np.sqrt(sum)
        return sum

    def subtract_rms(self, inv_rms1, inv_rms2):
        var_inv1 = np.square(inv_rms1)
        var_inv2 = np.square(inv_rms2)
        sum = var_inv1 - var_inv2 
        sum = np.sqrt(sum)
        return sum
    
    def coaddfeed(self, data, inv_rms):
        inv_var = np.square(inv_rms)
        coadded = np.sum(map * inv_var, axis = 0)
        coadded /= np.sum(inv_var, axis = 0)
        return coadded

    def coaddfeed_rms(self, inv_rms):
        inv_var = np.square(inv_rms)
        coadded = np.sum(inv_var, axis = 0)
        return coadded


if __name__ == "__main__":
    t = time.time()
    map = Atlas()
    print("Run time: ", time.time() - t, " sec")




    
