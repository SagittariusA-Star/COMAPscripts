import time
import sys
import getopt
import numpy as np
import h5py as h5
import textwrap
import ctypes

class Atlas:
    def __init__(self):
        self.jk_choices   = ["odde", "dayn", "half", "sdlb"]
        self.jk           = None
        self.jack         = False
        self.tool_choices   = ["coadd", "subtract", "dgradeXY", "dgradeZ", "dgradeXYZ",
                                                    "ugradeXY", "ugradeZ", "ugradeXYZ"]
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
        if self.infile1 != None and self.infile2 != None:
            if self.jack and ("jackknives" not in self.dfile1 or "jackknives" not in self.dfile2):
                print("One or both of the input files does not contain any jackknife information!")
                sys.exit()
        if not self.full and not self.beam and not self.jack:
            self.everything = True
            if self.infile1 != None and self.infile2 != None:
                if "jackknives" in self.dfile1 and "jackknives" in self.dfile2:
                    nhit_lst = [i for i in self.dfile1["jackknives"].keys() if "nhit_" in i]
                    self.jk =  [i.split("_")[1] for i in nhit_lst]
            else: 
                if "jackknives" in self.dfile1:
                    nhit_lst = [i for i in self.dfile1["jackknives"].keys() if "nhit_" in i]
                    self.jk =  [i.split("_")[1] for i in nhit_lst]
        self.operation()
        self.dfile1.close()
        if self.infile1 != None and self.infile2 != None:
            self.dfile2.close() 
        self.ofile.close()

    def input(self):
        if len(sys.argv) == 1:
            self.usage()

        try:
            opts, args = getopt.getopt(sys.argv[1:],"s:f:i:h:d:o:I:r:l:j:t:bw:a::F", ["sb=", "freq=", "infile1=", "help", "de=t",
                                                                                      "out=","infile2=","deepx", "deepy", "jk=", "tool=", "beam", "scale=","access=", "full"])
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
                conditionXY   = "dgradeXY" in arg.split(",") or "ugradeXY" in arg.split(",")
                conditionZ    = "dgradeZ" in arg.split(",") or "ugradeZ" in arg.split(",")
                conditionXYZ  = "dgradeXYZ" in arg.split(",") or "ugradeXYZ" in arg.split(",")
                
                if "coadd" in arg or "subtract" in arg:
                    if self.infile1 != None and self.infile2 == None:
                        print("To perform a coadd or subtraction two input files must be given!")
                        sys.exit()
                elif conditionXY and len(arg.split(",")) == 2:
                    if self.infile1 != None and self.infile2 != None:
                        print("Tool ugradeXY and dgradeXY are only supported for single input file!")
                        sys.exit()
                    self.tool, self.merge_numXY = arg.split(",")
                    self.merge_numXY = int(self.merge_numXY)
                    n_x = np.array(self.dfile1["n_x"])
                    n_y = np.array(self.dfile1["n_y"])
                    if n_x % self.merge_numXY != 0 or n_y % self.merge_numXY != 0: 
                        message = """\
                        Make sure that the pixel grid resolution of input map 
                        file is a multiple of the number of merging pixels!
                        """
                        print(textwrap.dedent(message))
                        sys.exit()
                elif conditionXY and len(arg.split(",")) != 2:
                    message = """\
                    To use ugradeXY or dgradeXY tool please provide a number of pixels to co-merge along each 
                    axis; e.g. -t dgrade,2 (don't forget the comma!!)!
                    """
                    print(textwrap.dedent(message))
                    sys.exit()
                
                elif conditionZ and len(arg.split(",")) == 2:
                    if self.infile1 != None and self.infile2 != None:
                        print("Tool ugradeZ and dgradeZ are only supported for single input file!")
                        sys.exit()
                    self.tool, self.merge_numZ = arg.split(",")
                    self.merge_numZ = int(self.merge_numZ)
                    n_z = self.dfile1["freq"].shape[1]
                    if n_z % self.merge_numZ != 0: 
                        message = """\
                        Make sure that the pixel grid resolution of input map file is a multiple
                        of the number of merging pixels!
                        """
                        print(textwrap.dedent(message))
                        sys.exit()
                elif conditionZ and len(arg.split(",")) != 2:
                    message = """\
                    To use ugradeZ or dgradeZ tool please provide a number of frequency channels to co-merge along 
                    each axis; e.g. -t dgrade,2 (don't forget the comma!!)!
                    """
                    print(textwrap.dedent(message))
                    sys.exit()
                
                elif conditionXYZ and len(arg.split(",")) == 3:
                    if self.infile1 != None and self.infile2 != None:
                        print("Tool ugradeXYZ and dgradeXYZ are only supported for single input file!")
                        sys.exit()
                    self.tool, self.merge_numXY, self.merge_numZ = arg.split(",")
                    self.merge_numXY, self.merge_numZ = int(self.merge_numXY), int(self.merge_numZ)
                    n_x = np.array(self.dfile1["n_x"])
                    n_y = np.array(self.dfile1["n_y"])
                    n_z = self.dfile1["freq"].shape[1]
                    conditionX = n_x % self.merge_numXY != 0
                    conditionY = n_y % self.merge_numXY != 0
                    conditionZ = n_z % self.merge_numZ  != 0
                    condition  = conditionX or conditionY or conditionZ
                    if condition: 
                        message = """\
                        Make sure that the pixel grid resolution and number of frequency channels of 
                        input map file is a multiple of the number of merging pixels!
                        """
                        print(textwrap.dedent(message))
                        sys.exit()
                elif conditionXYZ and len(arg.split(",")) != 3:
                    message = """
                    To use ugradeXYZ or dgradeXYZ tool please provide a number of pixels 
                    to co-merge along each axis as well as a number of frequency 
                    channels to co-merge; e.g. -t dgrade,2,2 (don't forget the commas!!)!
                    """
                    print(textwrap.dedent(message))
                    sys.exit()
                else:
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
        
    def writeMap(self, jackmode = None, write_the_rest = False, custom_data = None, custom_name = None):
        if custom_data != None and custom_name:
            if map_name in self.ofile and self.access == "a":
                """
                To overwrite existing dataset with different shape, the existing
                dataset must first be deleted.
                """   
                del self.ofile[map_name]
                self.ofile.create_dataset(custom_name, data = custom_data)
            else:
                self.ofile.create_dataset(custom_name, data = custom_data)

        elif jackmode != None:
            map_name    = "jackknives/map_" + jackmode
            nhit_name   = "jackknives/nhit_" + jackmode
            rms_name    = "jackknives/rms_" + jackmode
            if map_name in self.ofile and self.access == "a":
                """
                To overwrite existing dataset with different shape, the existing
                dataset must first be deleted.
                """   
                del self.ofile[map_name]
                del self.ofile[nhit_name]
                del self.ofile[rms_name]
                self.ofile.create_dataset(map_name, data = self.map)
                self.ofile.create_dataset(nhit_name, data = self.nhit)
                self.ofile.create_dataset(rms_name, data = self.rms)
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
                """
                To overwrite existing dataset with different shape, the existing
                dataset must first be deleted.
                """
                del self.ofile[map_name]
                del self.ofile[nhit_name]
                del self.ofile[rms_name]
                self.ofile.create_dataset(map_name, data = self.map)
                self.ofile.create_dataset(nhit_name, data = self.nhit)
                self.ofile.create_dataset(rms_name, data = self.rms)

            else:
                self.ofile.create_dataset(map_name, data = self.map)
                self.ofile.create_dataset(nhit_name, data = self.nhit)
                self.ofile.create_dataset(rms_name, data = self.rms)
        
        if write_the_rest:
            if self.infile1 != None and self.infile2 != None:
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
            else:   
                if self.tool == "dgradeXY" or self.tool == "ugradeXY":
                    self.merge_numZ = 1
                elif self.tool == "dgradeZ" or self.tool == "ugradeZ":
                    self.merge_numXY = 1
                
                condition1 = "x" in self.ofile and "y" in self.ofile 
                condition2 = "n_x" in self.ofile and "n_y" in self.ofile 
                condition3 = "nside" in self.ofile
                condition4 = "freq" in self.ofile
                condition  = condition1 and condition2 and condition3 and condition4

                if not condition and "dgrade" in self.tool:
                    x1, y1 = self.dfile1["x"][:], self.dfile1["y"][:]
                    x      = x1.reshape(int(len(x1) / self.merge_numXY), self.merge_numXY) 
                    y      = y1.reshape(int(len(y1) / self.merge_numXY), self.merge_numXY)
                    x      = np.mean(x, axis = 1)
                    y      = np.mean(y, axis = 1)
                    nside  = np.array(self.dfile1["nside"]) / self.merge_numXY                
                    freq   = self.dfile1["freq"][:]
                    freq   = freq.reshape(freq.shape[0], int(freq.shape[1] / self.merge_numZ), self.merge_numZ)
                    freq   = np.mean(freq, axis = 2)

                    self.ofile.create_dataset("x",      data = x)
                    self.ofile.create_dataset("y",      data = y)
                    self.ofile.create_dataset("n_x",    data = len(x))
                    self.ofile.create_dataset("n_y",    data = len(y))
                    self.ofile.create_dataset("nside",  data = nside)
                    self.ofile.create_dataset("freq",   data = freq)

                elif not condition and "ugrade" in self.tool:
                    x1, y1 = self.dfile1["x"][:], self.dfile1["y"][:]
                    x      = np.linspace(np.min(x1), np.max(x1), len(x1) * self.merge_numXY) 
                    y      = np.linspace(np.min(y1), np.max(y1), len(y1) * self.merge_numXY)
                    nside  = np.array(self.dfile1["nside"]) / self.merge_numXY                
                    freq1   = self.dfile1["freq"][:]
                    freq    = np.zeros((freq1.shape[0], freq1.shape[1] * self.merge_numZ))
                    for i in range(freq.shape[0]):
                        freq[i, :] = np.linspace(np.min(freq1[i, :]), 
                                                 np.max(freq1[i, :]), 
                                                 freq1.shape[1] * self.merge_numZ)
                
                    self.ofile.create_dataset("x",      data = x)
                    self.ofile.create_dataset("y",      data = y)
                    self.ofile.create_dataset("n_x",    data = len(x))
                    self.ofile.create_dataset("n_y",    data = len(y))
                    self.ofile.create_dataset("nside",  data = nside)
                    self.ofile.create_dataset("freq",   data = freq)

                data_not_to_copy = ["jackknives", "map", "map_beam", "nhit", 
                                    "nhit_beam", "rms", "rms_beam",
                                    "x", "y", "n_x", "n_y", "nside", "freq"]
                jk_data_not_to_copy = ["map_dayn",  "map_half",  "map_odde",  "map_sdlb",
                                       "nhit_dayn", "nhit_half", "nhit_odde", "nhit_sdlb",
                                       "rms_dayn",  "rms_half",  "rms_odde",  "rms_sdlb"]
                for name in self.dfile1.keys():
                    if name not in self.ofile.keys() and name not in data_not_to_copy:
                        self.ofile.create_dataset(name, data = self.dfile1[name])    
                
                    if "jackknives" in self.dfile1 and "jackknives" in self.ofile:
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
                            if len(self.map1.shape) == 6:
                                self.C_subtract6D(self.map1, self.nhit1, self.rms1,
                                                  self.map2, self.nhit2, self.rms2)  
                            
                            if len(self.map1.shape) == 5:
                                self.C_subtract5D(self.map1, self.nhit1, self.rms1,
                                                  self.map2, self.nhit2, self.rms2) 
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
                    self.C_subtract4D(self.map1, self.nhit1, self.rms1,
                                      self.map2, self.nhit2, self.rms2)  

                self.writeMap()
                self.beam = False
            
            if self.jack:
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

            if self.full:
                _beam = self.beam
                self.beam = False
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                self.map2, self.nhit2, self.rms2 = self.readMap(False)
                
                if self.tool == "coadd":
                    self.C_coadd5D(self.map1, self.nhit1, self.rms1,
                                    self.map2, self.nhit2, self.rms2)
    
                elif self.tool == "subtract": 
                    self.C_subtract5D(self.map1, self.nhit1, self.rms1,
                                      self.map2, self.nhit2, self.rms2)  
                        
                self.writeMap()
                self.beam = _beam
        
            if self.beam:
                _full = self.full
                self.full = False
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                self.map2, self.nhit2, self.rms2 = self.readMap(False)
                
                if self.tool == "coadd":
                    self.C_coadd4D(self.map1, self.nhit1, self.rms1,
                                    self.map2, self.nhit2, self.rms2)

                elif self.tool == "subtract": 
                    self.C_subtract4D(self.map1, self.nhit1, self.rms1,
                                      self.map2, self.nhit2, self.rms2)  
                        
                self.writeMap()
                self.full = _full
            self.writeMap(write_the_rest = True)

        if self.infile1 != None and self.infile2 == None:
            if self.everything:
                if "jackknives" in self.dfile1:
                    for jack in self.jk:
                        self.map1, self.nhit1, self.rms1 = self.readMap(True, jack)
                        if self.tool == "dgradeXY":
                            if len(self.map1.shape) == 6:
                                self.C_dgradeXY6D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)

                            elif len(self.map1.shape) == 5:
                                self.C_dgradeXY5D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)

                        elif self.tool == "dgradeZ":
                            if len(self.map1.shape) == 6:
                                self.C_dgradeZ6D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)
                            
                            elif len(self.map1.shape) == 5:
                                self.C_dgradeZ5D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)

                        elif self.tool == "dgradeXYZ":
                            if len(self.map1.shape) == 6:
                                self.C_dgradeXYZ6D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)
                            
                            elif len(self.map1.shape) == 5:
                                self.C_dgradeXYZ5D(self.map1, self.nhit1, self.rms1)
                                self.writeMap(jack)
                        
                        elif self.tool == "ugradeXY":
                            if len(self.map1.shape) == 6:
                                self.C_ugradeXY6D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                                self.C_ugradeXY6D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                                self.C_dgradeXY6D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                            
                            elif len(self.map1.shape) == 5:
                                self.C_ugradeXY5D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                                self.C_ugradeXY5D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                                self.C_dgradeXY5D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                            
                        elif self.tool == "ugradeZ":
                            if len(self.map1.shape) == 6:
                                self.C_ugradeZ6D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                                self.C_ugradeZ6D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                                self.C_ugradeZ6D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                            
                            elif len(self.map1.shape) == 5:
                                self.C_ugradeZ5D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                                self.C_ugradeZ5D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                                self.C_ugradeZ5D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)

                        elif self.tool == "ugradeXYZ":
                            if len(self.map1.shape) == 6:
                                self.C_ugradeXYZ6D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)
                                
                                self.C_ugradeXYZ6D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)
                                
                                self.C_ugradeXYZ6D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                            
                            elif len(self.map1.shape) == 5:
                                self.C_ugradeXYZ5D_float(self.map1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)
                                
                                self.C_ugradeXYZ5D_int(self.nhit1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)
                                
                                self.C_ugradeXYZ5D_float(self.rms1)
                                self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)

                self.full = True
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                
                if self.tool == "dgradeXY":
                    self.C_dgradeXY5D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()
                
                elif self.tool == "dgradeZ":
                    self.C_dgradeZ5D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()
                
                elif self.tool == "dgradeXYZ":
                    self.C_dgradeXYZ5D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()
                
                elif self.tool == "ugradeXY":
                    self.C_ugradeXY5D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")

                    self.C_ugradeXY5D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")

                    self.C_ugradeXY5D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")
                
                elif self.tool == "ugradeZ":
                    self.C_ugradeZ5D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")
                    
                    self.C_ugradeZ5D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")
                    
                    self.C_ugradeZ5D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")
                
                elif self.tool == "ugradeXYZ":
                    self.C_ugradeXYZ5D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")

                    self.C_ugradeXYZ5D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")

                    self.C_ugradeXYZ5D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")
                
                
                self.full = False
                self.beam = True
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                
                if self.tool == "dgradeXY":
                    self.C_dgradeXY4D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()

                elif self.tool == "dgradeZ":
                    self.C_dgradeZ4D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()

                elif self.tool == "dgradeXYZ":
                    self.C_dgradeXYZ4D(self.map1, self.nhit1, self.rms1)
                    self.writeMap()

                elif self.tool == "ugradeXY":
                    self.C_ugradeXY4D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")
                    
                    self.C_ugradeXY4D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")
                    
                    self.C_ugradeXY4D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")

                elif self.tool == "ugradeZ":
                    self.C_ugradeZ4D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")

                    self.C_ugradeZ4D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")

                    self.C_ugradeZ4D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")

                elif self.tool == "ugradeXYZ":
                    self.C_ugradeXYZ4D_float(self.map1)
                    self.writeMap(custom_data = self.map, custom_name = "map_")

                    self.C_ugradeXYZ4D_int(self.nhit1)
                    self.writeMap(custom_data = self.map, custom_name = "nhit_")

                    self.C_ugradeXYZ4D_float(self.rms1)
                    self.writeMap(custom_data = self.map, custom_name = "rms_")

                self.writeMap()
                self.beam = False
            if self.jack:
                for jack in self.jk:
                    self.map1, self.nhit1, self.rms1 = self.readMap(True, jack)
                    if self.tool == "dgradeXY":
                        if len(self.map1.shape) == 6:
                            self.C_dgradeXY6D(self.map1, self.nhit1, self.rms1)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_dgradeXY5D(self.map1, self.nhit1, self.rms1)
                        
                    elif self.tool == "dgradeZ":
                        if len(self.map1.shape) == 6:
                            self.C_dgradeZ6D(self.map1, self.nhit1, self.rms1)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_dgradeZ5D(self.map1, self.nhit1, self.rms1)

                    elif self.tool == "dgradeXYZ":
                        if len(self.map1.shape) == 6:
                            self.C_dgradeXYZ6D(self.map1, self.nhit1, self.rms1)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_dgradeXYZ5D(self.map1, self.nhit1, self.rms1)

                    elif self.tool == "ugradeXY":
                        if len(self.map1.shape) == 6:
                            self.C_ugradeXY6D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                            self.C_ugradeXY6D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                            self.C_dgradeXY6D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_ugradeXY5D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                            self.C_ugradeXY5D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                            self.C_dgradeXY5D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                        
                    elif self.tool == "ugradeZ":
                        if len(self.map1.shape) == 6:
                            self.C_ugradeZ6D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                            self.C_ugradeZ6D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                            self.C_ugradeZ6D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_ugradeZ5D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)

                            self.C_ugradeZ5D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)

                            self.C_ugradeZ5D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)

                    elif self.tool == "ugradeXYZ":
                        if len(self.map1.shape) == 6:
                            self.C_ugradeXYZ6D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)
                            
                            self.C_ugradeXYZ6D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)
                            
                            self.C_ugradeXYZ6D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)
                        
                        elif len(self.map1.shape) == 5:
                            self.C_ugradeXYZ5D_float(self.map1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/map_" + jack)
                            
                            self.C_ugradeXYZ5D_int(self.nhit1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/nhit_" + jack)
                            
                            self.C_ugradeXYZ5D_float(self.rms1)
                            self.writeMap(custom_data = self.map, custom_name = "jackknives/rms_" + jack)

                    self.writeMap(jack)

            if self.full:
                _beam = self.beam
                self.beam = False
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                
                if self.tool == "dgradeXY":
                    self.C_dgradeXY5D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "dgradeZ":
                    self.C_dgradeZ5D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "dgradeXYZ":
                    self.C_dgradeXYZ5D(self.map1, self.nhit1, self.rms1)

                elif self.tool == "ugradeXY":
                    self.C_ugradeXY5D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "ugradeZ":
                    self.C_ugradeZ5D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "ugradeXYZ":
                    self.C_ugradeXYZ5D(self.map1, self.nhit1, self.rms1)
            
                self.writeMap()
                self.beam = _beam
        
            if self.beam:
                _full = self.full
                self.full = False
                self.map1, self.nhit1, self.rms1 = self.readMap(True)
                
                if self.tool == "dgradeXY":
                    self.C_dgradeXY4D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "dgradeZ":
                    self.C_dgradeZ4D(self.map1, self.nhit1, self.rms1)

                elif self.tool == "dgradeXYZ":
                    self.C_dgradeXYZ4D(self.map1, self.nhit1, self.rms1)

                elif self.tool == "ugradeXY":
                    self.C_ugradeXY4D(self.map1, self.nhit1, self.rms1)
                
                elif self.tool == "ugradeZ":
                    self.C_ugradeZ4D(self.map1, self.nhit1, self.rms1)

                elif self.tool == "ugradeXYZ":
                    self.C_ugradeXYZ4D(self.map1, self.nhit1, self.rms1)

                self.writeMap()
                self.full = _full
            self.writeMap(write_the_rest = True)            

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

    def C_dgradeXY4D(self, map_h, nhit_h, rms_h):
            float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
            int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
            self.maputilslib.dgradeXY4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                                  float32_array4, int32_array4, float32_array4,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                  ctypes.c_int]
            n0, n1, n2, n3 = map_h.shape
            N2, N3 = int(n2 / self.merge_numXY), int(n3 / self.merge_numXY)
            
            self.map = np.zeros( (n0, n1, N2, N3), dtype = ctypes.c_float)
            self.nhit = np.zeros((n0, n1, N2, N3), dtype = ctypes.c_int)
            self.rms = np.zeros( (n0, n1, N2, N3), dtype = ctypes.c_float)

            self.maputilslib.dgradeXY4D(map_h,    nhit_h,     rms_h,
                                    self.map,   self.nhit,  self.rms,
                                    n0,         n1,         n2,
                                    n3,         N2,         N3,
                                    self.merge_numXY)

    def C_dgradeXY5D(self, map_h, nhit_h, rms_h):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.dgradeXY5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                              float32_array5, int32_array5, float32_array5,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3, n4 = map_h.shape
        N3, N4 = int(n3 / self.merge_numXY), int(n4 / self.merge_numXY)
        
        self.map = np.zeros( (n0, n1, n2, N3, N4), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, n2, N3, N4), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, n2, N3, N4), dtype = ctypes.c_float)

        self.maputilslib.dgradeXY5D(map_h,    nhit_h,     rms_h,
                                  self.map, self.nhit,  self.rms,
                                  n0,       n1,         n2,
                                  n3,       n4,         N3,
                                  N4,       self.merge_numXY)

    def C_dgradeXY6D(self, map_h, nhit_h, rms_h):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.dgradeXY6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                              float32_array6, int32_array6, float32_array6,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_h.shape
        N4, N5 = int(n4 / self.merge_numXY), int(n5 / self.merge_numXY)
        
        self.map = np.zeros( (n0, n1, n2, n3, N4, N5), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, n2, n3, N4, N5), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, n2, n3, N4, N5), dtype = ctypes.c_float)

        self.maputilslib.dgradeXY6D(map_h,        nhit_h,         rms_h,
                                  self.map,     self.nhit,      self.rms,
                                  n0,           n1,             n2,
                                  n3,           n4,             n5, 
                                  N4,           N5,             self.merge_numXY)

    def C_dgradeZ4D(self, map_h, nhit_h, rms_h):
            float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
            int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
            self.maputilslib.dgradeZ4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                                  float32_array4, int32_array4, float32_array4,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int]
            n0, n1, n2, n3  = map_h.shape
            N1              = int(n1 / self.merge_numZ)
            
            self.map = np.zeros( (n0, N1, n2, n3), dtype = ctypes.c_float)
            self.nhit = np.zeros((n0, N1, n2, n3), dtype = ctypes.c_int)
            self.rms = np.zeros( (n0, N1, n2, n3), dtype = ctypes.c_float)

            self.maputilslib.dgradeZ4D(map_h,    nhit_h,     rms_h,
                                       self.map,   self.nhit,  self.rms,
                                       n0,         n1,         n2,
                                       n3,         N1,         self.merge_numZ)

    def C_dgradeZ5D(self, map_h, nhit_h, rms_h):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.dgradeZ5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                              float32_array5, int32_array5, float32_array5,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int]
        n0, n1, n2, n3, n4 = map_h.shape
        N2 = int(n2 / self.merge_numZ)
        
        self.map = np.zeros( (n0, n1, N2, n3, n4), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, N2, n3, n4), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, N2, n3, n4), dtype = ctypes.c_float)

        self.maputilslib.dgradeZ5D(map_h,    nhit_h,     rms_h,
                                  self.map, self.nhit,  self.rms,
                                  n0,       n1,         n2,
                                  n3,       n4,         N2,
                                  self.merge_numZ)

    def C_dgradeZ6D(self, map_h, nhit_h, rms_h):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.dgradeZ6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                              float32_array6, int32_array6, float32_array6,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_h.shape
        N3 = int(n3 / self.merge_numZ)
        
        self.map = np.zeros( (n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, n2, N3, n4, n5), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)

        self.maputilslib.dgradeZ6D(map_h,        nhit_h,         rms_h,
                                  self.map,     self.nhit,      self.rms,
                                  n0,           n1,             n2,
                                  n3,           n4,             n5, 
                                  N3,           self.merge_numZ)
    
    def C_dgradeXYZ4D(self, map_h, nhit_h, rms_h):
            float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
            int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
            self.maputilslib.dgradeXYZ4D.argtypes = [float32_array4, int32_array4, float32_array4,
                                                  float32_array4, int32_array4, float32_array4,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                                  ctypes.c_int,   ctypes.c_int, ctypes.c_int]
            n0, n1, n2, n3 = map_h.shape
            N1, N2, N3 = int(n1 / self.merge_numZ), int(n2 / self.merge_numXY), int(n3 / self.merge_numXY)
            
            self.map = np.zeros( (n0, N1, N2, N3), dtype = ctypes.c_float)
            self.nhit = np.zeros((n0, N1, N2, N3), dtype = ctypes.c_int)
            self.rms = np.zeros( (n0, N1, N2, N3), dtype = ctypes.c_float)

            self.maputilslib.dgradeXYZ4D(map_h,    nhit_h,          rms_h,
                                         self.map, self.nhit,       self.rms,
                                         n0,       n1,              n2,
                                         n3,       N1,              N2,
                                         N3,       self.merge_numZ,  self.merge_numXY)

    def C_dgradeXYZ5D(self, map_h, nhit_h, rms_h):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.dgradeXYZ5D.argtypes = [float32_array5, int32_array5, float32_array5,
                                              float32_array5, int32_array5, float32_array5,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int]
        n0, n1, n2, n3, n4 = map_h.shape
        N2, N3, N4 = int(n2 / self.merge_numZ), int(n3 / self.merge_numXY), int(n4 / self.merge_numXY)
        
        self.map = np.zeros( (n0, n1, N2, N3, N4), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, N2, N3, N4), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, N2, N3, N4), dtype = ctypes.c_float)

        self.maputilslib.dgradeXYZ5D(map_h,    nhit_h,     rms_h,
                                  self.map, self.nhit,  self.rms,
                                  n0,       n1,         n2,
                                  n3,       n4,         N2,
                                  N3,       N4,         self.merge_numZ,
                                  self.merge_numXY)

    def C_dgradeXYZ6D(self, map_h, nhit_h, rms_h):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.dgradeXYZ6D.argtypes = [float32_array6, int32_array6, float32_array6,
                                              float32_array6, int32_array6, float32_array6,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int, ctypes.c_int,
                                              ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_h.shape
        N3, N4, N5 = int(n3 / self.merge_numZ), int(n4 / self.merge_numXY), int(n5 / self.merge_numXY)
        
        self.map = np.zeros( (n0, n1, n2, N3, N4, N5), dtype = ctypes.c_float)
        self.nhit = np.zeros((n0, n1, n2, N3, N4, N5), dtype = ctypes.c_int)
        self.rms = np.zeros( (n0, n1, n2, N3, N4, N5), dtype = ctypes.c_float)

        self.maputilslib.dgradeXYZ6D(map_h,         nhit_h,         rms_h,
                                  self.map,         self.nhit,      self.rms,
                                  n0,               n1,             n2,
                                  n3,               n4,             n5, 
                                  N3,               N4,             N5,
                                  self.merge_numZ,  self.merge_numXY)



    def C_ugradeXY4D_float(self, map_l):
        float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
        self.maputilslib.ugradeXY4D_float.argtypes = [float32_array4, float32_array4, ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int,   ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3 = map_l.shape
        N2, N3 = n2 * self.merge_numXY, n3 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, N2, N3), dtype = ctypes.c_float)

        self.maputilslib.ugradeXY4D_float(self.map,   map_l,      n0,         
                                    n1,         n2,         n3,         
                                    N2,         N3,         self.merge_numXY)

    def C_ugradeXY5D_float(self, map_l):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        self.maputilslib.ugradeXY5D_float.argtypes = [float32_array5,   float32_array5, ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int,   ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int,   ctypes.c_int,   
                                                        ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N3, N4 = n3 * self.merge_numXY, n4 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, N3, N4), dtype = ctypes.c_float)

        self.maputilslib.ugradeXY5D_float(self.map, map_l,        n0,       
                                        n1,         n2,         n3,       
                                        n4,         N3,         N4,       
                                        self.merge_numXY)

    def C_ugradeXY6D_float(self, map_l):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        self.maputilslib.ugradeXY6D_float.argtypes = [float32_array6,  float32_array6,    ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int,     ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int,     ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N4, N5 = n4 * self.merge_numXY, n5 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, n3, N4, N5), dtype = ctypes.c_float)

        self.maputilslib.ugradeXY6D_float(self.map,     map_l,      n0,           
                                        n1,             n2,         n3,           
                                        n4,             n5,         N4,           
                                        N5,             self.merge_numXY)

    def C_ugradeZ4D_float(self, map_l):
        float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
        self.maputilslib.ugradeZ4D_float.argtypes = [float32_array4, float32_array4,    ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int,       ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3  = map_l.shape
        N1              = n1 * self.merge_numZ
        
        self.map = np.zeros( (n0, N1, n2, n3), dtype = ctypes.c_float)

        self.maputilslib.ugradeZ4D_float(self.map,   map_l,     n0,         
                                        n1,          n2,        n3,         
                                        N1,         self.merge_numZ)

    def C_ugradeZ5D_float(self, map_l):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        self.maputilslib.ugradeZ5D_float.argtypes = [float32_array5, float32_array5,  ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N2 = n2 * self.merge_numZ
        
        self.map = np.zeros( (n0, n1, N2, n3, n4), dtype = ctypes.c_float)

        self.maputilslib.ugradeZ5D_float(self.map,  map_l,  n0,       
                                        n1,         n2,     n3,       
                                        n4,         N2,     self.merge_numZ)

    def C_ugradeZ6D_float(self, map_l):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        self.maputilslib.ugradeZ6D_float.argtypes = [float32_array6, float32_array6,  ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int,     ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int,     ctypes.c_int,   
                                                    ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N3 = n3 * self.merge_numZ
        
        self.map = np.zeros( (n0, n1, n2, N3, n4, n5), dtype = ctypes.c_float)

        self.maputilslib.ugradeZ6D_float(self.map,     map_l,     n0,           
                                        n1,             n2,     n3,           
                                        n4,             n5,     N3,           
                                        self.merge_numZ)
    
    def C_ugradeXYZ4D_float(self, map_l):
        float32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=4, flags="contiguous")
        self.maputilslib.ugradeXYZ4D_float.argtypes = [float32_array4, float32_array4,    ctypes.c_int,   
                                                ctypes.c_int,   ctypes.c_int,       ctypes.c_int,   
                                                ctypes.c_int,   ctypes.c_int,       ctypes.c_int,   
                                                ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3 = map_l.shape
        N1, N2, N3 = n1 * self.merge_numZ, n2 * self.merge_numXY, n3 * self.merge_numXY
        
        self.map = np.zeros( (n0, N1, N2, N3), dtype = ctypes.c_float)

        self.maputilslib.ugradeXYZ4D_float(self.map,    map_l,      n0,       
                                            n1,         n2,         n3,       
                                            N1,         N2,         N3,       
                                            self.merge_numZ,  self.merge_numXY)

    def C_ugradeXYZ5D_float(self, map_l):
        float32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=5, flags="contiguous")
        self.maputilslib.ugradeXYZ5D_float.argtypes = [float32_array5, float32_array5,    ctypes.c_int,   
                                                        ctypes.c_int,   ctypes.c_int,       ctypes.c_int,   
                                                        ctypes.c_int, ctypes.c_int,         ctypes.c_int,   
                                                        ctypes.c_int, ctypes.c_int,         ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N2, N3, N4 = n2 * self.merge_numZ, n3 * self.merge_numXY, n4 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, N2, N3, N4), dtype = ctypes.c_float)

        self.maputilslib.ugradeXYZ5D_float(self.map, map_l,                 n0,       
                                            n1,         n2,                 n3,       
                                            n4,         N2,                 N3,       
                                            N4,         self.merge_numZ,    self.merge_numXY)

    def C_ugradeXYZ6D_float(self, map_l):
        float32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_float, ndim=6, flags="contiguous")
        self.maputilslib.ugradeXYZ6D_float.argtypes = [float32_array6, float32_array6,    ctypes.c_int,   
                                                        ctypes.c_int, ctypes.c_int,       ctypes.c_int,   
                                                        ctypes.c_int, ctypes.c_int,       ctypes.c_int,   
                                                        ctypes.c_int, ctypes.c_int,       ctypes.c_int,   
                                                        ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N3, N4, N5 = n3 * self.merge_numZ, n4 * self.merge_numXY, n5 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, N3, N4, N5), dtype = ctypes.c_float)

        self.maputilslib.ugradeXYZ6D_float(self.map,         map_l,     n0,               
                                            n1,             n2,         n3,               
                                            n4,             n5,         N3,               
                                            N4,             N5,         self.merge_numZ,  
                                            self.merge_numXY)



    def C_ugradeXY4D_int(self, map_l):
        int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
        self.maputilslib.ugradeXY4D_int.argtypes = [int32_array4, int32_array4, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3 = map_l.shape
        N2, N3 = n2 * self.merge_numXY, n3 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, N2, N3), dtype = ctypes.c_int)

        self.maputilslib.ugradeXY4D_int(self.map,   map_l,      n0,         
                                    n1,         n2,         n3,         
                                    N2,         N3,         self.merge_numXY)

    def C_ugradeXY5D_int(self, map_l):
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.ugradeXY5D_int.argtypes = [int32_array5, int32_array5, ctypes.c_int,
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                    ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N3, N4 = n3 * self.merge_numXY, n4 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, N3, N4), dtype = ctypes.c_int)
        
        self.maputilslib.ugradeXY5D_int(self.map,   map_l,      n0,       
                                        n1,         n2,         n3,       
                                        n4,         N3,         N4,       
                                        self.merge_numXY)

    def C_ugradeXY6D_int(self, map_l):
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.ugradeXY6D_int.argtypes = [int32_array6, int32_array6,     ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int,     ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int,     ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N4, N5 = n4 * self.merge_numXY, n5 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, n3, N4, N5), dtype = ctypes.c_int)
        
        self.maputilslib.ugradeXY6D_int(self.map,     map_l,    n0,           
                                        n1,           n2,       n3,           
                                        n4,           n5,       N4,           
                                        N5,           self.merge_numXY)

    def C_ugradeZ4D_int(self, map_l):
        int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
        self.maputilslib.ugradeZ4D_int.argtypes = [int32_array4,    int32_array4,   ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int,   ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int]
        n0, n1, n2, n3  = map_l.shape
        N1              = n1 * self.merge_numZ
        
        self.map = np.zeros( (n0, N1, n2, n3), dtype = ctypes.c_int)

        self.maputilslib.ugradeZ4D_int(self.map,   map_l,   n0,         
                                        n1,        n2,      n3,         
                                        N1,        self.merge_numZ)

    def C_ugradeZ5D_int(self, map_l):
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.ugradeZ5D_int.argtypes = [int32_array5, int32_array5, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N2 = n2 * self.merge_numZ
        
        self.map = np.zeros( (n0, n1, N2, n3, n4), dtype = ctypes.c_int)

        self.maputilslib.ugradeZ5D_int(self.map,    map_l,      n0,       
                                        n1,         n2,         n3,       
                                        n4,         N2,         self.merge_numZ)

    def C_ugradeZ6D_int(self, map_l):
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.ugradeZ6D_int.argtypes = [int32_array6, int32_array6, ctypes.c_int,   
                                                   ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                   ctypes.c_int, ctypes.c_int, ctypes.c_int,   
                                                   ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N3 = n3 * self.merge_numZ
        
        self.map = np.zeros( (n0, n1, n2, N3, n4, n5), dtype = ctypes.c_int)

        self.maputilslib.ugradeZ6D_int(self.map,     map_l,     n0,           
                                        n1,             n2,     n3,           
                                        n4,             n5,     N3,           
                                        self.merge_numZ)
    
    def C_ugradeXYZ4D_int(self, map_l):
        int32_array4 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=4, flags="contiguous")
        self.maputilslib.ugradeXYZ4D_int.argtypes = [int32_array4, int32_array4,    ctypes.c_int,   
                                                    ctypes.c_int,   ctypes.c_int,   ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int,     ctypes.c_int,   
                                                    ctypes.c_int, ctypes.c_int]
        n0, n1, n2, n3 = map_l.shape
        N1, N2, N3 = n1 * self.merge_numZ, n2 * self.merge_numXY, n3 * self.merge_numXY
        
        self.map = np.zeros( (n0, N1, N2, N3), dtype = ctypes.c_int)

        self.maputilslib.ugradeXYZ4D_int(self.map,          map_l,       n0,       
                                         n1,                n2,          n3,       
                                         N1,                N2,          N3,       
                                         self.merge_numZ,   self.merge_numXY)

    def C_ugradeXYZ5D_int(self, map_l):
        int32_array5 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=5, flags="contiguous")
        self.maputilslib.ugradeXYZ5D_int.argtypes = [int32_array5, int32_array5, ctypes.c_int,   
                                                ctypes.c_int, ctypes.c_int,  ctypes.c_int,      
                                                ctypes.c_int, ctypes.c_int,  ctypes.c_int,   
                                                ctypes.c_int, ctypes.c_int,  ctypes.c_int]
        n0, n1, n2, n3, n4 = map_l.shape
        N2, N3, N4 = n2 * self.merge_numZ, n3 * self.merge_numXY, n4 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, N2, N3, N4), dtype = ctypes.c_int)

        self.maputilslib.ugradeXYZ5D_int(self.map, map_l,               n0,       
                                     n1,         n2,                n3,       
                                     n4,         N2,                N3,       
                                     N4,         self.merge_numZ,   self.merge_numXY)

    def C_ugradeXYZ6D_int(self, map_l):
        int32_array6 = np.ctypeslib.ndpointer(dtype=ctypes.c_int, ndim=6, flags="contiguous")
        self.maputilslib.ugradeXYZ6D_int.argtypes = [int32_array6,  int32_array6,    ctypes.c_int,   
                                                ctypes.c_int,   ctypes.c_int,    ctypes.c_int,   
                                                ctypes.c_int,   ctypes.c_int,    ctypes.c_int,   
                                                ctypes.c_int,   ctypes.c_int,    ctypes.c_int,   
                                                ctypes.c_int]
        n0, n1, n2, n3, n4, n5 = map_l.shape
        N3, N4, N5 = n3 * self.merge_numZ, n4 * self.merge_numXY, n5 * self.merge_numXY
        
        self.map = np.zeros( (n0, n1, n2, N3, N4, N5), dtype = ctypes.c_int)

        self.maputilslib.ugradeXYZ6D_int(self.map,        map_l,           n0,               
                                     n1,              n2,              n3,               
                                     n4,              n5,              N3,               
                                     N4,              N5,              self.merge_numZ,  
                                     self.merge_numXY)


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




    
