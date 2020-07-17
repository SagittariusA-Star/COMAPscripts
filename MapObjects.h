# ifndef MAPOBJECTS_H
# define MAPOBJECTS_H

# include <ctype.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include <string>
# include <vector>
# include <iostream>
# include <cstring>
# include <math.h>
# include "H5Cpp.h"
# include <armadillo>
#include <algorithm>    
# include <time.h>  


using std::cout;
using std::endl;
using std::vector;
# ifndef H5_NO_NAMESPACE
    using namespace H5;
# endif

class MapObjects {
    public:
        int opt;
        
        std::vector<int> sb_list    = {1,2,3,4};

        std::vector<int> feed_list = {1,2,3,4,5,6,7,8,9,10,
                                          11,12,13,14,15,16,17,18,19};

        std::vector<int> freq_list = {1,2,3,4,5,6,7,8,9,10,
                                          11,12,13,14,15,16,17,18,19,20,
                                          21, 22,23,24,25,26,27,28,29,30,
                                          31,32,33,34,35,36,37,38,39,40,
                                          41,42,43,44,45,46,47,48,49,50,
                                          51,52,53,54,55,56,57,58,58,60,
                                          61,62,63,64};
        
        int sb;
        int freq;
        int feed;

        int sim_numb = 0;
        int x_index;
        int y_index;
        int number_of_maps;
        int plot_len;

        char *infile;
        char *outfile;
        std::string plot = "all";
        std::string a;

        bool jupiter = false;
        bool deepx   = false;
        bool deepy   = false;

        double scale = 1;
        double rms_lim = 200000;

        int n;              // Length of char* lists

        double *freq_arr;
        double *x_arr;
        double *y_arr;
        double *map_arr;
        double *rms_arr;
        int    *hit_arr;

        double *dummy_arr;

        int Nsb;   
        int Nfeed; 
        int Nfreq; 
        int Nx;    
        int Ny;    
        
        MapObjects();
        ~MapObjects();
        void usage();
        void command_input(int argc, char *argv[]);
        void read_map(int argc, char *argv[]);
        

    private:
        char *split_str;
        
    
};



# endif