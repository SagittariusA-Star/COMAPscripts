# ifndef MAPOBJECTS_H
# define MAPOBJECTS_H

# include <ctype.h>
# include <stdio.h>
# include <stdlib.h>
# include <unistd.h>
# include <string>
# include <vector>
# include <iostream>
# include <cstring>
# include <math.h>
# include "H5Cpp.h"
# include <armadillo>
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
        std::vector<int> sb_list;
        std::vector<int> freq_list;
        std::vector<int> feed_list;
        int sb;
        int freq;
        int feed;
        int list_input;

        int sim_numb = 0;
        int x_index;
        int y_index;
        int number_of_maps;
        int plot_len;

        char *infile;
        char *outfile;
        std::string plot;
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

        
        MapObjects(int _argc, char *_argv[]);
        ~MapObjects();
        void usage();
        void command_input();
        void read_map();
        

    private:
        int argc;
        char *argv[];
    


};



# endif