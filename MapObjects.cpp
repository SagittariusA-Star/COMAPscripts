# include "MapObjects.h"

MapObjects::MapObjects(){
}

void MapObjects::usage(){
    printf("These are the supported operations:\n");
}

void MapObjects::command_input(int argc, char *argv[]){
    opterr = 0;

    while ((opt = getopt(argc, argv, "s:i:p:f:h:d:o:x:y:j:z:w:m:r:l:n:")) != -1)
        switch (opt)
        {
        case 's':
            sb_list.clear();
            split_str = strtok(optarg, "[],");
            while (split_str != NULL){
                sb = atoi(split_str);
                if (sb == 0){
                    printf("Please provide the number of the side band in 1-base, not 0-base!\n");
                    exit(EXIT_FAILURE);
                }
                else if (sb > 4){
                    printf("There are only four side bands!\n");
                    exit(EXIT_FAILURE);
                }
                sb_list.push_back(sb);
                split_str = strtok(NULL, ",");
            }                
            break;
        case 'i':
            infile = optarg;
            break;
        case 'p':
            plot = optarg;
            if (plot != "map" && plot != "rms" && 
                plot != "hit" && plot != "map/rms" &&
                plot != "feed" && plot != "var" &&
                plot != "all"){
                cout << "Make sure one of the valid map modes is chosen;";
                cout << " 'map', 'rms', 'hit',\n 'feed', 'var' or 'map/rms'." << endl;
                cout << "'all' will mean all the latter.\n" << endl;
                exit(EXIT_FAILURE);
            }
            break;
        case 'f':
            freq_list.clear();
            split_str = strtok(optarg, "[],");
            while (split_str != NULL){
                freq = atoi(split_str);
                if (freq == 0){
                    printf("Please provide the number of the frequency channels with 1-base, not 0-base!\n");
                    exit(EXIT_FAILURE);
                }
                else if (freq > 64){
                    printf("There are only 64 frequencies pr. side band!\n");
                    exit(EXIT_FAILURE);
                }
                freq_list.push_back(freq);
                split_str = strtok(NULL, ",");
            }   
            break;
        case 'h':
            usage();
            break;
        case 'd':
            feed_list.clear();
            split_str = strtok(optarg, "[],");
            while (split_str != NULL){
                feed = atoi(split_str);
                if (feed == 0){
                    printf("Please provide the number of the feed/detector with 1-base, not 0-base!\n");
                    exit(EXIT_FAILURE);
                }
                else if (feed > 20){
                    printf("There are only 20 feeds/detectors!\n");
                    exit(EXIT_FAILURE);
                }
                feed_list.push_back(feed);
                split_str = strtok(NULL, ",");
            }   
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'w':
            scale = atof(optarg);
            break;
        case 'm':
            rms_lim = atof(optarg);
            break;
        case 'l':
            deepx = true;
            y_index = atoi(optarg) - 1;
            break;
        case 'r':
            deepy = true;
            x_index = atoi(optarg) - 1;
            break;
        case 'n':
            number_of_maps = atoi(optarg);
            if (number_of_maps != 1 && number_of_maps != 2){
                printf("Only operations on one or two input maps are supported!\n");
                exit(EXIT_FAILURE);
            }
            break;
        case '?':
            printf("Unknown option: %c\n", optopt);
        default:
            abort();
        }
}

void MapObjects::read_map(int argc, char *argv[]){
    command_input(argc, argv);
    const H5std_string FILE_NAME(infile);
    const H5std_string DATASET_x("x");
    const H5std_string DATASET_y("y");

    H5File file( FILE_NAME, H5F_ACC_RDONLY);
    DataSet data_x = file.openDataSet( DATASET_x );
    DataSet data_y = file.openDataSet( DATASET_y );
    
    DataSpace dataspace_x = data_x.getSpace();
    DataSpace dataspace_y = data_y.getSpace();
    
    int RANK_x  = dataspace_x.getSimpleExtentNdims();
    int RANK_y  = dataspace_y.getSimpleExtentNdims();
    
    hsize_t dims_x[RANK_x];
    hsize_t dims_y[RANK_y];
    
    dataspace_x.getSimpleExtentDims(dims_x);
    dataspace_y.getSimpleExtentDims(dims_y);
    
    DataSpace mspace_x(RANK_x, dims_x);
    DataSpace mspace_y(RANK_y, dims_y);
    
    int Nx     = dims_x[0];
    int Ny     = dims_y[0];
    
    x_arr    = new double[dataspace_x.getSimpleExtentNpoints()];
    y_arr    = new double[dataspace_y.getSimpleExtentNpoints()];
    
    data_x.read(x_arr, PredType::NATIVE_DOUBLE, mspace_x, dataspace_x );
    data_y.read(y_arr, PredType::NATIVE_DOUBLE, mspace_y, dataspace_y );
    
    if (feed_list.size() == 19 && plot != "feed"){
        Nsb = sb_list.size();
        Nfreq = freq_list.size();
        /*
        const H5std_string DATASET_freq("freq");
        DataSet data_freq = file.openDataSet( DATASET_freq);
        DataSpace dataspace_freq = data_freq.getSpace();
        int RANK_freq  = dataspace_freq.getSimpleExtentNdims();
        hsize_t dims_freq[RANK_freq];
        dataspace_freq.getSimpleExtentDims(dims_freq);
        DataSpace mspace_freq(RANK_freq, dims_freq);
        
        hsize_t start[RANK_freq] = {(unsigned long long int) (sb_list[0] - 1), (unsigned long long int	) (freq_list[0] - 1)};
        hsize_t count[RANK_freq] = {(unsigned long long int) Nsb,
                                    (unsigned long long int) Nfreq};
        
        hsize_t start_out[RANK_freq] = {0, 0};
        hsize_t count_out[RANK_freq] = {(unsigned long long int) Nsb,
                                        (unsigned long long int) Nfreq};

        dataspace_freq.selectHyperslab(H5S_SELECT_SET, count, start);
        mspace_freq.selectHyperslab(H5S_SELECT_SET, count_out, start_out);
        hssize_t Npoints = mspace_freq.getSimpleExtentNpoints(); 
        freq_arr  = new double[Npoints];
        for (int i = 0; i < Npoints; i++){
            freq_arr[i] = 0;
        }
        data_freq.read(freq_arr, PredType::NATIVE_DOUBLE, mspace_freq, dataspace_freq); 
        data_freq.close();
        dataspace_freq.close();
        mspace_freq.close();
        */
        
        const H5std_string DATASET_map("map_beam");
        DataSet data_map = file.openDataSet( DATASET_map );
        DataSpace dataspace_map = data_map.getSpace();
        int RANK_map  = dataspace_map.getSimpleExtentNdims();
        hsize_t dims_map[RANK_map];
        dataspace_map.getSimpleExtentDims(dims_map);
        DataSpace mspace_map(RANK_map, dims_map);
    
        hsize_t start[RANK_map] = {(unsigned long long int) (sb_list[0] - 1), 
                                   (unsigned long long int) (freq_list[0] - 1), 
                                   0, 0};
        hsize_t count[RANK_map] = {(unsigned long long int) Nsb,
                                   (unsigned long long int) Nfreq,
                                   (unsigned long long int) Nx,
                                   (unsigned long long int) Ny};

        hsize_t start_out[RANK_map] = {0, 0, 0, 0};
        hsize_t count_out[RANK_map] = {(unsigned long long int) Nsb,
                                       (unsigned long long int) Nfreq,
                                       (unsigned long long int) Nx,
                                       (unsigned long long int) Ny};
        dataspace_map.selectHyperslab(H5S_SELECT_SET, count, start);
        mspace_map.selectHyperslab(H5S_SELECT_SET, count_out, start_out);
        map_arr  = new double[dataspace_map.getSimpleExtentNpoints()];
        data_map.read(map_arr, PredType::NATIVE_DOUBLE, mspace_map, dataspace_map); 
    }
    /*
    const H5std_string DATASET_freq("freq");
    const H5std_string DATASET_map("map_beam");
    const H5std_string DATASET_hit("nhit_beam");
    const H5std_string DATASET_rms("rms_beam");

    DataSet data_freq = file.openDataSet( DATASET_freq );
    DataSet data_map = file.openDataSet( DATASET_map );
    DataSet data_hit = file.openDataSet( DATASET_hit );
    DataSet data_rms = file.openDataSet( DATASET_rms );

    DataSpace dataspace_freq = data_freq.getSpace();
    DataSpace dataspace_map = data_map.getSpace();
    DataSpace dataspace_hit = data_hit.getSpace();
    DataSpace dataspace_rms = data_rms.getSpace();
    
    int RANK_freq  = dataspace_freq.getSimpleExtentNdims();
    int RANK_map  = dataspace_map.getSimpleExtentNdims();
    int RANK_hit  = dataspace_hit.getSimpleExtentNdims();
    int RANK_rms  = dataspace_rms.getSimpleExtentNdims();

    hsize_t dims_freq[RANK_freq];
    hsize_t dims_map[RANK_map];
    hsize_t dims_hit[RANK_hit];
    hsize_t dims_rms[RANK_rms];

    dataspace_freq.getSimpleExtentDims(dims_freq);
    dataspace_map.getSimpleExtentDims(dims_map);
    dataspace_hit.getSimpleExtentDims(dims_hit);
    dataspace_rms.getSimpleExtentDims(dims_rms);
    
    DataSpace mspace_freq(RANK_freq, dims_freq);
    DataSpace mspace_map(RANK_map, dims_map);
    DataSpace mspace_hit(RANK_hit, dims_hit);
    DataSpace mspace_rms(RANK_rms, dims_rms);

    int Nsb    = dims_freq[0];
    int Nfeed  = dims_map[0];
    int Nfreq  = dims_freq[1];

    freq_arr = new double[dataspace_freq.getSimpleExtentNpoints()];
    map_arr  = new double[dataspace_map.getSimpleExtentNpoints()];
    hit_arr  = new int[dataspace_hit.getSimpleExtentNpoints()];
    rms_arr  = new double[dataspace_rms.getSimpleExtentNpoints()];

    data_freq.read(freq_arr, PredType::NATIVE_DOUBLE, mspace_freq, dataspace_freq );
    data_hit.read(hit_arr, PredType::NATIVE_INT, mspace_hit, dataspace_hit );
    data_rms.read(rms_arr, PredType::NATIVE_DOUBLE, mspace_rms, dataspace_rms );
    data_map.read(map_arr, PredType::NATIVE_DOUBLE, mspace_map, dataspace_map );
    
    file.close();
    data_freq.close();
    data_x.close();
    data_y.close();
    data_map.close();
    data_hit.close();
    data_rms.close();
    dataspace_freq.close();
    dataspace_x.close();
    dataspace_y.close();
    dataspace_map.close();
    dataspace_hit.close();
    dataspace_rms.close();
    mspace_freq.close();
    mspace_x.close();
    mspace_y.close();
    mspace_map.close();
    mspace_hit.close();
    mspace_rms.close();
    */
}

MapObjects::~MapObjects()
{   
    delete [] x_arr;
    delete [] y_arr;
    delete [] map_arr;
    /*
    delete [] freq_arr;
    delete [] hit_arr;
    delete [] rms_arr;
    */
}