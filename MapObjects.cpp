# include "MapObjects.h"

MapObjects::MapObjects(int _argc, char *_argv[]){
    argc = _argc;
    for (int i = 0; i < argc; i++){
        argv[i] = _argv[i];
    }
}

void MapObjects::usage(){
    printf("These are the supported operations:\n");
}

void MapObjects::command_input(){
    opterr = 0;

    while ((opt = getopt(argc, argv, "s:i:p:f:h:d:o:x:y:j:z:w:m:r:l:n:")) != -1)
        switch (opt)
        {
        case 's':
            n = strlen(optarg);
            list_input = atoi(optarg);
            for (int i = n - 1; i >= 0; i--){
                sb = (list_input / (int) pow(10, i)) % 10;
                if (sb == 0){
                    printf("Please provide the number of the side band in 1-base, not 0-base!\n");
                    exit(EXIT_FAILURE);
                }
                else if (sb > 4){
                    printf("There are only four side bands!\n");
                    exit(EXIT_FAILURE);
                }
                sb_list.push_back(sb);
            }    
            break;
        case 'i':
            infile = optarg;
            break;
        case 'p':
            plot = optarg;
            if (plot != "map" && plot != "rms" && 
                plot != "hit" && plot != "map/rms"){
                printf("Make sure one of the valid map modes is chosen; 'map', 'rms', 'hit' or 'map/rms'.\n");
                exit(EXIT_FAILURE);
            }
            break;
        case 'f':
            freq = atoi(optarg);
            if (freq == 0){
                printf("Please provide the number of the frequency channels with 1-base, not 0-base!\n");
                exit(EXIT_FAILURE);
            }
            if (freq > 64){
                printf("There are only 64 frequencies pr. side band!\n");
                exit(EXIT_FAILURE);
            }
            
            break;
        case 'h':
            usage();
            break;
        case 'd':
            feed = atoi(optarg);
            if (feed == 0){
                printf("Please provide the number of the feed/detector with 1-base, not 0-base!\n");
                exit(EXIT_FAILURE);
            }
            if (feed > 20){
                printf("There are only 20 feeds/detectors!\n");
                exit(EXIT_FAILURE);
            }
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'j':
            jupiter == true;
            break;
        case 'z':
            sim_numb = atoi(optarg);
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

void MapObjects::read_map(){
    const H5std_string FILE_NAME("co6_good_map.h5");
    const H5std_string DATASET_freq( "freq" );
    const H5std_string DATASET_x( "x" );
    const H5std_string DATASET_y( "y" );
    const H5std_string DATASET_map( "map_beam" );
    const H5std_string DATASET_hit( "nhit_beam" );
    const H5std_string DATASET_rms( "rms_beam" );

    H5File file( FILE_NAME, H5F_ACC_RDONLY);

    DataSet data_freq = file.openDataSet( DATASET_freq );
    DataSet data_x = file.openDataSet( DATASET_x );
    DataSet data_y = file.openDataSet( DATASET_y );
    DataSet data_map = file.openDataSet( DATASET_map );
    DataSet data_hit = file.openDataSet( DATASET_hit );
    DataSet data_rms = file.openDataSet( DATASET_rms );

    DataSpace dataspace_freq = data_freq.getSpace();
    DataSpace dataspace_x = data_x.getSpace();
    DataSpace dataspace_y = data_y.getSpace();
    DataSpace dataspace_map = data_map.getSpace();
    DataSpace dataspace_hit = data_hit.getSpace();
    DataSpace dataspace_rms = data_rms.getSpace();
    
    int RANK_freq  = dataspace_freq.getSimpleExtentNdims();
    int RANK_x  = dataspace_x.getSimpleExtentNdims();
    int RANK_y  = dataspace_y.getSimpleExtentNdims();
    int RANK_map  = dataspace_map.getSimpleExtentNdims();
    int RANK_hit  = dataspace_hit.getSimpleExtentNdims();
    int RANK_rms  = dataspace_rms.getSimpleExtentNdims();

    hsize_t dims_freq[RANK_freq];
    hsize_t dims_x[RANK_x];
    hsize_t dims_y[RANK_y];
    hsize_t dims_map[RANK_map];
    hsize_t dims_hit[RANK_hit];
    hsize_t dims_rms[RANK_rms];

    dataspace_freq.getSimpleExtentDims(dims_freq);
    dataspace_x.getSimpleExtentDims(dims_x);
    dataspace_y.getSimpleExtentDims(dims_y);
    dataspace_map.getSimpleExtentDims(dims_map);
    dataspace_hit.getSimpleExtentDims(dims_hit);
    dataspace_rms.getSimpleExtentDims(dims_rms);
    
    DataSpace mspace_freq(RANK_freq, dims_freq);
    DataSpace mspace_x(RANK_x, dims_x);
    DataSpace mspace_y(RANK_y, dims_y);
    DataSpace mspace_map(RANK_map, dims_map);
    DataSpace mspace_hit(RANK_hit, dims_hit);
    DataSpace mspace_rms(RANK_rms, dims_rms);

    int Nsb    = dims_freq[0];
    int Nfeed  = dims_map[0];
    int Nfreq  = dims_freq[1];
    int Nx     = dims_x[0];
    int Ny     = dims_y[0];

    freq_arr = new double[dataspace_freq.getSimpleExtentNpoints()];
    x_arr    = new double[dataspace_x.getSimpleExtentNpoints()];
    y_arr    = new double[dataspace_y.getSimpleExtentNpoints()];
    map_arr  = new double[dataspace_map.getSimpleExtentNpoints()];
    hit_arr  = new int[dataspace_hit.getSimpleExtentNpoints()];
    rms_arr  = new double[dataspace_rms.getSimpleExtentNpoints()];

    data_freq.read(freq_arr, PredType::NATIVE_DOUBLE, mspace_freq, dataspace_freq );
    data_x.read(x_arr, PredType::NATIVE_DOUBLE, mspace_x, dataspace_x );
    data_y.read(y_arr, PredType::NATIVE_DOUBLE, mspace_y, dataspace_y );
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
}

MapObjects::~MapObjects()
{
    delete [] freq_arr;
    delete [] x_arr;
    delete [] y_arr;
    delete [] map_arr;
    delete [] hit_arr;
    delete [] rms_arr;
}