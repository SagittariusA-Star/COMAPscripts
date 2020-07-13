# include "MapObjects.h"
int main(int argc, char *argv[]){
MapObjects map = MapObjects(argc, argv);
map.read_map();
for (int i = 0; i < 4; i++){
    for (int j = 0; j < 64; j++){
        std::cout << " " << map.freq_arr[i * 64 + j] << " ";
    }
    std::cout << std::endl;
}

}

