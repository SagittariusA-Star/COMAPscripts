# include "MapObjects.h"

int main(int argc, char *argv[]){
    MapObjects map = MapObjects();
    map.read_map(argc, argv);
    
    /*
    for (int i = 0; i < 4;i++){
        for (int j = 0; j < 10;j++){
            cout << map.freq_arr[i * 64 + j] << " ";// << "," << i << "," << j << ")" << " ";
        }   
        cout << endl; 
    }
    cout << map.freq_arr[1 * 64 + 1] << endl;
    */
   for (int i = 0; i < 120; i++){
       cout << map.map_arr[(120*(64 * 1 + 29) + 49) + i] << endl;
   }
   cout << " hei " << endl;

}
