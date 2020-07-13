# ifdef MAPOBJECTS_H
# define MAPOBJECTS_H

# include "MapObjects.h"

class MapUtils {
    public: 
        MapUtils();
        MapUtils(MapObject map);
        MapUtils(MapObject map1, map2);
        void add();
        void weighted_sum();
        void subtract();
        void multiply();
        void divide();
        void perform_operation();

    private:

    
}





# endif