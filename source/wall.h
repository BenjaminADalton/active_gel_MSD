// Filament class. Includes nucleation if a MT reduces to zero length.
//
// Benjamin Dalton 08/10/2015

#ifndef wall_H
#define wall_H

#include <vector>
#include <array>
#include <math.h>

using namespace std;

class wall
{
    
    private:
    public:
    
    int x1_n,x2_n,z1_n,z2_n;
    
    vector<int> wl_x1,wl_x2,wl_z1,wl_z2;
    
    wall();
    void wall_init();
    void wall_x1_add(unsigned int fil_id);
    void wall_x2_add(unsigned int fil_id);
    void wall_z1_add(unsigned int fil_id);
    void wall_z2_add(unsigned int fil_id);
    
    void wall_2d_remove(unsigned int fil_id, unsigned int end_id);
    void wall_3d_remove(unsigned int fil_id, unsigned int end_id);
    
};

#endif