// Filament class. Includes nucleation if a MT reduces to zero length.
//
// Benjamin Dalton 08/10/2015

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>
#include <chrono>
#include <ctime>

using namespace std;

#include "utilities.h"
#include "filament.h"

void periodic_shift(int Nt, vector<vector<double>> & r_t, vector<vector<double>> & r_per_t, double L_x, double L_y, double L_z)
{    
    int t, x_shift = 0, y_shift = 0, z_shift = 0;
    double delta_x = 0.0, delta_y = 0.0, delta_z = 0.0;
    
    r_per_t[0][0] = r_t[0][0]; r_per_t[1][0] = r_t[1][0]; r_per_t[2][0] = r_t[2][0];
    
    for(t = 1; t < Nt; t++)
    {
        delta_x = r_t[0][t] - r_t[0][t-1]; delta_y = r_t[1][t] - r_t[1][t-1]; delta_z = r_t[2][t] - r_t[2][t-1];
        
        // calculate the periodic shift for the x-component
        if(delta_x >= L_x/2.0){delta_x = delta_x - L_x; x_shift = x_shift - 1;}
        if(delta_x <= -L_x/2.0){delta_x = delta_x + L_x; x_shift = x_shift + 1;}
        
        // calculate the periodic shift for the y-component
        if(delta_y >= L_y/2.0){delta_y = delta_y - L_y; y_shift = y_shift - 1;}
        if(delta_y <= -L_y/2.0){delta_y = delta_y + L_y; y_shift = y_shift + 1;}
        
        // calculate the periodic shift for the z-component
        if(delta_z >= L_z/2.0){delta_z = delta_z - L_z; z_shift = z_shift - 1;}
        if(delta_z <= -L_z/2.0){delta_z = delta_z + L_z; z_shift = z_shift + 1;}
        
        r_per_t[0][t] = (r_t[0][t] + x_shift*L_x);
        r_per_t[1][t] = (r_t[1][t] + y_shift*L_y);
        r_per_t[2][t] = (r_t[2][t] + z_shift*L_z);
  
    }
}

void periodic_scale(double tscale, int Nt, vector<vector<double>> & r_per_t)
{
    for(int t = 0; t < Nt; t++)
    {
        r_per_t[0][t] = tscale*r_per_t[0][t];
        r_per_t[1][t] = tscale*r_per_t[1][t];
        r_per_t[2][t] = tscale*r_per_t[2][t];
    }
}
