// Periodic Boundary Conditions and minimal image calculation

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "pbc.h"

using namespace std;

pbc :: pbc(){}

void pbc :: pbc_constructor(unsigned int n_dims, double L_x, double L_y, double L_z)
{
    
    vector<double> row_buffer;
    for(unsigned int i = 0; i < n_dims; i++)
    {
        ri_image.push_back(0.0);
        rj_image.push_back(0.0);
        box_vec.push_back(0.0);
    }
    
    if(n_dims == 2)
    {
        box_vec[0] = L_x;
        box_vec[1] = L_y;
    }
    if(n_dims == 3)
    {
        box_vec[0] = L_x;
        box_vec[1] = L_y;
        box_vec[2] = L_z;
    }

}

void pbc :: minimal_image(vector<double>& ri_com,vector<double>& rj_com,int n_dims)
{
    for(unsigned int i = 0; i < n_dims; i++)
    {
        ri_image[i] = ri_com[i];
        rj_image[i] = rj_com[i];
    
        if(ri_com[i]-rj_com[i] > box_vec[i]/2.0){ri_image[i] = ri_com[i]-box_vec[i];}
        else if(ri_com[i]-rj_com[i] < -box_vec[i]/2.0){ri_image[i] = ri_com[i] + box_vec[i];}
    }
}


void pbc :: pbc_calc(vector<double>& r_com,int n_dims, double pbc_index[3])
{
    for(unsigned int n = 0; n < n_dims; n++)
    {
        if(r_com[n] > box_vec[n]){r_com[n] = r_com[n]-box_vec[n]; pbc_index[n] = pbc_index[n] + 1;}
        if(r_com[n] < 0){r_com[n] = r_com[n] + box_vec[n];; pbc_index[n] = pbc_index[n] - 1;}
    }
}

void pbc :: pbc_y_calc(vector<double>& r_com,int n_dims, double pbc_index[3])
{
    if(r_com[1] > box_vec[1]){r_com[1] = r_com[1]-box_vec[1]; pbc_index[1] = pbc_index[1] + 1;}
    if(r_com[1] < 0){r_com[1] = r_com[1] + box_vec[1];; pbc_index[1] = pbc_index[1] - 1;}
}













