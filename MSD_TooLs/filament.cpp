#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "utilities.h"
#include "filament.h"

using namespace std;

filament :: filament(){}

void filament :: filament_init(int n_dim, int Nt_1, int Nt_2, int Nt_3)
{    
    int i,n;
    
    vector<double> buffer_1, buffer_2, buffer_3;
    
    for (i = 0; i < Nt_1; i++){buffer_1.push_back(0.0); t1.push_back(0.0); l_t1.push_back(0.0);}
    for (i = 0; i < Nt_2; i++){buffer_2.push_back(0.0); t2.push_back(0.0); l_t2.push_back(0.0);}
    for (i = 0; i < Nt_3; i++){buffer_3.push_back(0.0); t3.push_back(0.0); l_t3.push_back(0.0);}

    for (n = 0; n < n_dim; n++)
    {
        r_t1.push_back(buffer_1);
        r_t2.push_back(buffer_2);
        r_t3.push_back(buffer_3);
        
        r_per_t1.push_back(buffer_1);
        r_per_t2.push_back(buffer_2);
        r_per_t3.push_back(buffer_3);
        
        u_t1.push_back(buffer_1);
        u_t2.push_back(buffer_2);
        u_t3.push_back(buffer_3);
        
        u_per_t1.push_back(buffer_1);
        u_per_t2.push_back(buffer_2);
        u_per_t3.push_back(buffer_3);
    }
}









