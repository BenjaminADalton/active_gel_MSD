#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "utilities.h"
#include "filament.h"

using namespace std;

filament :: filament(){}

void filament :: filament_init(unsigned int n_dim)
{    
    for (unsigned int n = 0; n < n_dim; n++)
    {
        r_com.push_back(0.0);
        u_vec.push_back(0.0);
        r_init.push_back(0.0);
        u_init.push_back(0.0);
        r_accum.push_back(0.0);
        r_hist.push_back(0.0);
    }
    
    for (unsigned int n = 0; n < n_dim; n++){r_eqnm.push_back(0.0);}
    for (unsigned int n = 0; n < 3; n++){u_eqnm.push_back(0.0);}
    
    vector<double> buffer_u;
    
    if(n_dim == 2)
    {
        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        
        u_perp.push_back(buffer_u);
    }
    if(n_dim == 3)
    {

        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        buffer_u.push_back(0.0);
        
        u_perp.push_back(buffer_u);
        u_perp.push_back(buffer_u);
    }

    pbc_index[0] = 0.0;
    pbc_index[1] = 0.0;
    pbc_index[2] = 0.0;
    
    fil_id = 0;
    
}

void filament :: filament_drag_tao(double d,double kappa)
{
    
    log_ld = log(fil_length/d);
    
    para_poly = log_ld;
    perp_poly = log_ld;
    rot_poly  = log_ld;

    gamma_para = 2.0*kappa*fil_length/para_poly;
    gamma_perp = 4.0*kappa*fil_length/perp_poly;
    gamma_rot  = (kappa*fil_length*fil_length*fil_length)/(3.0*rot_poly);
    
}

void filament :: filament_drag_torre(double d,double kappa)
{
    
    log_ld     = log(fil_length/d);
    d_on_l     = d/fil_length;
    d_on_l_sqr = d_on_l*d_on_l;
    
    para_poly = log_ld - 0.207 + 0.980*d_on_l - 0.133*d_on_l_sqr;
    perp_poly = log_ld + 0.839 + 0.185*d_on_l + 0.233*d_on_l_sqr;
    rot_poly  = log_ld - 0.662 + 0.917*d_on_l - 0.050*d_on_l_sqr;
    
    gamma_para = 2.0*kappa*fil_length/para_poly;
    gamma_perp = 4.0*kappa*fil_length/perp_poly;
    gamma_rot  = (kappa*fil_length*fil_length*fil_length)/(3.0*rot_poly);
    
}







