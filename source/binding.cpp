// Calculating the probability for a pair binding (state-1 to state-2) event to a specific site on a neighbouring
// filament. Here is just a signe site calculation along 1 filament, the loop over all filaments is done in the main code
// main code
//
// See the crosslinker.cpp and crosslink.h files. Used together.
// cl_kon_i hold the probability of binding to the specific neighbour filament in the list (list updated periodically
// in nlist_cl) the sum over cl_kon_i gives the total probability. Individual values cl_kon_i are stored because they
// are subsequentially reused to choose which filament will be the pair. Once choosen, the epsilon_j_calc function
// finds the location on that filament where the second cross-linker domain will bind

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>
#include <chrono>
#include <ctime>

using namespace std;

#include "binding.h"
#include "filament.h"
#include "crosslinker.h"

binding :: binding(){}

void binding :: binding_init(int n_dims, double eta_sites)
{
    
    delta_lambda_j = 5.0e-08;
 
    rcli.resize(0);
    v_cli.resize(0);
    d_ij.resize(0);
    d_0.resize(0);
    r_image.resize(0);
    
    for(unsigned int i = 0; i < n_dims; i++)
    {
        rcli.push_back(0.0);
        v_cli.push_back(0.0);
        d_ij.push_back(0.0);
        d_0.push_back(0.0);
        r_image.push_back(0.0);
    }
}

// Calculates the probability of a cross-linker in state-1 undergoing a binding event to a second filament, thus
// creating a state-2 binding. Need to loop through all potential binding sites on all filaments that are in
// close enough proximaty to join in pair binding.
// (note: this is a very in-efficient inplementation since all sites on a filament are considered, even when the
// the alighnment is such that only a small section of the filament is within reach. Intended simplification,
// never implemented)
double binding :: kon_calc(filament& mt_i, filament& mt_j, vector<double>& rj_image, crosslinker& cl_k, int N_dim, double r_0)
{
    N_lambda_j = floor(mt_j.fil_length/delta_lambda_j);
    
    // use for the mid-point integration scheme
    lambda_j_mp = 0.0;
    
    a = 0.0, b = 0.0, r = 0.0, r_sqr = 0.0, q = 0.0, w = 0.0, s = 0.0, t = 0.0;
    d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0;
    k_on = 0.0, eps_r = 0.0, e_eps_r = 0.0;
    
    for (unsigned int n=0; n < N_dim; n++)
    {
        q = q + (mt_i.r_com[n] - rj_image[n])*(mt_i.r_com[n] - rj_image[n]);
        w = w + (mt_i.r_com[n] - rj_image[n])*mt_i.u_vec[n];
        s = s + (mt_i.r_com[n] - rj_image[n])*mt_j.u_vec[n];
        t = t + mt_i.u_vec[n]*mt_j.u_vec[n];
    }
    
    w = 2*cl_k.epsilon_i*w;
    s = 2*s;
    t = 2*cl_k.epsilon_i*t;
    
    a = q + w + cl_k.epsilon_i*cl_k.epsilon_i;
    b = s + t;
    
    for (unsigned int j=0; j < N_lambda_j; j++)
    {
        
        lambda_j_mp = -(mt_j.fil_length - delta_lambda_j)/2.0 + j*delta_lambda_j;
        
        d_ij_mag = 0.0;
        
        for(n=0; n<N_dim; n++)
        {
            d_ij[n] = (mt_i.r_com[n] + cl_k.epsilon_i*mt_i.u_vec[n]) - (rj_image[n] + lambda_j_mp*mt_j.u_vec[n]);
            d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
        }
        
        d_ij_mag = sqrt(d_ij_mag);
        
        if(d_ij_mag < 2.0e-7)
        {
            d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0;
            
            for(n=0; n<N_dim; n++)
            {
                d_0[n] = (r_0/d_ij_mag)*d_ij[n];
            }
            
            for(n=0; n<N_dim; n++)
            {
                d1 = d1 + (mt_i.r_com[n] - rj_image[n])*d_0[n];
                d2 = d2 + mt_i.u_vec[n]*d_0[n];
                d3 = d3 + d_0[n]*d_0[n];
                d4 = d4 + mt_j.u_vec[n]*d_0[n];
            }
    
            r_sqr       = (a - 2.0*d1 - 2.0*cl_k.epsilon_i*d2 + d3) - (b - 2.0*d4)*lambda_j_mp + lambda_j_mp*lambda_j_mp;
            r           = sqrt(r_sqr);
            eps_r       = cl_k.k_spr_kbT*(cl_k.cl_delta*r - 0.5*r*r);
            e_eps_r     = exp(eps_r);
            
            k_on        = k_on + e_eps_r*delta_lambda_j;
        }
    }
    
    if(k_on < 1.0e-12){k_on = 0;}
    
    return k_on;
}

// Only used once the main loop has determined that a pair binding event will definitely take place
double binding :: epsilon_j_calc(filament& mt_i, filament& mt_j, vector<double>& rj_image, crosslinker& cl_k, int N_dim, int m_id, double rand_in, double r_0)
{
    N_lambda_j = floor(mt_j.fil_length/delta_lambda_j);
    
    // use for the mid-point integration scheme
    lambda_j_mp = 0.0;
    lambda_0 = 0.0;
    
    a = 0.0, b = 0.0, r = 0.0, r_sqr = 0.0, q = 0.0, w = 0.0, s = 0.0, t = 0.0;
    k_on = 0.0, lambda_accum = 0.0, eps_r = 0.0, e_eps_r = 0.0;
    
    for (unsigned int n=0; n < N_dim; n++)
    {
        q = q + (mt_i.r_com[n] - rj_image[n])*(mt_i.r_com[n] - rj_image[n]);
        w = w + (mt_i.r_com[n] - rj_image[n])*mt_i.u_vec[n];
        s = s + (mt_i.r_com[n] - rj_image[n])*mt_j.u_vec[n];
        t = t + mt_i.u_vec[n]*mt_j.u_vec[n];
    }
    
    w = 2*cl_k.epsilon_i*w;
    s = 2*s;
    t = 2*cl_k.epsilon_i*t;
    
    a = q + w + cl_k.epsilon_i*cl_k.epsilon_i;
    b = s + t;
    
    // Perform loop to accumulate
    for (unsigned int j=0; j < N_lambda_j; j++)
    {
        lambda_j_mp = -(mt_j.fil_length - delta_lambda_j)/2.0 + j*delta_lambda_j;
        
        d_ij_mag = 0.0;
        
        for(n=0; n<N_dim; n++)
        {
            d_ij[n] = (mt_i.r_com[n] + cl_k.epsilon_i*mt_i.u_vec[n]) - (rj_image[n] + lambda_j_mp*mt_j.u_vec[n]);
            d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
        }
        
        d_ij_mag = sqrt(d_ij_mag);
        
        if(d_ij_mag < 2.0e-7)
        {
            d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0;
            
            for(n=0; n<N_dim; n++)
            {
                d_0[n] = (r_0/d_ij_mag)*d_ij[n];
            }
            
            for(n=0; n<N_dim; n++)
            {
                d1 = d1 + (mt_i.r_com[n] - rj_image[n])*d_0[n];
                d2 = d2 + mt_i.u_vec[n]*d_0[n];
                d3 = d3 + d_0[n]*d_0[n];
                d4 = d4 + mt_j.u_vec[n]*d_0[n];
            }
            
            r_sqr       = (a - 2.0*d1 - 2.0*cl_k.epsilon_i*d2 + d3) - (b - 2.0*d4)*lambda_j_mp + lambda_j_mp*lambda_j_mp;
            r           = sqrt(r_sqr);
            eps_r       = cl_k.k_spr_kbT*(cl_k.cl_delta*r - 0.5*r*r);
            e_eps_r     = exp(eps_r);
            
            k_on        = k_on + e_eps_r*delta_lambda_j;
        }
        
        if(rand_in < k_on/cl_k.cl_kon_i[m_id]){lambda_0 = lambda_j_mp; break;}
    }
    return lambda_0;
}

// a list of all filaments that are within interaction reach from a state-1 cross-linker. Need to loop through all
// potential filament to which a 2-state pair binding can attatch. Thus keep an update list, used for the
// kon_calc function
double binding :: nlist_cl(vector<double>& rj_com,double cl_r[3],vector<double>& ui_vec,double L_i,int N_dims)
{
    
    v_mag = 0.0;
    
    rcli = {0.0};
    v_cli = {0.0};
    
    lambda_min = 0.0;
    
    for(n = 0; n < N_dims; n++){rcli[n] = cl_r[n] - rj_com[n];}
    
    for(n = 0; n < N_dims; n++){lambda_min = lambda_min + rcli[n]*ui_vec[n];}
    
    if(lambda_min < -L_i/2)
    {
        lambda_0 = -L_i/2;
    }
    else if(lambda_min <= L_i/2 && lambda_min >= -L_i/2)
    {
        lambda_0 = lambda_min;
    }
    else if(lambda_min > L_i/2)
    {
        lambda_0 = L_i/2;
    }
    
    for(n=0; n < N_dims; n++)
    {
        v_cli[n] = rj_com[n] + lambda_0*ui_vec[n] - cl_r[n];
        v_mag   = v_mag + v_cli[n]*v_cli[n];
    }
    return sqrt(v_mag);
}



vector<double> binding :: image_calc_BC_cl(int BC_TYPE, int N_DIMS, vector<double>& ri_com, vector<double>& rj_com, vector<double>& box_vec)
{
    
    for(int n=0; n < N_DIMS; n++)
    {
        r_image[n] = ri_com[n];
    }
    if(BC_TYPE == 1)
    {
        for(int n=0; n < N_DIMS; n++)
        {
            if(ri_com[n] - rj_com[n] > box_vec[n]/2.0)
            {
                r_image[n] = ri_com[n]-box_vec[n];
            }
            else if(ri_com[n] - rj_com[n]< -box_vec[n]/2.0)
            {
                r_image[n] = ri_com[n] + box_vec[n];
            }
        }
    }
    if(BC_TYPE == 2)
    {
        if(ri_com[1] - rj_com[1] > box_vec[1]/2.0)
        {
            r_image[1] = ri_com[1]-box_vec[1];
        }
        else if(ri_com[1] - rj_com[1] < -box_vec[1]/2.0)
        {
            r_image[1] = ri_com[1] + box_vec[1];
        }
    }
    if(BC_TYPE == 4)
    {
        for(int n=0; n < 2; n++)
        {
            if(ri_com[n] - rj_com[n] > box_vec[n]/2.0)
            {
                r_image[n] = ri_com[n]-box_vec[n];
            }
            else if(ri_com[n] - rj_com[n]< -box_vec[n]/2.0)
            {
                r_image[n] = ri_com[n] + box_vec[n];
            }
        }
    }
    
    return r_image;
 
    for(int n=0; n < N_DIMS; n++)
    {
        r_image[n] = 0.0;
    }
    
}



