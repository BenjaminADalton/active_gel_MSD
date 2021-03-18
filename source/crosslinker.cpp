#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "utilities.h"
#include "crosslinker.h"
#include "filament.h"

using namespace std;

crosslinker :: crosslinker(){}

void crosslinker :: cl_init_start(double delta_in, double ksp_in)
{
    cl_delta = delta_in;
    k_spr    = ksp_in;
}

void crosslinker :: cl_init(unsigned int n_dim, unsigned int n_fils, double delta_in, double ksp_in)
{
    epsilon_i = 0.0;
    epsilon_j = 0.0;
    
    cl_delta = delta_in;
    k_spr = ksp_in;
    
    for(unsigned int i = 0; i < n_fils - 1; i++)
    {
        cl_neigh.push_back(0);
        cl_kon_i.push_back(0.0);
    }
}

void crosslinker :: neigh_append(unsigned int fil_id)
{
    cl_neigh.push_back(fil_id);
    cl_kon_i.push_back(0.0);
}

void crosslinker :: nlist_clean()
{
    
    n_neigh = 0;
    
    cl_neigh.resize(0);
    cl_kon_i.resize(0);
}

void crosslinker :: cl_walk_1_type_k(double fil_length_i, double motor_vel, double dt)
{
    epsilon_i = epsilon_i + motor_vel*dt;
    if(epsilon_i >= fil_length_i/2.0){state = 0; nlist_clean(); end_state = 0;}
}

void crosslinker :: cl_walk_1_type_d(double fil_length_i, double motor_vel, double dt)
{
    epsilon_i = epsilon_i + motor_vel*dt;
    
    if(epsilon_i <= -fil_length_i/2.0)
    {
        epsilon_i = -fil_length_i/2.0;
        
        if(dyn_detach == 1){state = 0; nlist_clean(); end_state = 0;}
    }
}

void crosslinker :: cl_walk_1_type_d1(double fil_length_i, double motor_vel, double dt) // only the j-filament domain walks
{
    if(epsilon_i <= -fil_length_i/2.0)
    {
        epsilon_i = -fil_length_i/2.0;
        
        if(dyn_detach == 1){state = 0; nlist_clean(); end_state = 0;}
    }
}

void crosslinker :: cl_walk_2_type_k(double fil_length_i,double fil_length_j, double motor_vel, double dt)
{
    // Force-velocity relationship. If f_para_i > F_stall: stall
    if(abs(f_para_i) <= F_stall)
    {
            epsilon_i = epsilon_i + motor_vel*dt*(1.0-abs(f_para_i/F_stall));
    }
    
    if(epsilon_i >= fil_length_i/2.0){end_state_i = 1;}
    
    // Force-velocity relationship. If f_para_j > F_stall: stall
    if(abs(f_para_j) <= F_stall)
    {
        epsilon_j = epsilon_j + motor_vel*dt*(1.0-abs(f_para_j/F_stall));
    }
    
    if(epsilon_j > fil_length_j/2.0){end_state_j = 1;}
}

void crosslinker :: cl_walk_2_type_d(double fil_length_i,double fil_length_j, double motor_vel, double dt)
{
    // Force-velocity relationship. If f_para_i > F_stall: stall
    if(abs(f_para_i) <= F_stall)
    {
        epsilon_i = epsilon_i + motor_vel*dt*(1.0-abs(f_para_i/F_stall));
    }
    
    if(epsilon_i <= -fil_length_i/2.0)
    {
        epsilon_i = -fil_length_i/2.0;
        
        if(dyn_detach == 1){end_state_i = 1;}
    }
    
    // Force-velocity relationship. If f_para_j > F_stall: stall
    if(abs(f_para_j) <= F_stall)
    {
        epsilon_j = epsilon_j + motor_vel*dt*(1.0-abs(f_para_j/F_stall));
    }
    
    if(epsilon_j < -fil_length_j/2.0)
    {
        epsilon_j = -fil_length_j/2.0;
        
        if(dyn_detach == 1){end_state_j = 1;}
    }

}

void crosslinker :: cl_walk_2_type_d1(double fil_length_i,double fil_length_j, double motor_vel, double dt) // only the j-filament domain walks
{
    if(epsilon_i <= -fil_length_i/2.0)
    {
        epsilon_i = -fil_length_i/2.0;
        
        if(dyn_detach == 1){end_state_i = 1;}
    }
    
    // Force-velocity relationship. If f_para_j > F_stall: stall
    if(abs(f_para_j) <= F_stall)
    {
        epsilon_j = epsilon_j + motor_vel*dt*(1.0-abs(f_para_j/F_stall));
    }
    
    if(epsilon_j < -fil_length_j/2.0)
    {
        epsilon_j = -fil_length_j/2.0;
        
        if(dyn_detach == 1){end_state_j = 1;}
    }
    
}

void crosslinker :: least_sep_calc_cl(vector<double>& ri_com,double cl_r[3],vector<double>& ui_vec,double L_i,int N_dims)
{
    
    vector<double> rcli(N_dims); // Need to added these as new to the .h file, to hold until the add_new is made
    
    double lambda_0;
    double lambda_min = 0.0;
    unsigned int i,j,n;
    
    for(n = 0; n < N_dims; n++){rcli[n] = cl_r[n] - ri_com[n];}
    
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
    
    double v_mag = 0.0;
    
    double* v_cli = rcli.data();//reusing memory of rcli
    
    for(n=0; n < N_dims; n++)
    {
        v_cli[n] = ri_com[n] + lambda_0*ui_vec[n] - cl_r[n];
        v_mag   = v_mag + v_cli[n]*v_cli[n];
    }
    
    v_cli_mag_new = sqrt(v_mag);
    lambda_0_new  = lambda_0;
    
    for(n = 0; n < N_dims; n++)
    {
        uv_cli_new[n] = v_cli[n]/v_cli_mag_new;
        v_cli_out_new[n] = v_cli[n];
    }
}

