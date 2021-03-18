// Benjamin Dalton 09/26/2017

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "wall_calc.h"

using namespace std;

wall_calc :: wall_calc(){}

void wall_calc :: wall_calc_init(int n_dims, double epsilon_in,double sig_shift_wall_in,double beta_sig_wall_in,double one_on_beta_wall_in)
{

    epsilon = epsilon_in;
    one_on_beta_wall = one_on_beta_wall_in;
    sig_shift_wall = sig_shift_wall_in;
    beta_sig_wall = beta_sig_wall_in;
    
    f_beta_wall  = -24.0*(epsilon/sig_shift_wall)*(2.0*pow(one_on_beta_wall,13.0) - pow(one_on_beta_wall,7.0));
    df_beta_wall = 24.0*(epsilon/sig_shift_wall)*(26.0*pow(sig_shift_wall,13.0)/pow(beta_sig_wall,14.0) - 7.0*pow(sig_shift_wall,7.0)/pow(beta_sig_wall,8.0));
    
    r_ij.resize(0);
    w_ij.resize(0);
    uw_ij_out.resize(0);
    ri_min.resize(0);
    rj_min.resize(0);
    r_min_a.resize(0);
    r_min_b.resize(0);

    for(unsigned int i = 0; i < n_dims; i++)
    {
        r_ij.push_back(0.0);
        w_ij.push_back(0.0);
        uw_ij_out.push_back(0.0);
        ri_min.push_back(0.0);
        rj_min.push_back(0.0);
        r_min_a.push_back(0.0);
        r_min_b.push_back(0.0);
    }
}

double wall_calc :: wall_wca_calc()
{
    F_i_wall   = 0.0;
    eps_on_sig = epsilon/sig_shift_wall;
    sig_on_wij = sig_shift_wall/w_ij_mag_out;
    
    if(w_ij_mag_out <= 1.0e-20 || isnan(w_ij_mag_out) == 1)
    {
        w_ij_mag_out = 1.0e-25;
    }
    
    if(w_ij_mag_out >= beta_sig_wall)
    {
        F_i_wall = -(24.0*eps_on_sig)*(2.0*pow(sig_on_wij,13.0) - pow(sig_on_wij,7.0));
    }
    else if(w_ij_mag_out < beta_sig_wall)
    {
        F_i_wall = f_beta_wall + df_beta_wall*(w_ij_mag_out - beta_sig_wall);
    }

    return F_i_wall;
}

void wall_calc :: wall_calc_class(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec,vector<double>& uj_vec, double L_i, double L_j, int N_dims, int id)
{
    
    lambda[0] = 0.0;
    lambda[1] = 0.0;
    
    rix_out = 2.0;
    
    uw_ij_out = {0.0};
    w_ij      = {0.0};
    r_ij      = {0.0};
    ri_min    = {0.0};
    rj_min    = {0.0};

    unsigned int i,j,n;
    
    for(i = 0; i < N_dims; i++){r_ij[i] = ri_com[i] - rj_com[i];}
    
    u_ij = 0.0; ui_rij_dp = 0.0; uj_rij_dp = 0.0;
    
    for(j=0; j < N_dims; j++){u_ij = u_ij + ui_vec[j]*uj_vec[j];}
    for(j=0; j < N_dims; j++){ui_rij_dp = ui_rij_dp + ui_vec[j]*r_ij[j];}
    for(j=0; j < N_dims; j++){uj_rij_dp = uj_rij_dp + uj_vec[j]*r_ij[j];}
    
    K = 1.0/(1.0-u_ij*u_ij);
    
    la1_0 = K*(-ui_rij_dp + u_ij*uj_rij_dp);
    la2_0 = K*(uj_rij_dp - u_ij*ui_rij_dp);
    
    abs_uij = abs(u_ij), abs_ui_rij_dp = abs(ui_rij_dp), abs_uj_rij_dp = abs(uj_rij_dp), abs_la1_0 = abs(la1_0), abs_la2_0 = abs(la2_0);
    
    la1_a = 0.0;
    la2_a = 0.0;
    la1_b = 0.0;
    la2_b = 0.0;
    
    r_min_a = {0.0};
    r_min_b = {0.0};
    abs_r_min_a = 0.0, abs_r_min_b = 0.0;
    
    if(abs_la1_0 < L_i/2.0 && abs_la2_0 < L_j/2.0)
    {
        lambda[0] = la1_0;
        lambda[1] = la2_0;
    }
    if(1-abs_uij < PREC3)
    {
        if ((L_i/(L_i+L_j))*abs_ui_rij_dp >= L_i/2.0)
        {
            if (uj_rij_dp >= 0.0 && uj_rij_dp <= PI/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= PI/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -PI/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -PI/2.0)
            {lambda[0] = -L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -PI/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= PI/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp >= 0.0 && uj_rij_dp <= PI/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -PI/2.0)
            {lambda[0] = -L_i/2.0; lambda[1] = L_j/2.0;}
        }
        else if((L_i/(L_i+L_j))*abs_ui_rij_dp < L_i/2.0)
        {
            lambda[0] = -(L_i/(L_i+L_j))*ui_rij_dp;
            lambda[1] = (L_j/(L_i+L_j))*uj_rij_dp;
        }
    }
    else if(abs_la1_0 > L_i/2.0 && abs_la2_0 > L_j/2.0)
    {
        // signbit boolean: if positive output 0 (no to negative), if negative output 1 (yes to negative)
        
        if(signbit(la1_0)==0){la1_a = L_i/2;}
        else if(signbit(la1_0)==1){la1_a = -L_i/2.0;}
        
        la2_a = la1_a*u_ij + uj_rij_dp;
        
        if(abs(la2_a) > L_j/2.0)
        {
            if(signbit(la2_a)==0){la2_a = L_j/2.0; } // pos
            else if(signbit(la2_a)==1){la2_a = -L_j/2.0; } // neg
        }
        
        for(j=0; j < N_dims; j++){r_min_a[j] = r_ij[j] + la1_a*ui_vec[j] - la2_a*uj_vec[j];}
        
        for(j=0; j < N_dims; j++){abs_r_min_a = abs_r_min_a + r_min_a[j]*r_min_a[j];}
        abs_r_min_a = sqrt(abs_r_min_a);
        
        
        if(signbit(la2_0)==0){la2_b = L_j/2;} else if(signbit(la2_0)==1){la2_b = -L_j/2.0;}
        la1_b = la2_b*u_ij - ui_rij_dp;
        
        if(abs(la1_b) > L_i/2.0)
        {
            if(signbit(la1_b)==0){la1_b = L_i/2.0;}
            else if(signbit(la1_b)==1){la1_b = -L_i/2.0;}
        }
        
        for(j=0; j < N_dims; j++){r_min_b[j] = r_ij[j] + la1_b*ui_vec[j] - la2_b*uj_vec[j];}
        for(j=0; j < N_dims; j++){abs_r_min_b = abs_r_min_b + r_min_b[j]*r_min_b[j];}
        abs_r_min_b = sqrt(abs_r_min_b);
        
        
        if (abs_r_min_a == abs_r_min_b){lambda[0] = la1_a; lambda[1] = la2_b;}
        else if(abs_r_min_a < abs_r_min_b){lambda[0] = la1_a; lambda[1] = la2_a;}
        else if (abs_r_min_b < abs_r_min_a){lambda[0] = la1_b; lambda[1] = la2_b;}
        
    }
    else if(abs_la1_0 <= L_i/2.0 && abs_la2_0 > L_j/2.0)
    {
        if(signbit(la2_0)==0){lambda[1] = L_j/2.0;} //no = pos
        else if(signbit(la2_0)==1){lambda[1] = -L_j/2.0;} //yes = neg
        
        lambda[0] = lambda[1]*u_ij - ui_rij_dp;
        
        if(abs(lambda[0]) > L_i/2.0)
        {
            if(signbit(lambda[0])==0){lambda[0] = L_i/2.0;}
            else if(signbit(lambda[0])==1){lambda[0] = -L_i/2.0;}
        }
    }
    else if(abs_la1_0 > L_i/2.0 && abs_la2_0 <= L_j/2.0)
    {
        if(signbit(la1_0)==0){lambda[0] = L_i/2.0;}
        else if(signbit(la1_0)==1){lambda[0] = -L_i/2.0;}
        
        lambda[1] = lambda[0]*u_ij + uj_rij_dp;
        
        if(abs(lambda[1]) > L_j/2.0)
        {
            if(signbit(lambda[1])==0){lambda[1] = L_j/2.0;}
            if(signbit(lambda[1])==1){lambda[1] = -L_j/2.0;}
        }
    }
    
    w_ij_mag = 0.0;
    
    for(n=0; n < N_dims; n++)
    {
        //r_ij_vec[n] = ri_com[n] - rj_com[n];
        ri_min[n]  = ri_com[n] + lambda[0]*ui_vec[n];
        rj_min[n]  = rj_com[n] + lambda[1]*uj_vec[n];
        w_ij[n]    = ri_min[n] - rj_min[n];
        w_ij_mag   = w_ij_mag + w_ij[n]*w_ij[n];
    }
    
    w_ij_mag = sqrt(w_ij_mag);
    
    w_ij_mag_out = w_ij_mag;
    
    
    for(n = 0; n < N_dims; n++)
    {
        uw_ij_out[n] = w_ij[n]/w_ij_mag;
    }
    
}
