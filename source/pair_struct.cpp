// Gear Predictor/corrector class.
//
// Benjamin Dalton 03/19/2016

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "pair_struct.h"

using namespace std;

pair_struct :: pair_struct(){}

void pair_struct :: pair_init(unsigned int n_dims)
{
    r_ij.resize(0);
    r_min_temp.resize(0);
    
    for(unsigned int i = 0; i < n_dims; i++)
    {
        r_ij.push_back(0.0);
        w_ij.push_back(0.0);
        r_min_temp.push_back(0.0);
    }
}

void pair_struct :: pair_add_new(unsigned int n_fils, unsigned int n_dims)
{
    
    // calculate the number of pairs for N_FILS
    n_pair = (n_fils*(n_fils-1))/2;
    
    vector<unsigned int> pair_buffer;
    pair_buffer.push_back(0.0);
    pair_buffer.push_back(0.0);
    pair_buffer.push_back(0.0);
    
    for (unsigned int i = 0; i < n_fils - 1; i++)
    {
        pair_buffer[0] = n_fils-1;
        pair_buffer[1] = i;
        pair_buffer[2] = 0;
        pair_list.push_back(pair_buffer);
    }
    
    vector<double> buffer_doub;
    buffer_doub.push_back(0.0);
    buffer_doub.push_back(0.0);
    
    vector<double> row_buffer;
    for(unsigned int i = 0; i < n_dims; i++)
    {
        row_buffer.push_back(0.0);
    }
    
    for(unsigned int i = 0; i < n_fils - 1; i++)
    {
        lambda_array.push_back(buffer_doub);
        uw_ij.push_back(row_buffer);
        w_ij_mag.push_back(0.0);
    }
}

void pair_struct :: nlist_init_state(unsigned int n_dims)
{
    w_ij_mag_single = 0;
    uw_ij_single.resize(0);
    lambda_single.resize(0);
    
    lambda_single.push_back(0.0);
    lambda_single.push_back(0.0);
    
    for(unsigned int i = 0; i < n_dims; i++)
    {
        uw_ij_single.push_back(0.0);
    }
}

void pair_struct :: outer_nlist_init(unsigned int n_dims)
{
    
    n_outer_pair = 0;
    outer_nlist.resize(0);
    
    w_ij_mag_single = 0;
    uw_ij_single.resize(0);
    lambda_single.resize(0);
    F_ij_vec.resize(0);
    
    lambda_single.push_back(0.0);
    lambda_single.push_back(0.0);
    
    for(unsigned int i = 0; i < n_dims; i++)
    {
        uw_ij_single.push_back(0.0);
    }
}

void pair_struct :: outer_nlist_build(unsigned int fil_i, unsigned int fil_j)
{
    n_outer_pair = n_outer_pair + 1;

    vector<unsigned int> id_buff;
    id_buff.push_back(fil_i);
    id_buff.push_back(fil_j);
    
    outer_nlist.push_back(id_buff);
    
}


void pair_struct :: inner_nlist_init(unsigned int n_dims)
{
    
    n_inner_pair = 0;
    w_ij_mag_single = 0.0;
    
    inner_nlist.resize(0);
    
    uw_ij_single.resize(0);
    lambda_single.resize(0);
    
    lambda_array.resize(0);
    uw_ij.resize(0);
    w_ij_mag.resize(0);
    F_ij_vec.resize(0);
    
    lambda_single.push_back(0.0);
    lambda_single.push_back(0.0);
    
    for(unsigned int i = 0; i < n_dims; i++)
    {
        uw_ij_single.push_back(0.0);
    }
}

void pair_struct :: inner_nlist_build(unsigned int fil_i, unsigned int fil_j)
{
    
    vector<unsigned int> id_buff; id_buff.resize(0);
    id_buff.push_back(fil_i);
    id_buff.push_back(fil_j);
    id_buff.push_back(n_inner_pair);
    
    n_inner_pair = n_inner_pair + 1;
 
    inner_nlist.push_back(id_buff);
 
    lambda_array.push_back(lambda_single);
    uw_ij.push_back(uw_ij_single);
    w_ij_mag.push_back(w_ij_mag_single);
    F_ij_vec.push_back(0.0);
    
}

void pair_struct :: inner_pair_constructor(unsigned int n_dims)
{
    
    lambda_array.resize(0);
    uw_ij.resize(0);
    w_ij_mag.resize(0);
    F_ij_vec.resize(0);
    
    vector<double> pair_buffer_doub;
    pair_buffer_doub.push_back(0.0);
    pair_buffer_doub.push_back(0.0);
    
    vector<double> row_buffer;
    for(unsigned int i = 0; i < n_dims; i++)
    {
        row_buffer.push_back(0.0);
    }
    
    for(unsigned int i = 0; i < n_inner_pair; i++)
    {
        lambda_array.push_back(pair_buffer_doub);
        uw_ij.push_back(row_buffer);
        w_ij_mag.push_back(0.0);
        F_ij_vec.push_back(0.0);
    }
}

void pair_struct :: least_sep_calc(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec ,vector<double>& uj_vec,double L_i, double L_j,int N_dims, unsigned int pair_id)

{
    
    r_ij = {0.0}; w_ij = {0.0}; r_min_temp = {0.0};
    
    lambda[0] = 0.0; lambda[1] = 0.0;
    abs_r_min_a = 0.0;
    abs_r_min_b = 0.0;
    
    unsigned int i,j;
    
    for(i = 0; i < N_dims; i++){r_ij[i] = ri_com[i] - rj_com[i];}

    u_ij      = 0.0;
    ui_rij_dp = 0.0;
    uj_rij_dp = 0.0;
    
    for(j=0; j < N_dims; j++){u_ij = u_ij + ui_vec[j]*uj_vec[j];}
    for(j=0; j < N_dims; j++){ui_rij_dp = ui_rij_dp + ui_vec[j]*r_ij[j];}
    for(j=0; j < N_dims; j++){uj_rij_dp = uj_rij_dp + uj_vec[j]*r_ij[j];}
    
    K = 1.0/(1.0-u_ij*u_ij);
    
    la1_0 = K*(-ui_rij_dp + u_ij*uj_rij_dp);
    la2_0 = K*(uj_rij_dp - u_ij*ui_rij_dp);
    
    abs_uij = abs(u_ij), abs_ui_rij_dp = abs(ui_rij_dp), abs_uj_rij_dp = abs(uj_rij_dp), abs_la1_0 = abs(la1_0), abs_la2_0 = abs(la2_0);
    
    la1_a = 0.0,la2_a = 0.0,la1_b = 0.0,la2_b = 0.0;
    
    if(abs_la1_0 < L_i/2.0 && abs_la2_0 < L_j/2.0)
    {
        lambda[0] = la1_0;
        lambda[1] = la2_0;
    }
    if(1-abs_uij < PREC)
    {
        if ((L_i/(L_i+L_j))*abs_ui_rij_dp >= L_i/2.0)
        {
            if (uj_rij_dp >= 0.0 && uj_rij_dp <= pi/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= pi/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -pi/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -pi/2.0)
            {lambda[0] = -L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -pi/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= pi/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp >= 0.0 && uj_rij_dp <= pi/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -pi/2.0)
            {lambda[0] = -L_i/2.0; lambda[1] = L_j/2.0;}
        }
        else if((L_i/(L_i+L_j))*abs_ui_rij_dp < L_i/2.0)
        {
            lambda[0] = -(L_i/(L_i+L_j))*ui_rij_dp;
            lambda[1] = (L_j/(L_i+L_j))*uj_rij_dp;
        }
    }
    else if(abs_la1_0 >= L_i/2.0 && abs_la2_0 >= L_j/2.0)
    {
        
        // signbit boolean: if positive output 0 (no to negative), if negative output 1 (yes to negative)
        if(signbit(la1_0)==0){la1_a = L_i/2;} else if(signbit(la1_0)==1){la1_a = -L_i/2.0;}
        la2_a = la1_a*u_ij + uj_rij_dp;
        
        if(abs(la2_a) > L_j/2.0)
        {
            if(signbit(la2_a)==0){la2_a = L_j/2.0; } // pos
            else if(signbit(la2_a)==1){la2_a = -L_j/2.0; } // neg
        }

        for(j=0; j < N_dims; j++){r_min_temp[j] = r_ij[j] + la1_a*ui_vec[j] - la2_a*uj_vec[j];}
        for(j=0; j < N_dims; j++){abs_r_min_a = abs_r_min_a + r_min_temp[j]*r_min_temp[j];}

        abs_r_min_a = sqrt(abs_r_min_a);
        
        if(signbit(la2_0)==0){la2_b = L_j/2;} else if(signbit(la2_0)==1){la2_b = -L_j/2.0;}
        la1_b = la2_b*u_ij - ui_rij_dp;
        
        if(abs(la1_b) > L_i/2.0)
        {
            if(signbit(la1_b)==0){la1_b = L_i/2.0;}
            else if(signbit(la1_b)==1){la1_b = -L_i/2.0;}
        }

        for(j=0; j < N_dims; j++){r_min_temp[j] = r_ij[j] + la1_b*ui_vec[j] - la2_b*uj_vec[j];}
        for(j=0; j < N_dims; j++){abs_r_min_b = abs_r_min_b + r_min_temp[j]*r_min_temp[j];}

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
    
    w_mag  = 0.0;
    ri_min = 0.0;
    rj_min = 0.0;
    
    for(j=0; j < N_dims; j++)
    {
        
        ri_min = ri_com[j] + lambda[0]*ui_vec[j];
        rj_min = rj_com[j] + lambda[1]*uj_vec[j];
        
        w_ij[j]   = ri_min - rj_min;
        w_mag     = w_mag + w_ij[j]*w_ij[j];
    }
    
    w_ij_mag[pair_id] = sqrt(w_mag);
    
    for(j = 0; j < N_dims; j++)
    {
        uw_ij[pair_id][j] = w_ij[j]/w_ij_mag[pair_id];
    }
    
    lambda_array[pair_id][0] = lambda[0];
    lambda_array[pair_id][1] = lambda[1];
}

void pair_struct :: least_sep_calc_single(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec ,vector<double>& uj_vec, double L_i, double L_j,int N_dims)

{
    r_ij = {0.0};
    w_ij = {0.0};
    
    r_min_temp = {0.0};
    
    lambda[0] = 0.0; lambda[1] = 0.0;
    abs_r_min_a = 0.0;
    abs_r_min_b = 0.0;
    
    unsigned int i,j;
    
    for(i = 0; i < N_dims; i++){r_ij[i] = ri_com[i] - rj_com[i];}
    
    u_ij = 0.0, ui_rij_dp = 0.0, uj_rij_dp = 0.0;
    
    for(j=0; j < N_dims; j++){u_ij = u_ij + ui_vec[j]*uj_vec[j];}
    for(j=0; j < N_dims; j++){ui_rij_dp = ui_rij_dp + ui_vec[j]*r_ij[j];}
    for(j=0; j < N_dims; j++){uj_rij_dp = uj_rij_dp + uj_vec[j]*r_ij[j];}
    
    K = 1.0/(1.0-u_ij*u_ij);
    
    la1_0 = K*(-ui_rij_dp + u_ij*uj_rij_dp);
    la2_0 = K*(uj_rij_dp - u_ij*ui_rij_dp);
    
    abs_uij = abs(u_ij), abs_ui_rij_dp = abs(ui_rij_dp), abs_uj_rij_dp = abs(uj_rij_dp), abs_la1_0 = abs(la1_0), abs_la2_0 = abs(la2_0);
    
    la1_a = 0.0,la2_a = 0.0,la1_b = 0.0,la2_b = 0.0;
    
    abs_r_min_a = 0.0, abs_r_min_b = 0.0;
    
    intercept_switch_2 = 0;
    
    if(abs_la1_0 <= L_i/2.0 && abs_la2_0 <= L_j/2.0)
    {
        lambda[0] = la1_0;
        lambda[1] = la2_0;
        
        intercept_switch_2 = 1;
    }
    if(1-abs_uij < PREC)
    {
        if ((L_i/(L_i+L_j))*abs_ui_rij_dp >= L_i/2.0)
        {
            if (uj_rij_dp >= 0.0 && uj_rij_dp <= pi/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= pi/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -pi/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -pi/2.0)
            {lambda[0] = -L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp <= 0.0 && uj_rij_dp >= -pi/2.0 && -ui_rij_dp >= 0.0 && -ui_rij_dp <= pi/2.0)
            {lambda[0] = L_i/2.0; lambda[1] = -L_j/2.0;}
            else if (uj_rij_dp >= 0.0 && uj_rij_dp <= pi/2.0 && -ui_rij_dp <= 0.0 && -ui_rij_dp >= -pi/2.0)
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
        if(signbit(la1_0)==0){la1_a = L_i/2;} else if(signbit(la1_0)==1){la1_a = -L_i/2.0;}
        la2_a = la1_a*u_ij + uj_rij_dp;
        
        if(abs(la2_a) > L_j/2.0)
        {
            if(signbit(la2_a)==0){la2_a = L_j/2.0; } // pos
            else if(signbit(la2_a)==1){la2_a = -L_j/2.0; } // neg
        }
        
        for(j=0; j < N_dims; j++){r_min_temp[j] = r_ij[j] + la1_a*ui_vec[j] - la2_a*uj_vec[j];}
        for(j=0; j < N_dims; j++){abs_r_min_a = abs_r_min_a + r_min_temp[j]*r_min_temp[j];}
        
        abs_r_min_a = sqrt(abs_r_min_a);
        
        if(signbit(la2_0)==0){la2_b = L_j/2;} else if(signbit(la2_0)==1){la2_b = -L_j/2.0;}
        la1_b = la2_b*u_ij - ui_rij_dp;
        
        if(abs(la1_b) > L_i/2.0)
        {
            if(signbit(la1_b)==0){la1_b = L_i/2.0;}
            else if(signbit(la1_b)==1){la1_b = -L_i/2.0;}
        }
        
        for(j=0; j < N_dims; j++){r_min_temp[j] = r_ij[j] + la1_b*ui_vec[j] - la2_b*uj_vec[j];}
        for(j=0; j < N_dims; j++){abs_r_min_b = abs_r_min_b + r_min_temp[j]*r_min_temp[j];}
        
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

    w_mag = 0.0;
    
    ri_min = 0.0;
    rj_min = 0.0;
    
    for(j=0; j < N_dims; j++)
    {
        ri_min = ri_com[j] + lambda[0]*ui_vec[j];
        rj_min = rj_com[j] + lambda[1]*uj_vec[j];
        
        w_ij[j] = ri_min - rj_min;
        w_mag   = w_mag + w_ij[j]*w_ij[j];
    }
    
    w_ij_mag_single = sqrt(w_mag);
    
    for(j = 0; j < N_dims; j++)
    {
        uw_ij_single[j] = w_ij[j]/w_ij_mag_single;
    }
    
    lambda_single[0] = lambda[0];
    lambda_single[1] = lambda[1];
    
}



