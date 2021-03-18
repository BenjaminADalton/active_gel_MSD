#ifndef crosslinker_H
#define crosslinker_H

using namespace std;

#include "utilities.h"
#include "filament.h"

class crosslinker
{
    
    private:
    public:
    
    unsigned int state = 0; // 0 = unbound, 1 = single bound, 2 = double bound
    unsigned int type  = 1; // 0 = passive, 1 = +end directed, 2 = -end directed
    
    double epsilon_i = 0.0, epsilon_j = 0.0;

    // motor parameters
    double cl_delta  = 0.0;
    double k_spr     = 0.0;
    double k_spr_kbT = 0.0;
    double k_on      = 0.0;
    double k_off     = 0.0;
    double F_stall   = 0.0;
    
    int dyn_detach = 0;
    
    int occ_time_on  = 0;
    int sing_time_on = 0;
    
    int end_state   = 0;
    int end_state_i = 0;
    int end_state_j = 0;
    
    double f_para_i = 0.0;
    double f_para_j = 0.0;
    
    double v_cli_mag_new = 0.0;
    double lambda_0_new  = 0.0;
    
    // set up filament proximate neighbors
    // Used in conjunction with the kon_calc, epsilon_j_calc and the binding.cpp functions.
    //
    unsigned int n_neigh = 0;
    vector<int> cl_neigh;
    vector<double> cl_kon_i;
    
    vector<double> uv_cli_new;
    vector<double> v_cli_out_new;
    
    crosslinker();

    unsigned int fil_i, fil_j;

    void nlist_clean();

    void cl_init(unsigned int n_dim, unsigned int n_fils, double delta_in, double ksp_in);
    void cl_init_start(double delta_in, double ksp_in);
    
    void neigh_append(unsigned int fil_id);
    
    void cl_walk_1_type_k(double fil_length_i, double motor_vel, double dt);
    void cl_walk_1_type_d(double fil_length_i, double motor_vel, double dt);
    void cl_walk_1_type_d1(double fil_length_i, double motor_vel, double dt);
    
    void cl_walk_2_type_k(double fil_length_i,double fil_length_j, double motor_vel, double dt);
    void cl_walk_2_type_d(double fil_length_i,double fil_length_j, double motor_vel, double dt);
    void cl_walk_2_type_d1(double fil_length_i,double fil_length_j, double motor_vel, double dt);
    

    void least_sep_calc_cl(vector<double>& ri_com,double cl_r[3],vector<double>& ui_vec,double L_i,int N_dims);

    
};

#endif