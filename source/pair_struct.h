// Filament class. Includes nucleation if a MT reduces to zero length.
//
// Benjamin Dalton 08/10/2015

#ifndef pair_struct_H
#define pair_struct_H

#define PI 3.141592653589793
#define PREC 5e-8

#include <vector>
#include <array>
#include <math.h>

using namespace std;

class pair_struct
{
    private:
    public:
    
    const double pi = 3.141592653589793;
    
    double w_ij_mag_single = 0.0;
    double abs_r_min_a = 0.0, abs_r_min_b = 0.0;
    double ri_min = 0.0, rj_min = 0.0;
    double u_ij = 0.0, ui_rij_dp = 0.0, uj_rij_dp = 0.0;
    double w_mag = 0.0, K = 0.0, la1_0 = 0.0, la2_0 = 0.0;
    double abs_uij = 0.0, abs_ui_rij_dp = 0.0, abs_uj_rij_dp = 0.0, abs_la1_0 = 0.0, abs_la2_0 = 0.0;
    double la1_a,la2_a,la1_b,la2_b;

    unsigned int n_pair, n_outer_pair, n_inner_pair;

    bool intercept_switch_1;
    bool intercept_switch_2;
    
    vector<double> r_ij, w_ij, r_min_temp;
    array<double,2> lambda;
    
    vector<double> w_ij_mag, F_ij_vec,lambda_single,uw_ij_single;
    vector< vector<double>> lambda_array,uw_ij;
    vector< vector<unsigned int>> pair_list, outer_nlist, inner_nlist;
    
    pair_struct();
    
    void pair_init(unsigned int n_dims);
    
    void pair_add_new(unsigned int n_fils, unsigned int n_dims);

    void nlist_init_state(unsigned int n_dims);
    void outer_nlist_init(unsigned int n_dims);
    void inner_nlist_init(unsigned int n_dims);
    void outer_nlist_build(unsigned int fil_i, unsigned int fil_j);
    void inner_nlist_build(unsigned int fil_i, unsigned int fil_j);
    
    void inner_pair_constructor(unsigned int n_dims);
    
    void least_sep_calc(vector<double>& ri_com,vector<double>& rj_com, vector<double>& ui_vec ,vector<double>& uj_vec,double L_i, double L_j,int N_dims,unsigned int pair_id);
    void least_sep_calc_single(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec ,vector<double>& uj_vec,double L_i, double L_j,int N_dims);
};

#endif