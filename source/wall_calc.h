// Filament class. Includes nucleation if a MT reduces to zero length.
//
// Benjamin Dalton 08/10/2015

#ifndef wall_calc_H
#define wall_calc_H

#define PI 3.141592653589793
#define PREC3 5e-7

#include <vector>
#include <array>
#include <math.h>

using namespace std;

class wall_calc
{
    
    private:
    public:
    
    double rix_out = 0.0;
    double w_ij_mag = 0.0;
    double w_ij_mag_out = 0.0;
    double u_ij = 0.0, ui_rij_dp = 0.0, uj_rij_dp = 0.0;
    double K = 0.0;
    double la1_0 = 0.0;
    double la2_0 = 0.0;
    double abs_uij = 0.0, abs_ui_rij_dp = 0.0, abs_uj_rij_dp = 0.0, abs_la1_0 = 0.0, abs_la2_0 = 0.0;
    double la1_a,la2_a,la1_b,la2_b;
    double abs_r_min_a = 0.0, abs_r_min_b = 0.0;

    double eps_on_sig = 0.0;
    double sig_on_wij = 0.0;
    
    double one_on_beta_wall;
    double sig_shift_wall;
    double F_i_wall = 0.0;
    double epsilon;
    double beta_sig_wall;
    double f_beta_wall = 0.0;
    double df_beta_wall = 0.0;
    
    vector<double> r_min_a , r_min_b;
    vector<double> w_ij, r_ij, ri_min, rj_min;
    
    wall_calc();
    
    array<double,2> lambda;
    
    vector<double> uw_ij_out;

    void wall_calc_init(int n_dims, double epsilon,double sig_shift_wall,double beta_sig_wall,double one_on_beta_wall);
    double wall_wca_calc();
    
    
    
    
    
    void wall_calc_class(vector<double>& ri_com,vector<double>& rj_com,vector<double>& ui_vec,vector<double>& uj_vec, double L_i, double L_j, int N_dims, int id);
    
};

#endif