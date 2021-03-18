// Utilities library for useful calculations
// 
// Benjamin Dalton 18/11/2015

#ifndef binding_H
#define binding_H

#include "filament.h"
#include "crosslinker.h"

using namespace std;

class binding
{
    private:
    public:
    
    double delta_lambda_j = 0.0;
    int N_lambda_j = 0;
    
    double lambda_j_mp = 0.0;
    
    double v_mag = 0.0;
    double d_ij_mag = 0.0;
    
    double a = 0.0, b = 0.0, r = 0.0, r_sqr = 0.0, q = 0.0, w = 0.0, s = 0.0, t = 0.0;
    double d1 = 0.0, d2 = 0.0, d3 = 0.0, d4 = 0.0;
    double k_on = 0.0, eps_r = 0.0, e_eps_r = 0.0;
    
    double lambda_0 = 0.0;
    double lambda_min = 0.0;
    double lambda_accum = 0.0;
    
    unsigned int i,j,n;
    
    vector<double> rcli, v_cli, d_ij, d_0, r_image;
    
    binding();
    
    void binding_init(int n_dims, double eta_sites);

    double kon_calc(filament& mt_i, filament& mt_j,  vector<double>& rj_image, crosslinker& cl_k, int N_dim, double r_0);
    double nlist_cl(vector<double>& ri_com, double cl_r[3], vector<double>& ui_vec,   double L_i,int N_dims);
    double epsilon_j_calc(filament& mt_i, filament& mt_j, vector<double>& rj_image, crosslinker& cl_k, int N_dim, int m_id, double rand_in, double r_0);
    vector<double> image_calc_BC_cl(int BC_TYPE, int N_DIMS, vector<double>& ri_com, vector<double>& rj_com, vector<double>& box_vec);
    
};


#endif