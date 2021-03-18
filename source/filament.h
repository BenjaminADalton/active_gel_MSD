#ifndef filament_H
#define filament_H

#include <random>

using namespace std;

class filament
{
    
    
    private:
    public:

    unsigned int pop_id = 1;
    unsigned int fil_id = 0;
    
    unsigned int fil_time_nuc = 0;
    
    int growth_phase = 0;
    int growth_stall = 0; // communicate with other functions when growth is stalled due to upper length
    
    double fil_length;
    double scaled_length = 0.0;
    double fil_min = 2.0e-6;
    double gamma_para, gamma_perp, gamma_rot;
    double eqnm_pass_1,eqnm_pass_2,eqnm_pass_3;
    
    double rand_para;
    double rand_perp_1;
    double rand_perp_2;
    double rand_rot_1;
    double rand_rot_2;
    
    double pbc_index [3];

    double para_poly,perp_poly,rot_poly,log_ld,d_on_l,d_on_l_sqr;
    
    double msd_accum_para = 0;
    double msd_accum_perp = 0;
    double msd_accum_r = 0;
    
    filament();
    
    vector<vector<double>> u_perp;
    
    vector<double> r_com,r_eqnm,u_vec,u_eqnm,r_init,u_init,r_accum,r_hist;

    void filament_init(unsigned int n_dim);
    void filament_drag_tao(double d,double kappa);
    void filament_drag_torre(double d,double kappa);
    
};

#endif
