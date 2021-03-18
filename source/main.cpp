
// Benjamin Dalton 04/10/2016

#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <string>
#include <ctime>
#include <typeinfo>
#include <algorithm>

#include "filament.h"
#include "crosslinker.h"
#include "utilities.h"
#include "binding.h"
#include "pair_struct.h"
#include "pbc.h"
#include "wall.h"
#include "wall_calc.h"

using namespace std;

#define T_STEP 20000000
#define N_DIMS 3            // Number of spatial dimensions (2 0r 3)
#define INIT 0              // initial configuration setup: 0 = new, 1 = read in configuration.txt
#define DRAG 0              // filament drag rule (0=Tao, 1=de la Torre)
#define COMP_FORCE 0        // compression force (used only for debugging)
#define TWEEZER 0           // compression force (used only for debugging)
#define BC_TYPE 4           // periodic boundar conditions (0,1,2,3 or 4 = see list below)
#define CL 1                // cross-linkers/motors (0 = off, 1 = on, 2 = just forces, skip binding/unbinding stochastics)
#define CL_INIT 0           // CL inital conditions (0 = all state 0, 1 = all state 1 and placed on MT) must be 0 if INIT = 1
#define CL_PLOT 0
#define SIJ_CALC 0          // mean square displacement and other correlations (0 = off, 1 = on)
#define DETATCH 2           // detatchment rules for motors at filament end points
#define CRYSTAL 0           // include fixed position filament structure (0 = off, 1 = on)

int main()
{
    ofstream log; log.open ("log.txt");
    
    time_t timev;
    time(&timev);
    
    // Set up clock properties:
    clock_t tic_full = clock();
    
    time_t now = time(0); tm *ltm = localtime(&now);
    
    int srand_int = 9 + (ltm->tm_hour)*(ltm->tm_min)*(ltm->tm_sec) + (ltm->tm_sec)*(ltm->tm_sec);
    
    char* time_date = ctime(&now); log << "Time and date of execution: " << time_date << endl;
    
    mt19937 generator (srand_int);
    mt19937 generator_util(2*srand_int+111);
    
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Random number generator seed
    cout << endl << "srand_int (random number seed): " << srand_int << endl;
    cout << "N_DIMS:    " << N_DIMS << endl;
    
    // number of filaments and crosslinkers
    int N_FILS    = 500;
    int N_CLS1    = 1000;
    int N_CLS2    = 0;
    int N_CLS_add = 0;
    
    // output for the MSD calculations
    int N_MSD_Calc = N_FILS;
    
    int N_CLS = N_CLS1 + N_CLS2;

    int CL_cryst_exclude = 1;       // prevent the 2nd CL population from interacting with CRYST
    int dyn_detach_type  = 0;       // 0=stay end bond, 1=detatch from end
    int dyn_single_loc   = 1;       // 0=anywhere, 1=only at end
    int dyn_single_walk  = 0;       // 0=both domains walk, 1=only j domain walks
    
    int N_CRYST          = 0;       // total number of pins
    int z_stack          = 0;       // must be integer divisable with N_CRYST
    int cryst_type       = 1;       // 0 = 2D, 1 = 3D
    double y_cryst       = 0.0e-6;  // place the c.o.m of pinning lattice in y

    // Used for outfile: number of CL's per lattice filament
    int cryst_nt = 0;
    int cryst_nt_array [1] = {0};
    
    // Choose the initial length distributions (0=uniform, 1=exponential)
    int L_dist = 0;
    double L_mean = 1.0e-6;
    double L_upper_limit = 3.0e-6;
    double L_lower_limit = 3.0e-7;
    
    // for the filament inital state setup
    double init_cutoff = 1.0;
    double init_length = 1.0e-6;
    double cryst_length = 20.0e-6;
    double bidisp_length = 1.0e-6;
    double upper_length = 200.0e-6;
    
    // box dimensions
    double L_x = 5.0e-6;
    double L_y = 5.0e-6;
    double L_z = 0.25e-6;
    
    // Kinesin-like or passive
    double k_0_3D    = 500; // (m^-1*s^-1*motor^-1)
    double k_off     = 0.033;
    
    // Dynein-like or passive
    double k_0_3D_2  = 500; // (m^-1*s^-1*motor^-1)
    double k_off_2   = 0.033;    // (s^-1)
    
    // Spring constants (3.0e-4N/m = 0.3pN/nm)
    double k_spr   = 1.0e-4;
    double k_spr_2 = 1.0e-4;

    // Stall force for kinesin and dynein
    double F_stall   = 5.0e-12;
    double F_stall_2 = 1.0e-12;
    
    double xz_nuc_frac = 1.0;
    double y_nuc_frac = 1.0;

    // velocity for the motor protein (0.0e-6 for passive)
    double motor_vel1 = 0.0e-6;
    double motor_vel2 = -0.0e-6;
    
    // (0) uy = 1, (1) random, (2) uy = -1, (3) 3D random,  (4) XWall random : note Cryst uy=1
    int rand_or_switch = 4;
    double rand_orient = 0.0;
    
    // depoly all filaments that turn in channel
    int uvec_stab = 0;
    
    // real time of each iteration
    const double dt = 5.0e-6;
    
    // print to standard out and txt files
    const unsigned int t_write = 10000;
    const unsigned int t_print = 500000;
    const unsigned int t_config    = 100000000;
    
    // time stepping parameters for the split time-scale trajectory outputs
    const unsigned int t_begin_write = 1000000;

    const unsigned int tscale_1_inc = 1;
    const unsigned int tscale_2_inc = 1000;
    const unsigned int tscale_3_inc = 100000;

    const unsigned int tscale_1_upper = 5000;   // tscale_1_upper/tscale_1_inc = n1_plot
    const unsigned int tscale_2_upper = 1000000; // tscale_2_upper/tscale_2_inc = n2_plot
    
    // trajectory_id files
    int t_write_2   = 10000000;
    int start_write = 150000000;
    
    // list update cycle lengths
    const unsigned int n_list_update       = 750;
    const unsigned int n_list_update_outer = 7500;
    const unsigned int wall_update         = 2000;
    
    // Distance cutoff for the inner and outer neighbour list method
    const double outer_nlist_cut = 2.5e-07;
    const double inner_nlist_cut = 0.75e-07;
    
    // time increments for the cross-linker event sampling
    const unsigned int cl_neigh_update      = 5000; // Must be the same as lowest of sample_inc_on1 and sample_inc_pair_on
    const unsigned int sample_inc_kn        = 5000;
    const unsigned int sample_inc_kcat      = 10000;
    const unsigned int sample_inc_on1       = 10000;
    const unsigned int sample_inc_off1      = 10000;
    const unsigned int sample_inc_pair_on   = 10000;
    const unsigned int sample_inc_pair_off  = 10000;
    
    // check filaments for the depolymerization state
    const unsigned int t_depoly = 5000;
    
    // parameters for fluid interactions - diffusion and drag
    const double k_B     = 1.381e-23;
    const double T       = 300;
    const double epsilon = k_B*T;
    const double kbt     = k_B*T;
    const double kTt_2   = 2.0*k_B*T*dt;
    const double pi      = 3.141592;
    const double eta     = 1.0;
    const double kappa   = pi*eta;
    const double lambda  = k_B*T/(2.0*pi*eta);

    // filament dimensions and WCA potential parameters
    const double d        = 2.5e-8; // Diameter, not sigma!!!!!
    double sig_shift      = 2.5e-8;
    double beta           = 0.9;
    double sig_shift_wall = 2.5e-8;
    double beta_wall      = 0.9;
    double sig_shift_min  = sig_shift*pow(2.0,1.0/6.0);
    double sig_min_wall   = sig_shift_wall*pow(2.0,1.0/6.0);
    
    double beta_sig         = beta*sig_shift;
    double one_on_beta      = sig_shift/beta_sig;
    double beta_sig_wall    = beta_wall*sig_shift_wall;
    double one_on_beta_wall = sig_shift_wall/beta_sig_wall;
    
    const double wall_cut  = 3.0;
    
    // cross-link/motor parameters
    double eta_sites   = 1.0*1.25e8;
    double k_0_2D      = 1.0e-18;         // (m^2/s)
    double sig_cl      = 0.3e-6;          // (m) min reach 2e-7m
    double k_spr_kbT   = k_spr/(k_B*T);
    double k_spr_kbT_2 = k_spr_2/(k_B*T);
    double cl_delta    = 1.38e-9;
    double cl_delta_2  = 1.38e-9;
    
    // use for spring force calculation
    double d_ij [3] = {0.0};
    double d_0_vec [3] = {0.0};
    double d_0         = 8.0e-8;
    double d_ij_mag    = 0.0;
    
    const double v_g = 0.25e-6; // 15um/min = 0.25e-6 m/s = 0.25e-6 um/us = 0.25 um/s
    const double v_s = 0.58e-6; // 35um/min

    double rand1 = 0.0;
    double k_on_dim = 0.0;
    
    // stall force u_vec projections
    double f_pi = 0.0;
    double f_pj = 0.0;
    
    // various parameters used and reused throughout
    int n, nt_cl1_1 = 0, nt_cl1_2 = 0, m_min=0,count=0;
    unsigned int i,j,k,l,m;
    double F_ij, F_i_wall, sig_wij, sig_wij_3, sig_wij_6, sig_wij_12, u_mod, umod, u_iw_dot, u_jw_dot;
    double sqrt2 = sqrt(2.0), sqrt2pi = sqrt(2.0*pi);
    
    // declare the filament and cross-linker objects
    vector<filament> mt(0);
    vector<crosslinker> cl(N_CLS+N_CLS_add);
    
    // basis vectors used in the stochastic dynamics
    double *e1_2d = new double[2]; e1_2d[0] = 1; e1_2d[1] = 0;
    double *e1_3d = new double[3]; e1_3d[0] = 1; e1_3d[1] = 0; e1_3d[2] = 0;
    double *e2_3d = new double[3]; e2_3d[0] = 0; e2_3d[1] = 1; e2_3d[2] = 0;
    
    // used for wall boundary conditions (don't need otherwise, but here they are anyway)
    vector<double> uwall1(3); uwall1[0] = 1.0; uwall1[1] = 0.0; uwall1[2] = 0.0;
    vector<double> uwall2(3); uwall2[0] = -1.0; uwall2[1] = 0.0; uwall2[2] = 0.0;
    vector<double> uwall3(3); uwall3[0] = 0.0; uwall3[1] = 0.0; uwall3[2] = 1.0;
    vector<double> uwall4(3); uwall4[0] = 0.0; uwall4[1] = 0.0; uwall4[2] = -1.0;
    
    // Use for the wall B.Cs 2D
    vector<double> ri_s2d(2);
    vector<double> uwall_2d(2);

    // Use for the wall B.Cs 3D
    vector<double> ri_s3d(3);
    vector<double> uwall_3d(3);
    double L_wall_2D = 0.0;
    double L_wall_3D = 20.0*L_y;
    
    double tot_fil_length = 0.0;
    double av_tot_fil_length = 0.0;
    double tot_scaled_length = 0.0;
    double cryst_fil_length = 0.0;
    double vol_fils  = 0.0;
    double vol_space = L_x*L_y*L_z;
    
    double gamma_array_2d[2][2]; double gamma_array_3d[3][3];
    
    double cl_check [3] = {0.0};
    
    // if(BC_TYPE == 1){pbc pbc1; pbc1.pbc_constructor(N_DIMS);}
    pbc pbc;
    pbc.pbc_constructor(N_DIMS, L_x, L_y, L_z);

    // declare wall objects
    wall wall;
    wall_calc wall_calc;
    wall_calc.wall_calc_init(N_DIMS,epsilon,sig_shift_wall,beta_sig_wall,one_on_beta_wall);

    binding binding;
    binding.binding_init(N_DIMS,eta_sites);

    double f_beta = 0.0;
    double df_beta = 0.0;

    double f_beta_wall = 0.0;
    double df_beta_wall = 0.0;
    
    double area = L_x*L_y; double volume = L_x*L_y*L_z;
    
    // initialise require class ans functions
    srand (srand_int);
    
    double u_id_dot = 0.0;
    double u_jd_dot = 0.0;
    
    // booking keeping tools for the cross-linker on-off events
    int cl_nst1 = 0;
    int cl_nst2 = 0;
    int *cl_st1_list = new int[N_CLS];   // list of all state 1 cls
    int *cl_st2_list = new int[N_CLS];   // list of all state 2 cls

    int filament_id_count = 0;
    
    for(i=0;i < N_CLS; i++)
    {
        cl[i].dyn_detach = dyn_detach_type; // 0=stay end bond, 1=detatch from end
        
        if(i < N_CLS1)
        {
            cl[i].type      = 1;
            cl[i].k_spr     = k_spr;
            cl[i].k_spr_kbT = k_spr_kbT;
            cl[i].cl_delta  = cl_delta;
            cl[i].k_on      = k_0_3D;
            cl[i].k_off     = k_off;
            cl[i].F_stall   = F_stall;
        }
        else if(i >= N_CLS1)
        {
            cl[i].type      = 2;
            cl[i].k_spr     = k_spr_2;
            cl[i].k_spr_kbT = k_spr_kbT_2;
            cl[i].k_on      = k_0_3D_2;
            cl[i].cl_delta  = cl_delta_2;
            cl[i].k_off     = k_off_2;
            cl[i].F_stall   = F_stall_2;
        }
    }
    
    // construct some vectors used throughout
    vector<double> box_vec,ri_image,rj_image,rcl_image;
    for(unsigned int n = 0; n < N_DIMS; n++)
    {
        ri_image.push_back(0.0);
        rj_image.push_back(0.0);
        rcl_image.push_back(0.0);
    }
    
    box_vec.push_back(L_x); box_vec.push_back(L_y); box_vec.push_back(L_z);
    
    // Create the order parameter tensor
    vector<double> buffer_S; vector<vector<double>> S_jk;
    for (j=0; j < N_DIMS; j++){buffer_S.push_back(0.0);}
    for (j=0; j < N_DIMS; j++){S_jk.push_back(buffer_S);}

    // Set up vectors for branching Gram Schmidt
    vector<double> u_perp_3D_1(3,0.0);
    vector<double> u_perp_3D_2(3,0.0);
    
    vector<double> network_fils;
    vector<double> network_fils_filtered;
    
    // constants used for the LJ-WCA Fij calculation
    f_beta  = -24.0*(epsilon/sig_shift)*(2.0*pow(one_on_beta,13.0) - pow(one_on_beta,7.0));
    df_beta = 24.0*(epsilon/sig_shift)*(26.0*pow(sig_shift,13.0)/pow(beta_sig,14.0) - 7.0*pow(sig_shift,7.0)/pow(beta_sig,8.0));
    
    // ************************* Set up pairs and create the inital state configuration *************************************************

    pair_struct pair;
    pair.pair_init(N_DIMS);
    pair.outer_nlist_init(N_DIMS);
    pair.inner_nlist_init(N_DIMS);
    pair.nlist_init_state(N_DIMS);

    if(INIT == 0) // Create new initial configuration
    {
        // create exponential length distribution and sort a list for longest first
        double *L_array = new double[N_FILS];
        if(L_dist == 1)
        {
            for(j=0; j < N_FILS; j++){L_array[j] = L_mean*std::log(1/dis(generator));}
            
            // need vector format for the sort (problem using vector above)
            std::vector<double> L_vector (L_array, L_array+N_FILS);
            std::sort (L_vector.begin(), L_vector.end());
            
            for(j=0; j < N_FILS; j++){L_array[j] = L_vector[N_FILS-j-1];}
            for(j=0; j < N_FILS; j++)
            {
                if(L_array[j] > L_upper_limit){L_array[j] = L_upper_limit;}
                if(L_array[j] < L_lower_limit){L_array[j] = L_lower_limit;}
            }
        }
        
        if(CRYSTAL == 1)
        {
            if(cryst_type == 0)
            {
                double shift_frac = 0.9;
                
                for(i=0; i < N_CRYST; i++)
                {
                    filament mt_new;
                    mt_new.filament_init(N_DIMS);
                    mt_new.fil_length = cryst_length;
                    
                    mt_new.r_com[0] = shift_frac*(L_x/(2.0*double(N_CRYST)) + double(i)*L_x/(double(N_CRYST)))+(L_x - shift_frac*L_x)/2.0;
                    
                    mt_new.r_com[1] = y_cryst;
                    
                    if(N_DIMS == 3){mt_new.r_com[2] = L_z/2.0;}
                    
                    mt_new.u_vec[1] = 1.0;
                    
                    mt_new.growth_phase = 3;
                    
                    mt.push_back(mt_new);
                    
                }
            }
            else if(cryst_type == 1)
            {
                for(k=0; k < z_stack; k++)
                {
                    for(i=0; i < N_CRYST/z_stack; i++)
                    {
                        filament mt_new;
                        mt_new.filament_init(N_DIMS);
                        mt_new.fil_length = cryst_length;
                        
                        mt_new.r_com[0] = (L_x/(2.0*double(N_CRYST/z_stack)) + double(z_stack*i)*L_x/(double(N_CRYST)));
                        
                        mt_new.r_com[1] = y_cryst;
                        
                        if(N_DIMS == 3){mt_new.r_com[2] = (k+1)*L_z/double(z_stack + 1);}
                        
                        mt_new.u_vec[1] = 1.0;
                        
                        mt_new.growth_phase = 3;
                        
                        mt.push_back(mt_new);
                        
                    }
                }
            }
        }
        
        int N_accum = N_CRYST;
        
        cout << endl << "Creating initial state:    " << endl;
        
        N_FILS = N_FILS + N_CRYST;
        
        int loop_int = 1;
        
        for(i = N_CRYST; i < N_FILS; i++)
        {
            filament mt_new;
            mt_new.filament_init(N_DIMS);
            if(L_dist == 0){mt_new.fil_length = init_length;}
            else if(L_dist == 1){mt_new.fil_length = L_array[i];}
            
            unsigned int nuc_block = 0;
            
            while(nuc_block == 0)
            {
                nuc_block = 1;
                loop_int = loop_int + 1;
                
                for (unsigned int n = 0; n < N_DIMS; n++)
                {
                    mt_new.r_com[n] = (box_vec[n]/xz_nuc_frac)*(uniform_random(i, i)-0.5) + box_vec[n]/2.0; mt_new.u_vec[n] = 0.0;
                }
                
                mt_new.r_com[1] = (box_vec[1]/y_nuc_frac)*(dis(generator)-0.5) + box_vec[1]/2.0;
                
                // choose either inital parallel or antiparallel initial orientations
                if(rand_or_switch == 0)
                {
                    rand_orient = 0.5;
                }
                else if(rand_or_switch == 1)
                {
                    rand_orient = dis(generator) - 0.5;
                }
                else if(rand_or_switch == 2)
                {
                    rand_orient = - 0.5;
                }
                
                if(rand_orient <= 0.0){mt_new.u_vec[1] = -1.0;}
                else if(rand_orient > 0.0){mt_new.u_vec[1] = 1.0;}
                
                if(rand_or_switch == 3)
                {
                    mt_new.u_vec[0] = dis(generator) - 0.5;
                    mt_new.u_vec[1] = dis(generator) - 0.5;
                    mt_new.u_vec[2] = dis(generator) - 0.5;
                    
                    u_mod = 0.0;
                    
                    for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt_new.u_vec[n]*mt_new.u_vec[n];}
                    u_mod = sqrt(u_mod);
                    for (unsigned int n = 0; n < N_DIMS; n++){mt_new.u_vec[n] = mt_new.u_vec[n]/u_mod;}
                    
                }
                if(rand_or_switch == 4)
                {
                    mt_new.u_vec[0] = dis(generator) - 0.5;
                    mt_new.u_vec[1] = dis(generator) - 0.5;
                    mt_new.u_vec[2] = 0.0;
                    
                    u_mod = 0.0;
                    
                    for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt_new.u_vec[n]*mt_new.u_vec[n];}
                    u_mod = sqrt(u_mod);
                    for (unsigned int n = 0; n < N_DIMS; n++){mt_new.u_vec[n] = mt_new.u_vec[n]/u_mod;}
                }
                
                if(N_accum > 0)
                {
                    for(unsigned int k = 0; k < N_accum; k++)
                    {
                        for(n=0; n<N_DIMS; n++){ri_image[n] = mt[k].r_com[n]; rj_image[n] = mt_new.r_com[n];}
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            if(mt[k].r_com[n] - mt_new.r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[k].r_com[n] -box_vec[n];}
                            else if(mt[k].r_com[n] - mt_new.r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[k].r_com[n] + box_vec[n];}
                        }
                        
                        pair.least_sep_calc_single(ri_image, rj_image, mt[k].u_vec,mt_new.u_vec, mt[k].fil_length,mt_new.fil_length,N_DIMS);
                        
                        if(pair.w_ij_mag_single <= init_cutoff*sig_shift_min){nuc_block = 0; break;}
                    }
                }
            }
            
            filament_id_count = filament_id_count + 1; mt_new.fil_id = filament_id_count;

            mt.push_back(mt_new);
            N_accum = N_accum + 1;
        }
    }
    
    u_mod = 0.0;
    for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[0].u_vec[n]*mt[0].u_vec[n];}
    u_mod = sqrt(u_mod);
    for (unsigned int n = 0; n < N_DIMS; n++){mt[0].u_vec[n] = mt[0].u_vec[n]/u_mod;}

    // ******************************************************************************************************************
    
    // initialise for the mean square displacement etc. correlations that require accumulations
    for(i=0;i<N_FILS;i++){for(n=0;n<N_DIMS;n++){mt[i].r_init[n] = mt[i].r_com[n]; mt[i].u_init[n] = mt[i].u_vec[n];};}
    
    // prepare the drag coefficients for each filament (update every time when polymerizing)
    if(DRAG == 0)
    {
        for (j=0; j < N_FILS; j++){mt[j].filament_drag_tao(d,kappa);}
    }
    else if(DRAG == 1)
    {
        for (j=0; j < N_FILS; j++){mt[j].filament_drag_torre(d,kappa);}
    }

    // open text files to be used for data writing
    ofstream trajectory;
    trajectory.open ("trajectory.txt");
    ofstream cl_trajectory;
    cl_trajectory.open ("trajectory_cl.txt");

    ofstream id_trajectory;
    id_trajectory.open ("trajectory_id.txt");
    ofstream id_orientation;
    id_orientation.open ("orientation_id.txt");
    ofstream network_trajectory;
    network_trajectory.open ("network_trajectory_id.txt");
    ofstream network_orientation;
    network_orientation.open ("network_orientation_id.txt");
    ofstream id_trajectory_com;
    id_trajectory_com.open ("trajectory_com_id.txt");
    ofstream minus_trajectory;
    minus_trajectory.open ("trajectory_plus.txt");
    ofstream pinning_stress;
    pinning_stress.open ("pinning_stress.txt");
    ofstream cl_nt;
    cl_nt.open ("trajectory_cl_nt.txt");
    ofstream cl_id;
    cl_id.open ("cl_pairs.txt");
    
    ofstream traj_tscale_1;
    traj_tscale_1.open ("traj_tscale_1.txt");
    ofstream traj_tscale_2;
    traj_tscale_2.open ("traj_tscale_2.txt");
    ofstream traj_tscale_3;
    traj_tscale_3.open ("traj_tscale_3.txt");
    
    ofstream fil_time;
    fil_time.open ("filament_lifetime.txt");
    ofstream occ_time;
    occ_time.open ("occ_time.txt");
    
    // printing to the log files.
    log << "Start parameters:" << endl;
    log << "\nPBC - type: " << BC_TYPE << "\nCLs - type: " << CL << endl;
    log << "\ntime steps:  " << T_STEP  << "\nt_write inc: " << t_write << "\ndt:          " << dt << endl;
    log << "\nN_DIMS: " << N_DIMS << "\nN_FILS: " << N_FILS  << "\nN_CL:  " << N_CLS << endl;
    log << "\nL_x: " << L_x  << "\nL_y: " << L_y << "\nL_z: " << L_z << endl;
    log << "\nInner update: " << n_list_update  << "\nOuter update: " << n_list_update_outer;
    log << "\nInner cut: " << inner_nlist_cut << "\nOuter cut: " << outer_nlist_cut << endl;

    log << "\nv_cl: " << motor_vel1 << endl<< endl;
    
    log << "CL parameters: " << endl;
    log << "eta_sites : " << eta_sites << endl;
    log << "k_off :     " << k_off << endl;
    log << "sig_cl :    " << sig_cl << endl;
    log << "k_spr :     " << k_spr << endl;
    log << "k_spr_kbT : " << k_spr_kbT << endl;
    log << "cl_delta :  " << cl_delta << endl;
    log << "F_stall :   " << F_stall << endl;
    log << "d_0 :       " << d_0 << endl << endl;
    
    for(j=0; j < N_FILS; j++){tot_fil_length = tot_fil_length + mt[j].fil_length;}
    for(j=0; j < N_CRYST; j++){cryst_fil_length = cryst_fil_length + mt[j].fil_length;}
    
    cout << "tot_fil_length    "  << tot_fil_length << endl;

    if(N_DIMS == 2){log << "P_0_2D : " << k_0_2D << endl;}
    if(N_DIMS == 2){log << "P_0_2D (system): " << (k_0_2D*eta_sites/area)*N_CLS*tot_fil_length*dt*sample_inc_on1 << endl;}
    if(N_DIMS == 2){log << "P_0_2D (per CL): " << (k_0_2D*eta_sites/area)*1*tot_fil_length*dt*sample_inc_on1 << endl << endl;}
    if(N_DIMS == 3){log << "P_0_3D : " << k_0_3D << endl;}
    if(N_DIMS == 3){log << "P_0_3D (system): " <<(k_0_3D*eta_sites/volume)*N_CLS*tot_fil_length*dt*sample_inc_on1 << endl;}
    if(N_DIMS == 3){log << "P_0_3D (per CL): " <<(k_0_3D*eta_sites/volume)*1*tot_fil_length*dt*sample_inc_on1 << endl << endl;}
    
    std::cout << "Start: number of CLs:   " << N_CLS << "    " << cl_nst1 << "    " << cl_nst2 << endl;
    
    unsigned int loop_int = 0;

    double length_acc = 0.0;
    
    // binding angle calculations
    double r_cl_ij_mag = 0.0;
    vector<double> r_cl_i(3,0.0);
    vector<double> r_cl_j(3,0.0);
    
    int CL0_loop, CL1_loop, CL2_loop;
    
    if(CL_cryst_exclude == 0) // 2nd population can interact with pinning lattice
    {
        CL0_loop = 0; CL1_loop = N_CLS; CL2_loop = 0;
    }
    else if(CL_cryst_exclude == 1) // prevent 2nd population from interacting with pinning lattice
    {
        CL0_loop = 0; CL1_loop = N_CLS1; CL2_loop = N_CLS;
    }

    cout << endl << "Start Main Loop:" << endl << endl;
    
    for(i=0;i<T_STEP;i++)
    {
        if(CL == 1)
        {
            // ********************  create state 1 bingings ********************
            //
            // We assume un-bound cross-links/motors are uniformly distributed in volume. Here we calculate
            // is an event occures where an unbound cross-linker/motor will be taken out of the environment
            // and located in state-1 somewhere on a filament. Location selection is from a uniform distribution
            // along the total length of all filaments. See cross-linker binding stochastics, 2nd paragraph
            
            if(i % sample_inc_on1 == 0)
            {
                // Initial filament binding
                tot_fil_length    = 0.0;
                
                av_tot_fil_length = double(N_FILS)*7.0e-6; // uniform weight for the dynien on-events
                for(j=0; j < N_FILS; j++){tot_fil_length = tot_fil_length + mt[j].fil_length;}
                
                // finding free location (0-state) to place new 1-state. Create MT neighbour list
                for(j=CL0_loop; j < CL1_loop; j++)
                {
                    if(cl[j].state == 0)
                    {
                        rand1 = dis(generator);
                        k_on_dim = 0.0;
                        
                        if(N_DIMS == 2){k_on_dim = (k_0_2D)*tot_fil_length*dt*sample_inc_on1;}
                        else if(N_DIMS == 3)
                        {
                            if(cl[j].type == 1){k_on_dim = (cl[j].k_on)*tot_fil_length*dt*sample_inc_on1;}
                            if(cl[j].type == 2){k_on_dim = (cl[j].k_on)*av_tot_fil_length*dt*sample_inc_on1;} // uniform weight for the dynien on-events
                        }
                        
                        // when binding will occur:
                        if(rand1 <= k_on_dim)
                        {
                            cl[j].n_neigh   = 0;
                            cl[j].nlist_clean();
                            cl[j].state     = 1;
                            
                            double rand2 = dis(generator);
                            double l_r2 = rand2*tot_fil_length;
                            double length_acc = 0;
                            
                            // calculate total filament location
                            if(cl[j].type == 1)
                            {
                                for(n=0; n < N_FILS; n++)
                                {
                                    length_acc = length_acc + mt[n].fil_length;
                                    
                                    if(l_r2 < length_acc)
                                    {
                                        cl[j].fil_i     = n;
                                        cl[j].epsilon_i = mt[n].fil_length/2.0 - (length_acc-l_r2);
                                        break;
                                    }
                                }
                            }
                            if(cl[j].type == 2)
                            {
                                for(n=0; n < N_FILS; n++)
                                {
                                 
                                    length_acc = length_acc + 7.0e-6; // uniform weight for the dynien on-events
                                    
                                    if(l_r2 < length_acc)
                                    {
                                        cl[j].fil_i     = n;
                                        cl[j].epsilon_i = -mt[cl[j].fil_i].fil_length/2.0;
                                        break;
                                    }
                                }
                            }
                            
                            cl[j].sing_time_on = i;
                            
                            cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                            
                            for(n=0; n<N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[j].fil_i].r_com[n] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[n];
                            }
                            
                            // once cl is placed in state-1 we need to create a list of all filaments with which it can potential interact
                            // and form a pair state binding
                            int cl_inc = 0;
                            
                            for (m=0; m < N_FILS; m++)
                            {
                                if(m != cl[j].fil_i)
                                {
                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[j].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[j].neigh_append(m);
                                        cl_inc = cl_inc + 1;
                                    }
                                }
                            }
                            cl[j].n_neigh = cl_inc;
                            
                            cl_nst1 = cl_nst1 + 1;
                        }
                    }
                }
                
                if(CL_cryst_exclude == 1) // only perform if 2nd population is to avoid pinning lattice
                {
                    for(j=CL1_loop; j < CL2_loop; j++)
                    {
                        if(cl[j].state == 0)
                        {
                            rand1 = dis(generator);
                            k_on_dim = 0.0;
                            
                            if(N_DIMS == 2){k_on_dim = (k_0_2D)*(tot_fil_length - cryst_fil_length)*dt*sample_inc_on1;}

                            else if(N_DIMS == 3){k_on_dim = (cl[j].k_on)*(double(N_FILS - N_CRYST)*7.0e-6)*dt*sample_inc_on1;}
                            
                            // when binding will occur:
                            if(rand1 <= k_on_dim)
                            {
                                cl[j].n_neigh   = 0;
                                cl[j].nlist_clean();
                                cl[j].state     = 1;
                                
                                double rand2 = dis(generator);

                                double l_r2 = rand2*double(N_FILS - N_CRYST)*7.0e-6;
                                double length_acc = 0;
                                
                                // calculate total filament location
                                for(n=N_CRYST; n < N_FILS; n++)
                                {
                                    length_acc = length_acc + 7.0e-6;
                                    
                                    if(l_r2 < length_acc)
                                    {
                                        cl[j].fil_i     = n;
                                        cl[j].epsilon_i = -mt[cl[j].fil_i].fil_length/2.0;
                                        break;
                                    }
                                }
                                
                                cl[j].sing_time_on = i;
                                
                                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                                
                                for(n=0; n<N_DIMS; n++)
                                {
                                    cl_check[n] = mt[cl[j].fil_i].r_com[n] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[n];
                                }
                                
                                // once cl is placed in state-1 we need to create a list of all filaments with which it can potential interact
                                // and form a pair state binding
                                int cl_inc = 0;
                                for (m=N_CRYST; m < N_FILS; m++)
                                {
                                    if(m != cl[j].fil_i)
                                    {

                                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[j].fil_i].r_com, box_vec);
                                        
                                        if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                        {
                                            cl[j].neigh_append(m);
                                            cl_inc = cl_inc + 1;
                                        }
                                    }
                                }
                                cl[j].n_neigh = cl_inc;
                                
                                cl_nst1 = cl_nst1 + 1;
                            }
                        }
                    }
                }
                
                // Update the cl1 and cl2 list and counts
                cl_nst1 = 0;
                cl_nst2 = 0;
                for(j=0; j < N_CLS; j++)
                {
                    if(cl[j].state == 1)
                    {
                        cl_st1_list[cl_nst1] = j;
                        cl_nst1 = cl_nst1 + 1;
                    }
                    else if(cl[j].state == 2)
                    {
                        cl_st2_list[cl_nst2] = j;
                        cl_nst2 = cl_nst2 + 1;
                    }
                }
            }
            
            if(i % sample_inc_off1 == 0)
            {
                // ******************** Un-binding from single filament ********************
                // loop though all state-1 filaments and sample whether it unbinds back into free diffusing state
                for(j=0; j < cl_nst1; j++)
                {
                    if(dis(generator) <= cl[cl_st1_list[j]].k_off*dt*sample_inc_off1)
                    {
                        cl[cl_st1_list[j]].state = 0;
                        cl[cl_st1_list[j]].nlist_clean();
                        
                    }
                }
                int cl_nst1 = 0;
                for(j=0; j < N_CLS; j++)
                {
                    if(cl[j].state == 1)
                    {
                        cl_st1_list[cl_nst1] = j;
                        cl_nst1 = cl_nst1 + 1;
                    }
                }
            }
            
            if(i % cl_neigh_update == 0)
            {
                for (j=0; j < cl_nst1; j++)
                {
                    int cl_inc = 0;
                    
                    cl[cl_st1_list[j]].nlist_clean();
                    
                    cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                    
                    for(n=0; n<N_DIMS; n++)
                    {
                        cl_check[n] = mt[cl[cl_st1_list[j]].fil_i].r_com[n] + cl[cl_st1_list[j]].epsilon_i*mt[cl[cl_st1_list[j]].fil_i].u_vec[n];
                    }
                    
                    // to exclude crystal for dynein
                    m_min = 0; if(CL_cryst_exclude == 1 && cl[cl_st1_list[j]].type == 2){m_min = N_CRYST;}
                    
                    for (m = m_min; m < N_FILS; m++)
                    {
                        // check all filament except the one that it is on
                        if(m != cl[cl_st1_list[j]].fil_i)
                        {
                            rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                            
                            if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                            {
                                cl[cl_st1_list[j]].neigh_append(m);
                                cl_inc = cl_inc + 1;
                            }
                        }
                    }
                    
                    cl[cl_st1_list[j]].n_neigh = cl_inc;
                }
            }
            
            //******************* create a state-2 binding ******************
            if(i % sample_inc_pair_on == 0)
            {
                for(j=0; j < cl_nst1; j++)
                {
                    double k_on = 0.0;
                    
                    // Make loop over all filaments that are close that the CL can extend to reach it
                    // cl[X].n_neigh list holds these ID's for state-1 cross-linkers.
                    for (m=0; m < cl[cl_st1_list[j]].n_neigh; m++)
                    {
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st1_list[j]].cl_neigh[m]].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                        
                        cl[cl_st1_list[j]].cl_kon_i[m] = binding.kon_calc(mt[cl[cl_st1_list[j]].fil_i], mt[cl[cl_st1_list[j]].cl_neigh[m]], rcl_image, cl[cl_st1_list[j]], N_DIMS, d_0);
                    
                        k_on = k_on + cl[cl_st1_list[j]].cl_kon_i[m];
                    }
                    
                    double rand_2state = dis(generator);
                    
                    // Use kon to determine if pair-binding will occur for cross-linker j (in state-1). See Eq 4 and 6 in document.
                    if(rand_2state < cl[cl_st1_list[j]].k_off*dt*eta_sites*sample_inc_pair_on*k_on)
                    {
                        // event will occur. use stored cl_kon_i for that cross-linker to find which filament form it's neighbor set
                        // will be involved in the pair.
                        double k_rand = dis(generator)*k_on;
                        double k_accum = 0.0;
                        
                        for(m=0; m < cl[cl_st1_list[j]].n_neigh; m++)
                        {
                            k_accum = k_accum + cl[cl_st1_list[j]].cl_kon_i[m];
                            
                            // The pair filament has been selected, now find the location on the filament to place second
                            // cross-linker head
                            if(k_rand < k_accum)
                            {
                                
                                cl[cl_st1_list[j]].state = 2;
                                cl[cl_st1_list[j]].occ_time_on = i;
                                cl[cl_st1_list[j]].fil_j = cl[cl_st1_list[j]].cl_neigh[m];
                                
                                // create cross-linker position vector (need to change later, inefficient double calculations)
                                
                                cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                                
                                for(n=0; n<N_DIMS; n++)
                                {
                                    cl_check[n] = mt[cl[cl_st1_list[j]].fil_i].r_com[n] + cl[cl_st1_list[j]].epsilon_i*mt[cl[cl_st1_list[j]].fil_i].u_vec[n];
                                }
                                
                                rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st1_list[j]].fil_j].r_com, mt[cl[cl_st1_list[j]].fil_i].r_com, box_vec);
                                
                                cl[cl_st1_list[j]].epsilon_j = binding.epsilon_j_calc(mt[cl[cl_st1_list[j]].fil_i], mt[cl[cl_st1_list[j]].fil_j], rcl_image, cl[cl_st1_list[j]], N_DIMS, m, dis(generator),d_0);
                                
                                cl[cl_st1_list[j]].nlist_clean();
                                
                                break;
                            }
                        }
                    }
                    
                } // for loop close
                
                cl_nst1 = 0;
                cl_nst2 = 0;
                
                for(m=0; m < N_CLS; m++)
                {
                    if(cl[m].state == 1)
                    {
                        cl_st1_list[cl_nst1] = m;
                        cl_nst1 = cl_nst1 + 1;
                    }
                    if(cl[m].state == 2)
                    {
                        cl_st2_list[cl_nst2] = m;
                        cl_nst2 = cl_nst2 + 1;
                    }
                }
            }
            
            //******************* remove a state-2 binding ******************
            if(i % sample_inc_pair_off == 0)
            {
                if(cl_nst2 > 0)
                {
                    for(j=0; j < cl_nst2; j++)
                    {
                        
                        rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st2_list[j]].fil_i].r_com, mt[cl[cl_st2_list[j]].fil_j].r_com, box_vec);
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            d_ij[n] = (rcl_image[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n]) - (mt[cl[cl_st2_list[j]].fil_j].r_com[n] + cl[cl_st2_list[j]].epsilon_j*mt[cl[cl_st2_list[j]].fil_j].u_vec[n]);
                        }
                        
                        d_ij_mag = 0.0;
                        
                        for(n=0; n<N_DIMS; n++)
                        {
                            d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
                        }
                        
                        d_ij_mag = sqrt(d_ij_mag);

                        double rand_2_1 = dis(generator);
                        
                        if(rand_2_1 < dt*sample_inc_pair_off*cl[cl_st2_list[j]].k_off*exp(cl[cl_st2_list[j]].k_spr_kbT*cl[cl_st2_list[j]].cl_delta*(d_ij_mag - d_0)))
                        {
                            
                            cl[cl_st2_list[j]].state = 1;
                            cl[cl_st2_list[j]].nlist_clean();
                            
                            occ_time << i - cl[cl_st2_list[j]].occ_time_on << endl;
                            cl[cl_st2_list[j]].occ_time_on = 0;
                            
                            cl[cl_st2_list[j]].sing_time_on = i;
                            
                            double bin_rand = dis(generator) - 0.5;
                            
                            if(dyn_single_loc == 1 && cl[cl_st2_list[j]].type==2){bin_rand = -1.0;}
                            
                            if(bin_rand <= 0.0)
                            {
                                cl[cl_st2_list[j]].epsilon_j = 0.0;
                                cl[cl_st2_list[j]].fil_j = 0;
                                cl[cl_st2_list[j]].state = 1;
                                cl[cl_st2_list[j]].nlist_clean();
                                cl[cl_st2_list[j]].end_state_j = 0;

                            }
                            else if(bin_rand > 0.0)
                            {
                                cl[cl_st2_list[j]].epsilon_i = cl[cl_st2_list[j]].epsilon_j;
                                cl[cl_st2_list[j]].fil_i = cl[cl_st2_list[j]].fil_j;


                                cl[cl_st2_list[j]].epsilon_j = 0.0;
                                cl[cl_st2_list[j]].fil_j = 0;
                                cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                                cl[cl_st2_list[j]].nlist_clean();
                                cl[cl_st2_list[j]].end_state_i = 0;
                                cl[cl_st2_list[j]].end_state_j = 0;
                            }
                            
                            cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                            
                            for(n=0; n < N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];
                            }
                            
                            // re-create nlist for the cl returned to state 1
                            int cl_inc_off = 0;
                            
                            m_min = 0;
                            
                            if(CL_cryst_exclude == 1 && cl[cl_st2_list[j]].type == 2){m_min = N_CRYST;}
                            
                            for (m = m_min; m < N_FILS; m++)
                            {
                                if(m != cl[cl_st2_list[j]].fil_i)
                                {

                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[cl_st2_list[j]].neigh_append(m);
                                        cl_inc_off = cl_inc_off + 1;
                                    }
                                }
                            }
                            cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                        }
                    }
                    
                    cl_nst1 = 0;
                    cl_nst2 = 0;
                    
                    for(m=0; m < N_CLS; m++)
                    {
                        if(cl[m].state == 1)
                        {
                            cl_st1_list[cl_nst1] = m;
                            cl_nst1 = cl_nst1 + 1;
                        }
                        if(cl[m].state == 2)
                        {
                            cl_st2_list[cl_nst2] = m;
                            cl_nst2 = cl_nst2 + 1;
                        }
                    }
                }
            }
            
            if(DETATCH == 2) // detatch both (fix_here1)
            {
                for (j=0; j < cl_nst1; j++)
                {
                    
                    if(cl[cl_st1_list[j]].type == 1){cl[cl_st1_list[j]].cl_walk_1_type_k(mt[cl[cl_st1_list[j]].fil_i].fil_length, motor_vel1, dt);}
                    if(cl[cl_st1_list[j]].type == 2 && dyn_single_walk == 0){cl[cl_st1_list[j]].cl_walk_1_type_d(mt[cl[cl_st1_list[j]].fil_i].fil_length, motor_vel2, dt);}
                    if(cl[cl_st1_list[j]].type == 2 && dyn_single_walk == 1){cl[cl_st1_list[j]].cl_walk_1_type_d1(mt[cl[cl_st1_list[j]].fil_i].fil_length, motor_vel2, dt);}
                    
                    if(cl[cl_st1_list[j]].end_state == 1)
                    {
                        cl[cl_st1_list[j]].state = 0;
                        cl[cl_st1_list[j]].nlist_clean();
                        cl[cl_st1_list[j]].end_state = 0;
                    }
                    
                }
                for (j=0; j < cl_nst2; j++)
                {
                    if(cl[cl_st2_list[j]].type == 1){cl[cl_st2_list[j]].cl_walk_2_type_k(mt[cl[cl_st2_list[j]].fil_i].fil_length,mt[cl[cl_st2_list[j]].fil_j].fil_length, motor_vel1, dt);}
                    if(cl[cl_st2_list[j]].type == 2 && dyn_single_walk == 0){cl[cl_st2_list[j]].cl_walk_2_type_d(mt[cl[cl_st2_list[j]].fil_i].fil_length,mt[cl[cl_st2_list[j]].fil_j].fil_length, motor_vel2, dt);}
                    if(cl[cl_st2_list[j]].type == 2 && dyn_single_walk == 1){cl[cl_st2_list[j]].cl_walk_2_type_d1(mt[cl[cl_st2_list[j]].fil_i].fil_length,mt[cl[cl_st2_list[j]].fil_j].fil_length, motor_vel2, dt);}
                }
                
                for (j=0; j < cl_nst2; j++)
                {
                    if(cl[cl_st2_list[j]].end_state_i == 1)
                    {
                        cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                        cl[cl_st2_list[j]].nlist_clean();
                        cl[cl_st2_list[j]].epsilon_i = cl[cl_st2_list[j]].epsilon_j;
                        cl[cl_st2_list[j]].fil_i = cl[cl_st2_list[j]].fil_j;
                        cl[cl_st2_list[j]].epsilon_j = 0.0;
                        cl[cl_st2_list[j]].fil_j = 0;
                        cl[cl_st2_list[j]].end_state_i = 0;
                        cl[cl_st2_list[j]].end_state_j = 0;
                        
                        occ_time << i - cl[cl_st2_list[j]].occ_time_on << endl;
                        cl[cl_st2_list[j]].occ_time_on = 0;
                        
                        if(CL == 1) // Note (999):
                        {
                            cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                            
                            for(n=0; n<N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];
                            }
                            
                            // re - create nlist for the cl returned to state 1
                            int cl_inc_off = 0;
                            
                            m_min = 0; if(CL_cryst_exclude == 1 && cl[cl_st2_list[j]].type == 2){m_min = N_CRYST;}
                            
                            for (m = m_min; m < N_FILS; m++)
                            {
                                if(m != cl[cl_st2_list[j]].fil_i)
                                {
                                    
                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[cl_st2_list[j]].neigh_append(m);
                                        cl_inc_off = cl_inc_off + 1;
                                    }
                                }
                            }
                            
                            cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                        }
                    }
                    if( cl[cl_st2_list[j]].end_state_j == 1)
                    {
                        cl[cl_st2_list[j]].state = 1; // NOTE: need to recalculate the neighbour list for this state-1 cl
                        cl[cl_st2_list[j]].nlist_clean();
                        cl[cl_st2_list[j]].epsilon_j = 0.0;
                        cl[cl_st2_list[j]].fil_j = 0;
                        cl[cl_st2_list[j]].end_state_j = 0;
                        
                        occ_time << i - cl[cl_st2_list[j]].occ_time_on << endl;
                        cl[cl_st2_list[j]].occ_time_on = 0;
                        
                        cl_check [0] = 0.0; cl_check [1] = 0.0; cl_check [2] = 0.0;
                        
                        if(CL == 1) // Note (999):
                        {
                            
                            for(n=0; n<N_DIMS; n++)
                            {
                                cl_check[n] = mt[cl[cl_st2_list[j]].fil_i].r_com[n] + cl[cl_st2_list[j]].epsilon_i*mt[cl[cl_st2_list[j]].fil_i].u_vec[n];
                            }
                            
                            // re - create nlist for the cl returned to state 1
                            int cl_inc_off = 0;
                            
                            m_min = 0; if(CL_cryst_exclude == 1 && cl[cl_st2_list[j]].type == 2){m_min = N_CRYST;}
                            
                            for (m = m_min; m < N_FILS; m++)
                            {
                                if(m != cl[cl_st2_list[j]].fil_i)
                                {
                                    
                                    rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[m].r_com, mt[cl[cl_st2_list[j]].fil_i].r_com, box_vec);
                                    
                                    if(binding.nlist_cl(rcl_image,cl_check,mt[m].u_vec,mt[m].fil_length,N_DIMS) <= sig_cl)
                                    {
                                        cl[cl_st2_list[j]].neigh_append(m);
                                        cl_inc_off = cl_inc_off + 1;
                                    }
                                }
                            }
                            cl[cl_st2_list[j]].n_neigh = cl_inc_off;
                            
                        }
                    }
                }
            }

            // Rebuild the cl state lists (Need to change from every time step)
            cl_nst1 = 0;
            cl_nst2 = 0;
            for(j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 1)
                {
                    cl_st1_list[cl_nst1] = j;
                    cl_nst1 = cl_nst1 + 1;
                }
                else if(cl[j].state == 2)
                {
                    cl_st2_list[cl_nst2] = j;
                    cl_nst2 = cl_nst2 + 1;
                }
            }
        }
        
        // numerically update the Eqns of Motion (2D or 3D). Include thermal noise contribution (Brownian-force)
        if(N_DIMS == 2)
        {
            for (j=N_CRYST; j < N_FILS; j++)
            {
                gram_schmidt_2d(mt[j].u_vec, mt[j].u_perp, e1_2d);
                stoch_2d_util(mt[j], kTt_2);

                u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
                u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}
            }
        }
        else if(N_DIMS == 3)
        {
            for (j = N_CRYST; j < N_FILS; j++)
            {
                gram_schmidt_3d(mt[j].u_vec, mt[j].u_perp, e1_3d, e2_3d);
                stoch_3d_util(mt[j], kTt_2);
                
                u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
                u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}
            }
            
        }
        
        // re-initialise the equations of motion
        for (j=N_CRYST; j < N_FILS; j++){for(n=0; n<N_DIMS; n++){mt[j].r_eqnm[n] = 0.0;} mt[j].u_eqnm[0] = 0.0; mt[j].u_eqnm[1] = 0.0; mt[j].u_eqnm[2] = 0.0;}
        
        if(BC_TYPE == 0) // No PBC and no Walls - infinite space
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        pair.least_sep_calc_single(mt[i].r_com, mt[j].r_com, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(i,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    pair.least_sep_calc_single(mt[pair.outer_nlist[k][0]].r_com, mt[pair.outer_nlist[k][1]].r_com, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                pair.inner_pair_constructor(N_DIMS);
            }
            for(unsigned int k = 0; k < pair.n_inner_pair; k++)
            {
                pair.least_sep_calc(mt[pair.inner_nlist[k][0]].r_com,mt[pair.inner_nlist[k][1]].r_com, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
            }
        }
        if(BC_TYPE == 1) // PBC in (x, y) for 2D and (x, y, z) for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[i].r_com[n];
                            rj_image[n] = mt[j].r_com[n];
                            
                            if(mt[i].r_com[n] - mt[j].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n]-box_vec[n];}
                            else if(mt[i].r_com[n] - mt[j].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n] + box_vec[n];}
                        }
                        pair.least_sep_calc_single(ri_image, rj_image, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(i,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    for(n=0;n<N_DIMS;n++)
                    {
                        ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.outer_nlist[k][1]].r_com[n];
                    }
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]-box_vec[n];}
                        else if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                
                // pair.inner_pair_constructor(N_DIMS);
                
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    for(unsigned int n = 0; n < N_DIMS; n++)
                    {
                        ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.inner_nlist[k][1]].r_com[n];
                        
                        if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] - box_vec[n];}
                        else if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] < - box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    pair.least_sep_calc(ri_image,rj_image, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        if(BC_TYPE == 2) // Walls in x with PBC in y for 2D and (x, z) with PBC in y for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS); pair.inner_nlist_init(N_DIMS);
                for (unsigned int k = 0; k < N_FILS-1; k++)
                {
                    for(unsigned int j = k+1; j < N_FILS; j++)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[k].r_com[n]; rj_image[n] = mt[j].r_com[n];
                        }
                        
                        if(mt[k].r_com[1] - mt[j].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[k].r_com[1] - box_vec[1];}
                        else if(mt[k].r_com[1] - mt[j].r_com[1] < -box_vec[1]/2.0){ri_image[1] = mt[k].r_com[1] + box_vec[1];}
                        
                        pair.least_sep_calc_single(ri_image, rj_image, mt[k].u_vec ,mt[j].u_vec, mt[k].fil_length,mt[j].fil_length,N_DIMS);
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(k,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    for(n=0;n<N_DIMS;n++)
                    {
                        ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.outer_nlist[k][1]].r_com[n];
                    }
                    
                    if(mt[pair.outer_nlist[k][0]].r_com[1] - mt[pair.outer_nlist[k][1]].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[pair.outer_nlist[k][0]].r_com[1]-box_vec[1];}
                    else if(mt[pair.outer_nlist[k][0]].r_com[1] - mt[pair.outer_nlist[k][1]].r_com[1] < -box_vec[1]/2.0){ri_image[1] = mt[pair.outer_nlist[k][0]].r_com[1] + box_vec[1];}
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    for(unsigned int n = 0; n < N_DIMS; n++)
                    {
                        ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n]; rj_image[n] = mt[pair.inner_nlist[k][1]].r_com[n];
                    }
                    
                    if(mt[pair.inner_nlist[k][0]].r_com[1] - mt[pair.inner_nlist[k][1]].r_com[1] > box_vec[1]/2.0){ri_image[1] = mt[pair.inner_nlist[k][0]].r_com[1] - box_vec[1];}
                    else if(mt[pair.inner_nlist[k][0]].r_com[1] - mt[pair.inner_nlist[k][1]].r_com[1] < - box_vec[1]/2.0){ri_image[1] = mt[pair.inner_nlist[k][0]].r_com[1] + box_vec[1];}
                    
                    pair.least_sep_calc(ri_image,rj_image, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
            
        }
        if(BC_TYPE == 3) // Walls in x with PBC in y for 2D and (x, z) open BC in y for 3D
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                
                for (unsigned int k = 0; k < N_FILS-1; k++)
                {
                    for(unsigned int j = k+1; j < N_FILS; j++)
                    {
                        
                        pair.least_sep_calc_single(mt[k].r_com, mt[j].r_com, mt[k].u_vec ,mt[j].u_vec, mt[k].fil_length,mt[j].fil_length,N_DIMS);
                        
                        if(k < N_CRYST && j < N_CRYST){pair.w_ij_mag_single = 1.0;}
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(k,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    pair.least_sep_calc_single(mt[pair.outer_nlist[k][0]].r_com, mt[pair.outer_nlist[k][1]].r_com, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                pair.inner_pair_constructor(N_DIMS);
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    pair.least_sep_calc(mt[pair.inner_nlist[k][0]].r_com,mt[pair.inner_nlist[k][1]].r_com, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        if(BC_TYPE == 4) // PBC in (x,y) with wall in z. USE FOR 3D ONLY!!!!
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS); pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            ri_image[n] = mt[i].r_com[n];
                            rj_image[n] = mt[j].r_com[n];
                        }
                        for(n=0; n<2; n++)
                        {
                            if(mt[i].r_com[n] - mt[j].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n]-box_vec[n];}
                            else if(mt[i].r_com[n] - mt[j].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[i].r_com[n] + box_vec[n];}
                        }
                        
                        pair.least_sep_calc_single(ri_image, rj_image, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        
                        if(pair.w_ij_mag_single <= outer_nlist_cut){pair.outer_nlist_build(i,j);}
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    for(n=0;n<N_DIMS;n++)
                    {
                        ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n];
                        rj_image[n] = mt[pair.outer_nlist[k][1]].r_com[n];
                    }
                    
                    for(n=0; n < 2; n++)
                    {
                        if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n]-box_vec[n];}
                        else if(mt[pair.outer_nlist[k][0]].r_com[n] - mt[pair.outer_nlist[k][1]].r_com[n] < -box_vec[n]/2.0){ri_image[n] = mt[pair.outer_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    
                    pair.least_sep_calc_single(ri_image, rj_image, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
            }
            {
                for(unsigned int k = 0; k < pair.n_inner_pair; k++)
                {
                    for(unsigned int n = 0; n < N_DIMS; n++)
                    {
                        ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n];
                        rj_image[n] = mt[pair.inner_nlist[k][1]].r_com[n];
                    }
                    for(unsigned int n = 0; n < 2; n++)
                    {
                        if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] > box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] - box_vec[n];}
                        else if(mt[pair.inner_nlist[k][0]].r_com[n] - mt[pair.inner_nlist[k][1]].r_com[n] < - box_vec[n]/2.0){ri_image[n] = mt[pair.inner_nlist[k][0]].r_com[n] + box_vec[n];}
                    }
                    pair.least_sep_calc(ri_image,rj_image, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
                }
            }
        }
        if(BC_TYPE == 5) // Open in (x,y) with wall in z. USE FOR 3D ONLY!!!!
        {
            if (i % n_list_update_outer == 0)
            {
                // Initialise the outer neighbour list
                pair.outer_nlist_init(N_DIMS);
                pair.inner_nlist_init(N_DIMS);
                for (unsigned int i = 0; i < N_FILS-1; i++)
                {
                    for(unsigned int j = i+1; j < N_FILS; j++)
                    {
                        
                        pair.least_sep_calc_single(mt[i].r_com, mt[j].r_com, mt[i].u_vec ,mt[j].u_vec, mt[i].fil_length,mt[j].fil_length,N_DIMS);
                        
                        if(pair.w_ij_mag_single <= outer_nlist_cut)
                        {
                            pair.outer_nlist_build(i,j);
                        }
                    }
                }
            }
            if (i % n_list_update == 0)
            {
                pair.inner_nlist_init(N_DIMS);
                
                for(unsigned int k = 0; k < pair.n_outer_pair; k++)
                {
                    pair.least_sep_calc_single(mt[pair.outer_nlist[k][0]].r_com, mt[pair.outer_nlist[k][1]].r_com, mt[pair.outer_nlist[k][0]].u_vec ,mt[pair.outer_nlist[k][1]].u_vec, mt[pair.outer_nlist[k][0]].fil_length,mt[pair.outer_nlist[k][1]].fil_length,N_DIMS);
                    if(pair.w_ij_mag_single <= inner_nlist_cut){pair.inner_nlist_build(pair.outer_nlist[k][0],pair.outer_nlist[k][1]);}
                }
                pair.inner_pair_constructor(N_DIMS);
            }
            for(unsigned int k = 0; k < pair.n_inner_pair; k++)
            {
                pair.least_sep_calc(mt[pair.inner_nlist[k][0]].r_com,mt[pair.inner_nlist[k][1]].r_com, mt[pair.inner_nlist[k][0]].u_vec ,mt[pair.inner_nlist[k][1]].u_vec, mt[pair.inner_nlist[k][0]].fil_length,mt[pair.inner_nlist[k][1]].fil_length,N_DIMS,k);
            }
        }
        
        //************************ Start force section of main loop **************************
        
        // Calculate force from "WCA potential" along all lines of least of interaction for the inner pair list
        for(unsigned int k = 0; k < pair.n_inner_pair; k++)
        {

            F_ij = 0.0;
            
            if(pair.w_ij_mag[k] <= sig_shift_min)
            {
                if(pair.w_ij_mag[k] > beta_sig)
                {
                    F_ij = -24.0*(epsilon/sig_shift)*(2.0*pow(sig_shift/pair.w_ij_mag[k],13.0) - pow(sig_shift/pair.w_ij_mag[k],7.0));
                }
                else if(pair.w_ij_mag[k] <= beta_sig)
                {
                    F_ij    = f_beta + df_beta*(pair.w_ij_mag[k] - beta_sig);
                }
                if(pair.w_ij_mag[k] < 1.0e-13)
                {
                    std :: cout << "Filaments intercepting: " << i << "    " << pair.w_ij_mag[k] << "\t" << pair.inner_nlist[k][0] << "\t" << pair.inner_nlist[k][1] << endl;
                    std :: cout << "Filament lengths:     " << mt[pair.inner_nlist[k][0]].fil_length << "\t" << mt[pair.inner_nlist[k][1]].fil_length << endl;
                    
                    if(N_DIMS == 2){abort();}
                    else if(N_DIMS == 3){F_ij = 0;}
                }
            }
            
            pair.F_ij_vec[k] = F_ij;
        }
        
        if(N_DIMS == 2)
        {
            // calculate steric interaction force for Eqns of Mot. needed for both position and orientation equations (steric-force)
            for(unsigned int k = 0; k < pair.n_inner_pair; k++)
            {
                
                mt[pair.inner_nlist[k][0]].r_eqnm[0] = mt[pair.inner_nlist[k][0]].r_eqnm[0] - pair.F_ij_vec[k]*pair.uw_ij[k][0];
                mt[pair.inner_nlist[k][1]].r_eqnm[0] = mt[pair.inner_nlist[k][1]].r_eqnm[0] + pair.F_ij_vec[k]*pair.uw_ij[k][0];
                
                mt[pair.inner_nlist[k][0]].r_eqnm[1] = mt[pair.inner_nlist[k][0]].r_eqnm[1] - pair.F_ij_vec[k]*pair.uw_ij[k][1];
                mt[pair.inner_nlist[k][1]].r_eqnm[1] = mt[pair.inner_nlist[k][1]].r_eqnm[1] + pair.F_ij_vec[k]*pair.uw_ij[k][1];
                
                u_iw_dot = mt[pair.inner_nlist[k][0]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][0]].u_vec[1]*pair.uw_ij[k][1];
                u_jw_dot = mt[pair.inner_nlist[k][1]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][1]].u_vec[1]*pair.uw_ij[k][1];
                
                
                mt[pair.inner_nlist[k][0]].u_eqnm[0] = mt[pair.inner_nlist[k][0]].u_eqnm[0] + pair.F_ij_vec[k]*pair.lambda_array[k][0]*(u_iw_dot*mt[pair.inner_nlist[k][0]].u_vec[0] - pair.uw_ij[k][0]);
                mt[pair.inner_nlist[k][1]].u_eqnm[0] = mt[pair.inner_nlist[k][1]].u_eqnm[0] - pair.F_ij_vec[k]*pair.lambda_array[k][1]*(u_jw_dot*mt[pair.inner_nlist[k][1]].u_vec[0] - pair.uw_ij[k][0]);
                
                mt[pair.inner_nlist[k][0]].u_eqnm[1] = mt[pair.inner_nlist[k][0]].u_eqnm[1] + pair.F_ij_vec[k]*pair.lambda_array[k][0]*(u_iw_dot*mt[pair.inner_nlist[k][0]].u_vec[1] - pair.uw_ij[k][1]);
                mt[pair.inner_nlist[k][1]].u_eqnm[1] = mt[pair.inner_nlist[k][1]].u_eqnm[1] - pair.F_ij_vec[k]*pair.lambda_array[k][1]*(u_jw_dot*mt[pair.inner_nlist[k][1]].u_vec[1] - pair.uw_ij[k][1]);
                
            }
        }
        else if(N_DIMS == 3)
        {
            // calculate steric interaction force for Eqns of Mot. needed for both position and orientation equations (steric-force)
            for(unsigned int k = 0; k < pair.n_inner_pair; k++)
            {
                for(n=0;n<N_DIMS;n++)
                {
                    mt[pair.inner_nlist[k][0]].r_eqnm[n] = mt[pair.inner_nlist[k][0]].r_eqnm[n] - pair.F_ij_vec[k]*pair.uw_ij[k][n];
                    mt[pair.inner_nlist[k][1]].r_eqnm[n] = mt[pair.inner_nlist[k][1]].r_eqnm[n] + pair.F_ij_vec[k]*pair.uw_ij[k][n];
                }
                
                u_iw_dot = mt[pair.inner_nlist[k][0]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][0]].u_vec[1]*pair.uw_ij[k][1] + mt[pair.inner_nlist[k][0]].u_vec[2]*pair.uw_ij[k][2];
                u_jw_dot = mt[pair.inner_nlist[k][1]].u_vec[0]*pair.uw_ij[k][0] + mt[pair.inner_nlist[k][1]].u_vec[1]*pair.uw_ij[k][1] + mt[pair.inner_nlist[k][1]].u_vec[2]*pair.uw_ij[k][2];
                
                for(n=0;n<N_DIMS;n++)
                {
                    mt[pair.inner_nlist[k][0]].u_eqnm[n] = mt[pair.inner_nlist[k][0]].u_eqnm[n] + pair.F_ij_vec[k]*pair.lambda_array[k][0]*(u_iw_dot*mt[pair.inner_nlist[k][0]].u_vec[n] - pair.uw_ij[k][n]);
                    mt[pair.inner_nlist[k][1]].u_eqnm[n] = mt[pair.inner_nlist[k][1]].u_eqnm[n] - pair.F_ij_vec[k]*pair.lambda_array[k][1]*(u_jw_dot*mt[pair.inner_nlist[k][1]].u_vec[n] - pair.uw_ij[k][n]);
                }
            }
        }
        
        // Calculate the wall boundary conditions (wall-force)
        if(BC_TYPE == 2 || BC_TYPE == 3)
        {
            if(N_DIMS == 2)
            {
                // Create a wall pair list for neighboring filaments (2D)
                if (i % wall_update == 0)
                {
                    wall.wall_init();
                    
                    for (j=N_CRYST; j < N_FILS; j++)
                    {
                        
                        ri_s2d[0]   = 0.0;
                        ri_s2d[1]   = mt[j].r_com[1];
                        uwall_2d[0] = 0.0;
                        uwall_2d[1] = 1.0;
                        
                        // L_wall_2D   = mt[j].fil_length*sqrt(1 - mt[j].u_vec[0]*mt[j].u_vec[0]);
                        L_wall_2D = L_y;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s2d, mt[j].u_vec,uwall_2d,mt[j].fil_length,L_wall_2D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= 1.0e-25 || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){cout << "warning wall x1 - " << j << endl;}
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x1_add(j);}
                        
                        ri_s2d[0]  = L_x;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s2d, mt[j].u_vec,uwall_2d,mt[j].fil_length,L_wall_2D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= 1.0e-25 || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){cout << "warning wall x2 - " << j << endl;}
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x2_add(j);}
                    }
                }
                // Calculate wall force for each listed wall-filament interaction (do for 2 walls in xy plane)
                {
                    for (j=0; j < wall.x1_n; j++)
                    {
                        ri_s2d[0]  = 0.0;
                        ri_s2d[1]  = mt[wall.wl_x1[j]].r_com[1];
                        uwall_2d[0] = 0.0;
                        uwall_2d[1] = 1.0;
                        F_i_wall    = 0.0;
                        
                        // L_wall_2D   = mt[wall.wl_x1[j]].fil_length*sqrt(1 - mt[wall.wl_x1[j]].u_vec[0]*mt[wall.wl_x1[j]].u_vec[0]);
                        L_wall_2D = L_y;
                        
                        wall_calc.wall_calc_class(mt[wall.wl_x1[j]].r_com, ri_s2d, mt[wall.wl_x1[j]].u_vec,uwall_2d,mt[wall.wl_x1[j]].fil_length,L_wall_2D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= sig_min_wall)
                        {
                            F_i_wall = wall_calc.wall_wca_calc();
                            
                            for(n=0;n<N_DIMS;n++)
                            {
                                mt[wall.wl_x1[j]].r_eqnm[n] = mt[wall.wl_x1[j]].r_eqnm[n] - F_i_wall*uwall1[n];
                            }
                            
                            u_iw_dot = mt[wall.wl_x1[j]].u_vec[0]*uwall1[0] + mt[wall.wl_x1[j]].u_vec[1]*uwall1[1];
                            
                            for(n=0;n<N_DIMS;n++)
                            {
                                mt[wall.wl_x1[j]].u_eqnm[n] = mt[wall.wl_x1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x1[j]].u_vec[n] - uwall1[n]);
                            }
                        }
                    }
                    for (j=0; j < wall.x2_n; j++)
                    {
                        ri_s2d[0]  = L_x;
                        ri_s2d[1]  = mt[wall.wl_x2[j]].r_com[1];
                        uwall_2d[0] = 0.0;
                        uwall_2d[1] = 1.0;
                        F_i_wall    = 0.0;
                        
                        // L_wall_2D   = mt[wall.wl_x2[j]].fil_length*sqrt(1 - mt[wall.wl_x2[j]].u_vec[0]*mt[wall.wl_x2[j]].u_vec[0]);
                        L_wall_2D = L_y;
                        
                        wall_calc.wall_calc_class(mt[wall.wl_x2[j]].r_com, ri_s2d, mt[wall.wl_x2[j]].u_vec,uwall_2d,mt[wall.wl_x2[j]].fil_length,L_wall_2D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= sig_min_wall)
                        {
                            F_i_wall = wall_calc.wall_wca_calc();
                            
                            for(n=0;n<N_DIMS;n++)
                            {
                                mt[wall.wl_x2[j]].r_eqnm[n] = mt[wall.wl_x2[j]].r_eqnm[n] - F_i_wall*uwall2[n];
                            }
                            
                            u_iw_dot = mt[wall.wl_x2[j]].u_vec[0]*uwall2[0] + mt[wall.wl_x2[j]].u_vec[1]*uwall2[1];
                            
                            for(n=0;n<N_DIMS;n++)
                            {
                                mt[wall.wl_x2[j]].u_eqnm[n] = mt[wall.wl_x2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x2[j]].u_vec[n] - uwall2[n]);
                            }
                        }
                    }
                }
            }
            
            if(N_DIMS == 3)
            {
                // Create a wall pair list for neighboring filaments (3D)
                if (i % wall_update == 0)
                {
                    wall.wall_init();
                    
                    for (j = N_CRYST; j < N_FILS; j++)
                    {
                        
                        ri_s3d[0]  = 0.0;
                        ri_s3d[1]  = mt[j].r_com[1];
                        ri_s3d[2]  = mt[j].r_com[2];
                        
                        // must normalise:
                        uwall_3d[0] = 0.0;
                        uwall_3d[1] = mt[j].u_vec[1];
                        uwall_3d[2] = mt[j].u_vec[2];
                        
                        double umod = uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2];
                        umod = sqrt(umod);
                        
                        uwall_3d[1] = uwall_3d[1]/umod;
                        uwall_3d[2] = uwall_3d[2]/umod;
                        
                        // L_wall_3D   = mt[j].fil_length*sqrt(1.0 - mt[j].u_vec[0]*mt[j].u_vec[0]);
                        L_wall_3D = 20.0*L_y;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec, uwall_3d, mt[j].fil_length, L_wall_3D, N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x1_add(j);}
                        
                        ri_s3d[0]  = L_x;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_x2_add(j);}
                        
                        ri_s3d[0]  = mt[j].r_com[0];
                        ri_s3d[1]  = mt[j].r_com[1];
                        ri_s3d[2]  = 0.0;
                        
                        // must normalise:
                        uwall_3d[0] = mt[j].u_vec[0];
                        uwall_3d[1] = mt[j].u_vec[1];
                        uwall_3d[2] = 0.0;
                        
                        umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1];
                        umod = sqrt(umod);
                        
                        uwall_3d[0] = uwall_3d[0]/umod;
                        uwall_3d[1] = uwall_3d[1]/umod;
                        
                        // L_wall_3D   = mt[j].fil_length*sqrt(1.0 - mt[j].u_vec[2]*mt[j].u_vec[2]);
                        L_wall_3D = 20.0*L_y;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_z1_add(j);}
                        
                        ri_s3d[2]  = L_z;
                        
                        wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                        
                        if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1){wall.wall_z2_add(j);}
                    }
                }
                
                // Calculate wall force for each listed wall-filament interaction (do for 2 walls in xy plane and 2 in yz plane)
                for (j=0; j < wall.x1_n; j++)
                {
                    // wall_here
                    
                    ri_s3d[0]  = 0.0;
                    ri_s3d[1]  = mt[wall.wl_x1[j]].r_com[1];
                    ri_s3d[2]  = mt[wall.wl_x1[j]].r_com[2];
                    
                    // must normalise:
                    uwall_3d[0] = 0.0;
                    uwall_3d[1] = mt[wall.wl_x1[j]].u_vec[1];
                    uwall_3d[2] = mt[wall.wl_x1[j]].u_vec[2];
                    
                    double umod = uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2];
                    umod = sqrt(umod);
                    
                    uwall_3d[1] = uwall_3d[1]/umod;
                    uwall_3d[2] = uwall_3d[2]/umod;
                    
                    L_wall_3D = 20.0*L_y;
                    
                    F_i_wall    = 0.0;
                    
                    wall_calc.wall_calc_class(mt[wall.wl_x1[j]].r_com, ri_s3d, mt[wall.wl_x1[j]].u_vec, uwall_3d, mt[wall.wl_x1[j]].fil_length, L_wall_3D, N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= sig_min_wall)
                    {
                        
                        F_i_wall = wall_calc.wall_wca_calc();
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_x1[j]].r_eqnm[n] = mt[wall.wl_x1[j]].r_eqnm[n] - F_i_wall*uwall1[n];
                        }
                        
                        u_iw_dot = mt[wall.wl_x1[j]].u_vec[0]*uwall1[0] + mt[wall.wl_x1[j]].u_vec[1]*uwall1[1] + mt[wall.wl_x1[j]].u_vec[2]*uwall1[2];
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_x1[j]].u_eqnm[n] = mt[wall.wl_x1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x1[j]].u_vec[n] - uwall1[n]);
                        }
                    }
                }
                for (j=0; j < wall.x2_n; j++)
                {
                    
                    ri_s3d[0]  = L_x;
                    ri_s3d[1]  = mt[wall.wl_x2[j]].r_com[1];
                    ri_s3d[2]  = mt[wall.wl_x2[j]].r_com[2];
                    
                    // must normalise:
                    uwall_3d[0] = 0.0;
                    uwall_3d[1] = mt[wall.wl_x2[j]].u_vec[1];
                    uwall_3d[2] = mt[wall.wl_x2[j]].u_vec[2];
                    
                    double umod = uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2];
                    umod = sqrt(umod);
                    
                    uwall_3d[1] = uwall_3d[1]/umod;
                    uwall_3d[2] = uwall_3d[2]/umod;
                    
                    // L_wall_3D = mt[wall.wl_x2[j]].fil_length*sqrt(1.0 - mt[wall.wl_x2[j]].u_vec[0]*mt[wall.wl_x2[j]].u_vec[0]);
                    L_wall_3D = 20.0*L_y;
                    
                    F_i_wall    = 0.0;
                    
                    wall_calc.wall_calc_class(mt[wall.wl_x2[j]].r_com, ri_s3d, mt[wall.wl_x2[j]].u_vec,uwall_3d,mt[wall.wl_x2[j]].fil_length,L_wall_3D,N_DIMS,j);
                    
                    
                    if(wall_calc.w_ij_mag_out <= sig_min_wall)
                    {
                        
                        F_i_wall = wall_calc.wall_wca_calc();
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_x2[j]].r_eqnm[n] = mt[wall.wl_x2[j]].r_eqnm[n] - F_i_wall*uwall2[n];
                        }
                        
                        u_iw_dot = mt[wall.wl_x2[j]].u_vec[0]*uwall2[0] + mt[wall.wl_x2[j]].u_vec[1]*uwall2[1] + mt[wall.wl_x2[j]].u_vec[2]*uwall2[2];
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_x2[j]].u_eqnm[n] = mt[wall.wl_x2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_x2[j]].u_vec[n] - uwall2[n]);
                        }
                    }
                }
                for (j=0; j < wall.z1_n; j++)
                {
                    
                    ri_s3d[0]  = mt[wall.wl_z1[j]].r_com[0];
                    ri_s3d[1]  = mt[wall.wl_z1[j]].r_com[1];
                    ri_s3d[2]  = 0.0;
                    
                    // must normalise:
                    uwall_3d[0] = mt[wall.wl_z1[j]].u_vec[0];
                    uwall_3d[1] = mt[wall.wl_z1[j]].u_vec[1];
                    uwall_3d[2] = 0.0;
                    
                    double umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1];
                    umod = sqrt(umod);
                    
                    uwall_3d[0] = uwall_3d[0]/umod;
                    uwall_3d[1] = uwall_3d[1]/umod;
                    
                    // L_wall_3D   = mt[wall.wl_z1[j]].fil_length*sqrt(1.0 - mt[wall.wl_z1[j]].u_vec[2]*mt[wall.wl_z1[j]].u_vec[2]);
                    L_wall_3D = 20.0*L_y;
                    
                    F_i_wall    = 0.0;
                    
                    wall_calc.wall_calc_class(mt[wall.wl_z1[j]].r_com, ri_s3d, mt[wall.wl_z1[j]].u_vec,uwall_3d,mt[wall.wl_z1[j]].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= sig_min_wall)
                    {

                        F_i_wall = wall_calc.wall_wca_calc();
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_z1[j]].r_eqnm[n] = mt[wall.wl_z1[j]].r_eqnm[n] - F_i_wall*uwall3[n];
                        }
                        
                        u_iw_dot = mt[wall.wl_z1[j]].u_vec[0]*uwall3[0] + mt[wall.wl_z1[j]].u_vec[1]*uwall3[1] + mt[wall.wl_z1[j]].u_vec[2]*uwall3[2];
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_z1[j]].u_eqnm[n] = mt[wall.wl_z1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z1[j]].u_vec[n] - uwall3[n]);
                        }
                    }
                }
                for (j=0; j < wall.z2_n; j++)
                {
                    ri_s3d[0]  = mt[wall.wl_z2[j]].r_com[0];
                    ri_s3d[1]  = mt[wall.wl_z2[j]].r_com[1];
                    ri_s3d[2]  = L_z;
                    
                    // must normalise:
                    uwall_3d[0] = mt[wall.wl_z2[j]].u_vec[0];
                    uwall_3d[1] = mt[wall.wl_z2[j]].u_vec[1];
                    uwall_3d[2] = 0.0;
                    
                    double umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1];
                    umod = sqrt(umod);
                    
                    uwall_3d[0] = uwall_3d[0]/umod;
                    uwall_3d[1] = uwall_3d[1]/umod;
                    
                    // L_wall_3D   = mt[wall.wl_z2[j]].fil_length*sqrt(1.0 - mt[wall.wl_z2[j]].u_vec[2]*mt[wall.wl_z2[j]].u_vec[2]);
                    L_wall_3D = 20.0*L_y;
                    
                    F_i_wall = 0.0;
                    
                    wall_calc.wall_calc_class(mt[wall.wl_z2[j]].r_com, ri_s3d, mt[wall.wl_z2[j]].u_vec,uwall_3d,mt[wall.wl_z2[j]].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= sig_min_wall)
                    {

                        F_i_wall = wall_calc.wall_wca_calc();
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_z2[j]].r_eqnm[n] = mt[wall.wl_z2[j]].r_eqnm[n] - F_i_wall*uwall4[n];
                        }
                        
                        u_iw_dot = mt[wall.wl_z2[j]].u_vec[0]*uwall4[0] + mt[wall.wl_z2[j]].u_vec[1]*uwall4[1] + mt[wall.wl_z2[j]].u_vec[2]*uwall4[2];
                        
                        for(n=0;n<N_DIMS;n++)
                        {
                            mt[wall.wl_z2[j]].u_eqnm[n] = mt[wall.wl_z2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z2[j]].u_vec[n] - uwall4[n]);
                        }
                    }
                }
            }
        }
        if(BC_TYPE == 4 || BC_TYPE == 5)
        {
            // Create a wall pair list for neighboring filaments (2D)
            if (i % wall_update == 0)
            {
                wall.wall_init();
                
                for (j=0; j < N_FILS; j++)
                {
                    
                    ri_s3d[0]  = mt[j].r_com[0];
                    ri_s3d[1]  = mt[j].r_com[1];
                    ri_s3d[2]  = 0.0;
                    
                    // must normalise:
                    uwall_3d[0] = mt[j].u_vec[0];
                    uwall_3d[1] = mt[j].u_vec[1];
                    uwall_3d[2] = 0.0;
                    
                    double umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1] + uwall_3d[2]*uwall_3d[2];
                    
                    umod = sqrt(umod);

                    uwall_3d[0] = uwall_3d[0]/umod;
                    uwall_3d[1] = uwall_3d[1]/umod;
                    uwall_3d[2] = uwall_3d[2]/umod;
                    
                    L_wall_3D = L_y;
                    

                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec, uwall_3d, mt[j].fil_length, L_wall_3D, N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1)
                    {
                        wall.wall_z1_add(j);
                    }
                    
                    ri_s3d[2]  = L_z;
                    
                    wall_calc.wall_calc_class(mt[j].r_com, ri_s3d, mt[j].u_vec,uwall_3d,mt[j].fil_length,L_wall_3D,N_DIMS,j);
                    
                    if(wall_calc.w_ij_mag_out <= wall_cut*sig_min_wall || wall_calc.w_ij_mag_out == 0.0 || isnan(wall_calc.w_ij_mag_out) == 1)
                    {
                        wall.wall_z2_add(j);
                    }
                }
            }
            for (j=0; j < wall.z1_n; j++)
            {
                ri_s3d[0]  = mt[wall.wl_z1[j]].r_com[0];
                ri_s3d[1]  = mt[wall.wl_z1[j]].r_com[1];
                ri_s3d[2]  = 0.0;
                
                // must normalise:
                uwall_3d[0] = mt[wall.wl_z1[j]].u_vec[0];
                uwall_3d[1] = mt[wall.wl_z1[j]].u_vec[1];
                uwall_3d[2] = 0.0;
                
                double umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1];
                umod = sqrt(umod);
                
                uwall_3d[0] = uwall_3d[0]/umod;
                uwall_3d[1] = uwall_3d[1]/umod;
                
                L_wall_3D = 20.0*L_y;
                
                F_i_wall = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_z1[j]].r_com, ri_s3d, mt[wall.wl_z1[j]].u_vec,uwall_3d,mt[wall.wl_z1[j]].fil_length,L_wall_3D,N_DIMS,j);
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {

                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0; n < N_DIMS; n++){mt[wall.wl_z1[j]].r_eqnm[n] = mt[wall.wl_z1[j]].r_eqnm[n] - F_i_wall*uwall3[n];}
                    
                    u_iw_dot = mt[wall.wl_z1[j]].u_vec[0]*uwall3[0] + mt[wall.wl_z1[j]].u_vec[1]*uwall3[1] + mt[wall.wl_z1[j]].u_vec[2]*uwall3[2];
                    
                    for(n=0; n < N_DIMS; n++)
                    {
                        mt[wall.wl_z1[j]].u_eqnm[n] = mt[wall.wl_z1[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z1[j]].u_vec[n] - uwall3[n]);
                    }
                }
            }
            for (j=0; j < wall.z2_n; j++)
            {
                
                ri_s3d[0]  = mt[wall.wl_z2[j]].r_com[0];
                ri_s3d[1]  = mt[wall.wl_z2[j]].r_com[1];
                ri_s3d[2]  = L_z;
                
                // must normalise:
                uwall_3d[0] = mt[wall.wl_z2[j]].u_vec[0];
                uwall_3d[1] = mt[wall.wl_z2[j]].u_vec[1];
                uwall_3d[2] = 0.0;
                
                double umod = uwall_3d[0]*uwall_3d[0] + uwall_3d[1]*uwall_3d[1];
                umod = sqrt(umod);
                
                uwall_3d[0] = uwall_3d[0]/umod;
                uwall_3d[1] = uwall_3d[1]/umod;
                
                // L_wall_3D   = mt[wall.wl_z2[j]].fil_length*sqrt(1.0 - mt[wall.wl_z2[j]].u_vec[2]*mt[wall.wl_z2[j]].u_vec[2]);
                L_wall_3D = 20.0*L_y;
                
                F_i_wall = 0.0;
                
                wall_calc.wall_calc_class(mt[wall.wl_z2[j]].r_com, ri_s3d, mt[wall.wl_z2[j]].u_vec,uwall_3d,mt[wall.wl_z2[j]].fil_length,L_wall_3D,N_DIMS,j);
                
                if(wall_calc.w_ij_mag_out <= sig_min_wall)
                {
                    
                    F_i_wall = wall_calc.wall_wca_calc();
                    
                    for(n=0;n<N_DIMS;n++){mt[wall.wl_z2[j]].r_eqnm[n] = mt[wall.wl_z2[j]].r_eqnm[n] - F_i_wall*uwall4[n];}
                    
                    u_iw_dot = mt[wall.wl_z2[j]].u_vec[0]*uwall4[0] + mt[wall.wl_z2[j]].u_vec[1]*uwall4[1] + mt[wall.wl_z2[j]].u_vec[2]*uwall4[2];
                    
                    for(n=0;n<N_DIMS;n++)
                    {
                        mt[wall.wl_z2[j]].u_eqnm[n] = mt[wall.wl_z2[j]].u_eqnm[n] + F_i_wall*wall_calc.lambda[0]*(u_iw_dot*mt[wall.wl_z2[j]].u_vec[n] - uwall4[n]);
                    }
                }
            }
        }
        
        // Add forces due to cross-linker/motor interactions to Eqns of Motion (motor-force)
        if(CL == 1 || CL == 2)
        {
            for(unsigned int k = 0; k < cl_nst2; k++)
            {

                rcl_image = binding.image_calc_BC_cl(BC_TYPE, N_DIMS, mt[cl[cl_st2_list[k]].fil_i].r_com, mt[cl[cl_st2_list[k]].fil_j].r_com, box_vec);
                
                d_ij_mag = 0.0;
                
                for(n=0; n<N_DIMS; n++)
                {
                    d_ij[n] = (rcl_image[n] + cl[cl_st2_list[k]].epsilon_i*mt[cl[cl_st2_list[k]].fil_i].u_vec[n]) - (mt[cl[cl_st2_list[k]].fil_j].r_com[n] + cl[cl_st2_list[k]].epsilon_j*mt[cl[cl_st2_list[k]].fil_j].u_vec[n]);
                    d_ij_mag = d_ij_mag + d_ij[n]*d_ij[n];
                }
                
                d_ij_mag = sqrt(d_ij_mag);
                
                for(n=0; n<N_DIMS; n++)
                {
                    d_0_vec[n] = (d_0/d_ij_mag)*d_ij[n];
                }
                
                for(n=0; n<N_DIMS; n++)
                {
                    mt[cl[cl_st2_list[k]].fil_i].r_eqnm[n] = mt[cl[cl_st2_list[k]].fil_i].r_eqnm[n] - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n]);
                    mt[cl[cl_st2_list[k]].fil_j].r_eqnm[n] = mt[cl[cl_st2_list[k]].fil_j].r_eqnm[n] + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n]);
                }
        
                f_pi = 0.0;
                f_pj = 0.0;
                
                if(cl[cl_st2_list[k]].type == 1 && d_ij_mag > d_0) // only calculate under extension
                {
                    
                    double uidijk_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*d_ij[0] + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*d_ij[1] + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*d_ij[2];
                    double ujdjik_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(-d_ij[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(-d_ij[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(-d_ij[2]);
                    
                    if(uidijk_dot > 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pi = f_pi - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_i].u_vec[n];
                        }
                    }
                    if(ujdjik_dot > 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pj = f_pj + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_j].u_vec[n];
                        }
                    }
                }
                if(cl[cl_st2_list[k]].type == 2 && d_ij_mag > d_0) // only calculate under extension
                {
                    
                    double uidijd_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*d_ij[0] + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*d_ij[1] + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*d_ij[2];
                    double ujdjid_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(-d_ij[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(-d_ij[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(-d_ij[2]);
                    
                    if(uidijd_dot < 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pi = f_pi - cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_i].u_vec[n];
                        }
                    }
                    if(ujdjid_dot < 0.0)
                    {
                        for(n=0; n<N_DIMS; n++)
                        {
                            f_pj = f_pj + cl[cl_st2_list[k]].k_spr*(d_ij[n]-d_0_vec[n])*mt[cl[cl_st2_list[k]].fil_j].u_vec[n];
                        }
                    }
                }
  
                u_id_dot = mt[cl[cl_st2_list[k]].fil_i].u_vec[0]*(d_ij[0]-d_0_vec[0]) + mt[cl[cl_st2_list[k]].fil_i].u_vec[1]*(d_ij[1]-d_0_vec[1]) + mt[cl[cl_st2_list[k]].fil_i].u_vec[2]*(d_ij[2]-d_0_vec[2]);
                u_jd_dot = mt[cl[cl_st2_list[k]].fil_j].u_vec[0]*(d_ij[0]-d_0_vec[0]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[1]*(d_ij[1]-d_0_vec[1]) + mt[cl[cl_st2_list[k]].fil_j].u_vec[2]*(d_ij[2]-d_0_vec[2]);
                
                
                for(n=0;n<N_DIMS;n++)
                {
                    mt[cl[cl_st2_list[k]].fil_i].u_eqnm[n] = mt[cl[cl_st2_list[k]].fil_i].u_eqnm[n] + cl[cl_st2_list[k]].k_spr*cl[cl_st2_list[k]].epsilon_i*(u_id_dot*mt[cl[cl_st2_list[k]].fil_i].u_vec[n] - (d_ij[n]-d_0_vec[n]));
                    mt[cl[cl_st2_list[k]].fil_j].u_eqnm[n] = mt[cl[cl_st2_list[k]].fil_j].u_eqnm[n] - cl[cl_st2_list[k]].k_spr*cl[cl_st2_list[k]].epsilon_j*(u_jd_dot*mt[cl[cl_st2_list[k]].fil_j].u_vec[n] - (d_ij[n]-d_0_vec[n]));
                }

                cl[cl_st2_list[k]].f_para_i = f_pi;
                cl[cl_st2_list[k]].f_para_j = f_pj;
            }
        }

        // numerically update the Eqns of Motion (2D or 3D). Include thermal noise contribution (Brownian-force)
        if(N_DIMS == 2)
        {
            for (j=N_CRYST; j < N_FILS; j++)
            {
                // Function call from utilities (Choose between the two)
                
                filament_drag_2d_util(mt[j].u_vec,mt[j].r_eqnm,mt[j].u_eqnm,gamma_array_2d,mt[j].gamma_para,mt[j].gamma_perp,mt[j].gamma_rot);
                
                mt[j].r_com[0] = mt[j].r_com[0] + dt*mt[j].r_eqnm[0];
                mt[j].r_com[1] = mt[j].r_com[1] + dt*mt[j].r_eqnm[1];
                
                mt[j].u_vec[0] = mt[j].u_vec[0] + dt*mt[j].u_eqnm[0];
                mt[j].u_vec[1] = mt[j].u_vec[1] + dt*mt[j].u_eqnm[1];
                
                u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
                u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}

            }
        }
        else if(N_DIMS == 3)
        {
            for (j = N_CRYST; j < N_FILS; j++)
            {
                
                filament_drag_3d_util(mt[j].u_vec,mt[j].r_eqnm,mt[j].u_eqnm,gamma_array_3d,mt[j].gamma_para,mt[j].gamma_perp,mt[j].gamma_rot);
                
                mt[j].r_com[0] = mt[j].r_com[0] + dt*mt[j].r_eqnm[0];
                mt[j].r_com[1] = mt[j].r_com[1] + dt*mt[j].r_eqnm[1];
                mt[j].r_com[2] = mt[j].r_com[2] + dt*mt[j].r_eqnm[2];
                
                mt[j].u_vec[0] = mt[j].u_vec[0] + dt*mt[j].u_eqnm[0];
                mt[j].u_vec[1] = mt[j].u_vec[1] + dt*mt[j].u_eqnm[1];
                mt[j].u_vec[2] = mt[j].u_vec[2] + dt*mt[j].u_eqnm[2];
                
                u_mod = 0.0; for (unsigned int n = 0; n < N_DIMS; n++){u_mod = u_mod + mt[j].u_vec[n]*mt[j].u_vec[n];}
                u_mod = sqrt(u_mod); for (unsigned int n = 0; n < N_DIMS; n++){mt[j].u_vec[n] = mt[j].u_vec[n]/u_mod;}
            }
        
        }
        
        // correct for PBCs if required (minimal image calculations)
        if(BC_TYPE == 1 )
        {
            for (j=0; j < N_FILS; j++)
            {
                pbc.pbc_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
            }
        }
        if(BC_TYPE == 2)
        {
            for (j=0; j < N_FILS; j++)
            {
                pbc.pbc_y_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
            }
        }
        if(BC_TYPE == 4)
        {
            for (j=0; j < N_FILS; j++)
            {
                // pbc.pbc_yz_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
                pbc.pbc_calc(mt[j].r_com,N_DIMS,mt[j].pbc_index);
            }
        }
        
        if (i % t_write == 0)
        {
            // Main trajectory file for filaments
            trajectory << i << "\t" << N_FILS << "\t";
            for (j=0; j < N_FILS; j++)
            {
                trajectory << mt[j].fil_length << "\t"  << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2]<< "\t";
            }

            trajectory << "\n";
            
            // output CL trajectories for visualisation. Only if required.
            if(CL_PLOT == 1)
            {
                cl_trajectory << i << "\t" << cl_nst1 << "\t" << cl_nst2 << "\t";
                
                if(N_DIMS == 2)
                {
                    for (j=0; j < N_CLS; j++)
                    {
                        if(cl[j].state == 1)
                        {
                            cl_trajectory << cl[j].state << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t";
                        }
                        if(cl[j].state == 2)
                        {
                            cl_trajectory << cl[j].state << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t";
                        }
                    }
                }
                else if(N_DIMS == 3)
                {
                    for (j=0; j < N_CLS; j++)
                    {
                        if(cl[j].state == 1)
                        {
                            cl_trajectory << cl[j].state << "\t" << cl[j].type << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t";
                        }
                        if(cl[j].state == 2)
                        {
                            cl_trajectory << cl[j].state << "\t" << cl[j].type << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[2] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[2] << "\t";
                        }
                    }
                }
                cl_trajectory << endl;
            }
            
            // output just the CL populations
            nt_cl1_1 = 0; nt_cl1_2 = 0;
            for (j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 2)
                {
                    if(cl[j].type == 1){nt_cl1_1 = nt_cl1_1 + 1;}
                    if(cl[j].type == 2){nt_cl1_2 = nt_cl1_2 + 1;}
                }
            }
            cl_nt << i << "\t" << cl_nst1 << "\t" << cl_nst2 << "\t" << nt_cl1_1 << "\t" << nt_cl1_2 << endl;
            
        }
        
        // split time scale output for detailed MSD calculations
        if (i >= t_begin_write && i < t_begin_write + tscale_1_upper)
        {
            if(i % tscale_1_inc == 0)
            {
                traj_tscale_1 << i << "\t" << N_MSD_Calc << "\t";

                for (j=0; j < N_MSD_Calc; j++)
                {
                    traj_tscale_1 << mt[j].fil_length << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1]  << "\t" << mt[j].u_vec[2] << "\t";
                }

                traj_tscale_1 << "\n";
            }
        }

        if (i >= t_begin_write && i < t_begin_write + tscale_2_upper)
        {
            if(i % tscale_2_inc == 0)
            {
                traj_tscale_2 << i << "\t" << N_MSD_Calc << "\t";

                for (j=0; j < N_MSD_Calc; j++)
                {
                    traj_tscale_2 << mt[j].fil_length << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1]  << "\t" << mt[j].u_vec[2] << "\t";
                }

                traj_tscale_2 << "\n";
            }
        }

        if (i >= t_begin_write)
        {
            if(i % tscale_3_inc == 0)
            {
                traj_tscale_3 << i << "\t" << N_MSD_Calc << "\t";

                for (j=0; j < N_MSD_Calc; j++)
                {
                    traj_tscale_3 << mt[j].fil_length << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1]  << "\t" << mt[j].u_vec[2] << "\t";
                }

                traj_tscale_3 << "\n";
            }
        }
        
        // Currently the following is unused. It contains advanced filament outputting where
        // filament details are written to file depending on ID and whether they belong
        // to cross-linked networks. Will likely not use, but do not want to delete for now,
        // since it is potentially very useful. Original prompt: if(i >= start_write)
        if(start_write == 9999)
        {
            if (i % t_write_2 == 0)
            {
                // write all filaments with IDs
                id_trajectory << double(i)*dt << "\t" << N_FILS - N_CRYST << "\t" << filament_id_count << "\t";
                id_orientation << double(i)*dt << "\t" << N_FILS - N_CRYST << "\t" << filament_id_count << "\t";
                if(N_DIMS == 3)
                {
                    for (j = N_CRYST; j < N_FILS; j++)
                    {
                        id_trajectory << mt[j].fil_id << "\t" << mt[j].r_com[0] - 0.5*mt[j].fil_length*mt[j].u_vec[0] << "\t" << mt[j].r_com[1] - 0.5*mt[j].fil_length*mt[j].u_vec[1] << "\t" << mt[j].r_com[2] - 0.5*mt[j].fil_length*mt[j].u_vec[2] << "\t";
                        id_orientation << mt[j].fil_id << "\t" << mt[j].fil_length << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2] << "\t";
                    }
                }
                id_trajectory  << "\n";
                id_orientation << "\n";

                // write only filaments connected to pinning field (No IDs)
                pinning_stress << double(i)*dt << "\t";
                for (j=0; j < N_CLS; j++)
                {
                    if(cl[j].state == 2)
                    {
                        if(cl[j].fil_i < N_CRYST && cl[j].fil_j >= N_CRYST )
                        {
                            pinning_stress << mt[cl[j].fil_i].r_com[0] << "\t" << mt[cl[j].fil_i].r_com[1] << "\t" << mt[cl[j].fil_i].r_com[2] << "\t" << mt[cl[j].fil_j].r_com[0] << "\t" << mt[cl[j].fil_j].r_com[1] << "\t" << mt[cl[j].fil_j].r_com[2] << "\t";
                            pinning_stress << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[2] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[2] << "\t";
                        }
                        if(cl[j].fil_i >= N_CRYST && cl[j].fil_j < N_CRYST )
                        {
                            pinning_stress << mt[cl[j].fil_i].r_com[0] << "\t" << mt[cl[j].fil_i].r_com[1] << "\t" << mt[cl[j].fil_i].r_com[2] << "\t" << mt[cl[j].fil_j].r_com[0] << "\t" << mt[cl[j].fil_j].r_com[1] << "\t" << mt[cl[j].fil_j].r_com[2] << "\t";
                            pinning_stress << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[2] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[2] << "\t";
                        }
                    }
                }
                pinning_stress <<  "\n";
                
                // Write out the filament IDs for each pair-bound CL
                cl_id << i << "\t" << N_FILS - N_CRYST << "\t";
                for(unsigned int k = 0; k < cl_nst2; k++)
                {
                    if(cl[cl_st2_list[k]].fil_i >= N_CRYST && cl[cl_st2_list[k]].fil_j >= N_CRYST)
                    {
                        cl_id << cl[cl_st2_list[k]].fil_i << "\t" << cl[cl_st2_list[k]].fil_j << "\t";
                    }
                }
                cl_id << endl;
                
                // Write the filament traj file in IDed format only for filaments invloved in a pair bond - newtworked filaments only (exclude pins)
                network_fils.resize(2*cl_nst2); count = 0;
                for(unsigned int k = 0; k < cl_nst2; k++)
                {
                    if(cl[cl_st2_list[k]].fil_i >= N_CRYST)
                    {
                        network_fils[count] = cl[cl_st2_list[k]].fil_i;
                        count = count + 1;
                    }
                    
                    if(cl[cl_st2_list[k]].fil_j >= N_CRYST)
                    {
                        network_fils[count] = cl[cl_st2_list[k]].fil_j;
                        count = count + 1;
                    }
                }
                
                std::sort (network_fils.begin(), network_fils.end());
                
                int Network_count = 0;
                network_fils_filtered.resize(0);
                if(cl_nst2 > 0)
                {
                    network_fils_filtered.push_back(network_fils[0]);
                    int last_fil = network_fils_filtered[0];
                    Network_count = 1;
                    for(unsigned int k = 1; k < 2*cl_nst2; k++)
                    {
                        if(network_fils[k] != last_fil)
                        {
                            network_fils_filtered.push_back(network_fils[k]);
                            last_fil = network_fils[k];
                            Network_count = Network_count + 1;
                        }
                    }
                }

                // write all filaments with IDs
                network_trajectory << double(i)*dt << "\t" << Network_count << "\t" << filament_id_count << "\t";
                network_orientation << double(i)*dt << "\t" << Network_count << "\t" << filament_id_count << "\t";
                if(N_DIMS == 3)
                {
                    for (j = 0; j < Network_count; j++)
                    {
                        network_trajectory << mt[network_fils_filtered[j]].fil_id << "\t" << mt[network_fils_filtered[j]].r_com[0] - 0.5*mt[network_fils_filtered[j]].fil_length*mt[network_fils_filtered[j]].u_vec[0] << "\t" << mt[network_fils_filtered[j]].r_com[1] - 0.5*mt[network_fils_filtered[j]].fil_length*mt[network_fils_filtered[j]].u_vec[1] << "\t" << mt[network_fils_filtered[j]].r_com[2] - 0.5*mt[network_fils_filtered[j]].fil_length*mt[network_fils_filtered[j]].u_vec[2] << "\t";
                        network_orientation << mt[network_fils_filtered[j]].fil_id << "\t" << mt[network_fils_filtered[j]].fil_length << "\t" << mt[network_fils_filtered[j]].u_vec[0] << "\t" << mt[network_fils_filtered[j]].u_vec[1] << "\t" << mt[network_fils_filtered[j]].u_vec[2] << "\t";
                    }
                }
                network_trajectory  << "\n";
                network_orientation << "\n";
                
                network_fils.resize(0);
                network_fils_filtered.resize(0);                
                
            // write the CL states out for stress calculation rather then for visualisation
            if(CL_PLOT == 2)
            {
                cl_trajectory << i << "\t" << cl_nst2 << "\t";
                
                if(N_DIMS == 3)
                {
                    for (j=0; j < N_CLS; j++)
                    {
                        if(cl[j].state == 2)
                        {
                            cl_trajectory << cl[j].fil_i  << "\t" << cl[j].fil_j << "\t" << cl[j].type << "\t" << mt[cl[j].fil_i].r_com[0] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[0] << "\t" << mt[cl[j].fil_i].r_com[1] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[1] << "\t" << mt[cl[j].fil_i].r_com[2] + cl[j].epsilon_i*mt[cl[j].fil_i].u_vec[2] << "\t" << mt[cl[j].fil_j].r_com[0] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[0] << "\t" << mt[cl[j].fil_j].r_com[1] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[1] << "\t" << mt[cl[j].fil_j].r_com[2] + cl[j].epsilon_j*mt[cl[j].fil_j].u_vec[2] << "\t";
                        }
                    }
                }
                cl_trajectory << endl;
            }
	}
        }
        
        if (i % t_print == 0)
        {
            nt_cl1_1 = 0; nt_cl1_2 = 0;
            
            for (j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 2)
                {
                    if(cl[j].type == 1){nt_cl1_1 = nt_cl1_1 + 1;}
                    if(cl[j].type == 2){nt_cl1_2 = nt_cl1_2 + 1;}
                }
            }
            
            clock_t toc_full = clock();
            
            double tic_toc_full = (double)(toc_full-tic_full);
            cout << endl;
            //                std::cout << "time step:    " << i << "    N_FILS:    " << N_FILS << "    time taken    " << tic_toc_full/CLOCKS_PER_SEC << "s" << endl;
            printf(  "\n%s%i%s%i %s%12.4g%s\n", "time step:    " , i,",    N_FILS:    ", N_FILS,",    Runtime: " ,tic_toc_full/CLOCKS_PER_SEC ,"sec");
            std::cout << "cl bound: " << cl_nst1 + cl_nst2 << ",    state 1: " << cl_nst1 << ",    state 2 (type 1): " << nt_cl1_1 << ",    state 2 (type 2): " << nt_cl1_2 << ",    total length: " << tot_fil_length << endl;
            cout << "inner pairs: " << pair.n_inner_pair << " outer pairs: " << pair.n_outer_pair << endl;
            cout << endl;
        }
        if (i % t_config == 0)
        {
            // Write the filament configuation details periodically
            ofstream config_write;
            config_write.open ("configuration.txt",ios::trunc);
            config_write << N_FILS << endl;
            if(N_DIMS == 2){for (j=0; j < N_FILS; j++){config_write << mt[j].fil_length << "\t" << mt[j].growth_phase << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1]<< "\n";}}
            else if(N_DIMS == 3){for (j=0; j < N_FILS; j++){config_write << mt[j].fil_length << "\t" << mt[j].growth_phase << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2]<< "\n";}}
            
            config_write.close();
            
            // Write the cross-linker details periodically
            ofstream config_write_cl_end;
            config_write_cl_end.open ("configuration_cl.txt",ios::trunc);
            config_write_cl_end << N_CLS << "\t" << cl_nst1 << "\t" << cl_nst2 << endl;
            for (j=0; j < N_CLS; j++)
            {
                if(cl[j].state == 0)
                {
                    config_write_cl_end << 0 << endl;
                }
                else if(cl[j].state == 1)
                {
                    config_write_cl_end << 1 << "\t" << cl[j].fil_i << "\t" << cl[j].epsilon_i << endl;
                }
                else if(cl[j].state == 2)
                {
                    config_write_cl_end << 2 << "\t" << cl[j].fil_i << "\t" << cl[j].epsilon_i << "\t" << cl[j].fil_j << "\t" << cl[j].epsilon_j << endl;
                }
            }
            config_write_cl_end.close();
        
        }
    }
    // end of main loop
    
    // write configuration and close files
    std::cout << "Number of filaments:   " << N_FILS << endl;
    std::cout << "Number of CLs:   " << N_CLS << "    " << cl_nst1 << "    " << cl_nst2 << endl;
    
    ofstream config_write_end;
    ofstream config_write_cl_end;

    
    // ##############################################################################################################################
    // Final configuration files
    
    config_write_end.open ("configuration.txt",ios::trunc);
    config_write_cl_end.open ("configuration_cl.txt",ios::trunc);
    
    config_write_end << N_FILS << endl;
    if(N_DIMS == 2){for (j=0; j < N_FILS; j++){config_write_end << mt[j].fil_length << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1]<< "\n";}}
    else if(N_DIMS == 3){for (j=0; j < N_FILS; j++){config_write_end << mt[j].fil_length << "\t" << mt[j].growth_phase << "\t" << mt[j].r_com[0] << "\t" << mt[j].r_com[1] << "\t" << mt[j].r_com[2] << "\t" << mt[j].u_vec[0] << "\t" << mt[j].u_vec[1] << "\t" << mt[j].u_vec[2]<< "\n";}}
    config_write_cl_end << N_CLS << "\t" << cl_nst1 << "\t" << cl_nst2 << endl;
    for (j=0; j < N_CLS; j++)
    {
        if(cl[j].state == 0)
        {
            config_write_cl_end << 0 << endl;
        }
        else if(cl[j].state == 1)
        {
            config_write_cl_end << 1 << "\t" << cl[j].fil_i << "\t" << cl[j].epsilon_i << endl;
        }
        else if(cl[j].state == 2)
        {
            config_write_cl_end << 2 << "\t" << cl[j].fil_i << "\t" << cl[j].epsilon_i << "\t" << cl[j].fil_j << "\t" << cl[j].epsilon_j << endl;
        }
    }

    config_write_end.close();
    config_write_cl_end.close();

    // ##############################################################################################################################

    trajectory.close();
    id_trajectory_com.close();

    cl_trajectory.close();
    id_trajectory.close();
    id_orientation.close();
    network_trajectory.close();
    network_orientation.close();
    pinning_stress.close();
    cl_id.close();
    cl_nt.close();
    occ_time.close();
    fil_time.close();
    
    traj_tscale_1.close();
    traj_tscale_2.close();
    traj_tscale_3.close();
    
    clock_t toc_full = clock();
        
    double tic_toc_full = (double)(toc_full-tic_full);
    cout << "\nTime taken: " << tic_toc_full/CLOCKS_PER_SEC << " (s)" << " , " << ((tic_toc_full/CLOCKS_PER_SEC))/3600.0 << " (hr)"<< endl << endl;

    log << "\n--------------------------------------------------------------- " << endl;
    log << "\nEnd of run stats: " << endl;
    log << "\nFinal N_FILS: " << N_FILS  << "\n\ncl bound: " << cl_nst1 + cl_nst2 << "\nstate 1:  " << cl_nst1 << "\nstate 2:  " << cl_nst2  << endl;
    
    log << "\nFinal time taken: " << tic_toc_full/CLOCKS_PER_SEC << " (s)" << " , " << ((tic_toc_full/CLOCKS_PER_SEC))/3600.0 << " (hr)"<< endl << endl;
    log.close();
    
    return 0;
}
