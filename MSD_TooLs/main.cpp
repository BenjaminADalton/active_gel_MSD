
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
#include "utilities.h"

using namespace std;

#define N_DIMS 3            // Number of spatial dimensions (2 or 3)

int main()
{
    // number of filaments and crosslinkers
    int N_FILS = 500;
    int N_Sets = 1;
    
    int N_Total = N_FILS*N_Sets;
    
    // box dimensions
    double L_x = 5.0e-6;
    double L_y = 5.0e-6;
    double L_z = 0.25e-6;
    
    // real time of each iteration (in seconds)
    double dt = 5.0e-6;
    
    // convert to micro-meters
    double tscale = 1.0e6;
    
    int t_inc_1 = 1;
    int t_inc_2 = 100;
    int t_inc_3 = 10000;
    
    double dt_1 = double(t_inc_1)*dt;
    double dt_2 = double(t_inc_2)*dt;
    double dt_3 = double(t_inc_3)*dt;

    int Nt_1 = 2000;   // tscale_1_upper/tscale_1_inc = n1_plot
    int Nt_2 = 5000;   // tscale_1_upper/tscale_1_inc = n1_plot
    int Nt_3 = 50;   // tscale_1_upper/tscale_1_inc = n1_plot
    
    int pad_length = 2;
    
    // various parameters used and reused throughout
    int i,j,n,t,time,nfil_in;
    
    // declare the filament and cross-linker objects
    vector<filament> mt(0);
    
    cout << endl << "Reading: " << N_Sets << " sets" << endl << endl;
    
    string set_number, padded_str, read_in_str;

    for(i = 0; i < N_Sets; i++)
    {
        // construct the filament classes
        for(j = 0; j < N_FILS; j++)
        {
            filament mt_new;
            mt_new.filament_init(N_DIMS,Nt_1,Nt_2,Nt_3);
            mt.push_back(mt_new);
        }
        
        // ####################### open trajectories files #######################
        set_number  = to_string(i+1);
        padded_str  = string(pad_length - set_number.length(), '0') + set_number;
        read_in_str = "../source/traj_tscale_1.txt";
        cout << "Reading: " + read_in_str << endl;
        
        ifstream t_read_1; t_read_1.open(read_in_str, ios::in);
        
        for(t = 0; t < Nt_1; t++)
        {
            t_read_1 >> time >> nfil_in;
            
//            cout << double(t)*dt_1 << endl;
            
            for(n = 0; n < N_FILS; n++)
            {
                mt[i*N_FILS + n].t1[t] = double(t)*dt_1;
                
                t_read_1 >> mt[i*N_FILS + n].l_t1[t];
                t_read_1 >> mt[i*N_FILS + n].r_t1[0][t]; t_read_1 >> mt[i*N_FILS + n].r_t1[1][t]; t_read_1 >> mt[i*N_FILS + n].r_t1[2][t];
                t_read_1 >> mt[i*N_FILS + n].u_t1[0][t]; t_read_1 >> mt[i*N_FILS + n].u_t1[1][t]; t_read_1 >> mt[i*N_FILS + n].u_t1[2][t];
            }
        }
        
        // ####################### open trajectories files #######################
        set_number  = to_string(i+1);
        padded_str  = string(pad_length - set_number.length(), '0') + set_number;
        read_in_str = "../source/traj_tscale_2.txt";
        cout << "Reading: " + read_in_str << endl;
        
        ifstream t_read_2; t_read_2.open(read_in_str, ios::in);
        
        for(t = 0; t < Nt_2; t++)
        {
            t_read_2 >> time >> nfil_in;
            
//            cout << double(t)*dt_2 << endl;
            
            for(n = 0; n < N_FILS; n++)
            {
                mt[i*N_FILS + n].t2[t] = double(t)*dt_2;
                
                t_read_2 >> mt[i*N_FILS + n].l_t2[t];
                t_read_2 >> mt[i*N_FILS + n].r_t2[0][t]; t_read_2 >> mt[i*N_FILS + n].r_t2[1][t]; t_read_2 >> mt[i*N_FILS + n].r_t2[2][t];
                t_read_2 >> mt[i*N_FILS + n].u_t2[0][t]; t_read_2 >> mt[i*N_FILS + n].u_t2[1][t]; t_read_2 >> mt[i*N_FILS + n].u_t2[2][t];
            }
        }
        
        // ####################### open trajectories files #######################
        set_number  = to_string(i+1);
        padded_str  = string(pad_length - set_number.length(), '0') + set_number;
        read_in_str = "../source/traj_tscale_3.txt";
        cout << "Reading: " + read_in_str << endl;
        
        ifstream t_read_3; t_read_3.open(read_in_str, ios::in);
        
        for(t = 0; t < Nt_3; t++)
        {
            t_read_3 >> time >> nfil_in;
            
//            cout << double(t)*dt_3 << endl;
            
            for(n = 0; n < N_FILS; n++)
            {
                mt[i*N_FILS + n].t3[t] = double(t)*dt_3;
                
                t_read_3 >> mt[i*N_FILS + n].l_t3[t];
                t_read_3 >> mt[i*N_FILS + n].r_t3[0][t]; t_read_3 >> mt[i*N_FILS + n].r_t3[1][t]; t_read_3 >> mt[i*N_FILS + n].r_t3[2][t];
                t_read_3 >> mt[i*N_FILS + n].u_t3[0][t]; t_read_3 >> mt[i*N_FILS + n].u_t3[1][t]; t_read_3 >> mt[i*N_FILS + n].u_t3[2][t];
            }
        }
    }
        
    // Calculate the periodic corrections for all stored trajectories
    for(n = 0; n < N_Total; n++)
    {
        periodic_shift(Nt_1, mt[n].r_t1, mt[n].r_per_t1, L_x, L_y, L_z);
        periodic_shift(Nt_2, mt[n].r_t2, mt[n].r_per_t2, L_x, L_y, L_z);
        periodic_shift(Nt_3, mt[n].r_t3, mt[n].r_per_t3, L_x, L_y, L_z);
    }
    
    // Calculate the periodic corrections for all stored trajectories
    for(n = 0; n < N_Total; n++)
    {
        periodic_scale(tscale, Nt_1, mt[n].r_per_t1);
        periodic_scale(tscale, Nt_2, mt[n].r_per_t2);
        periodic_scale(tscale, Nt_3, mt[n].r_per_t3);
    }

    
    // Calculate MSD for full set of filaments over three time scales
    vector<double> msd_1, msd_2, msd_3;
    vector<vector<double>> dr_sqr_1, dr_sqr_2, dr_sqr_3;

    msd_1.resize(Nt_1);
    msd_2.resize(Nt_2);
    msd_3.resize(Nt_3);
    dr_sqr_1.resize(N_Total);
    dr_sqr_2.resize(N_Total);
    dr_sqr_3.resize(N_Total);
    
    for(n=0; n < N_Total; n++){dr_sqr_1[n].resize(Nt_1); dr_sqr_2[n].resize(Nt_2); dr_sqr_3[n].resize(Nt_3);}
    
    for(n=0; n < N_Total; n++)
    {
        for(t = 0; t < Nt_1; t++)
        {
            dr_sqr_1[n][t] = pow(mt[n].r_per_t1[0][t] - mt[n].r_per_t1[0][0],2) + pow(mt[n].r_per_t1[1][t] - mt[n].r_per_t1[1][0],2);
            // dr_sqr_1[n][t] = 1e0*dr_sqr_1[n][t];
        }
        for(t = 0; t < Nt_2; t++)
        {
            dr_sqr_2[n][t] = pow(mt[n].r_per_t2[0][t] - mt[n].r_per_t2[0][0],2) + pow(mt[n].r_per_t2[1][t] - mt[n].r_per_t2[1][0],2);
            // dr_sqr_2[n][t] = 1e0*dr_sqr_2[n][t];
        }
        for(t = 0; t < Nt_3; t++)
        {
            dr_sqr_3[n][t] = pow(mt[n].r_per_t3[0][t] - mt[n].r_per_t3[0][0],2) + pow(mt[n].r_per_t3[1][t] - mt[n].r_per_t3[1][0],2);
            // dr_sqr_3[n][t] = 1e0*dr_sqr_3[n][t];
        }
    }
    
    for(t = 0; t < Nt_1; t++){for(n=0; n < N_Total; n++){msd_1[t] = msd_1[t] + dr_sqr_1[n][t];} msd_1[t] = msd_1[t]/N_Total;}
    for(t = 0; t < Nt_2; t++){for(n=0; n < N_Total; n++){msd_2[t] = msd_2[t] + dr_sqr_2[n][t];} msd_2[t] = msd_2[t]/N_Total;}
    for(t = 0; t < Nt_3; t++){for(n=0; n < N_Total; n++){msd_3[t] = msd_3[t] + dr_sqr_3[n][t];} msd_3[t] = msd_3[t]/N_Total;}
    
    cout << endl << "Total number of filaments: " << N_Total << endl << endl;
    // Calculate the periodic corrections for all trajectories
    
    // open text files to be used for data writing
    ofstream msd_t1; msd_t1.open("../msd_t1.txt");
    ofstream msd_t2; msd_t2.open("../msd_t2.txt");
    ofstream msd_t3; msd_t3.open("../msd_t3.txt");
    
    for(t = 0; t < Nt_1; t++){msd_t1 << mt[0].t1[t] << "\t" << msd_1[t] << endl;}
    for(t = 0; t < Nt_2; t++){msd_t2 << mt[0].t2[t] << "\t" << msd_2[t] << endl;}
    for(t = 0; t < Nt_3; t++){msd_t3 << mt[0].t3[t] << "\t" << msd_3[t] << endl;}
    
    msd_t1.close();
    msd_t2.close();
    msd_t3.close();

    return 0;
}
