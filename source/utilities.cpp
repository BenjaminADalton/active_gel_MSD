// Filament class. Includes nucleation if a MT reduces to zero length.
//
// Benjamin Dalton 08/10/2015

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>
#include <chrono>
#include <ctime>

using namespace std;

#include "utilities.h"
#include "filament.h"


time_t now    = time(0);
tm *ltm       = localtime(&now);
int srand_int = (ltm->tm_hour)*(ltm->tm_sec) + (ltm->tm_min)*(ltm->tm_sec);

mt19937 generator_uni(2*srand_int+ srand_int - 1);
mt19937 generator(srand_int);

uniform_real_distribution<double> dis_uni(0.0, 1.0);
uniform_real_distribution<double> dis_mer(0.0, 1.0);

void gram_schmidt_2d(vector<double>& u_vec, vector<vector<double>>& u_perp_2d, double e1[2])
{
    double e_dot_u = e1[0]*u_vec[0] + e1[1]*u_vec[1];
    
    double vx = e1[0] - e_dot_u*u_vec[0];
    double vy = e1[1] - e_dot_u*u_vec[1];

    double v_mag = sqrt(vx*vx + vy*vy);
    
    u_perp_2d[0][0] = vx/v_mag;
    u_perp_2d[0][1] = vy/v_mag;

}

void gram_schmidt_3d(vector<double>& u_vec, vector<vector<double>>& u_perp_3d, double e1[3], double e2[3])
{
    double e1_dot_u = e1[0]*u_vec[0] + e1[1]*u_vec[1] + e1[2]*u_vec[2];
    double vx = e1[0] - e1_dot_u*u_vec[0];
    double vy = e1[1] - e1_dot_u*u_vec[1];
    double vz = e1[2] - e1_dot_u*u_vec[2];
    
    double v_mag = sqrt(vx*vx + vy*vy + vz*vz);
    
    u_perp_3d[0][0]= vx/v_mag;
    u_perp_3d[0][1] = vy/v_mag;
    u_perp_3d[0][2] = vz/v_mag;
    
    u_perp_3d[1][0] = u_vec[1]*u_perp_3d[0][2]-u_vec[2]*u_perp_3d[0][1];
    u_perp_3d[1][1] = -u_vec[0]*u_perp_3d[0][2]+u_vec[2]*u_perp_3d[0][0];
    u_perp_3d[1][2] = u_vec[0]*u_perp_3d[0][1]-u_vec[1]*u_perp_3d[0][0];
    
}

void gram_schmidt_2d_2(vector<double>& u_vec, vector<double>& u_perp_2d, double e1[2])
{
    double e_dot_u = e1[0]*u_vec[0] + e1[1]*u_vec[1];
    double vx = e1[0] - e_dot_u*u_vec[0];
    double vy = e1[1] - e_dot_u*u_vec[1];
    
    double v_mag = sqrt(vx*vx + vy*vy);
    
    u_perp_2d[0] = vx/v_mag;
    u_perp_2d[1] = vy/v_mag;
    
}

void gram_schmidt_3d_2(vector<double>& u_vec, vector<double>& u_perp_3d_1, vector<double>& u_perp_3d_2, double e1[3], double e2[3])
{
    double e1_dot_u = e1[0]*u_vec[0] + e1[1]*u_vec[1] + e1[2]*u_vec[2];
    double vx = e1[0] - e1_dot_u*u_vec[0];
    double vy = e1[1] - e1_dot_u*u_vec[1];
    double vz = e1[2] - e1_dot_u*u_vec[2];
    
    double v_mag = sqrt(vx*vx + vy*vy + vz*vz);
    
    u_perp_3d_1[0]= vx/v_mag;
    u_perp_3d_1[1] = vy/v_mag;
    u_perp_3d_1[2] = vz/v_mag;
    
    u_perp_3d_2[0] = u_vec[1]*u_perp_3d_1[2]-u_vec[2]*u_perp_3d_1[1];
    u_perp_3d_2[1] = -u_vec[0]*u_perp_3d_1[2]+u_vec[2]*u_perp_3d_1[0];
    u_perp_3d_2[2] = u_vec[0]*u_perp_3d_1[1]-u_vec[1]*u_perp_3d_1[0];
    
}

void filament_drag_2d_util(vector<double>& u_vec, vector<double>& r_eqnm, vector<double>& u_eqnm, double gamma_array[2][2], double gamma_para, double gamma_perp, double gamma_rot)
{
    gamma_array[0][0] = 0.0;
    gamma_array[0][1] = 0.0;
    gamma_array[1][0] = 0.0;
    gamma_array[1][1] = 0.0;
    
    double eqnm_pass_1 = 0.0;
    double eqnm_pass_2 = 0.0;
    
    gamma_array[0][0] = u_vec[0]*u_vec[0]/gamma_para + (1-u_vec[0]*u_vec[0])/gamma_perp;
    gamma_array[0][1] = u_vec[0]*u_vec[1]/gamma_para - u_vec[0]*u_vec[1]/gamma_perp;
    gamma_array[1][0] = gamma_array[0][1];
    gamma_array[1][1] = u_vec[1]*u_vec[1]/gamma_para + (1-u_vec[1]*u_vec[1])/gamma_perp;
    
    eqnm_pass_1 = gamma_array[0][0]*r_eqnm[0] + gamma_array[0][1]*r_eqnm[1];
    eqnm_pass_2 = gamma_array[1][0]*r_eqnm[0] + gamma_array[1][1]*r_eqnm[1];
    
    r_eqnm[0] = eqnm_pass_1;
    r_eqnm[1] = eqnm_pass_2;
    
    u_eqnm[0] = u_eqnm[0]/gamma_rot;
    u_eqnm[1] = u_eqnm[1]/gamma_rot;
    
}

void filament_drag_3d_util(vector<double>& u_vec, vector<double>& r_eqnm, vector<double>& u_eqnm, double gamma_array[3][3], double gamma_para, double gamma_perp, double gamma_rot)
{
    gamma_array[0][0] = 0.0;
    gamma_array[0][1] = 0.0;
    gamma_array[0][2] = 0.0;
    gamma_array[1][0] = 0.0;
    gamma_array[1][1] = 0.0;
    gamma_array[1][2] = 0.0;
    gamma_array[2][0] = 0.0;
    gamma_array[2][1] = 0.0;
    gamma_array[2][2] = 0.0;
    
    double eqnm_pass_1 = 0.0;
    double eqnm_pass_2 = 0.0;
    double eqnm_pass_3 = 0.0;
    
    gamma_array[0][0] = u_vec[0]*u_vec[0]/gamma_para + (1-u_vec[0]*u_vec[0])/gamma_perp;
    gamma_array[0][1] = u_vec[0]*u_vec[1]/gamma_para - u_vec[0]*u_vec[1]/gamma_perp;
    gamma_array[0][2] = u_vec[0]*u_vec[2]/gamma_para - u_vec[0]*u_vec[2]/gamma_perp;

    gamma_array[1][0] = gamma_array[0][1];
    gamma_array[1][1] = u_vec[1]*u_vec[1]/gamma_para + (1-u_vec[1]*u_vec[1])/gamma_perp;
    gamma_array[1][2] = u_vec[1]*u_vec[2]/gamma_para - u_vec[1]*u_vec[2]/gamma_perp;

    gamma_array[2][0] = gamma_array[0][2];
    gamma_array[2][1] = gamma_array[1][2];
    gamma_array[2][2] = u_vec[2]*u_vec[2]/gamma_para + (1-u_vec[2]*u_vec[2])/gamma_perp;
    
    eqnm_pass_1 = gamma_array[0][0]*r_eqnm[0] + gamma_array[0][1]*r_eqnm[1] + gamma_array[0][2]*r_eqnm[2];
    eqnm_pass_2 = gamma_array[1][0]*r_eqnm[0] + gamma_array[1][1]*r_eqnm[1] + gamma_array[1][2]*r_eqnm[2];
    eqnm_pass_3 = gamma_array[2][0]*r_eqnm[0] + gamma_array[2][1]*r_eqnm[1] + gamma_array[2][2]*r_eqnm[2];
    
    r_eqnm[0] = eqnm_pass_1;
    r_eqnm[1] = eqnm_pass_2;
    r_eqnm[2] = eqnm_pass_3;

    u_eqnm[0] = u_eqnm[0]/gamma_rot;
    u_eqnm[1] = u_eqnm[1]/gamma_rot;
    u_eqnm[2] = u_eqnm[2]/gamma_rot;
    
}

void stoch_2d_util(filament& mt, const double kTt_2)
{

    double rand_para = rand_normal(0,sqrt(kTt_2/mt.gamma_para));
    double rand_perp = rand_normal(0,sqrt(kTt_2/mt.gamma_perp));
    double rand_rot = rand_normal(0,sqrt(kTt_2/mt.gamma_rot));

    mt.r_com[0] = mt.r_com[0] + rand_para*mt.u_vec[0] + rand_perp*mt.u_perp[0][0];
    mt.r_com[1] = mt.r_com[1] + rand_para*mt.u_vec[1] + rand_perp*mt.u_perp[0][1];

    mt.u_vec[0] = mt.u_vec[0] + rand_rot*mt.u_perp[0][0];
    mt.u_vec[1] = mt.u_vec[1] + rand_rot*mt.u_perp[0][1];

    double u_mod = sqrt(mt.u_vec[0]*mt.u_vec[0] + mt.u_vec[1]*mt.u_vec[1]);
    mt.u_vec[0] = mt.u_vec[0]/u_mod;
    mt.u_vec[1] = mt.u_vec[1]/u_mod;

}

void stoch_3d_util(filament& mt, const double kTt_2)
{
    double rand_para   = rand_normal(0,sqrt(kTt_2/mt.gamma_para));
    double rand_perp_1 = rand_normal(0,sqrt(kTt_2/mt.gamma_perp));
    double rand_perp_2 = rand_normal(0,sqrt(kTt_2/mt.gamma_perp));
    double rand_rot_1  = rand_normal(0,sqrt(kTt_2/mt.gamma_rot));
    double rand_rot_2  = rand_normal(0,sqrt(kTt_2/mt.gamma_rot));
    
    mt.r_com[0] = mt.r_com[0] + rand_para*mt.u_vec[0] + rand_perp_1*mt.u_perp[0][0] + rand_perp_2*mt.u_perp[1][0];
    mt.r_com[1] = mt.r_com[1] + rand_para*mt.u_vec[1] + rand_perp_1*mt.u_perp[0][1] + rand_perp_2*mt.u_perp[1][1];
    mt.r_com[2] = mt.r_com[2] + rand_para*mt.u_vec[2] + rand_perp_1*mt.u_perp[0][2] + rand_perp_2*mt.u_perp[1][2];
    
    mt.u_vec[0] = mt.u_vec[0] + rand_rot_1*mt.u_perp[0][0] + rand_rot_2*mt.u_perp[1][0];
    mt.u_vec[1] = mt.u_vec[1] + rand_rot_1*mt.u_perp[0][1] + rand_rot_2*mt.u_perp[1][1];
    mt.u_vec[2] = mt.u_vec[2] + rand_rot_1*mt.u_perp[0][2] + rand_rot_2*mt.u_perp[1][2];
    
}




double uniform_random(int index, int seed)
{
    double rnd_no = 0;

    unsigned int rand_max = 1e6;
    
    return (rand() % rand_max)/float(rand_max);
    
}

double uniform_random_Mersenne()
{
    return dis_uni(generator_uni);
}

// double normal_random(double mean, double std_dev, mt19937 generator_util_in)
double normal_random(double mean, double std_dev)
{

    normal_distribution<double> dist_norm(mean,std_dev);
    return dist_norm(generator);
    
}

// rand_normal - Box Muller method
//
// Takes two uniformly distributed random numbers (U1,U2) as inputs and produces
// two normally distributed numbers (n1,n2) as outputs.
//
// n1 = r*cos(2piU1), n2 = r*cos(2piU1), r = sqrt(-2ln(U2))
//
// only one normal random number as output so the second in stored and use
// on next call
//
double rand_normal(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {

            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            
            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double rand_normal_Mersenne(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*dis_mer(generator) - 1;
            y = 2.0*dis_mer(generator) - 1;
            
            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}


