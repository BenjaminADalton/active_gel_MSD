#ifndef filament_H
#define filament_H

#include <random>

using namespace std;

class filament
{
    private:
    public:
    
    filament();

    vector<double> t1, t2, t3, l_t1, l_t2, l_t3;
    
    vector<vector<double>> r_t1,r_t2,r_t3,u_t1,u_t2,u_t3;
    vector<vector<double>> r_per_t1,r_per_t2,r_per_t3,u_per_t1,u_per_t2,u_per_t3;
    
    void filament_init(int n_dim,int Nt_1,int Nt_2,int Nt_3);
    
};

#endif
