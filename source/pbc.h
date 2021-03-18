// Periodic Boundary Conditions and minimal image calculation
//
// Benjamin Dalton 07/07/2016

#ifndef pbc_H
#define pbc_H

using namespace std;

class pbc
{
    
    private:
    public:
    
    vector<double> box_vec,ri_image,rj_image;
    
    pbc();
    void pbc_constructor(unsigned int n_dims, double L_x, double L_y, double L_z);
    void minimal_image(vector<double>& ri_com,vector<double>& rj_com,int N_dims);
    void pbc_calc(vector<double>& r_com,int n_dims, double pbc_index[3]);
    void pbc_y_calc(vector<double>& r_com,int n_dims, double pbc_index[3]);
};

#endif