// Utilities library for useful calculations
// 
// Benjamin Dalton 18/11/2015

#ifndef utilities_H
#define utilities_H

#include "filament.h"

using namespace std;

void periodic_shift(int Nt, vector<vector<double>> & r_t, vector<vector<double>> & r_per_t, double L_x, double L_y, double L_z);
void periodic_scale(double tscale, int Nt, vector<vector<double>> & r_per_t);

#endif
