// Benjamin Dalton 09/26/2017

#include <iostream>
#include <fstream>
#include <random>
#include <cstdlib>
#include <math.h>

#include "wall.h"

using namespace std;

wall :: wall(){}

void wall :: wall_init()
{
    x1_n = 0;
    x2_n = 0;
    z1_n = 0;
    z2_n = 0;
    
    wl_x1.resize(0);
    wl_x2.resize(0);
    wl_z1.resize(0);
    wl_z2.resize(0);
}

void wall :: wall_x1_add(unsigned int fil_id)
{
    x1_n = x1_n + 1;
    wl_x1.push_back(fil_id);
}

void wall :: wall_x2_add(unsigned int fil_id)
{
    x2_n = x2_n + 1; wl_x2.push_back(fil_id);
}

void wall :: wall_z1_add(unsigned int fil_id)
{
    z1_n = z1_n + 1; wl_z1.push_back(fil_id);
}

void wall :: wall_z2_add(unsigned int fil_id)
{
    z2_n = z2_n + 1; wl_z2.push_back(fil_id);
}

void wall :: wall_2d_remove(unsigned int fil_id, unsigned int end_id)
{
    if(fil_id == end_id)
    {
        for(int j=0; j < x1_n; j++)
        {
            if(wl_x1[j] == fil_id)
            {
                wl_x1.erase(wl_x1.begin() + j);
                x1_n = x1_n - 1;
                break;
            }
        }
        for(int j=0; j < x2_n; j++)
        {
            if(wl_x2[j] == fil_id)
            {
                wl_x2.erase(wl_x2.begin() + j);
                x2_n = x2_n - 1;
                break;
            }
        }
    }
    else if(fil_id != end_id)
    {
        for(int j=0; j < x1_n; j++)
        {
            if(wl_x1[j] == end_id)
            {
                wl_x1[j] = fil_id;
                break;
            }
        }
        for(int j=0; j < x2_n; j++)
        {
            if(wl_x2[j] == end_id)
            {
                wl_x2[j] = fil_id;
                break;
            }
        }
    }
}

void wall :: wall_3d_remove(unsigned int fil_id, unsigned int end_id)
{

    for(int j=0; j < x1_n; j++)
    {
        if(wl_x1[j] == fil_id)
        {
            wl_x1[j] = fil_id;
            wl_x1.erase(wl_x1.begin() + j);
            x1_n = x1_n - 1;
            break;
        }
    }
    
    for(int j=0; j < x1_n; j++)
    {
        if(wl_x1[j] == end_id)
        {
            wl_x1[j] = fil_id;
            break;
        }
    }

    for(int j=0; j < x2_n; j++)
    {
        if(wl_x2[j] == fil_id)
        {
            wl_x2[j] = fil_id;
            wl_x2.erase(wl_x2.begin() + j);
            x2_n = x2_n - 1;
            break;
        }
    }
    
    for(int j=0; j < x2_n; j++)
    {
        if(wl_x2[j] == end_id)
        {
            wl_x2[j] = fil_id;
            break;
        }
    }
    

    for(int j=0; j < z1_n; j++)
    {
        if(wl_z1[j] == fil_id)
        {
            wl_z1[j] = fil_id;
            wl_z1.erase(wl_z1.begin() + j);
            z1_n = z1_n - 1;
            break;
        }
    }
    
    for(int j=0; j < z1_n; j++)
    {
        if(wl_z1[j] == end_id)
        {
            wl_z1[j] = fil_id;
            break;
        }
    }
    
    for(int j=0; j < z2_n; j++)
    {
        if(wl_z2[j] == fil_id)
        {
            wl_z2[j] = fil_id;
            wl_z2.erase(wl_z2.begin() + j);
            z2_n = z2_n - 1;
            break;
        }
    }
    
    for(int j=0; j < z2_n; j++)
    {
        if(wl_z2[j] == end_id)
        {
            wl_z2[j] = fil_id;
            break;
        }
    }
}