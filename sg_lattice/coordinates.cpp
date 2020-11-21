//
//  Created by Shae Machlus on 6/20/19.
//
//////////////////////////////////////////
//SIERPENSKI GASKET LATTICE COORDINATE GENERATOR
//////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "sg_connections.h"

using namespace std;

void sg_coords(int n, int s)
{

    int len = EdgeVerts(n, s)/s;
    
    ofstream coordfile;
    
    vector<double> x(3);
    vector<double> y(3);
    
    x[0] = 0.0; y[0] = 0.0;
    x[1] = 0.5; y[1] = sqrt(3)/2;
    x[2] = 1.0; y[2] = 0.0;
    
    int sizex = int(x.size());
    for(int k = 1; k <= n; k++)
    {
        for(int i = 0; i < sizex; i++)
        {
            x.push_back(x[i] + pow(2,k-1) * 0.5);
            y.push_back(y[i] + pow(2,k-1) * sqrt(3)/2);
            
            x.push_back(x[i] + pow(2,k-1) * 1);
            y.push_back(y[i]);
        }
        sizex = int(x.size());
    }
    
    //x triangles
    for(int l = 1; l < s; l++)
    {
        for(int i = 0; i < sizex; i++)
        {
            x.push_back(x[i] + l * len);
            y.push_back(y[i]);
        }
    }
    sizex = int(x.size());
    
    //y triangles
    for(int l = 1; l < s; l++)
    {
        for(int i = 0; i < sizex; i++)
        {
            x.push_back(x[i] + l * len/2);
            y.push_back(y[i] - len * l * sqrt(3)/2);
        }
    }
    
    //Output
    coordfile.open("mathematica_coords.txt");
    coordfile << "{";
    for(int j = 0; j <= x.size() - 1; j++)
    {
        if(j == x.size() - 1)
        {
            coordfile << "{" << x[j] << ", " << y[j] << "}}" << endl;
        }
        else
        {
            coordfile << "{" << x[j] << ", " << y[j] << "}, ";
        }
    }
    coordfile.close();
    
    return;
}
