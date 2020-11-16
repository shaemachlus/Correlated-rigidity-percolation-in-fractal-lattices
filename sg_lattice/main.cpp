//
// Created by Shae Machlus on 6/26/2019
//
/////////////////////////////////////////
//Input File for SG Lattice Pebble Game//
/////////////////////////////////////////
//

#include "siteRP.h"
#include "sg_connections.h"
#include "bond.h"
#include "sg_coordinates.h"
#include <time.h>

using namespace std;

int main()
{
    remove("mathematica_occsites_from_pebble.txt");
    remove("mathematica_lines_from_pebble.txt");
    remove("mathematica_coords.txt");
    remove("rclusterout.txt");
    int seed;
    srand(time(NULL));
    SiteRP a;
    a.OneTrialTest(0, 1);
    //a.PlotNetworkTest();
    sg_coords(a.n,a.s);
    return 0;
}
 
/*int main()
{
    remove("mathematica_lines.txt");
    
    //The number of iterations of the fractal
    int n = 1;
    
    //The side length of the unit cell in units of big triangles
    int s = 2;
    
    //The neighbor you want
    int d = 3;
    int nb;
    
    //The number of vertices in the lattice
    int tvc = TotalVerts(n, s);
    

    for(int dum = 0; dum < tvc; dum++)
    {
        nb = DoEverything(dum, n, s, tvc, d);
    }

    sg_coords(n, s);
    
    return 0;
}*/
