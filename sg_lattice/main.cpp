//
// Created by Shae Machlus on 6/26/2019
//
/////////////////////////////////////////
//Main file for running the pebble game
//on a lattice of SG's 
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
    // clear old output     
    remove("mathematica_occsites_from_pebble.txt");
    remove("mathematica_lines_from_pebble.txt");
    remove("mathematica_coords.txt");
    remove("rclusterout.txt");

    // the lattice is randomly seeded
    srand(time(NULL));

    // SiteRP object is created
    SiteRP a; 

    //a.OneTrialTest(0, 1);
    a.PlotNetworkTest();
    sg_coords(a.n,a.s);
    return 0;
}
