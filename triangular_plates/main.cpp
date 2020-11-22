//
// Created by Shae Machlus on 7/26/2019
//

#include "siteRP.h"
#include "bond.h"
#include <time.h>

int main()
{
    // clear old output 
    remove("bonds_from_pebble.txt");
    remove("occsites_from_pebble.txt");
    remove("rclusterout.txt");

    // the lattice is randomly seeded
    srand(int(time(NULL)));
    
    // SiteRP object is created 
    SiteRP a;
    
    a.OneTrialTest(0.0,1);
    //a.PlotNetworkTest();

    a.listalledges();

    return 0;
}
