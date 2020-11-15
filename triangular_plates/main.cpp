//
// Created by Shae Machlus on 7/26/2019
//

#include "siteRP.h"
#include "bond.h"
//#include <time.h>

int main()
{
    remove("bonds_from_pebble.txt");
    remove("occsites_from_pebble.txt");
    remove("rclusterout.txt");
    int seed = 50;
    //srand(int(time(NULL)));
    srand(seed);
    SiteRP a;
    a.OneTrialTest(0.0,1);
    a.listalledges();
    return 0;
}
