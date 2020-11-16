# Fractal-correlation-RP

The code contained in this repository was used to generate the data presented in the paper "Correlated rigidty percolation in fractal lattices."

Some of the code is currently unpolished in terms of comments. I will provide information on each file and make the comments more instructive later. Also, much of the code was not written by me. I will also try to make this more clear later. 

To compile a program that runs the pebble game algorithm on a lattice of Sierpinski Gaskets, travel to the sg_lattice directory and execute 
g++ -std=c++17 connections.cpp coordinates.cpp main.cpp SiteRP.cpp bond.cpp 

To compile a program that runs the pebble game algorithm on a lattice of triangular plates, travel to triangular_plates and execute 
g++ -std=c++17 bond.cpp main.cpp SiteRP.cpp 

More detail to come.




Note for both programs the following block of text was commented out, which differs from the original pebble game code. It has not seemed to be substantial for my purposes, but this may be important for future users.

// if rcluster_site has not stored the RigidIndex, store it
            /*if (std::find(rcluster_site[site_I].begin(),rcluster_site[site_I].end(),it->RigidIndex) == rcluster_site[site_I].end()){
                rcluster_site[site_I].push_back(it->RigidIndex);
            }
            if (std::find(rcluster_site[site_J].begin(),rcluster_site[site_J].end(),it->RigidIndex) == rcluster_site[site_J].end()){
                rcluster_site[site_J].push_back(it->RigidIndex);
            }*/
            //push the rigid indices to the rcluster_site in order to store the info

