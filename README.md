# Fractal-correlation-RP

The code contained in this repository was used to generate the data presented in the paper "Correlated rigidty percolation in fractal lattices."  

To compile a program that runs the pebble game algorithm on a lattice of Sierpinski Gaskets, travel to the sg_lattice directory and execute  
g++ -std=c++17 connections.cpp coordinates.cpp main.cpp SiteRP.cpp bond.cpp   

To compile a program that runs the pebble game algorithm on a lattice of upwards pointing triangular plates, travel to triangular_plates and execute  
g++ -std=c++17 bond.cpp main.cpp SiteRP.cpp  

Understanding output for SG RP:  
pc_and_phi_n2s7.txt is a two column text document which contains the critical packing and volume fractions for the lattice (the lattice of the code as is arbitrarily has  n=2 s=7)  
rclusterout.txt is a list of all the sites which belong to the largest rigid cluster  
mathematica_coords.txt is a list (in Mathematica style) of all the coordinates of the sites of the lattice, unoccupied or not  
mathematica_occsites_from_pebble.txt is a list of all the vertex numbers which are occupied  
mathematica_lines_from_pebble.txt is a two column list of all the sites which are connected to each other in the lattice  
mathematica_lines.txt is a two column list of all the sites which may in principle share a bond, whether or not they actually do from the pebble game  
        
Understanding output for triangle plate RP:  
n0s32macro_data.txt is a list of the critical packing fractions determined for a lattice of n=0 (i.e., a regular triangular lattice) and the arbitrarily chosen s=32.  
rclusterout.txt is again, a list of all the sites which belong to the largest rigid cluster
occsites_from_pebble.txt is a list of all the vertex numbers which are occupied  
bonds_from_pebble.txt is a two column list of all the sites which are connected to each other in the lattice  
