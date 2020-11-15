//
// Created by Shae Machlus on 6/26/19
//

#ifndef SG_H
#define SG_H

//Returns the total number of vertices in the latice
int TotalVerts(int N, int S);

//Returns the number of vertices along the bottom edge of the unit cell
int EdgeVerts(int N, int S);

//Checks to make sure the number of vertices is equal to the size of the lattice
void SizeCheck(std::vector<int> VERTICES, int N);

//
void SG_VertList(std::vector<int> &VERTICES, int N);

//Expands the list horizontally (right) to include a whole strip of the lattice
void SG_StripVertList(std::vector<int> &VERTICES, int S);

//Expands the list vertically (down) to include the entire lattice
void SG_LatticeVertList(std::vector<int> &VERTICES, int S);

//The element slots represent vertices, and contain the layer number at which that vertex can be found
void FillLayers(std::vector<int> &VERTICES, std::vector<int> &LAYERS);

//Determines what type of vertex dum is and records it in the vector type under dum's index
void TypeAssign(std::vector <int> &TYPE, int DUM, int TEST_LAYER, int SINGLE_TRI_INDEX, int VERT_COUNT, int VERT_COUNT_BELOW, int LL, int S);

//Returns how many vertices far right into a layer dum is
int Depth(std::vector<int> &LAYERS, int DUM, int DUM_LAYER, int SIZE);

//Defines important quantities and assigns neighbors
int DoEverything(int DUM, int N, int S, int SIZE, int D);

#endif
