//
//  one_open_sg.hpp
//  connections
//
//  Created by Shae Machlus on 7/5/19.
//  Copyright © 2019 Shae Machlus. All rights reserved.
//

#ifndef one_open_sg_hpp
#define one_open_sg_hpp
#include <vector>

//Returns the total number of vertices in the latice
int TotalVerts_osg(int N);

//Returns the number of vertices along the bottom edge of the unit cell
int EdgeVerts_osg(int N);

void SG_VertList_osg(std::vector<int> &VERTICES, int N);

//The element slots represent vertices, and contain the layer number at which that vertex can be found
void FillLayers_osg(std::vector<int> &VERTICES, std::vector<int> &LAYERS);

//Determines what type of vertex dum is and records it in the vector type under dum's index
void TypeAssign_osg(std::vector <int> &TYPE, int DUM, int TEST_LAYER, int SINGLE_TRI_INDEX, int VERT_COUNT, int VERT_COUNT_BELOW, int LL);

//Returns how many vertices far right into a layer dum is
int Depth_osg(std::vector<int> &LAYERS, int DUM, int DUM_LAYER, int SIZE);

//Defines important quantities and assigns neighbors
int DoEverything_osg(int DUM, int N, int SIZE, int D);

#endif
