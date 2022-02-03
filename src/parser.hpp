#ifndef PARSER_HPP
#define PARSER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>
#include <sstream>

#include "polygon.hpp"

using namespace std;

/*
	This file defines functions to read and write data from and to files
	Author: Debora Monego | Columbia University | 2021
    data structures:
    - vertices: list of vertices {v0, v1, ... vN} with positions (x0, y0, z0), (x1, y1, z1) ... (xN, yN, zN)

*/

/* void read_file(string filename)
{
	ifstream myfile;
	myfile.open(filename);
} */

/* Description: Read vertex indices for each polygon in the network from file
   Parameters: string with filename with indices - generated in initial configuration
   Returns: network_inidices: vector of vectors with indices for each vertex in each polygon in the network*/
std::vector<std::vector<int> > read_cell_indices(string networkfile);

// Build network from cell indices
std::vector<Polygon> build_network(std::vector<vector<int> > cell_indices, double A0);

// Read vertex coordinates
std::vector<vector<double> > read_vertices(string vertexfile);

// Read indices for edges between vertices
std::vector<vector<int> > read_edges(string edgefile);

// Write matrix to file
template <class myMatrix>
void write_matrix(std::vector<std::vector<myMatrix> > matrix, string outfile);

// Write vector to file
template <class myVector>
void write_vector(std::vector<myVector> vector, string outfile);

#endif