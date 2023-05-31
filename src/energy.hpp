#ifndef ENERGY_HPP
#define ENERGY_HPP

#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

#include "polygon.hpp"

using namespace std;

/*  Calculate elasticity energy
        data structures: vertices: vector of vectors with the coordinates of the polygon's vertices
                         network: vector with the indices for each cell in the network
                         ka: double with the modulus constant of the cell - resistance to area changes
                         L: vector with the simulation box dimensions
        returns: double with the elastic energy */
double get_energy_elasticity(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double ka, std::vector<double> L, std::vector<std::vector<int> > edges);

/* Calculate the adhesion energy - anisotropic Vertex Model
        data structures: vertices: vector of vectors with the coordinates of the polygon's vertices
                         edges: vector with a pair of indices for each vertex of the edge
                         Lambda: double with the adhesion energy per unit lenght
                         L: vector with the simulation box dimensions
        returns: double with the adhesion energy */
double get_energy_adhesion(std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, double Lambda, std::vector<double> L);

/* Calculate the contraction energy
        data structures: vertices: vector of vectors with the coordinates of the polygon's vertices
                         network: vector with the indices for each cell in the network
                         gamma: double with the elastic constant of the cell
                         L: vector with the simulation box dimensions
        returns: double with the contractility energy */
double get_energy_contraction(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double gamma, std::vector<double> L, std::vector<std::vector<int> > edges);

// Get total energy
double get_total_energy(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<std::vector<int> > edges, double ka, std::vector<double> L,
                        double Lambda, double gamma, double lambda_potts);

double get_energy_perimeter_potts(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<std::vector<int> > edges, double lambda_potts, std::vector<double> L);

double get_energy_adhesion_potts(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<double> L, std::vector<std::vector<int> > edges);

#endif