#ifndef FORCE_HPP
#define FORCE_HPP

#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

using namespace std;

// Get positions of the vertices in clockwise order
std::vector<double> get_clockwise(int index, std::vector<int> indices, std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);

// Get positions of the vertices in counter clockwise order
std::vector<double> get_counter_clockwise(int index, std::vector<int> indices, std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);

// Calculate force on vertex due to elasticity
std::vector<std::vector<double> > calc_force_elasticity(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double ka, std::vector<double> L, std::vector<std::vector<int> > edges);

// Calculate force due to contraction
std::vector<std::vector<double> > calc_force_contraction(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double gamma, std::vector<double> L, std::vector<std::vector<int> > edges);

// Calculate force due to adhesion
std::vector<std::vector<double> > calc_force_adhesion(std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, double Lambda, std::vector<double> L);

// Force to move vertices of polygons in a particular direction
std::vector<std::vector<double> > calc_force_motility(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double eta, double xi);

// Computes forces in the current configuration of the vertex model
std::vector<std::vector<double> > get_forces(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<std::vector<int> > edges,
                                             std::vector<double> L, double ka, double Lambda, double gamma, double eta, double xi);

// Employ forces to move vertices' positions
std::pair<std::vector<std::vector<double> >, std::vector<std::vector<int> > >  move_vertices(std::vector<std::vector<double> > vertices, std::vector<std::vector<double> > forces, std::vector<std::vector<int> > edges,
                                                                              std::vector<double> L, double delta_t, double lmin);

#endif