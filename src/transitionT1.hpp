#ifndef TRANSITIONT1_HPP
#define TRANSITIONT1_HPP

#include <math.h>
#include <string.h>
#include <vector>

#include "polygon.hpp"

/*
transitionT1.py - implements T1 transition for short bond lengths

4 polys involved in transition are 1-4 counter-clockwise order
Cells defined such that:
Cell 0: i4, i1, i2, i5
Cell 1: i3, i1, i4
Cell 2: i6, i2, i1, i3
Cell 3: i5, i2, i6

Edges defined such that:
Edge 0: i1 - i2
Edge 1: i1 - i3
Edge 2: i1 - i4
Edge 3: i2 - i1
Edge 4: i2 - i5
Edge 5: i2 - i6
Edge 6: i3 - i1 # reverse edges
Edge 7: i4 - i1
Edge 8: i5 - i2
Edge 9: i6 - i2
*/

// Get vertices associated with the transition
std::vector<int> get_vertex_indices(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids);

// Get cells and edges associated with short bond length
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices);

// Get cells and edges associated with left side transition
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T2_right(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices);

// Get cells and edges associated with rigth side transition
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_right(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices);

// Find 4 cells involved with 2 vertices
// Labeled cells 0 - 3 in counter-clockwise order
// Cell 0 and Cell 3 are neighbors
std::vector<int> get_4_cells(std::vector<Polygon> network, int i1, int i2);

// Perform T1 transition and check the energy change
// std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_transition(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<std::vector<int> > edges,
// 																			   double lx, double ly, double lmin, double ka, double Lambda, double gamma);
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_transition(std::vector<std::vector<double> > vertices, std::vector<Polygon> network,
																												   std::vector<std::vector<int> > edges, double lx, double ly, double lmin,
																												   double ka, double Lambda, double gamma, double ksep);

std::vector<std::vector<double> > set_T1_metastable(int i1, int i2, double lx, double ly, double ksep, std::vector<std::vector<double> > vertices);

// Set new cell indices and edges for T1 left transition
// std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > set_T1_left(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_l_data, std::vector<int> cell_ids,
// 																			 std::vector<std::vector<int> > edges, std::vector<int> vertex_indices);
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_left(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_l_data,
																												 std::vector<int> cell_ids, std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, double lx, double ly, double lmin, double ksep,
																												 std::vector<std::vector<double> > vertices);

// Set new cell indices and edges for T1 right transition
// std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > set_T1_right(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_r_data, std::vector<int> cell_ids,
// 																			  std::vector<std::vector<int> > edges, std::vector<int> vertex_indices);
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_right(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_r_data, std::vector<int> cell_ids,
																												  std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, double lmin, double lx, double ly, double ksep,
																												  std::vector<std::vector<double> > vertices);
#endif