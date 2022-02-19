#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <math.h>
#include <string>
#include <vector>
#include <algorithm>

//#include "vector.hpp"
using namespace std;

/* Calculates the distance of a point in respect to the periodic boundaries
    / data structures: v1 and v2 are the positions of vertices 1 and 2
    /                  L is a vector with the dimensions of the simulation box
    / returns: vector of doubles */
std::vector<double> pbc_diff(std::vector<double> v1, std::vector<double> v2, std::vector<double> L);

/* Calculates the area of the base of a polygon
    / Considers the base to be on the xy plane - To do: does this need to be generalized?
    / returns: double */
double get_area(std::vector<std::vector<double> > vertices);

/* Calculates the perimeter of a polygon
    / data structures: vertices: vector of vectors with the positions of the polygon's vertices
    / returns: double */
double get_perimeter(std::vector<std::vector<double> > vertices);

/* Calculates the Euclidian distance between two points with coordinates (x0,y0,z0) and (x1,y1,z1)
    / returns: double */
double get_euclidian_distance(double x0, double y0, double x1, double y1);

/* Sets a new separation between two vertices that went through a T1 transition according to the separation ratio ksep
    / returns: vector with new x0 and y0 coordinates */
std::vector<double> set_separation_transition(double x0, double y0, double x1, double y1, double ksep, double lmin, double angle);

/* Calculates the position of the center of a polygon
    / data structures: vertices: vector of vectors with the positions of the polygon's vertices
    / returns: a 3-element vector with the x y z coordinates for the center of the polygon */
std::vector<double> get_center(std::vector<std::vector<double> > vertices);

// Generates random angle theta 
double random_angle();

#endif