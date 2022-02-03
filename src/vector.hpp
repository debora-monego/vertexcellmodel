#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cassert>

using namespace std;

/* Returns a vector vzip that is the iteration of two other vectors v1 and v2, both with size n=n1=n2,
    / where the first element of vzip is the iteration of the first elements of v1 and v2 together,
    / the second element of vzip is the iteration of the second elements of v1 anv v2 together, and so on.
    / returns: a vector with n vectors */
std::vector<std::vector<double> > zip_vector(std::vector<double> v1, std::vector<double> v2);

std::vector<std::vector<std::vector<double> > > zip_matrix(std::vector<std::vector<double> > vertices1, std::vector<std::vector<double> > vertices2);

// Scales vector by a constant
std::vector<double> scale_vector(std::vector<double> vector, double s);

// Summation of vectors: creates a vector containing the sums of the elements of other vectors
std::vector<double> add_vectors(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> add_vectors(const std::vector<double> &v1, const std::vector<double> &v2, const std::vector<double> &v3, const std::vector<double> &v4);

// Subtraction of vectors: creates a vector containing the subtractions of the elements of two other vectors
std::vector<double> subtract_vectors(const std::vector<double>& v1, const std::vector<double>& v2);

// Calculates the dot product between two vectors v1 and v2
double get_dot_product(std::vector<double> v1, std::vector<double> v2);

// Calculates the dot product between a matrix m1 (2x2) and a vector v1 (2D)
std::vector<double> get_dot_product_matrix(std::vector<std::vector<double> > m1, std::vector<double> v1);

// Calculates the unit vector between two vectors
std::vector<double> get_unit_vector(std::vector<double> v1, std::vector<double> v2);

// Converts angle theta to unit vector 
std::vector<double> angle_2_vector(double theta);

// Converts vector to angle
double vector_2_angle(double x, double y);

// Generates random angle theta 
double random_angle(double min, double max);

// Calculates the sqrt of the square of each element in the vector
std::vector<double> absolute_vector(std::vector<double> vector);

#endif