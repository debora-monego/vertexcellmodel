#include <math.h>
#include <iostream>





// Calculates the cross product between two vectors v1 and v2
std::vector<double> get_cross_product(std::vector<double> v1, std::vector<double> v2)
{
    std::vector<double> cross_product{(v1[1] * v2[2] - v1[2] * v2[1]),
                                      (v1[2] * v2[0] - v1[0] * v2[2]),
                                      (v1[0] * v2[1] - v1[1] * v2[0])};

    return cross_product;
}




// Calculates the magnitude of a vector
double magnitude(std::vector<double> v)
{
    double magnitude = sqrt((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));

    return magnitude;
}

// Finds the angle between two vectors v1 and v2
double get_angle_vectors(std::vector<double> v1, std::vector<double> v2)
{
    double theta = get_dot_product(v1, v2);
    theta = theta / magnitude(v1) * magnitude(v2);

    return acos(theta);
}