#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>

//#include "vector.hpp"

using namespace std;

/*
	This file defines geometrical variables and operations
	Author: Debora Monego | Columbia University | 2021
    data structures:
    vertices: list of vertices {v0, v1, ... vN} with positions (x0, y0, z0), (x1, y1, z1) ... (xN, yN, zN) 
*/



/* To do: Calculate center of polyhedra */




/* Calculates the volume of a polyhedra
    / returns: double /
    double volume(std::vector<std::vector<double> > vertices) {
        double base_area = area(vertices);
        double height_perp;

        double volume = base_area * height_perp;
    } */




/* Calculates the difference between two angles theta1 and theta2 
double angle_diff(double theta1, double theta2)
{
    double theta = theta1 - theta2;
    return (theta - 2 * pi * floor((theta + pi) / (2 * pi)));
}

/* Finds the angle between the faces forming vertex v1 at point p1
    / returns
    // To do: Needs to be generalized to handle 3D as well 
double get_angle_points(std::vector<double> p1, std::vector<double> p2, std::vector<double> p3)
{
    double angle_radian = 0;
    double p12 = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) + (p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));
    double p13 = sqrt((p1[0] - p3[0]) * (p1[0] - p3[0]) + (p1[1] - p3[1]) * (p1[1] - p3[1]) + (p1[2] - p3[2]) * (p1[2] - p3[2]));
    double p23 = sqrt((p2[0] - p3[0]) * (p2[0] - p3[0]) + (p2[1] - p3[1]) * (p2[1] - p3[1]) + (p2[2] - p3[2]) * (p2[2] - p3[2]));

    if (p12 != 0 and p13 != 0)
    {
        angle_radian = acos((p12 * p12 + p13 * p13 - p23 * p23) / (2 * p12 * p13));
    }
    else
    {
        std::cout << "The lenth of one (or both) of the faces is equal to zero.\n";
    }

    return angle_radian;
}

/* Converts radian to degrees
    / returns: double 
double readian_2_degrees(double theta)
{
    return theta * (360) / (2 * pi);
}

/* Checks whether the order of the vertices is counter-clockwise
    returns: boolean 
bool is_counter_clockwise(std::vector<std::vector<double> > cell_polygon)
{
    double sumEdges = 0;
    int i;
    int n = cell_polygon.size();
    for (i = 0; i < n; i++)
    {
        double x0 = cell_polygon[0][0];
        double y0 = cell_polygon[0][1];
        double x2 = 0;
        double y2 = 0;
        if ((i + 1) != n)
        {
            x2 = cell_polygon[i + 1][0];
            y2 = cell_polygon[i + 1][1];
        }
        else if ((i + 1) == n)
        {
            x2 = x0;
            y2 = y0;
        }
        sumEdges += float(x2 - cell_polygon[i][0]) / float(y2 + cell_polygon[i][1]);
    }
    if (sumEdges > 0)
    {
        return true;
    }
    else if (sumEdges < 0)
    {
        return false;
    }
    else
    {
        return 0;
    }
} */
