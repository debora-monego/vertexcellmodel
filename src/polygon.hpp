#ifndef POLYGON_HPP
#define POLYGON_HPP

#include <math.h>
#include <string.h>
#include <vector>

using namespace std;

/*

Description of class that defines every polygon in the network

author: Debora Monego | Columbia University | 2021

Data structures:
vertices: list of all vertices in the network - <vector of vectores> [(x0,y0), (x1,y1), ..., (xNvertices,yNvertices)]
indices: list of indices mapping to vertices for every vertex in current polygon <list of integers> * counter-clockwise order
n_sides: number of sides in polygon for given cell
L: length of box * used to compute periodic boundary conditions

*/

class Polygon
{
public:
    int id;
    std::vector<int> indices;
    double A0;
    double P0;
    double J;
    double theta;

    void initialize(int read_id, std::vector<int> read_indices, double read_A0, double read_P0, double read_J, double read_theta)
    {
        id = read_id;
        indices = read_indices;
        A0 = read_A0;
        P0 = read_P0;
        J = read_J;
        theta = read_theta;
    }

    // Vertices with pbc
    // Returns: Matrix with the coordinates of vertices with periodic boundary conditions
    std::vector<vector<double> > get_polygon_vertices(std::vector<vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);

    // Polygon area
    // Returns: Double with polygon area
    double get_polygon_area(std::vector<vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);

    // Polygon perimeter
    // Returns: Double with polygon perimeter
    double get_polygon_perimeter(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);

    // Polygon center
    // Returns: Vector with coordinates of polygon center
    std::vector<double> get_polygon_center(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges);
};

#endif