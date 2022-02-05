#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include "MD.hpp"
#include "parser.hpp"
//#include "geometry.hpp"
//#include "vector.hpp"
#include "polygon.hpp"

using namespace std;

int main(int argc, char *argv[])
{
    // command line arguments for data files
    string vertex_file = argv[1];
    string edges_file = argv[2];
    string cells_file = argv[3];

    // Open output files to write data

    // Parameters
    // double lx = 9 * (2 / sqrt(3 * (sqrt(3))));
    // double ly = 4 * (2 / sqrt(sqrt(3)));
    double lx = 1.0;
    double ly = 1.0;

    std::vector<double> L{lx, ly};

    double ka = 1; //area force coefficint
    double A0 = 0.01; // current preferred area for polygon
    //double gamma = 0.04 * ka * A0; // hexagonal network
    double gamma = 1;
    // gamma = 0.1 * ka * A0 // soft network
    //double Lambda = 0.12 * ka * sqrt(pow(A0, 3)); // hexagonal network
    double Lambda = 0.5;
    // Lambda = -0.85 * ka * A0**(3/2) // soft network
    double lmin = 0.1; // cell rearangment treshold
    double ksep = 3.0;
    double delta_t = 0.001; // timestep
    double xi = 0.2; // motility coefficient
    double eta = 0.01; // noise scalling coefficient

    // maximum Time
    double T = 0.01;
    double pi = atan(1) * 4;

    // open vertices file
    std::vector<std::vector<double> > vertices;
    vertices = read_vertices(vertex_file);

    // open edges file
    std::vector<std::vector<int> > edges;
    edges = read_edges(edges_file);

    // open polygons (cells) file
    std::vector<std::vector<int> > vertex_indices;
    vertex_indices = read_cell_indices(cells_file);
    std::vector<Polygon> network;
    network = build_network(vertex_indices, A0);

    molecular_dynamics(vertices, edges, network, delta_t, lx, ly, T, ka, Lambda, gamma, eta, xi, lmin, ksep);

    return 0;
}