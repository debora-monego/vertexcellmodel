#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <chrono>

#include "MD.hpp"
#include "parser.hpp"
//#include "geometry.hpp"
//#include "vector.hpp"
#include "polygon.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
    srand(time(NULL));

    // command line arguments for data files
    string vertex_file = argv[1];
    string edges_file = argv[2];
    string cells_file = argv[3];

    // Parameters
    double lx = 10;
    double ly = 10;

    std::vector<double> L{lx, ly};

    double ka = 1;                 // area force coefficint
    double A0 = 1;               // current preferred area for polygon
    double gamma = 0.04 * ka * A0; // hexagonal network
    double Lambda = 0.12 * ka * sqrt(pow(A0, 3)); // hexagonal network
    double lmin = 0.07; // edge rearangment treshold
    double ksep = 1.5;
    double delta_t = 0.01; // timestep
    double xi = 0.2;       // motility coefficient
    double eta = 0.1;      // noise scalling coefficient

    double T = 0.02; // total simulation time
    bool T1_enabled = true; // enable T1 transitions

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

    auto start = high_resolution_clock::now();

    molecular_dynamics(vertices, edges, network, delta_t, L, T, ka, Lambda, gamma, eta, xi, lmin, ksep, T1_enabled);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken to run: "
         << duration.count() << " microseconds" << endl;

    return 0;
}