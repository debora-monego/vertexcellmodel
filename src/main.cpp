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
    double lx, ly; // simulation box dimensions
    sscanf(argv[4],"%lf",&lx);
    sscanf(argv[5],"%lf",&ly);
    std::vector<double> L{lx, ly};
    
    double ka, A0, P0, J;
    sscanf(argv[6],"%lf",&ka);
    sscanf(argv[7],"%lf",&A0);
    sscanf(argv[8],"%lf",&P0);
    sscanf(argv[9],"%lf",&J);
    
    double gamma_f, Lambda_f;
    sscanf(argv[17],"%lf",&gamma_f);
    sscanf(argv[18],"%lf",&Lambda_f);
    double gamma = gamma_f * ka * A0; // hexagonal network
    double Lambda = Lambda_f * ka * sqrt(pow(A0, 3)); // hexagonal network

    double lmin, ksep;
    sscanf(argv[10],"%lf",&lmin);
    sscanf(argv[11],"%lf",&ksep);

    double xi, eta;
    sscanf(argv[12],"%lf",&xi);
    sscanf(argv[13],"%lf",&eta);

    double delta_t, T;
    sscanf(argv[14],"%lf",&delta_t);
    sscanf(argv[15],"%lf",&T);

    bool T1_enabled = argv[16];

    string out_folder = argv[19];

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
    network = build_network(vertex_indices, A0, P0, J);

    auto start = high_resolution_clock::now();

    molecular_dynamics(vertices, edges, network, delta_t, L, T, ka, Lambda, gamma, eta, xi, J, lmin, ksep, T1_enabled, out_folder);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken to run: "
         << duration.count() << " microseconds" << endl;

    return 0;
}