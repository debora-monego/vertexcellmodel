#include "vector.hpp"
#include "geometry.hpp"
#include "polygon.hpp"
#include "energy.hpp"

using namespace std;

// Get total energy
double get_total_energy(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, std::vector<std::vector<int> > edges, double ka, std::vector<double> L,
                        double Lambda, double gamma)
{
    // Elasticity energy: e1 = sum_{J=1 to N} (ka/2) * (AJ - A0)^2
    // ka: effective modulus of cell (resistance to volume changes)
    // AJ: area of cell J; A0: preferred area of cell J 
    double e1 = get_energy_elasticity(vertices, network, ka, L);
    // Adhesion energy: e2 = sum{ij} Lambda * lij
    // Lambda: adhesion energy/unit length
    // lij: length of edge between vertex i and j
    double e2 = get_energy_adhesion(vertices, edges, Lambda, L);
    // Take into account double counting of edges
    e2 = e2 / 4;
    // Contraction energy: e3 = sum_{J=1 to N} (gamma/2) * LJ^2
    // gamma/2: elastic constant
    // LJ: perimeter of cell J
    double e3 = get_energy_contraction(vertices, network, gamma, L);

    return (e1 + e2 + e3);
}

// Calculate elasticity energy
double get_energy_elasticity(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double ka, std::vector<double> L)
{
    double e = 0;
    double area = 0;
    double A0 = 0;

    for (int i = 0; i < network.size(); i++)
    {
        Polygon cell = network[i];
        area = cell.get_polygon_area(vertices, L);
        A0 = cell.A0;
        e += (ka / 2) * pow((area - A0),2);
    }
    return e;
}

// Calculate the adhesion energy
double get_energy_adhesion(std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, double Lambda, std::vector<double> L)
{
    double e = 0;
    double dist = 0;

    for (int i = 0; i < edges.size(); i++)
    {
        int i1 = edges[i][0];
        int i2 = edges[i][1];
        std::vector<double> v1 = vertices[i1];
        std::vector<double> vertex2 = vertices[i2];
        std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L));
        dist = get_euclidian_distance(v1[0], v1[1], v2[0], v2[1]);
        e += Lambda * dist;
    }
    return e;
}

// Calculate the contraction energy
double get_energy_contraction(std::vector<std::vector<double> > vertices, std::vector<Polygon> network, double gamma, std::vector<double> L)
{
    double e = 0;
    double perimeter = 0;

    for (int i = 0; i < network.size(); i++)
    {
        Polygon cell = network[i];
        perimeter = cell.get_polygon_perimeter(vertices, L);
        e += ((gamma / 2) * pow(perimeter, 2));
    }
    return e;
}