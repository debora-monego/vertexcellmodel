#include <iostream>
#include <vector>
#include <string>
//#include <armadillo>

#include "geometry.hpp"
//#include "polygon.hpp"
#include "vector.hpp"
//#include "vector.hpp"
#include "parser.hpp"

//using namespace std;
//using namespace arma;

/*./cell vertices.txt allala

void main(dkjad;la{}){
    vertcies = file argv[1]
} */

int main()
{

    /*double i, j, k;

    std::vector<std::vector<double> > vertices{{1, 1, 0}, {1, 2, 0}, {2, 2, 0}, {2, 1, 0}};
    //arma::vec L{4,4,0};
    std::vector<double> L{4, 4, 0};

    Polygon poly;

    // Test of PBC_DIFF function from GEOMETRY.HPP
    //std::vector<double> periodic_diff;
    //periodic_diff = cells.pbc_diff(vertices[0], vertices[1], L);
    //for (double i : periodic_diff)
    //    std::cout << i << ' ' << '\n';

    

    // Test of GET_POLYOGON_AREA function from POLYGON.HPP
    double area;
    area = poly.get_polygon_area(vertices, L);
    std::cout << "area = " << area << '\n';

    // Test of GET_POLYGON_PERIMETER function from POLYGON.HPP
    double perimeter;
    perimeter = poly.get_polygon_perimeter(vertices, L);
    std::cout << "perimeter = " << perimeter << '\n';

    // Test of GET_POLYGON_CENTER function from POLYGON.HPP
    std::vector<double> center;
    center = poly.get_polygon_center(vertices, L);
    for (double i : center)
    {
        std::cout << i << ' ' << '\n';
    } */

    //Test of READ_POLYGON_INDICES function from PARSER.HPP
    string networkfile = "cell_indices.txt";
    std::vector<vector<int> > network_indices;
    network_indices = read_polygon_indices(networkfile);

    for (std::vector<int> i : network_indices)
    {
        for (int k : i)
        {
            std::cout << k << ' ';
        }
        std::cout << '\n';
    }

    // Test of BUILD_NETWORK function from PARSER.HPP
    std::vector<Polygon> network;
    network = build_network(network_indices, 3.0);

    /*std::vector<Polygon>::iterator it;
    for (it = network.begin(); it != network.end(); ++it)
    {
        for (int i : it->indices)
        {
            //for (int k : i){
            std::cout << i << '\n';
            // }
        }
        //std::cout << '\n';
    } 

    // Test of GET_POLYGON_VERTICES function from POLYGON.HPP

    std::vector<std::vector<double> > polygon_vertices;
    polygon_vertices = network[i].get_polygon_vertices(vertices, L);
    for (std::vector<double> i : polygon_vertices)
    {
        for (double k : i)
        {
            std::cout << k << ' ';
        }
        std::cout << '\n';
    } */

    // Test of READ_VERTICES function from PARSER.HPP
    string vertexfile = "vertices.txt";
    std::vector<vector<double> > vertices;
    vertices = read_vertices(vertexfile);

    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < vertices[i].size(); j++)
        {
            std::cout << vertices[i][j] << ' ';
        }
        std::cout << '\n';
    }

    /*// Test of READ_EDGES function from PARSER.HPP
    string edgefile = "edges.txt";
    std::vector<vector<int> > edges;
    edges = read_edges(edgefile);

    for (int i = 0; i < edges.size(); i++)
    {
        for (int j = 0; j < edges[i].size(); j++)
        {
            std::cout << edges[i][j] << ' ';
        }
        std::cout << '\n';
    }*/

    // Test of SCALE_VECTOR function from VECTOR.HPP
    std::vector<double> scaled_vector;
    scaled_vector = scale_vector(vertices[0], 4);

    for (double i : scaled_vector)
        std::cout << i << ' ' << '\n';

    string outputfile = "test.txt";
    write_vertices(vertices, outputfile);

    std::vector<double> summation_vector;
    summation_vector = add_vectors(vertices[0], vertices[1]);
    std::cout << vertices[1][1]; //<< ' ' << vertices[1] << '\n';
    for (double i : summation_vector)
        std::cout << i << ' ' << '\n';

    
    
/*
    std::vector<double> center_vertices;

    center_vertices = cells.center(vertices);

    for (double i: center_vertices)
        std::cout << i << ' ' << '\n';




    double magnitude;

    magnitude = cells.magnitude(vertices[0]);
    std::cout << magnitude << '\n';

    double angle;
    angle = cells.random_angle();
    std::cout << angle << '\n';

    double vect;
    vect = cells.angle_2_vector(angle);
    for (double i: vect)
        std::cout << i << ' ' << '\n'; 
    
    double x = cos(angle);
    double y = sin(angle);

    double new_angle = cells.vector_2_angle(y,x);
    std::cout << new_angle << '\n' << '\n';

    std::vector<vector<double> > uv;
    uv = cells.angle_2_vector(angle);
    for(i = 0; i < 3; i++ ) {
        for(j = 0; j < 3; j++) {
            std::cout << uv[i][j] << ' ' << '\n';
        }
        std::cout << '\n';
    }
    
    std::vector<double> unit_vector;
    unit_vector = cells.unit_vector(vertices[0], vertices[1]);
    for (double i: unit_vector)
        std::cout << i << ' ' << '\n';

    double theta = cells.angle_diff(angle, new_angle);
    std::cout << theta << '\n';

    double dot_product = cells.get_dot_product(vertices[0], vertices[1]);
    std::cout << "dot product:" << dot_product << '\n';
  
    std::vector<std::vector<std::vector<double> > > zip_test;
    zip_test = cells.zip_matrix(vertices, vertices);
    std::cout << "zip restult:\n";
    for (std::vector<std::vector<double> > i: zip_test){
        for (std::vector<double> j: i) {
            for (double k: j) {
            std::cout << k << ' ';
            }
        std::cout << ' ';
        }
    std::cout << '\n';
    }

    Test of READ_POLYGON_INDICES function from PARSER.HPP
    string indicesfile = "cell_indices.txt";
    std::vector<vector<int> > network_indices;
    network_indices = read_polygon_indices(indicesfile);

    for (std::vector<int> i: network_indices) {
        for (int k: i){
            std::cout << k << ' ';
        }
        std::cout << '\n';
    } */
}