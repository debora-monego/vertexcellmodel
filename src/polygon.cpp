#include "polygon.hpp"
#include "geometry.hpp"
#include "vector.hpp"

#include <iostream>

using namespace std;

// Return list of vertices with periodic boundary conditions
std::vector<std::vector<double> > Polygon::get_polygon_vertices(std::vector<std::vector<double> > vertices, std::vector<double> L)
{
    // array of x,y vertices in counter-clockwise order
    // moving vertices to maintain periodic boundaries
    std::vector<std::vector<double> > polygon_vertices;

    // align everything to previous vertex
    double x0 = vertices[indices[0]][0];
    double y0 = vertices[indices[0]][1];
    //double z0 = vertices[0][2];

    std::vector<double> v0;
    v0.insert(v0.end(), {x0, y0}); //, z0});
    std::vector<double> v_last(v0);

    for (int i = 0; i < indices.size(); i++)
    {
        double x = vertices[indices[i]][0];
        double y = vertices[indices[i]][1];
        // double z = vertices[i][2];

        std::vector<double> v;
        v.insert(v.end(), {x, y}); //, z});
        std::vector<double> v_next;
        v_next = pbc_diff(v, v_last, L);
        polygon_vertices.insert(polygon_vertices.end(), v_next);
        v_last = v_next;
    }
    return polygon_vertices;
}

// Return polygon area
double Polygon::get_polygon_area(std::vector<std::vector<double> > vertices, std::vector<double> L)
{
    std::vector<std::vector<double> > polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L);
    double area = get_area(polygon_vertices);
    return area;
}

// Return polygon perimeter
double Polygon::get_polygon_perimeter(std::vector<std::vector<double> > vertices, std::vector<double> L)
{
    std::vector<std::vector<double> > polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L);
    double perimeter = get_perimeter(polygon_vertices);
    return perimeter;
}

// Return polygon center
std::vector<double> Polygon::get_polygon_center(std::vector<std::vector<double> > vertices, std::vector<double> L)
{
    std::vector<std::vector<double> > polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L);
    std::vector<double> polygon_center = get_center(polygon_vertices);
    return polygon_center;
}