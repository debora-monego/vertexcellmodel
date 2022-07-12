#include "polygon.hpp"
#include "geometry.hpp"
#include "vector.hpp"

#include <vector>
#include <iostream>

using namespace std;

// Return list of vertices with periodic boundary conditions
std::vector<std::vector<double> > Polygon::get_polygon_vertices(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges)
{
    // array of x,y vertices in counter-clockwise order
    // moving vertices to maintain periodic boundaries
    std::vector<std::vector<double>> polygon_vertices;

    // align everything to previous vertex
    double x0 = vertices[indices[0]][0];
    double y0 = vertices[indices[0]][1];
    // double z0 = vertices[0][2];
    std::vector<int> q{0,0};

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

        if (i > 0)
        {
            for (int j = 0; j < edges.size(); j++)
            {
                if (edges[j][0] == indices[i - 1] && edges[j][1] == indices[i])
                {
                    std::vector<int>::const_iterator first = edges[j].begin() + 2;
                    std::vector<int>::const_iterator last = edges[j].begin() + 4;
                    std::vector<int> pbc(first, last);
                    q[0] = q[0] + pbc[0];
                    q[1] = q[1] + pbc[1];
                    // q.assign(first, last);
                    //  cout << i << '\t' << edges[j][0] << '\t' << q[0] << '\n';
                }
            }
        }

        //cout << q[0] << '\n';
        v_next = add_vectors(v_last, pbc_diff(v, v_last, L, q));
        polygon_vertices.insert(polygon_vertices.end(), v_next);
        v_last = v_next;
    }

    // for (std::vector<double> i : polygon_vertices)
    // {
    //     for (double j : i)
    //     {
    //         cout << j << '\t';
    //     }
    //     cout << '\n';
    // }
    return polygon_vertices;
}

// Return polygon area
double Polygon::get_polygon_area(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<std::vector<int> > edges)
{
    std::vector<std::vector<double>> polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L, edges);
    double area = get_area(polygon_vertices);
    return area;
}

// Return polygon perimeter
double Polygon::get_polygon_perimeter(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<vector<int> > edges)
{
    std::vector<std::vector<double>> polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L, edges);
    double perimeter = get_perimeter(polygon_vertices);
    return perimeter;
}

// Return polygon center
std::vector<double> Polygon::get_polygon_center(std::vector<std::vector<double> > vertices, std::vector<double> L, std::vector<vector<int> > edges)
{
    std::vector<std::vector<double>> polygon_vertices;
    polygon_vertices = get_polygon_vertices(vertices, L, edges);
    std::vector<double> polygon_center = get_center(polygon_vertices);
    return polygon_center;
}
