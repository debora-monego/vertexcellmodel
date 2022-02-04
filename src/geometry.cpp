#include "geometry.hpp"
#include "vector.hpp"

using namespace std;

double pi = atan(1) * 4;

std::vector<double> pbc_diff(std::vector<double> v1, std::vector<double> v2, std::vector<double> L)
{
    std::vector<double> periodic_diff;
    /*if (L[2] != 0)
    {
        double pbcx = fmod((v1[0] - v2[0] + L[0] / 2.), L[0]) - L[0] / 2;
        double pbcy = fmod((v1[1] - v2[1] + L[1] / 2.), L[1]) - L[1] / 2;
        double pbcz = fmod((v1[2] - v2[2] + L[2] / 2.), L[2]) - L[2] / 2;
        periodic_diff.insert(periodic_diff.end(), {pbcx, pbcy, pbcz});
    }
    else
    {*/
    double pbcx = fmod((v1[0] - v2[0] + (L[0] / 2.)), L[0]) - (L[0] / 2);
    double pbcy = fmod((v1[1] - v2[1] + (L[1] / 2.)), L[1]) - (L[1] / 2);
    periodic_diff.insert(periodic_diff.end(), {pbcx, pbcy});
    // }

    return periodic_diff;
}

double get_area(std::vector<std::vector<double> > vertices)
{
    int n_vertices = vertices.size();
    std::vector<std::vector<double> > vertices_rotate;
    for (int i = 0; i < n_vertices; i++)
    {
        vertices_rotate.push_back(vertices[i]);
    }
    std::rotate(vertices_rotate.begin(), vertices_rotate.begin() + n_vertices, vertices_rotate.end());

    std::vector<std::vector<std::vector<double> > > edges_coord;
    edges_coord = zip_matrix(vertices, vertices_rotate);

    double cross_product = 0;
    double x0, y0, x1, y1;

    for (int i = 0; i < n_vertices; i++)
    {
        x0 = edges_coord[i][0][0];
        y0 = edges_coord[i][0][1];
        // z0 = edges[i][0][2];
        if (i == (n_vertices - 1))
        {
            x1 = edges_coord[0][1][0];
            y1 = edges_coord[0][1][1];
            // z1 = edges[0][1][2];
        }
        else
        {
            x1 = edges_coord[i + 1][1][0];
            y1 = edges_coord[i + 1][1][1];
            // z1 = edges[i + 1][1][2];
        }
        cross_product += ((x0 * y1) - (x1 * y0));
    }

    return 0.5 * abs(cross_product);
}

// Calculates the Euclidian distance between two points with coordinates (x0,y0) and (x1,y1)
double get_euclidian_distance(double x0, double y0, double x1, double y1)
{
    double distance = sqrt(pow((x0 - x1), 2) + pow((y0 - y1), 2));
    return distance;
}

// Calculates the perimeter of a polygon
double get_perimeter(std::vector<std::vector<double> > vertices)
{

    int i;
    double x0, y0, x1, y1;
    int n = vertices.size();
    double perimeter = 0;

    for (i = 0; i < n; i++)
    {
        x0 = vertices[i][0];
        y0 = vertices[i][1];
        // z0 = vertices[i][2];

        if (i == (n - 1))
        {
            x1 = vertices[0][0];
            y1 = vertices[0][1];
            // z1 = vertices[0][2];
        }
        else
        {
            x1 = vertices[i + 1][0];
            y1 = vertices[i + 1][1];
            // z1 = vertices[i + 1][2];
        }
        double dist = get_euclidian_distance(x0, y0, x1, y1);
        perimeter += dist;
    }

    return perimeter;
}

// This function sets the new separation between two vertices after a T1 transition
// The new separation is equal to ksep * lmin
std::vector<double> set_separation_transition(double x0, double y0, double x1, double y1, double ksep, double lmin)
{
    double x_diff = x1 - x0;
    double y_diff = y1 - y0;

    double original_distance = get_euclidian_distance(x0, x1, y0, y1);
    double transition_distance = ksep * lmin;

    double new_x = x0;
    double new_y = y0;

    std::vector<double> set_separation{new_x, new_y};

    if (x_diff == 0)
    {
        new_x = x0;
        if (y_diff > 0)
        {
            new_y = y0 + original_distance;
        }
        else
        {
            new_y = y0 - original_distance;
        }
    }
    else if (y_diff == 0)
    {
        new_y = y0;
        if (x_diff > 0)
        {
            new_x = x1 + original_distance;
        }
        else
        {
            new_x = x1 - original_distance;
        }
    }
    else
    {
        double tan_angle = y_diff / x_diff;
        double angle = tan(tan_angle);
        new_x = round(x0 + transition_distance * cos(angle));
        new_y = round(y0 + transition_distance * sin(angle));
    }
    return set_separation;
}

// Calculates the position of the center of a polygon
std::vector<double> get_center(std::vector<std::vector<double> > vertices)
{

    int i;
    int n_vertices = vertices.size();
    double sumX = 0;
    double sumY = 0;
    // double sumZ = 0;

    // sum vertices positions in each component
    for (i = 0; i < n_vertices; i++)
    {
        sumX += vertices[i][0];
        sumY += vertices[i][1];
        // sumZ += vertices[i][2];
    }

    std::vector<double> center_polygon;
    center_polygon.insert(center_polygon.end(), {(sumX / n_vertices), (sumY / n_vertices)}); //, (sumZ / n_vertices)});

    return center_polygon;
}

// Generates random angle theta
double random_angle()
{
    srand(time(NULL));
    // general function for random number: double randNum = rand()%(max-min + 1) + min;
    // Angle theta between -pi and +pi
    double theta = fmod(rand(), ((2 * pi + 1) - pi));
    return theta;
}