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

double get_area(std::vector<std::vector<double>> vertices)
{
    int n_vertices = vertices.size();
    std::vector<std::vector<double>> vertices_rotate;
    for (int i = 0; i < n_vertices; i++)
    {
        vertices_rotate.push_back(vertices[i]);
    }
    std::rotate(vertices_rotate.begin(), vertices_rotate.begin() + n_vertices, vertices_rotate.end());

    std::vector<std::vector<std::vector<double>>> edges_coord;
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
double get_perimeter(std::vector<std::vector<double>> vertices)
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
std::vector<double> set_separation_transition(double x0, double y0, double x1, double y1, double ksep, double lmin, double angle)
{
    double original_separation = get_euclidian_distance(x0, y0, x1, y1);
    double new_separation = lmin * ksep;
    double m = new_separation / original_separation;

    std::vector<double> midpoint(2, 0);
    midpoint[0] = (x0 + x1) / 2;
    midpoint[1] = (y0 + y1) / 2;
    // cout << "x0 = " << x0 << " and x1 = " << x1 << '\n';
    // cout << "midpoint_x = " << midpoint[0];

    // Make the midpoint the origin
    std::vector<double> v0_mid(2, 0);
    std::vector<double> v1_mid(2, 0);
    v0_mid[0] = x0 - midpoint[0];
    v0_mid[1] = y0 - midpoint[1];
    v1_mid[0] = x1 - midpoint[0];
    v1_mid[1] = y1 - midpoint[1];
    //cout << "v0_mid = " << v0_mid[0] << " v1_mid = " << v1_mid[0] << '\n';

    std::vector<double> v0_rotated(2, 0);
    std::vector<double> v1_rotated(2, 0);
    v0_rotated[0] = cos(angle) * v0_mid[0] - sin(angle) * v0_mid[1];
    v0_rotated[1] = sin(angle) * v0_mid[0] + cos(angle) * v0_mid[1];
    v1_rotated[0] = cos(angle) * v1_mid[0] - sin(angle) * v1_mid[1];
    v1_rotated[1] = sin(angle) * v1_mid[0] + cos(angle) * v1_mid[1];
    // v0_rotated[0] = -v0_mid[1];
    // v0_rotated[1] = v0_mid[0];
    // v1_rotated[0] = -v1_mid[1];
    // v1_rotated[1] = v1_mid[0];
   // cout << "v0 rotated = " << v0_rotated[0] << " " << v0_rotated[1] << '\n';

    // Then add the midpoint coordinates to return to previous origin
    double x0_rotated = m*v0_rotated[0] + midpoint[0];
    double y0_rotated = m*v0_rotated[1] + midpoint[1];
    double x1_rotated = m*v1_rotated[0] + midpoint[0];
    double y1_rotated = m*v1_rotated[1] + midpoint[1];
    // double x0_rotated = v0_rotated[0] + midpoint[0];
    // double y0_rotated = v0_rotated[1] + midpoint[1];
    // double x1_rotated = v1_rotated[0] + midpoint[0];
    // double y1_rotated = v1_rotated[1] + midpoint[1];

    std::vector<double> set_separation{x0_rotated, y0_rotated, x1_rotated, y1_rotated};
    // std::vector<double> set_separation{new_x0, new_y0, new_x1, new_y1};
    // cout << x0 << "," << y0 << ","
    //      << "\t" << x1 << "," << y1 << '\n';
    // cout << x0_rotated << "," << y0_rotated << ","
    //      << "\t" << x1_rotated << "," << y1_rotated << '\n';
    return set_separation;
}

// //find the center
// cx = (x1+x2)/2;
// cy = (y1+y2)/2;

// //move the line to center on the origin
// x1-=cx; y1-=cy;
// x2-=cx; y2-=cy;

// //rotate both points
// xtemp = x1; ytemp = y1;
// x1=-ytemp; y1=xtemp;

// xtemp = x2; ytemp = y2;
// x2=-ytemp; y2=xtemp;

// //move the center point back to where it was
// x1+=cx; y1+=cy;
// x2+=cx; y2+=cy;
// {
//     double x_diff = x1 - x0;
//     double y_diff = y1 - y0;

//     double original_distance = get_euclidian_distance(x0, y0, x1, y1);
//     double transition_distance = ksep * lmin;

//     double new_x = x0;
//     double new_y = y0;

//     cout << "x_diff = " << x_diff << " y_diff = " << y_diff << '\n';

//     if (x_diff == 0)
//     {
//         new_x = x0;
//         if (y_diff > 0)
//         {
//             new_y = y0 + transition_distance;
//         }
//         else
//         {
//             new_y = y0 - transition_distance;
//         }
//     }
//     else if (y_diff == 0)
//     {
//         new_y = y0;
//         if (x_diff > 0)
//         {
//             new_x = x1 + transition_distance;
//         }
//         else
//         {
//             new_x = x1 - transition_distance;
//         }
//     }
//     else
//     {
//         double tan_angle = y_diff / x_diff;
//         double angle = tan(tan_angle);
//         new_x = round(x0 + transition_distance * (cos(angle)));
//         new_y = round(y0 + transition_distance * (sin(angle)));
//     }
//     std::vector<double> set_separation{new_x, new_y};
//     return set_separation;
// }

// Generates random angle theta
double random_angle()
{
    srand(time(NULL));
    // general function for random number: double randNum = rand()%(max-min + 1) + min;
    // Angle theta between -pi and +pi
    double theta = fmod(rand(), ((2 * pi + 1) - pi));
    return theta;
}

// Calculates the position of the center of a polygon
std::vector<double> get_center(std::vector<std::vector<double>> vertices)
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