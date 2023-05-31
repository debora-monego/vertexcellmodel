#include "vector.hpp"
#include "geometry.hpp"
#include "polygon.hpp"
#include "force.hpp"

using namespace std;

// Computes forces in the current configuration of the vertex model
std::vector<std::vector<double>> get_forces(std::vector<std::vector<double>> vertices, std::vector<Polygon> network, std::vector<std::vector<int>> edges,
                                            std::vector<double> L, double ka, double Lambda, double gamma, double eta, double xi, double J, double lambda_potts)
{
    std::vector<std::vector<double>> f1 = calc_force_elasticity(vertices, network, ka, L, edges);

    // // Vertex model
    // std::vector<std::vector<double> > f2 = calc_force_adhesion(vertices, edges, Lambda, L);
    // std::vector<std::vector<double> > f3 = calc_force_contraction(vertices, network, gamma, L, edges);

    // Potts model
    std::vector<std::vector<double> > f2 = calc_force_perimeter_potts(vertices, network, lambda_potts, L, edges);
    std::vector<std::vector<double> > f3 = calc_force_adhesion_potts(vertices, edges, network, J, L);

    // // Add motility force
    // std::vector<std::vector<double>> f4 = calc_force_motility(vertices, network, eta, xi);


    // Initialize total forces vector
    std::vector<std::vector<double>> result_force(f2.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(result_force[i].begin(), result_force[i].end(), 0);
    }

    for (int i = 0; i < vertices.size(); i++)
    {
        assert(f1.size() == f2.size());
        assert(f3.size() == f2.size());
        //assert(f4.size() == f1.size());

        // result_force[i] = scale_vector(add_vectors(f1[i], f2[i], f3[i], f4[i]), -1);
        result_force[i] = add_vectors(f1[i], f2[i]);
        result_force[i] = scale_vector(add_vectors(result_force[i], f3[i]), -1);
    }
    return result_force;
}

std::pair<std::vector<std::vector<double> >, std::vector<std::vector<int> > > move_vertices(std::vector<std::vector<double> > vertices, std::vector<std::vector<double> > forces, std::vector<std::vector<int> > edges,
                                                                                         std::vector<double> L, double delta_t, double lmin)
{
    // Initialize new vertices matrix
    std::vector<std::vector<double>> updated_vertices(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(updated_vertices[i].begin(), updated_vertices[i].end(), 0);
    }

    for (int i = 0; i < vertices.size(); i++)
    {
        assert(vertices.size() == forces.size());
        updated_vertices[i] = add_vectors(vertices[i], scale_vector(forces[i], delta_t));

        double diff_dist = get_euclidian_distance(vertices[i][0], vertices[i][1], updated_vertices[i][0], updated_vertices[i][1]);

        // Update vertices positions only if the vertex displacement is smaller than lmin/2
        // This step should prevent cells from intersecting each other. ref: https://doi.org/10.1016/j.pbiomolbio.2013.09.003
        if (diff_dist < (lmin / 2))
        {
            vertices[i] = updated_vertices[i];

            // Wrap around periodic boundaries and update periodicity vector q in edges
            if (vertices[i][0] < 0) // vertex moved across left boundary
            {
                // Wrap around to right
                vertices[i][0] = vertices[i][0] + L[0];
                for (int j = 0; j < edges.size(); j++)
                {
                    if (i == edges[j][0])
                    {
                        edges[j][2]++;
                    }
                    if (i == edges[j][1])
                    {
                        edges[j][2]--;
                    }
                }
            }
            else if (vertices[i][0] > L[0]) // vertex moved across right boundary
            {
                // Wrap around to left
                vertices[i][0] = vertices[i][0] - L[0];
                for (int j = 0; j < edges.size(); j++)
                {
                    if (i == edges[j][0])
                    {
                        edges[j][2]--;
                    }
                    if (i == edges[j][1])
                    {
                        edges[j][2]++;
                    }
                }
            }

            if (vertices[i][1] < 0) // vertex moved across bottom boundary
            {
                // Wrap around to top
                vertices[i][1] = vertices[i][1] + L[1];
                for (int j = 0; j < edges.size(); j++)
                {
                    if (i == edges[j][0])
                    {
                        edges[j][3]++;
                    }
                    if (i == edges[j][1])
                    {
                        edges[j][3]--;
                    }
                }
            }
            else if (vertices[i][1] > L[1]) // vertex moved across top boundary
            {
                // Wrap around to bottom
                vertices[i][1] = vertices[i][1] - L[1];
                for (int j = 0; j < edges.size(); j++)
                {
                    if (i == edges[j][0])
                    {
                        edges[j][3]--;
                    }
                    if (i == edges[j][1])
                    {
                        edges[j][3]++;
                    }
                }
            }
        }
    }

    std::pair<std::vector<std::vector<double>>, std::vector<vector<int>>> updated_connections(vertices, edges);

    return updated_connections;
}

std::vector<double> get_clockwise(int index, std::vector<int> indices, std::vector<std::vector<double>> vertices, std::vector<double> L, std::vector<std::vector<int>> edges)
{
    // Find position of vertex in indices list
    int position = -1;
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        if (indices[i] == index)
        {
            position = i;
            break;
        }
    }

    // Clockwise - left -> right
    // Wrap around to 0 if at the end of the list
    if (position == (indices.size() - 1))
    {
        position = 0;
    }
    else
    {
        position += 1;
    }

    // Calculate vertex position with respect to boundaries
    std::vector<double> v0 = vertices[index];
    std::vector<double> v = vertices[indices[position]];
    std::vector<int> q{0, 0};
    for (int j = 0; j < edges.size(); j++)
    {
        if (edges[j][0] == index && edges[j][1] == indices[position])
        {
            std::vector<int>::const_iterator first = edges[j].begin() + 2;
            std::vector<int>::const_iterator last = edges[j].begin() + 4;
            std::vector<int> pbc(first, last);
            q[0] = q[0] + pbc[0];
            q[1] = q[1] + pbc[1];
        }
    }
    std::vector<double> vc = add_vectors(v0, pbc_diff(v, v0, L, q));

    return vc;
}

std::vector<double> get_counter_clockwise(int index, std::vector<int> indices, std::vector<std::vector<double>> vertices, std::vector<double> L, std::vector<std::vector<int>> edges)
{
    // Find position of vertex in indices list
    int position = -1;
    for (std::size_t i = 0; i < indices.size(); ++i)
    {
        if (indices[i] == index)
        {
            position = i;
            break;
        }
    }

    // Counter-clockwise - rigth -> left
    // Wrap around to 0 if at the end of the list
    if (position == (0))
    {
        position = (indices.size() - 1);
    }
    else
    {
        position -= 1;
    }

    // Calculate vertex position with respect to boundaries
    std::vector<double> v0 = vertices[index];
    std::vector<double> v = vertices[indices[position]];
    std::vector<int> q{0, 0};
    for (int j = 0; j < edges.size(); j++)
    {
        if (edges[j][0] == index && edges[j][1] == indices[position])
        {
            std::vector<int>::const_iterator first = edges[j].begin() + 2;
            std::vector<int>::const_iterator last = edges[j].begin() + 4;
            std::vector<int> pbc(first, last);
            q[0] = q[0] + pbc[0];
            q[1] = q[1] + pbc[1];
        }
    }
    std::vector<double> vcc = add_vectors(v0, pbc_diff(v, v0, L, q));

    return vcc;
}

// Calculate force on vertex due to elasticity -- same for Potts and Vertex models
std::vector<std::vector<double>> calc_force_elasticity(std::vector<vector<double>> vertices, std::vector<Polygon> network, double ka, std::vector<double> L, std::vector<std::vector<int>> edges)
{
    // Initialize force associated with vertex
    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    // Iterate over vertices and calculate force
    for (int i = 0; i < vertices.size(); i++)
    {
        int this_vertex = i;

        // Find polygons with this vertex
        for (int j = 0; j < network.size(); j++)
        {
            Polygon cell = network[j];
            bool vertex_in_polygon = (std::find(cell.indices.begin(), cell.indices.end(), this_vertex) != cell.indices.end());
            if (vertex_in_polygon != 0)
            {
                std::vector<int> indices = cell.indices;
                // Get clockwise vector
                std::vector<double> vc = get_clockwise(this_vertex, indices, vertices, L, edges);

                // Get counter-clockwise vector
                std::vector<double> vcc = get_counter_clockwise(this_vertex, indices, vertices, L, edges);

                // Get the difference vector
                std::vector<double> diff = subtract_vectors(vc, vcc);

                // Coompute perpendicular vector to assure correct direction for force is being chosen, i.e. pointing towards vertex (2d)
                std::vector<std::vector<double>> perp_matrix(2, std::vector<double>(2));
                for (int k = 0; k < perp_matrix.size(); k++)
                {
                    std::fill(perp_matrix[k].begin(), perp_matrix[k].end(), 0);
                }
                perp_matrix[0][1] = 1;
                perp_matrix[1][0] = -1;

                std::vector<double> force = scale_vector((get_dot_product_matrix(perp_matrix, diff)), -0.5);

                // force contribution from this polygon is stored in force
                double coeff = ka * (cell.A0 - cell.get_polygon_area(vertices, L, edges));
                forces[i] = add_vectors(forces[i], (scale_vector(force, coeff)));
            }
        }
    }
    return forces;
}

// Calculate force perimeter -- Potts model
std::vector<std::vector<double>> calc_force_perimeter_potts(std::vector<vector<double>> vertices, std::vector<Polygon> network, double lambda_potts, std::vector<double> L, std::vector<std::vector<int>> edges)
{
    // Initialize force associated with vertex
    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    // Iterate over vertices and calculate force
    for (int i = 0; i < vertices.size(); i++)
    {
        int this_vertex = i;

        // Find polygons with this vertex
        for (int j = 0; j < network.size(); j++)
        {
            Polygon cell = network[j];
            bool vertex_in_polygon = (std::find(cell.indices.begin(), cell.indices.end(), this_vertex) != cell.indices.end());
            if (vertex_in_polygon != 0)
            {
                std::vector<int> indices = cell.indices;
                // Get clockwise vector
                std::vector<double> vc = get_clockwise(this_vertex, indices, vertices, L, edges);

                // Get counter-clockwise vector
                std::vector<double> vcc = get_counter_clockwise(this_vertex, indices, vertices, L, edges);

                // Get the difference vector
                std::vector<double> diff = subtract_vectors(vc, vcc);

                // Coompute perpendicular vector to assure correct direction for force is being chosen, i.e. pointing towards vertex (2d)
                std::vector<std::vector<double>> perp_matrix(2, std::vector<double>(2));
                for (int k = 0; k < perp_matrix.size(); k++)
                {
                    std::fill(perp_matrix[k].begin(), perp_matrix[k].end(), 0);
                }
                perp_matrix[0][1] = 1;
                perp_matrix[1][0] = -1;

                std::vector<double> force = scale_vector((get_dot_product_matrix(perp_matrix, diff)), -0.5);

                // force contribution from this polygon is stored in force
                double coeff = lambda_potts * (cell.P0 - cell.get_polygon_perimeter(vertices, L, edges));
                forces[i] = add_vectors(forces[i], (scale_vector(force, coeff)));
            }
        }
    }
    return forces;
}

// Calculate force adhesion -- Potts model
std::vector<std::vector<double>> calc_force_adhesion_potts(std::vector<std::vector<double>> vertices, std::vector<std::vector<int>> edges, std::vector<Polygon> network, double J, std::vector<double> L)
{

    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    for (int i = 0; i < edges.size(); i++)
    {
        int i1 = edges[i][0];
        int i2 = edges[i][1];
        std::vector<int>::const_iterator first = edges[i].begin() + 2;
        std::vector<int>::const_iterator last = edges[i].begin() + 4;
        std::vector<int> q(first, last);
        std::vector<double> v1 = vertices[i1];
        std::vector<double> vertex2 = vertices[i2];
        std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));
        std::vector<double> uv = get_unit_vector(v1, v2);
        forces[i1] = add_vectors(forces[i1], scale_vector(uv, J));
    }
    return forces;
}

// Calculate force due to contraction - Vertex model
std::vector<std::vector<double>> calc_force_contraction(std::vector<std::vector<double>> vertices, std::vector<Polygon> network, double gamma, std::vector<double> L, std::vector<std::vector<int>> edges)
{
    // Initialize force associated with vertex
    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    // Iterate over vertices and calculate force
    for (int i = 0; i < vertices.size(); i++)
    {
        int this_vertex = i;

        // Find polygons with this vertex
        for (int j = 0; j < network.size(); j++)
        {
            Polygon cell = network[j];
            bool vertex_in_polygon = (std::find(cell.indices.begin(), cell.indices.end(), this_vertex) != cell.indices.end());
            if (vertex_in_polygon != 0)
            {
                // Get clockwise vector
                std::vector<double> vc = get_clockwise(this_vertex, cell.indices, vertices, L, edges);
                std::vector<double> uvc = get_unit_vector(vertices[i], vc);

                // Get counter-clockwise vector
                std::vector<double> vcc = get_counter_clockwise(this_vertex, cell.indices, vertices, L, edges);
                std::vector<double> uvcc = get_unit_vector(vcc, vertices[i]);

                // Get the difference vector
                std::vector<double> diff = subtract_vectors(uvc, uvcc);

                // Coompute perpendicular vector to assure correct direction for force is being chosen, i.e. pointing towards vertex (2d)
                std::vector<std::vector<double>> perp_matrix(2, std::vector<double>(2));
                for (int k = 0; k < perp_matrix.size(); k++)
                {
                    std::fill(perp_matrix[k].begin(), perp_matrix[k].end(), 0);
                }
                perp_matrix[0][1] = 1;
                perp_matrix[1][0] = -1;

                std::vector<double> force = scale_vector((get_dot_product_matrix(perp_matrix, diff)), -0.5);

                // Get perimenter for this polygon
                double perimeter = cell.get_polygon_perimeter(vertices, L, edges);

                // Force contribution from this polygon is stored in force
                forces[i] = add_vectors(forces[i], scale_vector(force, (- gamma * perimeter)));
            }
        }
    }
    return forces;
}

// Calculate force due to adhesion - interfacial energy - Vertex model
std::vector<std::vector<double>> calc_force_adhesion(std::vector<std::vector<double>> vertices, std::vector<std::vector<int>> edges, double Lambda, std::vector<double> L)
{
    // Initialize force associated with vertex
    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    for (int i = 0; i < edges.size(); i++)
    {
        int i1 = edges[i][0];
        int i2 = edges[i][1];
        std::vector<int>::const_iterator first = edges[i].begin() + 2;
        std::vector<int>::const_iterator last = edges[i].begin() + 4;
        std::vector<int> q(first, last);
        std::vector<double> v1 = vertices[i1];
        std::vector<double> vertex2 = vertices[i2];
        std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));
        std::vector<double> uv = get_unit_vector(v1, v2);
        forces[i1] = add_vectors(forces[i1], scale_vector(uv, Lambda));
    }
    return forces;
}

// Force to move vertices of polygons in a particular direction
std::vector<std::vector<double>> calc_force_motility(std::vector<std::vector<double>> vertices, std::vector<Polygon> network, double eta, double xi)
{
    // Initialize force associated with vertex
    std::vector<std::vector<double>> forces(vertices.size(), std::vector<double>(2));
    for (int i = 0; i < vertices.size(); i++)
    {
        std::fill(forces[i].begin(), forces[i].end(), 0);
    }

    // Find neighbors for every polygon, i.e. two polygons that share a vertex
    std::vector<std::vector<double>> avg_angles(network.size(), std::vector<double>(2));
    for (int i = 0; i < network.size(); i++)
    {
        std::fill(avg_angles[i].begin(), avg_angles[i].end(), 0);
    }

    std::vector<double> neighbor_count(network.size());
    std::fill(neighbor_count.begin(), neighbor_count.end(), 1);

    for (int i = 0; i < network.size(); i++)
    {
        Polygon cell1 = network[i];
        avg_angles[i] = add_vectors(avg_angles[i], angle_2_vector(cell1.theta));
        for (int j = 0; j < 2; j++)
        {
            for (int m = 0; m < network.size(); m++)
            {
                Polygon cell2 = network[m];
                if (i != m)
                {
                    std::vector<int> a = cell1.indices;
                    std::vector<int> b = cell2.indices;
                    bool is_neighbor = false;
                    for (int k = 0; k < a.size(); k++)
                    {
                        if (is_neighbor)
                            break;
                        for (int n = 0; n < b.size(); n++)
                        {
                            if (a[k] == b[n])
                            {
                                is_neighbor = true;
                                break;
                            }
                        }
                    }
                    avg_angles[i] = add_vectors(avg_angles[i], angle_2_vector(cell2.theta));
                    neighbor_count[i] += 1;
                }
            }
        }
    }

    for (int i = 0; i < network.size(); i++)
    {
        Polygon cell = network[i];
        double pi = atan(1) * 4;

        // noise variable
        double noisex = random_angle(-pi, pi);
        double noisey = random_angle(-pi, pi);
        // double nz = random_angle(-pi, pi);
        std::vector<double> noise_random;
        noise_random.insert(noise_random.end(), {noisex, noisey}); //, nz});

        // Average over all unit vectos for angles
        std::vector<double> avg;
        avg = scale_vector(avg_angles[i], (1 / neighbor_count[i]));

        // Add this force direction for every vertex in current polygon
        for (int ind = 0; ind < cell.indices.size(); ind++)
        {
            int index = cell.indices[ind];
            forces[index] = add_vectors(forces[index], (scale_vector(add_vectors(avg, scale_vector(noise_random, eta)), xi)));
        }
        // theta = avg + eta * noise
        cell.theta = vector_2_angle((eta * noise_random[0] + avg[0]), (eta * noise_random[1] + avg[1]));
    }

    return forces;
}