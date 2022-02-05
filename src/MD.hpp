#include <math.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple>

#include <iostream>
#include <fstream>
#include <cstdlib>
//#include <filesystem>

#include "vector.hpp"
//#include "geometry.hpp"
#include "energy.hpp"
#include "force.hpp"
#include "transitionT1.hpp"
#include "parser.hpp"
#include "polygon.hpp"

using namespace std;

void molecular_dynamics(std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, std::vector<Polygon> network,
                        double delta_t, double lx, double ly, double T, double ka, double Lambda, double gamma, double eta, double xi, double lmin, double ksep)
{

    // Vertices coordinates
    ofstream out_vertexfile;
    out_vertexfile.open("./results/outvertex.txt");

    ofstream out_edgefile;
    out_edgefile.open("./results/outedges.txt");

    ofstream out_cellfile;
    out_cellfile.open("./results/outcells.txt");

    // Energy, force and log data
    ofstream out_logfile;
    out_logfile.open("./results/log.txt");
    out_logfile << "TOTAL ENERGY AND FORCE IN THE NETWORK AT EVERY TIME STEP\n";
    out_logfile << "TOTAL TIME STEPS = " << T / delta_t << '\n';
    out_logfile << "timestep" << std::setw(20) << "f_total" << std::setw(20) << "f_elasticity" << std::setw(20) << "f_contraction" << std::setw(20)
                << "f_adhesion" << std::setw(20) << "f_motility" << std::setw(20) << "e_total" << std::setw(20) << "e_elasticity" << std::setw(20)
                << "e_adhesion" << std::setw(20) << "e_contraction" << '\n';
    // out_logfile << "timestep" << std::setw(20) << "f_total" << std::setw(20) << "f_elasticity" << std::setw(20) << "f_contraction" << std::setw(20)
    //             << "f_adhesion" << std::setw(20) << "e_total" << std::setw(20) << "e_elasticity" << std::setw(20)
    //             << "e_adhesion" << std::setw(20) << "e_contraction" << '\n';

    std::vector<double> L{lx, ly};

    for (double t = 0; t < T; t = t + delta_t)
    {
        //////// Get energies for the network ////////

        // Get elasticity energy
        double e_elasticity = get_energy_elasticity(vertices, network, ka, L);

        // Get adhesion energy
        double e_adhesion = get_energy_adhesion(vertices, edges, Lambda, L);

        // Get contraction energy
        double e_contraction = get_energy_contraction(vertices, network, gamma, L);

        // Get total energy
        double e_total = get_total_energy(vertices, network, edges, ka, L, Lambda, gamma);

        //////// Get forces for network ////////

        // Get elasticity force
        // Initialize forces for every vertex
        std::vector<std::vector<double> > f_elasticity(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(f_elasticity[i].begin(), f_elasticity[i].end(), 0);
        }

        f_elasticity = calc_force_elasticity(vertices, network, ka, L);

        std::vector<std::vector<double> > abs_f_elasticity(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(abs_f_elasticity[i].begin(), abs_f_elasticity[i].end(), 0);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            abs_f_elasticity[i] = absolute_vector(f_elasticity[i]);
        }

        std::vector<double> f_elasticity_sum(vertices.size());

        for (int i = 0; i < vertices.size(); i++)
        {
            f_elasticity_sum[i] = std::accumulate(abs_f_elasticity[i].begin(), abs_f_elasticity[i].end(), 0.0);
        }
        double f_elasticity_total = std::accumulate(f_elasticity_sum.begin(), f_elasticity_sum.end(), 0.0);

        // Get contraction force
        std::vector<std::vector<double> > f_contraction(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(f_contraction[i].begin(), f_contraction[i].end(), 0);
        }

        f_contraction = calc_force_contraction(vertices, network, gamma, L);

        std::vector<std::vector<double> > abs_f_contraction(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(abs_f_contraction[i].begin(), abs_f_contraction[i].end(), 0);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            abs_f_contraction[i] = absolute_vector(f_contraction[i]);
        }

        std::vector<double> f_contraction_sum(vertices.size());

        for (int i = 0; i < vertices.size(); i++)
        {
            f_contraction_sum[i] = std::accumulate(abs_f_contraction[i].begin(), abs_f_contraction[i].end(), 0.0);
        }
        double f_contraction_total = std::accumulate(f_contraction_sum.begin(), f_contraction_sum.end(), 0.0);

        // Get adhesion force
        std::vector<std::vector<double> > f_adhesion(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(f_adhesion[i].begin(), f_adhesion[i].end(), 0);
        }

        f_adhesion = calc_force_adhesion(vertices, edges, Lambda, L);

        std::vector<std::vector<double> > abs_f_adhesion(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(abs_f_adhesion[i].begin(), abs_f_adhesion[i].end(), 0);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            abs_f_adhesion[i] = absolute_vector(f_adhesion[i]);
        }

        std::vector<double> f_adhesion_sum(vertices.size());

        for (int i = 0; i < vertices.size(); i++)
        {
            f_adhesion_sum[i] = std::accumulate(abs_f_adhesion[i].begin(), abs_f_adhesion[i].end(), 0.0);
        }
        double f_adhesion_total = std::accumulate(f_adhesion_sum.begin(), f_adhesion_sum.end(), 0.0);

        // Get motility force
        std::vector<std::vector<double> > f_motility(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(f_motility[i].begin(), f_motility[i].end(), 0);
        }

        f_motility = calc_force_motility(vertices, network, eta, xi);

        std::vector<std::vector<double> > abs_f_motility(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(abs_f_motility[i].begin(), abs_f_motility[i].end(), 0);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            abs_f_motility[i] = absolute_vector(f_motility[i]);
        }

        std::vector<double> f_motility_sum(vertices.size());

        for (int i = 0; i < vertices.size(); i++)
        {
            f_motility_sum[i] = std::accumulate(abs_f_motility[i].begin(), abs_f_motility[i].end(), 0.0);
        }
        double f_motility_total = std::accumulate(f_motility_sum.begin(), f_motility_sum.end(), 0.0);

        // Get total force
        // Initialize forces for every vertex
        std::vector<std::vector<double> > forces(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(forces[i].begin(), forces[i].end(), 0);
        }

        forces = get_forces(vertices, network, edges, lx, ly, ka, Lambda, gamma, eta, xi);

        std::vector<std::vector<double> > abs_forces(vertices.size(), std::vector<double>(2));
        for (int i = 0; i < vertices.size(); i++)
        {
            std::fill(abs_forces[i].begin(), abs_forces[i].end(), 0);
        }

        for (int i = 0; i < vertices.size(); i++)
        {
            abs_forces[i] = absolute_vector(forces[i]);
        }

        std::vector<double> f_sum(vertices.size());

        for (int i = 0; i < vertices.size(); i++)
        {
            f_sum[i] = std::accumulate(abs_forces[i].begin(), abs_forces[i].end(), 0.0);
        }
        double f_total = std::accumulate(f_sum.begin(), f_sum.end(), 0.0);

        // Print results for every time step in output file
        out_logfile << t << std::setw(20) << f_total << std::setw(20) << f_elasticity_total << std::setw(20) << f_contraction_total << std::setw(20)
                    << f_adhesion_total << std::setw(20) << f_motility_total << std::setw(20) << e_total << std::setw(20) << e_elasticity << std::setw(20)
                    << e_adhesion << std::setw(20) << e_contraction << '\n';
        // out_logfile << t << std::setw(20) << f_total << std::setw(20) << f_elasticity_total << std::setw(20) << f_contraction_total << std::setw(20)
        //             << f_adhesion_total << std::setw(20) << e_total << std::setw(20) << e_elasticity << std::setw(20)
        //             << e_adhesion << std::setw(20) << e_contraction << '\n';

        // Move vertices
        vertices = move_vertices(vertices, forces, lx, ly, delta_t, lmin);
        out_vertexfile << "timestep = " << t << '\n';
        out_vertexfile << "vx" << std::setw(10) << "vy" << '\n';

        // Check for T1 transitions
        std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data = T1_transition(vertices, network, edges, lx, ly, lmin, ka, Lambda, gamma, ksep);
        network = std::get<0>(T1_data);
        edges = std::get<1>(T1_data);
        vertices = std::get<2>(T1_data);
        out_edgefile << "timestep = " << t << '\n';
        out_edgefile << "edge1" << std::setw(10) << "edge2" << '\n';
        for (std::vector<int> i : edges)
        {
            for (int j : i)
            {
                out_edgefile << j << std::setw(10);
            }
            out_edgefile << '\n';
        }
        for (std::vector<double> i : vertices)
        {
            for (double j : i)
            {
                out_vertexfile << j << std::setw(10);
            }
            out_vertexfile << '\n';
        }
    }

    // Print log information
    out_logfile.close();

    std::ifstream auxfile;
    auxfile.open("log_aux.txt");

    std::ifstream ifile("log_aux.txt");
    std::ofstream ofile("./results/log.txt", std::ios::app);

    // check to see that the input file exists:
    if (!ifile.is_open())
    {
        cout << "Aux file couldn't be found, booo\n";
    }
    // check to see that the output file exists:
    else if (!ofile.is_open())
    {
        cout << "Log file couldn't be found, booo\n";
    }
    else
    {
        ofile << "\n\n";
        ofile << "List of the T1 transitions\n";
        ofile << "T" << std::setw(20) << "v_id1" << std::setw(20) << "v_id2" << std::setw(20) << "e_length\n";
        ofile << ifile.rdbuf();
        // then add more lines to the file if need be...
    }

    ofile.close();
    auxfile.close();
    // std::filesystem::remove(auxfile);

    // Print final configuration information: vertex coordinates
    ofstream final_vertex;
    final_vertex.open("./results/vertex_coord.txt");

    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < vertices[0].size(); j++)
        {
            if (j == vertices[0].size() - 1)
            {
                final_vertex << vertices[i][j];
            }
            else
            {
                final_vertex << vertices[i][j] << ' ';
            }
        }
        final_vertex << '\n';
    }

    final_vertex.close();

    // Print final configuration information: edges indices
    ofstream final_edge;
    final_edge.open("./results/edges_index.txt");

    for (int i = 0; i < edges.size(); i++)
    {
        for (int j = 0; j < edges[0].size(); j++)
        {
            if (j == edges[0].size() - 1)
            {
                final_edge << edges[i][j];
            }
            else
            {
                final_edge << edges[i][j] << '\t';
            }
        }
        final_edge << '\n';
    }

    final_edge.close();

    // Print final configuration information: cell indices in the network
    ofstream final_network;
    final_network.open("./results/network_cells.txt");

    for (int i = 0; i < network.size(); i++)
    {
        for (int j = 0; j < network[i].indices.size(); j++)
        {
            if (j == network[i].indices.size() - 1)
            {
                final_network << network[i].indices[j];
            }
            else
            {
                final_network << network[i].indices[j] << '\t';
            }
        }
        final_network << '\n';
    }

    final_network.close();

    // Delete aux file
    if (remove("log_aux.txt") != 0)
        perror("Error deleting transitions aux file");
    else
        puts("Transitions aux file successfully deleted");
}