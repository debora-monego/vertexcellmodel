#include <iostream>
#include <fstream>
#include <tuple>

#include "transitionT1.hpp"
#include "geometry.hpp"
#include "energy.hpp"
#include "vector.hpp"
#include "polygon.hpp"

using namespace std;

// Get vertices associated with the transition
std::vector<int> get_vertex_indices(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids)
{
	// Make a copy of polygons
	std::vector<Polygon> network_copy(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_copy[i] = network[cell_ids[i]];
	}

	// Define polygons
	Polygon cell_0 = network_copy[0];
	Polygon cell_1 = network_copy[1];
	Polygon cell_2 = network_copy[2];
	Polygon cell_3 = network_copy[3];

	// Find indices with respect to cell 1
	int position;
	std::vector<int>::iterator it;
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	// i3: cell 1 - vertex before i1 (clockwise)
	int i_cw;
	if (position == 0)
	{
		i_cw = cell_1.indices.size() - 1;
	}
	else
	{
		i_cw = position - 1;
	}
	int i3 = cell_1.indices[i_cw];
	// i4: cell 1 - vertex after i1
	int i_ccw;
	if (position == (cell_1.indices.size() - 1))
	{
		i_ccw = 0;
	}
	else
	{
		i_ccw = position + 1;
	}
	int i4 = cell_1.indices[i_ccw];

	// Find indeces with respect to cell 3
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	// i5: cell 3 - vertex before i2
	if (position == 0)
	{
		i_cw = cell_3.indices.size() - 1;
	}
	else
	{
		i_cw = position - 1;
	}
	int i5 = cell_3.indices[i_cw];

	// i6: cell 3 - vertex after i2
	if (position == (cell_3.indices.size() - 1))
	{
		i_ccw = 0;
	}
	else
	{
		i_ccw = position + 1;
	}
	int i6 = cell_3.indices[i_ccw];

	std::vector<int> vertex_indices;
	vertex_indices.insert(vertex_indices.end(), {i1, i2, i3, i4, i5, i6});

	return vertex_indices;
}

// Get polygons and edges associated with short bond length
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices, 
																	std::vector<std::vector<int> > edges)
{

	// Make a copy of polygons
	std::vector<Polygon> network_0(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_0[i] = network[cell_ids[i]];
	}

	// Define polygons
	Polygon cell_0 = network_0[0];
	Polygon cell_1 = network_0[1];
	Polygon cell_2 = network_0[2];
	Polygon cell_3 = network_0[3];

	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Initialize matrix of edges 
	std::vector<std::vector<int> > edges_0(10, std::vector<int>(4, 0));

	// Edge 0 : i1 - i2
	edges_0[0][0] = i1;
	edges_0[0][1] = i2;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i1 == edges[j][0] && i2 == edges[j][1])
		{
			edges_0[0][2] = edges[j][2];
			edges_0[0][3] = edges[j][3];
		}
	}

	// Edge 1 : i1 - i3
	edges_0[1][0] = i1;
	edges_0[1][1] = i3;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i1 == edges[j][0] && i3 == edges[j][1])
		{
			edges_0[1][2] = edges[j][2];
			edges_0[1][3] = edges[j][3];
		}
	}

	// Edge 2 : i1 - i4
	edges_0[2][0] = i1;
	edges_0[2][1] = i4;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i1 == edges[j][0] && i4 == edges[j][1])
		{
			edges_0[2][2] = edges[j][2];
			edges_0[2][3] = edges[j][3];
		}
	}

	// Edge 3 : i2 - i1
	edges_0[3][0] = i2;
	edges_0[3][1] = i1;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i2 == edges[j][0] && i1 == edges[j][1])
		{
			edges_0[3][2] = edges[j][2];
			edges_0[3][3] = edges[j][3];
		}
	}

	// Edge 4 : i2 - i5
	edges_0[4][0] = i2;
	edges_0[4][1] = i5;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i2 == edges[j][0] && i5 == edges[j][1])
		{
			edges_0[4][2] = edges[j][2];
			edges_0[4][3] = edges[j][3];
		}
	}

	// Edge 5 : i2 - i6
	edges_0[5][0] = i2;
	edges_0[5][1] = i6;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i2 == edges[j][0] && i6 == edges[j][1])
		{
			edges_0[5][2] = edges[j][2];
			edges_0[5][3] = edges[j][3];
		}
	}

	// Edge 6 : i3 - i2
	edges_0[6][0] = i3;
	edges_0[6][1] = i2;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i3 == edges[j][0] && i2 == edges[j][1])
		{
			edges_0[6][2] = edges[j][2];
			edges_0[6][3] = edges[j][3];
		}
	}

	// Edge 7 : i4 - i1
	edges_0[7][0] = i4;
	edges_0[7][1] = i1;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i4 == edges[j][0] && i1 == edges[j][1])
		{
			edges_0[7][2] = edges[j][2];
			edges_0[7][3] = edges[j][3];
		}
	}

	// Edge 8 : i5 - i2
	edges_0[8][0] = i5;
	edges_0[8][1] = i2;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i5 == edges[j][0] && i2 == edges[j][1])
		{
			edges_0[8][2] = edges[j][2];
			edges_0[8][3] = edges[j][3];
		}
	}

	// Edge 9 : i6 - i2
	edges_0[9][0] = i6;
	edges_0[9][1] = i2;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i6 == edges[j][0] && i2 == edges[j][1])
		{
			edges_0[9][2] = edges[j][2];
			edges_0[9][3] = edges[j][3];
		}
	}

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0_data(network_0, edges_0);

	return T1_0_data;
}

// Get cells and edges associated with cw side
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_cw(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices)
{

	// Make a copy of polygons
	std::vector<Polygon> network_cw(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_cw[i] = network[cell_ids[i]];
	}

	// Define polygons
	Polygon cell_0 = network_cw[0];
	Polygon cell_1 = network_cw[1];
	Polygon cell_2 = network_cw[2];
	Polygon cell_3 = network_cw[3];

	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Cell 0: remove i2
	int position = 0;
	std::vector<int>::iterator it;
	it = std::find(cell_0.indices.begin(), cell_0.indices.end(), i2);
	if (it != cell_0.indices.end())
	{
		position = std::distance(cell_0.indices.begin(), it);
	}
	cell_0.indices.erase(cell_0.indices.begin() + position);
	network_cw[0].indices = cell_0.indices;

	// Cell 1: insert i2 before i1
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	std::vector<int> c1_cw_vertex_indices(cell_1.indices.begin(), cell_1.indices.begin() + position);
	std::vector<int> c1_ccw_vertex_indices(cell_1.indices.begin() + position, cell_1.indices.end());
	c1_cw_vertex_indices.insert(c1_cw_vertex_indices.end(), i2);
	c1_cw_vertex_indices.insert(c1_cw_vertex_indices.end(), c1_ccw_vertex_indices.begin(), c1_ccw_vertex_indices.end());
	network_cw[1].indices = c1_cw_vertex_indices;

	// Cell 2: remove i1
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i1);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}

	cell_2.indices.erase(cell_2.indices.begin() + position);
	network_cw[2].indices = cell_2.indices;

	// Cell 3: insert i1 before i2
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	std::vector<int> c3_cw_vertex_indices(cell_3.indices.begin(), cell_3.indices.begin() + position);
	std::vector<int> c3_ccw_vertex_indices(cell_3.indices.begin() + position, cell_3.indices.end());
	c3_cw_vertex_indices.insert(c3_cw_vertex_indices.end(), i1);
	c3_cw_vertex_indices.insert(c3_cw_vertex_indices.end(), c3_ccw_vertex_indices.begin(), c3_ccw_vertex_indices.end());
	network_cw[3].indices = c3_cw_vertex_indices;

	// Edges
	// Initialize matrix of edges
	std::vector<std::vector<int> > edges_cw(10, std::vector<int>(4, 0));

	// Edge 0 : i1 - i2
	edges_cw[0][0] = i1;
	edges_cw[0][1] = i2;

	// Edge 1 : i2 - i3
	edges_cw[1][0] = i2;
	edges_cw[1][1] = i3;

	// Edge 2 : i1 - i4
	edges_cw[2][0] = i1;
	edges_cw[2][1] = i4;

	// Edge 3 : i2 - i1
	edges_cw[3][0] = i2;
	edges_cw[3][1] = i1;

	// Edge 4 : i1 - i5
	edges_cw[4][0] = i1;
	edges_cw[4][1] = i5;

	// Edge 5 : i2 - i6
	edges_cw[5][0] = i2;
	edges_cw[5][1] = i6;

	// Edge 6 : i3 - i2
	edges_cw[6][0] = i3;
	edges_cw[6][1] = i2;

	// Edge 7 : i4 - i1
	edges_cw[7][0] = i4;
	edges_cw[7][1] = i1;

	// Edge 8 : i5 - i1
	edges_cw[8][0] = i5;
	edges_cw[8][1] = i1;

	// Edge 9 : i6 - i2
	edges_cw[9][0] = i6;
	edges_cw[9][1] = i2;

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_cw_data(network_cw, edges_cw);

	return T1_cw_data;
}

// Get cells and edges associated with ccw side
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_ccw(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices)
{

	// Make a copy of polygons
	std::vector<Polygon> network_ccw(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_ccw[i] = network[cell_ids[i]];
	}

	// Define polygons
	Polygon cell_0 = network_ccw[0];
	Polygon cell_1 = network_ccw[1];
	Polygon cell_2 = network_ccw[2];
	Polygon cell_3 = network_ccw[3];

	// i1 = vertex_indices[0];
	// i2 = vertex_indices[1];
	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Cell 0: remove i1
	int position = 0;
	std::vector<int>::iterator it;
	it = std::find(cell_0.indices.begin(), cell_0.indices.end(), i1);
	if (it != cell_0.indices.end())
	{
		position = std::distance(cell_0.indices.begin(), it);
	}
	cell_0.indices.erase(cell_0.indices.begin() + position);
	network_ccw[0].indices = cell_0.indices;

	// Cell 1: insert i2 after i1
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	std::vector<int> c1_cw_vertex_indices(cell_1.indices.begin(), cell_1.indices.begin() + position + 1);
	std::vector<int> c1_ccw_vertex_indices(cell_1.indices.begin() + position + 1, cell_1.indices.end());
	c1_cw_vertex_indices.insert(c1_cw_vertex_indices.end(), i2);
	c1_cw_vertex_indices.insert(c1_cw_vertex_indices.end(), c1_ccw_vertex_indices.begin(), c1_ccw_vertex_indices.end());
	network_ccw[1].indices = c1_cw_vertex_indices;

	// Cell 2: remove i2
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i2);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}
	cell_2.indices.erase(cell_2.indices.begin() + position);
	network_ccw[2].indices = cell_2.indices;

	// Cell 3: insert i1 after i2
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	std::vector<int> c3_cw_vertex_indices(cell_3.indices.begin(), cell_3.indices.begin() + position + 1);
	std::vector<int> c3_ccw_vertex_indices(cell_3.indices.begin() + position + 1, cell_3.indices.end());
	c3_cw_vertex_indices.insert(c3_cw_vertex_indices.end(), i1);
	c3_cw_vertex_indices.insert(c3_cw_vertex_indices.end(), c3_ccw_vertex_indices.begin(), c3_ccw_vertex_indices.end());
	network_ccw[3].indices = c3_cw_vertex_indices;

	// Edges
	// Initialize matrix of edges
	std::vector<std::vector<int> > edges_ccw(10, std::vector<int>(4, 0));
	/*for (int i = 0; i < 10; i++)
	{
		std::fill(edges_ccw[i].begin(), edges_ccw[i].end(), 0);
	}*/

	// Edge 0 : i1 - i2
	edges_ccw[0][0] = i1;
	edges_ccw[0][1] = i2;

	// Edge 1 : i1 - i3
	edges_ccw[1][0] = i1;
	edges_ccw[1][1] = i3;

	// Edge 2 : i1 - i6
	edges_ccw[2][0] = i1;
	edges_ccw[2][1] = i6;

	// Edge 3 : i2 - i1
	edges_ccw[3][0] = i2;
	edges_ccw[3][1] = i1;

	// Edge 4 : i2 - i5
	edges_ccw[4][0] = i2;
	edges_ccw[4][1] = i5;

	// Edge 5 : i2 - i4
	edges_ccw[5][0] = i2;
	edges_ccw[5][1] = i4;

	// Edge 6 : i3 - i1
	edges_ccw[6][0] = i3;
	edges_ccw[6][1] = i1;

	// Edge 7 : i4 - i2
	edges_ccw[7][0] = i4;
	edges_ccw[7][1] = i2;

	// Edge 8 : i5 - i2
	edges_ccw[8][0] = i5;
	edges_ccw[8][1] = i2;

	// Edge 9 : i6 - i1
	edges_ccw[9][0] = i6;
	edges_ccw[9][1] = i1;

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_ccw_data(network_ccw, edges_ccw);

	return T1_ccw_data;
}

// Find 4 cells involved with 2 vertices
// Labeled cells 0 - 3 in counter-clockwise order
// Cell 0 and Cell 3 are neighbors
std::vector<int> get_4_cells(std::vector<Polygon> network, int i1, int i2)
{

	std::vector<int> cell_ids(4, -1);
	// std::fill(cell_ids.begin(), cell_ids.end(), -1); // Fill with -1 to cath errors later

	Polygon cell;

	// Cell 0 or Cell 2
	// Current neighboring polys
	// Cell 0 should have i1 before i2 in counter-clockwise order
	for (int i = 0; i < network.size(); i++)
	{
		cell = network[i];

		int position1 = 0;
		int position2 = 0;
		std::vector<int>::iterator it1;
		std::vector<int>::iterator it2;
		it1 = std::find(cell.indices.begin(), cell.indices.end(), i1);
		it2 = std::find(cell.indices.begin(), cell.indices.end(), i2);
		// cout << "entrou no for 1\n";

		if ((it1 != cell.indices.end()) && (it2 != cell.indices.end()))
		{
			position1 = std::distance(cell.indices.begin(), it1);
			position2 = std::distance(cell.indices.begin(), it2);
			// cout << "nao entrou no if 1\n";

			if ((position1 == 0) && (position2 == (cell.indices.size() - 1)))
			{
				position2 = -1;
			}
			if ((position2 == 0) && (position1 == (cell.indices.size() - 1)))
			{
				position1 = -1;
			}

			// if Cell 1: i1 is before i2
			if (position1 < position2)
			{
				cell_ids[0] = i;
			}
			// if Cell 3: i2 is before i1
			if (position2 < position1)
			{
				cell_ids[2] = i;
			}
		}
		// Cell 3
		if ((it2 != cell.indices.end()) && (it1 == cell.indices.end()))
		{
			cell_ids[3] = i;
		}
		// Cell 1
		if ((it1 != cell.indices.end()) && (it2 == cell.indices.end()))
		{
			cell_ids[1] = i;
		}
	}

	return cell_ids;
}

// Perform T1 transition and check the energy change
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_transition(std::vector<std::vector<double> > vertices, std::vector<Polygon> network,
																												   std::vector<std::vector<int> > edges, std::vector<double> L, double lmin,
																												   double ka, double Lambda, double gamma, double ksep)
{
	ofstream aux_logfile;
	aux_logfile.open("log_aux.txt", std::ios_base::app);

	std::vector<int> reverse;

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data;
	//std::vector<Polygon> net_tuple = std::get<0>(T1_data);

	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> edge = edges[i];
		int i1 = edge[0];
		int i2 = edge[1];
		std::vector<int> subs{i1, i2};
		std::vector<int>::const_iterator first = edge.begin() + 2;
        std::vector<int>::const_iterator last = edge.begin() + 4;
        std::vector<int> q(first, last);
		std::vector<double> v1 = vertices[i1];
		std::vector<double> vertex2 = vertices[i2];
		std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));

		double dist = get_euclidian_distance(v1[0], v1[1], v2[0], v2[1]);

		auto it = std::search(reverse.begin(), reverse.end(), subs.begin(), subs.end());
		if ((dist < lmin) && (it == reverse.end()))
		{
			std::vector<int> cell_ids = get_4_cells(network, i1, i2);
			std::vector<int>::iterator itr;
			itr = std::find(cell_ids.begin(), cell_ids.end(), -1);
			if (itr != cell_ids.end())
			{
				// cout << "not\n";
			}
			else
			{
				// cout << "entrou no else not\n";
				//  Find minimum configuration
				reverse.insert(reverse.end(), {i2, i1});

				// Find 6 indices for vertices involved in the transition
				std::vector<int> indices = get_vertex_indices(network, i1, i2, cell_ids);
				// cout << "i1 = " << indices[0] << " i2 = " << indices[1] << " i3 = " << indices[2] << " i4 = " << indices[3] << " i5 = " << indices[4] << " i6 = " << indices[5] << '\n';

				// Get original configuration
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0_data = T1_0(network, i1, i2, cell_ids, indices, edges);
				std::vector<Polygon> network_0 = T1_0_data.first;
				std::vector<std::vector<int> > edges_0 = T1_0_data.second;
				double E0 = get_total_energy(vertices, network_0, edges_0, ka, L, Lambda, gamma);
				cout << "E0 = " << E0 << "\n";

				// Get cw T1 transition
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_cw_data = T1_cw(network, i1, i2, cell_ids, indices);
				std::vector<Polygon> network_cw = T1_cw_data.first;
				std::vector<std::vector<int> > edges_cw = T1_cw_data.second;
				double E_cw = get_total_energy(vertices, network_cw, edges_cw, ka, L, Lambda, gamma);
				cout << "Ecw = " << E_cw << "\n";

				// Get ccw T1 transition
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_ccw_data = T1_ccw(network, i1, i2, cell_ids, indices);
				std::vector<Polygon> network_ccw = T1_ccw_data.first;
				std::vector<std::vector<int>> edges_ccw = T1_ccw_data.second;
				double E_ccw = get_total_energy(vertices, network_ccw, edges_ccw, ka, L, Lambda, gamma);
				cout << "Eccw\n" << E_ccw << "\n";

				std::vector<double> all_energies{E0, E_cw, E_ccw};
				// Get minimum
				auto min_energy = *std::min_element(all_energies.begin(), all_energies.end());
				std::vector<double>::iterator itd;

				itd = std::find(all_energies.begin(), all_energies.end(), min_energy);
				int min_index = -1;
				if (itd != all_energies.end())
				{
					min_index = std::distance(all_energies.begin(), itd);
					cout << "min_index = " << min_index << '\n';
				}

				// What to do?
				// Same configuration -> do nothing
				if (min_index == 0)
				{
					// std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_data(network, edges, vertices);
					T1_data = set_tuple_data(network, edges, vertices);
				}
				else if (min_index == 1)
				{
					T1_data = set_T1_cw(network, T1_cw_data, cell_ids, edges, indices, L, lmin, ksep, vertices);
					aux_logfile << "T1" << '\t' << i1 << '\t' << i2 << '\t' << dist << '\n';	
				}
				else
				{
					T1_data = set_T1_ccw(network, T1_ccw_data, cell_ids, edges, indices, L, lmin, ksep, vertices);
					aux_logfile << "T1" << '\t' << i1 << '\t' << i2 << '\t' << dist << '\n';
				}
			}
		}
	}
	for (std::vector<int> i : std::get<1>(T1_data))
	{
		for (int j : i)
		{
			cout << j << '\t';
		}
		cout << '\n';
	}
	return T1_data;
}

std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_tuple_data(std::vector<Polygon> network, std::vector<std::vector<int> > edges,
																													std::vector<std::vector<double> > vertices)
{
	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data_tuple(network, edges, vertices);
	return T1_data_tuple;
}

// Set new cell indices, vertex positions, and edges for T1 cw transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_cw(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_cw_data,
																												 std::vector<int> cell_ids, std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, std::vector<double> L,
																												 double lmin, double ksep, std::vector<std::vector<double> > vertices)
{
	std::vector<Polygon> network_cw = T1_cw_data.first;

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_cw = network_cw[i];
		network[cell_ids[i]].indices = cell_cw.indices;
	}

	int i1 = vertex_indices[0];
	int i2 = vertex_indices[1];
	int i3 = vertex_indices[2];
	int i5 = vertex_indices[4];

	cout << i1 << " " << i2 << " " << i3 << " " << i5 << "\n";
	
	// Set new edges
	
	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> edge = edges[i];
		// i1-i3 becomes i2-i3
		if (edge[0] == i1 && edge[1] == i3)
		{
			cout << "i1-i3 becomes i2-i3\n";
			cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i2;
			cout << "new = " << edges[i][0] << '\n';
		}
		// i2-i5 becomes i1-i5
		if (edge[0] == i2 && edge[1] == i5)
		{
			cout << "i2-i5 becomes i1-i5\n";
			cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i1;
			cout << "new = " << edges[i][0] << '\n';
		}
		// i3-i1 becomes i3-i2
		if (edge[0] == i3 && edge[1] == i1)
		{
			cout << "i3-i1 becomes i3-i2\n";
			cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i2;
			cout << "new = " << edges[i][1] << '\n';
		}
		// i5-i2 becomes i5-i1
		if (edge[0] == i5 && edge[1] == i2)
		{
			cout << "i5-i2 becomes i5-i1\n";
			cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i1;
			cout << "new = " << edges[i][1] << '\n';
		}
	}

	// Set new vertices positions so that i1 and i2 are separated by ksep * lmin
	double pi = atan(1) * 4;
	double angle = -pi / 2;

	std::vector<double> v1 = vertices[i1];
	std::vector<double> vertex2 = vertices[i2];
	std::vector<int> q{0, 0};
    for (int j = 0; j < edges.size(); j++)
    {
        if (edges[j][0] == i1 && edges[j][1] == i2)
        {
            std::vector<int>::const_iterator first = edges[j].begin() + 2;
            std::vector<int>::const_iterator last = edges[j].begin() + 4;
            std::vector<int> pbc(first, last);
            q[0] = q[0] + pbc[0];
            q[1] = q[1] + pbc[1];
        }
    }
	std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));

	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices[i1][0] = move_transition[0];
	vertices[i1][1] = move_transition[1];
	vertices[i2][0] = move_transition[2];
	vertices[i2][1] = move_transition[3];

	// Wrap around periodic boundaries and update periodicity vector q in edges
	// i1
	//FIXME
	if (vertices[i1][0] < 0)
	{
		// Wrap around to right
		vertices[i1][0] = vertices[i1][0] + L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][2]++;
			}
			if (i1 == edges[j][1])
			{
				edges[j][2]--;
			}
		}
	}
	else if (vertices[i1][0] > L[0])
	{
		// Wrap around to left
		vertices[i1][0] = vertices[i1][0] - L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][2]--;
			}
			if (i1 == edges[j][1])
			{
				edges[j][2]++;
			}
		}
	}

	if (vertices[i1][1] < 0) 
	{
		// Wrap around to top
		vertices[i1][1] = vertices[i1][1] + L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][3]++;
			}
			if (i1 == edges[j][1])
			{
				edges[j][3]--;
			}
		}
	}
	else if (vertices[i1][1] > L[1])
	{
		// Wrap around to bottom
		vertices[i1][1] = vertices[i1][1] - L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][3]--;
			}
			if (i1 == edges[j][1])
			{
				edges[j][3]++;
			}
		}
	}

	// i2
	if (vertices[i2][0] < 0)
	{
		// Wrap around to right
		vertices[i2][0] = vertices[i2][0] + L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2== edges[j][0])
			{
				edges[j][2]++;
			}
			if (i2 == edges[j][1])
			{
				edges[j][2]--;
			}
		}
	}
	else if (vertices[i2][0] > L[0])
	{
		// Wrap around to left
		vertices[i2][0] = vertices[i2][0] - L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][2]--;
			}
			if (i2 == edges[j][1])
			{
				edges[j][2]++;
			}
		}
	}

	if (vertices[i2][1] < 0) 
	{
		// Wrap around to top
		vertices[i2][1] = vertices[i2][1] + L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][3]++;
			}
			if (i2 == edges[j][1])
			{
				edges[j][3]--;
			}
		}
	}
	else if (vertices[i2][1] > L[1])
	{
		// Wrap around to bottom
		vertices[i2][1] = vertices[i2][1] - L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][3]--;
			}
			if (i2 == edges[j][1])
			{
				edges[j][3]++;
			}
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_CW_data(network, edges, vertices);
	return T1_CW_data;
}

// Set new cell indices and edges for T1 ccw transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_ccw(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_ccw_data, std::vector<int> cell_ids,
																												  std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, std::vector<double> L,
																												  double lmin, double ksep, std::vector<std::vector<double> > vertices)
{
	std::vector<Polygon> network_ccw = T1_ccw_data.first;

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_ccw = network_ccw[i];
		network[cell_ids[i]].indices = cell_ccw.indices;
	}

	// Set new edges
	int i1 = vertex_indices[0];
	int i2 = vertex_indices[1];
	int i4 = vertex_indices[3];
	int i6 = vertex_indices[5];

	for (int i = 0; i < edges.size(); i++)
	{
		// cout << "Set new edges entered\n";
		std::vector<int> edge = edges[i];
		// i1-i4 becomes i2-i4
		if (edge[0] == i1 && edge[1] == i4)
		{
			// cout << "i1-i4 becomes i2-i4\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i2;
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i2-i6 becomes i1-i6
		if (edge[0] == i2 && edge[1] == i6)
		{
			// cout << "i2-i6 becomes i1-i6\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i1;
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i4-i1 becomes i4-i2
		if (edge[0] == i4 && edge[1] == i1)
		{
			// cout << "i4-i1 becomes i4-i2\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i2;
			// cout << "new = " << edges[i][1] << '\n';
		}
		// i6-i2 becomes i6-i1
		if (edge[0] == i6 && edge[1] == i2)
		{
			// cout << "i6-i2 becomes i6-i1\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i1;
			// cout << "new = " << edges[i][1] << '\n';
		}
	}

	double pi = atan(1) * 4;
	double angle = pi / 2;

	std::vector<double> v1 = vertices[i1];
	std::vector<double> vertex2 = vertices[i2];
	std::vector<int> q{0, 0};
    for (int j = 0; j < edges.size(); j++)
    {
        if (edges[j][0] == i1 && edges[j][1] == i2)
        {
            std::vector<int>::const_iterator first = edges[j].begin() + 2;
            std::vector<int>::const_iterator last = edges[j].begin() + 4;
			std::vector<int> pbc(first, last);
            q[0] = q[0] + pbc[0];
            q[1] = q[1] + pbc[1];
        }
    }
	std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));

	// Set new vertices
	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices[i1][0] = move_transition[0];
	vertices[i1][1] = move_transition[1];
	vertices[i2][0] = move_transition[2];
	vertices[i2][1] = move_transition[3];

	// Wrap around periodic boundaries and update periodicity vector q in edges
	// i1
	if (vertices[i1][0] < 0)
	{
		// Wrap around to right
		vertices[i1][0] = vertices[i1][0] + L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][2]++;
			}
			if (i1 == edges[j][1])
			{
				edges[j][2]--;
			}
		}
	}
	else if (vertices[i1][0] > L[0])
	{
		// Wrap around to left
		vertices[i1][0] = vertices[i1][0] - L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][2]--;
			}
			if (i1 == edges[j][1])
			{
				edges[j][2]++;
			}
		}
	}

	if (vertices[i1][1] < 0) 
	{
		// Wrap around to top
		vertices[i1][1] = vertices[i1][1] + L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][3]++;
			}
			if (i1 == edges[j][1])
			{
				edges[j][3]--;
			}
		}
	}
	else if (vertices[i1][1] > L[1])
	{
		// Wrap around to bottom
		vertices[i1][1] = vertices[i1][1] - L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i1 == edges[j][0] || i1 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i1 == edges[j][0])
			{
				edges[j][3]--;
			}
			if (i1 == edges[j][1])
			{
				edges[j][3]++;
			}
		}
	}

	// i2
	if (vertices[i2][0] < 0)
	{
		// Wrap around to right
		vertices[i2][0] = vertices[i2][0] + L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2== edges[j][0])
			{
				edges[j][2]++;
			}
			if (i2 == edges[j][1])
			{
				edges[j][2]--;
			}
		}
	}
	else if (vertices[i2][0] > L[0])
	{
		// Wrap around to left
		vertices[i2][0] = vertices[i2][0] - L[0];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][2] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][2]--;
			}
			if (i2 == edges[j][1])
			{
				edges[j][2]++;
			}
		}
	}

	if (vertices[i2][1] < 0) 
	{
		// Wrap around to top
		vertices[i2][1] = vertices[i2][1] + L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][3]++;
			}
			if (i2 == edges[j][1])
			{
				edges[j][3]--;
			}
		}
	}
	else if (vertices[i2][1] > L[1])
	{
		// Wrap around to bottom
		vertices[i2][1] = vertices[i2][1] - L[1];
		for (int j = 0; j < edges.size(); j++){
			if (i2 == edges[j][0] || i2 == edges[j][1]){
				edges[j][3] = 0;
			}
		}
		for (int j = 0; j < edges.size(); j++)
		{
			if (i2 == edges[j][0])
			{
				edges[j][3]--;
			}
			if (i2 == edges[j][1])
			{
				edges[j][3]++;
			}
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_CCW_data(network, edges, vertices);
	return T1_CCW_data;
}