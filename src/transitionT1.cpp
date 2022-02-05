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
	// Polygon cell_0 = network_copy[cell_ids[0]];
	// Polygon cell_1 = network_copy[cell_ids[1]];
	// Polygon cell_2 = network_copy[cell_ids[2]];
	// Polygon cell_3 = network_copy[cell_ids[3]];
	Polygon cell_0 = network_copy[0];
	Polygon cell_1 = network_copy[1];
	Polygon cell_2 = network_copy[2];
	Polygon cell_3 = network_copy[3];

	// cout << cell_0.indices[0] << '\n'
	// 	 << cell_1.indices[0] << '\n'
	// 	 << cell_2.indices[0] << '\n'
	// 	 << cell_3.indices[0] << '\n';

	// Find indices with respect to cell 1
	int position;
	std::vector<int>::iterator it;
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	// i3: cell 1 - vertex before i1 (clockwise)
	int i_left;
	if (position == 0)
	{
		i_left = cell_1.indices.size() - 1;
	}
	else
	{
		i_left = position - 1;
	}
	int i3 = cell_1.indices[i_left];
	// i4: cell 1 - vertex after i1
	int i_right;
	if (position == (cell_1.indices.size() - 1))
	{
		i_right = 0;
	}
	else
	{
		i_right = position + 1;
	}
	int i4 = cell_1.indices[i_right];

	// Find indeces with respect to cell 3
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	// i5: cell 3 - vertex before i2
	if (position == 0)
	{
		i_left = cell_3.indices.size() - 1;
	}
	else
	{
		i_left = position - 1;
	}
	int i5 = cell_3.indices[i_left];

	// i6: cell 3 - vertex after i2
	if (position == (cell_3.indices.size() - 1))
	{
		i_right = 0;
	}
	else
	{
		i_right = position + 1;
	}
	int i6 = cell_3.indices[i_right];

	std::vector<int> vertex_indices;
	vertex_indices.insert(vertex_indices.end(), {i1, i2, i3, i4, i5, i6});

	return vertex_indices;
}

// Get polygons and edges associated with short bond length
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices)
{

	// Make a copy of polygons
	std::vector<Polygon> network_0(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_0[i] = network[cell_ids[i]];
	}

	// Define polygons
	// Polygon cell_0 = network_0[cell_ids[0]];
	// Polygon cell_1 = network_0[cell_ids[1]];
	// Polygon cell_2 = network_0[cell_ids[2]];
	// Polygon cell_3 = network_0[cell_ids[3]];
	Polygon cell_0 = network_0[0];
	Polygon cell_1 = network_0[1];
	Polygon cell_2 = network_0[2];
	Polygon cell_3 = network_0[3];

	// i1 = vertex_indices[0];
	// i2 = vertex_indices[1];
	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Initialize matrix of edges
	// std::vector<std::vector<int> > edges_0(edges);
	std::vector<std::vector<int> > edges_0(10, std::vector<int>(2, 0));
	/*for (int i = 0; i < 10; i++)
	{
		std::fill(edges_0[i].begin(), edges_0[i].end(), 0);
	}*/

	// Edge 0 : i1 - i2
	edges_0[0][0] = i1;
	edges_0[0][1] = i2;

	// Edge 1 : i1 - i3
	edges_0[1][0] = i1;
	edges_0[1][1] = i3;

	// Edge 2 : i1 - i4
	edges_0[2][0] = i1;
	edges_0[2][1] = i4;

	// Edge 3 : i2 - i1
	edges_0[3][0] = i2;
	edges_0[3][1] = i1;

	// Edge 4 : i2 - i5
	edges_0[4][0] = i2;
	edges_0[4][1] = i5;

	// Edge 5 : i2 - i6
	edges_0[5][0] = i2;
	edges_0[5][1] = i6;

	// Edge 6 : i3 - i2
	edges_0[6][0] = i3;
	edges_0[6][1] = i2;

	// Edge 7 : i4 - i1
	edges_0[7][0] = i4;
	edges_0[7][1] = i1;

	// Edge 8 : i5 - i2
	edges_0[8][0] = i5;
	edges_0[8][1] = i2;

	// Edge 9 : i6 - i2
	edges_0[9][0] = i6;
	edges_0[9][1] = i2;

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0_data(network_0, edges_0);

	return T1_0_data;
}

// Get cells and edges associated with left side
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_left(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices)
{

	// Make a copy of polygons
	std::vector<Polygon> network_l(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_l[i] = network[cell_ids[i]];
	}

	// Define polygons
	// Polygon cell_0 = network_l[cell_ids[0]];
	// Polygon cell_1 = network_l[cell_ids[1]];
	// Polygon cell_2 = network_l[cell_ids[2]];
	// Polygon cell_3 = network_l[cell_ids[3]];
	Polygon cell_0 = network_l[0];
	Polygon cell_1 = network_l[1];
	Polygon cell_2 = network_l[2];
	Polygon cell_3 = network_l[3];

	// i1 = vertex_indices[0];
	// i2 = vertex_indices[1];
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
	network_l[0].indices = cell_0.indices;

	// Cell 1: insert i2 before i1
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	std::vector<int> c1_left_vertex_indices(cell_1.indices.begin(), cell_1.indices.begin() + position);
	std::vector<int> c1_right_vertex_indices(cell_1.indices.begin() + position, cell_1.indices.end());
	c1_left_vertex_indices.insert(c1_left_vertex_indices.end(), i2);
	c1_left_vertex_indices.insert(c1_left_vertex_indices.end(), c1_right_vertex_indices.begin(), c1_right_vertex_indices.end());
	network_l[1].indices = c1_left_vertex_indices;

	// Cell 2: remove i1
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i1);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}

	cell_2.indices.erase(cell_2.indices.begin() + position);
	network_l[2].indices = cell_2.indices;

	// Cell 3: insert i1 before i2
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	std::vector<int> c3_left_vertex_indices(cell_3.indices.begin(), cell_3.indices.begin() + position);
	std::vector<int> c3_right_vertex_indices(cell_3.indices.begin() + position, cell_3.indices.end());
	c3_left_vertex_indices.insert(c3_left_vertex_indices.end(), i1);
	c3_left_vertex_indices.insert(c3_left_vertex_indices.end(), c3_right_vertex_indices.begin(), c3_right_vertex_indices.end());
	network_l[3].indices = c3_left_vertex_indices;

	// Edges
	// Initialize matrix of edges
	std::vector<std::vector<int> > edges_l(10, std::vector<int>(2, 0));
	/*for (int i = 0; i < 10; i++)
	{
		std::fill(edges_l[i].begin(), edges_l[i].end(), 0);
	}*/

	// Edge 0 : i1 - i2
	edges_l[0][0] = i1;
	edges_l[0][1] = i2;

	// Edge 1 : i2 - i3
	edges_l[1][0] = i2;
	edges_l[1][1] = i3;

	// Edge 2 : i1 - i4
	edges_l[2][0] = i1;
	edges_l[2][1] = i4;

	// Edge 3 : i2 - i1
	edges_l[3][0] = i2;
	edges_l[3][1] = i1;

	// Edge 4 : i1 - i5
	edges_l[4][0] = i1;
	edges_l[4][1] = i5;

	// Edge 5 : i2 - i6
	edges_l[5][0] = i2;
	edges_l[5][1] = i6;

	// Edge 6 : i3 - i2
	edges_l[6][0] = i3;
	edges_l[6][1] = i2;

	// Edge 7 : i4 - i1
	edges_l[7][0] = i4;
	edges_l[7][1] = i1;

	// Edge 8 : i5 - i1
	edges_l[8][0] = i5;
	edges_l[8][1] = i1;

	// Edge 9 : i6 - i2
	edges_l[9][0] = i6;
	edges_l[9][1] = i2;

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_l_data(network_l, edges_l);

	return T1_l_data;
}

// Get cells and edges associated with right side
std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_right(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices)
{

	// Make a copy of polygons
	std::vector<Polygon> network_r(4);
	for (int i = 0; i < cell_ids.size(); i++)
	{
		network_r[i] = network[cell_ids[i]];
	}

	// Define polygons
	Polygon cell_0 = network_r[0];
	Polygon cell_1 = network_r[1];
	Polygon cell_2 = network_r[2];
	Polygon cell_3 = network_r[3];

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
	network_r[0].indices = cell_0.indices;

	// Cell 1: insert i2 after i1
	it = std::find(cell_1.indices.begin(), cell_1.indices.end(), i1);
	if (it != cell_1.indices.end())
	{
		position = std::distance(cell_1.indices.begin(), it);
	}
	std::vector<int> c1_left_vertex_indices(cell_1.indices.begin(), cell_1.indices.begin() + position + 1);
	std::vector<int> c1_right_vertex_indices(cell_1.indices.begin() + position + 1, cell_1.indices.end());
	c1_left_vertex_indices.insert(c1_left_vertex_indices.end(), i2);
	c1_left_vertex_indices.insert(c1_left_vertex_indices.end(), c1_right_vertex_indices.begin(), c1_right_vertex_indices.end());
	network_r[1].indices = c1_left_vertex_indices;

	// Cell 2: remove i2
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i2);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}
	cell_2.indices.erase(cell_2.indices.begin() + position);
	network_r[2].indices = cell_2.indices;

	// Cell 3: insert i1 after i2
	it = std::find(cell_3.indices.begin(), cell_3.indices.end(), i2);
	if (it != cell_3.indices.end())
	{
		position = std::distance(cell_3.indices.begin(), it);
	}
	std::vector<int> c3_left_vertex_indices(cell_3.indices.begin(), cell_3.indices.begin() + position + 1);
	std::vector<int> c3_right_vertex_indices(cell_3.indices.begin() + position + 1, cell_3.indices.end());
	c3_left_vertex_indices.insert(c3_left_vertex_indices.end(), i1);
	c3_left_vertex_indices.insert(c3_left_vertex_indices.end(), c3_right_vertex_indices.begin(), c3_right_vertex_indices.end());
	network_r[3].indices = c3_left_vertex_indices;

	// Edges
	// Initialize matrix of edges
	std::vector<std::vector<int> > edges_r(10, std::vector<int>(2, 0));
	/*for (int i = 0; i < 10; i++)
	{
		std::fill(edges_r[i].begin(), edges_r[i].end(), 0);
	}*/

	// Edge 0 : i1 - i2
	edges_r[0][0] = i1;
	edges_r[0][1] = i2;

	// Edge 1 : i1 - i3
	edges_r[1][0] = i1;
	edges_r[1][1] = i3;

	// Edge 2 : i1 - i6
	edges_r[2][0] = i1;
	edges_r[2][1] = i6;

	// Edge 3 : i2 - i1
	edges_r[3][0] = i2;
	edges_r[3][1] = i1;

	// Edge 4 : i2 - i5
	edges_r[4][0] = i2;
	edges_r[4][1] = i5;

	// Edge 5 : i2 - i4
	edges_r[5][0] = i2;
	edges_r[5][1] = i4;

	// Edge 6 : i3 - i1
	edges_r[6][0] = i3;
	edges_r[6][1] = i1;

	// Edge 7 : i4 - i2
	edges_r[7][0] = i4;
	edges_r[7][1] = i2;

	// Edge 8 : i5 - i2
	edges_r[8][0] = i5;
	edges_r[8][1] = i2;

	// Edge 9 : i6 - i1
	edges_r[9][0] = i6;
	edges_r[9][1] = i1;

	std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_r_data(network_r, edges_r);

	return T1_r_data;
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
	// for (int i = 0; i < cell_ids.size(); i++)
	// {
	// 	cout << cell_ids[i] << '\n';
	// }
	// cout << "i1 = " << i1 << " i2 =" << i2 << "\n";
	return cell_ids;
}

// Perform T1 transition and check the energy change
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_transition(std::vector<std::vector<double> > vertices, std::vector<Polygon> network,
																												   std::vector<std::vector<int> > edges, double lx, double ly, double lmin,
																												   double ka, double Lambda, double gamma, double ksep)
{
	ofstream aux_logfile;
	aux_logfile.open("log_aux.txt", std::ios_base::app);

	std::vector<double> L{lx, ly};
	std::vector<int> reverse;

	// std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data(network, edges, vertices);
	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data; //(network, edges, vertices);

	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> edge = edges[i];
		int i1 = edge[0];
		int i2 = edge[1];
		std::vector<int> subs{i1, i2};

		std::vector<double> v1 = vertices[i1];
		std::vector<double> vertex2 = vertices[i2];
		std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L));

		double dist = get_euclidian_distance(v1[0], v1[1], v2[0], v2[1]);

		auto it = std::search(reverse.begin(), reverse.end(), subs.begin(), subs.end());
		if ((dist < lmin) && (it == reverse.end()))
		{
			// cout << "entrou 1 if\n";
			aux_logfile << "T1" << std::setw(20) << i1 << std::setw(20) << i2 << std::setw(30) << dist << '\n';
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
				cout << "i1 = " << indices[0] << " i2 = " << indices[1] << " i3 = " << indices[2] << " i4 = " << indices[3] << " i5 = " << indices[4] << " i6 = " << indices[5] << '\n';

				// Get original configuration
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_0_data = T1_0(network, i1, i2, cell_ids, indices);
				std::vector<Polygon> network_0 = T1_0_data.first;
				std::vector<std::vector<int> > edges_0 = T1_0_data.second;
				double E0 = get_total_energy(vertices, network_0, edges_0, ka, L, Lambda, gamma);

				// Get left T1 transition
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_l_data = T1_left(network, i1, i2, cell_ids, indices);
				std::vector<Polygon> network_l = T1_l_data.first;
				std::vector<std::vector<int> > edges_l = T1_l_data.second;
				double E_left = get_total_energy(vertices, network_l, edges_l, ka, L, Lambda, gamma);

				// Get right T1 transition
				std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_r_data = T1_right(network, i1, i2, cell_ids, indices);
				std::vector<Polygon> network_r = T1_r_data.first;
				std::vector<std::vector<int> > edges_r = T1_r_data.second;
				double E_right = get_total_energy(vertices, network_r, edges_r, ka, L, Lambda, gamma);

				std::vector<double> all_energies{E0, E_left, E_right};
				// Get minimum
				auto min_energy = *std::min_element(all_energies.begin(), all_energies.end());
				std::vector<double>::iterator itd;

				itd = std::find(all_energies.begin(), all_energies.end(), min_energy);
				int min_index = -1;
				if (itd != all_energies.end())
				{
					min_index = std::distance(all_energies.begin(), itd);
				}

				// What to do?
				// Same configuration -> do nothing
				if (min_index == 0)
				{
					// cout << "same configuration\n";
				}
				else if (min_index == 1)
				{
					T1_data = set_T1_left(network, T1_l_data, cell_ids, edges, indices, lx, ly, lmin, ksep, vertices);
				}
				else
				{
					T1_data = set_T1_right(network, T1_r_data, cell_ids, edges, indices, lx, ly, lmin, ksep, vertices);
				}
			}
		}
	}
	return T1_data;
}

// Set new cell indices, vertex positions, and edges for T1 left transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_left(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_l_data,
																												 std::vector<int> cell_ids, std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, double lx, double ly,
																												 double lmin, double ksep, std::vector<std::vector<double> > vertices)
{
	std::vector<Polygon> network_l = T1_l_data.first;

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_left = network_l[i];
		network[cell_ids[i]].indices = cell_left.indices;
	}

	// Set new edges
	int i1 = vertex_indices[0];
	int i2 = vertex_indices[1];
	int i3 = vertex_indices[2];
	int i5 = vertex_indices[4];
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

	std::vector<double> L{lx, ly};

	std::vector<double> v1 = vertices[i1];
	std::vector<double> vertex2 = vertices[i2];
	std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L));

	// Set new vertices positions so that i1 and i2 are separated by ksep * lmin
	std::vector<double> move_v1_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin);
	vertices[i1][0] = move_v1_transition[0];
	vertices[i1][1] = move_v1_transition[1];

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_left_data(network, edges, vertices);
	return T1_left_data;
}

// Set new cell indices and edges for T1 right transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > set_T1_right(std::vector<Polygon> network, std::pair<std::vector<Polygon>, std::vector<std::vector<int> > > T1_r_data, std::vector<int> cell_ids,
																												  std::vector<std::vector<int> > edges, std::vector<int> vertex_indices, double lx, double ly,
																												  double lmin, double ksep, std::vector<std::vector<double> > vertices)
{
	std::vector<Polygon> network_r = T1_r_data.first;

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_right = network_r[i];
		network[cell_ids[i]].indices = cell_right.indices;
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
			cout << "i1-i4 becomes i2-i4\n";
			cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i2;
			cout << "new = " << edges[i][0] << '\n';
		}
		// i2-i6 becomes i1-i6
		if (edge[0] == i2 && edge[1] == i6)
		{
			cout << "i2-i6 becomes i1-i6\n";
			cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i1;
			cout << "new = " << edges[i][0] << '\n';
		}
		// i4-i1 becomes i4-i2
		if (edge[0] == i4 && edge[1] == i1)
		{
			cout << "i4-i1 becomes i4-i2\n";
			cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i2;
			cout << "new = " << edges[i][1] << '\n';
		}
		// i6-i2 becomes i6-i1
		if (edge[0] == i6 && edge[1] == i2)
		{
			cout << "i6-i2 becomes i6-i1\n";
			cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i1;
			cout << "new = " << edges[i][1] << '\n';
		}
	}

	std::vector<double> L{lx, ly};

	std::vector<double> v1 = vertices[i1];
	std::vector<double> vertex2 = vertices[i2];
	std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L));

	// Set new vertices
	std::vector<double> move_v1_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin);
	vertices[i1][0] = move_v1_transition[0];
	vertices[i1][1] = move_v1_transition[1];

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_right_data(network, edges, vertices);
	return T1_right_data;
}