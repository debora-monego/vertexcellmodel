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
// Example below. Note that cells do not need to be all hexagons necessarily, as all 6 indices of interest will be
// directly connected to i1 and i2, i.e., the edge going through T1
//      / \
//  i3 | 1 | i4
//    / \ / \     i1 (center)
//   | 2 | 0 |
//	  \ / \ /     i2 (center)
//  i6 | 3 | i5
//      \ /
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

	// Find indices with respect to cell 3
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

// Get polygons and edges associated with original short bond length
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
	std::vector<std::vector<int>> edges_0(10, std::vector<int>(4, 0));

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

	// Edge 6 : i3 - i1
	edges_0[6][0] = i3;
	edges_0[6][1] = i1;
	for (int j = 0; j < edges.size(); j++)
	{
		if (i3 == edges[j][0] && i1 == edges[j][1])
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
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_cw(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices,
																										std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, std::vector<double> L,
																										double ksep, double lmin)
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
	// Update edges
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

	// Vertices
	// Initialize copy of vertices
	std::vector<std::vector<double> > vertices_cw(vertices.size(), std::vector<double>(2));
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices_cw[i] = vertices[i];
	}
	//cout << vertices_cw.size() << "\n";
	// for (std::vector<double> i : vertices_cw)
	// {
	// 	for (double j : i)
	// 	{
	// 		cout << j << '\t';
	// 	}
	// 	cout << '\n';
	// }

	// Set new vertices positions so that i1 and i2 are separated by ksep * lmin
	double pi = atan(1) * 4;
	double angle = -pi / 2;

	std::vector<double> v1 = vertices_cw[i1];
	std::vector<double> vertex2 = vertices_cw[i2];
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

	// cout << "set separation cw\n";
	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices_cw[i1][0] = move_transition[0];
	vertices_cw[i1][1] = move_transition[1];
	vertices_cw[i2][0] = move_transition[2];
	vertices_cw[i2][1] = move_transition[3];

	// std::cout << vertices_cw[i1][0] << " " << vertices_cw[i1][1] << " " << vertices_cw[i2][0] << " " << vertices_cw[i2][1] << "\n";

	// Wrap around periodic boundaries and update periodicity vector q in edges

	int qx = q[0];
	int qy = q[1];
	for (int j = 0; j < edges_cw.size(); j++)
	{
		if ((edges_cw[j][0] == i1 && edges_cw[j][1] == i2) || (edges_cw[j][0] == i2 && edges_cw[j][1] == i1)) 
		{
			if ((vertices_cw[edges_cw[j][0]][0] > 0 && vertices_cw[edges_cw[j][1]][0] < 0) || (vertices_cw[edges_cw[j][0]][0] > L[0] && vertices_cw[edges_cw[j][1]][0] > 0))
			{
				edges_cw[j][2]--;
				if (qx == 0){
					qx = edges_cw[j][2];
				}
				
			}
			else if ((vertices_cw[edges_cw[j][0]][0] < 0 && vertices_cw[edges_cw[j][1]][0] > 0) || (vertices_cw[edges_cw[j][0]][0] > 0 && vertices_cw[edges_cw[j][1]][0] > L[0]))
			{
				edges_cw[j][2]++;
				if (qx == 0){
					qx = edges_cw[j][2];
				}
			}
			if ((vertices_cw[edges_cw[j][0]][1] > 0 && vertices_cw[edges_cw[j][1]][1] < 0) || (vertices_cw[edges_cw[j][0]][1] > L[1] && vertices_cw[edges_cw[j][1]][1] > 0))
			{
				edges_cw[j][3]--;
				if (qy == 0){
					qy = edges_cw[j][3];
				}
			}
			else if ((vertices_cw[edges_cw[j][0]][1] < 0 && vertices_cw[edges_cw[j][1]][1] > 0) || (vertices_cw[edges_cw[j][0]][1] > 0 && vertices_cw[edges_cw[j][1]][1] > L[1]))
			{
				edges_cw[j][3]++;
				if (qy == 0){
					qy = edges_cw[j][3];
				}
			}
		} else if (edges_cw[j][0] == i1 || edges_cw[j][0] == i2)
		{
			cout << "entered else if\n";
			//cout << "test" << edges_cw[j][0] << " " << edges_cw[j][1] << "\n";
			std::vector<double> vertex2_cw = vertices_cw[edges_cw[j][1]];
			std::vector<int> q_cw{0, 0};
			for (int k = 0; k < edges.size(); k++)
			{
				if ((edges[k][0] == i1 || edges[k][0] == i2) && edges[k][1] == edges_cw[j][1])
				{
					std::vector<int>::const_iterator first = edges[k].begin() + 2;
					std::vector<int>::const_iterator last = edges[k].begin() + 4;
					std::vector<int> pbc_cw(first, last);
					q_cw[0] = q_cw[0] + pbc_cw[0] + qx;
					q_cw[1] = q_cw[1] + pbc_cw[1] + qy;
				}
			}
			std::vector<double> v2_cw = add_vectors(vertices_cw[edges_cw[j][0]], pbc_diff(vertex2_cw, vertices_cw[edges_cw[j][0]], L, q_cw));
			// cout << "coord = " << v2_cw[0] << " " << v2_cw[1] << "\n";

			if ((vertices_cw[edges_cw[j][0]][0] > 0 && v2_cw[0] < 0) || (vertices_cw[edges_cw[j][0]][0] > L[0] && v2_cw[0] > 0))
			{
				edges_cw[j][2]--;
			}
			else if ((vertices_cw[edges_cw[j][0]][0] < 0 && v2_cw[0] > 0) || (vertices_cw[edges_cw[j][0]][0] > 0 && v2_cw[0] > L[0]))
			{
				edges_cw[j][2]++;
			}
			if ((vertices_cw[edges_cw[j][0]][1] > 0 && v2_cw[1] < 0) || (vertices_cw[edges_cw[j][0]][1] > L[1] && v2_cw[1] > 0))
			{
				edges_cw[j][3]--;
			}
			else if ((vertices_cw[edges_cw[j][0]][1] < 0 && v2_cw[1] > 0) || (vertices_cw[edges_cw[j][0]][1] > 0 && v2_cw[1] > L[1]))
			{
				cout<< "entered right if\n";
				edges_cw[j][3]++;
			}
		} else if (edges_cw[j][1] == i1 || edges_cw[j][1] == i2)
		{
			std::vector<double> vertex2_cw = vertices_cw[edges_cw[j][0]];
			std::vector<int> q_cw{0, 0};
			for (int k = 0; k < edges.size(); k++)
			{
				if (edges[k][0] == edges_cw[j][0] && (edges[k][1] == i1 || edges[k][1] == i2))
				{
					std::vector<int>::const_iterator first = edges[k].begin() + 2;
					std::vector<int>::const_iterator last = edges[k].begin() + 4;
					std::vector<int> pbc_cw(first, last);
					q_cw[0] = q_cw[0] + pbc_cw[0] + qx;
					q_cw[1] = q_cw[1] + pbc_cw[1] + qy;
				}
			}
			std::vector<double> v2_cw = add_vectors(vertices_cw[edges_cw[j][1]], pbc_diff(vertex2_cw, vertices_cw[edges_cw[j][1]], L, q_cw));

			if ((vertices_cw[edges_cw[j][1]][0] > 0 && v2_cw[0] < 0) || (vertices_cw[edges_cw[j][1]][0] > L[0] && v2_cw[0] > 0))
			{
				edges_cw[j][2]++;
			}
			else if ((vertices_cw[edges_cw[j][1]][0] < 0 && v2_cw[0] > 0) || (vertices_cw[edges_cw[j][1]][0] > 0 && v2_cw[0] > L[0]))
			{
				edges_cw[j][2]--;
			}
			if ((vertices_cw[edges_cw[j][1]][1] > 0 && v2_cw[1] < 0) || (vertices_cw[edges_cw[j][1]][1] > L[1] && v2_cw[1] > 0))
			{
				edges_cw[j][3]++;
			}
			else if ((vertices_cw[edges_cw[j][1]][1] < 0 && v2_cw[1] > 0) || (vertices_cw[edges_cw[j][1]][1] > 0 && v2_cw[1] > L[1]))
			{
				edges_cw[j][3]--;
			}
		}
	}

	// Update position of vertices i1 and i2 with pbc
	if (vertices_cw[i1][0] < 0)
	{
		vertices_cw[i1][0] = vertices_cw[i1][0] + L[0];
	}
	else if (vertices_cw[i1][0] > L[0])
	{
		vertices_cw[i1][0] = vertices_cw[i1][0] - L[0];
	}
	if (vertices_cw[i1][1] < 0)
	{
		vertices_cw[i1][1] = vertices_cw[i1][1] + L[1];
	}
	else if (vertices_cw[i1][1] > L[1])
	{
		vertices_cw[i1][1] = vertices_cw[i1][1] - L[1];
	}
	if (vertices_cw[i2][0] < 0)
	{
		vertices_cw[i2][0] = vertices_cw[i2][0] + L[0];
	}
	else if (vertices_cw[i2][0] > L[0])
	{
		vertices_cw[i2][0] = vertices_cw[i2][0] - L[0];
	}
	if (vertices_cw[i2][1] < 0)
	{
		vertices_cw[i2][1] = vertices_cw[i2][1] + L[1];
	}
	else if (vertices_cw[i2][1] > L[1])
	{
		vertices_cw[i2][1] = vertices_cw[i2][1] - L[1];
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_cw_data(network_cw, edges_cw, vertices_cw);

	return T1_cw_data;
}

// Get cells and edges associated with ccw side
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_ccw(std::vector<Polygon> network, int i1, int i2, std::vector<int> cell_ids, std::vector<int> vertex_indices,
																										std::vector<std::vector<double> > vertices, std::vector<std::vector<int> > edges, std::vector<double> L,
																										double ksep, double lmin)
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
	std::vector<std::vector<int>> edges_ccw(10, std::vector<int>(4, 0));
	// Update edges
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

	// Vertices
	// Initialize copy of vertices
	std::vector<std::vector<double> > vertices_ccw(vertices.size(), std::vector<double>(2));
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices_ccw[i] = vertices[i];
	}

	// Set new vertices positions so that i1 and i2 are separated by ksep * lmin
	double pi = atan(1) * 4;
	double angle = pi / 2;

	std::vector<double> v1 = vertices_ccw[i1];
	std::vector<double> vertex2 = vertices_ccw[i2];
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

	// cout << "set separation ccw\n";
	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices_ccw[i1][0] = move_transition[0];
	vertices_ccw[i1][1] = move_transition[1];
	vertices_ccw[i2][0] = move_transition[2];
	vertices_ccw[i2][1] = move_transition[3];

	// Wrap around periodic boundaries and update periodicity vector q in edges 
	//FIXME maybe I have to change this so that the vertices go to the same point (i1=i2) and then separate into the new distance
	// i1
	int qx = q[0];
	int qy = q[1];
	for (int j = 0; j < edges_ccw.size(); j++)
	{
		cout << "entered first loop \n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n" << vertices_ccw[edges_ccw[j][0]][0] << " " << vertices_ccw[edges_ccw[j][1]][0];
		if ((edges_ccw[j][0] == i1 && edges_ccw[j][1] == i2) || (edges_ccw[j][0] == i2 && edges_ccw[j][1] == i1))
		{
			if ((vertices_ccw[edges_ccw[j][0]][0] > 0 && vertices_ccw[edges_ccw[j][1]][0] < 0) || (vertices_ccw[edges_ccw[j][0]][0] > L[0] && vertices_ccw[edges_ccw[j][1]][0] > 0))
			{
				edges_ccw[j][2]--;
				if (qx == 0){
					qx = edges_ccw[j][2];
				}
			}
			else if ((vertices_ccw[edges_ccw[j][0]][0] < 0 && vertices_ccw[edges_ccw[j][1]][0] > 0) || (vertices_ccw[edges_ccw[j][0]][0] > 0 && vertices_ccw[edges_ccw[j][1]][0] > L[0]))
			{
				edges_ccw[j][2]++;
				if (qx == 0){
					qx = edges_ccw[j][2];
				}
			}
			if ((vertices_ccw[edges_ccw[j][0]][1] > 0 && vertices_ccw[edges_ccw[j][1]][1] < 0) || (vertices_ccw[edges_ccw[j][0]][1] > L[1] && vertices_ccw[edges_ccw[j][1]][1] > 0))
			{
				edges_ccw[j][3]--;
				if (qy == 0){
					qy = edges_ccw[j][3];
				}
			}
			else if ((vertices_ccw[edges_ccw[j][0]][1] < 0 && vertices_ccw[edges_ccw[j][1]][1] > 0) || (vertices_ccw[edges_ccw[j][0]][1] > 0 && vertices_ccw[edges_ccw[j][1]][1] > L[1]))
			{
				edges_ccw[j][3]++;
				if (qy == 0){
					qy = edges_ccw[j][3];
				}
			}
			cout << "q " << edges_ccw[j][2] << " " << edges_ccw[j][3] << "\n";
		} else if (edges_ccw[j][0] == i1 || edges_ccw[j][0] == i2)
		{
			cout << "entered second loop\n";
			std::vector<double> vertex2_ccw = vertices_ccw[edges_ccw[j][1]];
			std::vector<int> q_ccw{0, 0};
			for (int k = 0; k < edges.size(); k++)
			{
				if ((edges[k][0] == i1 || edges[k][0] == i2) && edges[k][1] == edges_ccw[j][1])
				{
					std::vector<int>::const_iterator first = edges[k].begin() + 2;
					std::vector<int>::const_iterator last = edges[k].begin() + 4;
					std::vector<int> pbc_ccw(first, last);
					q_ccw[0] = q_ccw[0] + pbc_ccw[0] + qx;
					q_ccw[1] = q_ccw[1] + pbc_ccw[1] + qy;
				}
			}
			std::vector<double> v2_ccw = add_vectors(vertices_ccw[edges_ccw[j][0]], pbc_diff(vertex2_ccw, vertices_ccw[edges_ccw[j][0]], L, q_ccw));
			cout << "v1 coord = " << vertices_ccw[edges_ccw[j][0]][0] << " " << vertices_ccw[edges_ccw[j][0]][1] << " v2 coord = " << v2_ccw[0] << " " << v2_ccw[1] << "\n";

			if ((vertices_ccw[edges_ccw[j][0]][0] > 0 && v2_ccw[0] < 0) || (vertices_ccw[edges_ccw[j][0]][0] > L[0] && v2_ccw[0] > 0))
			{
				cout << "entered cww loop 2\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][2]--;
			}
			else if ((vertices_ccw[edges_ccw[j][0]][0] < 0 && v2_ccw[0] > 0) || (vertices_ccw[edges_ccw[j][0]][0] > 0 && v2_ccw[0] > L[0]))
			{
				cout << "entered cww loop 3\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][2]++;
			}
			if ((vertices_ccw[edges_ccw[j][0]][1] > 0 && v2_ccw[1] < 0) || (vertices_ccw[edges_ccw[j][0]][1] > L[1] && v2_ccw[1] > 0))
			{
				cout << "entered cww loop 4\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][3]--;
			}
			else if ((vertices_ccw[edges_ccw[j][0]][1] < 0 && v2_ccw[1] > 0) || (vertices_ccw[edges_ccw[j][0]][1] > 0 && v2_ccw[1] > L[1]))
			{
				cout << "entered cww loop 5\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][3]++;
			}
			cout << "q " << edges_ccw[j][2] << " " << edges_ccw[j][3] << "\n";
		} else if (edges_ccw[j][1] == i1 || edges_ccw[j][1] == i2)
		{
			std::vector<double> vertex2_ccw = vertices_ccw[edges_ccw[j][0]];
			std::vector<int> q_ccw{0, 0};
			for (int k = 0; k < edges.size(); k++)
			{
				if (edges[k][0] == edges_ccw[j][0] && (edges[k][1] == i1 || edges[k][1] == i2))
				{
					std::vector<int>::const_iterator first = edges[k].begin() + 2;
					std::vector<int>::const_iterator last = edges[k].begin() + 4;
					std::vector<int> pbc_ccw(first, last);
					q_ccw[0] = q_ccw[0] + pbc_ccw[0] + qx;
					q_ccw[1] = q_ccw[1] + pbc_ccw[1] + qy;
					cout << pbc_ccw[0] << " qx = " << q[0] << " " << pbc_ccw[1] << " qy = " << q[1] << "\n";
				}
			}
			std::vector<double> v2_ccw = add_vectors(vertices_ccw[edges_ccw[j][1]], pbc_diff(vertex2_ccw, vertices_ccw[edges_ccw[j][1]], L, q_ccw));
			cout << "v1 coord = " << vertices_ccw[edges_ccw[j][1]][0] << " " << vertices_ccw[edges_ccw[j][1]][1] << " v2 coord = " << v2_ccw[0] << " " << v2_ccw[1] << "\n";
			cout << "q " << edges_ccw[j][2] << " " << edges_ccw[j][3] << "\n";

			if ((vertices_ccw[edges_ccw[j][1]][0] > 0 && v2_ccw[0] < 0) || (vertices_ccw[edges_ccw[j][1]][0] > L[0] && v2_ccw[0] > 0))
			{
				cout << "entered cww loop 20\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][2]++;
			}
			else if ((vertices_ccw[edges_ccw[j][1]][0] < 0 && v2_ccw[0] > 0) || (vertices_ccw[edges_ccw[j][1]][0] > 0 && v2_ccw[0] > L[0]))
			{
				cout << "entered cww loop 30\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][2]--;
			}
			if ((vertices_ccw[edges_ccw[j][1]][1] > 0 && v2_ccw[1] < 0) || (vertices_ccw[edges_ccw[j][1]][1] > L[1] && v2_ccw[1] > 0))
			{
				cout << "entered cww loop 40\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][3]++;
			}
			else if ((vertices_ccw[edges_ccw[j][1]][1] < 0 && v2_ccw[1] > 0) || (vertices_ccw[edges_ccw[j][1]][1] > 0 && v2_ccw[1] > L[1]))
			{
				cout << "entered cww loop 50\n" << "e0 = " << edges_ccw[j][0] << " e1 = " << edges_ccw[j][1] << "\n";
				edges_ccw[j][3]--;
			}
		}
	}

	// Update position of vertices i1 and i2 with pbc
	if (vertices_ccw[i1][0] < 0)
	{
		vertices_ccw[i1][0] = vertices_ccw[i1][0] + L[0];
	}
	else if (vertices_ccw[i1][0] > L[0])
	{
		vertices_ccw[i1][0] = vertices_ccw[i1][0] - L[0];
	}
	if (vertices_ccw[i1][1] < 0)
	{
		vertices_ccw[i1][1] = vertices_ccw[i1][1] + L[1];
	}
	else if (vertices_ccw[i1][1] > L[1])
	{
		vertices_ccw[i1][1] = vertices_ccw[i1][1] - L[1];
	}
	if (vertices_ccw[i2][0] < 0)
	{
		vertices_ccw[i2][0] = vertices_ccw[i2][0] + L[0];
	}
	else if (vertices_ccw[i2][0] > L[0])
	{
		vertices_ccw[i2][0] = vertices_ccw[i2][0] - L[0];
	}
	if (vertices_ccw[i2][1] < 0)
	{
		vertices_ccw[i2][1] = vertices_ccw[i2][1] + L[1];
	}
	else if (vertices_ccw[i2][1] > L[1])
	{
		vertices_ccw[i2][1] = vertices_ccw[i2][1] - L[1];
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_ccw_data(network_ccw, edges_ccw, vertices_ccw);

	return T1_ccw_data;
}

// Perform T1 transition and check the energy change
std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_transition(std::vector<std::vector<double>> vertices, std::vector<Polygon> network,
																												std::vector<std::vector<int>> edges, std::vector<double> L, double lmin,
																												double ka, double Lambda, double gamma, double ksep)
{
	ofstream aux_logfile;
	aux_logfile.open("log_aux.txt", std::ios_base::app);

	std::vector<int> reverse;

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_data;
	// std::vector<Polygon> net_tuple = std::get<0>(T1_data);

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
				std::pair<std::vector<Polygon>, std::vector<std::vector<int>>> T1_0_data = T1_0(network, i1, i2, cell_ids, indices, edges);
				std::vector<Polygon> network_0 = T1_0_data.first;
				std::vector<std::vector<int>> edges_0 = T1_0_data.second;
				double E0 = get_total_energy(vertices, network_0, edges_0, ka, L, Lambda, gamma);
				cout << "E0 = " << E0 << "\n";

				// Get cw T1 transition
				std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_cw_data = T1_cw(network, i1, i2, cell_ids, indices, vertices, edges, L, ksep, lmin);
				std::vector<Polygon> network_cw = std::get<0>(T1_cw_data);
				std::vector<std::vector<int> > edges_cw = std::get<1>(T1_cw_data);
				std::vector<std::vector<double> > vertices_cw = std::get<2>(T1_cw_data);
				double E_cw = get_total_energy(vertices_cw, network_cw, edges_cw, ka, L, Lambda, gamma);
				cout << "Ecw = " << E_cw << "\n";

				// Get ccw T1 transition
				std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_ccw_data = T1_ccw(network, i1, i2, cell_ids, indices, vertices, edges, L, ksep, lmin);
				std::vector<Polygon> network_ccw = std::get<0>(T1_ccw_data);
				std::vector<std::vector<int> > edges_ccw = std::get<1>(T1_ccw_data);
				std::vector<std::vector<double> > vertices_ccw = std::get<2>(T1_ccw_data);
				double E_ccw = get_total_energy(vertices_ccw, network_ccw, edges_ccw, ka, L, Lambda, gamma);
				cout << "Eccw = "  << E_ccw << "\n";

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
					cout << "transition happened here\n";
				}
				else
				{
					T1_data = set_T1_ccw(network, T1_ccw_data, cell_ids, edges, indices, L, lmin, ksep, vertices);
					aux_logfile << "T1" << '\t' << i1 << '\t' << i2 << '\t' << dist << '\n';
					cout << "transition happened here\n";
				}

				network = std::get<0>(T1_data);
				edges = std::get<1>(T1_data);
				vertices = std::get<2>(T1_data); 
			}
		}
	}
	return T1_data;
}

std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> set_tuple_data(std::vector<Polygon> network, std::vector<std::vector<int>> edges,
																												 std::vector<std::vector<double>> vertices)
{
	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_data_tuple(network, edges, vertices);
	return T1_data_tuple;
}

// Set new cell indices, vertex positions, and edges for T1 cw transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> set_T1_cw(std::vector<Polygon> network, std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_cw_data,
																											std::vector<int> cell_ids, std::vector<std::vector<int>> edges, std::vector<int> vertex_indices, std::vector<double> L,
																											double lmin, double ksep, std::vector<std::vector<double>> vertices)
{
	std::vector<Polygon> network_cw = std::get<0>(T1_cw_data);
	std::vector<std::vector<int> > edges_cw = std::get<1>(T1_cw_data);
	std::vector<std::vector<double> > vertices_cw = std::get<2>(T1_cw_data);

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_cw = network_cw[i];
		network[cell_ids[i]].indices = cell_cw.indices;
	}

	// Set new vertices
	vertices = vertices_cw;

	int i1 = vertex_indices[0];
	int i2 = vertex_indices[1];
	int i3 = vertex_indices[2];
	int i5 = vertex_indices[4];

	// Set new edges
	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> edge = edges[i];
		if (edge[0] == i1 && edge[1] == i2)
		{
			edges[i][2] = edges_cw[0][2];
			edges[i][3] = edges_cw[0][3];
		}
		if (edge[0] == i2 && edge[1] == i1)
		{
			edges[i][2] = edges_cw[3][2];
			edges[i][3] = edges_cw[3][3];
		}
		// i1-i3 becomes i2-i3
		if (edge[0] == i1 && edge[1] == i3)
		{
			// cout << "i1-i3 becomes i2-i3\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i2;
			edges[i][2] = edges_cw[1][2];
			edges[i][3] = edges_cw[1][3];
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i2-i5 becomes i1-i5
		if (edge[0] == i2 && edge[1] == i5)
		{
			// cout << "i2-i5 becomes i1-i5\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i1;
			edges[i][2] = edges_cw[4][2];
			edges[i][3] = edges_cw[4][3];
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i3-i1 becomes i3-i2
		if (edge[0] == i3 && edge[1] == i1)
		{
			// cout << "i3-i1 becomes i3-i2\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i2;
			edges[i][2] = edges_cw[6][2];
			edges[i][3] = edges_cw[6][3];
			// cout << "new = " << edges[i][1] << '\n';
		}
		// i5-i2 becomes i5-i1
		if (edge[0] == i5 && edge[1] == i2)
		{
			// cout << "i5-i2 becomes i5-i1\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i1;
			edges[i][2] = edges_cw[8][2];
			edges[i][3] = edges_cw[8][3];
			// cout << "new = " << edges[i][1] << '\n';
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_CW_data(network, edges, vertices);
	return T1_CW_data;
}

// Set new cell indices and edges for T1 ccw transition
std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> set_T1_ccw(std::vector<Polygon> network, std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_ccw_data, std::vector<int> cell_ids,
																											 std::vector<std::vector<int>> edges, std::vector<int> vertex_indices, std::vector<double> L,
																											 double lmin, double ksep, std::vector<std::vector<double>> vertices)
{
	std::vector<Polygon> network_ccw = std::get<0>(T1_ccw_data);
	std::vector<std::vector<int> > edges_ccw = std::get<1>(T1_ccw_data);
	std::vector<std::vector<double> > vertices_ccw = std::get<2>(T1_ccw_data);

	// Set new cell indices
	for (int i = 0; i < cell_ids.size(); i++)
	{
		Polygon cell_ccw = network_ccw[i];
		network[cell_ids[i]].indices = cell_ccw.indices;
	}

	// Set new vertices
	vertices = vertices_ccw;

	// Set new edges
	int i1 = vertex_indices[0];
	int i2 = vertex_indices[1];
	int i4 = vertex_indices[3];
	int i6 = vertex_indices[5];

	for (int i = 0; i < edges.size(); i++)
	{
		// cout << "Set new edges entered\n";
		std::vector<int> edge = edges[i];
		if (edge[0] == i1 && edge[1] == i2)
		{
			edges[i][2] = edges_ccw[0][2];
			edges[i][3] = edges_ccw[0][3];
		}
		if (edge[0] == i2 && edge[1] == i1)
		{
			edges[i][2] = edges_ccw[3][2];
			edges[i][3] = edges_ccw[3][3];
		}
		// i1-i4 becomes i2-i4
		if (edge[0] == i1 && edge[1] == i4)
		{
			// cout << "i1-i4 becomes i2-i4\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i2;
			edges[i][2] = edges_ccw[5][2];
			edges[i][3] = edges_ccw[5][3];
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i2-i6 becomes i1-i6
		if (edge[0] == i2 && edge[1] == i6)
		{
			// cout << "i2-i6 becomes i1-i6\n";
			// cout << "original = " << edges[i][0] << '\n';
			edges[i][0] = i1;
			edges[i][2] = edges_ccw[2][2];
			edges[i][3] = edges_ccw[2][3];
			// cout << "new = " << edges[i][0] << '\n';
		}
		// i4-i1 becomes i4-i2
		if (edge[0] == i4 && edge[1] == i1)
		{
			// cout << "i4-i1 becomes i4-i2\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i2;
			edges[i][2] = edges_ccw[7][2];
			edges[i][3] = edges_ccw[7][3];
			// cout << "new = " << edges[i][1] << '\n';
		}
		// i6-i2 becomes i6-i1
		if (edge[0] == i6 && edge[1] == i2)
		{
			// cout << "i6-i2 becomes i6-i1\n";
			// cout << "original = " << edges[i][1] << '\n';
			edges[i][1] = i1;
			edges[i][2] = edges_ccw[9][2];
			edges[i][3] = edges_ccw[9][3];
			// cout << "new = " << edges[i][1] << '\n';
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_CCW_data(network, edges, vertices);
	return T1_CCW_data;
}