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

		if ((it1 != cell.indices.end()) && (it2 != cell.indices.end()))
		{
			position1 = std::distance(cell.indices.begin(), it1);
			position2 = std::distance(cell.indices.begin(), it2);

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

// /*
// Brings the two vertices that are close enough to go through a T1 transition to the same position,
// i.e., they are now one vertice connected to the other 4 vertices of interest
// data that will be updated: vertices, edges, polygons
std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_intermediate(std::vector<Polygon> network, int i1, int i2,
																												  std::vector<std::vector<double> > vertices, std::vector<int> vertex_indices,
																												  std::vector<std::vector<int> > edges, std::vector<double> L, int T1_type)
{
	// Define polygons
	Polygon cell_0 = network[0];
	Polygon cell_1 = network[1];
	Polygon cell_2 = network[2];
	Polygon cell_3 = network[3];

	if (T1_type == 0)
	{ // CW
	// Cell 0: remove i2
	int position = 0;
	std::vector<int>::iterator it;
	it = std::find(cell_0.indices.begin(), cell_0.indices.end(), i2);
	if (it != cell_0.indices.end())
	{
		position = std::distance(cell_0.indices.begin(), it);
	}
	cell_0.indices.erase(cell_0.indices.begin() + position);
	network[0].indices = cell_0.indices;

	// Cell 2: remove i1
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i1);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}
	cell_2.indices.erase(cell_2.indices.begin() + position);
	network[2].indices = cell_2.indices;
	}
	else if (T1_type == 1)
	{ // CCW
	// Cell 0: remove i1
	int position = 0;
	std::vector<int>::iterator it;
	it = std::find(cell_0.indices.begin(), cell_0.indices.end(), i1);
	if (it != cell_0.indices.end())
	{
		position = std::distance(cell_0.indices.begin(), it);
	}
	cell_0.indices.erase(cell_0.indices.begin() + position);
	network[0].indices = cell_0.indices;

	// Cell 2: remove i2
	it = std::find(cell_2.indices.begin(), cell_2.indices.end(), i2);
	if (it != cell_2.indices.end())
	{
		position = std::distance(cell_2.indices.begin(), it);
	}
	cell_2.indices.erase(cell_2.indices.begin() + position);
	network[2].indices = cell_2.indices;
	}

	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Vertices
	// Initialize copy of vertices
	std::vector<std::vector<double>> vertices_int(vertices.size(), std::vector<double>(2));
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices_int[i] = vertices[i];
	}

	std::vector<double> v1 = vertices[i1];
	std::vector<double> vertex2 = vertices[i2];
	std::vector<int> q{0, 0};
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i1 && edges[i][1] == i2)
		{
			std::vector<int>::const_iterator first = edges[i].begin() + 2;
			std::vector<int>::const_iterator last = edges[i].begin() + 4;
			std::vector<int> pbc(first, last);
			q[0] = q[0] + pbc[0];
			q[1] = q[1] + pbc[1];
		}
	}
	std::vector<double> v2 = add_vectors(v1, pbc_diff(vertex2, v1, L, q));

	// Updated positions of i1 and i2 to be equal to the midpoint between them
	for (int i = 0; i < vertices.size(); i++)
	{
		if (i == i1 || i == i2)
		{
			vertices_int[i][0] = (v1[0] + v2[0]) / 2;
			vertices_int[i][1] = (v1[1] + v2[1]) / 2;

			if (vertices_int[i][0] < 0)
			{
				vertices_int[i][0] = vertices_int[i][0] + L[0];
				for (int j = 0; edges.size(); j++)
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
			else if (vertices_int[i][0] > L[0])
			{
				vertices_int[i][0] = vertices_int[i][0] - L[0];
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
			if (vertices_int[i][1] < 0)
			{
				vertices_int[i][1] = vertices_int[i][1] + L[1];
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
			else if (vertices_int[i][1] > L[1])
			{
				vertices_int[i][1] = vertices_int[i][1] - L[1];
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

	// Edges
	// Initialize matrix of edges
	std::vector<std::vector<int>> edges_int(18, std::vector<int>(4, 0));
	// Update edges
	// Edge 0 : i1 - i2
	edges_int[0][0] = i1;
	edges_int[0][1] = i2;
	edges_int[0][2] = 0;
	edges_int[0][3] = 0;

	// Edge 1 : i1 - i3
	edges_int[1][0] = i1;
	edges_int[1][1] = i3;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i1 && edges[i][1] == i3)
		{
			edges_int[1][2] = edges[i][2];
			edges_int[1][3] = edges[i][3];
		}
	}

	// Edge 2 : i1 - i4
	edges_int[2][0] = i1;
	edges_int[2][1] = i4;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i1 && edges[i][1] == i4)
		{
			edges_int[2][2] = edges[i][2];
			edges_int[2][3] = edges[i][3];
		}
	}

	// Edge 3 : i2 - i1
	edges_int[3][0] = i2;
	edges_int[3][1] = i1;
	edges_int[3][2] = 0;
	edges_int[3][3] = 0;

	// Edge 4 : i2 - i5
	edges_int[4][0] = i2;
	edges_int[4][1] = i5;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i2 && edges[i][1] == i5)
		{
			edges_int[4][2] = edges[i][2];
			edges_int[4][3] = edges[i][3];
		}
	}

	// Edge 5 : i2 - i6
	edges_int[5][0] = i2;
	edges_int[5][1] = i6;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i2 && edges[i][1] == i6)
		{
			edges_int[5][2] = edges[i][2];
			edges_int[5][3] = edges[i][3];
		}
	}

	// Edge 6 : i3 - i1
	edges_int[6][0] = i3;
	edges_int[6][1] = i1;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i3 && edges[i][1] == i1)
		{
			edges_int[6][2] = edges[i][2];
			edges_int[6][3] = edges[i][3];
		}
	}

	// Edge 7 : i4 - i1
	edges_int[7][0] = i4;
	edges_int[7][1] = i1;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i4 && edges[i][1] == i1)
		{
			edges_int[7][2] = edges[i][2];
			edges_int[7][3] = edges[i][3];
		}
	}

	// Edge 8 : i5 - i2
	edges_int[8][0] = i5;
	edges_int[8][1] = i2;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i5 && edges[i][1] == i2)
		{
			edges_int[8][2] = edges[i][2];
			edges_int[8][3] = edges[i][3];
		}
	}

	// Edge 9 : i6 - i2
	edges_int[9][0] = i6;
	edges_int[9][1] = i2;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i6 && edges[i][1] == i2)
		{
			edges_int[9][2] = edges[i][2];
			edges_int[9][3] = edges[i][3];
		}
	}

	// Edge 10 : i1 - i5
	edges_int[10][0] = i1;
	edges_int[10][1] = i5;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i2 && edges[i][1] == i5)
		{
			edges_int[10][2] = edges[i][2];
			edges_int[10][3] = edges[i][3];
		}
	}	

	// Edge 11 : i5 - i1
	edges_int[10][0] = i5;
	edges_int[10][1] = i1;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i5 && edges[i][1] == i2)
		{
			edges_int[11][2] = edges[i][2];
			edges_int[11][3] = edges[i][3];
		}
	}

	// Edge 12 : i1 - i6
	edges_int[12][0] = i1;
	edges_int[12][1] = i6;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i2 && edges[i][1] == i6)
		{
			edges_int[12][2] = edges[i][2];
			edges_int[12][3] = edges[i][3];
		}
	}	
	
	// Edge 13: i6 - i1
	edges_int[13][0] = i6;
	edges_int[13][1] = i1;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i6 && edges[i][1] == i2)
		{
			edges_int[13][2] = edges[i][2];
			edges_int[13][3] = edges[i][3];
		}
	}

	// Edge 14: i2 - i3
	edges_int[14][0] = i2;
	edges_int[14][1] = i3;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i1 && edges[i][1] == i3)
		{
			edges_int[14][2] = edges[i][2];
			edges_int[14][3] = edges[i][3];
		}
	}

	// Edge 15: i3 - i2
	edges_int[15][0] = i3;
	edges_int[15][1] = i2;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i3 && edges[i][1] == i1)
		{
			edges_int[15][2] = edges[i][2];
			edges_int[15][3] = edges[i][3];
		}
	}

	// Edge 16: i2 - i4
	edges_int[16][0] = i2;
	edges_int[16][1] = i4;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i1 && edges[i][1] == i4)
		{
			edges_int[16][2] = edges[i][2];
			edges_int[16][3] = edges[i][3];
		}
	}

	// Edge 17: i4 - i2
	edges_int[17][0] = i4;
	edges_int[17][1] = i2;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i][0] == i4 && edges[i][1] == i1)
		{
			edges_int[17][2] = edges[i][2];
			edges_int[17][3] = edges[i][3];
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_int_data(network, edges_int, vertices_int);

	return T1_int_data;
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
	int T1_type = 0;

	// Create intermediate state in which i1=i2
	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_int = T1_intermediate(network, i1, i2, vertices, vertex_indices, edges, L, T1_type);
	std::vector<Polygon> network_int = std::get<0>(T1_int);
	std::vector<std::vector<int> > edges_int= std::get<1>(T1_int);
	std::vector<std::vector<double> > vertices_int = std::get<2>(T1_int);

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
	edges_cw[1][2] = edges_int[14][2];
	edges_cw[1][3] = edges_int[14][3];

	// Edge 2 : i1 - i4
	edges_cw[2][0] = i1;
	edges_cw[2][1] = i4;
	edges_cw[2][2] = edges_int[2][2];
	edges_cw[2][3] = edges_int[2][3];

	// Edge 3 : i2 - i1
	edges_cw[3][0] = i2;
	edges_cw[3][1] = i1;

	// Edge 4 : i1 - i5
	edges_cw[4][0] = i1;
	edges_cw[4][1] = i5;
	edges_cw[4][2] = edges_int[10][2];
	edges_cw[4][3] = edges_int[10][3];

	// Edge 5 : i2 - i6
	edges_cw[5][0] = i2;
	edges_cw[5][1] = i6;
	edges_cw[5][2] = edges_int[5][2];
	edges_cw[5][3] = edges_int[5][3];

	// Edge 6 : i3 - i2
	edges_cw[6][0] = i3;
	edges_cw[6][1] = i2;
	edges_cw[6][2] = edges_int[15][2];
	edges_cw[6][3] = edges_int[15][3];

	// Edge 7 : i4 - i1
	edges_cw[7][0] = i4;
	edges_cw[7][1] = i1;
	edges_cw[7][2] = edges_int[7][2];
	edges_cw[7][3] = edges_int[7][3];

	// Edge 8 : i5 - i1
	edges_cw[8][0] = i5;
	edges_cw[8][1] = i1;
	edges_cw[8][2] = edges_int[11][2];
	edges_cw[8][3] = edges_int[11][3];

	// Edge 9 : i6 - i2
	edges_cw[9][0] = i6;
	edges_cw[9][1] = i2;
	edges_cw[9][2] = edges_int[9][2];
	edges_cw[9][3] = edges_int[9][3];

	// Vertices
	// Initialize copy of vertices
	std::vector<std::vector<double> > vertices_cw(vertices.size(), std::vector<double>(2));
	for (int i = 0; i < vertices.size(); i++)
	{
		vertices_cw[i] = vertices[i];
	}

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

	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices_cw[i1][0] = move_transition[0];
	vertices_cw[i1][1] = move_transition[1];
	vertices_cw[i2][0] = move_transition[2];
	vertices_cw[i2][1] = move_transition[3];

	// Wrap around periodic boundaries and update periodicity vector q in edges
	for (int j = 0; j < edges_cw.size(); j++)
	{
		if ((edges_cw[j][0] == i1 && edges_cw[j][1] == i2) || (edges_cw[j][0] == i2 && edges_cw[j][1] == i1)) 
		{
			if ((vertices_cw[edges_cw[j][0]][0] > 0 && vertices_cw[edges_cw[j][1]][0] < 0) || (vertices_cw[edges_cw[j][0]][0] > L[0] && vertices_cw[edges_cw[j][1]][0] > 0))
			{
				edges_cw[j][2]--;
			}
			else if ((vertices_cw[edges_cw[j][0]][0] < 0 && vertices_cw[edges_cw[j][1]][0] > 0) || (vertices_cw[edges_cw[j][0]][0] > 0 && vertices_cw[edges_cw[j][1]][0] > L[0]))
			{
				edges_cw[j][2]++;
			}
			if ((vertices_cw[edges_cw[j][0]][1] > 0 && vertices_cw[edges_cw[j][1]][1] < 0) || (vertices_cw[edges_cw[j][0]][1] > L[1] && vertices_cw[edges_cw[j][1]][1] > 0))
			{
				edges_cw[j][3]--;
			}
			else if ((vertices_cw[edges_cw[j][0]][1] < 0 && vertices_cw[edges_cw[j][1]][1] > 0) || (vertices_cw[edges_cw[j][0]][1] > 0 && vertices_cw[edges_cw[j][1]][1] > L[1]))
			{
				edges_cw[j][3]++;
			}
		} else if (edges_cw[j][0] == i1 || edges_cw[j][0] == i2)
		{
			std::vector<double> vertex2_int = vertices_int[edges_cw[j][0]];
			std::vector<double> vertex2_cw = vertices_cw[edges_cw[j][0]];
			if ((vertex2_int[0] > 0 && vertex2_cw[0] < 0) || (vertex2_int[0] > L[0] && vertex2_cw[0] > 0))
			{
				edges_cw[j][2]++;
			}
			else if ((vertex2_int[0] < 0 && vertex2_cw[0] > 0) || (vertex2_int[0] > 0 && vertex2_cw[0] > L[0]))
			{
				edges_cw[j][2]--;
			}
			if ((vertex2_int[1] > 0 && vertex2_cw[1] < 0) && (vertex2_int[1] > L[1] && vertex2_cw[1] > 0))
			{
				edges_cw[j][3]++;
			}
			else if ((vertex2_int[1] < 0 && vertex2_cw[1] > 0) || (vertex2_int[1] > 0 && vertex2_cw[1] > L[1]))
			{
				edges_cw[j][3]--;
			}
		} else if (edges_cw[j][1] == i1 || edges_cw[j][1] == i2)
		{
			std::vector<double> vertex1_int = vertices_int[edges_cw[j][1]];
			std::vector<double> vertex1_cw = vertices_cw[edges_cw[j][1]];
			if ((vertex1_int[0] > 0 && vertex1_cw[0] < 0) || (vertex1_int[0] > L[0] && vertex1_cw[0] > 0))
			{
				edges_cw[j][2]--;
			}
			else if ((vertex1_int[0] < 0 && vertex1_cw[0] > 0) || (vertex1_int[0] > 0 && vertex1_cw[0] > L[0]))
			{
				edges_cw[j][2]++;
			}
			if ((vertex1_int[1] > 0 && vertex1_cw[1] < 0) && (vertex1_int[1] > L[1] && vertex1_cw[1] > 0))
			{
				edges_cw[j][3]--;
			}
			else if ((vertex1_int[1] < 0 && vertex1_cw[1] > 0) || (vertex1_int[1] > 0 && vertex1_cw[1] > L[1]))
			{
				edges_cw[j][3]++;
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
	int T1_type = 1;

	// Create intermediate state in which i1=i2
	std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_int = T1_intermediate(network, i1, i2, vertices, vertex_indices, edges, L, T1_type);

	std::vector<Polygon> network_int = std::get<0>(T1_int);
	std::vector<std::vector<int> > edges_int= std::get<1>(T1_int);
	std::vector<std::vector<double> > vertices_int = std::get<2>(T1_int);

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
	edges_ccw[1][2] = edges_int[1][2];
	edges_ccw[1][3] = edges_int[1][3];

	// Edge 2 : i1 - i6
	edges_ccw[2][0] = i1;
	edges_ccw[2][1] = i6;
	edges_ccw[2][2] = edges_int[12][2];
	edges_ccw[2][3] = edges_int[12][3];

	// Edge 3 : i2 - i1
	edges_ccw[3][0] = i2;
	edges_ccw[3][1] = i1;

	// Edge 4 : i2 - i5
	edges_ccw[4][0] = i2;
	edges_ccw[4][1] = i5;
	edges_ccw[4][2] = edges_int[4][2];
	edges_ccw[4][3] = edges_int[4][3];

	// Edge 5 : i2 - i4
	edges_ccw[5][0] = i2;
	edges_ccw[5][1] = i4;
	edges_ccw[5][2] = edges_int[16][2];
	edges_ccw[5][3] = edges_int[16][3];

	// Edge 6 : i3 - i1
	edges_ccw[6][0] = i3;
	edges_ccw[6][1] = i1;
	edges_ccw[6][2] = edges_int[6][2];
	edges_ccw[6][3] = edges_int[6][3];

	// Edge 7 : i4 - i2
	edges_ccw[7][0] = i4;
	edges_ccw[7][1] = i2;
	edges_ccw[7][2] = edges_int[17][2];
	edges_ccw[7][3] = edges_int[17][3];

	// Edge 8 : i5 - i2
	edges_ccw[8][0] = i5;
	edges_ccw[8][1] = i2;
	edges_ccw[8][2] = edges_int[8][2];
	edges_ccw[8][3] = edges_int[8][3];

	// Edge 9 : i6 - i1
	edges_ccw[9][0] = i6;
	edges_ccw[9][1] = i1;
	edges_ccw[9][2] = edges_int[13][2];
	edges_ccw[9][3] = edges_int[13][3];

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

	std::vector<double> move_transition = set_separation_transition(v1[0], v1[1], v2[0], v2[1], ksep, lmin, angle);
	vertices_ccw[i1][0] = move_transition[0];
	vertices_ccw[i1][1] = move_transition[1];
	vertices_ccw[i2][0] = move_transition[2];
	vertices_ccw[i2][1] = move_transition[3];

	// Wrap around periodic boundaries and update periodicity vector q in edges 
	for (int j = 0; j < edges_ccw.size(); j++)
	{
		if ((edges_ccw[j][0] == i1 && edges_ccw[j][1] == i2) || (edges_ccw[j][0] == i2 && edges_ccw[j][1] == i1))
		{
			if ((vertices_ccw[edges_ccw[j][0]][0] > 0 && vertices_ccw[edges_ccw[j][1]][0] < 0) || (vertices_ccw[edges_ccw[j][0]][0] > L[0] && vertices_ccw[edges_ccw[j][1]][0] > 0))
			{
				edges_ccw[j][2]--;
			}
			else if ((vertices_ccw[edges_ccw[j][0]][0] < 0 && vertices_ccw[edges_ccw[j][1]][0] > 0) || (vertices_ccw[edges_ccw[j][0]][0] > 0 && vertices_ccw[edges_ccw[j][1]][0] > L[0]))
			{
				edges_ccw[j][2]++;
			}
			if ((vertices_ccw[edges_ccw[j][0]][1] > 0 && vertices_ccw[edges_ccw[j][1]][1] < 0) || (vertices_ccw[edges_ccw[j][0]][1] > L[1] && vertices_ccw[edges_ccw[j][1]][1] > 0))
			{
				edges_ccw[j][3]--;
			}
			else if ((vertices_ccw[edges_ccw[j][0]][1] < 0 && vertices_ccw[edges_ccw[j][1]][1] > 0) || (vertices_ccw[edges_ccw[j][0]][1] > 0 && vertices_ccw[edges_ccw[j][1]][1] > L[1]))
			{
				edges_ccw[j][3]++;
			}
		} else if (edges_ccw[j][0] == i1 || edges_ccw[j][0] == i2)
		{
			std::vector<double> vertex2_int = vertices_int[edges_ccw[j][0]];
			std::vector<double> vertex2_ccw = vertices_ccw[edges_ccw[j][0]];
			if ((vertex2_int[0] > 0 && vertex2_ccw[0] < 0) || (vertex2_int[0] > L[0] && vertex2_ccw[0] > 0))
			{
				edges_ccw[j][2]++;
			}
			else if ((vertex2_int[0] < 0 && vertex2_ccw[0] > 0) || (vertex2_int[0] > 0 && vertex2_ccw[0] > L[0]))
			{
				edges_ccw[j][2]--;
			}
			if ((vertex2_int[1] > 0 && vertex2_ccw[1] < 0) && (vertex2_int[1] > L[1] && vertex2_ccw[1] > 0))
			{
				edges_ccw[j][3]++;
			}
			else if ((vertex2_int[1] < 0 && vertex2_ccw[1] > 0) || (vertex2_int[1] > 0 && vertex2_ccw[1] > L[1]))
			{
				edges_ccw[j][3]--;
			}
		} else if (edges_ccw[j][1] == i1 || edges_ccw[j][1] == i2)
		{
			std::vector<double> vertex1_int = vertices_int[edges_ccw[j][1]];
			std::vector<double> vertex1_ccw = vertices_ccw[edges_ccw[j][1]];
			if ((vertex1_int[0] > 0 && vertex1_ccw[0] < 0) || (vertex1_int[0] > L[0] && vertex1_ccw[0] > 0))
			{
				edges_ccw[j][2]--;
			}
			else if ((vertex1_int[0] < 0 && vertex1_ccw[0] > 0) || (vertex1_int[0] > 0 && vertex1_ccw[0] > L[0]))
			{
				edges_ccw[j][2]++;
			}
			if ((vertex1_int[1] > 0 && vertex1_ccw[1] < 0) && (vertex1_int[1] > L[1] && vertex1_ccw[1] > 0))
			{
				edges_ccw[j][3]--;
			}
			else if ((vertex1_int[1] < 0 && vertex1_ccw[1] > 0) || (vertex1_int[1] > 0 && vertex1_ccw[1] > L[1]))
			{
				edges_ccw[j][3]++;
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
																												double ka, double Lambda, double gamma, double ksep, double lambda_potts)
{
	ofstream aux_logfile;
	aux_logfile.open("log_aux.txt", std::ios_base::app);

	std::vector<int> reverse;

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_data;
	T1_data = set_tuple_data(network, edges, vertices);

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
			}
			else
			{
				//  Find minimum configuration
				reverse.insert(reverse.end(), {i2, i1});

				// Find 6 indices for vertices involved in the transition
				std::vector<int> indices = get_vertex_indices(network, i1, i2, cell_ids);

				// Get original configuration
				std::pair<std::vector<Polygon>, std::vector<std::vector<int>>> T1_0_data = T1_0(network, i1, i2, cell_ids, indices, edges);
				std::vector<Polygon> network_0 = T1_0_data.first;
				std::vector<std::vector<int>> edges_0 = T1_0_data.second;
				double E0 = get_total_energy(vertices, network_0, edges_0, ka, L, Lambda, gamma, lambda_potts);

				// Get cw T1 transition
				std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_cw_data = T1_cw(network, i1, i2, cell_ids, indices, vertices, edges, L, ksep, lmin);
				std::vector<Polygon> network_cw = std::get<0>(T1_cw_data);
				std::vector<std::vector<int> > edges_cw = std::get<1>(T1_cw_data);
				std::vector<std::vector<double> > vertices_cw = std::get<2>(T1_cw_data);
				double E_cw = get_total_energy(vertices_cw, network_cw, edges_cw, ka, L, Lambda, gamma, lambda_potts);

				// Get ccw T1 transition
				std::tuple<std::vector<Polygon>, std::vector<std::vector<int> >, std::vector<std::vector<double> > > T1_ccw_data = T1_ccw(network, i1, i2, cell_ids, indices, vertices, edges, L, ksep, lmin);
				std::vector<Polygon> network_ccw = std::get<0>(T1_ccw_data);
				std::vector<std::vector<int> > edges_ccw = std::get<1>(T1_ccw_data);
				std::vector<std::vector<double> > vertices_ccw = std::get<2>(T1_ccw_data);
				double E_ccw = get_total_energy(vertices_ccw, network_ccw, edges_ccw, ka, L, Lambda, gamma, lambda_potts);

				std::vector<double> all_energies{E0, E_cw, E_ccw};
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
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	// Set new edges
	for (int i = 0; i < edges.size(); i++)
	{	
		std::vector<int> edge = edges[i];
		// i1-i2 and i2-i1 stay the same
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
		// i1-i4 and i4-i1 stay the same
		if (edge[0] == i1 && edge[1] == i4)
		{
			edges[i][2] = edges_cw[2][2];
			edges[i][3] = edges_cw[2][3];
		}
		if (edge[0] == i4 && edge[1] == i1)
		{
			edges[i][2] = edges_cw[7][2];
			edges[i][3] = edges_cw[7][3];
		}
		// i2-i6 and i6-i2 stay the same
		if (edge[0] == i2 && edge[1] == i6)
		{
			edges[i][2] = edges_cw[5][2];
			edges[i][3] = edges_cw[5][3];
		}
		if (edge[0] == i6 && edge[1] == i2)
		{
			edges[i][2] = edges_cw[9][2];
			edges[i][3] = edges_cw[9][3];	
		}
		// i1-i3 becomes i2-i3
		if (edge[0] == i1 && edge[1] == i3)
		{
			edges[i][0] = i2;
			edges[i][2] = edges_cw[1][2];
			edges[i][3] = edges_cw[1][3];
		}
		// i2-i5 becomes i1-i5
		if (edge[0] == i2 && edge[1] == i5)
		{
			edges[i][0] = i1;
			edges[i][2] = edges_cw[4][2];
			edges[i][3] = edges_cw[4][3];
		}
		// i3-i1 becomes i3-i2
		if (edge[0] == i3 && edge[1] == i1)
		{
			edges[i][1] = i2;
			edges[i][2] = edges_cw[6][2];
			edges[i][3] = edges_cw[6][3];
		}
		// i5-i2 becomes i5-i1
		if (edge[0] == i5 && edge[1] == i2)
		{
			edges[i][1] = i1;
			edges[i][2] = edges_cw[8][2];
			edges[i][3] = edges_cw[8][3];
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
	int i3 = vertex_indices[2];
	int i4 = vertex_indices[3];
	int i5 = vertex_indices[4];
	int i6 = vertex_indices[5];

	for (int i = 0; i < edges.size(); i++)
	{
		std::vector<int> edge = edges[i];
		// i1 -i2 and i2-i1 stay the same
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
		// i1 -i3 and i3-i1 stay the same
		if (edge[0] == i1 && edge[1] == i3)
		{
			edges[i][2] = edges_ccw[1][2];
			edges[i][3] = edges_ccw[1][3];
		}
		if (edge[0] == i3 && edge[1] == i1)
		{
			edges[i][2] = edges_ccw[6][2];
			edges[i][3] = edges_ccw[6][3];
		}
		// i2 -i5 and i5-i2 stay the same
		if (edge[0] == i2 && edge[1] == i5)
		{
			edges[i][2] = edges_ccw[4][2];
			edges[i][3] = edges_ccw[4][3];
		}
		if (edge[0] == i5 && edge[1] == i2)
		{
			edges[i][2] = edges_ccw[8][2];
			edges[i][3] = edges_ccw[8][3];
		}
		// i1-i4 becomes i2-i4
		if (edge[0] == i1 && edge[1] == i4)
		{
			edges[i][0] = i2;
			edges[i][2] = edges_ccw[5][2];
			edges[i][3] = edges_ccw[5][3];
		}
		// i2-i6 becomes i1-i6
		if (edge[0] == i2 && edge[1] == i6)
		{
			edges[i][0] = i1;
			edges[i][2] = edges_ccw[2][2];
			edges[i][3] = edges_ccw[2][3];
		}
		// i4-i1 becomes i4-i2
		if (edge[0] == i4 && edge[1] == i1)
		{
			edges[i][1] = i2;
			edges[i][2] = edges_ccw[7][2];
			edges[i][3] = edges_ccw[7][3];
		}
		// i6-i2 becomes i6-i1
		if (edge[0] == i6 && edge[1] == i2)
		{
			edges[i][1] = i1;
			edges[i][2] = edges_ccw[9][2];
			edges[i][3] = edges_ccw[9][3];
		}
	}

	std::tuple<std::vector<Polygon>, std::vector<std::vector<int>>, std::vector<std::vector<double>>> T1_CCW_data(network, edges, vertices);
	
	return T1_CCW_data;
}