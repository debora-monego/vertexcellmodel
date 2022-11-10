#include "parser.hpp"
#include "polygon.hpp"
#include "geometry.hpp"

using namespace std;

// Read vertex indices for each polygon in the network from file
std::vector<std::vector<int> > read_cell_indices(string networkfile)
{

	ifstream myfile;
	myfile.open(networkfile);
	std::string polygon;

	std::vector<std::vector<int> > vertex_indices;
	std::vector<int> cell_indices;

	while (std::getline(myfile, polygon))
	{
		cell_indices.clear();

		string delimiter = "\t";
		int start = 0;
		int end = polygon.find(delimiter);

		while (end != -1)
		{
			cell_indices.push_back(stoi(polygon.substr(start, end - start)));
			start = end + delimiter.size();
			end = polygon.find(delimiter, start);
		}
		cell_indices.push_back(stoi(polygon.substr(start, end - start)));

		vertex_indices.push_back(cell_indices);
	}
	return vertex_indices;
}

// Build network from cell indices
std::vector<Polygon> build_network(std::vector<std::vector<int> > vertex_indices, double A0, double P0, double J)
{
	std::vector<Polygon> network;
	for (int i = 0; i < vertex_indices.size(); i++)
	{

		double pi = atan(1) * 4;
		std::vector<int> indices = vertex_indices[i];
		double theta = random_angle(-pi, pi);
		Polygon cell;
		cell.initialize(i, indices, A0, P0, J, theta);
		network.push_back(cell);
	}
	return network;
}

// Read vertex coordinates
std::vector<std::vector<double> > read_vertices(string vertexfile)
{
	ifstream myfile;
	myfile.open(vertexfile);

	std::vector<std::vector<double> > vertices;

	if (myfile)
	{
		std::string line;

		while (std::getline(myfile, line))
		{

			vertices.push_back(std::vector<double>());

			// Break down the row into column values
			std::stringstream split(line);
			double value;

			while (split >> value)
				vertices.back().push_back(value);
		}
	}
	return vertices;
}

// Read indices for edges between vertices
std::vector<std::vector<int> > read_edges(string edgefile)
{
	ifstream myfile;
	myfile.open(edgefile);

	std::vector<std::vector<int> > edges;

	if (myfile)
	{
		std::string line;

		while (std::getline(myfile, line))
		{

			edges.push_back(std::vector<int>());

			// Break down the row into column values
			std::stringstream split(line);
			int value;

			while (split >> value)
				edges.back().push_back(value);
		}
	}
	return edges;
}

// Write matrix to file
template <class myMatrix>
void write_matrix(std::vector<std::vector<myMatrix> > matrix, string outfile)
{
	ofstream outputfile;
	outputfile.open(outfile);

	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix[0].size(); j++)
		{
			outputfile << matrix[i][j] << ' ';
		}
		outputfile << '\n';
	}

	outputfile.close();
}

// Write vector to file
template <class myVector>
void write_vector(std::vector<myVector> vector, string outfile)
{
	ofstream outputfile;
	outputfile.open(outfile);

	for (int i = 0; i < vector.size(); i++)
	{
		outputfile << vector[i] << '\n';
	}
	outputfile.close();
}
