#!/usr/bin/python
import numpy as np 
import matplotlib.pyplot as plt
import sys
#from periodic_voronoi import L

lx = float(sys.argv[2])
ly = float(sys.argv[3])

L = np.array([lx,ly])

# Difference with respect to periodic boundaries
# def periodic_diff(v1, v2, L):
# 	return ((v1 - v2 + L/2.) % L) - L/2.
def periodic_diff(v1, v2, L, q1):
    return (v1 - v2 + q1 * L)


def periodic_plot(x1,y1,x2,y2,color,L,q1):

	v1 = np.array((x1,y1))
	v2 = np.array((x2,y2))
	v2 = v1 + periodic_diff(v2, v1, L, q1)
	x2,y2 = v2
	plt.plot([x1,x2], [y1,y2], c=color)

	# v2 = np.array((x2,y2))
	# v1 = v2 + periodic_diff(v1, v2, L)
	# x1,y1 = v1
	# plt.plot([x1,x2], [y1,y2], c=color)

	return 


def plot_points(vertices):
	for x,y in vertices:
		plt.scatter(x, y, color="r")
	return


def plot_edges(vertices, edges, L):
	for i1,i2,i3,i4 in edges:
		if i1 != -1 and i2 != -1:
			x1,y1 = vertices[i1]
			x2,y2 = vertices[i2]
			q1 = np.array(i3,i4)
			periodic_plot(x1,y1,x2,y2,"k",L,q1)
	return



def plot_polygons(vertices, polys, edges, L):

	for poly in polys:
		for i,index in enumerate(poly):

			if index == poly[-1]:
				if index != -1 and poly[0] != -1:
					x1,y1 = vertices[index]
					x2,y2 = vertices[poly[0]]
					indexv2 = poly[0]
			else:
				if index != -1 and poly[i+1] != -1:
					x1,y1 = vertices[index]
					x2,y2 = vertices[poly[i+1]]
					indexv2 = poly[i+1]

			for edge in edges:
				i1 = edge[0]
				i2 = edge[1]
				if i1 == index and i2 == indexv2:
					q1 = edge[2:]

			periodic_plot(x1,y1,x2,y2,"k",L,q1)
            
	return

def save_plot(outfile, L):

	# remove tick marks
	frame = plt.gca()
	frame.axes.get_xaxis().set_ticks([])
	frame.axes.get_yaxis().set_ticks([])

	plt.axis([0,L[0],0,L[1]])

	# save and close plot
	plt.savefig(outfile)
	plt.close()
	return



def read_poly(file):
	indices = []
	f = open(file)
	for line in f:
		cell_indices = []
		linesplit = line.strip().split("\t")
		for i in linesplit:
			cell_indices.append(int(i))
		indices.append(cell_indices)
	f.close()
	return indices



# filename for plot
outfile = sys.argv[1]

vertex_file = "periodic_vertices.txt"
edge_file = "periodic_edges.txt"
poly_file = "periodic_polygons.txt"


vertices = np.loadtxt(vertex_file)
edges = np.loadtxt(edge_file, dtype=int)
polys = read_poly(poly_file)


plot_points(vertices)
plot_polygons(vertices, polys, edges, L)
#plot_edges(vertices, edges, L)

save_plot(outfile,L)