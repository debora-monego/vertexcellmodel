#!/usr/bin/python
import numpy as np
from scipy.spatial import Voronoi
import sys
from tile import *
from periodic import *

np.random.seed([34123])

N = int(sys.argv[1])
lx = float(sys.argv[2])
ly = float(sys.argv[3])

L = np.array([lx,ly])

# generate random points
points = np.random.uniform(0, lx, size=(N, 2))

# tile points
point_tile = tile_points(points, N, L)

# generate voronoi diagram
vor = Voronoi(point_tile)

# get vertices
# create an index map from old vertex index to new vertex index
tile_vertices = vor.vertices
vertices, index_map = get_vertices(tile_vertices, L)

# get edges
tile_edges = vor.ridge_vertices
edges = get_edges(tile_vertices, tile_edges, index_map, L)

# get polygons
# new index map finds the points that are just outside the box
# but will correspond to the indices of polygons built in tile
tile_vertices = vor.vertices
index_map = get_new_index_map(tile_vertices, vertices, index_map, L)

# polygons in tiled vertices
regions = vor.regions
polygons = get_polygons(regions, index_map, vertices, L, edges)

# write vertices
np.savetxt("periodic_vertices.txt", vertices)

# write edges
f = open("periodic_edges.txt", "w+")
for edge in edges:
    f.write("%d\t%d\t%d\t%d\n" % (edge[0], edge[1], edge[2], edge[3]))
f.close()

# f = open("periodic_qindex.txt", "w+")
# for pbc in q:
#     f.write("%d\t%d\n" % (pbc[0], pbc[1]))
# f.close()

# write polygons
f = open("periodic_polygons.txt", "w+")
for poly in polygons:
    if len(poly) != 0:
        for index in poly:
            if index == poly[0]:
                f.write("%d" % index)
            else:
                f.write("\t%d" % index)
        f.write("\n")
f.close()

