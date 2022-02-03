import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, pi, sin, cos, acos, atan2, floor

def periodic_diff(v1, v2, L):
	return ((v1 - v2 + L/2.) % L) - L/2.

lx = 9 * (2 / (3 * (3**0.5)))**0.5
ly = 4 * (2 / (3**0.5))**0.5
L = np.array((lx,ly))

with open('../data2/cell_indices.txt') as f:
#with open('network_cells.txt') as f:
    cell_indices = [[int(a) for a in line.split()] for line in f]
    # print(cell_indices)

with open('../data2/vertices.txt') as f:
#with open('vertex_coord.txt') as f:
    # x, y = [float(a) for a in next(f).split()]
    vertices = [[float(a) for a in line.split()] for line in f]
    vertex_indices = {v: k for v, k in enumerate(vertices)}
    # print(vertex_indices)

with open('../data2/edges.txt') as f:
#with open('edges_index.txt') as f:
    edges = [[int(a) for a in line.split()] for line in f]
    # print(edges)

plt.cla()
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for x,y in vertices:
    ax.scatter(x, y, c="m", marker=".", s=50)

for i in cell_indices:
    indices = i
    for j, index in enumerate(indices):
        x1, y1 = vertices[index]
        if j == len(indices) - 1:
            x2, y2 = vertices[indices[0]]
        else:
            x2, y2 = vertices[indices[j+1]]
        v1 = np.array((x1,y1))
        v2 = np.array((x2,y2))
        v2 = v1 + periodic_diff(v2, v1, L)
        x2, y2 = v2
        ax.plot([x1,x2], [y1,y2], c="black")
        
        v2 = np.array((x2,y2))
        v1 = v2 + periodic_diff(v1, v2, L)
        x1, y1 = v1
        ax.plot([x1,x2], [y1,y2], c="black")


G = nx.Graph()

#G.add_edges_from(edges)

options = {
    "font_size": 11,
    "node_size": 1,
    "node_color": "white",
    "edgecolors": "black",
}

#nx.draw_networkx(G, vertex_indices, **options)

# Set margins for the axes so that nodes aren't clipped
ax = plt.gca()
ax.margins(0.20)
plt.axis("off")
plt.show()
