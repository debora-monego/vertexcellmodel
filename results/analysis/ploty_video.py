#!/usr/bin/python
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import sys
from matplotlib.animation import FuncAnimation

lx = float(sys.argv[2])
ly = float(sys.argv[3])
L = np.array([lx, ly])

# Difference with respect to periodic boundaries


def periodic_diff(v1, v2, L, q1):
    return (v1 - v2 + q1 * L)


def periodic_plot(x1, y1, x2, y2, color,L,q1):
    #L = np.array([lx,ly])
    v1 = np.array((x1, y1))
    v2 = np.array((x2, y2))
    v2 = v1 + periodic_diff(v2, v1, L, q1)
    x2, y2 = v2
    plt.plot([x1, x2], [y1, y2], c=color)

    return


def plot_points(vertices):
    for x, y in vertices:
        plt.scatter(x, y, color="r")
    return


def plot_edges(vertices, edges, L):
    for i1, i2, i3, i4 in edges:
        if i1 != -1 and i2 != -1:
            x1, y1 = vertices[i1]
            x2, y2 = vertices[i2]
            q1 = np.array(i3, i4)
            periodic_plot(x1, y1, x2, y2, "k", L, q1)
    return


def plot_polygons(vertices, polys):

    for poly in polys:
        for i, index in enumerate(poly):

            if index == poly[-1]:
                if index != -1 and poly[0] != -1:
                    x1, y1 = vertices[index]
                    x2, y2 = vertices[poly[0]]
                    indexv2 = poly[0]
            else:
                if index != -1 and poly[i+1] != -1:
                    x1, y1 = vertices[index]
                    x2, y2 = vertices[poly[i+1]]
                    indexv2 = poly[i+1]

            for edge in edges:
                i1 = edge[0]
                i2 = edge[1]
                if i1 == index and i2 == indexv2:
                    q1 = edge[2:]

            periodic_plot(x1, y1, x2, y2, "k", L, q1)
    return


def save_plot(outfile, L):

    # remove tick marks
    frame = plt.gca()
    frame.axes.get_xaxis().set_ticks([])
    frame.axes.get_yaxis().set_ticks([])

    plt.axis([0, L[0], 0, L[1]])

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

# fig = plt.figure()
# ax = plt.axes(xlim=(0, lx), ylim=(0, ly))


vertex_file = "../vertex_coord.txt"
edge_file = "../edges_index.txt"
poly_file = "../network_cells.txt"


vertices = np.loadtxt(vertex_file)
edges = np.loadtxt(edge_file, dtype=int)
polys = read_poly(poly_file)



fig, ax = plt.subplots()
xdata, ydata0, ydata1 = [], [], []
vert, = plt.plot([], [], 'r', animated=True)
pol, = plt.plot([], [], 'b', animated=True)
f = np.linspace(0, 5, 500)

def init():
    ax.set_xlim(0, lx)
    ax.set_ylim(0, ly)
    vert.set_data(xdata,ydata0)
    pol.set_data(xdata,ydata1)
    return vert, pol

def update(frame):
    xdata.append(frame)
    ydata0.append(vertices)
    ydata1.append(polys)
    vert.set_data(xdata, ydata0) #sets may works

    pol.set_data(xdata, ydata1)
    return vert, pol,

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=2.5, repeat=False)
plt.show()

plot_points(vertices)
plot_polygons(vertices, polys)


save_plot(outfile,L)
