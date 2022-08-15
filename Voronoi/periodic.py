#!/usr/bin/python
import numpy as np

# Difference with respect to periodic boundaries


def periodic_diff(v1, v2, L, q1):
    return (v1 - v2 + q1 * L)


# get vertices inside original bounds
# box length is determined by L
def get_vertices(tile_vertices, L):
    vertices = []
    index_map = {}

    for i, (x, y) in enumerate(tile_vertices):
        if x >= 0 and x <= L[0]:
            if y >= 0 and y <= L[1]:
                curr_indx = len(vertices)
                index_map[i] = curr_indx
                vertices.append((x, y))
    return np.array(vertices), index_map


# get edges with respect to periodic boundaries
# vertices and edges are from voronoi with tiled coordinates

def get_edges(vertices, edges, index_map, L):
    e = []
    #q = []
    for edge in edges:
        i1 = edge[0]
        i2 = edge[1]

        # Case 1: i1 and i2 in index map
        # Append indices for vertex list wrt periodic bounds
        if i1 in index_map and i2 in index_map:
            e.append((index_map[i1], index_map[i2], 0, 0))
            # q.append((0,0))
            e.append((index_map[i2], index_map[i1], 0, 0))
            # q.append((0,0))

        # Case 2: neither in index map
        # Do nothing

        # Case 3: i1 or i2 in index map
        # 		  but not both
        # Because of periodicity, this edge needs to have a
        # connecting edge that wraps around
        if i1 in index_map and i2 not in index_map:
            # find the vertex that it "wraps around" to
            # v = vertex in plane
            v = np.array(vertices[i1])
            # v1 = vertex out of plane
            v1 = np.array(vertices[i2])
            qx = 0
            qy = 0

            # print v1
            if v1[0] < 0:
                v1[0] = L[0] + v1[0]
                qx -= 1

            if v1[0] > L[0]:
                v1[0] = v1[0] - L[0]
                qx += 1

            if v1[1] < 0:
                v1[1] = L[1] + v1[1]
                qy -= 1

            if v1[1] > L[1]:
                v1[1] = v1[1] - L[1]
                qy += 1

            # # # find index of this vertex
            for key in index_map:
                if abs(vertices[key][0] - v1[0]) < 10**-6:
                    e.append((index_map[i1], index_map[key], qx, qy))
                    #q.append((qx, qy))

        # include this IF you want duplicate edges...
        if i1 not in index_map and i2 in index_map:
            # find the vertex that it "wraps around" to
            # v = vertex in plane
            v = np.array(vertices[i2])
            # v1 = vertex out of plane
            v1 = np.array(vertices[i1])
            qx = 0
            qy = 0

            # print v1
            if v1[0] < 0:
                v1[0] = L[0] + v1[0]
                qx -= 1

            if v1[0] > L[0]:
                v1[0] = v1[0] - L[0]
                qx += 1

            if v1[1] < 0:
                v1[1] = L[1] + v1[1]
                qy -= 1

            if v1[1] > L[1]:
                v1[1] = v1[1] - L[1]
                qy += 1

            # # find index of this vertex
            for key in index_map:
                if abs(vertices[key][0] - v1[0]) < 10**-6:
                    e.append((index_map[i2], index_map[key], qx, qy))
                    #q.append((qx, qy))

    return e  # , q


# add indices mapping to new vertices outside of bounds
def get_new_index_map(vertices, v, index_map, L):

    for i, (x, y) in enumerate(vertices):
        x1 = x
        y1 = y

        if x < 0 and x > -L[0]:
            x1 = x + L[0]

        if x > L[0] and x < (2 * L[0]):
            x1 = x - L[0]

        if y < 0 and y > -L[1]:
            y1 = y + L[1]

        if y > L[1] and y < (2 * L[1]):
            y1 = y - L[1]

        # look up new x,y in list
        for j, (x2, y2) in enumerate(v):
            if abs(x1 - x2) < 10**-6:
                if abs(y1 - y2) < 10**-6:
                    index_map[i] = j

    return index_map


def check_counter(vertices, poly, L, edges):
    sumEdges = 0
    q = [0,0]
    pbc = [0,0]

    for i, index in enumerate(poly):

        x1, y1 = vertices[index] + q * L

        if i == len(poly) - 1:
            x2, y2 = vertices[poly[0]]
            indexv2 = poly[0]
        else:
            x2, y2 = vertices[poly[i+1]]
            indexv2 = poly[i+1]

        v1 = np.array([x1, y1])
        v2 = np.array([x2, y2])

        for edge in edges:
            i1 = edge[0]
            i2 = edge[1]
            if i1 == index and i2 == indexv2:
                pbc = edge[2:]
                q[0] = q[0] + pbc[0]
                q[1] = q[1] + pbc[1]

        v_next = v1 + periodic_diff(v2, v1, L, q)
        print(v1, v_next)
        x2, y2 = v_next
        # if qy != 0:
        #     sumEdges = sumEdges - ((x2 - x1) * (abs(y2) + abs(y1)))
        # else:
        sumEdges += (x2 - x1) * (y2 + y1)
        

    # print(sumqx)
    # print(sumqy)
    # if sumqy != 0:
    #      sumEdges = -sumEdges
    # print(sumEdges)

    if sumEdges < 0:
        return True

    if sumEdges >= 0:
        return False

    return True


# check if lists are the same order, but not necessarily
# indexed the same
def same_list(list1, list2):

    i0 = list1[0]

    i02 = -1
    for i, index in enumerate(list2):
        if i0 == index:
            i02 = i

    if i02 == -1:
        return False

    if i02 != -1:
        i2 = i02
        for i in range(0, len(list1)):

            if list1[i] != list2[i2]:
                return False

            i2 += 1
            if i2 >= len(list2):
                i2 = 0

        return True


# get polygons
def get_polygons(regions, index_map, vertices, L, edges):

    polygons = []

    for region in regions:
        count = 0
        poly = []
        for index in region:
            if index in index_map:
                count += 1
                poly.append(index_map[index])
        if count == len(region) and count != 0:
            polygons.append(poly)

    # assure all polygons in counter-clockwise order
    # This is not working yet because need to use periodic difference

    for i, poly in enumerate(polygons):
        cc = check_counter(vertices, poly, L, edges)
        print(poly)
        print(cc)
        if cc == False:
            rev_poly = []
            for index in reversed(poly):
                rev_poly.append(index)
            polygons[i] = rev_poly
            print(polygons[i])

    # remove duplicates
    new_polygons = []
    for poly in polygons:
        add = True
        for poly2 in new_polygons:
            if len(poly) == len(poly2):
                # check if they are the same ordered list
                same = same_list(poly, poly2)
                if same == True:
                    add = False
        if add == True:
            new_polygons.append(poly)

    return new_polygons
