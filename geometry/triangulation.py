# -*- coding: utf8 -*- 
"""
Constrained delaunay implementation based on 

Domiter, V. and Žalik, B.
Sweep‐line algorithm for constrained Delaunay triangulation

(this implementation is not complete)
"""

import numpy as np
from math import cos

# Distance between points A and B
def distance(A,B):
    n = len(A)
    assert len(B) == n
    return sum((A[i]-B[i])**2 for i in xrange(n))**0.5

# Cosine of angle ABC
def cosine(A,B,C):
    a,b,c = distance(B,C), distance(A,C), distance(A,B)
    return (a*a+c*c-b*b)/(2*a*c)

# Cartesian coordinates of the point whose barycentric coordinates
# with respect to the triangle ABC are [p,q,r]
def barycentric(A,B,C,p,q,r):
    n = len(A)
    assert len(B) == len(C) == n
    s = p+q+r
    p, q, r = p/s, q/s, r/s
    return tuple([p*A[i]+q*B[i]+r*C[i] for i in xrange(n)])

# Cartesian coordinates of the point whose trilinear coordinates
# with respect to the triangle ABC are [alpha,beta,gamma]
def trilinear(A,B,C,alpha,beta,gamma):
    a = distance(B,C)
    b = distance(A,C)
    c = distance(A,B)
    return barycentric(A,B,C,a*alpha,b*beta,c*gamma)
               
# Cartesian coordinates of the circumcenter of triangle ABC
def circuminfo(A,B,C):
    cosA = cosine(C,A,B)
    cosB = cosine(A,B,C)
    cosC = cosine(B,C,A)
    cc = trilinear(A,B,C,cosA,cosB,cosC)
    # returns circumcenter and circumradius
    return cc, distance(cc, A)

# Check if the points lie in counter-clockwise order or not
def iscounterclockwise(a, b, c):
    A = np.array(list(pts[a]))
    B = np.array(list(pts[b]))
    C = np.array(list(pts[c]))
    return np.cross(B-A, C-B)>0

def isintersects(edge1, edge2):
    A = np.array(list(pts[edge1[0]]))
    B = np.array(list(pts[edge1[1]]))
    C = np.array(list(pts[edge2[0]]))
    D = np.array(list(pts[edge2[1]]))

    E = B-A
    F = D-C
    P = np.array([-E[1], E[0]])
    h = ((A-C).dot(P))/F.dot(P)
    return (h>=0 and h<1)

#user input data - points and constraining edges
pts = [(0, 0),
       (10, 0),
       (10, 10),
       (20, 10),
       (20, 20),
       (10, 17),
       (9, 30),
       (5, 25),
       (6, 15),
       (10, 12),
       (0, 5)]
edges = [(4,1), (6,7)]

pts = np.array(pts, dtype=[('x', float), ('y', float)])
edges = np.array(edges, dtype=int)
edges_lookup = {}
## Initialization (sec. 3.3)

# sort by y, then x
order = np.argsort(pts, order=('y', 'x'))
pts = pts[order]
# update display edges to match new point order
invorder = np.argsort(order)
edges = invorder[edges]

# make artificial points P-1 and P-2
xmax = pts['x'].max()
xmin = pts['x'].min()
ymax = pts['y'].max()
ymin = pts['y'].min()
xa = (xmax-xmin) * 0.3
ya = (ymax-ymin) * 0.3
p1 = (pts['x'].min() - xa, pts['y'].min() - ya)
# error in the equation in the paper, should be x_max+delta_x, not -delta_x
p2 = (pts['x'].max() + xa, pts['y'].min() - ya)

# prepend artificial points to point list
newpts = np.empty(len(pts)+2, dtype=pts.dtype)
newpts[0] = p1
newpts[1] = p2
newpts[2:] = pts
pts = newpts
edges += 2

# find topmost point in each edge
tops = set(edges.max(axis=1))

# inintialize sweep front
front = [0, 2, 1]
tris = []


# visual debugging: draw edges, front, triangles
import pyqtgraph as pg
import time
app = pg.mkQApp()
win = pg.plot()
gpts = pts.view(float).reshape(len(pts), 2)
graph = pg.GraphItem(pos=gpts, adj=edges, pen={'width': 3, 'color':(0,100,0)})
win.addItem(graph)
front_line = pg.PlotCurveItem(pen={'width': 2, 'dash': [5, 5], 'color': 'y'})
win.addItem(front_line)
tri_shapes = []
def draw_state():
    front_pts = pts[np.array(front)]
    front_line.setData(front_pts['x'], front_pts['y'])
    for i in range(100):  # sleep ~1 sec, but keep ui responsive
        app.processEvents()
        time.sleep(0.01)
def draw_tri(tri):
    tpts = pts[np.array(tri)]
    path = pg.arrayToQPath(tpts['x'], tpts['y'])
    shape = pg.QtGui.QGraphicsPathItem(path)
    shape.setPen(pg.mkPen(255,255,255,100))
    shape.setBrush(pg.mkBrush(0,255,255,50))
    win.addItem(shape)
    tri_shapes.append(shape)
draw_state()

## Legalize recursively - incomplete
def legalize((f00, f11, p)):
    print "Legalizing points = {}, {}, {}".format(f00, f11, p)
    a = pts[f00]
    b = pts[f11]
    c = pts[p]
    cc, cr = circuminfo(a, b, c)
    for point in pts:
        if (point==a or point==b or point==c):
            continue
        elif distance(cc, point) < cr:
           print "Illegal point"
           print point
    return (f00, f11, p)


## Sweep (sec. 3.4)

def add_tri(f0, f1, p):
    # todo: legalize!
    global front, tris
    if iscounterclockwise(front[f0], front[f1], p):
        edges_lookup[(front[f0], front[f1])] = p
        edges_lookup[(front[f1], p)] = front[f0]
        edges_lookup[(p, front[f0])] = front[f1]
    else:
        edges_lookup[(front[f1], front[f0])] = p
        edges_lookup[(p, front[f1])] = front[f0]
        edges_lookup[(front[f0], p)] = front[f1]
    tri = legalize((front[f0], front[f1], p))
    tris.append(tri)
    front.insert(f1, p)
    print "Inserted at position f1={} value p={}".format(f1,p)
    print "front = {}".format(front)
    # visual debugging
    draw_tri(tris[-1])

for i in range(3, len(pts)):
    pi = pts[i]
    # First, triangulate from front to new point
    # This applies to both "point events" (3.4.1) and "edge events" (3.4.2).

    # get index along front that intersects pts[i]
    l = 0
    while pts[front[l+1]]['x'] <= pi['x']:
        l += 1
    pl = pts[front[l]]
    pr = pts[front[l+1]]
    if pi['x'] > pl['x']:  # "(i) middle case"
        # Add a single triangle connecting pi,pl,pr
        add_tri(l, l+1, i)
    else: # "(ii) left case"
        ps = pts[l-1]
        # Add triangles connecting pi,pl,ps and pi,pl,pr
        add_tri(l, l+1, i)
        add_tri(l-1, l, i)
        front.pop(l+1)
        
    # todo: Continue adding triangles to smooth out front
    # (use heuristics shown in figs. 9, 10)
                    
    if i in tops: # this is an "edge event" (sec. 3.4.2)
        print "Locate first intersected triangle"
        # Locate the other endpoint
        found = False
        for e in edges:
            if i in e:
                found = True
                break
        if not found:
            print "Other end point not located"
            continue
        # (i) locate intersected triangles
        """
        If the first intersected triangle contains the top point,
        then start traversal there. Also, if an intersected triangle
        contains the top point, then it has to be the first intersected
        triangleself.
        """
        print edges_lookup
        vals = edges_lookup.values()
        intersects = False
        for value in vals:
            if value==i:
                current_edge = edges_lookup.keys()[vals.index(i)]
                if isintersects(current_edge, e):
                    intersects = True
                    break

        if intersects:
            print current_edge+(i,)
            # now recursively check and remove all intersecting triangles
            pass
        else :
            # perform other intersection tests
            pass

        # (ii) remove intersected triangles
        # (iii) triangluate empty areas

    draw_state()

## Finalize (sec. 3.5)

# (i) Remove all triangles that include at least one artificial point
# (ii) Add bordering triangles to fill hull
#print edges_lookup
raw_input()