from numpy import *
from dolfin import *

# Mobius strip parameters
# Number of times you'd like your strip to twist (only even numbers allowed!)
num_twists = 4
# The width of the strip
width = 0.5

# Mesh parameters
# Number of nodes along the length of the strip
nl = 300
# Number of nodes along the width of the strip (>= 2)
nw = 14

# Generate suitable ranges for parameters
u_range = arange(nl, dtype='d')/(nl)*2*pi
v_range = arange(nw, dtype='d')/(nw - 1.0)*width

# Create the mesh
mesh = Mesh()
editor = MeshEditor()
editor.open(mesh, "triangle", 2, 3)

# Create a list to store the vertices
editor.init_vertices(nl*nw)

# Populate the list of vertices
j = 0
for u in u_range:
    for v in v_range:
        editor.add_vertex(j,  cos(u) + v*cos(num_twists*u/2.0)*cos(u), \
                              sin(u) + v*cos(num_twists*u/2.0)*sin(u), \
                              v*sin(num_twists*u/2.0))
        j = j + 1

# Create a list to store the cells
editor.init_cells(nl*(nw - 1)*2)

# Populate the list of cells
k = 0
for i in range(nl - 1):
    for j in range(nw - 1):
        editor.add_cell(k    , i*nw + j, (i + 1)*nw + j + 1, i*nw + j + 1) 
        editor.add_cell(k + 1, i*nw + j, (i + 1)*nw + j    , (i + 1)*nw + j + 1)
        k = k + 2
# Close the geometry
for j in range(nw - 1):
    editor.add_cell(k    , (nl - 1)*nw + j, j + 1, (nl - 1)*nw + j + 1) 
    editor.add_cell(k + 1, (nl - 1)*nw + j, j    , j + 1)
    k = k + 2

editor.close()

File('moebius.xml.gz') << mesh
