from __future__ import print_function, division
from fenics import *
from mshr import *

ysize = 85 # length of plates in y
wall_distance = 30 # distance to the simulation box wall

simulation_box = Box(Point(-wall_distance, - ((85/2)+wall_distance), -(16.5+wall_distance)),
             Point(100, ((85/2)+wall_distance), 16.5 + wall_distance)
             )

pcb = (Box(Point(0, -85/2, -1.5/2), Point(70, 85/2, 1.5/2)) -
       Box(Point(25, -50/2, -1.5/2), Point(70, 50/2, 1.5/2))
       )

domain = simulation_box - pcb

mesh = generate_mesh(domain, 32)

solution_file = File("solution.pvd")

# Define the boundaries
x_min = 'near(-30, x[1], x[2])' 
x_max = 'near(100, x[1], x[2])'

y_min = 'near(x[0], -85/2, x[2])'
y_max = 'near(x[0], 85/2, x[2])'

z_min = 'near(x[0], x[1], -(16.5+30))'
z_max = 'near(x[0], x[1], (16.5+30))'

walls = [x_min, x_max, y_min, y_max, z_min, z_max]

pcb_boundary = 'on_boundary && x[0]> 0-1 && x[0]< 70+1  && x[1]> -85/2-1 && x[1]< 85/2 && x[2]> -1.5/2-1 && x[2] < 1.5/2+1'

V = FunctionSpace(mesh, 'Lagrange', 1)
bc_box = [DirichletBC(V, Constant(0.0), wall) for wall in walls]
bc_pcb = DirichletBC(V, Constant(700.0), pcb_boundary)

bcs = bc_box.append(bc_pcb)

u = TrialFunction(V)
v = TestFunction(V)
a  = inner(grad(u), grad(v))*dx
L = v*dx


u = Function(V)
solve(a == L, u, bcs)

solution_file << u

