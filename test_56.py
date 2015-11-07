from dolfin import *
from mshr import *

domain = Rectangle(Point(0, 0), Point(1, 1))
mesh = generate_mesh(domain, 20)

h = CellVolume(mesh)

S = FunctionSpace(mesh, 'DG', 0)
# One way to get sizes
size0 = project(h, S)

# Now by element-by-element way
size1 = Function(S)
size_1_vec = size1.vector()

dofmap = S.dofmap()
elements = CellFunction('size_t', mesh, 0)
dK = Measure('dx', subdomain_id=1, domain=mesh, subdomain_data=elements)
form = Constant(1)*dK

for cell in cells(mesh):
    # Set to zero
    elements.set_all(0)
    # Light up me
    elements[cell] = 1
    # Integrate
    size = assemble(form)
    # Assign to dof of this cell
    dof = dofmap.cell_dofs(cell.index())[0]
    size_1_vec[dof] = size

# Check that they are the same
size0.vector().axpy(-1, size_1_vec)
print ':)' if size0.vector().norm('linf') < DOLFIN_EPS else ':('
