from fenics import *

class material_parameters_rank4(Expression):
    def __init__(self, cell_function, params, domain):
        self.cell_function = cell_function
        self.params = params
        self._domain = domain
    def eval_cell(self, values, x, cell):
        domain_id = self.cell_function[cell.index] - 1
        values = self.params[domain_id]
    def value_shape(self):
        return (3,3,3,3)

xl,yl,zl = 500.,100.,100.
mesh = BoxMesh(Point(0, 0, 0), Point(xl, yl, zl), 25, 5, 5)
cells = CellFunction('size_t', mesh)
facets = FacetFunction('size_t', mesh)
N = FacetNormal(mesh)

Space = VectorFunctionSpace(mesh, 'P', 1)
dA = Measure('ds', domain=mesh, subdomain_data=facets)
dV = Measure('dx', domain=mesh, subdomain_data=cells)

left = CompiledSubDomain('near(x[0], 0) && on_boundary')
right = CompiledSubDomain('near(x[0], length) && on_boundary', length=xl)
facets.set_all(0)
right.mark(facets, 1)
left.mark(facets, 2)
bc1 = DirichletBC(Space.sub(1), 1.0, facets, 1)
bc2 = DirichletBC(Space, Constant((0.0, 0.0, 0.0)), facets, 2)
bc = [bc1,bc2]

i, j, k, l = indices(4)
delta = Identity(3)
C_b = as_tensor(10E10*delta[i,j]*delta[k,l] + 1E10*delta[i,k]*delta[j,l] + 1E10*delta[i,l]*delta[j,k], (i,j,k,l))

#mat1 = CompiledSubDomain('x[0] > half',half=xl/2.)
#cells.set_all(1)
#mat1.mark(cells,2)

DG0 = FunctionSpace(mesh, 'DG', 0)
chi = interpolate(Expression('x[0] > half ? 2 : 1', half=xl/2.), DG0)

def tensor_conditional(predicate, tvalue, fvalue):
    assert tvalue.ufl_shape == fvalue.ufl_shape

    # Scalar
    rank = len(tvalue.ufl_shape)
    if rank == 0: 
        return conditional(predicate, tvalue, fvalue)
    # Vector
    elif rank == 1:
        dim = tvalue.ufl_shape[0]
        conds = [conditional(predicate, tvalue[i], fvalue[i]) for i in range(dim)]
    # Matrices
    elif rank == 2:
        nrows, ncols = tvalue.ufl_shape
        conds = [[conditional(predicate, tvalue[i, j], fvalue[i, j]) for j in range(ncols)]
                 for i in range(nrows)]
    # Generalize later
    else:
        raise ValueError

    return as_tensor(conds)

#C =  material_parameters_rank4(cells, [C_b, C_b], domain=mesh)
C = tensor_conditional(eq(chi, 2), C_b, C_b)
# C=C_b

du = TrialFunction(Space)
del_u = TestFunction(Space)
u = Function(Space)

eps = sym(grad(u))
tau = as_tensor(C[j,i,k,l]*eps[k,l] , (j,i))

F_u =  tau[j,i]*del_u[i].dx(j) *(dV(1)+dV(2))

Form = F_u 
Gain = derivative(Form, u, du)

solve(Form == 0, u, bc, J=Gain, \
    solver_parameters={"newton_solver":{"linear_solver": "mumps", "relative_tolerance": 1e-3} }, \
    form_compiler_parameters={"cpp_optimize": True, "representation": "quadrature", "quadrature_degree": 2}  )
