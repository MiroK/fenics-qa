# Example 1:
from dolfin import *
import numpy as np

class Linear_Elliptic_Example1(object):
    '''
    This is the Example_Linear_Elliptic_Example1 class
    '''
def __init__(self):
    self.mesh = mesh
    self.boundary = DirichletBoundary()
    self.boundary_value = Linear_Elliptic_Example1_solution()
    self.solution = Linear_Elliptic_Example1_solution()
    self.source =Linear_Elliptic_Example1_source()
    return None

def __str__(self):
    return 'Linear_Elliptic_Example1'

#define mesh
mesh = UnitSquareMesh(10,10, "crossed")
boundaries = FacetFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)

aa = Constant(100.0)
c = Constant(1.0)

# compute the exact solution and the correspunding source function
class Linear_Elliptic_Example1_solution(Expression):
    def __init__(self):
        return None
    def eval(self, values, x):
        values[0] = (1./aa)*((x[1]-x[0])*(1-x[0]-x[1])*(x[0]-0.50)**2\
            *(x[1]-0.5)**2)

class Linear_Elliptic_Example1_source(Expression):
    def __init__(self):
        return None
    def eval(self, values, x):
        values[0] = (1./aa)*((x[1]-x[0])*(1-x[0]-x[1])*(x[0]-0.50)**2\
            *(x[1]-0.5)**2) + 4*(1-x[0]-x[1])*(x[0]-0.5)*(x[1]-0.5)**2\
            + 4*(x[1]-x[0])*(x[0]-0.5)*(x[1]-0.5)**2-2*(x[1]-x[0])\
            *(1-x[0]-x[1])*(x[1]-0.5)**2-4*(1-x[0]-x[1])*(x[0]-0.5)**2\
            *(x[1]-0.5)+4*(x[1]-x[0])*(x[0]-0.5)**2*(x[1]-0.5)\
            -2*(x[1]-x[0])*(1-x[0]-x[1])*(x[0]-0.5)**2

# aa* grad(u)
class Linear_Elliptic_Example1_Gradiant(Expression):
    def __init__(self):
        return None

    def eval(self, values, x):
        values[0] = -(1.-x[0]-x[1])*(x[0]-0.5)**2*(x[1]-0.5)**2-(x[1]-x[0])*(x[0]-0.5)**2\
            *(x[1]-0.5)**2 + 2.*(x[1]-x[0])*(1-x[0]-x[1])*(x[0]-0.5)*(x[1]-0.5)**2
        values[1] = (1.-x[0]-x[1])*(x[0]-0.5)**2*(x[1]-0.5)**2-(x[1]-x[0])*(x[0]-0.5)**2\
            *(x[1]-0.5)**2 + 2.*(x[1]-x[0])*(1-x[0]-x[1])*(x[0]-0.5)**2*(x[1]-0.5)

    def value_shape(self):
        return (2,)

# Define function G such that G \cdot n = u_n
class BoundarySource(Expression):
    def __init__(self, mesh):
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.cell_normal()
        #cell.normal(ufc_cell.local_facet)
        g = Linear_Elliptic_Example1_Gradiant()
        values[0] = g[0]*n[0] + g[1]*n[1]

u_n = BoundarySource(mesh)


# Define boundary conditions

class DirichletBoundary_top(SubDomain):
    def inside(self ,x, on_boundary):
        return x[1] < DOLFIN_EPS

class DirichletBoundary_bottom(SubDomain):
    def inside(self ,x, on_boundary):
        return x[1] > 1.0 - DOLFIN_EPS

Gamma_D1 = DirichletBoundary_top()
Gamma_D1.mark(boundaries, 2)
Gamma_D2 = DirichletBoundary_top()
Gamma_D2.mark(boundaries, 4)

class NeumannBoundary_left(SubDomain):
    def inside(self, x, ob_boundary):
        return x[0] < DOLFIN_EPS

class NeumannBoundary_right(SubDomain):
    def inside(self, x, ob_boundary):
        return x[0] > 1.0 - DOLFIN_EPS

Gamma_N1 = NeumannBoundary_left()
Gamma_N1.mark(boundaries, 1)
Gamma_N2 = NeumannBoundary_right()
Gamma_N2.mark(boundaries, 3)

# Problem
V = FunctionSpace(mesh, "DG", 3)
source = Linear_Elliptic_Example1_source()
alpha = Constant(40)
h_T = CellSize(mesh)
h_T_avg = (h_T('+') + h_T('-'))/2.0
n = FacetNormal(mesh)

u = TrialFunction(V)
v = TestFunction(V)

a = dot(aa*grad(u),grad(v))*dx- dot(avg(aa*grad(u)),n('+'))*jump(v)*dS\
        + jump(u)*dot(avg(aa*grad(v)),n('+'))*dS\
        + alpha('+')/h_T_avg*jump(u)*jump(v)*dS

L = source*v*dx + u_n*v*ds(1)# + u_n*v*ds(3)

u_top = Expression('(x[0]*(1.-x[0])*(x[0]-0.5)*(x[0]-0.5))/400.')
u_bottom = Expression('(x[0]*(x[0]-1.)*(x[0]-0.5)*(x[0]-0.5))/400.')

# Boundary conditions
bcs =[DirichletBC(V, u_top, Gamma_D1),
    DirichletBC(V, u_bottom, Gamma_D2)]

# Compute solution
u_h = Function(V)

solve(a == L, u_h, bcs)
