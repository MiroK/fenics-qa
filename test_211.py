from dolfin import *

mesh = UnitCubeMesh(10,10,10) 

V1 = VectorFunctionSpace(mesh,'P',2) 

class Conductor(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[2], (0.4, 0.6)) and between(x[1], (0.4, 0.6)) and between(x[0], (0.4, 0.6)))

conductor = Conductor()

domain_marker =  CellFunction("size_t", mesh) 
domain_marker.set_all(0)
conductor.mark(domain_marker, 1)
dx = Measure('dx', domain=mesh, subdomain_data=domain_marker)

f1 = Expression(('sin(x[0])','sin(x[1])','sin(x[2])'), degree=4)
f2 = Expression(('sin(x[2])','sin(x[0])','sin(x[1])'), degree=4)

class Source(Expression):
    def __init__(self,domain_marker, f1, f2, **kwargs):
        self.domain_marker = domain_marker
        self.f1 = f1
        self.f2 = f2

    def eval_cell(self, values, x, cell):
        if self.domain_marker[cell.index] == 0:
            values[:] =  self.f1(x)
        else:
            values[:] =  self.f2(x)

    def value_shape(self): return (3, )

f = Source(domain_marker,f1,f2, degree=4)  

plot(interpolate(f,V1),interactive = True)
