from dolfin import *
#from dolfin_adjoint import *
set_log_level(INFO)

vel_angular = 1
omega = Constant((0.0, 0.0, vel_angular))
rho = 1.0

rect = Rectangle(0.0,0.0,1.0,1.0)
mesh = UnitSquareMesh(10, 10)

U = VectorFunctionSpace(mesh, "CG", 3) # velocity function space    -U
P = FunctionSpace(mesh, "CG", 1)       # pressure function space    -P
W = MixedFunctionSpace([U, P])         # mixed Taylor-Hood function space

def BC_1(x, on_boundary):
    return on_boundary and \
        near(x[0],0.0) or near(x[1],1.0) or near(x[0],1.0)

def BC_2(x, on_boundary):
    return on_boundary and \
        near(x[1],0.0)

Gamma1 = DirichletBC(W.sub(1),0.0,BC_1)
Gamma2 = DirichletBC(W.sub(1),10.0,BC_2)

bcs = [Gamma1, Gamma2]

def forward(rho):
    w = Function(W)
    (u, p) = split(w)
    (v, q) = TestFunctions(W)

    #F = (2.0*rho*inner(cross(omega,u),v)*dx)


    F=(2.0*rho*inner(as_vector([omega[2]*u[1], omega[2]*u[0]]) , v)*dx )

    solve(F == 0, w, bcs)
    return w

if __name__ == "__main__":
    w   = forward(rho)
    (u, p) = split(w)
