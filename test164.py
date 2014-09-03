from dolfin import *

if has_linear_algebra_backend("Epetra"):
    parameters["linear_algebra_backend"] = "Epetra"

class East(SubDomain):
    def inside(self, x , on_boundary):
        return near(x[0], 1.0)

class West(SubDomain):
    def inside(self, x , on_boundary):
        return near(x[0], 0.0)

class North(SubDomain):
    def inside(self, x , on_boundary):
        return near(x[1], 1.0)

class South(SubDomain):
    def inside(self, x , on_boundary):
        return near(x[1], 0.0)

if __name__ == '__main__':

    mesh = UnitSquareMesh(64, 64)
    File("mesh.pvd") << mesh

    V = FunctionSpace(mesh, "Lagrange", 1)

    bcs = [DirichletBC(V, Constant(150), North()), \
           DirichletBC(V, Constant(250), East()), \
           DirichletBC(V, Constant(100), South()), \
           DirichletBC(V, Constant(50), West())]

    u = Function(V)
    v = TestFunction(V)
    F = inner(grad(u), grad(v))*dx

    # Compute solution
    solve(F == 0, u, bcs, solver_parameters={"newton_solver":
                                            {"relative_tolerance": 1e-6}})

    # Plot solution and solution gradient
    plot(u, title="Solution")
    File('fausto.xdmf') << u
    interactive()
