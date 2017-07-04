from dolfin import *

def my_adaptive_solver(u_e, f, mesh):
    V = FunctionSpace(mesh, "Lagrange", 1)

    bc = DirichletBC(V, u_e, DomainBoundary())

    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v))*dx()
    L = f*v*dx

    # Define function for the solution
    u = Function(V)
    # Define goal functional (quantity of interest)
    M = u*dx()
    tol = 1.e-4

    problem = LinearVariationalProblem(a, L, u, bc)
    solver = AdaptiveLinearVariationalSolver(problem, M)
    solver.parameters['max_iterations'] = 20
    solver.solve(tol)
    #solver.summary()
    #data = solver.adaptive_data()

    u_h_list = []; u_h_list.append(u.root_node())
    num_dofs_list = []; num_dofs_list.append(u.root_node().function_space().dim())
    temp = u.root_node()
    while temp.has_child():
        V = temp.child().function_space()
        u_h_list.append(Function(V, temp.child().vector().copy()))
        num_dofs_list.append(V.dim())
        temp = temp.child()
##############################
    error = [errornorm(u_e, u_h_list[i], norm_type='H10') for i in range(len(u_h_list))]
    h = [num_dofs_list[i]**(-1.0/2) for i in range(len(num_dofs_list))]
    rate = []
    from math import log as ln
    for i in range(len(u_h_list) - 1):
        rate.append(ln(error[i]/error[i+1]) / ln(h[i]/h[i+1]))
    print rate
##############################

    return u_h_list, num_dofs_list


if __name__ == '__main__':
    import pdb
    pdb.set_trace()
    mesh = UnitSquareMesh(2, 2)
    u_e = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
    f = Constant(-6.0)

    u_h_list, num_dofs_list = my_adaptive_solver(u_e, f, mesh)

##############################
    error = [errornorm(u_e, u_h_list[i], norm_type='H10') for i in range(len(u_h_list))]
##############################
