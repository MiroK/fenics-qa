from dolfin import *
import math
import sympy as sym

T = 1.0
theta = 0.5

def derive_rhs():
    x,y,t, pi = sym.symbols(["x[0]", "x[1]", "t", "pi"])
    u_ana_sympy = t**3*sym.cos(x*pi)*sym.cos(y*pi)

    f_sympy = sym.diff(u_ana_sympy, t) - sym.diff(u_ana_sympy, x, 2) - sym.diff(u_ana_sympy, y, 2)
    print "u = ", u_ana_sympy
    print "f = ", f_sympy
    return u_ana_sympy, f_sympy

u_ana_sympy,f_sympy = derive_rhs()

def get_error(n):
    nt = 2**n
    dt = Constant(T/float(nt))

    mesh = UnitSquareMesh(2**n, 2**n)

    u_ana = Expression(sym.ccode(u_ana_sympy), t = 0, cell=mesh.ufl_cell(),
            degree=4)
    f = Expression(sym.ccode(f_sympy), t = 0, degree=4)   


    V = FunctionSpace(mesh, "CG", 1)
    ut = TrialFunction(V)
    v = TestFunction(V)

    u0 = Function(V)
    u1 = Function(V)

    n = FacetNormal(mesh)
    u_theta = ut*theta + u0*(1 - theta)
    dvdt = (ut - u0)/dt

    F = inner(dvdt, v)*dx + inner(grad(u_theta), grad(v))*dx - inner(f,v)*dx#\
        #-Constant(theta)*inner(dot(grad(u_ana), n), v)*ds\
        #-Constant(1-theta)*inner(dot(grad(u0), n), v)*ds

    a = lhs(F)
    b = rhs(F)

    bcs = []
    #bcs = DirichletBC(V, u_ana, "on_boundary")

    t = 0
    errors = []

    u0.interpolate(u_ana)
    for i in range(nt -1):
        f.t = t + theta*float(dt)

        t += float(dt)  
        u_ana.t = t

        solve(a == b, u1, bcs)

        u0.assign(u1)
        
        u_ana.t = t
        errors.append(errornorm(u_ana, u1))     
    plot(u_ana, mesh=mesh)
    plot(u1)
    return float(dt)*sum(errors)

def convergence_order(errors, base = 2):
    orders = [0.0] * (len(errors) - 1)
    for i in range(len(errors) - 1):
        if errors[i+1] ==0: errors[i] = 0.0
        else:
            ratio = errors[i]/errors[i+1]
            if ratio == 0: orders[i] = 0
            else:
                orders[i] = math.log(ratio, base)
    return orders

def run_test(): 
    errors = []
    for n in range(2,6):
        errors.append(get_error(n))

    conv_rates = convergence_order(errors)

    print "Error     Rate"
    for e,r in zip(errors,[0] + conv_rates):
        print "{:.3e} {:.3e}".format(e,r)
run_test()

interactive()
