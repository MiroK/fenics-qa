import dolfin
import numpy
import pylab
from dolfin import TrialFunction, TestFunction, inner, dx, Function, solve, assemble
from dolfin import PETScMatrix, PETScVector, SLEPcEigenSolver, assemble_system

# define an exact stream function
psi_exact_str = 'x[1]<=pi ? epsilon*cos(x[0])-(1.0/(cosh((x[1]-0.5*pi)/delta)*cosh((x[1]-0.5*pi)/delta)))/delta : epsilon*cos(x[0]) + (1.0/(cosh((1.5*pi-x[1])/delta)*cosh((1.5*pi-x[1])/delta)))/delta'   
epsilon = 0.05
delta = numpy.pi/15.0

psi_exact = dolfin.Expression(psi_exact_str,epsilon=epsilon,delta=delta)

# define the number of elements and polynomial degree of basis functions
n_elements = [5,10,20,40,80]

pp = numpy.arange(1,3)

# define the mesh
xBounds = numpy.array([0.,2.0*numpy.pi])
yBounds = numpy.array([0.,2.0*numpy.pi])

div_max = numpy.zeros([len(pp),len(n_elements)])

def project(expr, space):
    u, v = TrialFunction(space), TestFunction(space)
    a = inner(u, v)*dx
    L = inner(expr, v)*dx
    A, b = PETScMatrix(), PETScVector()
    assemble_system(a, L, A_tensor=A, b_tensor=b)

    uh = Function(space)
    x = uh.vector()
    solve(A, x, b, 'lu')

    lmax = SLEPcEigenSolver(A)
    lmax.parameters["spectrum"] = "largest magnitude"
    lmax.parameters["problem_type"] = "hermitian"
    lmax.solve(2)
    lmax = max([lmax.get_eigenpair(i)[0] for i in range(lmax.get_number_converged())])

    lmin = SLEPcEigenSolver(A)
    lmin.parameters["spectrum"] = "smallest magnitude"
    lmin.parameters["problem_type"] = "hermitian"
    lmin.solve(2)
    lmin = max([lmin.get_eigenpair(i)[0] for i in range(lmin.get_number_converged())])

    print space.dim(), 'Cond number', lmax/lmin



    
    return uh

for kp,p in enumerate(pp):
    for kn,n in enumerate(n_elements):
        print 'p = ' + str(p) + '   k = ' + str(n)
        mesh = dolfin.RectangleMesh(dolfin.Point(xBounds[0],yBounds[0],0.0),dolfin.Point(xBounds[1],yBounds[1],0.0),n,n,'crossed')

        # define the function spaces
        U = dolfin.FunctionSpace(mesh,'RT',p)      # velocity
        PSI = dolfin.FunctionSpace(mesh,'CG',p)    # stream function
        P = dolfin.FunctionSpace(mesh,'DG',p-1)    # divergence of velocity


        # compute the finite element approximation of the analytical stream function
        psi_h = project(psi_exact,PSI)

        # compute the velocity 
        u_h = project(dolfin.curl(psi_h), U)

        # compute the divergence 
        p_h = project(dolfin.div(u_h),P )


        # get the maximum value of the divergence
        div_max[kp,kn] = numpy.abs(p_h.vector().array()).max()

print div_max

pylab.semilogy(pp,div_max)
pylab.legend(n_elements,loc=4)
pylab.xlabel('p')
pylab.show()
