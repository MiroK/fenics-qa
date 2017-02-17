from dolfin import *
from mshr import *
import numpy as np


mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, 'CG', 2)
u = interpolate(Expression(('x[0]', 'x[0]*x[1]'), degree=2), V)

# ------------------
# Parameters
# ------------------
#set_log_level(ERROR) # log level
#parameters.parse()   # read paramaters from command line
# set some dolfin specific parameters
parameters["form_compiler"]["optimize"]     = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
# Geometry dimensions
L = 2 ; H = 2
# Material constants
E, nu = 210000.0, 0.3
mu, lamda = E/(2.0*(1.0+nu)), E*nu/((1.0+nu)*(1.0-2.0*nu))
# Loading
ut = 1. # reference value for the loading (imposed displacement)
load_min = 0. # load multiplier min value
load_max = 1. # load multiplier max value
load_steps = 50 # number of time steps
# Numerical parameters of the alternate minimization
maxiter = 300 
toll = 1.e-8
#mesh = Mesh("Mode_II_0pt005.xml")
ndim = 2 # get number of space dimensions
#-------------------
# Useful definitions
#-------------------
# zero and unit vectors
zero_v = Constant((0.,)*ndim)
e1 = Constant([1.,0.])
e2 = Constant([0.,1.])
# Strain
def eps(v):
    return sym(grad(v))    
# Eigenvalues are roots of characteristic polynomial
def eig_plus(A): 
    return (tr(A) + sqrt(tr(A)**2-4*det(A)))/2
def eig_minus(A): 
    return (tr(A) - sqrt(tr(A)**2-4*det(A)))/2
# Diagonal matrix with positive and negative eigenvalues as the diagonal elements    
def diag_eig(A):
    lambdap1 = 0.5*(eig_plus(A)+ abs(eig_plus(A)))
    lambdap2 = 0.5*(eig_minus(A) + abs(eig_minus(A)))
    lambdan1 = 0.5*(eig_plus(A) - abs(eig_plus(A)))
    lambdan2 = 0.5*(eig_minus(A) - abs(eig_minus(A)))
    matp = as_matrix(((lambdap1,0),(0,lambdap2)))
    matn = as_matrix(((lambdan1,0),(0,lambdan2)))
    return (matp,matn)    
# Eigenvectors of the matrix arranged in columns of matrix     
def eig_vecmat(A):
     lambdap = eig_plus(A)
     lambdan = eig_minus(A)
     a = A[0,0]
     b = A[0,1]
     c = A[1,0]
     d = A[1,1]
     v11 = lambdap - d
     v12 = lambdan - d
     nv11 = sqrt(v11**2 + c**2)
     nv12 = sqrt(v12**2 + c**2)
     a1 = v11/nv11
     b1 = v12/nv12
     c1 = c/nv11
     d1 = c/nv12
     v21 = lambdap - a
     v22 = lambdan - a
     nv21 = sqrt(v21**2 + b**2)
     nv22 = sqrt(v22**2 + b**2)
     A1 = b/nv21
     B1 = b/nv22
     C1 = v21/nv21
     D1 = v22/nv22
     tol = 1.e-10
     Eigvecmat = conditional(gt(abs(c), tol) ,as_matrix(((a1,b1),(c1,d1))), conditional(gt(abs(b), tol),as_matrix(((A1,B1),(C1,D1))),Identity(2)))
     return Eigvecmat
def eps_split(A):
    P = eig_vecmat(A)
    (diag_eigp,diag_eign) = diag_eig(A)
    epsp = dot(P,dot(diag_eigp,P.T))
    epsn = dot(P,dot(diag_eign,P.T))
    return (epsp,epsn)  

def psi(A,v):
    tol = 1.e-10
    epst = eps(v)
    (epsp,epsn) = eps_split(epst)
    H = conditional(gt(tr(A),tol),1,0)
    psip = 0.5*lamda*(tr(A)*H)**2 + mu*inner(epsp,epsp)
    psin = 0.5*lamda*(tr(A)-tr(A)*H)**2 + mu*inner(epsn,epsn)
    psit = 0.5*lamda*(tr(A)*H)**2 + mu*inner(epst,epst)
    return (psip,psin,psit)

def sigma_0(A,v):
    epst = eps(v)
    (epsp,epsn) = eps_split(epst)
    tol = 1.e-10
    H = conditional(gt(tr(A),tol),1,0)
    sigmap = 2.0*mu*(epsp) + lamda*(tr(A)*H)*Identity(ndim)
    sigman = 2.0*mu*(epsn) + lamda*(tr(A)*H)*Identity(ndim)
    sigmat = 2.0*mu*(epst) + lamda*(tr(epst))*Identity(ndim)
    return (sigmap,sigman,sigmat)
#----------------------------------------------------------------------------
# Define boundary sets for boundary conditions
#----------------------------------------------------------------------------
class Left(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-5 # tolerance for coordinate comparisons
        return on_boundary and abs(x[0] - 0) < tol         
class Right(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1E-5 # tolerance for coordinate comparisons
        return on_boundary and abs(x[0] - 2) < tol
# Initialize sub-domain instances
left = Left() 
right = Right()
# define meshfunction to identify boundaries by numbers
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1) # mark left as 1
right.mark(boundaries, 2) # mark right as 2
# Define new measure including boundary naming 
ds = Measure("ds")[boundaries] # left: ds(1), right: ds(2)    
# Define function spaces
V_u = VectorFunctionSpace(mesh, "CG", 1)

# Define the function, test and trial fields
u, du, v = Function(V_u), TrialFunction(V_u), TestFunction(V_u)
# Define boundary condiition                                  
zero = Constant(0.0)
zero2 = Constant((0.0,0.0))
# Dirichlet boundary condition for a traction test boundary
u_L = zero2
U_roll = zero
u_R = Expression(("0.0","0.02*t"),t = 0.,degree=1) 
# Dispalcement boundary conditions
bc1 = DirichletBC(V_u, u_L, boundaries, 1) 
bc2 = DirichletBC(V_u, u_R, boundaries, 2) 
bc_u = [bc1, bc2]
# Energy functions
strain  = eps(u)
(psip,psin,psit) = psi(strain,u)
(sp,sn,st) = sigma_0(strain,u)
elastic_energy = (psip + psin)*dx 
# Weak form 
E_u = derivative(elastic_energy,u,v)
E_u_u = derivative(E_u,u,du)

# Writing tangent problems in term of test and trial functions for matrix assembly
E_du = replace(E_u,{u:du})

# Variational problem for the displacement
problem_u = NonlinearVariationalProblem(E_u, u, bc_u, E_u_u)

