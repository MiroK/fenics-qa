from dolfin import *
import numpy as np
 
import sys

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
 
discret=2
mesh =  UnitCubeMesh(discret,discret,discret)
 
dolfin.plot(mesh, "3D mesh")
 
FS = FunctionSpace(mesh, "CG", 2)
 
W = MixedFunctionSpace([FS, FS, FS, FS, FS, FS, FS, FS, FS, FS,
                        FS, FS, FS, FS, FS, FS, FS, FS, FS, FS,
                        FS, FS, FS, FS, FS, FS, FS, FS, FS, FS])
 
(u0, u1, u2, v000, v001, v002, v010, v011, v012, v020, v021, v022, v100, v101, v102, v110, v111, v112, v120, v121, v122, v200, v201, v202, v210, v211, v212, v220, v221, v222) = TrialFunctions(W)
(du0, du1, du2, dv000, dv001, dv002, dv010, dv011, dv012, dv020, dv021, dv022, dv100, dv101, dv102, dv110, dv111, dv112, dv120, dv121, dv122, dv200, dv201, dv202, dv210, dv211, dv212, dv220, dv221, dv222) = TestFunctions(W)



# Lower face
def DB_0(x, on_boundary):
    return x[2] < DOLFIN_EPS and on_boundary
 
# Point on upper face
TOL = 1e-3
class Pinpoint(SubDomain):
    def __init__(self, coords):
        self.coords = np.array(coords)
        SubDomain.__init__(self)
    def move(self, coords):
        self.coords[:] = np.array(coords)
    def inside(self, x, on_boundary):
        return np.linalg.norm(x-self.coords) < TOL
 
pinpoint = Pinpoint([0.5, 0.5, 1.0])
 
# Neumann boundary: everything but the lower face
 
class Neumann(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary or \
        x[0] > 1.0 - DOLFIN_EPS and on_boundary or \
        x[1] < DOLFIN_EPS and on_boundary or \
        x[1] > 1.0 - DOLFIN_EPS and on_boundary  or \
        x[2] > 1.0 - DOLFIN_EPS and on_boundary
 
subdomains = MeshFunction("double", mesh)
neumann = Neumann()
neumann.mark(subdomains, 1)
 
# Displacement boundary condition on Dirichlet boundary (lower face and pinpoint)
u0_0=Expression('0.0')
u1_0=Expression('0.0')
u2_0=Expression('0.0')
u0_1=Expression('0.0')
u1_1=Expression('0.0')
u2_1=Expression('-0.01')
 
bc = [
DirichletBC(W.sub(0), u0_0, DB_0),
DirichletBC(W.sub(1), u1_0, DB_0),
DirichletBC(W.sub(2), u2_0, DB_0),
DirichletBC(W.sub(0), u0_1, pinpoint, 'pointwise'),
DirichletBC(W.sub(1), u1_1, pinpoint, 'pointwise'),
DirichletBC(W.sub(2), u2_1, pinpoint, 'pointwise')]
 
# Isotropic stiffness tetradic and hexadic
delta=[[1,0,0],[0,1,0],[0,0,1]]
lambada=100.0
mu=100.0
c1=0.0
c2=0.0
c3=0.0
c4=10.0
c5=0.0
 
C4=[[[[
    lambada*delta[i][j]*delta[k][l]+mu*(delta[i][k]*delta[j][l]+delta[i][l]*delta[j][k])
    for i in range(3)]
    for j in range(3)]
    for k in range(3)]
    for l in range(3)]
   
C6=[[[[[[
    c1 * ( delta[i][j]*delta[k][l]*delta[m][n] + delta[i][j]*delta[k][m]*delta[l][n] + delta[i][k]*delta[j][n]*delta[l][m] + delta[i][n]*delta[k][j]*delta[l][m] ) +
    c2 * ( delta[i][k]*delta[j][l]*delta[m][n] + delta[i][l]*delta[j][k]*delta[m][n] + delta[i][k]*delta[j][m]*delta[l][n] + delta[i][m]*delta[k][j]*delta[l][n] ) +
    c3 * ( delta[i][l]*delta[j][n]*delta[k][m] + delta[i][n]*delta[k][m]*delta[j][l] + delta[i][m]*delta[k][l]*delta[j][n] + delta[i][n]*delta[k][l]*delta[m][j] ) +
    c4 * ( delta[i][m]*delta[k][n]*delta[j][l] + delta[i][l]*delta[j][m]*delta[k][n] ) +
    c5 *   delta[i][j]*delta[k][n]*delta[l][m]
    for i in range(3)]
    for j in range(3)]
    for k in range(3)]
    for l in range(3)]
    for m in range(3)]
    for n in range(3)]
 
# Generation of equation string
 
eq1 = "("
for i in range(3):
    for j in range(3):
        for k in range(3):
            eq1 = eq1 + ("+v"+"%g%g%g" %(i,j,k) +"*dv"+"%g%g%g" %(i,j,k))
            eq1 = eq1 + ("+u%g.dx(%g)" %(i,j) +"*dv%g%g%g.dx(%g)" %(i,j,k,k))
            for l in range(3):
                if(C4[i][j][k][l]!=0.0):
                    eq1 = eq1 + "+C4[%g][%g][%g][%g]*v%g%g%g*du%g" %(i,j,k,l,k,l,j,i)
                for m in range(3):
                    for n in range(3):
                        if(C6[i][j][k][l][m][n]!=0.0):
                            eq1 = eq1 + "+C6[%g][%g][%g][%g][%g][%g]*v%g%g%g.dx(%g)*du%g.dx(%g)" %(i,j,k,l,m,n,l,m,n,j,i,k)
 
eq1 = eq1+")*dx"
 
n = FacetNormal(mesh)
 
eq2 = "("
 
for i in range(3):
    for j in range(3):
        for k in range(3):
            eq2 = eq2 + "-u%g.dx(%g)*n[%g]*dv%g%g%g" %(i,j,k,i,j,k)
            for l in range(3):
                for m in range(3):
                    for o in range(3):
                        if(C6[i][j][k][l][m][o]!=0.0):
                            eq2 = eq2 + "-C6[%g][%g][%g][%g][%g][%g]*v%g%g%g.dx(%g)*n[%g]*du%g" %(i,j,k,l,m,o,l,m,o,j,k,i)
 
eq2 = eq2 + ")*ds(1)"

print 'Booom!'
 
a=eval(eq1+"+"+eq2)
 
z = Constant(0.0)
L = z*(du0+du1+du2+dv000+dv001+dv002+dv010+dv011+dv012+dv020+dv021+dv022+dv100+dv101+dv102+dv110+dv111+dv112+dv120+dv121+dv122+dv200+dv201+dv202+dv210+dv211+dv212+dv220+dv221+dv222)*dx
 
w = Function(W)
 
A, b = assemble_system(a, L, bc) 
 
#solve(a == L, w, bc)
