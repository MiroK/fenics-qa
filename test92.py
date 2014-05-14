from dolfin import *
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, cg, bicg, bicgstab, minres, gmres

parameters.linear_algebra_backend = "uBLAS" 

def test(mesh, normalize=False, scipy_solver='lu'):
  V = FunctionSpace(mesh, 'Lagrange', 1)
  bc= DirichletBC(V, Constant(0.0), DomainBoundary())

  u  =TrialFunction(V)
  v = TestFunction(V)

  f = Constant(1.0)
  a = inner(grad(u), grad(v))*dx
  L = f*v*dx

  A, b = assemble_system(a, L, bc)

  if normalize:
    max = A.array().max()
    A /= max
    b /= max

  u = Function(V)
  solve(A, u.vector(), b)
  plot(u, interactive=True, title='dolfin')

  x = A.array()
  print "Value of A are [%g, %g]" %(x.min(), x.max())
  print "Num of entries in A larger theb 100", np.where(x > 1E2)[0].sum()
  print "Number of mesh cells", mesh.num_cells()

  # see about scipy
  rows, cols, values = A.data()
  B = csr_matrix((values, cols, rows))
  d = b.array().T

  if scipy_solver == 'cg':
    v_, info = cg(B, d)              
  elif scipy_solver == 'bicg':
    v_, info = bicg(B, d)              
  elif scipy_solver == 'bicgstab':
    v_, info = bicgstab(B, d)              
  elif scipy_solver == 'gmres':
    v_, info = gmres(B, d)              
  elif scipy_solver == 'minres':
    v_, info = minres(B, d)              
  else:
    v_ = spsolve(B, d)              

  v = Function(V)
  v.vector()[:] = v_
  plot(v, interactive=True, title='scipy')
  try:
    print "info", info, 'v_max', v_.max()
  except:
    pass

#------------------------------------------------------------------------------

def get_jacobian(mesh):
  V = FunctionSpace(mesh, 'DG', 0)
  j = Function(V)
  dofmap = V.dofmap()
  
  cell_f = CellFunction('bool', mesh, False)

  for cell in cells(mesh):
    j.vector()[dofmap.cell_dofs(cell.index())] = cell.volume()

  min_ = j.vector().min()
  max_ = j.vector().max()

  for cell in cells(mesh):
    if near(cell.volume(), min_):
      cell_f[cell] = True

  #plot(j, interactive=True)
  #plot(cell_f, interactive=True)

  return min_, max_

#------------------------------------------------------------------------------

domain = Circle(0, 0, 10)
mesh0 = Mesh(domain, 10)

mesh1 = UnitSquareMesh(13, 13)
mesh2 = CircleMesh(Point(0., 0.), 1, 0.1)

names = ['circle', 'square', 'unitcircle']
meshes = [mesh0, mesh1, mesh2]

scipy_solvers = ['lu', 'cg', 'bicg', 'bicgstab', 'gmres', 'minres']

for name, mesh in zip(names, meshes):
  print name, get_jacobian(mesh)
  #for scipy_solver in scipy_solvers[0:1]:
  #    print name, scipy_solver
  #    test(mesh, normalize=False, scipy_solver=scipy_solver)
  #    print



