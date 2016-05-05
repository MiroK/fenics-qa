import numpy
from dolfin import *

def get_eig(hes):
  mesh = hes.function_space().mesh()
  [eigL,eigR] = numpy.linalg.eig(hes.vector().array().reshape([mesh.num_vertices(),2,2]))
  eig = Function(VectorFunctionSpace(mesh,'CG',1))
  eig.vector().set_local(eigL.flatten())
  return eig

mesh = RectangleMesh(Point(0.,0.), Point(2.,1.),20,10)
u = interpolate(Expression(("cos(x[0])", "sin(x[1]+x[0])")),VectorFunctionSpace(mesh,'CG',2))
hes = project(grad(u),TensorFunctionSpace(mesh,'CG',1))
eig = get_eig(hes)
plot(eig[0])
plot(Expression('-sin(x[0])'), mesh=eig.function_space().mesh())
interactive()
