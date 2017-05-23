from dolfin import *                                                               

mesh = UnitIntervalMesh(3)                                                         
w = FunctionSpace(mesh, 'P', 1)                                                    
u = TrialFunction(w)                                                              
v = TestFunction(w)                                                                

def gamma(x): return x

x = MeshCoordinates(mesh)                                                          
f = Constant(1.0, cell = mesh.ufl_cell()) * gamma(x[0])                               
R = assemble(f*v.dx(0)*dx(domain = mesh))                                  

print(R.array())   
