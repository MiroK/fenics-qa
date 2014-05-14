from dolfin import *

mesh = UnitIntervalMesh(2)   # 2d ok, 3d? meshsize?

V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

a1 = u*v*dx

A1 = assemble(a1,
              form_compiler_parameters={'quadrature_rule': 'vertex',
                                        'quadrature_degree': 1}
              )
