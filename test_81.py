import dolfin as dl

mesh = dl.UnitSquareMesh(10,10)
V1 = dl.VectorFunctionSpace(mesh, "CG", 2)
V2 = dl.FunctionSpace(mesh, "CG", 1)
V3 = dl.FunctionSpace(mesh, "CG", 1)
V = dl.MixedFunctionSpace([V1,V2,V3])

f1 = dl.Expression( ("A + B*x[1]*(1-x[1])", "0."), A=1., B=1.)
f2 = dl.Expression( "A*x[0]", A=1.)
f3 = dl.Expression( "A", A = -3.)

F1, F2, F3 = [dl.interpolate(fi, Vi) for fi, Vi in zip((f1, f2, f3), (V1, V2, V3))]

assigner = dl.FunctionAssigner(V, [V1, V2, V3])
state = dl.Function(V)
assigner.assign(state, [F1, F2, F3])

print "|| state ||^2_L^(2) = ", dl.assemble( dl.inner(state, state)*dl.dx)

dl.plot(state.sub(0))
dl.plot(state.sub(1)) 
dl.plot(state.sub(2))

dl.interactive()

