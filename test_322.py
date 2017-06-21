import dolfin as dl

mesh = dl.UnitCubeMesh(1,1,1)

vectorP2Element = dl.VectorElement("P", mesh.ufl_cell(), 2) 
scalarP1Element = dl.FiniteElement("P", mesh.ufl_cell(), 1)
TaylorHoodV  = dl.FunctionSpace(mesh, dl.MixedElement([vectorP2Element,scalarP1Element]))
parameterV = dl.VectorFunctionSpace(mesh,"R",0,7)  

y = dl.Function(TaylorHoodV)
ytest = dl.TestFunction(TaylorHoodV)
ytrial = dl.TrialFunction(TaylorHoodV)
u,p = dl.split(y)

a = dl.Function(parameterV)
atest = dl.TestFunction(parameterV)
atrial = dl.TrialFunction(parameterV)
#Model Parameters
m0 = 543.09
m1 = 4.73
m2 = 4.57
m3 = 5.27
m4 = 4.49
m5 = 4.49
m6 = 12.02 
m = (m0,m1,m2,m3,m4,m5,m6)
parameters = dl.interpolate(dl.Constant(m),parameterV)

#Model
d    = len(u)
I    = dl.Identity(d)            # Identity tensor
F    = I + dl.grad(u)            # Deformation gradient
C    = F.T*F                     # Right Cauchy-Green tensor
J    = dl.det(F)                 # Jacobian of deformation gradient
Cbar = C*(J**(-2.0/3.0))           # Deviatoric Right Cauchy-Green tensor
E    = dl.Constant(0.5)*(Cbar-I) # Deviatoric Green strain

c0,c1,c2,c3,c4,c5,c6 = dl.split(parameters)

Q = c1*(E[0,0]**2)+c2*(E[1,1]**2)+c3*(E[2,2]**2) + c4*(dl.Constant(4.0)*(E[1,2]**2)) #+c5*(dl.Constant(4.0)*(E[0,1]**2))#)#+c6*(dl.Constant(4.0)*(E[0,2]**2))
Energy = dl.Constant(0.5)*c0*(dl.exp(Q)-dl.Constant(1.0))
IncompressibilityConstraint = (J-dl.Constant(1.0))
Lagrangian = p*IncompressibilityConstraint*dl.dx  + Energy*dl.dx
yhat = dl.interpolate(dl.Constant((1,1,1,1)),TaylorHoodV)
residual_form = dl.derivative(Lagrangian,y,yhat) 
residual_form_a = dl.derivative(residual_form, a,atest)
residual_form_y = dl.derivative(residual_form, y,ytest)

residual_form_ya = dl.assemble(dl.derivative(residual_form_y, a,atrial))
residual_form_ay = dl.assemble(dl.derivative(residual_form_a, y,ytrial))
residual_form_aa = dl.assemble(dl.derivative(residual_form_a, a,atrial))
residual_form_yy = dl.assemble(dl.derivative(residual_form_y, y,ytrial))
