from dolfin import *

width, height = 1, 0.5

mesh = RectangleMesh(0, 0, width, height, 8, 4, "right")

transverse_order, axial_order = 2, 1

V_N = FunctionSpace(mesh, "Nedelec 1st kind H(curl)", transverse_order)
V_L = FunctionSpace(mesh, "Lagrange", axial_order)

combined_space = V_N * V_L

(N_i, L_i) = TestFunctions(combined_space)
(N_j, L_j) = TrialFunctions(combined_space)

class HalfLoadedDielectric(Expression):
    def eval(self, values, x):
        if x[1] < 0.25:
            values[0] = 4.0
        else:
            values[0] = 1.0;

e_r = HalfLoadedDielectric()
one_over_u_r = Expression("1.0") # no magnetic material

k_o_squared = Expression("4*pi*pi*nu*nu*1e-4/9", nu=0.0) # nu is given in MHz

def curl_t(w):
    return Dx(w[1],0)-Dx(w[0],1)

s_tt = one_over_u_r*dot(curl_t(N_i), curl_t(N_j))
t_tt = e_r*dot(N_i, N_j)

s_zz = one_over_u_r*dot(grad(L_i), grad(L_j))
t_zz = e_r*L_i*L_j

b_tt = one_over_u_r*dot(N_i, N_j)
b_tz = one_over_u_r*dot(N_i, grad(L_j))
b_zt = one_over_u_r*dot(grad(L_i), N_j)

a_tt = s_tt - k_o_squared*t_tt
b_zz = s_zz - k_o_squared*t_zz

class ElectricWalls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

zero = Constant((0.0, 0.0, 0.0)) # a vector constant value

dirichlet_bc = DirichletBC(combined_space, zero, ElectricWalls())

# Now assemble the matrices
TO DO
# Form the system
TO DO
# Apply the dirichlet_bc
TO DO
# plot the tranverse and longitudinal part

f = Function(combined_space, e)
(transverse, axial) = f.split()
plot(transverse)
plot(axial)
interactive()
