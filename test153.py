from dolfin import *
import numpy as np

domain = Circle(0., 0., 1) + Rectangle(-1.5, -3, 1.5, 0)
mesh = Mesh(domain, 1)

tol = 1E-14
def inside_arch(x, on_boundary):
    on_circle = x[1] > 0 - DOLFIN_EPS and abs(np.hypot(x[0], x[1]) - 1) < tol
    return on_circle and on_boundary

Arch = AutoSubDomain(inside_arch)

check_midpoint_false = FacetFunction('size_t', mesh, 0)
Arch.mark(check_midpoint_false, 1, False)

check_midpoint_true = FacetFunction('size_t', mesh, 0)
Arch.mark(check_midpoint_true, 1, True)

plot(check_midpoint_false, title='False')
plot(check_midpoint_true, title='True')
interactive()
