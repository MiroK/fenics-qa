from dolfin import *
from sympy.utilities.codegen import ccode
from sympy import symbols
import sympy as sp

x, y, t, w = symbols('x[0], x[1], t, w')

# Sympy expression
r = sp.sqrt((x-2)**2 + (y-2)**2)
f = sp.cos(w*sp.pi*r*t)/(sp.pi*w*r)
# Derivative
df = f.diff(t, 1)
print df
# Convert to Dolfin expression. Note that args t, w are given values
df = Expression(ccode(df).replace('M_PI', 'pi'), t=1, w=1)

# Demo
mesh = RectangleMesh(Point(-1, -1), Point(1, 1), 20, 20)
V = FunctionSpace(mesh, 'CG', 1)
dfV = interpolate(df, V)
p = plot(dfV, elevate=0., rescale=False)

dt=1E-1
t = 1
while t < 5:
    t += dt
    df.t = t
    dfV = interpolate(df, V)
    p.plot(dfV)
interactive()
