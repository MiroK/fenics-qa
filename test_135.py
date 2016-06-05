from dolfin import *

mesh = UnitSquareMesh(10, 10)
# All is tagged as 'red' first
domains = CellFunction('size_t', mesh, 5)
# Let some shape mark its domain
uh = CompiledSubDomain('x[1] > 0.5')
uh.mark(domains, 1)

plot(domains)
interactive()
