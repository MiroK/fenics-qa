from dolfin import *


mesh = UnitSquareMesh(200, 200)

ellipse = CompiledSubDomain("(x[0] - 0.5)*(x[0] - 0.5)/0.25+\
                             (x[1] - 0.5)*(x[1] - 0.5)/0.125 < 1")

all_cells = CellFunction('size_t', mesh, 0)
ellipse.mark(all_cells, 1)
ellipse_cells = SubsetIterator(all_cells, 1)

plot(all_cells)
interactive()
