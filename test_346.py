from mshr import Rectangle, generate_mesh
from random import choice
import numpy as np
from dolfin import *


def test(n):
    domain = Rectangle(Point(0, 0), Point(2, 3))
    mesh = generate_mesh(domain, n)

    # Cell-cell connectivity will be defined over facet
    tdim = mesh.topology().dim()
    mesh.init(tdim, tdim-1)
    mesh.init(tdim-1, tdim)

    c2f = mesh.topology()(tdim, tdim-1)
    f2c = mesh.topology()(tdim-1, tdim)

    c2c = lambda c: set(np.hstack([f2c(f) for f in c2f(c)]))

    def cells_in_radius(c0, mesh, cell_connectivity, radius):
        x0 = c0.midpoint()
        c0 = c0.index()
        N = set([c0])              # cells in the patch
        C = cell_connectivity(c0)  # where to look elements of the patch
        C.remove(c0)

        while C:
            c = C.pop()
            x = Cell(mesh, c).midpoint()
            if x.distance(x0) < radius:
                N.add(c)
                new = cell_connectivity(c) - N
                C.update(new)
        return N
     
    f = CellFunction('size_t', mesh, 0)

    # c0 = Cell(mesh, choice(range(mesh.num_cells())))
    # for c in  cells_in_radius(c0, mesh, c2c, 0.4):
    #     f[int(c)] = 1
    # f[c0] = 2
    # plot(f, interactive=True)

    t = Timer('dt')
    c0 = choice(range(mesh.num_cells()))
    c0 = Cell(mesh, c0)
    t.start()
    print '\tPatch size', len(cells_in_radius(c0, mesh, c2c, 0.4))
    dt = t.stop()
    
    return mesh.num_cells(), dt

# ------------------------------------------------------------------------------

if __name__ == '__main__':
    size0, dt0 = -1, -1
    for n in 4, 8, 16, 32, 64, 128, 256, 512:
        size, dt = test(n)

        if size0 > 0:
            rate = ln(dt/dt0)/ln(size/size0)
        else:
            rate = -1
        print size, dt, rate
        size0, dt0 = size, dt
