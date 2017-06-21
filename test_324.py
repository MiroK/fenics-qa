from dolfin import *
from petsc4py import PETSc
import numpy as np


def diff_coeff_matrix(V, points):
    mesh = V.mesh()

    element = V.element()
    basis_values = np.zeros(element.space_dimension())

    nrows = len(points)
    ncols = V.dim()
    mat = PETSc.Mat().createAIJ(size=(nrows, ncols))
    mat.setUp() 
    mat.assemblyBegin()

    dmap = V.dofmap()
    tree = mesh.bounding_box_tree()
    for row, pt in enumerate(points):
        cell_id = tree.compute_first_entity_collision(Point(*pt))

        # Point missing
        if cell_id >= mesh.num_cells(): continue

        cell = Cell(mesh, cell_id)
        vertex_coordinates = cell.get_vertex_coordinates()
        cell_orientation = cell.orientation()

        element.evaluate_basis_all(basis_values, pt, vertex_coordinates, cell_orientation)

        row_indices = [row]
        col_indices = dmap.cell_dofs(cell_id)

        mat.setValues(row_indices, col_indices, basis_values, PETSc.InsertMode.INSERT_VALUES)
    mat.assemblyEnd()
    return PETScMatrix(mat)
    
# ----------------------------------------------------------------------------------------

if __name__ == '__main__':
    mesh = UnitSquareMesh(3, 3)
    V = FunctionSpace(mesh, 'CG', 8)

    # The matrix should have an identity property if points are
    # dof coordinates
    points = V.tabulate_dof_coordinates().reshape((V.dim(), -1))
    # Let's test it 
    I = diff_coeff_matrix(V, points)
    assert I.size(0) == I.size(1)

    x = as_backend_type(I).mat().createVecRight()
    y = x.copy()
    x.setRandom()

    x, y = map(PETScVector, (x, y))
   
    print x.norm('l2'), y.norm('l2')
    I.mult(x, y)
    print (x-y).norm('l2')
