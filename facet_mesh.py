from dolfin import *
import numpy as np


def facet_mesh(gamma_f, gamma_id=1):
    '''
    Build mesh out of facets marked to gamma_id by facet function gamma_f.
    '''
    # Background mesh
    mesh0 = gamma_f.mesh()
    gdim0 = mesh0.geometry().dim()
    tdim0 = mesh0.topology().dim()

    # Need facet function
    tdim = gamma_f.dim()
    assert tdim == tdim0 - 1

    # Make sure we have cell-vertex connectivity
    mesh0.init(tdim0, 0)
    # And the facet-cell connectivity
    mesh0.init(tdim, tdim0)
    # Finally vertices to build up cells of facet mesh
    mesh0_vertices = mesh0.coordinates().reshape((-1, gdim0))

    vertex_map = {}
    counter = 0
    # Vertices and cells (of vertices in their new numbering) of the facet mesh
    mesh_vertices, mesh_cells = [], []
    # i-th item of the map is the facet index of the mesh0 which is now the cell
    # of mesh
    cell_map = []
    for cell, facet in enumerate(SubsetIterator(gamma_f, gamma_id)):
        cell_map.append(facet.index())

        # Assign to vertex a new id and add it to vertex coordinates if
        # it is not already there
        vertices = facet.entities(0)
        for vertex in vertices:
            if vertex not in vertex_map:
                vertex_map[vertex] = counter
                counter += 1
                mesh_vertices.append(mesh0_vertices[vertex, :])
        # Vertices making up the cell in new numbering
        mesh_cells.append([vertex_map[vertex] for vertex in vertices])
   
    mesh_vertices = np.array(mesh_vertices)

    # Must have at least a cell...
    assert all(size > 0 for size in map(len, (mesh_cells, mesh_vertices)))

    # Build the mesh
    mesh = Mesh()
    editor = MeshEditor()
    editor.open(mesh, tdim, mesh_vertices.shape[1])
    editor.init_vertices(len(mesh_vertices))
    editor.init_cells(len(mesh_cells))

    # Add vertices
    for index, v in enumerate(mesh_vertices): editor.add_vertex(index, v)

    # Add cells
    for index, vs in enumerate(mesh_cells): editor.add_cell(index, *vs)

    editor.close()
    mesh.order()

    return mesh, cell_map

# ----------------------------------------------------------------------------

if __name__ == '__main__':

    mesh0 = UnitCubeMesh(10, 10, 10)
    gamma_f = FacetFunction('size_t', mesh0, 0)
    DomainBoundary().mark(gamma_f, 1)
    CompiledSubDomain('near(x[2], 0.)').mark(gamma_f, 0)
    CompiledSubDomain('near(x[2], 1.)').mark(gamma_f, 0)

    mesh, cell_map = facet_mesh(gamma_f)

    # 'Volume' a.k.a surface of the 'cube'
    assert abs(sum(cell.volume() for cell in cells(mesh))-4) < 1E-13

    # Same but involve the function space
    V = FunctionSpace(mesh, 'CG', 1)
    c = interpolate(Constant(1), V)
    assert abs(assemble(c*dx)-4) < 1E-13
