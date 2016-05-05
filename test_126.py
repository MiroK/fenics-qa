from dolfin import *

mesh = UnitSquareMesh(10, 10)                                                   
V = FunctionSpace(mesh, 'CG', 1)  
u = Function(V)                                                                 

vec = u.vector()
values = vec.get_local()                                            

dofmap = V.dofmap()                                                             
my_first, my_last = dofmap.ownership_range()                # global

# Let's build manually function x**2 + y**2
visited = []

# 'Handle' API change of tabulate coordinates
if dolfin_version().split('.')[1] == '7':
    x = V.tabulate_dof_coordinates().reshape((-1, 2))
else:
    x = V.dofmap().tabulate_all_coordinates(mesh)

for cell in cells(mesh):                                                        
    dofs = dofmap.cell_dofs(cell.index())                   # local                                    
    for dof in dofs:                                                            
        if not dof in visited:
            global_dof = dofmap.local_to_global_index(dof)  # global
            if my_first <= global_dof < my_last:                  
                visited.append(dof)                                                      
                values[dof] = sum(x[dof]**2)

vec.set_local(values)
vec.apply('insert')

# Check
u0 = interpolate(Expression('x[0]*x[0]+x[1]*x[1]', degree=2), V)
u0.vector().axpy(-1, vec)

error = u0.vector().norm('linf')
if MPI.rank(mpi_comm_world()) == 0: print error
