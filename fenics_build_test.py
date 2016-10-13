RED = '\033[1;37;31m%s\033[0m'
GREEN = '\033[1;37;32m%s\033[0m'

count = 0
def print_red(s): global count; count = count +1; print RED % s

# Imports
for module in ('dolfin', 'petsc4py', 'slepc4py', 'mpi4py'):
    code = '''
try:
    import %s
except ImportError:
    print_red('Failed to import %s')
    ''' % (module, module)
    exec(code)

# Function space construction
try:
    mesh = dolfin.UnitSquareMesh(10, 10)
    V = dolfin.FunctionSpace(mesh, 'CG', 1)
    status = True
except:
    print_red('Failed to construct function space')
    status = False

# Simple form
if status:
    try:
        u, v = dolfin.TrialFunction(V), dolfin.TestFunction(V)
        A = dolfin.assemble(dolfin.inner(u, v)*dolfin.dx)
    except:
        print_red('Failed to assemble matrix')

# Instant
try:
    f = dolfin.Expression('x[0]+x[1]', degree=1)
except:
    print_red('Failed to compile Expression')

# See if we have petsc LA
try:
    la = dolfin.linear_algebra_backends()
    assert 'PETSc' in la
except AssertionError:
    print_red('No PETSc in LA backends. Only %s are available' % ',' .join(la.keys()))

if count == 0: print GREEN % 'All is well with FEniCS stack'

import sys
sys.exit(0)
