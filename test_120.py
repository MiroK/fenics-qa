from dolfin import *

# This is to check the values i'm using are OK
from scipy import special as sp
print sp.kv( 0.5, 0.0002 )

# The FunctionSpace seems to be irrelevant
V = FunctionSpace(  UnitIntervalMesh( 22 ), "CG", 1 )

# If I comment this line out I get a warning message from PETsC! -------------
fund = Function( V )

# Read code from file.
c_file = open( "tmp.cpp" , 'r' )  

# We can compile the expression without a problem
xp = Expression( c_file.read() )
