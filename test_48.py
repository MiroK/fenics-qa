from dolfin import *
from mshr import *

# Define Domain

diam    = 5.0
mainBox = Rectangle(Point(0,0),Point(diam,diam))

geo = mainBox

# Prepare domain labels

inletRadius     = 0.5
INmid       = Point(1.25,diam)
OUTmid      = Point(3.75,0.0)

def dist(mid,x):
    dx = mid[0] - x[0]
    dy = mid[1] - x[1]
    return sqrt(dx*dx + dy*dy)

def inlet(x):   
    if ( near(x[1],diam) and  dist(INmid,x) < inletRadius+DOLFIN_EPS ): 
        return True
    else:
        return False

def outlet(x):
    if ( near(x[1],0.0) and  dist(OUTmid,x) < inletRadius+DOLFIN_EPS ): 
        return True
    else:
        return False


def noslip(x):
    if ( not outlet(x) and not inlet(x) ) \
    and (   near(x[0],0.0) or near(x[0],diam)\
        or near(x[1],0.0) or near(x[1],diam) ):
        return True
    else:
        return False

# Build mesh

mesh3d = generate_mesh(geo, 10)

# Boundary description
bndry   = FacetFunction("size_t", mesh3d)

noV     = AutoSubDomain(lambda x,bndry: bndry and noslip(x) )
inletV      = AutoSubDomain(lambda x,bndry: bndry and inlet(x) )
outletV     = AutoSubDomain(lambda x,bndry: bndry and outlet(x) )

# noV     .mark(bndry,1)
inletV      .mark(bndry,2)
outletV     .mark(bndry,3)

plot(bndry,interactive=True)
