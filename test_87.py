from dolfin import *
from mshr import *

ll0 = Point(0, 0)
ur0 = Point(-3, -3)

ll1 = Point(-1, -1)
ur1 = Point(-2, -2)

domain0 = Rectangle(ur0, ll0)
domain1 = Rectangle(ur1, ll1)
domain = domain0 + domain1
#domain.set_subdomain(1, domain0)
#domain.set_subdomain(2, domain1)

mesh_OK = generate_mesh(domain, 16)
plot(mesh_OK)


interactive()
