from dolfin import *

mesh = Mesh('testing.xml.gz')
mvc = MeshValueCollection('size_t', mesh, 'testing_facets.xml.gz')
d = mvc.dim()

default_value = 0
mf = MeshFunction('size_t', mesh, d, default_value)
values = mf.array()

D = mesh.topology().dim()
mesh.init(D, d)
connectivity = mesh.topology()(D, d)

mvc_values = mvc.values()
for ci_lei, value in mvc_values.iteritems():
    cell_index, local_entity_index = ci_lei
    entity_index = connectivity(cell_index)[local_entity_index]
    values[entity_index] = value

plot(mf, interactive=True)
plot(mf, range_min=1., range_max=3., interactive=True)
