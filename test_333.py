from dolfin import *
import subprocess

domain = '''
SIZE = 0.01;
Point(1) = {0, 0, 0, SIZE};
Point(2) = {1, 0, 0, SIZE};
Point(3) = {1, 1, 0, SIZE};
Point(4) = {0, 1, 0, SIZE};
// 2d domain
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Cracks
// 1
Point(5) = {0.25, 0.25, 0, SIZE};
Point(6) = {0.5, 0.5, 0, SIZE};
Line(7) = {5, 6};
Physical Line(1) = {7};
Line{7} In Surface{6};

// 2
Point(7) = {0.6, 0.7, 0, SIZE};
Line(8) = {6, 7};
Physical Line(2) = {8};
Line{8} In Surface{6};

Point(8) = {0.8, 0.2, 0, SIZE};
Line(9) = {6, 8};
Physical Line(3) = {9};
Line{9} In Surface{6};

Physical Surface(1) = {6};
'''

with open('domain.geo', 'w') as f: f.write(domain)

subprocess.call(['gmsh -2 domain.geo'], shell=True)
subprocess.call(['dolfin-convert domain.msh domain.xml'], shell=True)

mesh = Mesh('domain.xml')
facet_f = MeshFunction('size_t', mesh, 'domain_facet_region.xml')

for tag in (1, 2, 3):
    print sum(1 for _ in SubsetIterator(facet_f, tag)), 'edges marked as', tag

plot(facet_f)
interactive()
