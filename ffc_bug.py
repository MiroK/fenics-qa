import numpy
from FIAT.polynomial_set import mis
from FIAT.reference_element import default_simplex
from FIAT.quadrature import make_quadrature

order = 1
quadrature = make_quadrature(default_simplex(1), order)

from FIAT.lagrange import Lagrange

degree = 1
element = Lagrange(default_simplex(1), degree)

vertices = [n.get_point_dict().keys()[0] for n in element.dual.get_nodes()]

quadpts = numpy.array(quadrature.get_points(), dtype=numpy.float64)
quadwts = numpy.array(quadrature.get_weights(), dtype=numpy.float64)
numQuadPts = len(quadpts)
evals = element.get_nodal_basis().tabulate(quadrature.get_points(), 1)
basis = numpy.array(evals[mis(1, 0)[0]], dtype=numpy.float64).transpose()
numBasis = element.get_nodal_basis().get_num_members()
basisDeriv = numpy.array([evals[alpha] for alpha in mis(1, 1)], dtype=numpy.float64).transpose()

print "order: %d" % order
print "degree: %d" % degree
print "numQuadPts: %d" % numQuadPts
print "basis:" 
print basis
print "basisDeriv:"
print basisDeriv
