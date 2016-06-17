import ufl
import dolfin.cpp as cpp
import numpy

# Following https://github.com/FEniCS/dolfin/blob/master/site-packages/dolfin/functions/constant.py

class MyTensorConstant(ufl.Coefficient, cpp.Constant):

    def __init__(self, value, cell, name=None, symmetry=None):

        cell = ufl.as_cell(cell)
        ufl_domain = None

        array = numpy.array(value)
        rank = len(array.shape)
        floats = list(map(float, array.flat))

        assert rank == 2

        ufl_element = ufl.TensorElement("Real", cell, 0, shape=array.shape,
                                        symmetry=symmetry)
        cpp.Constant.__init__(self, list(array.shape), floats)

        ufl_function_space = ufl.FunctionSpace(ufl_domain, ufl_element)
        ufl.Coefficient.__init__(self, ufl_function_space, count=self.id())

        name = name or "f_%d" % self.count()
        self.rename(name, "a Constant")

    def cell(self):
        return self.ufl_element().cell()

    def __float__(self):
        # Overriding UFL operator in this particular case.
        if self.ufl_shape:
            raise TypeError("Cannot convert nonscalar constant to float.")
        return cpp.Constant.__float__(self)

    def __str__(self):
        "x.__str__() <==> print(x)"
        return self.name()

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    
    mesh = UnitSquareMesh(10, 10)
    V = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    c_values = ((1, 0), (0, 2))
    c1 = MyTensorConstant(c_values, cell=mesh.ufl_cell())
    c2 = MyTensorConstant(c_values, cell=mesh.ufl_cell(), symmetry={(0, 1): (1, 0)})

    a1 = inner(dot(c1, grad(u)), grad(v))*dx
    a2 = inner(dot(c2, grad(u)), grad(v))*dx

    print a1.signature() == a2.signature()
