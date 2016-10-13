import dolfin
import numpy as np

mesh = dolfin.UnitCubeMesh(1, 1, 1)

dx = dolfin.Measure("dx", domain=mesh)

code2 = '''

namespace dolfin
{
  
  std::vector<std::size_t> shape = {3,3,3,3};

class MyFunc2 : public Expression
{

  mutable Array<double> UX;

public:

  std::shared_ptr<MeshFunction<std::size_t> > cell_data;
  std::shared_ptr<const Function> U;

  MyFunc2() : Expression(shape), UX(2)
  {
  }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    for(int i=0; i<81;i++){
        values[i] = i;
    }
  }
};
}'''


Tensor4 = dolfin.Expression(code2, degree=1)
assert Tensor4.ufl_shape == (3, 3, 3, 3)
assert max(np.abs(Tensor4(0, 0, 0) - np.arange(81))) < 1E-14

from dolfin import *

I = Identity(3)
foo = Tensor4[i, j, k, l]*I[k, l]
