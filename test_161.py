code = '''
class MyFoo : public Expression
{
public:

  std::shared_ptr<Array<double> >* foo;

  MyFoo() : Expression() { }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    //assert(foo);
    //std::cout << "foo = " << foo->str(1) << std::endl;
    values[0] = 0.0;
  }
};'''

from dolfin import *
import numpy as np

mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 1)

my_foo = dolfin.Expression(cppcode=code, element=V.ufl_element())

foo = np.random.rand(3)
my_foo.foo = foo

values = np.array([100.])
x = np.array([0., 0.3])
my_foo.eval(values, x)

print values
