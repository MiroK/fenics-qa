from dolfin import *

my_expression ='''
class test : public Expression
{
public:

  test() : Expression() { }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    f->eval(values, x);
  }

  std::shared_ptr<const Function> f; // DOLFIN 1.4.0
  //boost::shared_ptr<const Function> f; // DOLFIN 1.3.0
};
'''

mesh = RectangleMesh(-1, -1, 1, 1, 40, 40)

V = FunctionSpace(mesh, 'CG', 1)

f = interpolate(Expression('sin(2*pi*(x[0]*x[0] + x[1]*x[1]))'), V)

g = Expression(my_expression)
g.f = f

plot(g, mesh=mesh)
interactive()
