from dolfin import *

class MyExpression(Expression):
    def __init__(self, a=0, b=0):
        self.a = a
        self.b = b

    def eval(self, value, x):
        dx = x[0] - self.a
        dy = x[1] - self.b
        value[0] = exp(-(dx*dx + dy*dy)/0.02)
        value[1] = exp(-(dx*dx + dy*dy)/0.02)

    def value_shape(self):
        return (2,)

cpp_code = '''
class MyCppExpression : public Expression
{
public:
  MyCppExpression() : Expression(2), a(0), b(0) { }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    const double dx = x[0] - a;
    const double dy = x[1] - b;
    values[0] = exp(-(dx*dx + dy*dy)/0.02);
    values[1] = exp(-(dx*dx + dy*dy)/0.02);
  }
public:
    double a, b;
};'''

mesh = UnitSquareMesh(20, 20)

f = MyExpression(a=0.5, b=0.5)
f.a = 0.5
f.b = 0.25

f_cpp = Expression(cpp_code)
f_cpp.a = 0.5
f_cpp.b = 0.25

# Compare the two
V = VectorFunctionSpace(mesh, 'CG', 1)
F = interpolate(f, V).vector()
F_cpp = interpolate(f_cpp, V).vector()
F.axpy(-1, F_cpp)
print '|F-F_cpp|', F.norm('l2')
