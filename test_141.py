import dolfin

mesh = dolfin.UnitSquareMesh(1,1)
dX = dolfin.dx(mesh)

fe = dolfin.FiniteElement(
    family="Quadrature",
    cell=mesh.ufl_cell(),
    degree=1,
    quad_scheme="default")

cppExprCode='''
namespace dolfin
{

class CppExpr : public Expression
{
public:

    CppExpr(): Expression(0)
    {
    }

    void eval(Array<double>& values, const Array<double>& position) const
    {
        std::cout << "position = " << position << std::endl;
        values[0] = 1.;
        std::cout << "values = " << values << std::endl;
    }
};

}'''
cppExpr = dolfin.Expression(cppExprCode, element=fe)
dolfin.assemble(cppExpr * dX)
