#include <dolfin.h>
#include "moebius.h"

using namespace dolfin;

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return on_boundary; }
};

int main()
{
  Mesh mesh("moebius.xml.gz");
  moebius::FunctionSpace V(mesh);

  Constant u0(0.0);
  DirichletBoundary boundary;
  DirichletBC bc(V, u0, boundary);

  moebius::BilinearForm a(V, V);
  moebius::LinearForm L(V);

  std::cout << V.element()->geometric_dimension() << std::endl;
  std::cout << V.element()->topological_dimension() << std::endl;

  Constant f(1.0);
  L.f = f;

  Function u(V);
  solve(a == L, u, bc);
  
  plot(u);
  interactive();
  return 0;
}
