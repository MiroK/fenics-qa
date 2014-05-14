from dolfin import *
import random

code = '''
#include <boost/math/special_functions/bessel.hpp>
using boost::math::cyl_bessel_i;
using boost::math::cyl_bessel_j;
using boost::math::cyl_bessel_k;
using boost::math::cyl_neumann;

namespace dolfin {
    class MyFun : public Expression
    {
        public:
            MyFun(): Expression() {};
        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = cyl_bessel_j(0,x[0]);
        }
    };
}'''

f=Expression(code)
for N in [3, 4, 5, 6, 7]:
	t = Timer('bessel')
	t.start()
	for x in range(int(10**N)):
		f(random.random())
	t.stop()
	print N, timing('bessel')
