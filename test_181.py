from dolfin import Expression

ccode = '(t-x[0])*(((t-x[0]) < 0) ? 0 : 1)*(((t-x[0]) > tmax) ? 0 : 1)'
foo = Expression(ccode, tmax=20, t=0.0)

for t in (0.1, 0.2, 0.3):
    foo.t = t
    print [foo(x) for x in (-1, -2, -3)]
