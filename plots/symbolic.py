import sympy
from sympy import *

x1, x2, dv1, dv2, v1, v2, l1, l2, k1, k2, m1, m2 = symbols('x1 x2 v1 dv2 dv1 v2 l1 l2 k1 k2 m1 m2')
system = [(((-k1 * (x1 - l1)) + (k2 * (x2 - x1 - l2))) / m1),
 (-k2 * (x2 - x1 - l2)) / m2]
stablePoints = solve(system,[x1,x2])

