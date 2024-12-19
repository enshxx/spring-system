import sympy
from sympy import *

# https://math.stackexchange.com/questions/3470910/how-to-fit-ordinary-differential-equations-to-empirical-data
# theta = { k1, k2, m1, m2, l1, l2} = { theta5, ..., theta10}
# x = {x_1, x_2, x_3, x_4} = { v1, x1, v2, x2}
# h1 = x_1, h2 = x_2, h3 = x_3, h4 = x_4
#  x_1(0) = theta1, x_2(0) = theta2, x_3(0) = theta3, x_4(0) = theta4
# {y1, y2, y3, y4} = {v1, x1, v2, x2}
x10, x20, l1, l2, k1, k2, m1, m2 = symbols('x10, x20, l1 l2 k1 k2 m1 m2', positive=True)
dt, x1m, x2m, v1m, v2m, x1, x2, v1, v2, x1n_p, x1n, x2n_p, x2n, v1n_p_1, v1n, v2n_p_1, v2n, v10, v20 = \
  symbols('dt, x1m, x2m, v1m, v2m, x1, x2, v1, v2, x1n_p, x1n, x2n_p, x2n, v1n_p_1, v1n, v2n_p_1, v2n, v10, v20')
dv1_dt = (((-k1 * (x1 - l1)) + (k2 * (x2 - x1 - l2))) / m1)
dv2_dt = (-k2 * (x2 - x1 - l2)) / m2

Params = Matrix([v10, x10, v20, x20, k1, k2, l1, l2, m1, m2])
x_1, x_2, x_3, x_4 = symbols('x_1, x_2, x_3, x_4')
State = Matrix([v1, x1, v2, x2])
StateX = Matrix([x_1, x_2, x_3, x_4])
theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10 = \
    symbols('theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10')
Theta = Matrix([theta1,theta2, theta3, theta4, theta5, theta6, theta7, theta8, theta9, theta10])
subs_state_dict = { State[i]:StateX[i] for i in range(State.shape[0]) } 
subs_params_dict = { Params[i]:Theta[i] for i in range(Params.shape[0]) }
F = Matrix([ \
    dv1_dt.subs(subs_state_dict).subs(subs_params_dict),\
    x_1,\
    dv2_dt.subs(subs_state_dict).subs(subs_params_dict),\
    x_3])
FdiffState = F.jacobian(StateX)
FdiffTheta = F.jacobian(Theta)
H = Matrix( [x_1, x_2, x_3, x_4] )
Y = H
YdiffX = H.jacobian(StateX)
# equations of dynamic system
# Xdot = F
# Y = H
YdiffTheta = H.jacobian(Theta)
# components of dX/dtheta

sdt11, sdt12, sdt13, sdt14, sdt15, sdt16, sdt17, sdt18, sdt19, sdt110 \
    = symbols('sdt11, sdt12, sdt13, sdt14, sdt15, sdt16, sdt17, sdt18, sdt19, sdt110')
sdt21, sdt22, sdt23, sdt24, sdt25, sdt26, sdt27, sdt28, sdt29, sdt210 \
    = symbols('sdt21, sdt22, sdt23, sdt24, sdt25, sdt26, sdt27, sdt28, sdt29, sdt210')
sdt31, sdt32, sdt33, sdt34, sdt35, sdt36, sdt37, sdt38, sdt39, sdt310 \
    = symbols('sdt31, sdt32, sdt33, sdt34, sdt35, sdt36, sdt37, sdt38, sdt39, sdt310')
sdt41, sdt42, sdt43, sdt44, sdt45, sdt46, sdt47, sdt48, sdt49, sdt410 \
    = symbols('sdt41, sdt42, sdt43, sdt44, sdt45, sdt46, sdt47, sdt48, sdt49, sdt410')

StateDiffTheta = Matrix(\
    [[sdt11, sdt12, sdt13, sdt14, sdt15, sdt16, sdt17, sdt18, sdt19, sdt110], \
    [sdt21, sdt22, sdt23, sdt24, sdt25, sdt26, sdt27, sdt28, sdt29, sdt210], \
    [sdt31, sdt32, sdt33, sdt34, sdt35, sdt36, sdt37, sdt38, sdt39, sdt310], \
    [sdt41, sdt42, sdt43, sdt44, sdt45, sdt46, sdt47, sdt48, sdt49, sdt410]] \
)

StateDiffThetaDot = FdiffState*StateDiffTheta + FdiffTheta
# print(StateDiffThetaDot)
# YDiffTheta = YdiffX*StateDiffTheta + 