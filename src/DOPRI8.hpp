#ifndef     DOPRI8_H
#define     DOPRI8_H

#include "function.hpp"


void DOPRI8(
    std::array<double, 4> &array,
    Params const& p,
    int steps,
    double h,
    std::vector<std::array<double, 4> > &stepPoints) {
    stepPoints.resize(steps);

    // Таблица Бутчера
    static double a[14][13] = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)1 / 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)1 / 48, (double)1 / 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)1 / 32, 0, (double)3 / 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)5 / 16, 0, (double)-75 / 64, (double)75 / 64, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)3 / 80, 0, 0, (double)3 / 16, (double)3 / 20, 0, 0, 0, 0, 0, 0, 0, 0},
        {(double)29443841 / 614563906, 0, 0, (double)77736538 / 692538347, (double)-28693883 / 1125000000, (double)23124283 / 1800000000, 0, 0, 0, 0, 0, 0, 0},
        {(double)16016141 / 946692911, 0, 0, (double)61564180 / 158732637, (double)22789713 / 633445777, (double)545815736 / 2771057229, (double)-180193667 / 1043307555, 0, 0, 0, 0, 0, 0},
        {(double)39632708 / 573591083, 0, 0, (double)-433636366 / 683701615, (double)-421739975 / 2616292301, (double)100302831 / 723423059, (double)790204164 / 839813087, (double)800635310 / 3783071287, 0, 0, 0, 0, 0},
        {(double)246121993 / 1340847787, 0, 0, (double)-37695042795 / 15268766246, (double)-309121744 / 1061227803, (double)-12992083 / 490766935, (double)6005943493 / 2108947869, (double)393006217 / 1396673457, (double)123872331 / 1001029789, 0, 0, 0, 0},
        {(double)-1028468189 / 846180014, 0, 0, (double)8478235783 / 508512852, (double)1311729495 / 1432422823, (double)-10304129995 / 1701304382, (double)-48777925059 / 3047939560, (double)15336726248 / 1032824649, (double)-45442868181 / 3398467696, (double)3065993473 / 597172653, 0, 0, 0},
        {(double)185892177 / 718116043, 0, 0, (double)-3185094517 / 667107341, (double)-477755414 / 1098053517, (double)-703635378 / 230739211, (double)5731566787 / 1027545527, (double)5232866602 / 850066563, (double)-4093664535 / 808688257, (double)3962137247 / 1805957418, (double)65686358 / 487910083, 0, 0},
        {(double)403863854 / 491063109, 0, 0, (double)-5068492393 / 434740067, (double)-411421997 / 543043805, (double)652783627 / 914296604, (double)11173962825 / 925320556, (double)-13158990841 / 6184727034, (double)3936647629 / 1978049680, (double)-160528059 / 685178525, (double)248638103 / 1413531060, 0, 0}
    };

    double b_coefs[] = {(double)14005451 / 335480064, 0, 0, 0, 0, (double)-59238493 / 1068277825, (double)181606767 / 758867731, (double)561292985 / 797845732, (double)-1041891430 / 1371343529, (double)760417239 / 1151165299, (double)118820643 / 751138087, (double)-528747749 / 2220607170, (double)1 / 4};

    std::array<std::array<double, 4>, 13> k = {};
    // setting initial data
    array[X1] = p[X10];
    array[X2] = p[X20];
    array[V1] = p[V10];
    array[V2] = p[V20];

    do
    {
        k[0] = f(array, p);
        std::array<double, 4> sums = {0};

        for (int i = 1; i <= 12; i++)
        {
            sums.fill(0.0);

            for (int j = 0; j < i; j++)
            {
                for (int ind = 0; ind < 4; ++ind) { 
                    sums[ind] += a[i][j] * k[ind][j];
                }
            }
            for (int ind = 0; ind < 4; ++ind) { 
                sums[ind] = array[ind] + sums[ind] * h;
            }

            for (int ind = 0; ind < 4; ++ind) { 
                k[i] = f(sums, p);
            }
        }

        sums.fill(0.0);

        for (int i = 0; i <= 12; i++)
        {
            for (int ind = 0; ind < 4; ++ind) { 
                sums[ind] += b_coefs[i] * k[ind][i];
            }
        }
        for (int ind = 0; ind < 4; ++ind) { 
            array[ind] += sums[ind] * h;
        }

        for (int i = 0; i < 4; ++i) {
            stepPoints[stepPoints.size()-steps][i] = array[i];
        }

    } while (--steps != 0);
    
}

#endif