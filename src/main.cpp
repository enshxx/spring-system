#include <iostream>
#include <array>


const int X1 = 0;
const int X2 = 1;
const int V1 = 2;
const int V2 = 3;

double h = 0.001;
int steps = 30000;
double k1, k2, m1, m2, l1, l2;



std::array<double, 4> f(std::array<double, 4> state){
    std::array<double, 4> derivative;
    derivative[X1] = state[V1];
    derivative[X2] = state[V2];
    derivative[V1] = ((-k1 * (state[X1] - l1)) + (k2 * (state[X2] - state[X1] - l2))) / m1;
    derivative[V2] = (-k2 * (state[X2] - state[X1] - l2)) / m2;
    return derivative;
}

void DOPRI8(std::array<double, 4> &array, int steps, double h) {

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

    do
    {
        k[0] = f(array);
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
                k[i] = f(sums);
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

    } while (--steps != 0);
    
}

int main(int argc, char **argv)
{   
    std::array<double, 4> array;
    if (argc < 9)
    {
        fprintf(stderr, "Usage: prog <x1> <x2> <l1> <l2> <k1> <k2> <m1> <m2>\n");
        return -1;
    }
    if (sscanf(argv[1], "%lf", &array[X1]) < 1 || sscanf(argv[2], "%lf", &array[X2]) < 1 || sscanf(argv[3], "%lf", &l1) < 1 || sscanf(argv[4], "%lf", &l2) < 1 || sscanf(argv[5], "%lf", &k1) < 1 || sscanf(argv[6], "%lf", &k2) < 1 || sscanf(argv[7], "%lf", &m1) < 1 || sscanf(argv[8], "%lf", &m2) < 1)
    {
        fprintf(stderr, "Usage: prog <x1> <x2> <l1> <l2> <k1> <k2> <m1> <m2>\n");
        return -1;
    }
    array[V1] = 0;
    array[V2] = 0;
    double l1 = array[X1];
    double l2 = array[X2];
    DOPRI8(array, steps, h);
    for (int i = 0; i < 2; i++) {
        printf("%lf ", array[i]);
    }
    printf("\n");
    return 0;
}
