
#include "solver.h"

using Eigen::VectorXd;
using Eigen::Matrix;

using Vector5d = Matrix<double, 5, 1>;

double f1(double x1, double x2, double x3, double x4, double x5) {
    return 4*x1*x1 + x2*x2 + 3*x3*x3 + 8*x4*x4 + x5*x5;
}

double f2(double x1, double x2, double x3, double x4, double x5) {
    return 2*x1*x2 + 6*x2*x4 + 3*x3*x5 + 3*x1*x5;
}

double f3(double x1, double x2, double x3, double x4, double x5) {
    return x1*x4 + 5*x3*x4 + x2*x3 + 7*x5*x5;
}

double f4(double x1, double x2, double x3, double x4, double x5) {
    return 3*x2*x2 + 4*x5*x4 + x1*x3;
}

double f5(double x1, double x2, double x3, double x4, double x5) {
    return x1*x1 + 8*x3*x4 + 3*x2*x5;
}

Vector5d F(Vector5d const& x) {
    Vector5d y;

    y[0] = f1(x(0), x(1), x(2), x(3), x(4));
    y[1] = f2(x(0), x(1), x(2), x(3), x(4));
    y[2] = f3(x(0), x(1), x(2), x(3), x(4));
    y[3] = f4(x(0), x(1), x(2), x(3), x(4));
    y[4] = f5(x(0), x(1), x(2), x(3), x(4));

    return y;
}

void test_system1() {
    Vector5d rhs;
    rhs << 33, 15, 74, 14, 17;

    VectorXd x0 = Vector5d::Random() * 5;
    Solver::NR_solve(F, rhs, x0, 1e-6, 100, true);
}

int main() {
    test_system1();
}