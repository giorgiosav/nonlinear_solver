
#include <cassert>
#include <iostream>
#include "solver.h"

namespace Solver {
    using std::cout;

    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    namespace {
        MatrixXd jacobian(int ncols, int nrows, NonLinearSystem const& F, VectorXd const& x,
                          double epsilon = 1e-6) {
            assert(ncols == x.size());
            MatrixXd J(nrows, ncols);
            const double epsilon2 = 2 * epsilon;

            VectorXd xe = x;
            for (int i = 0; i < ncols; ++i) {
                double xei_prev = xe(i);

                xe(i) += epsilon;
                VectorXd F_plus_epsilon = F(xe);

                xe(i) = xei_prev - epsilon;
                VectorXd F_minus_epsilon = F(xe);

                J.col(i) = (F_plus_epsilon - F_minus_epsilon) / epsilon2;
            }
            return J;
        }
    }

    NRSolveStatus NR_solve(NonLinearSystem const& F, VectorXd const& rhs,
                      VectorXd& x0, double epsilon, 
                      int max_iterations, bool verbose) {


        VectorXd& x = x0;

        if (verbose) {
            cout << "\tx = " << x.transpose() << '\n';
            cout << "\trhs = " << rhs.transpose() << '\n';
        }
        for (int i = 0; i < max_iterations; ++i) {
            VectorXd Fx = F(x);
            if (verbose)
                cout << "\tF(x) = " << Fx.transpose() << '\n';
            MatrixXd J = jacobian(x.size(), Fx.size(), F, x);
            // cout << "J:\n = " << J << '\n';
            VectorXd delta_x = J.colPivHouseholderQr().solve(-Fx);

            x += delta_x;
            // double norm = delta_x.norm();
            double norm = (Fx - rhs).norm();
            
            if (verbose) {
                cout << i << ": norm = " << norm << '\n';
                cout << "\tx = " << x.transpose() << '\n';
                cout << "\tdelta_x = " << delta_x.transpose() << '\n';
            }

            if (norm <= epsilon) {
                if (verbose) {
                    cout << "Converged after " << i + 1 << " iterations.\n";
                    cout << "Solution is: " << x.transpose() << '\n';
                    cout << "F(x) = " << F(x).transpose() << '\n';
                    cout << "rhs = " << rhs.transpose() << '\n';
                }
                return NRSolveStatus::CONVERGED;
            }
        }
        
        if (verbose) cout << "NR did not converge\n";
        return NRSolveStatus::NON_CONVERGED;

    }
}