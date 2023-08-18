
#ifndef SOLVER_H
#define SOLVER_H

#include <functional>
#include <Eigen/Dense>

namespace Solver {

    using Eigen::VectorXd;
    using NonLinearSystem = std::function<VectorXd(VectorXd const&)>;

    enum class NRSolveStatus {
        CONVERGED,
        NON_CONVERGED
    };

    NRSolveStatus NR_solve(NonLinearSystem const& F, VectorXd const& rhs,
                      VectorXd& guess, double epsilon, 
                      int max_iterations, bool verbose = false);

}

#endif // SOLVER_H
