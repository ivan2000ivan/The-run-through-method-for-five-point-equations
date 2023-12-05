#ifndef MATRIXSOLVER_H
#define MATRIXSOLVER_H

#include <vector>

class MatrixSolver
{
public:
    MatrixSolver();
    static void solve5PointRun(const std::vector<double> &a, const std::vector<double> &b,
                               const std::vector<double> &c, const std::vector<double> &d,
                               const std::vector<double> &e, const std::vector<double> &f,
                               std::vector<double> &x);
};

#endif // MATRIXSOLVER_H
