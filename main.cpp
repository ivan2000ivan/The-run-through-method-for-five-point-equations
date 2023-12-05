#include <iostream>
#include <vector>
#include <algorithm>
#include "matrixsolver.h"
#include <fstream>

// Тест для алгоритма немонотонной пятиточечной прогонки
// [Самарский А. А., Николаев Е. С. Методы решения сеточных уравнений. – 1978., стр 102]
void test()
{
    std::vector<double> a(11, 0), b(11, 0), c(11, 0), d(11, 0), e(11, 0), f(11, 0), y(11, 0);
    f.at(0) = 2;

    c.at(0) = 1;
    d.at(0) = 1;
    e.at(0) = 2;

    b.at(1) = 1;
    c.at(1) = 1;
    d.at(1) = 1;
    e.at(1) = 1;

    for (int i = 2; i <= 8; ++i)
    {
        a.at(i) = 1;
        b.at(i) = 1;
        c.at(i) = 2;
        d.at(i) = 1;
        e.at(i) = 1;
    }

    a.at(10) = 2;
    b.at(10) = 1;
    c.at(10) = 1;

    a.at(9) = 1;
    b.at(9) = 1;
    c.at(9) = 1;
    d.at(9) = 1;
    MatrixSolver::solve5PointRun(a, b, c, d, e, f, y);
    std::cout << "Test result: " << std::endl;
    for (auto val : y) std::cout << val << std::endl;
}

int main()
{
    const double L = 0.8; // [m]
    const double U = 0.008; // [m / s]
    const double H1 = 0.04, H2 = 0.02; // [m]
    const double E = 2E+11; // [Па]
    const double delta = 0.01; // [м]
    const double I = delta * delta * delta / 12; // [м^3]
    const double mu = 0.001; // [Па * с]
    auto Q = [](double U, double H) {return U * H;};
    const double Q1 = Q(U, H1), Q2 = Q(U, H2);
    auto g1 = [&](double x)
    {
        return 6 * mu * (Q1 / std::pow(H1, 2) + Q2 / std::pow(H2, 2)) * (L - x);
    };
    auto g2 = [&](double x)
    {
        return 36 * mu * (Q1 / std::pow(H1, 4) + Q2 / std::pow(H2, 4)) * (L - x);
    };
    auto g3 = [&](double x)
    {
        return 12 * mu * (Q1 / std::pow(H1, 3) - Q2 / std::pow(H2, 3)) * (L - x);
    };

    const int N = 103; // N = 1 + N_cells + 1
    const double h = L / (N - 3);
    std::vector<double> x(N), y(N);
    x.at(0) = - h;
    for (int i = 1; i < N; ++i) //mesh
    {
        x.at(i) = x.at(i - 1) + h;
    }

    std::vector<double> a(N), b(N), c(N), d(N), e(N), f(N);
    const double g0 = E * I;
    const double h2 = std::pow(h, 2), h4 = std::pow(h, 4);
    for (int i = 2; i < N - 2; ++i)
    {
        const double g1_x = g1(x.at(i));
        const double g2_x = g2(x.at(i));
        const double g3_x = g3(x.at(i));
        a.at(i) = g0 / h4;
        b.at(i) = g0 * 4 / h4 + g1_x / h2;
        c.at(i) = g0 * 6 / h4 + g1_x * 2 / h2 + g2_x;
        d.at(i) = g0 * 4 / h4 + g1_x / h2;
        e.at(i) = g0 / h4;
        f.at(i) = g3_x;
    }
    c.at(0) = -1;
    d.at(0) = 0;
    e.at(0) = 1;

    b.at(1) = 0;
    c.at(1) = 1;
    d.at(1) = 0;
    e.at(1) = 0;

    a.at(N - 2) = 1;
    b.at(N - 2) = 3;
    c.at(N - 2) = 3;
    d.at(N - 2) = 1;

    a.at(N - 1) = 1;
    b.at(N - 1) = 2;
    c.at(N - 1) = 1;

    MatrixSolver::solve5PointRun(a, b, c, d, e, f, y);

    for (auto val : y) std::cout << val << std::endl;
    std::cout << "Max w: " << *std::min_element(y.begin(), y.end()) << std::endl;
    //test();

    std::ofstream out("D:\\study\\NIR 3 Semester\\lab_2\\program\\program\\result.plt");
    out << "VARIABLES = \"x\", \"w\"" << std::endl;
    for (int i = 0; i < x.size(); ++i)
    {
        out << x.at(i) << " " << std::abs(y.at(i)) << std::endl;
    }
    return 0;
}
