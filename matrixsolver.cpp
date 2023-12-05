#include "matrixsolver.h"

MatrixSolver::MatrixSolver()
{

}

// Алгоритм немонотонной пятиточечной прогонки
// [Самарский А. А., Николаев Е. С. Методы решения сеточных уравнений. – 1978., стр 101]
void MatrixSolver::solve5PointRun(const std::vector<double> &a, const std::vector<double> &b,
                                  const std::vector<double> &c, const std::vector<double> &d,
                                  const std::vector<double> &e, const std::vector<double> &f,
                                  std::vector<double> &y)
{
    const int N = a.size() - 1;
    const int size = a.size();
    std::vector<double> alpha(size), beta(size),
            gamma(size + 1), delta(size);
    std::vector<int> kappa(size), nu(size), teta(size);

    double C = c.at(0), D = d.at(0), B = b.at(1),
            Q = c.at(1), S = a.at(2), T = b.at(2),
            R = 0, A = a.at(3), F = f.at(0),
            PHI = f.at(1), G = f.at(2), H = f.at(3);
    kappa.at(0) = 0;
    nu.at(0) = 1;
    for (int i = 0; i <= N - 2; ++i)
    {
        if (std::abs(C) >= std::abs(D) && std::abs(C) >= std::abs(e.at(i)))
        {

            alpha.at(i + 1) = D / C;
            beta.at(i + 1) = e.at(i) / C;
            gamma.at(i + 1) = F / C;

            C = Q - B * alpha.at(i + 1);
            D = d.at(i + 1) - B * beta.at(i + 1);
            F = PHI + B * gamma.at(i + 1);

            B = T - S * alpha.at(i + 1);
            Q = c.at(i + 2) - S * beta.at(i + 1);
            PHI = G - S * gamma.at(i + 1);

            if (i != N - 2)
            {
                S = A - R * alpha.at(i + 1);
                T = b.at(i + 3) - R * beta.at(i + 1);
                G = H + R * gamma.at(i + 1);
            }

            if (i < N - 3)
            {
                R = 0;
                A = a.at(i + 4);
                H = f.at(i + 4);
            }

            teta.at(i + 1) = kappa.at(i);
            kappa.at(i + 1) = nu.at(i);
            nu.at(i + 1) = i + 2;

        }
        else if (std::abs(D) > std::abs(C) && std::abs(D) >= std::abs(e.at(i)))
        {

            alpha.at(i + 1) = C / D;
            beta.at(i + 1) = -e.at(i) / D;
            gamma.at(i + 1) = -F / D;

            C = Q * alpha.at(i + 1) - B;
            D = Q * beta.at(i + 1) + d.at(i + 1);
            F = PHI - Q * gamma.at(i + 1);

            B = T * alpha.at(i + 1) - S;
            Q = T * beta.at(i + 1) + c.at(i + 2);
            PHI = T * gamma.at(i + 1) + G;


            if (i != N - 2)
            {
                S = A * alpha.at(i + 1) - R;
                T = A * beta.at(i + 1) + b.at(i + 3);
                G = H - A * gamma.at(i + 1);
            }

            if (i < N - 3)
            {
                R = 0;
                A = a.at(i + 4);
                H = f.at(i + 4);
            }

            teta.at(i + 1) = nu.at(i);
            kappa.at(i + 1) = kappa.at(i);
            nu.at(i + 1) = i + 2;

        }
        else if (std::abs(e.at(i)) > C && std::abs(e.at(i)) > std::abs(D))
        {

            alpha.at(i + 1) = D / e.at(i);
            beta.at(i + 1) = C / e.at(i);
            gamma.at(i + 1) = F / e.at(i);

            C = Q - d.at(i + 1) * alpha.at(i + 1);
            D = B - d.at(i + 1) * beta.at(i + 1);
            F = PHI + d.at(i + 1) * gamma.at(i + 1);

            B = T - c.at(i + 2) * alpha.at(i + 1);
            Q = S - c.at(i + 2) * beta.at(i + 1);
            PHI = G - c.at(i + 2) * gamma.at(i + 1);


            if (i != N - 2)
            {
                S = A - b.at(i + 3) * alpha.at(i + 1);
                T = R - b.at(i + 3) * beta.at(i + 1);
                G = H + b.at(i + 3) * gamma.at(i + 1);
            }

            if (i < N - 3)
            {
                R = -a.at(i + 4) * alpha.at(i + 1);
                A = -a.at(i + 4) * beta.at(i + 1);
                H = f.at(i + 4) - a.at(i + 4) * gamma.at(i + 1);
            }

            teta.at(i + 1) = i + 2;
            kappa.at(i + 1) = nu.at(i);
            nu.at(i + 1) = kappa.at(i);

        }
    }

    if (std::abs(C) >= std::abs(D))
    {
        alpha.at(N) = D / C;
        gamma.at(N) = F / C;
        gamma.at(N + 1) = (PHI + B * gamma.at(N)) / (Q - B * alpha.at(N));
        teta.at(N) = kappa.at(N - 1);
        kappa.at(N) = nu.at(N - 1);
    }
    else if (std::abs(D) > std::abs(C))
    {
        alpha.at(N) = C / D;
        gamma.at(N) = -F / D;
        gamma.at(N + 1) = (PHI - Q * gamma.at(N)) / (Q * alpha.at(N) - B);
        teta.at(N) = nu.at(N - 1);
        kappa.at(N) = kappa.at(N - 1);
    }

    int m = teta.at(N), n = kappa.at(N), k = 0;
    y.at(n) = gamma.at(N + 1);
    y.at(m) = alpha.at(N) * y.at(n) + gamma.at(N);
    for (int i = N - 2; i >= 0; --i)
    {
        m = teta.at(i + 1);
        n = kappa.at(i + 1);
        k = nu.at(i + 1);
        y.at(m) = alpha.at(i + 1) * y.at(n) - beta.at(i + 1) * y.at(k) + gamma.at(i + 1);
    }
}
