#include <Rcpp.h>
using namespace Rcpp;

LogicalVector nextCombination(LogicalVector comb)
{
    for (int i = 0; i < comb.size(); i++)
    {
        comb[i] = !comb[i];
        if (comb[i])
            break;
    }
    return comb;
}

// [[Rcpp::export]]
double logDensityBestWorst(NumericMatrix e_u, NumericVector weights)
{
    int n_choices = e_u.ncol();
    int n_phi = n_choices - 2;
    int last_index = n_choices - 1;
    int n = e_u.nrow();
    LogicalVector phi(n_choices - 2);
    NumericVector p(n);

    do {
        int sgn = 1 - 2 * (sum(phi) % 2);
        for (int i = 0; i < n; i++)
        {
            double phi_dot_omega = 0;
            for (int j = 0; j < n_phi; j++)
                if (phi[j])
                    phi_dot_omega += e_u(i, j + 1);
                p[i] += sgn / (1 + phi_dot_omega / e_u(i, last_index));
        }
        phi = nextCombination(phi);
    } while (is_true(any(phi)));

    double result = 0;
    for (int i = 0; i < n; i++)
    {
        double p_best = 0;
        for (int j = 0; j < n_choices; j++)
            p_best += e_u(i, j);
        result += (log(e_u(i, 0) / p_best) + log(p[i])) * weights[i];
    }
    return result;
}

// // [[Rcpp::export]]
// NumericVector gradientBestWorst(NumericMatrix e_u, NumericMatrix x)
// {
//
// }


