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

NumericVector densitiesP(NumericMatrix e_u)
{
    int n = e_u.nrow();
    int n_choices = e_u.ncol();
    NumericVector result(n);

    int n_phi = n_choices - 2;
    int last_index = n_choices - 1;
    LogicalVector phi(n_phi);

    do {
        int sgn = 1 - 2 * (sum(phi) % 2);
        for (int i = 0; i < n; i++)
        {
            double phi_dot_omega = 0;
            for (int j = 0; j < n_phi; j++)
            {
                if (phi[j])
                    phi_dot_omega += e_u(i, j + 1);
            }
            result[i] += sgn / (1 + phi_dot_omega / e_u(i, last_index));
        }
        phi = nextCombination(phi);
    } while (is_true(any(phi)));

    return result;
}

// [[Rcpp::export]]
NumericVector logDensitiesBestWorst(NumericMatrix e_u, NumericVector weights)
{
    int n = e_u.nrow();
    int n_choices = e_u.ncol();
    NumericVector result(n);
    if (n_choices > 2)
    {
        NumericVector densities_p = densitiesP(e_u);
        for (int i = 0; i < n; i++)
        {
            double total = 0;
            for (int j = 0; j < n_choices; j++)
                total += e_u(i, j);
            result[i] = (log(densities_p[i]) + log(e_u(i, 0)) - log(total)) * weights[i];
        }
    }
    else if (n_choices == 2)
    {
        for (int i = 0; i < n; i++)
            result[i] = (log(e_u(i, 0)) - log(e_u(i, 0) + e_u(i, 1))) * weights[i];
    }
    return result;
}

// [[Rcpp::export]]
double logDensityBestWorst(NumericMatrix e_u, NumericVector weights)
{
    return sum(logDensitiesBestWorst(e_u, weights));
}

// [[Rcpp::export]]
NumericVector gradientBestWorst(NumericMatrix e_u, IntegerMatrix x, NumericVector weights, int n_pars)
{
    int n = e_u.nrow();
    int n_choices = e_u.ncol();
    NumericVector result(n_pars);

    if (n_choices < 2)
        return result;

    // Gradient component for the best selection
    for (int i = 0; i < n; i++)
    {
        double non_best = 0;
        for (int j = 1; j < n_choices; j++)
            non_best += e_u(i, j);
        double total = e_u(i, 0) + non_best;
        if (x(i, 0) > 0)
            result[x(i, 0) - 1] += non_best * weights[i] / total;
        for (int j = 1; j < n_choices; j++)
            if (x(i, j) > 0)
                result[x(i, j) - 1] += -e_u(i, j) * weights[i] / total;
    }

    if (n_choices > 3)
    {
        int n_phi = n_choices - 3;
        int last_index = n_choices - 1;
        LogicalVector phi(n_phi);
        int n_selected = n_choices - 2;
        NumericMatrix dP_dOmega_selected(n, n_selected);

        // Calculate dP/dOmega for the items that are neither best or worst
        // by looping over all combinations
        do {
            int sgn = 1 - 2 * (sum(phi) % 2);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n_selected; j++)
                {
                    double phi_dot_omega = 0;
                    for (int k = 0; k < n_phi; k++)
                    {
                        if (phi[k])
                        {
                            if (k < j)
                                phi_dot_omega += e_u(i, k + 1);
                            else
                                phi_dot_omega += e_u(i, k + 2);
                        }
                    }
                    double tmp = e_u(i, j + 1) + e_u(i, last_index) + phi_dot_omega;
                    dP_dOmega_selected(i, j) += sgn * (e_u(i, last_index) / tmp) / tmp;
                }
            }
            phi = nextCombination(phi);
        } while (is_true(any(phi)));

        NumericVector densities_p = densitiesP(e_u);
        for (int i = 0; i < n; i++)
        {
            // Gradient for items that are neither best or worst
            for (int c = 1; c < last_index; c++)
            {
                if (x(i, c) > 0)
                    result[x(i, c) - 1] += (dP_dOmega_selected(i, c - 1) * e_u(i, c)) * weights[i] / densities_p[i];
            }

            // Gradient for worst item
            if (x(i, last_index) > 0)
            {
                double dP_dOmega_times_dn = 0;
                for (int l = 0; l < n_selected; l++)
                    dP_dOmega_times_dn += e_u(i, l + 1) * dP_dOmega_selected(i, l);
                result[x(i, last_index) - 1] += -dP_dOmega_times_dn * weights[i] / densities_p[i];
            }
        }
    }
    else if (n_choices == 3)
    {
        NumericVector densities_p = densitiesP(e_u);
        for (int i = 0; i < n; i++)
        {
            double tmp = e_u(i, 1) + e_u(i, 2);
            double gr = (e_u(i, 1) / tmp) * (e_u(i, 2) / tmp) * weights[i] / densities_p[i];
            if (x(i, 1) > 0)
                result[x(i, 1) - 1] += gr;
            if (x(i, 2) > 0)
                result[x(i, 2) - 1] += -gr;
        }
    }

    return result;
}


