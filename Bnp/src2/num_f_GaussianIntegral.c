/* ======================================================
   FILENAME:  num_f_GaussIntegral.c

   PURPOSE:   Gaussian quadrature : integration (a >0 to b>0) or (a to infinity) or (b to infinity)
                                with weight function exp(-x*x/2)/(2*Pi)^0.5   (NUMC p158)
   ====================================================== */

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_proba.h"
#include "num_h_tridiagQLi.h"
#include "utallhdr.h"
#include "utconst.h"

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*order eigen values and corresponding eigenvectors in descending order*/

static void eigsrt(double d[], double** v, int n)

/*Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
(§11.1) or tqli (§11.3), this routine sorts the eigenvalues into descending order, and rearranges
the columns of v correspondingly. The method is straight insertion.*/

{
    int    k, j, i;
    double p;

    for (i = 1; i < n; i++)
    {
        p = d[k = i];
        for (j = i + 1; j <= n; j++)
            if (d[j] >= p)
                p = d[k = j];
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; j++)
            {
                p       = v[j][i];
                v[j][i] = v[j][k];
                v[j][k] = p;
            }
        }
    }
}

static void tqli(double d[], double e[], int n, double** z)

/*QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 §11.2. On
input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns
the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix,
with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix
that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].*/

{
    double pythag(double a, double b);
    int    m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;
    for (i = 2; i <= n; i++)
        e[i - 1] = e[i];
    e[n] = 0.0;
    for (l = 1; l <= n; l++)
    {
        iter = 0;
        do
        {
            for (m = l; m <= n - 1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if ((fabs(e[m]) + dd) == dd)
                    break;
            }
            if (m != l)
            {
                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p     = 0.0;
                for (i = m - 1; i >= l; i--)
                {
                    f        = s * e[i];
                    b        = c * e[i];
                    e[i + 1] = (r = pythag(f, g));
                    if (r == 0.0)
                    {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s        = f / r;
                    c        = g / r;
                    g        = d[i + 1] - p;
                    r        = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g        = c * r - b;

                    for (k = 1; k <= n; k++)
                    {
                        f           = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i]     = c * z[k][i] - s * f;
                    }
                }
                if (r == 0.0 && i >= l)
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
}

Err GaussianIntegral(
    double  start,  // start < end
    double  end,
    int     Is_start_Infinity,  //=1 if integral from -Infinity to b
    int     Is_end_Infinity,    //=1 if integral from a to Infinity
    int     n,                  // n has to be greater than or equal to 1
    double* x,
    double* w)

{
    Err     err = NULL;
    int     i, j, looptmp;
    double *moments, *a, *b, **sig, **z, rec1, rec2;

    moments = dvector(1, 2 * n);

    // computation of the canonic moments

    if (Is_start_Infinity == 1)
    {
        if (Is_end_Infinity == 1)
        {
            moments[1] = 1;
            moments[2] = 0;
            rec1       = 0;
            rec2       = 0;
        }
        else
        {
            moments[1] = norm_accurate(end);
            moments[2] = INV_SQRT_TWO_PI * (-exp(-0.5 * end * end));
            rec1       = 0;
            rec2       = INV_SQRT_TWO_PI * exp(-0.5 * end * end);
        }
    }
    else
    {
        if (Is_end_Infinity == 1)
        {
            moments[1] = 1 - norm_accurate(start);
            moments[2] = INV_SQRT_TWO_PI * (exp(-start * start * 0.5));
            rec1       = INV_SQRT_TWO_PI * exp(-0.5 * start * start);
            rec2       = 0;
        }
        else
        {
            moments[1] = norm_accurate(end) - norm_accurate(start);
            moments[2] = INV_SQRT_TWO_PI * (exp(-start * start * 0.5) - exp(-0.5 * end * end));
            rec1       = INV_SQRT_TWO_PI * exp(-0.5 * start * start);
            rec2       = INV_SQRT_TWO_PI * exp(-0.5 * end * end);
        }
    }
    // recurrence relation for the canonic moments
    for (i = 3; i <= 2 * n; i++)
        moments[i] = (i - 2) * moments[i - 2] + rec1 * pow(start, i - 2) - rec2 * pow(end, i - 2);

    /*computation of the recurrence coefficients to construct the orthonomial polynomials
            (NUM_C p 159) with coefficients alpha and beta =0 */

    sig     = matrix(1, 2 * n + 1, 1, 2 * n + 1);
    looptmp = 2 * n;
    a       = dvector(1, n);
    b       = dvector(1, n);

    for (i = 3; i <= looptmp; i++)
        sig[1][i] = 0.0;
    looptmp++;

    for (i = 2; i <= looptmp; i++)
        sig[2][i] = moments[i - 1];
    a[1] = moments[2] / moments[1];
    b[1] = 0.0;

    for (i = 3; i <= n + 1; i++)
    {
        looptmp = 2 * n - i + 3;
        for (j = i; j <= looptmp; j++)
        {
            sig[i][j] = sig[i - 1][j + 1] + (0 - a[i - 2]) * sig[i - 1][j] -
                        b[i - 2] * sig[i - 2][j] + 0 * sig[i - 1][j - 1];
        }
        a[i - 1] = 0 + sig[i][i + 1] / sig[i][i] - sig[i - 1][i] / sig[i - 1][i - 1];
        b[i - 1] = sig[i][i] / sig[i - 1][i - 1];
    }

    /*computation of abcissas and weights for integration
            (NUM_C p 157) */

    z = matrix(1, n, 1, n);
    for (i = 1; i <= n; i++)
    {
        if (i != 1)
            b[i] = sqrt(b[i]);
        for (j = 1; j <= n; j++)
            z[i][j] = (double)(i == j);
    }

    tqli(a, b, n, z);
    eigsrt(a, z, n);
    for (i = 1; i <= n; i++)
    {
        x[i] = a[i];
        w[i] = moments[1] * z[1][i] * z[1][i];
    }

    free_matrix(z, 1, n, 1, n);
    free_dvector(moments, 1, 2 * n);
    free_dvector(a, 1, n);
    free_dvector(b, 1, n);
    free_matrix(sig, 1, 2 * n + 1, 1, 2 * n + 1);

    return err;
}

void lanczos(
    double start_point, double end_point, int n, int k, double* x, double* w, int weight_function)

{
    int      i, j;
    double **z, *b, *y, *q, *q_prev, *v, *alpha, *beta, *weight, norm_b = 0.0;

    y      = dvector(1, n);
    weight = dvector(1, n);
    b      = dvector(1, n);
    q      = dvector(1, n);
    q_prev = dvector(1, n);
    v      = dvector(1, n);
    alpha  = dvector(1, k);
    beta   = dvector(1, k);
    z      = matrix(1, k, 1, k);

    gauleg(start_point, end_point, y, weight, n);

    for (i = 1; i <= n; i++)
    {
        switch (weight_function)
        {
        case 1:
            b[i] = sqrt(
                weight[i] * (1 + (0.5 * 0.5 + 0.5 - y[i] * y[i]) /
                                     (exp(2.0 * log(fabs(0.5 * 0.5 + 0.5 - y[i] * y[i]))) +
                                      exp(2.0 * log(fabs((2 * 0.5 + 1.0) * y[i]))))));
            break;
        case 2:
            b[i] = sqrt(weight[i] * exp(-y[i] * y[i]) / (y[i] * y[i]));
            break;
        default:
            b[i] = sqrt(weight[i]);
            break;
        }

        norm_b += b[i] * b[i];
        q[i] = b[i];
    }

    norm_b = sqrt(norm_b);
    for (i = 1; i <= n; i++)
        q[i] = q[i] / norm_b;

    for (i = 1; i <= k; i++)
    {
        for (j = 1; j <= n; j++)
            v[j] = y[j] * q[j] - beta[i] * q_prev[j];
        for (j = 1; j <= n; j++)
            alpha[i] += v[j] * q[j];
        for (j = 1; j <= n; j++)
            v[j] = v[j] - alpha[i] * q[j];
        for (j = 1; j <= n; j++)
            beta[((i + 1) > k ? 1 : i + 1)] += v[j] * v[j];
        beta[i + 1] = sqrt(beta[i + 1]);
        for (j = 1; j <= n; j++)
        {
            q_prev[j] = q[j];
            q[j]      = v[j] / beta[((i + 1) > k ? 1 : i + 1)];
        }
    }

    for (i = 1; i <= k; i++)
    {
        for (j = 1; j <= k; j++)
            z[i][j] = (double)(i == j);
    }

    tqli(alpha, beta, k, z);
    eigsrt(alpha, z, k);

    for (i = 1; i <= k; i++)
    {
        x[i] = alpha[i];
        w[i] = norm_b * norm_b * z[1][i] * z[1][i];
    }

    free_dvector(y, 1, n);
    free_dvector(b, 1, n);
    free_dvector(q, 1, n);
    free_dvector(q_prev, 1, n);
    free_dvector(v, 1, n);
    free_dvector(weight, 1, n);
    free_dvector(alpha, 1, k);
    free_matrix(z, 1, k, 1, k);
}
