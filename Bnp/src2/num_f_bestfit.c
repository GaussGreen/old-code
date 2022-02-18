/* ===========================================================
   FILENAME :    num_f_bestfit.c

   PURPOSE:      multi-dimensional best-fit algorithm
   =========================================================== */

#include "math.h"
#include "num_h_allhdr.h"
#include "utallhdr.h"

/*	quadratic best fit on model:
        y[i] = a + <grad,x[i]> + tx[i]*hess*x[i]
        where i = 1..n
        and x has dimension dim

        inputs:
        x[1..n][1..dim]
        y[1..n]
        n
        dim

        outputs (must be allocated):
        a
        grad[1..dim]
        hess[1..dim[1..dim]
        min[1..dim] (estimated minimum)
*/

#define ROTATE(a, i, j, k, l)        \
    g       = a[i][j];               \
    h       = a[k][l];               \
    a[i][j] = g - s * (h + g * tau); \
    a[k][l] = h + s * (g - h * tau);

static Err jacobi(double** a, long n, double d[], double** v, long* nrot)
{
    long   j, iq, ip, i;
    double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

    b = vector(1, n);
    z = vector(1, n);
    for (ip = 1; ip <= n; ip++)
    {
        for (iq = 1; iq <= n; iq++)
            v[ip][iq] = 0.0;
        v[ip][ip] = 1.0;
    }
    for (ip = 1; ip <= n; ip++)
    {
        b[ip] = d[ip] = a[ip][ip];
        z[ip]         = 0.0;
    }
    *nrot = 0;
    for (i = 1; i <= 50; i++)
    {
        sm = 0.0;
        for (ip = 1; ip <= n - 1; ip++)
        {
            for (iq = ip + 1; iq <= n; iq++)
                sm += fabs(a[ip][iq]);
        }
        if (sm == 0.0)
        {
            free_vector(z, 1, n);
            free_vector(b, 1, n);
            return NULL;
        }
        if (i < 4)
            tresh = 0.2 * sm / (n * n);
        else
            tresh = 0.0;
        for (ip = 1; ip <= n - 1; ip++)
        {
            for (iq = ip + 1; iq <= n; iq++)
            {
                g = 100.0 * fabs(a[ip][iq]);
                if (i > 4 && (double)(fabs(d[ip]) + g) == (double)fabs(d[ip]) &&
                    (double)(fabs(d[iq]) + g) == (double)fabs(d[iq]))
                    a[ip][iq] = 0.0;
                else if (fabs(a[ip][iq]) > tresh)
                {
                    h = d[iq] - d[ip];
                    if ((double)(fabs(h) + g) == (double)fabs(h))
                        t = (a[ip][iq]) / h;
                    else
                    {
                        theta = 0.5 * h / (a[ip][iq]);
                        t     = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        if (theta < 0.0)
                            t = -t;
                    }
                    c   = 1.0 / sqrt(1 + t * t);
                    s   = t * c;
                    tau = s / (1.0 + c);
                    h   = t * a[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    d[ip] -= h;
                    d[iq] += h;
                    a[ip][iq] = 0.0;
                    for (j = 1; j <= ip - 1; j++)
                    {
                        ROTATE(a, j, ip, j, iq)
                    }
                    for (j = ip + 1; j <= iq - 1; j++)
                    {
                        ROTATE(a, ip, j, j, iq)
                    }
                    for (j = iq + 1; j <= n; j++)
                    {
                        ROTATE(a, ip, j, iq, j)
                    }
                    for (j = 1; j <= n; j++)
                    {
                        ROTATE(v, j, ip, j, iq)
                    }
                    ++(*nrot);
                }
            }
        }
        for (ip = 1; ip <= n; ip++)
        {
            b[ip] += z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
        }
    }
    return serror("Too many iterations in routine jacobi");
}
#undef ROTATE

Err quadr_best_fit(
    double** x, double* y, long n, long dim, double* a, double* grad, double** hess, double* min)
{
    /* declaration and allocation */

    long i, j, k, l, nrot, def, nvar = (dim + 1) * (dim + 2) / 2, *indx;

    double **X = dmatrix(1, n, 1, nvar), **Y = dmatrix(1, n, 1, 1), **transX = NULL,
           **transXX = NULL, **inv_transXX = NULL, **transXY = NULL, **inv_transXX_transXY = NULL,
           **inv_hess = NULL, **grad_matrix = dmatrix(1, dim, 1, 1), **inv_hess_grad = NULL,
           *eigen_val = dvector(1, dim), **eigen_vec = dmatrix(1, dim, 1, dim);

    Err err;

    /*	transform the problem into unidimensional one
            variables are 1, x1, ..., xdim, x1^2, x2*x1, x2^2, x3*x1, x3*x2, x3^2, ..., xdim^2 */

    /*	values are transfered into X[1..n][1..nvar] */

    for (i = 1; i <= n; i++)
    {
        Y[i][1] = y[i];

        X[i][1] = 1.0;

        for (j = 1; j <= dim; j++)
        {
            X[i][j + 1] = x[i][j];
        }

        l = 0;
        for (j = 1; j <= dim; j++)
        {
            for (k = 1; k <= j; k++)
            {
                l++;
                X[i][dim + 1 + l] = x[i][j] * x[i][k];
            }
        }
    }

    /*	find best fit using formula P = inv (trans (X) * X) * (trans (X) * Y) */

    transX = transpose_matrix(X, 1, n, 1, nvar);
    if (!transX)
        return serror("error in transpose_matrix");

    transXX = product_matrix(transX, 1, nvar, 1, n, X, 1, n, 1, nvar);
    if (!transXX)
        return serror("error in product_matrix");

    transXY = product_matrix(transX, 1, nvar, 1, n, Y, 1, n, 1, 1);
    if (!transXY)
        return serror("error in product_matrix");

    inv_transXX = inverse_matrix(transXX, 1, nvar);
    if (!inv_transXX)
        return serror("error in inverse_matrix");

    inv_transXX_transXY = product_matrix(inv_transXX, 1, nvar, 1, nvar, transXY, 1, nvar, 1, 1);
    if (!inv_transXX_transXY)
        return serror("error in inverse_matrix");

    /* interpret results in terms of a, grad, hess with y[i] = a + <grad,x[i]> + tx[i]*hess*x[i] */

    *a = inv_transXX_transXY[1][1];

    for (j = 1; j <= dim; j++)
    {
        grad[j] = inv_transXX_transXY[j + 1][1];
    }

    l = 0;
    for (j = 1; j <= dim; j++)
    {
        for (k = 1; k <= j; k++)
        {
            l++;
            hess[j][k] = hess[k][j] = inv_transXX_transXY[dim + 1 + l][1];
        }
    }

    /* check if hess is definite positive and put the sign in def */

    err = jacobi(hess, dim, eigen_val, eigen_vec, &nrot);
    if (err)
    {
        return err;
    }

    def = +1;
    for (j = 1; j <= dim; j++)
    {
        if (eigen_val[j] <= 0.0)
            def = -1;
    }

    /* if hess is definite positive, find minimum with formula min = - 0.5 * inv (hess) * grad */

    if (def > 0)
    {
        smessage(" hessian is definite positive");
        smessage(" returning quadratic minimum...");
        inv_hess = inverse_matrix(hess, 1, dim);
        if (!inv_hess)
            return serror("error in inverse_matrix");

        for (i = 1; i <= dim; i++)
            grad_matrix[i][1] = grad[i];

        inv_hess_grad = product_matrix(inv_hess, 1, dim, 1, dim, grad_matrix, 1, dim, 1, 1);
        if (!inv_hess_grad)
            return serror("error in product_matrix");

        for (i = 1; i <= dim; i++)
            min[i] = -0.5 * inv_hess_grad[i][1];

        free_dmatrix(X, 1, n, 1, nvar);
        free_dmatrix(Y, 1, n, 1, 1);
        free_dmatrix(transX, 1, nvar, 1, n);
        free_dmatrix(transXX, 1, nvar, 1, nvar);
        free_dmatrix(transXY, 1, nvar, 1, 1);
        free_dmatrix(inv_transXX, 1, nvar, 1, nvar);
        free_dmatrix(inv_transXX_transXY, 1, nvar, 1, 1);
        free_dmatrix(grad_matrix, 1, dim, 1, 1);
        free_dmatrix(inv_hess, 1, dim, 1, dim);
        free_dmatrix(inv_hess_grad, 1, dim, 1, 1);
        free_dvector(eigen_val, 1, dim);
        free_dmatrix(eigen_vec, 1, dim, 1, dim);
    }

    /* otherwise give best point */
    else
    {
        smessage(" hessian is NOT definite positive");
        smessage(" returning best point...");
        indx = lngvector(1, n);
        indexx(n, y, indx);
        memcpy(&(min[1]), &(x[indx[1]][1]), dim * sizeof(double));

        free_dmatrix(X, 1, n, 1, nvar);
        free_dmatrix(Y, 1, n, 1, 1);
        free_dmatrix(transX, 1, nvar, 1, n);
        free_dmatrix(transXX, 1, nvar, 1, nvar);
        free_dmatrix(transXY, 1, nvar, 1, 1);
        free_dmatrix(inv_transXX, 1, nvar, 1, nvar);
        free_dmatrix(inv_transXX_transXY, 1, nvar, 1, 1);
        free_lngvector(indx, 1, n);
        free_dvector(eigen_val, 1, dim);
        free_dmatrix(eigen_vec, 1, dim, 1, dim);
    }

    return NULL;
}
