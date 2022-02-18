/* ==================================================================
   FILE_NAME:	num_f_sobenberg.c

   PURPOSE:     minimisation of a function using SobenBerg method

   DESCRIPTION: fits the equation funcs (param, data[i]) = target[i]
                    with weighted least squares error function
                                using sobenberg algorithm
        (cf. optimal calibration document)

   MAIN FUNCTION:  Err sobenberg:
        double data[1..ndata]
        double target[1..ndata]
        double weight[1..ndata]
        double param[1..nparam]
        double min_param[1..nparam] (lower bounds for param)
        double max_param[1..nparam] (upper bounds for param)
        long niter (LM iterations)
        long tot_pts (initial sobol sample)
        long bst_pts (final sobol sample)
        long nclust (number of clusters)
        Err funcs (double data, double param[], double *return_value, int nparam)
        Err dfuncs (double data, double param[], double *return_value, int nparam]
   ================================================================== */

#include "math.h"
#include "num_h_allhdr.h"
#include "utallhdr.h"

#define CLUST_ALGO_NITER 500
#define BIG 1.0e+16
#define BP *1.0e-4

#define alloc_everything                                                    \
    double**  pts             = dmatrix(1, tot_pts, 1, nparam);             \
    double*   values          = dvector(1, tot_pts);                        \
    long*     index           = lngvector(1, tot_pts);                      \
    double**  bpts            = dmatrix(1, bst_pts, 1, nparam);             \
    double*   bvalues         = dvector(1, bst_pts);                        \
    double**  kernel          = dmatrix(1, nclust, 1, nparam);              \
    long*     npts_in_clust   = lngvector(1, nclust);                       \
    long*     kernel_of_point = lngvector(1, bst_pts);                      \
    double*** x               = f3tensor(1, nclust, 1, bst_pts, 1, nparam); \
    double**  y               = dmatrix(1, nclust, 1, bst_pts);             \
    double*   grad            = dvector(1, nparam);                         \
    double**  hess            = dmatrix(1, nparam, 1, nparam);              \
    double**  min             = dmatrix(1, nclust, 1, nparam);              \
    double*   minval          = dvector(1, nclust);                         \
    long*     indexc          = lngvector(1, nclust);

#define free_everything                                     \
    {                                                       \
        free_dmatrix(pts, 1, tot_pts, 1, nparam);           \
        free_dvector(values, 1, tot_pts);                   \
        free_lngvector(index, 1, tot_pts);                  \
        free_dmatrix(bpts, 1, bst_pts, 1, nparam);          \
        free_dvector(bvalues, 1, bst_pts);                  \
        free_dmatrix(kernel, 1, nclust, 1, nparam);         \
        free_lngvector(npts_in_clust, 1, nclust);           \
        free_lngvector(kernel_of_point, 1, bst_pts);        \
        free_f3tensor(x, 1, nclust, 1, bst_pts, 1, nparam); \
        free_dmatrix(y, 1, nclust, 1, bst_pts);             \
        free_dvector(grad, 1, nparam);                      \
        free_dmatrix(hess, 1, nparam, 1, nparam);           \
        free_dmatrix(min, 1, nclust, 1, nparam);            \
        free_dvector(minval, 1, nclust);                    \
        free_lngvector(indexc, 1, nclust);                  \
    }

static Err error_function(
    double* data,
    double* target,
    double* weight,
    long    ndata,
    double* param_0_1,
    double* min_param,
    double* max_param,
    long    nparam,
    Err (*funcs)(double, double[], double*, int),
    double* chisq)
{
    long   i;
    double funcs_value, *param = dvector(1, nparam);
    Err    err;

    for (i = 1; i <= nparam; i++)
    {
        param[i] = min_param[i] + (max_param[i] - min_param[i]) * param_0_1[i];
    }

    *chisq = 0;
    for (i = 1; i <= ndata; i++)
    {
        err = funcs(data[i], param, &funcs_value, nparam);
        if (err)
        {
            free_dvector(param, 1, nparam);
            return err;
        }
        *chisq += ((funcs_value - target[i]) / weight[i]) * ((funcs_value - target[i]) / weight[i]);
    }

    free_dvector(param, 1, nparam);
    return NULL;
}

/* ---------------------------------------------------------------- */
Err sobenberg(
    double* data,
    double* target,
    double* weight,
    long    ndata,
    double* param,
    double* min_param,
    double* max_param,
    long    nparam,
    long    tot_pts,
    long    bst_pts,
    long    nclust,
    char*   rsc_mth,
    long    niter,
    Err (*funcs)(double, double[], double*, int),
    Err (*dfuncs)(double, double[], double*, double[], int),
    double* chisq)
{
    /* Declarations */

    long   i, j, k, status = 0, *indexn;
    double a, dist, temp;
    Err    err;

    /* Allocations AND declarations */
    alloc_everything;

    /* Warning : about to start */

    smessage("Starting Sobenberg algorithm");
    smessage("");

    /* Sobol sampling of the initial space of parameters */

    smessage("Stage 1: Sobol sampling");
    smessage("");

    err = sobol_init(1, tot_pts, 1, 1, 1, nparam);
    if (err)
    {
        free_everything;
        return err;
    }

    for (i = 1; i <= tot_pts; i++)
    {
        err = sobol_vector(&pts[i][1], 0, nparam - 1);
        if (err)
        {
            free_everything;
            return err;
        }

        err = error_function(
            data, target, weight, ndata, pts[i], min_param, max_param, nparam, funcs, &(values[i]));
        if (err)
        {
            free_everything;
            return err;
        }
        smessage(" sobol point: %d", (int)i);
        smessage(" criteria: %.2f", 100 * sqrt(values[i]));
    }
    smessage("");
    err = sobol_free();
    if (err)
    {
        free_everything;
        return err;
    }

    /* Sort the points */
    smessage("Stage 2: Ranking");
    smessage("");

    err = indexx(tot_pts, values, index);
    if (err)
    {
        free_everything;
        return err;
    }

    for (i = 1; i <= bst_pts; i++)
    {
        memcpy(&(bpts[i][1]), &(pts[index[i]][1]), nparam * sizeof(double));
        err = error_function(
            data,
            target,
            weight,
            ndata,
            bpts[i],
            min_param,
            max_param,
            nparam,
            funcs,
            &(bvalues[i]));
        if (err)
        {
            free_everything;
            return err;
        }
    }

    /* clustering */

    smessage("Stage 3: Clustering");
    smessage("");

    find_initial_bars(bpts, kernel, bst_pts, nclust, nparam);

    for (i = 1; ((i <= CLUST_ALGO_NITER) || (status == 0)); i++)
    {
        status = clustering(bpts, kernel, bst_pts, nclust, nparam);
    }

    memset(&(npts_in_clust[1]), 0, nclust * sizeof(long));

    for (i = 1; i <= bst_pts; i++)
    {
        dist = BIG;
        for (j = 1; j <= nclust; j++)
        {
            temp = dist_sqr(bpts[i], kernel[j], nparam);
            if (temp < dist)
            {
                dist               = temp;
                kernel_of_point[i] = j;
            }
        }
        npts_in_clust[kernel_of_point[i]]++;
        memcpy(
            &(x[kernel_of_point[i]][npts_in_clust[kernel_of_point[i]]][1]),
            &(bpts[i][1]),
            nparam * sizeof(double));
        y[kernel_of_point[i]][npts_in_clust[kernel_of_point[i]]] = bvalues[i];
    }

    /* research of the starting point in each cluster */

    smessage("Stage 4: Starting point research in each cluster");
    smessage("");

    for (i = 1; i <= nclust; i++)
    {
        smessage(" cluster: %d", (int)i);
        smessage(" reseach method: %s", rsc_mth);

        /* quadratic best fit */
        if (!strcmp(rsc_mth, "BESTFIT"))
        {
            err = quadr_best_fit(x[i], y[i], npts_in_clust[i], nparam, &a, grad, hess, min[i]);
            if (err)
            {
                free_everything;
                return err;
            }
        }

        /* baricenter */
        else if (!strcmp(rsc_mth, "BARC"))
        {
            for (j = 1; j <= nparam; j++)
            {
                min[i][j] = 0.0;
                for (k = 1; k <= npts_in_clust[i]; k++)
                {
                    min[i][j] += x[i][k][j];
                }
                min[i][j] /= npts_in_clust[i];
            }
        }

        /* best point of the cluster */
        else if (!strcmp(rsc_mth, "BESTPOINT"))
        {
            indexn = lngvector(1, npts_in_clust[i]);
            err    = indexx(npts_in_clust[i], y[i], indexn);
            if (err)
            {
                free_everything;
                return err;
            }

            memcpy(&(min[i][1]), &(x[i][indexn[1]][1]), nparam * sizeof(double));
            free_lngvector(indexn, 1, npts_in_clust[i]);
        }

        else
        {
            free_everything;
            return serror("Unknown research method: %s", rsc_mth);
        }

        err = error_function(
            data, target, weight, ndata, min[i], min_param, max_param, nparam, funcs, &(minval[i]));
        if (err)
        {
            free_everything;
            return err;
        }
        smessage(" criteria cluster %d: %.2f", (int)i, 100 * sqrt(minval[i]));
        smessage("");
    }
    smessage("");

    /* process levenberg from each min */

    smessage("Stage 5: Levenberg processes");
    smessage("");

    for (i = 1; i <= nclust; i++)
    {
        smessage(" process: %d", (int)i);
        smessage("");

        for (j = 1; j <= nparam; j++)
        {
            min[i][j] = min_param[j] + (max_param[j] - min_param[j]) * min[i][j];
        }

        err = levenberg_marquardt(
            data, target, weight, ndata, min[i], nparam, niter, dfuncs, &(minval[i]));
        if (err)
        {
            free_everything;
            return err;
        }
    }

    err = indexx(nclust, minval, indexc);
    if (err)
    {
        free_everything;
        return err;
    }

    memcpy(&(param[1]), &(min[indexc[1]][1]), nparam * sizeof(double));
    *chisq = minval[indexc[1]];
    smessage(" terminal criteria: %.2f", 100 * sqrt(*chisq));
    smessage("");

    /* free memory */

    free_everything;

    return NULL;
} /* END Err sobenberg(...) */
