/*--------------------------------------------------------------
        FILE: LGM1FCalibToProduct.c
        PURPOSE: LGM1f calibration to any product
        AUTHOR: Dimitri Mayevski
        DATE: 11/07/2002
  --------------------------------------------------------------*/

#include "LGM1FCalibToProduct.h"

#include "math.h"
#include "num_h_zbrent.h"
#include "srt_h_all.h"

typedef struct
{
    long     d;
    double   t, exp_2lamdt, phi;
    double **dfs, *gammas;
    int *    nmat, nmkts;
    double   df;
    double   target;
} LGMExData;

typedef struct
{
    int         mkt_idx;
    SrtProduct* product;
    double      phi_prev;
    LGMExData*  ex_data;
    double      lam;
    int         nx;
    double      qto_fudge;
} TrySigmaData;

#define MAXNX 1000
#define MAXCPN 600

static Err try_sigma(double sig, double* diff, void* static_data)
{
    Err           err = NULL;
    TrySigmaData* p   = (TrySigmaData*)static_data;
    double        phi, std, sum, gam, payoff, qshft, c;
    double        x[MAXNX], w[MAXNX];
    double        dfs_d[MAXCPN], dfs_f[MAXCPN];
    double*       dfs[2] = {dfs_d, dfs_f};
    int           i, j, k;

    if (p->nx > MAXNX)
        return serror("nx should be no more than %d", MAXNX);
    if (p->ex_data->nmkts > 2)
        return serror("Maximum two markets supported");

    /* Calculate phi */
    phi = p->phi_prev * p->ex_data->exp_2lamdt +
          sig * sig / (2.0 * p->lam) * (1.0 - p->ex_data->exp_2lamdt);
    std = sqrt(phi);

    /* Quanto shift */
    qshft = -p->qto_fudge * sqrt(phi * p->ex_data->t);

    /* Calculate Gauss-Legendre abscissas */
    gauleg(qshft - 4.0 * std, qshft + 4.0 * std, x - 1, w - 1, p->nx);

    /* Copy forward DFs for non-calibrated markets */
    for (j = 0; j < p->ex_data->nmkts; j++)
        if (j != p->mkt_idx)
            memcpy(dfs[j], p->ex_data->dfs[j], p->ex_data->nmat[j] * sizeof(double));

    /* Integrate */
    for (i = 0, sum = 0.0; i < p->nx; i++)
    {
        /* Calculate dfs */
        for (k = 0; k < p->ex_data->nmat[p->mkt_idx]; k++)
        {
            gam = p->ex_data->gammas[k];
            dfs[p->mkt_idx][k] =
                p->ex_data->dfs[p->mkt_idx][k] * exp(-gam * (x[i] + 0.5 * gam * phi));
        }
        /* Calculate payoff at x[i] */
        err = p->product->Payoff(p->product, p->ex_data->d, p->ex_data->d, dfs, &payoff);
        if (err)
            return err;

        c = x[i] - qshft;
        sum += w[i] * INV_SQRT_TWO_PI * exp(-0.5 * c * c / phi) / std * payoff;
    }

    *diff = sum * p->ex_data->df - p->ex_data->target;
    return NULL;
}

static Err calib_sigma(
    int         mkt_idx,
    SrtProduct* product,
    int         isig,
    LGMExData*  ex_data,
    double      lam,
    double*     sig,
    int         nx,
    double      qto_fudge)
{
    Err          err = NULL;
    TrySigmaData data;
    double       x1, x2, f1, f2;
    double       tol = 1.0e-5;

    /* Fill in static data */
    data.mkt_idx   = mkt_idx;
    data.product   = product;
    data.phi_prev  = (isig == 0 ? 0.0 : ex_data[isig - 1].phi);
    data.ex_data   = &(ex_data[isig]);
    data.lam       = lam;
    data.nx        = nx;
    data.qto_fudge = qto_fudge;

    /* Take first guess from sig */
    x1  = 0.001; /*0.9 * sig[isig];*/
    x2  = 0.03;  /*1.1 * sig[isig];*/
    err = try_sigma(x1, &f1, &data);
    if (err)
        return err;
    err = try_sigma(x2, &f2, &data);
    if (err)
        return err;

    /* Bracket the root */
    /*	err = num_f_zbrac(try_sigma, &x1, &x2, &f1, &f2, &data);
            if (err) return err;*/

    /* Calibrate */
    err = num_f_zbrent(try_sigma, x1, x2, f1, f2, tol, &data, &(sig[isig]));
    if (err)
        return err;

    /* Calculate phi at ex_data[isig] */
    ex_data[isig].phi = data.phi_prev * ex_data[isig].exp_2lamdt +
                        sig[isig] * sig[isig] / (2.0 * lam) * (1.0 - ex_data[isig].exp_2lamdt);

    return NULL;
}

Err srt_f_lgm1f_calib_to_product(
    char**      yc_names,
    int         n_mkts,
    int         mkt_idx,
    SrtProduct* product,
    int         nex,
    long*       ex_dates,
    double      lam,
    double*     sig,
    int         nx,
    double      qto_fudge)
{
    Err         err = NULL;
    int         i, j, k;
    LGMExData*  ex_data = NULL;
    SrtCurvePtr yc_ptr;
    long        today;
    long**      dates = NULL;
    int         n_unds, *n_dfs = NULL;
    double      t1, t2;

    yc_ptr = lookup_curve(yc_names[0]);
    if (!yc_ptr)
        return serror("Yield curve %s not found", yc_names[0]);
    today = get_today_from_curve(yc_ptr);

    ex_data = (LGMExData*)calloc(nex, sizeof(LGMExData));
    if (!ex_data)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }
    memset(ex_data, 0, nex * sizeof(LGMExData));

    /* Fill in ex_data */
    for (i = 0, t1 = 0.0; i < nex; i++)
    {
        ex_data[i].d = ex_dates[i];
        ex_data[i].t = t2     = (ex_dates[i] - today) * YEARS_IN_DAY;
        ex_data[i].exp_2lamdt = exp(-2.0 * lam * (t2 - t1));
        ex_data[i].df         = swp_f_df(today, ex_dates[i], yc_names[0]); /* domestic */
        if (!product->RequestDfDates)                                      /* if no DFs required */
        {
            err = product->Payoff(product, today, ex_dates[i], NULL, &(ex_data[i].target));
            if (err)
                goto FREE_RETURN;
        }
        else /* otherwise ask the product */
        {
            err = product->RequestDfDates(product, ex_dates[i], &n_unds, &n_dfs, &dates);
            if (err)
                goto FREE_RETURN;
            if (n_unds > n_mkts)
            {
                err = serror("Foreign market not specified");
                goto FREE_RETURN;
            }

            ex_data[i].nmkts = n_unds;
            ex_data[i].dfs   = (double**)calloc(n_unds, sizeof(double*));
            memset(ex_data[i].dfs, 0, n_unds * sizeof(double*));
            ex_data[i].nmat = (int*)calloc(n_unds, sizeof(int));
            memcpy(ex_data[i].nmat, n_dfs, n_unds * sizeof(int));

            /* Get DFs needed to calculate target */
            for (j = 0; j < n_unds; j++)
            {
                ex_data[i].dfs[j] = (double*)calloc(n_dfs[j], sizeof(double));
                for (k = 0; k < n_dfs[j]; k++)
                    ex_data[i].dfs[j][k] = swp_f_df(today, dates[j][k], yc_names[j]);

                if (j == mkt_idx)
                {
                    ex_data[i].gammas = (double*)calloc(n_dfs[j], sizeof(double));
                    for (k = 0; k < n_dfs[j]; k++)
                    {
                        ex_data[i].gammas[k] =
                            (1.0 - exp(-lam * (dates[j][k] - ex_dates[i]) * YEARS_IN_DAY)) / lam;
                    }
                }
            }

            /* Calculate target */
            err =
                product->Payoff(product, today, ex_dates[i], ex_data[i].dfs, &(ex_data[i].target));
            if (err)
                goto FREE_RETURN;

            /* Now replace dfs by forward dfs */
            for (j = 0; j < n_unds; j++)
                for (k = 0; k < n_dfs[j]; k++)
                    ex_data[i].dfs[j][k] = swp_f_df(ex_dates[i], dates[j][k], yc_names[j]);

            /* Free dates and n_dfs */
            for (j = 0; j < n_unds; j++)
                free(dates[j]);
            free(dates);
            dates = NULL;
            free(n_dfs);
            n_dfs = NULL;

        } /* if (!product->RequestDfDates) */
        t1 = t2;
    } /* for (i=0, t1=0.0; i < nex; i++) */

    /* At this point ex_data is fully initialized except the phi field */
    /* Now calibrate all sigmas */

    for (i = 0; i < nex; i++)
    {
        err = calib_sigma(mkt_idx, product, i, ex_data, lam, sig, nx, qto_fudge);
        if (err)
            goto FREE_RETURN;
    }

FREE_RETURN:

    if (dates)
        for (j = 0; j < n_unds; j++)
            free(dates[j]);
    free(dates);
    free(n_dfs);

    if (ex_data)
        for (i = 0; i < nex; i++)
        {
            if (ex_data[i].dfs)
                for (j = 0; j < ex_data[i].nmkts; j++)
                    free(ex_data[i].dfs[j]);
            free(ex_data[i].dfs);
            free(ex_data[i].nmat);
            free(ex_data[i].gammas);
        }
    free(ex_data);

    return err;
}
