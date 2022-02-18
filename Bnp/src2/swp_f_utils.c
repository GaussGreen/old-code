/*--------------------------------------------------------------
        FILE: swp_f_utils.c
        PURPOSE: Useful swap utilities
        AUTHOR: Dimitri Mayevski
        DATE: 3/10/2002
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "swp_h_utils.h"

#define MAX_DATES 400

Err SwapElements_Init(
    long     start,
    long     end,
    char*    freqStr,
    char*    basisStr,
    char*    refrateStr,
    char*    recpayStr,
    double   fixed,
    int      elements,
    long     today,
    double** times,
    double** cflows,
    int*     len)
{
    Err             err = NULL;
    long            dates[MAX_DATES];
    double          cpns[MAX_DATES];
    long            idx[MAX_DATES];
    SrtCompounding  freq;
    SrtBasisCode    basis;
    SrtReceiverType rec_pay;
    int             i, j, k = 0;
    long            start_date, end_date, theo_date, act_date, temp_date, di, dj;
    double          sign, tmp;

    if (!elements)
        return serror("No swap elements specified");

    err = interp_rec_pay(recpayStr, &rec_pay);
    if (err)
        return err;

    if (rec_pay == SRT_RECEIVER)
        sign = 1.0;
    else
        sign = -1.0;

    start_date = bus_date_method(start, MODIFIED_SUCCEEDING);
    end_date   = bus_date_method(end, MODIFIED_SUCCEEDING);
    if (start_date >= end_date)
        return serror("Start >= End");

    if (elements & SE_FIXED)
    {
        err = interp_compounding(freqStr, &freq);
        if (err)
            return err;

        err = interp_basis(basisStr, &basis);
        if (err)
            return err;

        theo_date = end;
        act_date  = end_date;

        while (act_date > start_date)
        {
            theo_date = add_unit(theo_date, -12 / freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
            temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
            if (temp_date < start_date)
                temp_date = start_date;
            dates[k] = act_date;
            cpns[k]  = sign * fixed * coverage(temp_date, act_date, basis);
            k++;
            act_date = temp_date;
        }
    }
    if (elements & SE_SPREADS)
    {
        err = swp_f_get_ref_rate_details(refrateStr, &basis, &freq);
        if (err)
            return err;

        theo_date = end;
        act_date  = end_date;

        while (act_date > start_date)
        {
            theo_date = add_unit(theo_date, -12 / freq, SRT_MONTH, NO_BUSDAY_CONVENTION);
            temp_date = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
            if (temp_date < start_date)
                temp_date = start_date;
            dates[k] = act_date;
            tmp      = swp_f_spread(temp_date, act_date, refrateStr);
            if (tmp == SRT_SPREAD_ERROR)
                return serror("Error: swp_f_spread");
            cpns[k] = -sign * tmp * coverage(temp_date, act_date, basis);
            k++;
            act_date = temp_date;
        }
    }
    if (elements & SE_NOTIONALS)
    {
        dates[k] = start_date;
        cpns[k]  = -sign;
        k++;
        dates[k] = end_date;
        cpns[k]  = sign;
        k++;
    }

    err = indexx_ll(dates, idx, k);
    if (err)
        return err;

    j  = -1;
    dj = 0;
    for (i = 0; i < k; i++)
    {
        di = dates[idx[i]];
        if (di != dj)
        {
            j++;
            dates[idx[j]] = dj = di;
            cpns[idx[j]]       = cpns[idx[i]];
        }
        else
            cpns[idx[j]] += cpns[idx[i]];
    }
    *len    = j + 1;
    *times  = (double*)calloc(*len, sizeof(double));
    *cflows = (double*)calloc(*len, sizeof(double));
    if (!*times || !*cflows)
        return serror("Memory failure");

    for (i = 0; i < *len; i++)
    {
        (*times)[i]  = (dates[idx[i]] - today) * YEARS_IN_DAY;
        (*cflows)[i] = cpns[idx[i]];
    }
    return NULL;
}

Err CashFlows_Init(SCashFlows* g, int nex)
{
    g->nex  = nex;
    g->ex   = (double*)calloc(nex, sizeof(double));
    g->nmat = (int*)calloc(nex, sizeof(int));
    if (!g->ex || !g->nmat)
        return serror("Memory failure");
    g->mat = (double**)calloc(nex, sizeof(double*));
    if (g->mat)
        memset(g->mat, 0, nex * sizeof(double*));
    else
        return serror("Memory failure");
    g->cpn = (double**)calloc(nex, sizeof(double*));
    if (g->cpn)
        memset(g->cpn, 0, nex * sizeof(double*));
    else
        return serror("Memory failure");
    return NULL;
}

Err CashFlows_Free(SCashFlows* g)
{
    int i;

    if (g->mat)
        for (i = 0; i < g->nex; i++)
            free(g->mat[i]);
    if (g->cpn)
        for (i = 0; i < g->nex; i++)
            free(g->cpn[i]);
    free(g->mat);
    free(g->cpn);
    free(g->nmat);
    free(g->ex);
    memset(g, 0, sizeof(SCashFlows));

    return NULL;
}
