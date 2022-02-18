/* ------------------------------------------------------------------------
   FILENAME:  	srt_f_ts_eq_fct.c

   PURPOSE:     Function to work with an EQ TermStruct:
   ------------------------------------------------------------------------ */
#include "SrtUndUtils.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_fwdcurve.h"
#include "srt_h_repocurve.h"
#include "srt_h_ts_eq.h"
#include "srtaccess.h"

double** sig_matrix;
#define MAX_ITER 10
#define VOL_TOL 1e-4

static SrtListAtom* next_element_with_name(SrtListAtom* lst, char* name, int skip)
{
    /* Check that the element exists else nothing to do */
    if (!lst)
        return NULL;

    /* If the skip flag is set to "NO" and points to the right element: return it */
    if (!skip && !strcmp(lst->element->name, name))
        return lst;

    /* Go from next to next until the right name is found or end of list */
    lst = lst->next;
    while (lst && strcmp(lst->element->name, name) != 0)
        lst = lst->next;

    /* Return the element (even if NULL) */
    return lst;

} /* END SrtListAtom *next_element_with_name (... ) */
/* ------------------------------------------------------------------------ */

/* Finds the Value of the loval volatility at one given point in time */

double find_eq_sig(double time, TermStruct* l)
{
    SrtLst* ls;

    ls = l->head;

    while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;

    if (ls != NULL)
    {
        return ((EquityTermStructVal*)ls->element->val.pval)->sig;
    }
    else
    {
        return ((EquityTermStructVal*)l->tail->element->val.pval)->sig;
    }
}

/* ---------------------------------------------------------------------------- */

double eq_cum_vol_func(double time, SRT_Boolean sigsq, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               val, cum_vol_time = 0.0, prev_time;

    ls = l->head;

    if (l->head == l->tail)
    {
        val = ((EquityTermStructVal*)ls->element->val.pval)->sig;
        if (sigsq == SRT_YES)
            val *= val;
        return (val * time);
    }
    /* first case : time is before the first sigma bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        prev_time = 0.0;
        val       = ((EquityTermStructVal*)ls->element->val.pval)->sig;
    }
    /* general case */
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            val   = tsval->sig;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            val   = ((EquityTermStructVal*)ls->element->val.pval)->sig;
        }

        prev_time    = tsval->time;
        cum_vol_time = (sigsq == SRT_YES ? tsval->int_sig2_dt : tsval->int_sig_dt);
    }
    if (sigsq == SRT_YES)
        val *= val;
    cum_vol_time += val * (time - prev_time);

    return cum_vol_time;
}

/* Stochastic Rate & Volatility Gamma Smile Model functions */

double find_eq_omega(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               omega;

    ls = l->head;

    if (l->head == l->tail)
    {
        omega = ((EquityTermStructVal*)ls->element->val.pval)->omega;
        return omega;
    }
    /* first case : time is before the first omega bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        omega = ((EquityTermStructVal*)ls->element->val.pval)->omega;
        return omega;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            omega = tsval->omega;
            return omega;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            omega = ((EquityTermStructVal*)ls->element->val.pval)->omega;
            return omega;
        }
    }
}

double find_eq_beta(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               beta;

    ls = l->head;

    if (l->head == l->tail)
    {
        beta = ((EquityTermStructVal*)ls->element->val.pval)->beta;
        return beta;
    }
    /* first case : time is before the first beta bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        beta = ((EquityTermStructVal*)ls->element->val.pval)->beta;
        return beta;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            beta  = tsval->beta;
            return beta;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            beta  = ((EquityTermStructVal*)ls->element->val.pval)->beta;
            return beta;
        }
    }
}

double find_eq_gamma(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               gamma;

    ls = l->head;

    if (l->head == l->tail)
    {
        gamma = ((EquityTermStructVal*)ls->element->val.pval)->gamma;
        return gamma;
    }
    /* first case : time is before the first gamma bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        gamma = ((EquityTermStructVal*)ls->element->val.pval)->gamma;
        return gamma;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            gamma = tsval->gamma;
            return gamma;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            gamma = ((EquityTermStructVal*)ls->element->val.pval)->gamma;
            return gamma;
        }
    }
}

double find_eq_basevol(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               basevol;

    ls = l->head;

    if (l->head == l->tail)
    {
        basevol = ((EquityTermStructVal*)ls->element->val.pval)->basevol;
        return basevol;
    }
    /* first case : time is before the first basevol bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        basevol = ((EquityTermStructVal*)ls->element->val.pval)->basevol;
        return basevol;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval   = (EquityTermStructVal*)l->tail->element->val.pval;
            basevol = tsval->basevol;
            return basevol;
        }
        else
        {
            tsval   = (EquityTermStructVal*)ls->previous->element->val.pval;
            basevol = ((EquityTermStructVal*)ls->element->val.pval)->basevol;
            return basevol;
        }
    }
}

double find_eq_voldrift(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               voldrift;

    ls = l->head;

    if (l->head == l->tail)
    {
        voldrift = ((EquityTermStructVal*)ls->element->val.pval)->voldrift;
        return voldrift;
    }
    /* first case : time is before the first voldrift bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        voldrift = ((EquityTermStructVal*)ls->element->val.pval)->voldrift;
        return voldrift;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval    = (EquityTermStructVal*)l->tail->element->val.pval;
            voldrift = tsval->voldrift;
            return voldrift;
        }
        else
        {
            tsval    = (EquityTermStructVal*)ls->previous->element->val.pval;
            voldrift = ((EquityTermStructVal*)ls->element->val.pval)->voldrift;
            return voldrift;
        }
    }
}

double find_eq_vovol(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               vovol;

    ls = l->head;

    if (l->head == l->tail)
    {
        vovol = ((EquityTermStructVal*)ls->element->val.pval)->vovol;
        return vovol;
    }
    /* first case : time is before the first vovol bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        vovol = ((EquityTermStructVal*)ls->element->val.pval)->vovol;
        return vovol;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            vovol = tsval->vovol;
            return vovol;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            vovol = ((EquityTermStructVal*)ls->element->val.pval)->vovol;
            return vovol;
        }
    }
}

double find_eq_rho(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal* tsval;
    double               rho;

    ls = l->head;

    if (l->head == l->tail)
    {
        rho = ((EquityTermStructVal*)ls->element->val.pval)->rho_spot_vol;
        return rho;
    }
    /* first case : time is before the first rho bucket time */
    if (((EquityTermStructVal*)ls->element->val.pval)->time >= time)
    {
        rho = ((EquityTermStructVal*)ls->element->val.pval)->rho_spot_vol;
        return rho;
    }
    else
    {
        while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
            ls = ls->next;

        if (ls == NULL)
        {
            tsval = (EquityTermStructVal*)l->tail->element->val.pval;
            rho   = tsval->rho_spot_vol;
            return rho;
        }
        else
        {
            tsval = (EquityTermStructVal*)ls->previous->element->val.pval;
            rho   = ((EquityTermStructVal*)ls->element->val.pval)->rho_spot_vol;
            return rho;
        }
    }
}

/* Function to compute : int(s = 0, s = t) (rho(s)*sig_eq(s)*sig_ir(s)/F(s)) */
double V_ir_e_func(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal *tsval, *tsval_p;
    double               ir_sig, tau, F;
    double               res;

    ls = l->head;

    while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (EquityTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (EquityTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    ir_sig = tsval->StochRate_Ts.sig;
    tau    = tsval->StochRate_Ts.tau;
    F      = tsval->StochRate_Ts.F;

    res = tsval->V_ir_e + tsval->rho_e_ir * tsval->sig * (ir_sig / F) * (exp(time / tau) - 1) * tau;

    return res;
}
/* Function to compute : int(s = 0, s = t) (rho(s)*sig_eq(s)*sig_ir(s)*Psi(s)/F(s)) */
double W_ir_e_func(double time, TermStruct* l)
{
    SrtLst*              ls;
    EquityTermStructVal *tsval, *tsval_p;
    double               ir_sig, tau, Psi, F;
    double               res;

    ls = l->head;

    while ((ls != NULL) && (((EquityTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;
    if (ls == NULL)
        ls = l->tail;

    tsval = (EquityTermStructVal*)ls->element->val.pval;

    if (ls != l->head)
    {
        tsval_p = (EquityTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    ir_sig = tsval->StochRate_Ts.sig;
    tau    = tsval->StochRate_Ts.tau;
    F      = tsval->StochRate_Ts.F;
    Psi    = tsval->StochRate_Ts.Psi;

    res = tsval->W_ir_e +
          tsval->rho_e_ir * tsval->sig * ir_sig * (Psi / F) * (exp(time / tau) - 1) * tau +
          tsval->rho_e_ir * tsval->sig * ir_sig * tau * ((exp(time / tau) - 1) * tau - time);

    return res;
}

/* Function to get the volatility of a zero coupon bond */
Err srt_f_zc_bond_vol(
    double expiry,
    char*  und_name,

    double* zc_bond_vol)
{
    Err         err = NULL;
    SrtUndPtr   und;
    TermStruct* ts;

    /* get the underlying and its term structure */
    und = lookup_und(und_name);
    if (!und)
        return serror("can not get underlying");

    err = get_underlying_ts(und, &ts);
    if (err)
        return serror("can not get term structure");

    (*zc_bond_vol) = 2 * (Psi_func(expiry, ts) * L_func(expiry, ts) - O_func(expiry, ts)) / expiry;

    return NULL;
}

Err one_expiry_calib_func(
    char*  model,
    double spot,
    char*  und_name,
    char*  disc_und_name,
    char*  ccy,
    char*  dvd_name,
    char*  repoName,

    long     n_dividend,
    double** dividend,
    long     n_repo,
    double** repo,

    double expiry,
    int    ex_index,
    int    n_ex,

    double new_sig,

    SrtUndPtr* new_und,

    double* return_value)
{
    Err         err = NULL;
    SrtUndPtr   disc_und;
    SrtUndPtr   und;
    Date        date_today;
    TermStruct *ts, *disc_und_ts;
    char*       dividend_curve_name;
    char*       repo_curve_name;

    /* get today from the underlying */
    und = lookup_und(und_name);
    if (!und)
        return serror("can not find undelying");
    date_today = get_today_from_underlying(und);

    /* get the dividend curve name */
    dividend_curve_name = (String)malloc(strlen(dvd_name) + 1);
    strcpy(dividend_curve_name, dvd_name);

    /* get the repo curve name */
    repo_curve_name = (String)malloc(strlen(repoName) + 1);
    strcpy(repo_curve_name, repoName);

    /* destroy the underlying */
    err = srt_f_destroy_und(und_name);
    if (err)
        return err;

    /* thus re-initialise it */

    sig_matrix[1][ex_index - 1] = new_sig;

    err = SrtInitDividends(date_today, dividend_curve_name, ccy, n_dividend, 2, dividend);

    if (err)
        return err;

    err = SrtInitRepos(date_today, repo_curve_name, ccy, n_repo, 2, repo);

    if (err)
        return err;

    err = SrtInitEQUnd(
        und_name,
        model,
        spot,
        disc_und_name,
        dividend_curve_name,
        repo_curve_name,
        n_ex,
        2,
        sig_matrix,

        1.0,
        0.0,
        0.0,

        0.0,
        0.0,
        0.0);

    if (err)
        return err;

    (*new_und) = lookup_und(und_name);
    if (!(*new_und))
        return serror("can not find underlying");

    err = get_underlying_ts((*new_und), &ts);
    if (err)
        return serror("can not get term struct");

    disc_und = lookup_und(disc_und_name);
    if (!disc_und)
        return serror("can find underlying");

    err = get_underlying_ts(disc_und, &disc_und_ts);
    if (err)
        return serror("can not get term struct");

    (*return_value) =
        eq_cum_vol_func(expiry, SRT_YES, ts) -
        2 * (W_ir_e_func(expiry, ts) - Psi_func(expiry, disc_und_ts) * V_ir_e_func(expiry, ts));

    (*return_value) /= expiry;

    srt_free(dividend_curve_name);
    srt_free(repo_curve_name);

    return NULL;
}

/* Function to get the forward of an equity underlying */
Err srt_f_eq_forward(
    Date      date,
    SrtUndPtr und,

    char* yc_name,

    double* forward)
{
    Err    err = NULL;
    Date   date_today;
    double df_date;
    double time;
    double spot;

    String      dvd_crv_name, repo_crv_name;
    SrtCurvePtr dvd_crv, repo_crv;

    /* Get the date today from the underlying and the time to maturity */
    date_today = get_today_from_underlying(und);
    time       = (date - date_today) * YEARS_IN_DAY;

    /* Get the discount factor */
    df_date = swp_f_df(date_today, date, yc_name);

    /* Get the dividend & repo curve */
    dvd_crv_name = get_dividend_name_from_underlying(und);
    dvd_crv      = lookup_curve(dvd_crv_name);

    repo_crv_name = get_repo_name_from_underlying(und);
    repo_crv      = lookup_curve(repo_crv_name);

    /* Get the spot value */
    spot = get_spot_from_eqund(und);

    /* compute the forward */
    (*forward) = spot * srt_f_forward_from_fwdcrv(time, dvd_crv, repo_crv) / df_date;

    return err;
}

/* Function to get the implied volaility of an equity underlying */
Err srt_f_eq_implied_vol(double expiry, char* und_name, double* black_vol)
{
    Err          err = NULL;
    double       spot_eq_cum_vol;
    double       zc_bond_vol, covar;
    double       V_ir_e, W_ir_e;
    SrtModelType model_type;
    SrtUndPtr    und, ir_und;
    TermStruct * ts, *disc_und_ts;
    char*        disc_und_name;

    und = lookup_und(und_name);
    if (!und)
        return serror("can not get underying");

    model_type = get_mdltype_from_eqund(und);

    err = get_underlying_ts(und, &ts);
    if (err)
        return serror("can not get the equity term structure ");

    /* check if the model is a stoch rate one */
    if ((model_type != EQ_STOCH_RATES) && (model_type != EQ_STOCH_RATES_SRVGS))
    {
        spot_eq_cum_vol = eq_cum_vol_func(expiry, SRT_YES, ts);
        if (expiry > 0.0)
            *black_vol = sqrt(spot_eq_cum_vol / expiry);
        else
            *black_vol = find_eq_sig(0.0, ts);

        return NULL;
    }

    /* get the interest rate underlying term structure */
    disc_und_name = get_discname_from_underlying(und);
    ir_und        = lookup_und(disc_und_name);
    if (!ir_und)
        return serror(" can not get underlying from name");

    err = get_underlying_ts(ir_und, &disc_und_ts);
    if (err)
        return serror("can not get interest rate term structure");

    if (expiry <= 0.0)
    {
        *black_vol = find_eq_sig(0.0, ts);
        return NULL;
    }

    /* compute the spot equity volatility */
    spot_eq_cum_vol = 1.0 / expiry * eq_cum_vol_func(expiry, SRT_YES, ts);

    /* compute the zero-coupon bond volatility */
    zc_bond_vol = 2 *
                  (Psi_func(expiry, disc_und_ts) * L_func(expiry, disc_und_ts) -
                   O_func(expiry, disc_und_ts)) /
                  expiry;

    /* compute the covariance between the zc bond and the spot equity */
    V_ir_e = V_ir_e_func(expiry, ts);
    W_ir_e = W_ir_e_func(expiry, ts);
    covar  = (-Psi_func(expiry, disc_und_ts) * V_ir_e + W_ir_e) / expiry;

    /*compute the equity black & scholes volatility */
    (*black_vol) = sqrt(spot_eq_cum_vol + zc_bond_vol - 2 * covar);

    return NULL;
}

/* function to calibrate the local volatility of an equity under in a stoch rates environment */
Err srt_f_eq_calib(
    long n_ex,

    double* dates,
    double* vols,

    char*  disc_und_name,
    char** und_name)
{
    Err        err = NULL;
    SrtUndPtr  und;
    Date       date_today;
    long       k, l;
    SrtMdlType model_type;
    SrtMdlDim  model_dim;
    SrtUndPtr  new_und, disc_und;
    double*    newton_target;
    double     zc_bond_vol;
    double     spot;
    double     forward;

    double       nstop, a[3], b[3];
    char *       dvd_name, *yc_name, *repo_name;
    SrtCurvePtr  dvd_crv, repo_crv, yc_crv;
    double **    dvd_matrix, **repo_matrix;
    SrtListAtom* top;

    SrtRepoObj *repo_obj, *dvd_obj;
    char*       ccy;
    long        n_dvd, n_repo;

    /* get the sort underlying ptr and its model type */
    und = lookup_und((*und_name));
    if (!und)
        return serror("can not get underlying");
    model_type = get_mdltype_from_eqund(und);

    /* verify that the underlying is of EQ_STOCH_RATES or BLACK_SCHOLES type */
    if ((model_type != EQ_STOCH_RATES) && (model_type != BLACK_SCHOLES))
        return serror("underlying must be of EQ_STOCH_RATES or BLACK_SCHOLES type");

    /* get the model name */
    err = get_underlying_mdldim(und, &model_dim);
    if (err)
        return err;

    /* get the groth name the spot and the date today*/
    dvd_name = get_dividend_name_from_underlying(und);
    dvd_crv  = lookup_curve(dvd_name);
    dvd_obj  = get_dvdobj_from_dvdcrv(dvd_crv);

    top   = next_element_with_name(dvd_obj->head, "DVD", 0);
    n_dvd = 0;
    while (top)
    {
        n_dvd += 1;
        top = next_element_with_name(top, "DVD", 1);
    }

    dvd_matrix = dmatrix(0, 1, 0, n_dvd - 1);

    top = next_element_with_name(dvd_obj->head, "DVD", 0);
    for (l = 1; l <= n_dvd; l++)
    {
        dvd_matrix[0][l - 1] = top->element->key;
        dvd_matrix[1][l - 1] = top->element->val.dval;

        top = next_element_with_name(top, "DVD", 1);
    }

    repo_name = get_repo_name_from_underlying(und);
    repo_crv  = lookup_curve(repo_name);
    repo_obj  = get_repoobj_from_repocrv(repo_crv);

    top    = next_element_with_name(repo_obj->head, "REPO", 0);
    n_repo = 0;
    while (top)
    {
        n_repo += 1;
        top = next_element_with_name(top, "REPO", 1);
    }

    repo_matrix = dmatrix(0, 1, 0, n_repo - 1);

    top = next_element_with_name(repo_obj->head, "REPO", 0);
    for (l = 1; l <= n_repo; l++)
    {
        repo_matrix[0][l - 1] = top->element->key;
        repo_matrix[1][l - 1] = top->element->val.dval;

        top = next_element_with_name(top, "REPO", 1);
    }

    spot       = get_spot_from_eqund(und);
    date_today = get_today_from_underlying(und);

    sig_matrix = dmatrix(0, 1, 0, n_ex - 1);

    /* get the yc name and the ccy */
    disc_und = lookup_und(disc_und_name);
    if (!disc_und)
        return serror("can not get disc underlying");
    yc_name = get_discname_from_underlying(disc_und);

    yc_crv = lookup_curve(yc_name);
    if (yc_crv == NULL)
        return serror("fatal: can not get yield curve");
    ccy = get_curve_ccy(yc_crv);

    /* init. the sig curve */
    for (l = 1; l <= n_ex; l++)
    {
        sig_matrix[0][l - 1] = (double)dates[l - 1];
        sig_matrix[1][l - 1] = vols[l - 1];
    }

    /* define the newton target */
    newton_target = dvector(1, n_ex);

    for (l = 1; l <= n_ex; l++)
    {
        err = srt_f_zc_bond_vol(
            (dates[l - 1] - date_today) * YEARS_IN_DAY, disc_und_name, &zc_bond_vol);

        if (err)
            return err;

        newton_target[l] = (vols[l - 1] * vols[l - 1] - zc_bond_vol);
    }
    for (l = 1; l <= n_ex; l++)
    {
        nstop = 0.0;
        k     = 0;

        a[0] = vols[l - 1];
        if (l > 1)
        {
            dvd_name  = get_dividend_name_from_underlying(new_und);
            repo_name = get_repo_name_from_underlying(new_und);
        }

        err = one_expiry_calib_func(
            "EQ_STOCH_RATES",
            spot,
            *und_name,
            disc_und_name,
            ccy,
            dvd_name,
            repo_name,

            n_dvd,
            dvd_matrix,
            n_repo,
            repo_matrix,

            (dates[l - 1] - date_today) * YEARS_IN_DAY,
            l,
            n_ex,

            a[0],

            &new_und,
            &b[0]);
        if (err)
            return err;

        err = srt_f_eq_forward((long)dates[l - 1], new_und, yc_name, &forward);
        if (err)
            return err;

        a[1]      = vols[l - 1] + 1.0 / 100;
        dvd_name  = get_dividend_name_from_underlying(new_und);
        repo_name = get_repo_name_from_underlying(new_und);

        err = one_expiry_calib_func(
            "EQ_STOCH_RATES",
            spot,
            *und_name,
            disc_und_name,
            ccy,
            dvd_name,
            repo_name,

            n_dvd,
            dvd_matrix,
            n_repo,
            repo_matrix,

            (dates[l - 1] - date_today) * YEARS_IN_DAY,
            l,
            n_ex,

            a[1],

            &new_und,
            &b[1]);
        if (err)
            return err;

        a[2]      = vols[l - 1] + 2.0 / 100;
        dvd_name  = get_dividend_name_from_underlying(new_und);
        repo_name = get_repo_name_from_underlying(new_und);

        err = one_expiry_calib_func(
            "EQ_STOCH_RATES",
            spot,
            *und_name,
            disc_und_name,
            ccy,
            dvd_name,
            repo_name,

            n_dvd,
            dvd_matrix,
            n_repo,
            repo_matrix,

            (dates[l - 1] - date_today) * YEARS_IN_DAY,
            l,
            n_ex,

            a[2],

            &new_und,
            &b[2]);

        if (err)
            return err;

        while (nstop < 1 && k < MAX_ITER)
        {
            newton(newton_target[l], 2.0, a, b, &nstop);

            dvd_name  = get_dividend_name_from_underlying(new_und);
            repo_name = get_repo_name_from_underlying(new_und);

            err = one_expiry_calib_func(
                "EQ_STOCH_RATES",
                spot,
                *und_name,
                disc_und_name,
                ccy,
                dvd_name,
                repo_name,

                n_dvd,
                dvd_matrix,
                n_repo,
                repo_matrix,

                (dates[l - 1] - date_today) * YEARS_IN_DAY,
                l,
                n_ex,

                a[2],

                &new_und,
                &b[2]);

            if (err)
                return err;
            k++;
        }
    }

    if (sig_matrix)
        free_dmatrix(sig_matrix, 0, 1, 0, n_ex - 1);
    sig_matrix = NULL;
    if (dvd_matrix)
        free_dmatrix(dvd_matrix, 0, 1, 0, n_dvd - 1);
    dvd_matrix = NULL;
    if (repo_matrix)
        free_dmatrix(repo_matrix, 0, 1, 0, n_repo - 1);
    repo_matrix = NULL;
    if (newton_target)
        free_dvector(newton_target, 1, n_ex);
    newton_target = NULL;

    return err;
}

#undef MAX_ITER
#undef VOL_TOL
