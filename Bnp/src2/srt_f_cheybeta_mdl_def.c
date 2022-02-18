/* ==============================================================================

   FILE NAME:      srt_f_cheybeta_mdl_def.c

   OBJECT:         local vol Cheyette model definition

  =============================================================================== */

#include "math.h"
#include "srt_h_all.h"

/*	Read terms structures	*/

static SrtLst* get_next_sigma_date(SrtLst* ts)
{
    SrtLst* ls;
    ls = ts;
    while ((ls != NULL) && (((IrmTermStructVal*)ls->element->val.pval)->val_origin != SIGMA_DATE) &&
           (((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE))
    {
        ls = ls->next;
    }

    return ls;
}

static SrtLst* get_next_tau_date(SrtLst* ts)
{
    SrtLst* ls;

    ls = ts;
    while ((ls != NULL) && (((IrmTermStructVal*)ls->element->val.pval)->val_origin != TAU_DATE) &&
           (((IrmTermStructVal*)ls->element->val.pval)->val_origin != BOTH_DATE))
    {
        ls = ls->next;
    }

    return ls;
}

/*	Initialise, i.e. set pointers to NULL	*/
void chey_beta_mdl_init(CHEYBETA_MDL* mdl)
{
    mdl->yc_name      = NULL;
    mdl->num_sigma    = 0;
    mdl->sigma_times  = NULL;
    mdl->sigma        = NULL;
    mdl->num_lambda   = 0;
    mdl->lambda_times = NULL;
    mdl->lambda       = NULL;
    mdl->num_beta     = 0;
    mdl->beta_times   = NULL;
    mdl->beta         = NULL;
}

/*	Free	*/
void chey_beta_mdl_free(CHEYBETA_MDL* mdl)
{
    if (mdl->num_sigma)
    {
        free(mdl->sigma_times);
        free(mdl->sigma);
        mdl->num_sigma   = 0;
        mdl->sigma_times = NULL;
        mdl->sigma       = NULL;
    }

    if (mdl->num_lambda)
    {
        free(mdl->lambda_times);
        free(mdl->lambda);
        mdl->num_lambda   = 0;
        mdl->lambda_times = NULL;
        mdl->lambda       = NULL;
    }

    if (mdl->num_beta)
    {
        free(mdl->beta_times);
        free(mdl->beta);
        mdl->num_beta   = 0;
        mdl->beta_times = NULL;
        mdl->beta       = NULL;
    }
}

/*	Build from description	*/
void chey_beta_mdl_build(
    CHEYBETA_MDL* mdl,
    /*	Name of the yield curve	*/
    char* yc_name,
    /*	Number of sigma points	*/
    int num_sigma,
    /*	Sigma times	*/
    double* sigma_times,
    /*	Sigmas	*/
    double* sigma,
    /*	Number of tau points	*/
    int num_lambda,
    /*	Lambda times	*/
    double* lambda_times,
    /*	Lambdas	*/
    double* lambda,
    /*	Beta	*/
    double beta)
{
    SrtCurvePtr yldcrv = lookup_curve(yc_name);
    mdl->today_date    = (Date)get_clcndate_from_curve(yldcrv);

    mdl->yc_name = yc_name;

    mdl->num_sigma   = num_sigma;
    mdl->sigma_times = (double*)calloc(num_sigma, sizeof(double));
    memcpy(mdl->sigma_times, sigma_times, num_sigma * sizeof(double));
    mdl->sigma = (double*)calloc(num_sigma, sizeof(double));
    memcpy(mdl->sigma, sigma, num_sigma * sizeof(double));

    mdl->num_lambda   = num_lambda;
    mdl->lambda_times = (double*)calloc(num_lambda, sizeof(double));
    memcpy(mdl->lambda_times, lambda_times, num_lambda * sizeof(double));
    mdl->lambda = (double*)calloc(num_lambda, sizeof(double));
    memcpy(mdl->lambda, lambda, num_lambda * sizeof(double));

    mdl->num_beta      = 1;
    mdl->beta_times    = (double*)calloc(1, sizeof(double));
    mdl->beta_times[0] = 1.0;
    mdl->beta          = (double*)calloc(1, sizeof(double));
    mdl->beta[0]       = beta;
}

/*	Build from underlying	*/
void chey_beta_mdl_build_from_und(CHEYBETA_MDL* mdl, SrtUndPtr und)
{
    /*	Get Cheyette beta model info from underlying	*/
    SrtIrDesc*  ird = (SrtIrDesc*)(und->spec_desc);
    TermStruct* ts  = (TermStruct*)(ird->ts);

    /*	To read from the term structures	*/
    SrtLst*           ls;
    IrmTermStructVal* tsv;
    int               i;

    /*	Get today	*/
    mdl->today_date = get_today_from_underlying(und);

    /*	Get yc name	*/
    mdl->yc_name = ird->yc_name;

    /*	Calculate number of sigmas	*/
    mdl->num_sigma = 0;
    ls             = get_next_sigma_date(ts->head);
    while (ls)
    {
        mdl->num_sigma++;
        ls = ls->next;
        ls = get_next_sigma_date(ls);
    }

    /*	Allocate and fill sigmas and betas	*/
    mdl->num_beta = mdl->num_sigma;

    mdl->sigma_times = (double*)calloc(mdl->num_sigma, sizeof(double));
    mdl->sigma       = (double*)calloc(mdl->num_sigma, sizeof(double));
    mdl->beta_times  = (double*)calloc(mdl->num_sigma, sizeof(double));
    mdl->beta        = (double*)calloc(mdl->num_sigma, sizeof(double));

    for (i = 0, ls = get_next_sigma_date(ts->head); i < mdl->num_sigma;
         i++, ls = ls->next, ls = get_next_sigma_date(ls))
    {
        tsv                 = (IrmTermStructVal*)(ls->element->val.pval);
        mdl->sigma_times[i] = mdl->beta_times[i] = tsv->time;
        mdl->sigma[i]                            = tsv->sig;
        mdl->beta[i]                             = tsv->beta;
    }

    /*	Calculate number of lambdas	*/
    mdl->num_lambda = 0;
    ls              = get_next_tau_date(ts->head);
    while (ls)
    {
        mdl->num_lambda++;
        ls = ls->next;
        ls = get_next_tau_date(ls);
    }

    /*	Allocate and fill lambdas	*/
    mdl->lambda_times = (double*)calloc(mdl->num_lambda, sizeof(double));
    mdl->lambda       = (double*)calloc(mdl->num_lambda, sizeof(double));

    for (i = 0, ls = get_next_tau_date(ts->head); i < mdl->num_lambda;
         i++, ls = ls->next, ls = get_next_tau_date(ls))
    {
        tsv                  = (IrmTermStructVal*)(ls->element->val.pval);
        mdl->lambda_times[i] = tsv->time;
        mdl->lambda[i]       = 1.0 / tsv->tau;
    }
}

/*	Fill the diffusion parameters structure at a given date	*/
void chey_beta_mdl_param(
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Times	*/
    double t1,
    double t2,
    /*	Result	*/
    CHEYBETA_PARAM* param)
{
    int i;

    /*	Search relevant sigma index (such that sigma_times[i-1] < t1 <= sigma_times[i] */
    i = 0;
    while (i < mdl->num_sigma - 1 && mdl->sigma_times[i] < t1)
    {
        i++;
    }
    /*	Fill sigma and beta information	*/
    param->sigma = mdl->sigma[i];
    param->beta  = mdl->beta[i];

    /*	Search relevant lambda index (such that lambda_times[i-1] < t1 <= lambda_times[i] */
    i = 0;
    while (i < mdl->num_lambda - 1 && mdl->lambda_times[i] < t1)
    {
        i++;
    }
    /*	Fill lambda information	*/
    param->lambda = mdl->lambda[i];

    param->dt = t2 - t1;
}

/*	Get normal instantaneous variance at t	*/
double chey_beta_mdl_norm_var_at_t(
    /*	Param	*/
    CHEYBETA_PARAM* param,
    /*	Statevars	*/
    double x,
    double phi,
    /*	Forward rate	*/
    double f,
    /*	Max var	*/
    double maxvar,
    /*	Min var	*/
    double minvar)
{
    /*	Cheyette beta formula	*/
    double var = param->sigma * param->sigma * pow(fabs(f + x), 2 * param->beta);

    /*	Bounds	*/
    if (var > maxvar)
    {
        var = maxvar;
    }

    if (var < minvar)
    {
        var = minvar;
    }

    return var;
}

/*	Get drift, actually Et2 [ Xt2 - Xt1 / Ft1 ]	*/
double chey_beta_mdl_drift(
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Param	*/
    CHEYBETA_PARAM* param,
    /*	Statevars	*/
    double x,
    double phi)
{
    return (phi - param->lambda * x) * param->dt;
}

/*	Get var, actually Vt2 [ Xt2 - Xt1 / Ft1 ]	*/
double chey_beta_mdl_var(
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Param	*/
    CHEYBETA_PARAM* param,
    /*	Statevars	*/
    double x,
    double phi,
    /*	Norm var	*/
    double norm_var)
{
    return norm_var * param->dt;
}

/*	Get forward phi	*/
double chey_beta_mdl_phi(
    /*	Model	*/
    CHEYBETA_MDL* mdl,
    /*	Param	*/
    CHEYBETA_PARAM* param,
    /*	Statevars	*/
    double x,
    double phi,
    /*	Output from chey_beta_mdl_var	*/
    double var)
{
    return phi + var - 2 * param->lambda * phi * param->dt;
}
