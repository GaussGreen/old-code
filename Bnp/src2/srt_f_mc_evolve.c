/* ====================================================================================
   FILENAME:  srt_f_mc_evolve.c

   PURPOSE:   description of all the functions used to evolve IR or Stock/FX models
              in a Monte-Carlo type discretisation

   DESCRIPTION: all the functions have the same type of name, and the same inputs
   ==================================================================================== */
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_betaetaM.h"
#include "srt_h_cheybeta2fdynamics.h"
#include "srt_h_cheybetadynamics.h"
#include "srt_h_mc_evolve.h"

#define MAX_STDEV 3.0
#define MIN_SPOT EPS

Err monte_CHE_1f_Euler_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    double       drift;
    double       variance;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    step = top;

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Initialization */

        sam_get(sam[0], und_index, PHI)        = 0.0;
        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, SHORT_RATE);
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_market */
        }

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];
            /* The Euler step to get the next short rate */

            drift = (sam_get(sam[i], und_index, PHI) -
                     tminf->ev.onef.lambda * sam_get(sam[i], und_index, STATEVAR) +
                     tminf->quanto_adjustment * tminf->ev.onef.sig *
                         sam_get(sam[i], und_index, SHORT_RATE)) *
                    step->delta_t;

            variance = tminf->ev.onef.sig * sam_get(sam[i], und_index, SHORT_RATE) *
                       step->sqrt_delta_t * rndm_vec[i + 1];

            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + drift + variance;

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            /* The finite difference step to get the next cum. vol, (A.19) */

            sam_get(sam[i + 1], und_index, PHI) =
                sam_get(sam[i], und_index, PHI) * (1 - 2 * tminf->ev.onef.lambda * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    tminf->ev.onef.sig2 * step->delta_t;

            /* The step to get 1/rolling money-market account,
                  product in (A.31)  (using renormalisation techniques)*/

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            {
                sam[i + 1].numeraire = exp(sam_get(sam[i], und_index, STATEVAR) * step->delta_t) *
                                       sam[i].numeraire; /* pv_money_mkt */
            }
            i++;
            step  = step->next;
            tminf = nxttminf;
        } /* end while step loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST ";
    }

    return 0;

} /* end monte_CHE_Euler_evolve fn definition */

/*----------------------------------------------------------------------------*/

/* function " monte_CHE_1f_Milshtein_evolve " describes the temporal evolution of
   the short rate "sam[i].short_rate" using STD. GAUSSIAN random numbers
   generated elsewhere, (and uses the MILSHTEIN approximation
   to the underlying SDE), and the cumulative volatility
   statistic  "sam[i].phi"  and the pathwise money-market
   account "sam[i].numeraire"

   Inputs:
   *rndm_vec : pointer to VECTOR containing sample path of
                Monte Carlo random vars, each row corresponds to one path.
   *step : pointer to linked list of structures which contains all preprocessed
            step (which is path-independent) utilized in the evolution equations.
   *sam : pointer to VECTOR of structure which the function fills with
           sample path of the short_rate, v_square(called phi), and
           the money_market account corresponding to THIS path of the
           short_rate.


*/

Err monte_CHE_1f_Milshtein_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    double       drift;
    double       variance;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    step  = top;
    tminf = step->tminf[und_index];

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Initialization */
        sam_get(sam[0], und_index, PHI)        = 0.0;
        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, SHORT_RATE);
        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0; /* pv money market */

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];
            /* The Milshtein step to get the next short rate */

            drift = (sam_get(sam[i], und_index, PHI) -
                     tminf->ev.onef.lambda * sam_get(sam[i], und_index, STATEVAR) +
                     tminf->quanto_adjustment * tminf->ev.onef.sig *
                         sam_get(sam[i], und_index, SHORT_RATE)) *
                        step->delta_t +
                    0.5 * tminf->ev.onef.sig2 * sam_get(sam[i], und_index, SHORT_RATE) *
                        step->delta_t * (rndm_vec[i + 1] * rndm_vec[i + 1] - 1);

            variance = tminf->ev.onef.sig * sam_get(sam[i], und_index, SHORT_RATE) *
                       step->sqrt_delta_t * rndm_vec[i + 1];

            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + drift + variance;

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            /* The finite difference step to get the next cum. vol, (A.19) */
            sam_get(sam[i + 1], und_index, PHI) =
                sam_get(sam[i], und_index, PHI) * (1 - 2 * tminf->ev.onef.lambda * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    tminf->ev.onef.sig2 * step->delta_t;

            /* The step to get 1/rolling money-market account, product in (A.31)*/
            if (und_index == sam[i + 1].numeraire_index)
                sam[i + 1].numeraire = exp(sam_get(sam[i], und_index, STATEVAR) * step->delta_t) *
                                       sam[i].numeraire; /* pv_money_mkt */

            i++;
            step  = step->next;
            tminf = nxttminf;

        } /* end while step->next!= NULL loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST ";
    }

    return 0;

} /* end monte_CHE_Milshtein_evolve fn definition */

/*----------------------------------------------------------------------------*/

Err monte_CHE_2f_Euler_evolve(double** rndm_mat, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    double       g[2]; /*Correlated Random numbers */
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    double       alpha, beta;

    step = top;

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Initialization */
        /* Please note that all along this function, r(t) = f(0,t) + X1 + X2 */

        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, X1) = sam_get(sam[0], und_index, X2) = 0.0;
        sam_get(sam[0], und_index, PHI1) = sam_get(sam[0], und_index, PHI2) = 0.0;
        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0; /* pv_money_mkt */

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            alpha = sqrt((1 + tminf->correl_x) / 2);
            beta  = sqrt((1 - tminf->correl_x) / 2);

            g[0] = alpha * rndm_mat[0][i + 1] + beta * rndm_mat[1][i + 1];

            /* The Euler step to get the next short rate */

            sam_get(sam[i + 1], und_index, X1) =
                sam_get(sam[i], und_index, X1) +
                (sam_get(sam[i], und_index, PHI1) + sam_get(sam[i], und_index, CROSSPHI) -
                 tminf->ev.twof[0].lambda * sam_get(sam[i], und_index, X1)) *
                    step->delta_t -
                tminf->ev.twof[0].sig * sam_get(sam[i], und_index, SHORT_RATE) *
                    step->sqrt_delta_t * g[0];

            g[1] = alpha * rndm_mat[0][i + 1] - beta * rndm_mat[1][i + 1];

            sam_get(sam[i + 1], und_index, X2) =
                sam_get(sam[i], und_index, X2) +
                (sam_get(sam[i], und_index, PHI2) + sam_get(sam[i], und_index, CROSSPHI) -
                 tminf->ev.twof[1].lambda * sam_get(sam[i], und_index, X2)) *
                    step->delta_t -
                tminf->ev.twof[1].sig * sam_get(sam[i], und_index, SHORT_RATE) *
                    step->sqrt_delta_t * g[1];

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(nxttminf->fwd_sam, und_index, F_0_t) + sam_get(sam[i + 1], und_index, X1) +
                sam_get(sam[i + 1], und_index, X2);

            /* The finite difference step to get the next cum. vol, */

            sam_get(sam[i + 1], und_index, PHI1) =
                sam_get(sam[i], und_index, PHI1) *
                    (1 - 2 * tminf->ev.twof[0].lambda * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    tminf->ev.twof[0].sig2 * step->delta_t;

            sam_get(sam[i + 1], und_index, PHI2) =
                sam_get(sam[i], und_index, PHI2) *
                    (1 - 2 * tminf->ev.twof[1].lambda * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    tminf->ev.twof[1].sig2 * step->delta_t;

            sam_get(sam[i + 1], und_index, CROSSPHI) =
                sam_get(sam[i], und_index, CROSSPHI) *
                    (1 - (tminf->ev.twof[0].lambda + tminf->ev.twof[1].lambda) * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    tminf->ev.twof[0].sig * tminf->ev.twof[1].sig * tminf->correl_x * step->delta_t;

            /* The step to get 1/rolling money-market account,
                  product in (A.31)  (using renormalisation techniques)*/

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            {
                sam[i + 1].numeraire =
                    sam[i].numeraire * exp(step->delta_t * (sam_get(sam[i], und_index, X1) +
                                                            sam_get(sam[i], und_index, X2)));
                /* pv_money_mkt */
            }
            i++;
            step  = step->next;
            tminf = nxttminf;
        } /* end while step loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST ";
    }
    return 0;

} /* END monte_CHE_2f_Euler_evolve definition */

/*----------------------------------------------------------------------------*/

/* LGM 1 F JUMPING - Last modifification done by R BENICHOU 10/12/98 */
Err monte_LGM_1f_Jumping_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    double       yld;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    step = top;

    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Initialization */

        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0; /* pv_money_mkt */

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            {
                sam_get(sam[i + 1], und_index, STATEVAR) =
                    sam_get(sam[i], und_index, STATEVAR) * nxttminf->ev.onef.F / tminf->ev.onef.F +
                    nxttminf->ev.onef.F * tminf->rf.onef.G *
                        (nxttminf->ev.onef.Psi - tminf->ev.onef.Psi) +
                    rndm_vec[i + 1] * nxttminf->rf.onef.stdev_x;
            }

            else if (und_index != sam[i + 1].numeraire_index)
            {
                sam_get(sam[i + 1], und_index, STATEVAR) =
                    sam_get(sam[i], und_index, STATEVAR) * nxttminf->ev.onef.F / tminf->ev.onef.F +
                    nxttminf->ev.onef.F * tminf->rf.onef.G *
                        (nxttminf->ev.onef.Psi - tminf->ev.onef.Psi) +
                    nxttminf->fxquanto_adjustment + nxttminf->ffdquanto_adjustment +
                    nxttminf->sfdquanto_adjustment + rndm_vec[i + 1] * nxttminf->rf.onef.stdev_x;
            }

            /* Rebuilds the Short Rate 	: in this case : f(0,t) + X */
            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            /* Reconstruction of the Numeraire if Domestic Underlying */
            if (und_index == sam[i + 1].numeraire_index)
            {
                /* Rebuilds the Zero Coupon from this date to the next one */
                Y_T_at_t_compute(1, &sam[i], &tminf->yp, &yld, und_index, ONE_FAC, LGM);

                /* Increment the Numeraire (with the renormalisation with f(0,t) */
                sam[i + 1].numeraire =
                    sam[i].numeraire *
                    exp(yld - sam_get(tminf->fwd_sam, und_index, F_0_t) * step->delta_t);
            }

            i++;
            tminf = nxttminf;
            step  = step->next;
        } /* end while(step->next != NULL) */
    }
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* END monte_LGM_1f_evolve */

/*----------------------------------------------------------------------------*/

Err monte_vasicek_1f_euler_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    long         i;

    step = top;
    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        sam_get(sam[0], und_index, SHORT_RATE) = tminf->vasicek_init_cond;

        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0;
        sam_get(sam[0], und_index, BT) = 1.0;

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i], und_index, SHORT_RATE) -
                tminf->ev.onef.lambda * sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t +
                tminf->ev.onef.lambda * tminf->mean_rev_level * step->delta_t +
                tminf->quanto_adjustment * tminf->ev.onef.sig * step->delta_t +
                tminf->ev.onef.sig * step->sqrt_delta_t * rndm_vec[i + 1];

            if (und_index == sam[i + 1].numeraire_index)
                sam[i + 1].numeraire =
                    sam[i].numeraire * exp(sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t);

            sam_get(sam[i + 1], und_index, BT) =
                sam_get(sam[i], und_index, BT) *
                exp(sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t);

            i++;
            tminf = nxttminf;
            step  = step->next;
        }
    }
    else
    {
        return "step == NULL at top of LIST ";
    }

    return NULL;

} /*End of monte_vasicek_1f_euler_evolve */

/* This is for LGM without jumping numeraire 	*/
Err monte_LGM_1f_Euler_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;

    step = top;
    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Initialization */
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);
        /* Here STATEVAR is r(t) - f(0,t) */
        sam_get(sam[0], und_index, STATEVAR) = 0.0;
        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0;

        sam_get(sam[0], und_index, BT) = 1.0;

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + nxttminf->rf.onef.H - tminf->rf.onef.H -
                tminf->ev.onef.lambda * step->delta_t *
                    (sam_get(sam[i], und_index, STATEVAR) - tminf->rf.onef.H) +
                tminf->quanto_adjustment * tminf->ev.onef.sig * step->delta_t +
                tminf->ev.onef.sig * step->sqrt_delta_t * rndm_vec[i + 1];

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            /* The step to get 1/rolling money-market account,
            product in (A.31)  (using renormalisation techniques)*/
            {
                sam[i + 1].numeraire = sam[i].numeraire * exp(sam_get(sam[i], und_index, STATEVAR) *
                                                              step->delta_t); /* pv_money_mkt */
            }

            sam_get(sam[i + 1], und_index, BT) =
                sam_get(sam[i], und_index, BT) *
                exp(sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t);

            i++;
            tminf = nxttminf;
            step  = step->next;
        } /* end while(step->next != NULL) */
    }
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* end monte_LGM_1f_Euler_evolve */

/*----------------------------------------------------------------------------*/

/* This is for NEWLGM without jumping numeraire 	*/

Err monte_NEWLGM_1f_Euler_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    double       vol_sqdt;

    step = top;
    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Initialization */
        /* we don't need anymore the SHORT_RATE !!!! we set it to STATEVAR */

        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(sam[0], und_index, STATEVAR);
        sam_get(sam[0], und_index, PHI)        = 0.0;

        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0;

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            vol_sqdt = sqrt(nxttminf->rf.onef.G - tminf->rf.onef.G);

            sam_get(sam[i + 1], und_index, STATEVAR) = sam_get(sam[i], und_index, STATEVAR) +
                                                       tminf->quanto_adjustment * vol_sqdt +
                                                       vol_sqdt * rndm_vec[i + 1];

            sam_get(sam[i + 1], und_index, PHI)        = nxttminf->rf.onef.G;
            sam_get(sam[i + 1], und_index, SHORT_RATE) = sam_get(sam[i + 1], und_index, STATEVAR);

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            /* The step to get 1/rolling money-market account,
            product in (A.31)  (using renormalisation techniques)*/
            {
                sam[i + 1].numeraire =
                    exp(-sam_get(sam[i + 1], und_index, STATEVAR) * nxttminf->rf.onef.H +
                        nxttminf->rf.onef.H * nxttminf->rf.onef.H *
                            sam_get(sam[i + 1], und_index, PHI)); /* pv_money_mkt */
            }

            i++;
            tminf = nxttminf;
            step  = step->next;
        } /* end while(step->next != NULL) */
    }
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* end monte_NEWLGM_1f_Euler_evolve */

/* This is for NEWLGM without jumping numeraire 	*/

Err monte_NEWCHEYBETA_1f_Euler_evolve(
    double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    double       vol_sqdt;
    double       temp_var;

    step = top;
    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Initialization */
        /* we don't need anymore the SHORT_RATE !!!! we set it to STATEVAR */

        sam_get(sam[0], und_index, STATEVAR) = temp_var = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE)          = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, PHI)                 = 0.0;

        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0;

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            vol_sqdt = sqrt(nxttminf->rf.onef.G - tminf->rf.onef.G) *
                       pow(fabs(sam_get(sam[i], und_index, SHORT_RATE)), tminf->ev.onef.beta);

            /*	* pow(fabs(sam_get(sam[i],und_index,STATEVAR)),tminf->ev.onef.beta);*/

            sam_get(sam[i + 1], und_index, STATEVAR) = temp_var + vol_sqdt * rndm_vec[i + 1];

            sam_get(sam[i + 1], und_index, PHI) =
                sam_get(sam[i], und_index, PHI) + vol_sqdt * vol_sqdt;

            temp_var = sam_get(sam[i + 1], und_index, STATEVAR);

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(nxttminf->fwd_sam, und_index, F_0_t) -
                nxttminf->ev.onef.F *
                    (temp_var - nxttminf->rf.onef.H * sam_get(sam[i + 1], und_index, PHI));

            if (und_index == sam[i + 1].numeraire_index) /*i.e. if Domestic Underlying */
            /* The step to get 1/rolling money-market account,
            product in (A.31)  (using renormalisation techniques)*/
            {
                sam[i + 1].numeraire =
                    exp(-temp_var * nxttminf->rf.onef.H +
                        nxttminf->rf.onef.H * nxttminf->rf.onef.H *
                            sam_get(sam[i + 1], und_index, PHI)); /* pv_money_mkt */
            }

            i++;
            tminf = nxttminf;
            step  = step->next;
        } /* end while(step->next != NULL) */
    }
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* end monte_NEWCHEYBETA_1f_Euler_evolve */

/*----------------------------------------------------------------------------*/

Err monte_LGM_2f_evolve(double** rndm_mat, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    double       alpha, beta;
    double       g[2]; /*Correlated Random numbers */
    double       yld;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;

    step = top;

    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Initialization */

        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, X1) = sam_get(sam[0], und_index, X2) = 0.0;
        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0; /* pv_money_mkt */

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            alpha = sqrt((1 + nxttminf->rf.twof[0][1].stdev_x) / 2);
            beta  = sqrt((1 - nxttminf->rf.twof[0][1].stdev_x) / 2);

            g[0] = alpha * rndm_mat[0][i + 1] + beta * rndm_mat[1][i + 1];

            sam_get(sam[i + 1], und_index, X1) =
                sam_get(sam[i], und_index, X1) * nxttminf->ev.twof[0].F / tminf->ev.twof[0].F +
                nxttminf->ev.twof[0].F * tminf->rf.twof[0][0].G *
                    (nxttminf->ev.twof[0].Psi - tminf->ev.twof[0].Psi) +
                nxttminf->ev.twof[0].F * tminf->rf.twof[0][1].G *
                    (nxttminf->ev.twof[1].Psi - tminf->ev.twof[1].Psi) +
                g[0] * nxttminf->rf.twof[0][0].stdev_x;

            g[1] = alpha * rndm_mat[0][i + 1] - beta * rndm_mat[1][i + 1];

            sam_get(sam[i + 1], und_index, X2) =
                sam_get(sam[i], und_index, X2) * nxttminf->ev.twof[1].F / tminf->ev.twof[1].F +
                nxttminf->ev.twof[1].F * tminf->rf.twof[1][1].G *
                    (nxttminf->ev.twof[1].Psi - tminf->ev.twof[1].Psi) +
                nxttminf->ev.twof[1].F * tminf->rf.twof[1][0].G *
                    (nxttminf->ev.twof[0].Psi - tminf->ev.twof[0].Psi) +
                g[1] * nxttminf->rf.twof[1][1].stdev_x;
            /*
            No quanto adjustment for the moment: might be very tricky....
            */

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(nxttminf->fwd_sam, und_index, F_0_t) + sam_get(sam[i + 1], und_index, X1) +
                sam_get(sam[i + 1], und_index, X2);

            /* Here, the Numeraire is the Jumping Numeraire : multiply current by B(Ti,Ti+1) | F(Ti)
             */
            Y_T_at_t_compute(1, &sam[i], &tminf->yp, &yld, und_index, TWO_FAC, LGM);

            if (und_index == sam[i + 1].numeraire_index) /* i.e. if domestic underlying */
            {
                /* Here again, do the Numeraire calculation taking the Renormalisation into account
                 */
                sam[i + 1].numeraire =
                    exp(yld - sam_get(tminf->fwd_sam, und_index, F_0_t) * step->delta_t) *
                    sam[i].numeraire; /* pv_money_mkt */
            }

            i++;
            tminf = nxttminf;
            step  = step->next;
        } /* end while(step->next != NULL) */
    }     /* end  if(step != NULL)*/
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* end monte_LGM_2f_evolve */

/* ----------------------------------------------------------------------------*/
/* OVE for POWER model */
/* ----------------------------------------------------------------------------
   Evolve dx = -lam * x * dt + sig * dWt
   The conditionnal expectations and variances of X(T)|X(t) are perfectly
   well known:
               E[ X_T | X_t ] = F(T)/F(t) * X_t
                           V[ X_T | X_t ] = F(T) * F(T) * [G(T) -G(t)]
   we can therefore use a "Jumping" technique (the numeraire is the same )
   ---------------------------------------------------------------------------- */

Err monte_BETAETA_1f_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    double       stddev;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    double tmp1;

    step = top;

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Initialization */
        sam_get(sam[0], und_index, STATEVAR) = 0.0;
        if (und_index == sam[0].numeraire_index) /* i.e. if domestic underlying */
        {
            sam[0].numeraire = 1.0; /* numeraire N(t)  */
        }

        /* Loop on time steps */

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];
            /* The step to get the State Variable at next time step */

            stddev = tminf->ev.onef.sig *
                     exp(log(fabs(1.0 + tminf->the_beta * sam_get(sam[i], und_index, STATEVAR))) *
                         tminf->eta) *
                     step->sqrt_delta_t;

            /*
                                    if   (1.0 + tminf->beta * sam_get(sam[i],und_index,STATEVAR) <
               0.) sprintf(msg,"1+beta*x has gone too small %f", (1.0 + tminf->beta *
               sam_get(sam[i],und_index,STATEVAR)));


                                    if   (1.0 + tminf->beta * sam_get(sam[i],und_index,STATEVAR) >=
               0.) sprintf(msg,"1+beta*x is ok %f", (1.0 + tminf->beta *
               sam_get(sam[i],und_index,STATEVAR)));

                                    smessage(msg);
            */
            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + stddev * rndm_vec[i + 1];

            /* The step to get the Numeraire at next step */
            if (und_index == sam[i + 1].numeraire_index) /* i.e. if Domestic Underlying */
            {
                /*	h_t_x = x_power_fct( sam_get(sam[i+1],und_index,STATEVAR)*
                   nxttminf->h_t_power, nxttminf->x_power );

                        sam[i+1].numeraire = exp ( nxttminf->A_t + h_t_x) ;
                */

                tmp1 = sam_get(sam[i + 1], und_index, STATEVAR);
                /*	sam[i+1].numeraire = exp ( nxttminf->A_t +
                                nxttminf->lambda_t * sam_get(sam[i+1],und_index,STATEVAR)) ; */
                sam[i + 1].numeraire = exp(-sam_get(sam[i + 1], und_index, STATEVAR));
            }

            i++;
            step  = step->next;
            tminf = nxttminf;

        } /* end while step loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST ";
    }
    return 0;

} /* end monte_BETAETA_1f_evolve fn definition */

/*----------------------------------------------------------------------------*/

Err monte_LGM_1f_STOCHVOL_evolve(
    double** rndm_mat, SrtStpPtr top, SrtSample* sam, double* Pt_by_Zt, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    double       drift;
    double       stdev;
    double       g[2];
    double       alpha;
    double       beta;

    step = top;

    if (step != NULL)
    {
        tminf = (SrtIRMTmInf*)step->tminf[und_index];
        i     = 0;

        /* Use alpha and beta to get symetrical CORRELATED values for the Brownians */
        alpha = sqrt((1 + tminf->rho) / 2);
        beta  = sqrt((1 - tminf->rho) / 2);

        /* Initialization */
        /* Here STATEVAR is r(t) - f(0,t) : r(t) = X(t) + f(0,t)  */
        sam_get(sam[0], und_index, PHI)        = 0.0;
        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);

        sam_get(sam[0], und_index, SIGMA) = tminf->ev.onef.sig;

        if (und_index == sam[0].numeraire_index)
            sam[0].numeraire = 1.0;

        while (step->next != NULL)
        {
            g[0] = (alpha * rndm_mat[0][i + 1] + beta * rndm_mat[1][i + 1]);
            g[1] = (alpha * rndm_mat[0][i + 1] - beta * rndm_mat[1][i + 1]);

            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            drift = (sam_get(sam[i], und_index, PHI) -
                     tminf->ev.onef.lambda * sam_get(sam[i], und_index, STATEVAR)) *
                    step->delta_t;

            stdev = sam_get(sam[i], und_index, SIGMA) * step->sqrt_delta_t;

            /* The finite difference step to get the next state var */
            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + drift + stdev * g[0];

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            /* The finite difference step to get the next cum. vol */
            sam_get(sam[i + 1], und_index, PHI) =
                sam_get(sam[i], und_index, PHI) * (1 - 2 * tminf->ev.onef.lambda * step->delta_t) +
                sam_get(sam[i], und_index, SIGMA) * sam_get(sam[i], und_index, SIGMA) *
                    step->delta_t;

            /* The finite difference step to get the next vol */
            sam_get(sam[i + 1], und_index, SIGMA) = sam_get(sam[i], und_index, SIGMA) *
                                                    nxttminf->ev.onef.sig / tminf->ev.onef.sig *
                                                    exp(-0.5 * tminf->vovol_sqr * step->delta_t +
                                                        tminf->vovol * g[1] * step->sqrt_delta_t);

            /* Volatility is bounded in order to avoid getting explosive rates*/
            sam_get(sam[i + 1], und_index, SIGMA) = DMIN(
                sam_get(sam[i + 1], und_index, SIGMA),
                nxttminf->ev.onef.sig * exp(2 * tminf->vovol));

            /* If domestic underlying, accrue the numeraire : money market */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire =
                    exp(sam_get(sam[i], und_index, STATEVAR) * step->delta_t) * sam[i].numeraire;
            }
            i++;

            tminf = nxttminf;
            step  = step->next;
        }
    } /* end if (step != NULL ) */
    else
    {
        return "step==NULL at top of LIST ";
    }

    return NULL;

} /* END of monte_LGM_1f_STOCHVOL_evolve */

/*----------------------------------------------------------------------------*/

Err monte_CHE_1f_STOCHVOL_evolve(
    double** rndm_mat, SrtStpPtr top, SrtSample* sam, double* Pt_by_Zt, int und_index)
{
    long         i;
    SrtStpPtr    step;
    double       drift;
    double       variance;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    /* AA added certain parameters */

    double g[2];
    double alpha;
    double beta;

    step = top;

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Use alpha and beta as symetrical CORRELATED values for the Brownians */
        alpha = sqrt((1 + tminf->rho) / 2);
        beta  = sqrt((1 - tminf->rho) / 2);

        /* Initialization */
        sam_get(sam[0], und_index, PHI)        = 0.0;
        sam_get(sam[0], und_index, STATEVAR)   = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE) = sam_get(tminf->fwd_sam, und_index, F_0_t);

        /* AA added this.... */
        sam_get(sam[0], und_index, SIGMA) = tminf->ev.onef.sig;

        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
        }

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            g[0] = (alpha * rndm_mat[0][i + 1] + beta * rndm_mat[1][i + 1]);
            g[1] = (alpha * rndm_mat[0][i + 1] - beta * rndm_mat[1][i + 1]);

            drift = (sam_get(sam[i], und_index, PHI) -
                     tminf->ev.onef.lambda * sam_get(sam[i], und_index, STATEVAR)) *
                    step->delta_t;

            variance = sam_get(sam[i], und_index, SIGMA) * sam_get(sam[i], und_index, SHORT_RATE) *
                       step->sqrt_delta_t * g[0];

            sam_get(sam[i + 1], und_index, SIGMA) = sam_get(sam[i], und_index, SIGMA) *
                                                    nxttminf->ev.onef.sig / tminf->ev.onef.sig *
                                                    exp(-0.5 * tminf->vovol_sqr * step->delta_t +
                                                        tminf->vovol * g[1] * step->sqrt_delta_t);

            /* AA - Volatility is bound in order to avoid getting explosive rates*/
            sam_get(sam[i + 1], und_index, SIGMA) = DMIN(
                sam_get(sam[i + 1], und_index, SIGMA),
                nxttminf->ev.onef.sig * exp(3 * tminf->vovol));

            /* The Euler step to get the next short rate */

            sam_get(sam[i + 1], und_index, STATEVAR) =
                sam_get(sam[i], und_index, STATEVAR) + drift + variance;

            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                sam_get(sam[i + 1], und_index, STATEVAR) +
                sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            /* The finite difference step to get the next cum. vol, (A.19) */

            sam_get(sam[i + 1], und_index, PHI) =
                sam_get(sam[i], und_index, PHI) * (1 - 2 * tminf->ev.onef.lambda * step->delta_t) +
                sam_get(sam[i], und_index, SHORT_RATE) * sam_get(sam[i], und_index, SHORT_RATE) *
                    sam_get(sam[i], und_index, SIGMA) * sam_get(sam[i], und_index, SIGMA) *
                    step->delta_t;

            /* If domestic underlying, accrue the pv money market as numeraire */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire =
                    exp(sam_get(sam[i], und_index, STATEVAR) * step->delta_t) * sam[i].numeraire;
            }
            i++;
            step = step->next;

            tminf = nxttminf;

        } /* end while step loop */

    } /* end if ( step!= NULL ) */
    else
    {
        return "step==NULL at top of LIST ";
    }

    return 0;

} /* END monte_CHE_1f_STOCHVOL_evolve */

/*----------------------------------------------------------------------------*/
Err monte_CHEYBETA_1f_STOCHVOL_evolve(
    double** rndm_mat, SrtStpPtr top, SrtSample* sam, double* Pt_by_Zt, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    SrtIRMTmInf* nxttminf;

    /* AA added certain parameters */

    double g[2];
    double alpha;
    double beta;
    double statevar;
    double volsqdt;
    double cheybeta;
    double volvar;
    double phi;
    int    numeraire = 0;

    step = top;

    if (step != NULL)
    {
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        /* Use alpha and beta as symetrical CORRELATED values for the Brownians */
        alpha = sqrt((1 + tminf->rho) / 2);
        beta  = sqrt((1 - tminf->rho) / 2);

        /* Initialization */
        sam_get(sam[0], und_index, PHI) = phi = 0.0;
        sam_get(sam[0], und_index, STATEVAR) = statevar = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE)          = sam_get(tminf->fwd_sam, und_index, F_0_t);
        volvar                                          = tminf->ev.onef.sig;

        /* AA added this.... */
        cheybeta = tminf->ev.onef.beta;

        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
            numeraire        = 1;
        }

        sam_get(sam[0], und_index, BT) = 1.0;

        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            g[0] = (alpha * rndm_mat[0][i + 1] + beta * rndm_mat[1][i + 1]);
            g[1] = (alpha * rndm_mat[0][i + 1] - beta * rndm_mat[1][i + 1]);

            volsqdt = volvar * step->sqrt_delta_t *
                      pow(fabs(sam_get(sam[i], und_index, SHORT_RATE)), tminf->ev.onef.beta);

            volsqdt = DMIN(volsqdt, 0.3 * step->sqrt_delta_t);

            statevar += (phi - tminf->ev.onef.lambda * statevar) * step->delta_t + volsqdt * g[0];

            phi = phi + volsqdt * volsqdt - 2 * tminf->ev.onef.lambda * phi * step->delta_t;

            sam_get(sam[i + 1], und_index, STATEVAR) = statevar;
            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                statevar + sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            sam_get(sam[i + 1], und_index, PHI) = phi;

            volvar *=
                (nxttminf->ev.onef.sig / tminf->ev.onef.sig) *
                exp(-tminf->vovol_sqr * step->delta_t + tminf->vovol * g[1] * step->sqrt_delta_t);

            /* AA - Volatility is bound in order to avoid getting explosive rates*/
            /*	sam_get(sam[i+1],und_index,SIGMA) =
                            DMIN(sam_get(sam[i+1],und_index,SIGMA),nxttminf->ev.onef.sig*
                            exp(tminf->vovol*step->sqrt_delta_t*sqrt((double)
               i+1)-tminf->vovol_sqr*step->delta_t*(i+1)));

            */
            /* If domestic underlying, accrue the pv money market as numeraire */
            if (numeraire)
            {
                sam[i + 1].numeraire = exp(statevar * step->delta_t) * sam[i].numeraire;
            }

            sam_get(sam[i + 1], und_index, BT) =
                sam_get(sam[i], und_index, BT) *
                exp(sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t);
            i++;
            step = step->next;

            tminf = nxttminf;

        } /* end while step loop */

    } /* end if ( step!= NULL ) */
    else
    {
        return "step==NULL at top of LIST ";
    }

    return 0;

} /* END monte_CHEYBETA_1f_STOCHVOL_evolve */

/*----------------------------------------------------------------------------*/

Err monte_BLACKSCHOLES_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtLogTmInf *tminf, *nxttminf;

    step = top;

    if (step != NULL)
    {
        /* Initialisation */
        tminf                            = step->tminf[und_index];
        sam_get(sam[0], und_index, SPOT) = tminf->init_fwd_val;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_mkt */
        }
        i = 0;
        /* Loop on all the time steps */
        while (step->next != NULL)
        {
            nxttminf = step->next->tminf[und_index];
            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) * nxttminf->init_fwd_val / tminf->init_fwd_val *
                exp((-tminf->int_sig2_dt / 2) +
                    (rndm_vec[i + 1] + tminf->quanto_adjustment) * tminf->sqrt_int_sig2_dt);

            /* Sets pv_money_mkt */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }

            i++;
            tminf = nxttminf;
            step  = step->next;
        }
    }
    else
    {
        return "step==NULL at top of LIST ";
    }
    return NULL;

} /* END monte_BLACKSCHOLES_evolve */

/* ---------------------------------------------------------------------------- */
/*----------------------------------------------------------------------------*/

Err monte_NORMALBS_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtLogTmInf *tminf, *nxttminf;

    step = top;

    if (step != NULL)
    {
        /* Initialisation */
        tminf                            = step->tminf[und_index];
        sam_get(sam[0], und_index, SPOT) = tminf->init_fwd_val;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_mkt */
        }
        i = 0;

        /* Loop on all the time steps */
        while (step->next != NULL)
        {
            nxttminf = step->next->tminf[und_index];
            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) * (1.0 + tminf->drift) +
                (rndm_vec[i + 1] + tminf->quanto_adjustment) * tminf->sqrt_int_sig2_dt;

            /* Sets pv_money_mkt if underlying is domestic */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }
            i++;
            tminf = nxttminf;
            step  = step->next;
        }
    }
    else
    {
        return "step==NULL at top of LIST ";
    }
    return NULL;

} /* END monte_NORMALBS_evolve */

/* ---------------------------------------------------------------------------- */
Err monte_EQ_STOCH_RATES_evolve(
    double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index, int dom_index)
{
    long         i;
    SrtStpPtr    step;
    SrtLogTmInf *tminf, *nxttminf;

    step = top;
    if (step != NULL)
    {
        tminf                            = step->tminf[und_index];
        sam_get(sam[0], und_index, SPOT) = tminf->init_fwd_val;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_mkt */
        }
        i = 0;

        while (step->next != NULL)
        {
            nxttminf = step->next->tminf[und_index];
            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) * nxttminf->init_fwd_val / tminf->init_fwd_val *
                exp(+(sam_get(sam[i], dom_index, STATEVAR)) * step->delta_t -
                    (0.5 * tminf->int_sig2_dt) + (rndm_vec[i + 1]) * tminf->sqrt_int_sig2_dt);
            /* Sets pv_money_mkt if underlying is domestic */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }
            i++;
            tminf = nxttminf;
            step  = step->next;
        }
    }
    else
    {
        return "INFO == NULL at top of LIST";
    }

    return NULL;
}

Err monte_EQ_STOCH_RATES_SRVGS(
    double*    rndm_vec_spot,
    double*    rndm_vec_vol,
    SrtStpPtr  top,
    SrtSample* sam,
    int        und_index,
    int        dom_index)
{
    long         i;
    SrtStpPtr    step;
    SrtLogTmInf *tminf, *nxttminf;
    double       spot, omega, beta, gamma, basevol, voldrift, vovol, nxt_time, g[2];
    double       int_sig2_dt, sqrt_int_sig2_dt, stdev;

    step = top;
    if (step != NULL)
    {
        tminf = step->tminf[und_index];
        spot = sam_get(sam[0], und_index, SPOT) = tminf->init_fwd_val;
        sam_get(sam[0], und_index, SIGMA)       = tminf->basevol;

        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_mkt */
        }
        i = 0;

        while (step->next != NULL)
        {
            /* correlate the brownian motions */
            g[0] = rndm_vec_spot[i + 1];
            g[1] = tminf->rho * rndm_vec_spot[i + 1] +
                   sqrt(1 - tminf->rho * tminf->rho) * rndm_vec_vol[i + 1];

            /* get the next step time information */
            nxttminf = step->next->tminf[und_index];

            omega    = tminf->omega;
            beta     = tminf->beta;
            gamma    = tminf->gamma;
            basevol  = tminf->basevol;
            voldrift = tminf->voldrift;
            vovol    = tminf->vovol;
            nxt_time = step->next->time;

            /* compute sqrt_int_sig2_dt & int_sig2_dt */
            stdev = log(sam_get(sam[i], und_index, SPOT) / spot) / (basevol * sqrt(nxt_time));
            stdev = omega * tanh(stdev / omega);

            /* compute the new SIGMA */
            sam_get(sam[i + 1], und_index, SIGMA) =
                sam_get(sam[i], und_index, SIGMA) + (nxttminf->basevol - tminf->basevol) +
                voldrift * step->delta_t + vovol * g[1] * step->sqrt_delta_t;

            sqrt_int_sig2_dt = sam_get(sam[i + 1], und_index, SIGMA) *
                               (1 + beta * stdev + gamma * stdev * stdev) * sqrt(step->delta_t);
            int_sig2_dt = sqrt_int_sig2_dt * sqrt_int_sig2_dt;

            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) * nxttminf->init_fwd_val / tminf->init_fwd_val *
                exp(+(sam_get(sam[i], dom_index, STATEVAR)) * step->delta_t - 0.5 * int_sig2_dt +
                    g[0] * sqrt_int_sig2_dt);

            /* Sets pv_money_mkt if underlying is domestic */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }
            i++;
            tminf = nxttminf;
            step  = step->next;
        }
    }
    else
    {
        return "INFO == NULL at top of LIST";
    }

    return NULL;
}

Err monte_FX_STOCH_RATES_evolve(
    double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index, int dom_index, int for_index)
{
    long         i;
    SrtStpPtr    step;
    SrtLogTmInf* tminf;

    step = top;
    if (step != NULL)
    {
        tminf                            = step->tminf[und_index];
        sam_get(sam[0], und_index, SPOT) = tminf->init_fwd_val;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0; /* pv_money_mkt */
        }
        i = 0;

        while (step->next != NULL)
        {
            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) *
                exp(+(sam_get(sam[i], dom_index, SHORT_RATE)) * step->delta_t -
                    (sam_get(sam[i], for_index, SHORT_RATE)) * step->delta_t -
                    (0.5 * tminf->int_sig2_dt) + (rndm_vec[i + 1]) * tminf->sqrt_int_sig2_dt);

            /* Sets pv_money_mkt if underlying is domestic */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }
            i++;
            step  = step->next;
            tminf = step->tminf[und_index];
        }
    }
    else
    {
        return "INFO == NULL at top of LIST";
    }

    return NULL;
}

/* ---------------------------------------------------------------------------- */
/* After having called LGM with Jumping Numeraire */

Err monte_FX_STOCH_RATES_Jumping_evolve(
    double* rndm_vect, SrtStpPtr top, SrtSample* sam, int und_index, int dom_index, int for_index)

{
    SrtStpPtr   step;
    SrtFXTmInf *tminf, *nexttminf;
    int         i;
    // FILE			*fp0,*fp1,*fp2;

    //	fp0 = fopen("D:\\WORKAREA\\DEBUG\\fx_stoch_rates_debug_rdm_vect.txt", "w");
    //	fp1 = fopen("D:\\WORKAREA\\DEBUG\\fx_stoch_rates_debug_dom_state_var.txt", "w");
    //	fp2 = fopen("D:\\WORKAREA\\DEBUG\\fx_stoch_rates_debug_for_state_var.txt", "w");

    step = top;
    if (step != NULL)
    {
        tminf                            = (SrtFXTmInf*)step->tminf[und_index];
        sam_get(sam[0], und_index, SPOT) = tminf->init_spot;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
        }
        i = 0;
        while (step->next != NULL)
        {
            tminf     = (SrtFXTmInf*)step->tminf[und_index];
            nexttminf = (SrtFXTmInf*)step->next->tminf[und_index];
            /* Rebuilds the Spot  */

            sam_get(sam[i + 1], und_index, SPOT) =
                sam_get(sam[i], und_index, SPOT) *
                exp(+sam_get(sam[i], dom_index, STATEVAR) * tminf->dom_lambda_dt) *
                exp(-sam_get(sam[i], for_index, STATEVAR) * tminf->for_lambda_dt) *
                exp(tminf->mean_ln_fx - 0.5 * tminf->int_sigx2_dt) *
                exp(sqrt(tminf->var_ln_fx) * rndm_vect[i + 1]);

            //			fprintf(fp0,"%.16f\n",rndm_vect[i+1]);
            //			fprintf(fp1,"%.16f\n",sam_get(sam[i+1],dom_index,STATEVAR));
            //			fprintf(fp2,"%.16f\n",sam_get(sam[i+1],for_index,STATEVAR));

            /* Move on to the next step */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire = 1.0;
            }

            i++;
            step = step->next;
        }
    }

    else
    {
        return "INFO == NULL at top of LIST";
    }

    //	fclose(fp0);
    //	fclose(fp1);
    //	fclose(fp2);

    return NULL;

} /* END Err monte_FX_STOCH_RATES_Jumping_evolve(...) */

/* ---------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
   Set sam[t].pv_money_mkt at each time step to 1 for a deterministic
   underlying.
   --------------------------------------------------------------------------- */
Err monte_DETERMINISTIC_evolve(SrtStpPtr top, SrtSample* sam)
{
    long      i = 0;
    SrtStpPtr step;

    step = top;

    if (step != NULL)
    {
        while (step != NULL)
        {
            sam[i].numeraire = 1.0; /* pv_money_mkt */
            step             = step->next;
            i++;
        }
    }
    else
        return "step==NULL at top of LIST ";

    return NULL;

} /* monte_DETERMINISTIC_evolve */

/*----------------------------------------------------------------------------*/

/* The two following functions simulate General Cheyette evolutions, in a one-factor
world (function 1) and two-factor world (function). They use the local volatility
functions defined in srt_f_ts_init.c and srt_f_ts_init_2f.c, which return the local
volatility as a function of t and r, according to the model. */

Err monte_CHEYBETA_1f_Euler_evolve(double* rndm_vec, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf *tminf, *nxttminf;
    double       beta, volsqdt, statevar, phi;
    int          numeraire = 0;

    step  = top;
    tminf = (SrtIRMTmInf*)step->tminf[und_index];
    beta  = tminf->ev.onef.beta;

    if (step != NULL)
    {
        /* Initialization */
        i     = 0;
        tminf = (SrtIRMTmInf*)step->tminf[und_index];

        sam_get(sam[0], und_index, PHI) = phi = 0.0;
        sam_get(sam[0], und_index, STATEVAR) = statevar = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE)          = sam_get(tminf->fwd_sam, und_index, F_0_t);
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
            numeraire        = 1;
        }

        sam_get(sam[0], und_index, BT) = 1.0;
        while (step->next != NULL)
        {
            nxttminf = (SrtIRMTmInf*)step->next->tminf[und_index];

            volsqdt = tminf->ev.onef.sig * step->sqrt_delta_t *
                      pow(fabs(sam_get(sam[i], und_index, SHORT_RATE)), tminf->ev.onef.beta);

            statevar += (phi - tminf->ev.onef.lambda * statevar) * step->delta_t +
                        volsqdt * rndm_vec[i + 1];

            phi = phi + volsqdt * volsqdt - 2 * phi * tminf->ev.onef.lambda * step->delta_t;

            sam_get(sam[i + 1], und_index, STATEVAR) = statevar;
            sam_get(sam[i + 1], und_index, SHORT_RATE) =
                statevar + sam_get(nxttminf->fwd_sam, und_index, F_0_t);

            sam_get(sam[i + 1], und_index, PHI) = phi;

            if (numeraire)
            {
                sam[i + 1].numeraire = exp(statevar * step->delta_t) * sam[i].numeraire;
            }

            sam_get(sam[i + 1], und_index, BT) =
                sam_get(sam[i], und_index, BT) *
                exp(sam_get(sam[i], und_index, SHORT_RATE) * step->delta_t);

            i++;
            step  = step->next;
            tminf = nxttminf;
        } /* end while step loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST";
    }

    return NULL;
} /* end monte_CHEYBETA_1f_Euler_evolve fn definition */

/*----------------------------------------------------------------------------*/

Err monte_CHEYBETA_2f_Euler_evolve(double** rndm_mat, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    Err          err;

    step = top;
    if (step != NULL)
    {
        /* Initialization */
        i                              = 0;
        tminf                          = (SrtIRMTmInf*)step->tminf[und_index];
        sam_get(sam[0], und_index, X1) = sam_get(sam[0], und_index, X2) = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE)   = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, PHI1)         = sam_get(sam[0], und_index, PHI2) =
            sam_get(sam[0], und_index, CROSSPHI) = 0.0;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
        }

        while (step->next != NULL)
        {
            /* Evolve the ith Sample and store the result in the next one */
            err = srt_f_CheyBeta2f_evolve_sample(
                step, &sam[i], und_index, rndm_mat[0][i + 1], rndm_mat[1][i + 1], &sam[i + 1]);

            /* Update the Numeraire (the Bank account unit)(minus the forward rate...)  */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire =
                    exp((sam_get(sam[i], und_index, X1) + sam_get(sam[i], und_index, X2)) *
                        step->delta_t) *
                    sam[i].numeraire;
            }

            /* Move on to the next step */
            i++;
            step = step->next;

        } /* end while step loop */

    } /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST";
    }

    return NULL;
} /* end monte_CHEYBETA_2f_Euler_evolve fn definition */

/*----------------------------------------------------------------------------*/

Err monte_MIXEDBETA_2f_Euler_evolve(double** rndm_mat, SrtStpPtr top, SrtSample* sam, int und_index)
{
    long         i;
    SrtStpPtr    step;
    SrtIRMTmInf* tminf;
    Err          err;

    step = top;
    if (step != NULL)
    {
        /* Initialization */
        i                              = 0;
        tminf                          = (SrtIRMTmInf*)step->tminf[und_index];
        sam_get(sam[0], und_index, X1) = sam_get(sam[0], und_index, X2) = 0.0;
        sam_get(sam[0], und_index, SHORT_RATE)   = sam_get(tminf->fwd_sam, und_index, F_0_t);
        sam_get(sam[0], und_index, PHI1)         = sam_get(sam[0], und_index, PHI2) =
            sam_get(sam[0], und_index, CROSSPHI) = 0.0;
        if (und_index == sam[0].numeraire_index)
        {
            sam[0].numeraire = 1.0;
        }

        while (step->next != NULL)
        {
            /* Evolve the ith Sample and store the result in the next one */
            err = srt_f_CheyBeta2f_evolve_sample(
                step, &sam[i], und_index, rndm_mat[0][i + 1], rndm_mat[1][i + 1], &sam[i + 1]);

            /* Update the Numeraire (the Bank account unit)(minus the forward rate...)  */
            if (und_index == sam[i + 1].numeraire_index)
            {
                sam[i + 1].numeraire =
                    exp((sam_get(sam[i], und_index, X1) + sam_get(sam[i], und_index, X2)) *
                        step->delta_t) *
                    sam[i].numeraire;
            }

            i++;
            step = step->next;
        } /* end while step loop */
    }     /* end if != NULL */
    else
    {
        return "step==NULL at top of LIST";
    }

    return NULL;

} /* end monte_MIXEDBETA_2f_Euler_evolve fn definition */

/*----------------------------------------------------------------------------*/

#undef MAX_STDEV
#undef MIN_SPOT