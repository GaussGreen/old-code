/* =================================================================================

   FILENAME: srt_f_cheybeta2fdynamics.c

   PURPOSE:  Give all the Equations required to descretise the model:
                                - drift of r1 and r2
                                - expectations of r1 and r2 at t+1 knowing t
                                - local volatilities of r1 and r2
                                - varaince of r1 and r2 at t+1 knowing t
                                (and same for phi's and cross-phi)
   ================================================================================= */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_cheybeta2fdynamics.h"

/* ---------------------------------------------------------------------------------- */
/* The drifts of the short rate components ( drift = Phi + CrossPhi - Lam * X )  */

Err srt_f_CheyBeta2f_shortrate_drifts(
    SrtIRMTmInf* tminf, SrtSample* cur_sample, long und_index, double* drift_x1, double* drift_x2)
{
    Err err = NULL;

    *drift_x1 = samptr_get(cur_sample, und_index, PHI1) +
                samptr_get(cur_sample, und_index, CROSSPHI) -
                tminf->ev.twof[0].lambda * samptr_get(cur_sample, und_index, X1);

    *drift_x2 = samptr_get(cur_sample, 0, PHI2) + samptr_get(cur_sample, 0, CROSSPHI) -
                tminf->ev.twof[1].lambda * samptr_get(cur_sample, und_index, X2);

    /* Success message */
    return err;

} /* END Err srt_f_CheyBeta2f_shortrate_drifts(...) */

/* ----------------------------------------------------------------------------------- */

/* The Expected values of the short rate components ( E(Xt+dt|Xt) = Xt + drift * dt )  */

Err srt_f_CheyBeta2f_shortrate_forwards(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    long         und_index,
    double       delta_t,
    double*      forward_x1,
    double*      forward_x2)
{
    Err err = NULL;

    *forward_x1 =
        samptr_get(cur_sample, und_index, X1) +
        (samptr_get(cur_sample, und_index, PHI1) + samptr_get(cur_sample, und_index, CROSSPHI) -
         tminf->ev.twof[0].lambda * samptr_get(cur_sample, und_index, X1)) *
            delta_t;

    *forward_x2 =
        samptr_get(cur_sample, und_index, X2) +
        (samptr_get(cur_sample, und_index, PHI2) + samptr_get(cur_sample, und_index, CROSSPHI) -
         tminf->ev.twof[1].lambda * samptr_get(cur_sample, und_index, X2)) *
            delta_t;

    /* Success message */
    return err;
} /* END Err srt_f_CheyBeta2f_shortrate_forwards(...) */

/* ----------------------------------------------------------------------------------- */

/* The Local Vols of the short rate components ( vol = sigma * (r ^ Beta) )  */

Err srt_f_CheyBeta2f_shortrate_localvols(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    int          und_index,
    double*      localvol_x1,
    double*      localvol_x2,
    double*      correl_x1x2)
{
    double sigma1;
    double sigma2;
    double short_rate;
    double beta1;
    double beta2;
    Err    err = NULL;

    /* Extracts the sigma parameter and the short rate */
    sigma1 = tminf->ev.twof[0].sig;
    sigma2 = tminf->ev.twof[1].sig;

    short_rate = samptr_get(cur_sample, und_index, SHORT_RATE);

    /* Compute the local volatility */
    beta1        = tminf->ev.twof[0].beta;
    beta2        = tminf->ev.twof[1].beta;
    *localvol_x1 = sigma1 * pow(fabs(short_rate), beta1);
    *localvol_x2 = sigma2 * pow(fabs(short_rate), beta2);

    /* Get the correlation from the TmInf */
    *correl_x1x2 = tminf->correl_x;

    /* Returns a success string */
    return err;

} /* END Err srt_f_CheyBeta2f_shortrate_localvols (...) */

/* ----------------------------------------------------------------------------------- */

/* The Variances/Covariance of the short rate components ( Var[Xt+dt] = vol * vol * dt) */

Err srt_f_CheyBeta2f_shortrate_variances(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    int          und_index,
    double       delta_t,
    double*      variance_x1,
    double*      variance_x2,
    double*      covariance_x1x2)
{
    double localvol_x1;
    double localvol_x2;
    double correl_x1x2;
    Err    err = NULL;

    /* Compute the local volatilities and the local correlation */
    err = srt_f_CheyBeta2f_shortrate_localvols(
        tminf, cur_sample, und_index, &localvol_x1, &localvol_x2, &correl_x1x2);
    if (err)
        return err;

    /* Compute the variances */
    *variance_x1 = localvol_x1 * localvol_x1 * delta_t;
    *variance_x2 = localvol_x2 * localvol_x2 * delta_t;

    /* The Covariance... */
    *covariance_x1x2 = correl_x1x2 * localvol_x1 * localvol_x2 * delta_t;

    /* Returns a success string */
    return err;

} /* Err srt_f_CheyBeta2f_shortrate_variances (...) */

/* ------------------------------------------------------------------------------- */

/* The drifts of the Phi's variables( drift = sig * sig - (lam + lam) * Phi )      */

Err srt_f_CheyBeta2f_phi_drifts(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    long         und_index,
    double*      drift_phi1,
    double*      drift_phi2,
    double*      drift_phi12)
{
    double localvol_x1;
    double localvol_x2;
    double correl_x1x2;
    Err    err = NULL;

    /* Compute the local volatilities */
    err = srt_f_CheyBeta2f_shortrate_localvols(
        tminf, cur_sample, und_index, &localvol_x1, &localvol_x2, &correl_x1x2);
    if (err)
        return err;

    /* Compute the drifts for the straight Phi's */
    *drift_phi1 = localvol_x1 * localvol_x1 -
                  2 * tminf->ev.twof[0].lambda * samptr_get(cur_sample, und_index, PHI1);

    *drift_phi2 = localvol_x2 * localvol_x2 -
                  2 * tminf->ev.twof[1].lambda * samptr_get(cur_sample, und_index, PHI2);

    /* Compute the drifts for the Cross Phi (with the correlation between the Brownians) */
    *drift_phi12 = correl_x1x2 * localvol_x1 * localvol_x2 -
                   (tminf->ev.twof[0].lambda + tminf->ev.twof[1].lambda) *
                       samptr_get(cur_sample, und_index, CROSSPHI);

    /* Success message */
    return err;

} /* END Err srt_f_CheyBeta2f_phi_drifts(...) */

/* ------------------------------------------------------------------------------------ */
/* The forwards of the Phi's variables( E(Phi t+1) = Phi t + Drift * delta_t )      */

Err srt_f_CheyBeta2f_phi_forwards(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    long         und_index,
    double       delta_t,
    double*      forward_phi1,
    double*      forward_phi2,
    double*      forward_phi12)
{
    double variance_x1;
    double variance_x2;
    double covariance_x1x2;
    Err    err = NULL;

    /* Compute the variances for the delta_t time step  */
    err = srt_f_CheyBeta2f_shortrate_variances(
        tminf, cur_sample, und_index, delta_t, &variance_x1, &variance_x2, &covariance_x1x2);

    /* Compute the forwards for the straight Phi's */
    *forward_phi1 =
        samptr_get(cur_sample, und_index, PHI1) * (1.0 - 2.0 * tminf->ev.twof[0].lambda * delta_t) +
        variance_x1;

    *forward_phi2 =
        samptr_get(cur_sample, und_index, PHI2) * (1.0 - 2.0 * tminf->ev.twof[1].lambda * delta_t) +
        variance_x2;

    /* Compute the forward for the Cross Phi */
    *forward_phi12 = samptr_get(cur_sample, und_index, CROSSPHI) *
                         (1.0 - (tminf->ev.twof[0].lambda + tminf->ev.twof[1].lambda) * delta_t) +
                     covariance_x1x2;

    /* Success message */
    return err;

} /* END Err srt_f_CheyBeta2f_phi_forwards(...) */

/* --------------------------------------------------------------------------------- */
/* Populates fwd_sample with the forwardvalues of the state varaibles, computed
   taking cur_sample as a starting point from cur_step to the next one */

Err srt_f_CheyBeta2f_forward_sample(
    SrtStpPtr cur_step, SrtSample* cur_sample, SrtSample* fwd_sample, long und_index)
{
    SrtIRMTmInf* cur_tminf;
    SrtIRMTmInf* next_tminf;
    double       forward_x1;
    double       forward_x2;
    double       forward_phi1;
    double       forward_phi2;
    double       forward_phi12;
    Err          err = NULL;

    /* If there is not following time step: there is no point in doing this */
    if (!cur_step->next)
    {
        return NULL;
    }

    /* Extract the TmInfo's for this step and the next one */
    cur_tminf  = cur_step->tminf[und_index];
    next_tminf = cur_step->next->tminf[und_index];

    /* Compute the forwards of the short rate components */
    err = srt_f_CheyBeta2f_shortrate_forwards(
        cur_tminf, cur_sample, und_index, cur_step->delta_t, &forward_x1, &forward_x2);
    if (err)
        return err;

    /* Compute the forwards of the Phi varaibles */
    err = srt_f_CheyBeta2f_phi_forwards(
        cur_tminf,
        cur_sample,
        und_index,
        cur_step->delta_t,
        &forward_phi1,
        &forward_phi2,
        &forward_phi12);

    /* Repopulate the Sample at the next Step with these Forward values  */
    samptr_get(fwd_sample, und_index, X1)       = forward_x1;
    samptr_get(fwd_sample, und_index, X2)       = forward_x2;
    samptr_get(fwd_sample, und_index, PHI1)     = forward_phi1;
    samptr_get(fwd_sample, und_index, PHI2)     = forward_phi2;
    samptr_get(fwd_sample, und_index, CROSSPHI) = forward_phi12;

    /* Reconstruct the Short Rate from F(0,t) stored in the forward sam of next step */
    samptr_get(fwd_sample, und_index, SHORT_RATE) = sam_get(next_tminf->fwd_sam, und_index, F_0_t) +
                                                    samptr_get(fwd_sample, und_index, X1) +
                                                    samptr_get(fwd_sample, und_index, X2);

    /* Return a success message */
    return NULL;

} /* END srt_f_CheyBeta2f_forward_sample(...) */

/* --------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------- */

/* Evolves a Sample from Step to Step->Next given two uncorrelated Gaussian
   Random Numbers */

Err srt_f_CheyBeta2f_evolve_sample(
    SrtStpPtr  cur_step,
    SrtSample* cur_sample,
    long       und_index,
    double     rand1,
    double     rand2,
    SrtSample* next_sample)
{
    SrtIRMTmInf* cur_tminf;
    SrtIRMTmInf* next_tminf;
    double       coeff1;
    double       coeff2;
    double       g1;
    double       g2;
    double       fwd_x1;
    double       fwd_x2;
    double       vol_x1;
    double       vol_x2;
    double       corr_x1x2;
    double       fwd_phi1;
    double       fwd_phi2;
    double       fwd_phi12;
    Err          err = NULL;

    /* If there is not following time step: there is no point in doing this */
    if (!cur_step->next)
    {
        return NULL;
    }

    /* Extract the TmInfo's for this step and the next one */
    cur_tminf  = cur_step->tminf[und_index];
    next_tminf = cur_step->next->tminf[und_index];

    /* Compute the forwards of the short rate components */
    err = srt_f_CheyBeta2f_shortrate_forwards(
        cur_tminf, cur_sample, und_index, cur_step->delta_t, &fwd_x1, &fwd_x2);
    if (err)
        return err;

    /* Compute the local volatilities of the short rate components and their correlation */
    err = srt_f_CheyBeta2f_shortrate_localvols(
        cur_tminf, cur_sample, und_index, &vol_x1, &vol_x2, &corr_x1x2);
    if (err)
        return err;

    /* Compute the forwards of the Phi variables */
    err = srt_f_CheyBeta2f_phi_forwards(
        cur_tminf, cur_sample, und_index, cur_step->delta_t, &fwd_phi1, &fwd_phi2, &fwd_phi12);

    /* Compute the correlated Random numbers from the correlation */
    coeff1 = sqrt((1.0 + corr_x1x2) / 2);
    coeff2 = sqrt((1.0 - corr_x1x2) / 2);
    g1     = coeff1 * rand1 + coeff2 * rand2;
    g2     = coeff1 * rand1 - coeff2 * rand2;

    /* Populate the Next Sample with the new values : cur + drift  + random  */
    samptr_get(next_sample, und_index, X1)       = fwd_x1 + g1 * vol_x1 * cur_step->sqrt_delta_t;
    samptr_get(next_sample, und_index, X2)       = fwd_x2 + g2 * vol_x2 * cur_step->sqrt_delta_t;
    samptr_get(next_sample, und_index, PHI1)     = fwd_phi1;
    samptr_get(next_sample, und_index, PHI2)     = fwd_phi2;
    samptr_get(next_sample, und_index, CROSSPHI) = fwd_phi12;

    /* Reconstruct the Short Rate from F(0,t) stored in the forward sam of next step */
    samptr_get(next_sample, und_index, SHORT_RATE) =
        sam_get(next_tminf->fwd_sam, und_index, F_0_t) + samptr_get(next_sample, und_index, X1) +
        samptr_get(next_sample, und_index, X2);

    /* Return a success message */
    return NULL;

} /* END srt_f_CheyBeta2f_evolve_sample(...) */

/* --------------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------------- */
/*           FOR A GENERIC CHEYETTE BETA 2 FACTOR MODEL                              */
/* --------------------------------------------------------------------------------- */
/* The Local Vols of the short rate components ( vol = sigma * (r ^ Beta) )  */

Err srt_f_GenericCheyBeta2f_shortrate_localvols(
    SrtIRMTmInf* tminf,
    SrtSample*   cur_sample,
    int          und_index,
    SrtMdlType   mdl_type,
    double*      localvol_x1,
    double*      localvol_x2)
{
    double sigma1;
    double sigma2;
    double short_rate;
    double beta1;
    double beta2;
    Err    err = NULL;

    /* Extracts the sigma parameter and the short rate */
    sigma1 = tminf->ev.twof[0].sig;
    sigma2 = tminf->ev.twof[1].sig;

    short_rate = samptr_get(cur_sample, und_index, SHORT_RATE);

    /* Compute the local volatility */
    switch (mdl_type)
    {
    case LGM:
        *localvol_x1 = sigma1;
        *localvol_x2 = sigma2;
        break;
    case CHEY:
        *localvol_x1 = sigma1 * short_rate;
        *localvol_x2 = sigma2 * short_rate;
        break;
    case CHEY_BETA:
    case MIXED_BETA:
        beta1        = tminf->ev.twof[0].beta;
        beta2        = tminf->ev.twof[1].beta;
        *localvol_x1 = sigma1 * pow(fabs(short_rate), beta1);
        *localvol_x2 = sigma2 * pow(fabs(short_rate), beta2);
        break;
    default:
        return serror("Unknown model in srt_f_CheyBeta2f_shortrate_localvols");
        break;
    }

    /* Returns a success string */
    return err;

} /* END Err srt_f_GenericCheyBeta2f_shortrate_localvols (...) */

/* ----------------------------------------------------------------------------------- */