/* =============================================================================
   FILENAME: srt_f_cheybetadynamics.c

   PURPOSE:  Give all the functions required to descretise the Cheyette Beta 
             model (based on the diffusion equation):
				- drift of r
				- expectation of r at t+1 knowing t
				- local volatility of r
				- varaince of r at t+1 knowing t (and same for phi)
   ================================================================================= */
#include "srt_h_all.h"
#include "srt_h_cheybetadynamics.h"
#include "math.h"

/* ----------------------------------------------------------------------------- */

/* Find the local vol in a one factor Cheyette model */

Err srt_f_CheyBeta_local_vol (
	  SrtMdlType   mdl_type, 
	  SrtIRMTmInf  *tminf,
	  SrtSample    *cur_sample,
	  int          und_index,
	  double       *local_vol)
{
double	      sigma;
double        beta;
double        short_rate;
	
/* Extracts the volatility and the short rate */
	sigma = tminf->ev.onef.sig;
	short_rate = samptr_get(cur_sample, und_index, SHORT_RATE);

	switch (mdl_type)
	{
		case LGM:
			*local_vol = sigma;
		break;
		case CHEY:
			*local_vol = sigma * short_rate;
		break;
		case CHEY_BETA:
			beta = tminf->ev.onef.beta;
			*local_vol = sigma * pow (fabs (short_rate), beta);
		break;
		default:
			return serror ("Unknown model in srt_f_one_fac_find_local_vol");
		break;
	}

/* Returns a success string */	
	return NULL;
}

/* ------------------------------------------------------------------------------- */

/* Evolves deterministically the short rate and Phi, and stores it in drift_sam */

Err  srt_f_CheyBeta_drift_at_sam		
					 (
		 SrtMdlType     mdl_type,
		 SrtStpPtr		stp, 
		 SrtSample		*cur_sam, 
		 SrtSample		*drift_sam, 
		 int			und_index
					 )

{
double       vol;
double       mu;
double       dphi;  
SrtIRMTmInf  *tminf;
SrtIRMTmInf  *nxttminf;
Err          err             = NULL;


	if (!stp->next)
	{
	  return serror("No more time step to compute the drift");
	}

	tminf    = stp->tminf[und_index];
	nxttminf = stp->next->tminf[und_index];


/* The STATEVAR is X = r - f(0,t) */
	samptr_get (cur_sam, 0, STATEVAR) = samptr_get (cur_sam, 0, SHORT_RATE)
									  - sam_get (tminf->fwd_sam, 0, F_0_t);
/* Drift of the STATEVAR */
	mu = - tminf->ev.onef.lambda * samptr_get (cur_sam, 0, STATEVAR)
	   + samptr_get (cur_sam, 0 ,PHI);  
	mu   *= stp->delta_t;

/* Gets the local volatility */
	err = srt_f_CheyBeta_local_vol( mdl_type, tminf, cur_sam, und_index, &vol);
	if (err)
		return err;

/* Drift of PHI */
	dphi = vol * vol - 2.0 * tminf->ev.onef.lambda * samptr_get (cur_sam,0,PHI);
	dphi *= stp->delta_t;

/* Update the new sample */
	samptr_get (drift_sam, 0, STATEVAR) = samptr_get (cur_sam, 0, STATEVAR) + mu;
	samptr_get (drift_sam, 0, SHORT_RATE) = sam_get (nxttminf->fwd_sam, 0, F_0_t)
										  + samptr_get (drift_sam, 0, STATEVAR);

	samptr_get (drift_sam, 0, PHI)  = samptr_get (cur_sam, 0, PHI) + dphi;

	return NULL;
}

/* ------------------------------------------------------------------------------- */

/* Evolves deterministically the short rate and Phi, and stores it in drift_sam */

Err  srt_f_CheyBeta_var_at_sam		
					 (
		 SrtMdlType     mdl_type,
		 SrtStpPtr		stp, 
		 SrtSample		*cur_sam, 
		 double         *var_at_sam, 
		 int			und_index
					 )
{
double        vol;
SrtIRMTmInf   *tminf;
Err           err = NULL;

/* Get the relevant time information for that step */
  tminf = (SrtIRMTmInf *)stp->tminf[und_index];

/* Gets the local volatility */
	err = srt_f_CheyBeta_local_vol( mdl_type, tminf, cur_sam, und_index, &vol);
	if (err)
		return err;

/* The variance is the local volatility times the time step */
	*var_at_sam = vol * vol * stp->delta_t;
  
/* Return a success string */
	return NULL;
}

/* ------------------------------------------------------------------------------- */
