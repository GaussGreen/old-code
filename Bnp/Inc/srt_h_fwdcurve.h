/* =======================================================================================
   
   FILENAME :       srt_h_fwdcurve.h
	 
   PURPOSE:         functions to access  information from the SrtCurvePtr object stored 
				    in the double linked list, refering to a SrtFwdDesc, containing a
					SrtFwdObj
   
   ======================================================================================= */
   
#ifndef SRT_H_FWDCURVE_H
#define SRT_H_FWDCURVE_H

/* ----------------------------------------------------------------------------------
    FUNCTIONS TO OPERATE ON A SrtCurvePtr THAT REFERS TO A FORWARD CURVE
   ---------------------------------------------------------------------------------- */

/* Computes the drift from start_time till end_time implied by the forward curve */
double srt_f_drift_from_fwdcrv(double start_time, double end_time, SrtCurvePtr crv);

/* Computes the implied disc (inverse of forward ratio) from start_time till end_time */
double srt_f_disc_from_fwdcrv(double start_time, double end_time, SrtCurvePtr crv);

/* Compute the forward from 0.0 to start_time */
double srt_f_forward_from_fwdcrv(double start_time, SrtCurvePtr dvdcrv,SrtCurvePtr repocrv);


/* ----------------------------------------------------------------------------
   THE FUNCTION TO CREATE, INITIALISE AND ADD A SrtFwdObj TO THE SrtCurveList
             (referred to through a SrtFwdDesc and a SrtCrvePtr ) 
   ---------------------------------------------------------------------------- */
Err srt_f_init_DividendCurve(	
			Date       today,
			char      *ccy,
			double   **forward_curve, 
			int        ncols,
			int        nrows,
			String     fwd_name);


#endif

