/* =================================================================================

   FILENAME: srt_h_cheybeta2fdynamics.h

   PURPOSE:  Give all the Equations required to descretise the model:
				- drift of r1 and r2
				- expectations of r1 and r2 at t+1 knowing t
				- local volatilities of r1 and r2
				- varaince of r1 and r2 at t+1 knowing t
				(and same for phi's and cross-phi)
   ================================================================================= */
   
#ifndef SRT_H_CHEYBETA2FDYNAMICS_H
#define SRT_H_CHEYBETA2FDYNAMICS_H

/* The drifts of the short rate components ( drift = Phi + CrossPhi - Lam * X )  */
Err srt_f_CheyBeta2f_shortrate_drifts(
			SrtIRMTmInf    *tminf,
			SrtSample      *cur_sample,
			long            und_index,
			double         *drift_x1,
			double         *drift_x2);

/* The Expected values of the short rate components ( E(Xt+dt|Xt) = Xt + drift * dt )  */
Err srt_f_CheyBeta2f_shortrate_forwards(
			SrtIRMTmInf    *tminf,
			SrtSample      *cur_sample,
			long            und_index,
			double          delta_t,
			double         *forward_x1,
			double         *forward_x2);

/* The Local Vols of the short rate components ( vol = sigma * (r ^ Beta) )  */
Err srt_f_CheyBeta2f_shortrate_localvols (
		  SrtIRMTmInf  *tminf,
		  SrtSample    *cur_sample,
		  int          und_index,
		  double       *localvol_x1,
		  double       *localvol_x2,
		  double       *correl_x1x2);

/* The Variances/Covariance of the short rate components ( Var[Xt+dt] = vol * vol * dt) */
Err srt_f_CheyBeta2f_shortrate_variances (
		  SrtIRMTmInf  *tminf,
		  SrtSample    *cur_sample,
		  int          und_index,
		  double       delta_t,
		  double       *variance_x1,
		  double       *variance_x2,
		  double       *covariance_x1x2);


/* The drifts of the Phi's variables( drift = sig * sig - (lam + lam) * Phi )      */
Err srt_f_CheyBeta2f_phi_drifts(
			SrtIRMTmInf    *tminf,
			SrtSample      *cur_sample,
			long            und_index,
			double         *drift_phi1,
			double         *drift_phi2,
			double         *drift_phi12);

/* The forwards of the Phi's variables( E(Phi t+1) = Phi t + Drift * delta_t )      */
Err srt_f_CheyBeta2f_phi_forwards(
			SrtIRMTmInf    *tminf,
			SrtSample      *cur_sample,
			long            und_index,
			double          delta_t,
			double         *forward_phi1,
			double         *forward_phi2,
			double         *forward_phi12);

/* Populates fwd_sample with the forwardvalues of the state varaibles, computed
   taking cur_sample as a starting point from cur_step to the next one */
Err srt_f_CheyBeta2f_forward_sample(
			SrtStpPtr       cur_step,
			SrtSample      *cur_sample,
			SrtSample      *fwd_sample, 
			long            und_index);

/* Evolves a Sample from Step to Step->Next given two uncorrelated Gaussian
   Random Numbers */
Err srt_f_CheyBeta2f_evolve_sample(
			SrtStpPtr       cur_step,
			SrtSample      *cur_sample,
			long            und_index,
			double          rand1,
			double          rand2,
			SrtSample      *next_sample);



/* --------------------------------------------------------------------------------- */
/*          FOR A GENERIC CHEYETTE BETA 2 FACTOR MODEL                              */
/* --------------------------------------------------------------------------------- */
/* The Local Vols of the short rate components ( vol = sigma * (r ^ Beta) )  */
Err srt_f_GenericCheyBeta2f_shortrate_localvols (
		  SrtIRMTmInf  *tminf,
		  SrtSample    *cur_sample,
		  int           und_index,
		  SrtMdlType    mdl_type,
		  double       *localvol_x1,
		  double       *localvol_x2);
#endif