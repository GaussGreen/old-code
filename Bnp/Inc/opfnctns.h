/*******************************************************************************
**                                                                          
**              Copyright (c) 1993 PARIBAS Capital Markets Group            
**                                                                         
********************************************************************************
**                                                                          
**      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins              
**      SUB_SYSTEM:     OPT     Option Tools                                
**                                                                          
**      MODULE NAME:    OPTION_TOOLS                                        
**                                                                          
**      PURPOSE:        Import Header				            
**									    
**	INCLUDED BY:	option_tools				    
**			srt_f_bond_trin_compute                             
**			proba_tools
**			capfloor_compute                             
**                                                                          
**      AUTHORS:        Guillaume AMBLARD,                                  
**                      Remy KLAMMERS                                       
**                      Eric AULD                                           
**                      Alexandra POITE                                     
**
**	MODIFIED:	Ken LINTON 		MAY 94
**			Olivier VAN EYSEREN	JUN 95
**			Jasbir MALHI		JUL 95
**                                                                          
**      DATE:           1st October, 1992                                   
**                                                                          
**      DESCRIPTION:    Header for the library of general related option    
**			functions used in the SORT program applications     
*******************************************************************************/

#include        "utallhdr.h"
#include        <NUM_H_ALLHDR.H>

/* ========================================================================== */

#ifndef		OPFNCTNS_H
#define		OPFNCTNS_H

/* ========================================================================== */


/*******************************************************************************
**                      Macros, Typedefs and Constants                      
*******************************************************************************/

#define 	VOL_SHIFT 	    (double)	0.001
#define		VOL_ADD			(double)	0.001
#define		VOL_INIT		(double)	0.1
#define 	NULL_VOL 		(double) 	1.0e-12
#define 	STRIKE_SHIFT 	(double) 	0.01
#define 	ALPHA_SHIFT     (double)	0.01
#define 	BETA_SHIFT      (double)	0.01
#define 	RHO_SHIFT	    (double)	0.01
#define 	ZETA_SHIFT	    (double)	0.01
#define 	DELTA_SHIFT	    (double)	0.0001
#define 	SABRVOL_SHIFT	(double)	0.01
#define 	VANNA_SHIFT		(double)	0.01

#define 	MAX_ITER 					100
#define 	PREM_TOL 		(double) 	0.0000005
#define 	PREM_TOL_PCT 				0.001

#define 	MAX_RUNG 					100
#define 	MAX_STEP					1000
#define		DEFAULT_NUM_STEPS			400
#define		DEFAULT_NUM_TERMS			5
#define 	MAX_NUM_DATES 				50

#define		INFINITY					1E10

/******************************************************************************/

/* ========================================================================== */


/* global variable for options tools */
typedef struct{
	double vol_add; 
	double vol_shift;
	double sabrvol_shift;
	double strike_shift;
	double delta_shift;
	double vanna_shift;
	double theta_shift;
	double alpha_shift;
	double beta_shift;
	double rho_shift;
	double zeta_shift;
}GlobVarOptionTool;

extern GlobVarOptionTool GVOPT;

/*******************************************************************************
*                                                                           
* PURPOSE      	: iniatializes the shift parameters for greek computation
*                                                                                                          
*******************************************************************************/


Err srt_f_defaultgreekshift(void);
Err srt_f_setstrikeshift(double strike_shift);
double srt_f_getstrikeshift(void);
Err srt_f_setdeltashift(double delta_shift);
double srt_f_getdeltashift(void);

Err srt_f_setthetashift(double theta_shift);
double srt_f_getthetashift(void);

double srt_f_getvegashift(void);
Err srt_f_setvegashift(double vol_add);
double srt_f_getvannashift(void);
Err srt_f_setvannashift(double rates_shift);
double srt_f_getalphashift(void);
Err srt_f_setalphashift(double alpha_shift);
Err srt_f_setbetashift(double beta_shift);
double srt_f_getbetashift(void);
Err srt_f_setrhoshift(double rho_shift);
double srt_f_getrhoshift(void);
double srt_f_getsabrvolshift(void);
Err srt_f_setsabrvolshift(double vol_add);
Err srt_f_setzetashift(double zeta_shift);
double srt_f_getzetashift(void);

/*========================================================================== 

				PART 1
			STANDARD OPTIONS

==========================================================================*/

/* =============================================================================
  FUNCTION     : optmerton(...)  
 ============================================================================ */
Err optmerton(
			  double Fwd,
			  double Strike,
			  double sigma,
			  double U1,
			  double lambda1,
			  double U2,
			  double lambda2,
			  double maturity,
			  char *Logornorm,
			  SrtGreekType Greek,
			  double *answer
			  );
/* =============================================================================
  FUNCTION     : optmertonpremium(...)  
 ============================================================================ */
double optmertonpremium(
			  double dFwd,
			  double Strike,
			  double sigma,
			  double U1,
			  double lambda1,
			  double U2,
			  double lambda2,
			  double T,
			  char *Logornorm
			  );

/* =============================================================================
  FUNCTION     : optcalibmerton(...)  
 ============================================================================ */
 Err optcalibmerton(
					  double dFwd_loc,
					  int n_strikes_loc,
	                  double *strikes_loc,
					  double *market_vols,
					  double *param,
					  double maturity,
					  char   *Logornorm,
					  double *chisq, 
					  long   *ia,
					  double *calibvols
					  );
 
 
 
 /* =============================================================================
  FUNCTION     : optmertonsmile(...)  
 ============================================================================ */

 Err optmertonsmile (double dFwd,
					int n_strikes,
					double *strikes,
					double *parameters,
					double maturity,
					char   *Logornorm,
					double *impvols);

/* =============================================================================
  FUNCTION     : srt_f_optblksch(...)  
 ============================================================================ */

double 	srt_f_optblksch(	
		double 			fwd_price, 
		double 			k, 
		double 			vol, 
		double 			mat,
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		);

double srt_f_optblksch_accurate(	
		double 			fwd_price,
		double 			strike,
		double 			vol,
		double 			mat, 
		double 			disc, 
		SrtCallPutType 	call_put, 
		SrtGreekType	greek
		);

/* A function that checks the Black-Scholes inputs */
Err srt_f_optchkblksch(
		double 			fwd_price,
		double 			strike,
		double 			vol,
		double 			mat);

/* A Higher level interface function to call Black Scholes */
Err OptBlkSch(	
		double 				fwd_price,
		double 				strike,
		double 				vol,
		double 				mat, 
		double 	  			disc, 
		char               *call_put_str, 
		char               *greek_str,
		double             *result
		);


/* =============================================================================
  FUNCTION     : srt_f_optimpvol(...)
 ============================================================================ */
Err 		srt_f_optimpvol(		
		double 				premium, 
		double 				fwd_price, 
		double 				strike,
		double 				mat, 
		double 				disc, 
		SrtCallPutType		call_put, 
		SrtDiffusionType	lognormal_normal,
		double              *implied_vol
		);

Err 		srt_f_optimpvol_accurate(		
		double 				premium, 
		double 				fwd_price, 
		double 				strike,
		double 				mat, 
		double 				disc, 
		SrtCallPutType		call_put, 
		SrtDiffusionType	lognormal_normal,
		double              *implied_vol
		);

/* =============================================================================
  FUNCTION     : srt_f_optblkschbeta(...) 
 ============================================================================ */
double 	srt_f_optblkschbeta(
		double			fwd_price,	
		double 			strike, 
		double 			betavol,
		double 			mat,  
		double 			beta,
		double			disc,
		SrtCallPutType  call_put, 
		SrtGreekType    greek);

/* =============================================================================
  FUNCTION     : srt_f_optblkschbetaquick(...) 
 ============================================================================ */
double 	srt_f_optblkschbetaquick(
		double			forward,	
		double 			strike, 
		double 			betavol,
		double 			mat,  
		double 			beta,
		double			disc,
		SrtCallPutType  call_put, 
		SrtGreekType    greek);

/* =============================================================================
  FUNCTION     : srt_f_optblkschbetastochquick(...) 
 ============================================================================ */

double srt_f_optblkschbetastochquick(
		double            forward,							
		double            strike,
		double            maturity,
		double            betavol,
		double			  alpha,
		double			  beta,
		double            rho,
		double            disc,
		SrtDiffusionType  typeinput,
		SrtDiffusionType  log_norm, 
		SrtCallPutType    call_put, 
		SrtGreekType      greek);

/* =============================================================================
  FUNCTION     : srt_f_optimpvolbeta(...)
 ============================================================================ */
double 	srt_f_optimpvolbeta	(
			double 	premium,  
			double 	fwd_price, 
			double 	strike,
			double 	mat, 
			double  beta,
			double 	disc, 
			SrtCallPutType 		call_put
			);

/* =============================================================================
  FUNCTION     : srt_f_optblkvoltobetavol(...)
 ============================================================================ */
Err srt_f_optblkvoltobetavol(
		double            forward,							
		double            strike, 
		double            bsvol, 
		double            mat,  
		double            beta,
		double           *betavol);

/* =============================================================================
  FUNCTION     : srt_f_optblkvoltobetavol(...)
 ============================================================================ */
Err srt_f_optbetavoltoblkvol(
		double            forward,							
		double            strike, 
		double            betavol, 
		double            mat,  
		double            beta,
		double           *bsvol);

/* =============================================================================
  FUNCTION     : srt_f_optbetastochvoltoblkvol(...)
 ============================================================================ */
Err srt_f_optbetastochvoltoblkvol(
		double            forward,							
		double            strike, 
		double            betavol,
		double			  vovol,
		double			  rho,
		double            mat,  
		double            beta,
		double           *bsvol);
/* =============================================================================
  FUNCTION     : srt_f_optbetastochvoltoblkvolmeanreversion(...)
 ============================================================================ */
Err srt_f_optbetastochvoltoblkvolmeanreversion(
		double            forward,							
		double            strike, 
		double            betavol,
		double			  vovol,
		double			  rho, 
		double            tau,
		double            mr,
		double            mat,  
		double            beta,
		double           *bsvol);
/* =============================================================================
  FUNCTION     : srt_f_optbetastochvoltoblkvoldrift(...)
 ============================================================================ */
Err srt_f_optbetastochvoltoblkvoldrift(
		double            forward,							
		double            strike, 
		double            betavol,
		double			  vovol,
		double			  rho,
		double            mat,  
		double            beta,
		double           *bsvol);
/* =============================================================================
  FUNCTION     : srt_f_optbetastochvoltoblkvolpremium(...)
 ============================================================================ */

Err srt_f_optnewbetastochvoltoblkpremium(
		double            F,							
		double            K, 
		double            sigma,
		double            beta,
		double			  vovol,
		double			  rho,
		double            mat,  
		double           *premium);
/* =============================================================================
  FUNCTION     : srt_f_optsarbvol(...) & utilitties
 ============================================================================ */

Err	srt_f_optsarbvol(
		double            forward,							
		double            strike, 
		double            mat, 
		double            volinput,  
		double            alpha,
		double			  beta,
		double			  rho,
		SrtDiffusionType  input,
		SrtDiffusionType  output,
		double            *vol);

Err srt_f_optsabrvolq(
        double            forward,							
		double            strike, 
		double            maturity, 
		double            volFwd,  
		double            alpha,
		double            beta,
		double            corrFwdVol,
        double            volFx, 
        double            corrFwdFx, 
        double            corrVolFx, 
		SrtDiffusionType  output,
		double            *forward_adj,
		double            *vol_adj,
		double            *alpha_adj,
        double            *rho_adj);

Err		srt_f_optsabr_mr_vol(
		double            forward,
		double            strike, 
		double            mat, 
		double            volinput,  
		double            alpha,
		double			  beta,
		double			  rho,
		double			  lambda,
		double			  *numericalparams,
		int				  nparams,
		SrtDiffusionType  input,
		SrtDiffusionType  output,
		double            *vol);

Err		srt_f_cutvol(
		double            forward,
		double            strike, 
		double            mat, 
		double            volinput,  
		double            alpha,
		double			  betaordelay,
		double			  rho,
		double			  lambda,
		double			  *numericalparams,
		int				  nparams,
		SrtDiffusionType  input,
		SrtDiffusionType  output,
		int				  nstd,
		double			  lowbetacut,
		double			  highbetacut,
		Err				  (*bsvol)(	double F, double K, double T, double sigma, double alpha, double betaordelay, double rho, double lambda, double *numericalparams, int nparams, SrtDiffusionType  typeinput, SrtDiffusionType  typeoutput, double *voloutput),
		double            *vol);


Err srt_f_optblkvolATMtobetavolStochVol(
		double            forward,							
		double            strike, 
		double            bsvol, 
		double            mat, 
		double            alpha,
		double            beta,
		double            rho,
		double           *betavol);

/* =============================================================================
  FUNCTION     : srt_f_optimpstr(...)        
 ============================================================================ */

double 		srt_f_optimpstr(		
		double 				premium, 
		double 				fwd_price, 
		double 				vol,
		double 				mat, 
		double 				disc, 
		SrtCallPutType		call_put,
		SrtDiffusionType	lognormal_normal
		);


/* =============================================================================
  FUNCTION     : srt_f_optblkd12(...)       
 ============================================================================ */

void 	srt_f_optblkd12(			
		double 			fwd_price,			
		double 			strike,		
		double 			vol,	
		double 			mat, 	
		double 			*d1,
		double 			*d2
		);


/* =============================================================================
  FUNCTION     : srt_f_optblknrm(...)
 ============================================================================ */

double 		srt_f_optblknrm(
		double 			fwd_price, 
		double 			strike,
		double 			vol, 
		double 			mat, 
		double 			df,
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		);

double 		srt_f_optblknrm_accurate(
		double 			fwd_price, 
		double 			strike,
		double 			vol, 
		double 			mat, 
		double 			df,
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		);


/* A function that checks the Normal Black-Scholes inputs */
Err srt_f_optchkblknrm(
		double 			vol,
		double 			mat);

/* A Higher level interface function to call Normal Black Scholes with checks*/
Err OptBlkNrm(	
		double 				fwd_price,
		double 				strike,
		double 				vol,
		double 				mat, 
		double 	  			disc, 
		char               *call_put_str, 
		char               *greek_str,
		double             *result
		);
		
/* =============================================================================
  FUNCTION     : srt_f_optameopt(...)                                     
 ============================================================================ */

double 	srt_f_optameopt(
		double 			spot, 
		double			forward, 
		double 			strike, 
		double 			disc_fact, 
		double 			vol, 
		double 			opt_mat, 
		SrtCallPutType		call_put, 
		int 			step_num,
		SrtGreekType		greek);


/* =============================================================================
  FUNCTION     : srt_f_optamevol(...)                                     
 ============================================================================ */

Err 		srt_f_optamevol(
		double 			premium, 
		double			spot,
		double			forward, 
		double 			strike, 
		double 			disc_fact, 
		double 			opt_mat, 
		SrtCallPutType		call_put, 
		int 			step_num,
		double          *implied_vol
		);



/* ========================================================================== 

				PART 2
			   BARRIER OPTIONS

  ========================================================================== */


/* =============================================================================
  FUNCTION     : srt_f_optexting(...)                                      
 ============================================================================ */
 
double 		srt_f_optexting(
		double 			fwd_price, 
		double 			spot, 
		double 			strike, 
		double 			barrier, 
		double 			vol, 
		double 			mat, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType 	down_up, 
		SrtGreekType		greek
		);


/* =============================================================================
  FUNCTION     : srt_f_optprtext(...)                                 
 ============================================================================ */

double 	srt_f_optprtext(
		double 			fwd2, 
		double 			fwd1, 
		double 			spot, 
		double 			strike,
		double 			barrier, 
		double 			vol1, 
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType	down_up,
		SrtGreekType		greek);


/* =============================================================================
  FUNCTION     : srt_f_optfwdext(...)                                 
 ============================================================================ */

double 	srt_f_optfwdext(
		double 			fwd1, 
		double 			fwd2, 
		double 			strike,
		double 			barrier, 
		double 			vol1, 
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType	down_up,
		SrtGreekType		greek); 

/* =============================================================================
  FUNCTION     : srt_f_optlgtabl(...)                                       
 ============================================================================ */

double 		srt_f_optlgtabl(
		double 			fwd_price, 
		double 			spot, 
		double 			strike, 	
		double 			barrier, 
		double 			vol, 
		double 			mat, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType	down_up, 
		SrtGreekType		greek);


/* =============================================================================
  FUNCTION     : srt_f_optprtlgt(...)                                  
 ============================================================================ */

double 	srt_f_optprtlgt(
		double 			fwd1, 
		double 			fwd2, 
		double 			spot, 
		double 			strike,
		double 			barrier, 
		double 			vol1,
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType	down_up,
		SrtGreekType		greek);


/* =============================================================================
  FUNCTION     : srt_f_optfwdlgt(...)                                  
 ============================================================================ */

double 	srt_f_optfwdlgt(
		double 			fwd1, 
		double 			fwd2, 
		double 			strike,
		double 			barrier, 
		double 			vol1,
		double			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType	down_up,
		SrtGreekType		greek
		);


/* ========================================================================== 
   srt_f_optdblbar(...)
   ========================================================================== */

double 	srt_f_optdblbar(
		double 			spot,          
		double 			forward,          
		double 			strike, 
		double 			b_do, 
		double 			b_up, 
		double 			vol, 
		double 			opt_mat,
		double 			disc, 
		SrtCallPutType		call_put, 
		int			nb_term,
		SrtGreekType		greek
		);


/* ========================================================================== 
   srt_f_optprtdbl(...)
   ========================================================================== */

double srt_f_optprtdbl(
		double 			fwd1,          
		double 			fwd2,          
		double 			spot,          
		double 			strike, 
		double 			b_do, 
		double 			b_up, 
		double 			vol1,  
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtCallPutType		call_put, 
		int 			nb_term,
		SrtGreekType		greek
		);


/* ========================================================================== 
   srt_f_optprtdbl(...)
   ========================================================================== */

double srt_f_optfwddbl(
		double 	fwd1,         
		double 	fwd2,         
		double 	spot,         
		double 	strike, 
		double 	b_do, 
		double 	b_up, 
		double 	vol1, 
		double 	vol2, 
		double 	mat1,  
		double 	mat2,    
		double 	disc, 
		SrtCallPutType 	call_put, 
		int 	nb_term,
		SrtGreekType 	greek
		);                   


/* ========================================================================== 
   lightable_swap_fct(...)
   ========================================================================== */

double lightable_swap_fct	(
		int			num_path,
		double			y20,
		double			sy2,
		double			y30,
		double			sy3,
		double			rho,
		int			compd,
		int 			nfp,
		double			coupon,
		double			spread,
		double			t,
		double			df,
		int			condition	
		);

/* ========================================================================== 

				PART 3
			 TWO UNDERLYING OPTIONS

  ========================================================================== */


/* =============================================================================
  FUNCTION     : srt_f_optsprnum(...)              	    
 ============================================================================ */

double 	srt_f_optsprnum(
		double 			fwdy, 
		double 			fwdx,
		double 			sigy,
		double 			sigx,
		double 			rho,
		double 			disc,
		double 			a,
		double 			b,
		double 			strike,
		double 			mat,
		SrtCallPutType		call_put, 
		SrtGreekType		greek
		);

/* =============================================================================
  FUNCTION     : srt_f_optworstbest(...)              	    
 ============================================================================ */
double srt_f_optworstbest	(double				FirstForward,
							 double				SecondForward,
							double				FirstStrike,
							double				SecondStrike,
							double				Gearing,
							double				FirstBSVol,
							double				SecondBSVol,
							double				Rho,
							double				Df,
							double				Maturity,
							SrtCallPutType 		Call_PutOne,
							SrtCallPutType 		Call_PutTwo,
							SrtBestWorstType	Best_Worst,
							SrtGreekType 		Greek);
							
							
/* =============================================================================
  FUNCTION     : srt_f_optextspr(...)              	    
 ============================================================================ */

double 	srt_f_optextspr(
		double 			fwdy, 
		double 			fwdx,
		double 			sigy,
		double 			sigx,
		double 			rho,
		double 			disc,
		double 			a,
		double          ap,
		double 			b,
		double          bp,
		double 			k,
		double 			t,
		SrtMinmaxType	min_max,
		SrtCallPutType	call_put,
		SrtGreekType	greek
		);

/* =============================================================================
  FUNCTION     : srt_f_opttfteur(...)
 ============================================================================ */
 
double 	srt_f_opttfteur(	
		double 			fwd1, 
		double 			fwd2, 
		double 			strike,
		double 			a, 
		double 			b, 
		double 			vol1, 
		double 			vol2, 
		double 			rho12,
		double 			mat,
		double 			disc_factor,
		SrtCallPutType		call_put, 
		int 			step,
		SrtGreekType		greek
		);
             
/* =============================================================================
  FUNCTION     : srt_f_opttftame(...)       
 ============================================================================ */
 
double 	srt_f_opttftame(	
		double 			fwd1, 
		double 			fwd2, 
		double 			spot1,
		double 			spot2, 
		double 			strike,
		double 			a, 
		double 			b, 
		double 			vol1, 
		double 			vol2, 
		double 			rho12,
		double 			mat,
		double 			disc_factor,
		SrtCallPutType		call_put, 
		int 			step,
		SrtGreekType		greek
		);


/* =============================================================================
  FUNCTION     : srt_f_nono_bar(...)
 ============================================================================ */

double srt_f_nono_bar (double		spot1,
			double		spot2,
			double		fwd1,
			double		fwd2,
			double		strike,
			double		mat,
			double		vol1,
			double		vol2,
			double		rho,
			double		barrier[4],
			double		disc_fac,
			long		n_steps,
			SrtCallPutType	call_put,
			SrtGreekType	greek ) ;
    
/* =============================================================================
  FUNCTION     : srt_f_opteurscdout(...)
 ============================================================================ */

double 	srt_f_opteurscdout(
		double 			fwdx, 
		double 			fwdy, 
		double 			strike,
		double 			barrier, 
		double 			sigx, 
		double 			sigy, 
		double 			rho, 
		double 			matx,
		double			maty,
   		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType 		down_up,
		SrtGreekType		greek
		);

/* =============================================================================
  FUNCTION     : srt_f_optscdout(...)
 ============================================================================ */

double 	srt_f_optscdout(
		double 			fwdx, 
		double 			fwdy, 
		double 			spoty, 
		double 			strike,
		double 			barrier, 
		double 			sigx, 
		double 			sigy, 
		double 			rho, 
		double 			mat1,
		double			mat2,
   		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType 		down_up,
		SrtGreekType		greek
		);

/* =============================================================================
  FUNCTION     : srt_f_optscdinn(...)                                         
 ============================================================================ */

double 	srt_f_optscdinn(
		double 			fwdx, 
		double 			fwdy, 
		double 			spoty, 
		double 			strike,
		double 			barrier, 
		double 			sigx, 
		double 			sigy, 
		double 			rho, 
		double 			mat1,
		double 			mat2,
   		double 			disc, 
		SrtCallPutType		call_put, 
		SrtBarrierType 		down_up,
		SrtGreekType		greek
		);

/* =============================================================================
  FUNCTION     : srt_f_optdbscud(...)
 ============================================================================ */

double 	srt_f_optdbscud	(
		double 			fwdx, 
		double 			fwdy, 
		double 			spotx, 
		double 			spoty, 
		double 			strike,
		double 			barrierx, 
		double 			barriery, 
		double 			sigx, 
		double 			sigy, 
		double 			rho, 
		double 			mat,
   		double 			disc, 
		SrtCallPutType		call_put, 
		int 			scud_type1,
		int 			scud_type2
		);

/* =============================================================================
  FUNCTION     : srt_f_opttriscud(...)
****************************************************************************** 
**   		ATTENTION: X AND Y ARE SWAPPED FROM ADDIN IN 2020
**			   JUST A CASE OF NOTATION             
 ============================================================================ */

double 	srt_f_opttriscud	(
		double 			fwdy, 
		double 			fwdx, 
		double 			spoty, 
		double 			spotx, 
		double 			strikey,   
		double 			barriery_up, 
		double 			barriery_dwn, 
		double 			barrierx,   
		double 			sigy,       
		double 			sigx,       
		double 			rho, 
		double 			mat,
   		double 	 		disc, 
		SrtCallPutType		call_put,   /* ca_pu in y */     
		int 			scud_type, /* up_do in x */     
		int 			nterms,
		SrtGreekType		greek
		);


/* =============================================================================
* FUNCTION     	: srt_f_optscdspr(14)                                         
 ============================================================================ */


double 	srt_f_optscdspr	(
		double 			fwdx, 
		double 			fwdy, 
		double 			fwdz, 
		double 			spoty, 
		double 			spotz, 
		double 			strike,
		double 			barrier, 
		double 			sigx, 
		double 			sigy, 
		double 			sigz, 
		double 			rhoxy, 
		double 			rhoyz, 
		double 			rhoxz, 
		double 			matx,
		double 			matyz,
   		double 			disc, 
		SrtCallPutType		call_put,
		int 			scud_type, 	/* DOWN or UP*/
		int 			ext_lit, 	/* EXTING or LIGHT*/
		SrtGreekType		greek
		); 


/* =============================================================================
  FUNCTION     	: srt_f_optmltspd(...)                                         
 ============================================================================ */
double  srt_f_optmltspd (
                        double  fwdx,
                        double  fwdy,
                        double  fwdz,
                        double  sigx,
                        double  sigy,
                        double  sigz,
                        double  rhoxy,
						double 	rhoxz,
						double  rhoyz,
						double  a_x,
                        double  b_y,
						double 	c_z,
                        double  strike,
                        double  mat,
                        double  disc,
                        SrtCallPutType         call_put,
						long paths,
                        SrtGreekType    greek
                        );


/* =============================================================================
  FUNCTION     	: srt_f_EuroBasketSmile(...)                                         
 ============================================================================ */

Err	  srt_f_EuroBasketSmile (
						int		nFwds,	
						double  *Fwds,
						int		nWeights,
						double	*Weights,
						double	strike,
						int		nVols,
						double	*Vols,
						int		nBeta,
						double	*Beta,
						int		nAlpha,
						double	*Alpha,
						int		nRho,
						double	*Rho,
                        double  **Correlation,
                        double  mat,
                        double  disc,
                        SrtCallPutType         call_put,
						long paths,
						long points,
                        SrtGreekType    greek,
                        double  *premium,
						double  **ppdSpreadVol);

/* =============================================================================
  FUNCTION PAYOFF 	: MAX(0,w1*S1 + w2*S2 + ... + wN*SN )  with Student copulas                                        
 ============================================================================ */

Err srt_f_EuroBasketSmileCpl (
								int		nFwds,	
								double  *Fwds,
								int		nWeights,
								double	*Weights,
								double	strike,
								int		nVols,
								double	*Vols,
								int		nBeta,
								double	*Beta,
								int		nAlpha,
								double	*Alpha,
								int	   nRho,
								double *Rho,
								double **Correlation,
								double mat,
								double disc,
								SrtCallPutType         call_put,
								long nPaths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								double **ppdSpreadVol,
								int n_conv,
								SrtDiffusionType  TypeInput,
								SrtMCSamType  MCType 
								);

/* Basket Copula but with Mixed of shifted log as underlying dynamics */
/* ****************************************************************** */
Err srt_f_Basket_MixedSL_Copula (
								double			dMaturity,
												
								int				iNbFwd,
								double			*dShiftedFwdDown,
								double			*dShiftedFwdUp,
								double			*dProba,
								double			*dShift,
								double			*dVolDown,
								double			*dVolUp,
												
								double			**dCorrelation,
												
								double			*dWeights,
								int				iNbStrike,
								double			*dStrike,
																
								SrtCallPutType  eCallPut,
								long			lNbPaths,
								long			lNbPoints,
								double			dNbStd,
								int				iStudDegree,        // this is the only difference: input the degree of the 
																	// Student t Copula. 0 Corresponds to the gaussian copula								
								SrtMCSamType	MCType,
								int				iNConv,
								double			*dPremium);



/* =====================================================================================================
  FUNCTION PAYOFF	:  MAX(w1+S1, w2+S2, ..., wN+SN) + MarginFix if  MAX(w1+S1, w2+S2, ..., wN+SN) > strike
				               CouponFix						 if  MAX(w1+S1, w2+S2, ..., wN+SN) < strike                                                                                    
				    with Student copulas                                        
 ====================================================================================================== */

Err srt_f_BestOfDigitalSmileCpl (
								int		nFwds,	
								double  *Fwds,
								int		nWeights,
								double	*Weights,
								double	strike,
								int		nVols,
								double	*Vols,
								int		nBeta,
								double	*Beta,
								int		nAlpha,
								double	*Alpha,
								int	   nRho,
								double *Rho,
								double **Correlation,
								double mat,
								double disc,
								SrtCallPutType         call_put,
								long nPaths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								double **ppdSpreadVol,
								double dCoupon,
								double dMargin,
								int n_conv,
								SrtDiffusionType  TypeInput,
								SrtMCSamType  MCType 
								);

/* =====================================================================================================
  FUNCTION PAYOFF	:  MAX (0,MAX(w1+S1, w2+S2, ..., wN+SN)-strike)                                                                                   
				                  with Student copulas                                        
 ====================================================================================================== */

Err srt_f_BestOfSmileCpl (
								int		nFwds,	
								double  *Fwds,
								int		nWeights,
								double	*Weights,
								double	strike,
								int		nVols,
								double	*Vols,
								int		nBeta,
								double	*Beta,
								int		nAlpha,
								double	*Alpha,
								int	   nRho,
								double *Rho,
								double **Correlation,
								double mat,
								double disc,
								SrtCallPutType         call_put,
								long nPaths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								double **ppdSpreadVol,
								int n_conv,
								SrtDiffusionType  TypeInput,
								SrtMCSamType  MCType 
								);

/* =====================================================================================================
  FUNCTION PAYOFF	:	S1^w1*S2^w2* ...*Sn^wn + MarginFix		if  S1^w1*S2^w2* ...*Sn^wn > strike
				        CouponFix								if  S1^w1*S2^w2* ...*Sn^wn < strike                                                                                    
						with Student copulas                                        
 ====================================================================================================== */

Err srt_f_WeightedPdtSmileCpl (
								int		nFwds,	
								double  *Fwds,
								int		nWeights,
								double	*Weights,
								double	strike,
								int		nVols,
								double	*Vols,
								int		nBeta,
								double	*Beta,
								int		nAlpha,
								double	*Alpha,
								int	   nRho,
								double *Rho,
								double **Correlation,
								double mat,
								double disc,
								SrtCallPutType         call_put,
								long nPaths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								double **ppdSpreadVol,
								double dCoupon,
								double dMargin,
								int n_conv,
								SrtDiffusionType  TypeInput,
								SrtMCSamType  MCType 
								);


/* ====================================================================================
  FUNCTION     	: srt_f_CompundBasketCpl(...)  Compound option on a basket with copulas                                       
 =================================================================================== */

Err srt_f_CompoundBasketCpl (
								int		nFwds,	
								double  *Fwds,
								int		nWeights,
								double	*Weights,

								int		nStrikes,
								double *Strikes,

								int		nVols,
								double	*Vols,
								int		nBeta,
								double	*Beta,
								int		nAlpha,
								double	*Alpha,
								int	   nRho,
								double *Rho,
								double **Correlation,
								double mat,
								double disc,
								SrtCallPutType         call_put,
								long paths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								int n_conv,
								SrtDiffusionType  TypeInput,
								SrtDiffusionType  log_norm 
								);

/* ====================================================================================
  FUNCTION     	: srt_f_CompundBasketCpl(...)  Callable option on a basket with copulas                                       
 =================================================================================== */

Err srt_f_CallableBasketCpl (
								 double **InFwds,
								double	*InWeights,
								double *InStrikes,
								double **InVols,
								double	*InBeta,
								double	*InAlpha,
								double *InRho,
								double	*InLvl,
								double **InCorr,
								double mat,
								double disc,

								SrtCallPutType         call_put,

								long paths,
								long points,
								int StudDegree,        // this is the only difference: input the degree of the 
								                       // Student t Copula. 0 Corresponds to the gaussian copula
								SrtGreekType greek,
								double *premium,
								int n_conv,
								const long nFwdsRows,
								const long nFwdsCols,
								SrtDiffusionType  TypeInput,
								SrtDiffusionType  log_norm 
								);


/* =============================================================================
  FUNCTION     : srt_f_optmultpl(...)                                        
 ============================================================================ */

double 		srt_f_optmultpl(
		double 			fwd1, 
		double 			strike1, 
		double 			fwd2, 
		double 			strike2, 
		double 			vol1, 
		double 			vol2, 
		double 			corr, 
		double 			mat, 
		double 			disc,
		SrtCallPutType		call_put1, 
		SrtCallPutType		call_put2, 
		SrtGreekType		greek
		);



/* ========================================================================== 

				PART 4
			PATH DEPENDENT OPTIONS

  ========================================================================== */


/* =============================================================================
  FUNCTION     : srt_f_optikosmo(...)
 ============================================================================ */

double 	srt_f_optikosmo(
		double 			fwd_price1, 
		double 			fwd_price2, 
		double 			spot_price,
		double 			barrier,
		double 			hist_extreme, 
		double 			vol1, 
		double 			vol2, 
		double 			mat1, 
		double 			mat2,
		double 			disc,
		SrtCallPutType		call_put, 
		int 			num_rungs,
		double 			rungs[MAX_RUNG], 
		double 			strikes[MAX_RUNG],
		int 			step_num,
		SrtGreekType	greek
		);

/* =============================================================================
* FUNCTION     	: srt_f_optratcht()                                       
 ============================================================================ */

double 	srt_f_optratcht(
		double 			fwd1, 
		double 			fwd2, 
		double 			spot,
		double 			hist_ext, 
		double 			vol1,
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc,
		SrtCallPutType		call_put, 
		SrtGreekType		greek,
		double 			rungs[MAX_RUNG], 
		double 			strikes[MAX_RUNG]
		);


/* ============================================================================
   srt_f_optonetou(...)
   ==========================================================================*/

double	srt_f_optonetou	(	
		double 			fwd_price, 
		double 			spot, 
	   	double 			strike, 
	  	double 			barrier, 
		double 			vol, 
		double 			mat, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtGreekType		greek 
		);

/* ============================================================================
   srt_f_optrebate(...)
   ==========================================================================*/

double 	srt_f_optrebate	(
		double 			fwd_price, 
		double 			spot, 
		double 			barrier, 
		double 			rebate, 
		double 			vol, 
		double 			mat, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtGreekType		greek 
		);

/* ============================================================================
   srt_f_optladder(...)
   ==========================================================================*/

double  srt_f_optladder	(	
		double 			fwd_price, 
		double 			spot, 
		double 			strike, 
		double 			hist_extr, 
		double 			rungs[MAX_RUNG], 
		double 			vol, 
		double 			mat, 
		double 			disc, 
		SrtCallPutType		call_put, 
		SrtLadderType		grad_best, 
		SrtGreekType		greek);

/* ============================================================================
   srt_f_optlokbak(...)
   ==========================================================================*/

double 	srt_f_optlokbak (	
		double 			fwd_price, 
		double 			spot, 
		double 			hist_extr,
		double 			vol, 
		double 			mat, 
		double 			disc,
		SrtCallPutType		call_put, 
		SrtGreekType		greek 
		);

/* ============================================================================
   srt_f_optpaloba(...)
   ==========================================================================*/
                        
double 	srt_f_optpaloba	( 	
		double 			fwd_price1, 
		double 			fwd_price2,
		double 			spot,
		double 			hist_extr,
		double 			vol1, 
		double 			vol2, 
		double 			mat1, 
		double 			mat2,
		double 			disc,
		SrtCallPutType		call_put, 
		SrtGreekType		greek 
		);
                        

/* ============================================================================
   srt_f_optextrema(...)
   ==========================================================================*/

double 	srt_f_optextrema (	
		double 			fwd_price, 
		double 			spot, 
		double 			strike, 
		double 			hist_extr,
		double 			vol, 
		double 			mat, 
		double 			disc,
		SrtCallPutType		call_put, 
		SrtGreekType		greek 
		);


/* =========================================================================
  FUNCTION     : srt_f_optlifopt(...)               
   ========================================================================= */

double srt_f_optlifopt(	
		double 			fwd, 
		double	 		spot, 
		double 			barrier, 
		double 			vol, 
		double	 		mat, 
		double 			disc, 
		SrtBarrierType 	down_up,
		SrtGreekType		greek);

/* ============================================================================
   FUNCTION: srt_f_optfwdlife
   ==========================================================================*/

double srt_f_optfwdlife(	
		double 			fwd1, 
		double 			fwd2, 
		double 			spot, 
		double 			barrier, 
		double 			vol1, 
		double 			vol2, 
		double 			mat1, 
		double 			mat2, 
		double 			disc1,
		double			disc2,
		SrtBarrierType 	down_up,
		SrtGreekType		greek) ;


/* ============================================================================
   FUNCTION: srt_f_optlifout
   ==========================================================================*/

double srt_f_optlifout(	
		double 			fwd, 
		double			spot,
		double			b_life,
		double			b_out,
		double			vol,
		double			mat,
		double			disc,
		int			nterms,
		SrtGreekType		greek);

Err srt_f_asian( int nforwards, double* forwards_date, double* forwards, 
				  int nfixing, double* fixing_date,
				  int nvol, double* vol_date, double* volat,
				  double spot, double maturity, double strike, 
				  int npast_fix, double avg_cur, double disc, 
			      SrtCallPutType 	call_put,
				  SrtGreekType greek,
				  double *answer);

/* ===========================================================================

				PART 5
			 DIGITAL OPTIONS

  ============================================================================*/


/* =============================================================================
  FUNCTION     : srt_f_opteurdig(...)                                      
 ============================================================================ */

double srt_f_opteurdig(	
		double 			fwd, 
		double			barrier,
		double 			vol, 
		double 			mat,
		double 			disc, 
		SrtBarrierType		below_above,
		SrtGreekType		greek
		);

/* =============================================================================
  FUNCTION     : srt_f_optamedig(...)                                      
 ============================================================================ */

double 	srt_f_optamedig(	
		double 			fwd,
		double 			spot,
		double			barrier,
		double 			vol,
		double 			mat, 
		double 			disc, 
		SrtBarrierType		below_above,
		SrtMinmaxType 		min_max,
		SrtGreekType		greek
			);

/* ========================================================================== 
   srt_f_optfwdamdig(...)
   ========================================================================== */

double 	srt_f_optfwdamdig(	
		double 			fwd_price1,
		double 			fwd_price2,
		double 			strike,
		double 			vol1,
		double 			vol2,
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtBarrierType		below_above,
		SrtMinmaxType 		min_max,
		SrtGreekType		greek
		);
                                   
/* ========================================================================== 
   srt_f_opttpamdi(...)                                      
   ========================================================================== */

double 	srt_f_opttpamdi  	(
		double 	fwd1, 
		double 	fwd2, 
		double 	barrier, 
		double 	vol1, 
		double 	vol2, 
		double 	mat1, 
		double 	mat2, 
		double 	disc1, 
		double 	disc2, 
		SrtBarrierType 	below_above,
		SrtGreekType 	greek
		); 
 

/* ========================================================================== 
   srt_f_optfirsec(...)                                      
   ========================================================================== */

double 	srt_f_optfirsec(	
		double 			fwd, 
		double 			spot,  
		double			fir_bar, 
		double			sec_bar, 
		double 			vol, 
		double 			mat,
		double 			disc,    
		int			nterms,
		SrtGreekType		greek
		);


/* ========================================================================== 
   srt_f_opteursla(...)
   ========================================================================== */

double 		srt_f_opteursla(	
		double 			fwd_price1,
		double 			fwd_price2,
		double 			barrier,
		double 			vol1,
		double 			vol2,  
		double 			mat1, 
		double 			mat2, 
		double 			disc,
		SrtBarrierType		below_above,
		SrtGreekType		greek 
		);
                                 
/* ========================================================================== 
   srt_f_optamesla(...)
   ========================================================================== */

double 		srt_f_optamesla(	
		double 			fwd_price1,
		double 			fwd_price2,
		double 			strike,
		double 			vol1,
		double 			vol2,    
		double 			mat1, 
		double 			mat2, 
		double 			disc, 
		SrtBarrierType		below_above,
		SrtGreekType		greek 
		);
                                   

/* ========================================================================== 
   srt_f_optslalom(...)
   ========================================================================== */

double srt_f_optslalom(
		double 			forward, 
		double 			spot, 
		double 			barrier, 
		double 			vol, 
		double 			disc_fact, 
		int 			today, 
		int 			nb_dates, 
		int 			date[MAX_NUM_DATES], 
		int 			up_do[MAX_NUM_DATES], 
		int 			step_num,
		SrtGreekType	greek
		);

/* ========================================================================== 
   srt_f_optrainbow(...)
   ========================================================================== */

double 	srt_f_optrainbow(	double	T,		/* maturity			*/
							double	df,		/* discount factor	*/
							double	K,
							double	nX,
							double	X0,
							double	sigX,
							double	nY,
							double	Y0,
							double	sigY,
							double	rho);

/* ========================================================================== 
   srt_f_optmovbar(...)
   ========================================================================== */

double srt_f_optmovbar(
		double spot,
		double fwd1,
		double fwd2,
		double bar1,
		double bar2,
		double sigma1,
		double sigma2,
		double mat1,
		double mat2, 
		double disc_fact,
		SrtBarrierType 	below_above,
		SrtGreekType	greek
		);


/* ========================================================================== 

				PART 6
		  	  COMPOUND OPTIONS

  ========================================================================== */

/* =============================================================================
  FUNCTION     : srt_f_opt_extcombar(...)                                  
 ============================================================================ */

double 	srt_f_opt_extcombar(     
		double 			fwd1, 
		double 			fwd2, 
		double 			spot, 
		double 			strike1,
		double 			strike2,
		double 			barrier1,
		double 			barrier2, 
		double 			sig1,      
		double 			sig2,     
		double 			mat1, 
		double 			mat2, 
		double 			disc1, 
		double 			disc2,     
		SrtCallPutType 	call_put1,
		SrtCallPutType 	call_put2, 
		SrtBarrierType	scud_type1, 
		SrtBarrierType  scud_type2,
		SrtGreekType 	return_type );


/* ========================================================================== 

				PART 7
		  DISTRIBUTIONS and PROBABILITIES

  ========================================================================== */




/* =============================================================================
  FUNCTION     : srt_f_optnrmdis(...) 	    
 ============================================================================ */

double srt_f_optnrmdis(
		double			value
		);


/* =============================================================================
  FUNCTION     : srt_f_optbivdis(...) 	    
 ============================================================================ */

double srt_f_optbivdis(
		double			variable_x,
		double			variable_y,
		double			rho
		);

/* ========================================================================== 
   trivar(...)
   ========================================================================== */

double 	trivar( double 			a1, 
		double 			a2, 
		double 			a3, 
		double 			r12, 
		double 			r23, 
		double 			r13);


/* ========================================================================== 
   dbarrier_prob(...)
   This function computes the following probabilities for a drifted
   normal Brownian motion:
				Prob[ X(T) > k; min X(t) > b_do ; max X(t) < b_up ]
				Prob[ X(T) < k; min X(t) > b_do ; max X(t) < b_up ]
   ========================================================================== */


double dbarrier_prob( 	
			double 		k, 
			double 		b_up,
			double		b_do,
			double 		mu, 
			double 		vol, 
			double 		mat,
			SrtCallPutType	call_put, 
			int 		nb_terms);


/* ========================================================================== 
   stay_inside_fct(...)
   ========================================================================== */

double stay_inside_fct( 
			double 	spot, 
			double 	fwd, 
			double 	up, 
			double 	low, 
			double 	mat, 
			double 	vol, 
			int 	n_terms, 
			int 	measure_flag);

/* ========================================================================== 
   stay_inside_fwd_fct(...)
   ========================================================================== */

double stay_inside_fwd_fct( 
			double 	spot, 
			double 	fwd1, 
			double 	fwd2, 
			double 	up, 
			double 	low, 
			double 	t1, 
			double 	t2, 
			double 	vol0, 
			double 	vol, 
			int 	n, 
			int 	ff);

/* ========================================================================== 
   fct_gold_stay_inside(...)
   ========================================================================== */

double fct_gold_stay_inside(	
			double 	x, 
			va_list argptr);

/* ========================================================================== 
   fct_gold_stay_inside_fwd(...)
   ========================================================================== */

double fct_gold_stay_inside_fwd(	
			double 	x, 
			va_list argptr);

/* ========================================================================== 
   hit_one_barrier_fct(...)
   ========================================================================== */

double hit_one_barrier_fct(
			double 	spot, 
			double 	fwd, 
			double 	b_yes, 
			double 	b_no, 
			double 	mat, 
			double 	vol, 
			int 	n, 
			int 	ff);


/* ========================================================================== 
   fct_gold_hit_one_barrier(...)
   ========================================================================== */

double fct_gold_hit_one_barrier(
			double 	x, 
			va_list argptr);


/* ========================================================================== 
   dens_four_n(...)
   ========================================================================== */

double dens_four_n(	int 	n,
			double 	x,
			double 	a, 
			double 	b,
			double 	coeff1, 
			double 	coeff2);



/* ========================================================================== 
   dens_gauss_n(...)
   ========================================================================== */

double dens_gauss_n(
		double 	var,
		double 	x1_n, 
		double 	x2_n, 
		double 	mu,
		double 	vol,
		double 	mat);

/* ========================================================================== 
   part_dbar_fn(...)
   ========================================================================== */

double part_dbar_fn(
		double 	xx, 
		double 	up, 
		double 	dwn, 
		double 	mu1,
		double 	mu2,   
		double 	vol1,  
		double 	vol2, 
		double 	mat1,    
		double 	mat2,      
		int 	nb_term  
		);             


/* ========================================================================== 
prob_ExpectedDays:  calculates the expected occupation time (in days) for a corridor
   ========================================================================== */
double prob_ExpectedDays( double x1, double x0, double* dBarrier, double dVol, double dStartTime, double dLength, 
							int nStep, double* dParameters, int BBflag );


/* ==========================================================================

   				Utilities

   ========================================================================== */
 
/* =============================================================================
  FUNCTION     : comb(...) 	    
 ============================================================================ */

double 	comb(
		double 	n, 
		double 	n_steps,
		double 	p,
		int 	ii
		);                  


/* =============================================================================
  FUNCTION     : smooth(...) 	    
 ============================================================================ */

double 	smooth (
		double 	y[],
		int 	n, 
		double 	pts
		);

/* ========================================================================== */

/* ============================================================================
   trivar(...)
   ==========================================================================*/

double trivar( 	
		double 			a1, 
		double 			a2, 
		double 			a3, 
		double 			r12, 
		double 			r23, 
		double 			r13 );

/* ============================================================================
   Asian Price for IRD, Moment Matching (order 2 and 3)
   ==========================================================================*/

Err srt_f_AsianCapPrice(double				*dStartDates, 
						double				*dCmsForwards,
						double				*dCmsVols,
						double				*dCorrels,
						int					iCms,
						double				dStrike, 
						SrtCallPutType 		SrtCallPut,
						SrtGreekType		SrtGreek,
						SrtDiffusionType	SrtVolType,
						int					iMethod,
						double				*dAnswer);

double srt_f_GammaInvOption (double 			dStrike,
							 double 			dM1,
							 double				dM2,
							 SrtCallPutType 	call_put, 
							 SrtGreekType		greek);

double srt_f_GammaInvOpt(double 			strike,
						 double 			a,
						 double				b,
						 SrtCallPutType 	call_put, 
						 SrtGreekType		greek);

double srt_f_GammaOption(double 			dStrike,
						 double 			dM1,
						 double				dM2,
						 SrtCallPutType 	call_put, 
						 SrtGreekType		greek);

double srt_f_GammaOpt(double 			strike,
					  double 			a,
					  double			b,
					  SrtCallPutType 	call_put, 
					  SrtGreekType		greek);

/* ======================================================================== */

/* ==========================================================================

   				New SABR functions

   ========================================================================== */

#include "opsabr.h"
#include "opbasket.h"

#endif
