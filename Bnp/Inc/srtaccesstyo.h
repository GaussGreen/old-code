/*	SrtAccessTyo
	Files from Tokyo for SrtAccess	*/

#ifndef SRTACCESSTYO_H
#define	SRTACCESSTYO_H

#include "srt_h_all.h"
#include "grf_h_all.h"

#include "LGMSVPDE.h"
#include "MCEBOptimisation.h"
#include "Fx3FBetaDLMUtil.h"

char *SrtGrfnFxSabrAdi(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val);

char *SrtGrfnFxSabrAdi2(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		*is_up_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val);

Err FxSabrSmile(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res);

Err FxSabrSLSmile(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res);

Err FxSabrQuadSmile(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  long		mat_date,
   				  double	*strike,
				  int		nb_strike,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	**res);

Err	FxSabrAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  double	barrier,
				  int		is_up,		/* 1 Up, 0: Down */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res);

char *SrtGrfnLGM2Fpde(
				  char		*underlying,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  int		nstept,
				  int		nstepx,
				  int		*nprod,
				  double	**prod_val);

char *SrtGrfnLGM2FTaupde(
				  char		*underlying,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  int		nstept,
				  int		nstepx,
				  int		*nb_prod,
				  double	**prod_val);

char *SrtGrfnLGM2FTaupde2(
				  char		*underlying,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  int		nstept,
				  int		nstepx,
				  int		*nb_prod,
				  double	**prod_val);

char *SrtGrfnLGM2FMC(
				  char		*lgmund,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,	
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  long		num_paths,  
				  int		do_pecs,
				  int		do_jump,
				  double	***prod_val);

// LGM2F MC with lambda term structure
char *SrtGrfnLGM2FMClambda(
				  char		*lgmund,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  long		num_paths,  
				  int		do_pecs,
				  int		do_jump,
				  double	***prod_val);

char *SrtGrfnLGM2FMCEBlambda(
				  char			*lgmund,
	              int			numeventdates, 
				  long			*eventdates,
				  int			*nUsedEventDates,
				  int			*optimise,
				  double		*fwd_iv,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,
				  int			is_end_of_day_fixing,
				  int			is_end_of_day_payment,
				  long			num_paths,  
				  int			do_pecs,
				  int			do_jump,
				  double		***prod_val);
// END of LGM2F MC with lambda term structure

char *SrtGrfnLGM2FMCEB(
				  char			*lgmund,
	              int			numeventdates, 
				  long			*eventdates,				  
				  int			*nUsedEventDates,
				  int			*optimise,
				  double		*fwd_iv,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,
				  int			is_end_of_day_fixing,
				  int			is_end_of_day_payment,
				  long			num_paths,  
				  int			do_pecs,
				  int			do_jump,
				  double		***prod_val);

char *SrtGrfnLGMSVpde(
				  char		*underlying,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,	
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,

				  LGMSVParam Params,				 
				  int		nstept,
				  int		nstepx,
				  int		nstepeps,
				  int		nstepphi,				  

				  int		*nb_prod,
				  double	**prod_val);

char *SrtGrfnLGMSVMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,
					
					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnLGMSVMC_1F(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnLGMSVMC_2F(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					int			*nUsedEventDates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					double		*fwd_iv,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnLGMSVMCCV(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					long		tableauRows,
					long		tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,

					LGMSVParam	Params,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val);

char *SrtGrfnLGMSVFFT(
					/* Model */
					char		*underlying,

					/* Product */
					int			numeventdates, 
					long		*eventdates,
					long		tableauRows,
					long		tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,	
					
					int			is_end_of_day_fixing,
					int			is_end_of_day_payment,
				  
					/* Parameter of grids */
					int			iNbPhi,			/* Number of Phi : Should be a power of two */
					int			iNbft,			/* Number of ft : Should be a power of two */	
					double		iNbSigmaPhiGridLeft,
					double		iNbSigmaPhiGridRight,
					double		iNbSigmaftLeft,
					double		iNbSigmaftRight,
									
					double		iRatioPhi,
					double		iRatioFt,
					int			iPriorityFreqPhi,
					int			iPriorityFreqFt,

					/* outputs */
					int			*nb_prod,
					double		**prod_val);

char *SrtGrfn3DFXTree(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  double	*barrier,
				  int		*bar_col,
				  long		num_stp, 
				  int		*num_prod, 
				  int		discount,
				  double	**prod_val);

/*
char *SrtGrfn3DFXTree_thread(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  long		num_stp, 
				  int		*num_prod, 
				  int		discount,
				  double	**prod_val);
*/

char *SrtGrfn3DFXMc(
				  char		*und3dfx,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  int		is_end_of_day_fixing,
				  int		is_end_of_day_payment,
				  long		num_paths,
				  int		do_pecs,
				  double	***prod_val);

char *SrtGrfn3DFXBetaTree(
						  char		*und3dfx,
						  double	alpha,
						  double	beta,
						  int		numeventdates, 
						  long		*eventdates,
						  long		tableauRows,
						  long		tableauCols,
						  char		***tableauStrings,
						  int		**tableauMask,
						  long		auxWidth,
						  long		*auxLen,
						  double	**aux,
						  int		is_end_of_day_fixing,
						  int		is_end_of_day_payment,
						  double	*barrier, 
						  int		*bar_col,
						  long		num_stp, 
						  int		*num_prod, 
						  int		discount,
						  double	**prod_val);

/*
char *SrtGrfn3DFXBetaTree_thread(
						  char		*und3dfx,
						  double	alpha,
						  double	beta,
						  int		numeventdates, 
						  long		*eventdates,
						  long		tableauRows,
						  long		tableauCols,
						  char		***tableauStrings,
						  int		**tableauMask,
						  long		auxWidth,
						  long		*auxLen,
						  double	**aux,
						  double	*barrier, 
						  int		*bar_col,
						  long		num_stp, 
						  int		*num_prod, 
						  int		discount,
						  double	**prod_val);
*/

char *SrtGrfn3DFXAlphaBetaTree(
							  char		*und3dfx,
							  double	alpha,
							  double	beta,
							  int		numeventdates, 
							  long		*eventdates,
							  long		tableauRows,
							  long		tableauCols,
							  char		***tableauStrings,
							  int		**tableauMask,
							  long		auxWidth,
							  long		*auxLen,
							  double	**aux,
							  int		is_end_of_day_fixing,
							  int		is_end_of_day_payment,
							  double	*barrier,
							  int		*bar_col,
							  long		num_stp, 
							  int		*num_prod, 
							  int		discount,
							  double	**prod_val);

/*
char *SrtGrfn3DFXAlphaBetaTree_thread(
							  char		*und3dfx,
							  double	alpha,
							  double	beta,
							  int		numeventdates, 
							  long		*eventdates,
							  long		tableauRows,
							  long		tableauCols,
							  char		***tableauStrings,
							  int		**tableauMask,
							  long		auxWidth,
							  long		*auxLen,
							  double	**aux,
							  double	*barrier,
							  int		*bar_col,
							  long		num_stp, 
							  int		*num_prod, 
							  int		discount,
							  double	**prod_val);

*/

char *SrtGrfn3DFXAlphaBetaMcTree(
						  char		*und3dfx,
						  double	alpha,
						  double	lambda,
						  double	beta,
						  long		nbPaths,
						  int		numeventdates, 
						  long		*eventdates,
						  long		tableauRows,
						  long		tableauCols,
						  char		***tableauStrings,
						  int		**tableauMask,
						  long		auxWidth,
						  long		*auxLen,
						  double	**aux,
						  double	*barrier, 
						  int		*bar_col,
						  long		num_stp, 
						  int		*num_prod, 
						  int		discount,
						  double	***prod_val);

char *SrtGrfn3DFXBetaMc(
				  char		*und3dfx,
				  double	(*vol_ln_func_Beta) (double	t,
												double	S,
												double	F,
												double	*param),
				  double	*param,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  long		num_paths,
				  double	max_time,
				  int		do_pecs,
				  double	***prod_val);

char *SrtGrfn3DFXQuadTree(
				  char		*und3dfx,
				  double	beta,
				  double	gamma,
				  double	sig0,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  long		num_stp, 
				  int		*num_prod, 
				  int		discount,
				  double	*vols,
				  double	*vols_time,
				  int		nbVols,
				  double	**prod_val);

char *SrtGrfn3DFXQuadMc(
				  char		*und3dfx,
				  double	beta,
				  double	gamma,
				  double	sig0,
				  int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  long		num_paths,  
				  long		num_stp,
				  int		do_pecs,
				  double	***prod_val);

char *SrtMultiGrfn3DFXMc(
				  char		*underlying,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  long		num_paths,  
				  int		do_pecs,
				  double	***prod_val);

char *SrtMultiGrfn3DFXMcEB(
				  char			*underlying,
	              int			numeventdates, 
				  long			*eventdates,
				  int			*optimise,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,				  
				  long			num_paths,  
				  int			do_pecs,
				  double		***prod_val);

char *SrtGrfn5DFXMcTest(
				  char		*und3dfx,
				  double	**correl,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,	
				  long		num_dates,
				  long		num_paths,
				  int		do_pecs,
				  double	***prod_val);

char *SrtGrfn5DFXMc(
				  char		*und3dfx,
				  double	**correl,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		*tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,				  
				  long		num_paths,
				  int		do_pecs,
				  double	***prod_val);

char *SrtGrfn5DFXMcEB(
				  char			*und3dfx,
				  double		**correl,				  
	              int			numeventdates, 
				  long			*eventdates,
				  int			*nUsedEventDates,
				  int			*optimise,
				  double		*fwd_iv,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,				  
				  long			num_paths,
				  int			do_pecs,
				  double		***prod_val);

char *SrtGrfnFxSabrSLAdi(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val);

char *SrtGrfnFxSabrQuadAdi(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              int		numeventdates, 
				  long		*eventdates,
                  long		tableauRows,
				  long		tableauCols,
				  char		***tableauStrings,
				  int		**tableauMask,
				  long		auxWidth,
				  long		*auxLen,
				  double	**aux,
				  double	*barrier,
				  int		*bar_col,
				  int		*is_cont_bar,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  int		*nprod,
				  double	**prod_val);

Err FxSabrSLAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		is_digital,	/* 1: digital payoff, 0, regular option payoff */
				  int		is_american,
				  double	barrier_up,
				  double	barrier_down,
				  double	rebate_up,
				  double	rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks);

Err FxSabrQuadAdiKO(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */
				  double	maturity,
   				  double	strike,
				  int		is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		is_digital,	/* 1: digital payoff, 0, regular option payoff */
				  int		is_american,
				  double	barrier_up,
				  double	barrier_down,
				  double	rebate_up,
				  double	rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks);

Err FxSabrQuadAdi_MultiKO(
				  char		*underlying,
				  double	alpha,
				  double	gamma,
				  double	beta,
				  double	rho,
				  double	lambda,
				  double	floormu,
	              /*	Product data */	
				  int		nb_product,
				  double	*notional,
				  long		*maturity,
				  long		*settlmt,
   				  double	*strike,
				  int		*is_call,	/* 1 Call, 0: Put */
				  int		is_ko,		/* 1 KO, 0: KI */
				  int		*is_cvx,		/* 1 use 1 / Fx, 0 use Fx */
				  int		*is_digital,	/* 1: digital payoff, 0, regular option payoff */				  
				  double	barrier_up,
				  double	barrier_down,
				  double	*rebate_up,
				  double	*rebate_down,
				  int		nstp,
				  int		nstpfx,
				  int		nstpvol,
				  double	*res,
				  double	**greeks);


char *SrtGrfn3DFXMcEB(
				  char			*und3dfx,
	              int			numeventdates, 
				  long			*eventdates,
				  int			*nUsedEventDates,
				  int			*optimise,
				  double		*fwd_iv,
				  MCEBPARAMS	params,
				  long			*resRows,
                  long			tableauRows,
				  long			*tableauCols,
				  char			***tableauStrings,
				  int			**tableauMask,
				  long			auxWidth,
				  long			*auxLen,
				  double		**aux,
				  int			is_end_of_day_fixing,
				  int			is_end_of_day_payment,
				  long			num_paths,  
				  int			do_pecs,
				  double		***prod_val);

#endif
