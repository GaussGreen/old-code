/**************************************************************************************/
/*********************** Convolver ****************************************************/

/* Uses convolution to evaluate deals within a LGM framework **************************/
/* Outputs:
	*answer = value of the deal,
	xExBdry[0,1,...,nEx-1] are the x values of the exercise boundary at each exercise
								if multiple exercise points, gives the nearest to the money
*/
/* Convolver requires three sets of arguments
	1) the number of exercises and the zeta's at the exercise points
	2a) information about the market (including today, a yield curve, and a calibrated
			LGM term structure
	2b) information about caplet/swaption vols (for deals whose payoffs include options)
	3) information about the deal, and a function to calculate the payoff for all
			the exercise dates 
The arguments of the payoff function are:
	(*payofffunc)(double payoff[][], long nx, double x[], double reduction[] 
					LGMDealType dealtype, void *dealPtr, 
					Date EvalDate, String ycName, LGM_TS tsPtr)
 */

LGMErr Midat2DPayoff(double *payoff, double *swap, double X, double Y,void *dealPtr, 
					 Date tNow, long jEx,String ycName,
					 double *Zeta1,double *Zeta2,double *Zeta12,
					 double *H1,double *H2,double gamma,
					 LGMErr	(*GetVol)(Date, Date, double, SRT_Boolean, double*),	
				     LGMErr	(*GetBeta)(Date, Date, double*));	

LGMErr LGMAutoCal2DTree(
/* info about convolutions */
		long		nEx,			/* number of exercises */
		Date		*tEx,
		double		*zeta11,			/* [0, 1, ..., nEx-1], values of zeta11 at the exercise dates */
		double		*zeta12,			/* [0, 1, ..., nEx-1], values of zeta12 at the exercise dates */
		double		*zeta22,			/* [0, 1, ..., nEx-1], values of zeta22 at the exercise dates */
		double		*G1,				/* [0, 1, ..., nEx-1], values of G1 at the exercise dates */
		double		*G2,				/* [0, 1, ..., nEx-1], values of G2 at the exercise dates */
		double      gamma,
	/* info about today's discount curve and swaption/cap vols */ 
		Date		EvalDate,		/* tNow for deal evaluation */
		String		ycName,			/* yield curve name */
		LGMErr		(*GetVol)(Date, Date, double, SRT_Boolean, double*),	/* swaption/cap vols */
		LGMErr		(*GetBeta)(Date, Date, double*),		/* swaption/cap exponents (beta) */
	/* information about the deal */
		void		*dealPtr,
		LGMErr		(*payofffunc)(),
	/* output */
		double		*answer,		/* value of the deal */
		double		**xExBdry);		/* array of exercise points (in x) */

LGMErr LGMAutoCal2DVega(
		int			nScenarii,		/* number of scenarii different than spot */
		long		nEx,			/* number of exercises */
		Date		*tEx,
		double		**zeta11,		/* for each scenario, values of */
		double		**zeta12,		/* zeta11 at the exercise dates */
		double		**zeta22,			
		double		**G1,				
		double		**G2,
		double		gamma,
	/* info about today's discount curve and swaption/cap vols */ 
		Date		EvalDate,		/* tNow for deal evaluation */
		String		ycName,			/* yield curve name */
		LGMErr		(*GetVol)(Date, Date, double, SRT_Boolean, double*),	/* swaption/cap vols */
		LGMErr		(*GetBeta)(Date, Date, double*),		/* swaption/cap exponents (beta) */
	/* information about the deal */
		void		*dealPtr,
		LGMErr		(*payofffunc)(),
	/* output */
		double		*answer,		/* value of the deal */
		double		**xExBdry);		/* array of exercise points (in x) */

/* ========= END OF FILE =================================================== */
