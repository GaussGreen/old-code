/*******************************************************************************
*			                      Include Files                               
*******************************************************************************/
#include "math.h"
#include		<SPFNCTNS.H>
#include        "utallhdr.h"
#include        <OPFNCTNS.H>

/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_sprrf_tok_anal()                               	  
*                                                                         
* PURPOSE     	: compute price of RRF autoreset in the normal model
*				  analytically with approx price = disc * proba * (fwd + mrg) * cvg
*                                                                         
* DESCRIPTION  	: XX 
*
* CALLS			: norm in gen_math.c                            
*                                                                         
* PARAMETERS  	: today_long 	- today date (in days)
*				: fixing_long	- fixing date (in days)
*				: pay_long		- pay date (in days)
*				: CMSfixing		- forward fixing rate CMS adjusted for pay date
*				: CMSpay		- forward pay rate CMS adjusted for pay date
*				: cvgpay		- coverage of payment
*				: discpay		- discount factor to payment
*				: margin		- margin on payment
*								____________________________________________
*						
*				: K1			- lower bound strike
*				: K2			- upper bound strike
*             	: slopevol		- vol the slope between today and fixing
*             	: fwdvol		- fwd vol of Fpay between fixing and pay  
*				: in_out_flag	- proba to be in (SRT_IN) or out (SRT_OUT)                       	  
*				: 
*                                                                         
* RETURNS      	: *ans			- price
*                                                                         
*******************************************************************************/
     
Err srt_f_sprrf_tok_anal		(long today_long, long fixing_long, long pay_long,
								double CMSfixing, double CMSpay, 
								double cvgpay, double discpay, double margin, 
								double K1, double K2,
								double slopevol, double fwdvol, 
								long in_out_flag, double *ans)
{
	
	double slopevar, fwdvar;
	double CMSslope = CMSpay - CMSfixing;

	slopevar = slopevol * slopevol * (fixing_long - today_long) * YEARS_IN_DAY;
	fwdvar = fwdvol * fwdvol * (pay_long - fixing_long) * YEARS_IN_DAY;

	if (slopevar + fwdvar > 0)
	{

		*ans = 0.0;

		*ans += norm ((K2 - CMSslope) / (sqrt (slopevar + fwdvar)));
		*ans -= norm ((-K1 - CMSslope) / (sqrt (slopevar + fwdvar)));
	}
	else
		*ans = ((CMSslope < K2) && (CMSslope > -K1))? 1.0: 0.0;

	if (in_out_flag == SRT_OUT)
		*ans = 1.0 - (*ans);
	
	*ans *= (CMSpay + margin) * cvgpay * discpay;

	return NULL;
}

/******************************************************************************/

/*******************************************************************************
*                                                                         
* FUNCTION     	: srt_f_sprrf_tok_mc()                               	  
*                                                                         
* PURPOSE     	: compute price of RRF autoreset in the normal model
*				  exactly using NSOBOL sobol monte carlo paths
*                                                                         
* DESCRIPTION  	: XX 
*
* CALLS			: norm in gen_math.c                            
*                                                                         
* PARAMETERS  	: today_long 	- today date (in days)
*				: fixing_long	- fixing date (in days)
*				: pay_long		- pay date (in days)
*				: CMSfixing		- forward fixing rate CMS adjusted for pay date
*				: CMSpay		- forward pay rate CMS adjusted for pay date
*				: cvgpay		- coverage of payment
*				: discpay		- discount factor to payment
*				: margin		- margin on payment
*								____________________________________________
*						
*				: K1			- lower bound strike
*				: K2			- upper bound strike
*				: CMSfvol		- vol of Ffixing between today and fixing
*				: CMSpvol		- vol of Fpay between today and fixing
*             	: rhofp			- correl of Ffixing and Fpay between today and fixing
*             	: fwdvol		- fwd vol of Fpay between fixing and pay  
*				: in_out_flag	- proba to be in (SRT_IN) or out (SRT_OUT)                       	  
*				: 
*                                                                         
* RETURNS      	: *ans			- price
*                                                                         
*******************************************************************************/

Err srt_f_sprrf_tok_mc			(long today_long, long fixing_long, long pay_long,
								double CMSfixing, double CMSpay, 
								double cvgpay, double discpay, double margin,
								double K1, double K2,
								double CMSfvol, double CMSpvol, double rhofp, 
								double fwdvol, 
								long in_out_flag, double *ans)
{
	double slopevol 
		= sqrt (CMSfvol * CMSfvol + CMSpvol * CMSpvol - 2 * rhofp * CMSfvol * CMSpvol);
	double slopevar, fwdvar;
	double CMSslope = CMSpay - CMSfixing;
	double **brownians, correl_brownian;
	double f1, f2;
	double payoff;
	long i;

	Err err;

	slopevar = slopevol * slopevol * (fixing_long - today_long) * YEARS_IN_DAY;
	fwdvar = fwdvol * fwdvol * (pay_long - fixing_long) * YEARS_IN_DAY;

	if (slopevar + fwdvar > 0)
	{

		if (!(brownians = dmatrix (1, 1, 1, 3)))
			return serror ("error in srt_f_sprrf_tok_mc: memory allocation error");

		if (err = sobol_init (1,1, 1, 1, 1, 3))
			return serror ("error in srt_f_sprrf_tok_mc: sobol initialisation error");

		*ans = 0.0;

		for (i=1; i<=NSOBOL; i++)
		{
			if (err = sobol_matrix (brownians, 1, 1, 1, 3))
				return serror ("error in srt_f_sprrf_tok_mc: sobol error");
		
			f1 = CMSfixing 
				+ CMSfvol * sqrt (YEARS_IN_DAY * (fixing_long - today_long))
				* brownians[1][1];

			correl_brownian = rhofp * brownians[1][1] + sqrt (1 - rhofp * rhofp) * brownians[1][2];
	
			f2 = CMSpay
				+ CMSpvol * sqrt (YEARS_IN_DAY * (fixing_long - today_long))
				* correl_brownian;
		
			f2 += fwdvol * sqrt (YEARS_IN_DAY * (pay_long - fixing_long)) * brownians[1][3];
				
			payoff = ((f2 > f1 - K1) && (f2 < f1 + K2))? 1.0: 0.0;

			if (in_out_flag == SRT_OUT)
				payoff = 1.0 - payoff;

			payoff *= (f2 + margin);

			*ans += payoff;
		}

		*ans *= cvgpay * discpay / NSOBOL;
	
		if (err = sobol_free())
			return serror ("error in srt_f_sprrf_tok_mc: sobol free error");

		free_dmatrix (brownians, 1, 1, 1, 3);
	}
	else
	{
		payoff = ((CMSslope < K2) && (CMSslope > -K1))? 1.0: 0.0;

		if (in_out_flag == SRT_OUT)
				payoff = 1.0 - payoff;
		
		*ans = payoff * CMSpay * cvgpay * discpay;
	}

	return NULL;
}

/******************************************************************************/
