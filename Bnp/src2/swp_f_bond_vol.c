/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SWP     SWAP, Fixed Income Addins                     */
/*      SUB_SYSTEM:     BDT     BOND TOOLS                                    */
/*                                                                            */
/*      MODULE NAME:    SWP_F_BOND_VOL                                        */
/*                                                                            */
/*      PURPOSE:        EQUIVALENT BLACK-SCHOLES VOLATILITY                   */
/*						FROM BETA VOLATILITY								  */
/*                                                                            */
/*      AUTHORS:		FIRST					                        	  */
/*                                                                            */
/*      DATE:                                                                 */
/*                                                                            */
/*      VERSION:        01                                                    */
/*                                                                            */
/*      DESCRIPTION:    XX                                                    */
/*                                                                            */
/*      FUNCTIONS USED: XXX_X_XXXXXXXXX                                       */
/*              Must include all imported function call made by the module    */
/*                                                                            */
/*      PARAMETERS:     <not applicable>                                      */
/*                                                                            */
/*      RETURNS:                                                              */
/*                                                                            */
/*      DATA ACCESSED:                                                        */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/***  include files  **********************************************************/
#include "math.h"   
#include "num_h_allhdr.h"             
#include "swp_h_all.h"             

/*** declare static variables  ************************************************/
  

/***  functions  **************************************************************/

Err swp_f_bond_vol(	Date 		fut, 
					SwapDP 		p, 
					double 		coupon, 
					double 		clean_price,
					double		strike,
					double		volatility,
					Date		today,
					double		alpha,
					double		beta,
					double		rho,
					String 		swap, 
					String 		repo, 
					double 		redemption, 
					double		first_coupon,
					double		*implied_vol) 
{
double fwd_clean_price;
double PsiY, PsiYY, PsiYYY;
double yield, yieldmean, yieldstar;
double Tex; /* maturity of the option in number of years */
double pmean;

Err err = NULL;
SwapDP pf;
DateList list;      
double fut_first_coupon;
SimpleSwap s;


	/* Computation forward dirty price and dirty strike */
	pf = p;   /*Mainly to have pf.end = p.end :
			this alllows us to define a proper BKWD direction
			for the forward bond dates and then get the exact
			forward accrued interests */
	pf.start = fut;
	pf.first_full_fixing = fut;
	pf.direction = BKWD;               
	list = SwapDP_to_DateList(&p,NO_BUSDAY_CONVENTION);
	fut_first_coupon = (fut<list.date[1])?first_coupon:coupon;

	fwd_clean_price = fwd_clean_price_fct(fut, p, coupon, clean_price, swap, repo, first_coupon);
					
	pmean = (fwd_clean_price + strike) / 2.;

	/* Computation of yield, yieldmean and yieldstar */
	yield     = yield_fct(pf,coupon, fwd_clean_price, redemption, fut_first_coupon);
	yieldmean = yield_fct(pf,coupon, pmean, redemption, fut_first_coupon);
	yieldstar = yield_fct(pf,coupon, strike, redemption, fut_first_coupon);

	/* Computation of PsiY, PsiYY, PsiYYY */
	s = make_SimpleSwap(&pf,coupon,0.,1.,pf.start,BOND);
	if(err=duration_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			s.sdp.basis_code,
			s.sdp.compd,
			yieldmean,
			&PsiY))
		return err;
  	if(err=convexity_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			s.sdp.basis_code,
			s.sdp.compd,
			yieldmean,
			&PsiYY))
		return err;
  	if(err=thirdmoment_range(
			s.dl.date, 
			s.cpn.d,
			s.cpn.len,
			s.dl.type,
			s.dl.prev,
			s.sdp.basis_code,
			s.sdp.compd,
			yieldmean,
			&PsiYYY))
		return err;
	free_inSimpleSwap(&s);
	
	/* OK now what is the implied volatility ? */
	Tex = (add_unit(fut, -p.spot_lag, SRT_BDAY, NO_BUSDAY_CONVENTION) - today) / 365.0;

	err = swp_f_bond_vol_yldbeta_to_pricebs(fwd_clean_price, strike, 
											yield, yieldmean, yieldstar,
											- PsiY, PsiYY, - PsiYYY, 
											volatility, Tex, 
											alpha, beta, rho,
											implied_vol);
	
	return err;
}


Err swp_f_bond_vol_yldbeta_to_pricebs(	double 		fwd_clean_price,
										double		clean_strike,
										double		yield,
										double		yieldmean,
										double		yieldstar,
										double		PsiY,
										double		PsiYY,
										double		PsiYYY,
										double		volatility,
										double		optmat,
										double		alpha,
										double		beta,
										double		rho,  
										double		*implied_vol) 
{
double pmean, gamma1, gamma2, z, xz;
double integrand;
double ratio_z_xz, ratio_integrand;
double CFactor;

	/* Computation of pmean */
	pmean = (fwd_clean_price + clean_strike) / 2.;

	/* Computation of z, xz, gamma1, gamma2 */
	z = alpha / volatility * (fwd_clean_price - clean_strike) / PsiY / pow(yieldmean,beta);
	
	if ((1.-rho) < DBL_EPSILON*100.)
		xz = 1.;
	else
		xz = log((sqrt(1. - 2.*rho*z + z*z) + z - rho) / (1.-rho));
	
	if (fabs(z) < DBL_EPSILON*100.)
		ratio_z_xz = 1.;
	else
		ratio_z_xz = z / xz;

	gamma1 = PsiYY / (PsiY*PsiY) 
			+ beta / (PsiY*yieldmean);
	
	gamma2 = PsiYYY / pow(PsiY,3) 
			- pow(PsiYY,2) / pow(PsiY,4) 
			+ PsiYY / pow(PsiY,3)*beta/yieldmean 
			+ beta*(beta-1)/pow(yieldmean,2)/pow(PsiY,2); 
	
	/* OK now what is the implied volatility ? */
	if ((1.-beta) < DBL_EPSILON*100.)
		integrand = log(yield / yieldstar);
	else
		integrand = (pow(yield,1.-beta) - pow(yieldstar, 1.-beta)) / (1.-beta);

	if (fabs(integrand) < DBL_EPSILON*100.)
		ratio_integrand = PsiY * pow(yield, beta) / fwd_clean_price;
	else
		ratio_integrand = log(fwd_clean_price / clean_strike) / integrand; 

	CFactor = (2.*gamma2 - gamma1*gamma1 + 1./pmean/pmean)/24.
					+ 1./4.*rho*gamma1*alpha/volatility/PsiY/pow(yieldmean,beta)
					+ (2.-3.*rho*rho)/24.*alpha*alpha/pow(volatility*PsiY*pow(yieldmean,beta),2) ;

	(*implied_vol) = volatility * ratio_integrand * ratio_z_xz;
	(*implied_vol) *= (1.+ CFactor * pow(volatility * PsiY * pow(yieldmean,beta),2) * optmat);
	
	/* when yield is going up, price is going down */
	/* we have a transformation from a yield vol to a price vol */
	/* (so we have to take the negative value of the preceding volatility) */
	(*implied_vol) *= -1;

	return NULL;
}


Err swp_f_bond_vol_yldbeta_to_yldbs(	double 		fwd_clean_price,
										double		clean_strike,
										double		yield,
										double		yieldmean,
										double		yieldstar,
										double		PsiY,
										double		PsiYY,
										double		PsiYYY,
										double		volatility,
										double		optmat,
										double		alpha,
										double		beta,
										double		rho,  
										double		*implied_vol) 
{
double pmean, gamma1, gamma2, z, xz;
double ratio_z_xz, ratio_integrand;

	/* Computation of pmean */
	pmean = (fwd_clean_price + clean_strike) / 2.;

	/* Computation of z, xz, gamma1, gamma2 */
	z = alpha / volatility * (fwd_clean_price - clean_strike) / PsiY / pow(yieldmean,beta);
	
	if ((1.-rho) < DBL_EPSILON*100.)  // threshold has been increased to 100*epsilon to avoid numerical problems
		                              // on prices that are very close to ATM
		xz = 1.;
	else
		xz = log((sqrt(1. - 2.*rho*z + z*z) + z - rho) / (1.-rho));
	
	if (fabs(z) < DBL_EPSILON*100.)
		ratio_z_xz = 1.;
	else
		ratio_z_xz = z / xz;
	
	if (fabs(yield - yieldstar) < DBL_EPSILON*100.)
		ratio_integrand = 1. / pow(yield, 1.-beta);
	else
	{
		if ((1.-beta) < DBL_EPSILON*100.)
			ratio_integrand = 1.;
		else
			ratio_integrand = log(yield / yieldstar) / ((pow(yield,1.-beta) - pow(yieldstar, 1.-beta)) / (1.-beta));
	}

	gamma1 = PsiYY / (PsiY*PsiY) 
			+ beta / (PsiY*yieldmean);
	
	gamma2 = PsiYYY / pow(PsiY,3) 
			- pow(PsiYY,2) / pow(PsiY,4) 
			+ PsiYY / pow(PsiY,3)*beta/yieldmean 
			+ beta*(beta-1)/pow(yieldmean,2)/pow(PsiY,2); 
	
	/* OK now what is the implied volatility ? */
	*implied_vol = volatility * ratio_integrand * ratio_z_xz
				* (1.+ 
				  ( pow((1. - beta),2)/24.*pow(volatility,2)*pow(yieldmean,2.*(beta -1.))
					+ 1./4.*rho*gamma1*alpha*volatility*PsiY*pow(yieldmean,beta)
					+ (2.-3.*rho*rho)/24.*alpha*alpha ) * optmat);

	return NULL;
}

