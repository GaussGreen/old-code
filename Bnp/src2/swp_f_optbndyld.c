/*
 * module: swp_f_optbndyld.c
 * author: E.Auld Nov 95
 * this module has two public functions:
 * swp_f_optbndyld: value bond option numerically assuming yield is normally
 *     or lognormally distributed
 * swp_f_bndcms: compute bond convexity adjustment necessary given above model
 *
 * see Tokyo Swaps Pricing Methodologies
 *     E.Auld and A.Tsuchiya
 */
 

/*******************************************************************************
**  Include Files   
*******************************************************************************/


#include "math.h"
#include  <swp_h_all.h"
#include  <num_h_allhdr.h"
#include  <swp_h_optbndyld.h"

#define OPTBNDYLDMAXITER		20
#define OPTBNDYLDTOLERANCE		(1.0e-7)
#define OPTBNDYLDSIMPSONTOL		(1.0e-8)
#define OPTBNDYLDMAXINTEGLIMIT	50
#define OPTBNDYLDMINSTRIKE		(.0001)
#define MAXYLD					(0.6)
#define NUMGLOBALBONDS			2

#undef DBGBNDCMS 

/*************************************************************

NI stands for numerical integration
use static global variables because we call numerical integrating
routine sm_qsimp.
******************************************************************/

static  double				NI_sig;
static  double				NI_a;
static  double				NI_k;
static  double				NI_logy0;
static  double				NI_ext_adj1;
static  double				NI_ext_adj2;
static  int					NI_ext_flg;
static  SrtDiffusionType	NI_lognrm_nrm;
/* to represent bonds for numerical integration functions */
static DateList				glob_dl[NUMGLOBALBONDS];
static SwapDP				glob_sdp[NUMGLOBALBONDS];
static double				NI_cpn[NUMGLOBALBONDS];
static double				NI_cf[NUMGLOBALBONDS];
static FILE	*				fptr;


/*
 * local function to convert yield to price
 * use local function to avoid recomputing dates
 * every time
 */

static double local_price_fct (double irr, int i) 
{
double   price, disc1;
double   dr;
static double redemption = 1.0;

	disc1 = 1/(1+irr/(double)glob_sdp[i].compd);
  
	if(glob_sdp[i].basis_code==BASIS_ACT_USD)  
	{ 
		dr=day_count_date(glob_dl[i].date[0],glob_dl[i].date[1],glob_sdp[i].basis_code);
		dr=dr/day_count_date(glob_dl[i].prev,glob_dl[i].date[1],glob_sdp[i].basis_code);    
		dr=dr/glob_sdp[i].compd;
	}
	else
	{
		dr=coverage(glob_dl[i].date[0],glob_dl[i].date[1],glob_sdp[i].basis_code);
	}
  
	if(disc1 != 1.0)
	{
		price = (NI_cpn[i]/(double)glob_sdp[i].compd)*
				(pow(disc1,glob_dl[i].len-(glob_dl[i].type==BROKEN))-disc1)
				/(disc1-1)      
				+ redemption*pow(disc1,glob_dl[i].len-1-(glob_dl[i].type==BROKEN));
	}
	else
	{ 
		price = (NI_cpn[i]/(double)glob_sdp[i].compd)
				*(glob_dl[i].len-1-(glob_dl[i].type==BROKEN))
				+ redemption;
  
	}
  
	if(glob_dl[i].type == BROKEN)
	{
		price =	(price+NI_cpn[i]/(double)glob_sdp[i].compd)
				* pow(disc1,glob_sdp[i].compd*dr) ;
	}

	return(price);   
}  


/*
 * function to numerically integrate
 * = value of bond option contingent on log of yield at maturity
 * times density
 * sometimes we do extinguishable as well
 */
static double bond_opt_integrand (double logy)
{  
double yield;
double payoff;  
double price;
double dens;

  yield = NI_logy0 + NI_sig*logy;

  if(NI_lognrm_nrm == SRT_LOGNORMAL)
  {
    yield = exp(yield);
  }

  price = local_price_fct(yield,0);

  /*
   * we don't need to take max because we are only integrating in
   * in the area where there is positive payoff
   */

  payoff = NI_a*(price - NI_k);
  dens = gauss(logy);
  
  if(NI_ext_flg != 0)
  {
    dens -= INV_SQRT_TWO_PI*exp(NI_ext_adj2-(logy-NI_ext_adj1)*(logy-NI_ext_adj1)/ 2.0);
  }
  
  payoff *= dens;
  
  return payoff;
}

/******************************************************************************/

               

/*-------------------------------------------------------------------------*/

/* PUBLIC FUNCTIONS IN THIS MODULE: */

/*
 * swp_f_bndcms
 * give adjusted bond yield
 * to value option on bond yield
 * Works by iteratively calling swp_f_optbndyld
 * valueing a call bond option with zero price strike until
 * we get back the right number
 */

SrtErr swp_f_bndcms  (  
  Date				today,
  Date				fixdt,
  Date				paydt,
  SwapDP *			sdp,
  double			cpn,
  double			vol,
  double			fwd,
  SrtDiffusionType	lognrm_nrm,
  double *			answer
  )
{
SrtErr			err;
double			strike = 0;
SrtCallPutType	call_put = SRT_CALL;
SrtGreekType	greek = PREMIUM;
int niter, ext_flg = 0;
double x0,x1,y0,y1,dydx;
double fwd_pr;
double ext_barrier = 1000;
double mat;

	#ifdef DBGBNDCMS
		fprintf(fptr,"bndcms...\n");
	#endif

	/*
	* need to compute end date because bond code cannot handle nfp
	*/

	sdp->end = add_unit(sdp->start,12/(int)sdp->compd * sdp->nfp,
						SRT_MONTH, NO_BUSDAY_CONVENTION);

	/*             
	* find yield equivalent of fwd if it is a fwd price
	*/
	if(fwd > MAXYLD)
	{
		fwd_pr = fwd;
		fwd = yield_fct(*sdp,cpn,fwd,1.0,cpn);
	}
	else
	{
		fwd_pr = clean_price_fct(*sdp,cpn, fwd, 1.0, cpn);
	}


	/*
	* newton rhapson iteration
	* want to find x such that swp_f_optbndyld(x) = y(x) = 
	*/  
	x0 = fwd;
	x1 = fwd + .001;

	err = swp_f_optbndyld(today,fixdt,sdp,cpn,strike,vol,x0,ext_barrier,x0,call_put,
							ext_flg,greek,lognrm_nrm,&y0);
	if(err) return err;
	
	y0 -= fwd_pr;

	err = swp_f_optbndyld(today,fixdt,sdp,cpn,strike,vol,x1,ext_barrier,x1,call_put,
							ext_flg,greek,lognrm_nrm,&y1);
	if(err) return err;
  
	y1 -= fwd_pr;                         

	niter = 0;
	while(niter < OPTBNDYLDMAXITER && fabs(y1) > OPTBNDYLDTOLERANCE)
	{
		#ifdef DBGBNDCMS
			fprintf(fptr,"OPTBNDCMS fpr %lf niter %d x0 %lf y0 %lf x1 %lf y1 %lf\n",
					fwd_pr,niter,x0,y0,x1,y1);
		#endif

		dydx = (y1 - y0)/(x1-x0);
		y0 = y1;
		x0 = x1;
		x1 = x1 - y1/dydx;
		
		err = swp_f_optbndyld(today,fixdt,sdp,cpn,strike,vol,x1,ext_barrier,x1,call_put,
								ext_flg,greek,lognrm_nrm,&y1);
		if(err) return err;
    
		y1 -= fwd_pr;
		
		niter++;
	}

	/* 
	* if the CMS is payed late, there will be another convexity effect
	* due to correlation between the zero coupon yield to paydt and the
	* yield of the bond.  These adjustments are just "Sven Effect" delayed
	* digital adjustments.
	*/
	if(paydt > fixdt)
	{
	    mat = (fixdt - today) * YEARS_IN_DAY;
	    vol *= sqrt(mat);
	    mat = (paydt - fixdt) * YEARS_IN_DAY;
	
		if(lognrm_nrm == SRT_LOGNORMAL)
		{
			x1 = x1*exp(-x1*vol*vol*mat);
		} 
		else
		{
			x1 = (x1 - vol * vol * mat);
		}
	}     

	*answer = x1;
	
	return NULL;
}


/*
 * swp_f_optbndyld
 * returns price of bond option if yield was 
 * normal/lognormal
 * doesn't understand anything about accrued interest
 * everything is on a notional of 1
 * maturity of option is fixdt
 * in case of extinguishable, only prices down and out bond options
 * in this case KO if spot_yld is ABOVE extinguishing level
 * FWD should already by adjusted for any yield/delayed payment effects.
 * value returned by this formula must be multiplied by discount factor to 
 * payment date.
 */
SrtErr swp_f_optbndyld  (  
  Date				today,
  Date				fixdt,
  SwapDP *			sdp,
  double			cpn,
  double			strike,
  double			vol,
  double			fwd,
  double			ext_barrier,
  double			spot_yld,
  SrtCallPutType	call_put,
  int				ext_flg,
  SrtGreekType		greek,
  SrtDiffusionType	lognrm_nrm,
  double *			answer
  )
{ 
double (*f)(double);
double lower, upper, log_kyld, log_byld;
double mat;

	/*
	* need to compute end date because bond code cannot handle nfp
	*/

	sdp->end = add_unit(sdp->start,12/(int)sdp->compd * sdp->nfp,
						SRT_MONTH,NO_BUSDAY_CONVENTION);
	glob_sdp[0] = *sdp;

	/*
	* find yield equivalents if inputs are prices
	*/
	if(fwd > MAXYLD)
	{
	    fwd = yield_fct(glob_sdp[0],cpn,fwd,1.0,cpn);
	}
	
	if((ext_flg != 0) && (spot_yld > MAXYLD))
	{
	    spot_yld = yield_fct(glob_sdp[0],cpn,spot_yld,1.0,cpn);
	}
	
	if((ext_flg != 0) && (ext_barrier > MAXYLD))
	{
		ext_barrier = yield_fct(glob_sdp[0],cpn,ext_barrier,1.0,cpn);
	}
  
	if((ext_flg != 0) && (ext_barrier <= spot_yld))
	{
		return serror("already hit barrier");
	}
  
	if(fixdt > sdp->start)
	{
	    return serror("bad fix date");
	}

	/*
	 * global variables for numerical integration
	*/

	mat = (fixdt-today)*YEARS_IN_DAY;
	NI_lognrm_nrm = lognrm_nrm;
	NI_a = 1.0 - 2.0 * (double)call_put;
	NI_sig = vol * sqrt(mat);
  

	if(lognrm_nrm == SRT_LOGNORMAL)
	{
	    NI_logy0 = log(fwd) - 0.5 * NI_sig * NI_sig;
	}
	else
	{
	    NI_logy0 = fwd;
	}
	
	NI_k = strike;
	NI_cpn[0] = cpn;
	
	if(ext_flg != 0)
	{
		if(lognrm_nrm == SRT_LOGNORMAL)
		{
			NI_ext_adj1 = 2 * log(ext_barrier/spot_yld) / NI_sig;
			NI_ext_adj2 = (2.0 * (NI_logy0 - log(spot_yld)) / 
							 (NI_sig * NI_sig) * 
							 log(ext_barrier/spot_yld));
			NI_ext_flg = ext_flg;
		}
		else
		{
			NI_ext_adj1 = 2 * (ext_barrier-spot_yld) / NI_sig;
			NI_ext_adj2 = (2.0 * (NI_logy0 - spot_yld) / 
							(NI_sig * NI_sig) * 
							(ext_barrier - spot_yld));
			NI_ext_flg = ext_flg;
		}
	}
	else
	{
	    NI_ext_flg = 0;
	}


	/*
	* log yield equivalent of strike and barrier
	*/
	if(strike > OPTBNDYLDMINSTRIKE)
	{
	    log_kyld = yield_fct(glob_sdp[0],cpn,strike,1.0,cpn);
	}
	else
	{
	    log_kyld = 1000.0;
	}
	
	if(lognrm_nrm == SRT_LOGNORMAL)
	{
	    log_kyld = (log(log_kyld) - NI_logy0)/NI_sig;
	}
	else
	{
	    log_kyld = (log_kyld - NI_logy0)/NI_sig;
	}

	if(ext_flg != 0)
	{
	    if(lognrm_nrm == SRT_LOGNORMAL)
	    {
			log_byld = (log(ext_barrier) - NI_logy0)/NI_sig;
		}
		else
		{
			log_byld = (ext_barrier - NI_logy0)/NI_sig;
		}
	}

	/*
	* get the dates so we don't have to get them again and again
	* during the integration
	*/
	glob_dl[0] = SwapDP_to_DateList(&glob_sdp[0],NO_BUSDAY_CONVENTION);
	
	f = bond_opt_integrand;
	
	/*
	* limits of integration
	*/

	if((call_put == SRT_CALL) && (ext_flg == 0))
	{
	    lower = DMIN(-10.0,-20.0+log_kyld);
	    upper = log_kyld;
	}
	else if((call_put == SRT_CALL) && (ext_flg == 1))
	{
	    lower = DMIN(-10.0,-20.0+log_kyld);
	    upper = DMIN(log_byld,log_kyld);
	}
	else if((call_put == SRT_PUT) && (ext_flg == 0))
	{
	    lower = log_kyld;
	    upper = DMAX(10.0,20.0+log_kyld);
	}
	else
	{
	    lower = log_kyld;
	    upper = log_byld;
	}

	if(lower < upper)
	{
	    lower = DMAX(lower, -1 * OPTBNDYLDMAXINTEGLIMIT);
	    upper = DMIN(upper, OPTBNDYLDMAXINTEGLIMIT);
	    *answer = sm_qsimp(f,lower,upper,OPTBNDYLDSIMPSONTOL);
	}
	else
	{
	    *answer = 0.0;
	}

	srt_free(glob_dl[0].date);
	
	if(greek != PREMIUM )
		return serror("greek not yet implemented");

	return NULL;
}