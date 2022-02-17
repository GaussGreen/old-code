 /****************************************************************************/
/* Purpose */
/****************************************************************************/

/****************************************************************************/
/* EXTERNAL ROUTINES USED */
/****************************************************************************/
/*

*/
/****************************************************************************/
/* Headers */
/****************************************************************************/
#include "srt_h_all.h"
#include "opfnctns.h"
#include "srt_h_lgmUStypes.h"
#include "srt_h_lgmUSprotos.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmeval.h"
#include "math.h"

/****************************************************************************/
/* Definition */
/****************************************************************************/

#define SRT_MIN_REDUC_FACT 1.10


/****************************************************************************/
/* Prototypes of private functions */
/****************************************************************************/

/* static counter_getvolfudgetimeswap_new; */		/* for debug */

/* calculates the integer part and remainder of x */
static long IntPart(double x, double* remainder);

/* Interpolates value at k+dk from an array of values for convolver */
static double InterpfromArray(long k, double dk, long nk, double *TheArray);

/* same as before, but extrapolates flat */
static double InterpfromArrayTS(long k, double dk, long nk, double *TheArray);

/* Computes integration weights for convolver */
static void GenConvWeights(double *w, double h, long m);

/* Calculates the integration weights for the kink correction */
static void GetCorrection(double* w, double x);




LGMErr GetVol_SABRValues(
						long	tNow,
						long	value_date,
						long	start,
						long	end,
						SrtCallTimeSwapPtr deal,
						double	*alpha, double *beta, double *sigma, double* rho);

LGMErr GetVol_FudgeTimeSwap_new(
						long	tNow,
						long	value_date,
						long	start,
						long	end,
						double	strike,
						double	fra,
						SrtCallTimeSwapPtr deal,
						double	*result,
						double  alpha,
						double  beta,
						double  sigma,
						double  rho);


/****************************************************************************/
/**** Comments ****/
/* Routines to be completed later		*/
/* LGMErr ConvEvalSimMidAtwBarrier()	*/

/**** HOLIDAYS ***/
/* Wherever we need to use the add_unit function, 
the currency code ccy is available */

/* Multiple integration that computes the price of one option  */
void NewConvolverTimeSwap(
							long	nEx,	/* number of exercises */
							long	nx,		/* number of values for the state variable */
							long	n0,
							long	nz,
							double	dx,
							double	h,
							long	m,
							long	killkinks,
							double	*reduction,	/*[0..nEx-1] */
							double	*weights,	/* [0..nx] */
							double	*varr, /*[0..nEx-1] */
							double	**payofftable,
							double	*answer
							)
{
	long integershift, knew, ipos;
	double	*Qplus = NULL;
	double	*Qarr = NULL;
	double	*differ = NULL;
	double	*payoff = NULL;
	double a1, a2, a3, exactkink, xkink, y, theta, kshift;
	double alpha, gamma, Qinterp, Pinterp;
	long	i, j, k;
	long kk, ifoundkink, icor;
	double ratio;
	double Cweight[4];
	double	corr;

	/* Memory allocation */
	Qplus = dvector(0,nx);
	Qarr = dvector(0,nx);
	differ = dvector(0,nx);
	
	/* Initialization of the solution */
	for (k=0; k<=nx; k++)									/* Value is zero after last exercise has passed */
		Qarr[k]=0;
	
	/* Backward induction */
	for (j=(nEx-1); j>=0; j--)								/* For each exercise date (starting from the last) */
	{	
		if (j<nEx-1 && reduction[j]<0.90*reduction[j+1])	/* reduce grid size if needed */
		{	
			ratio = reduction[j]/reduction[j+1];
			for (knew=0; knew<=nx; knew++)					/* transfer results */
				Qplus[knew] = Qarr[knew];
			for (knew=0; knew<=nx; knew++)					/* interpolate to get Qarr on new grid */
			{	
				k = IntPart(ratio*((double)(knew-n0)) + (double)n0, &gamma);
				Qarr[knew] = InterpfromArrayTS(k, gamma, nx, Qplus);
			}
		}
		payoff = payofftable[j];

		for (k=0; k<=nx; k++)
		{	Qplus[k] = 0;
			differ[k] = Qarr[k] - payoff[k];
		}
		if (varr[j]>0)										/* integrate */
		{	alpha = h*varr[j]/(dx*reduction[j]);			/* Qplus[k] = sum of weight[i]*max{Qarr[k*], payoff[k*]} */			
			for (i=0; i<=nz; i++)							/*  at k* = k + (i-m)*alpha */
			{	kshift = alpha * ((double) (i - m));		/* set up interpolation */
				integershift = IntPart(kshift, &gamma);
				for (k=0; k<=nx; k++)
				{	knew = k+integershift;					/* interpolate to get value at k+integershift+gamma */
					Qinterp = InterpfromArray(knew, gamma, nx, Qarr);
					Pinterp = InterpfromArray(knew, gamma, nx, payoff);
					Qplus[k] = Qplus[k] + weights[i]*max(Qinterp, Pinterp);	
				}												
			}
		}
		else
		{	for (k=0; k<=nx; k++)
				Qplus[k] = max(Qarr[k], payoff[k]);
		}

	/* Correct values for any kinks */
/* should use continue statements to straighten out this nest of vipers */
		ifoundkink = 0;													/* Search to find each crossing (kink) */
		for (kk=1; kk<=nx-2; kk++)			
		{	if ((differ[kk]>0 && differ[kk+1]>0) ||
				(differ[kk]<0 && differ[kk+1]<0) ||
				(differ[kk-1]>0 && differ[kk]==0 && differ[kk+1]>0) ||
				(differ[kk-1]<0 && differ[kk]==0 && differ[kk+1]<0) )
				continue; /* No kink go to next point */
			
			a1 = differ[kk];										/* otherwise there's a kink between kk, kk+1 */								
			a2 = differ[kk+1] - differ[kk];							/* locate it exactly */						
			a3 = 0.25*(differ[kk+2] + differ[kk-1] - differ[kk] - differ[kk+1]);
			a2 = a2 - a3;

			if (fabs(a3/a2)<1.e-06)								
				exactkink = -(a1/a2)*(1.0 + a1*a3/(a2*a2));
			else if (a2>0)
				exactkink = (-a2 + sqrt(a2*a2-4.0*a1*a3))/(2.0*a3);
			else
				exactkink = (-a2 - sqrt(a2*a2-4.0*a1*a3))/(2.0*a3);	
			exactkink = exactkink + (double)kk;						/* kink is at exactkink */

			xkink = (exactkink - (double)n0)*dx*reduction[j];		/* store kink for output */
			ifoundkink = 1;

			if (killkinks == 0 || varr[j] <=0 || alpha <= 1.e-4)
				continue; /* continue to the next point unless we need to kill the kink */

			for (k=0; k<=nx; k++)								/* Correct Qplus[k] for kink */
			{	y = (exactkink-(double)k)/alpha + (double)m;
				i = IntPart(y, &theta);

				if (i>=2 && i<=nz-3)		
				{	GetCorrection(Cweight, theta);
					corr = 0.0;
					kshift = ((double)k) + alpha*((double)(i-1-m));
			  		for (icor=0; icor<=3; icor++)
					{	knew = IntPart(kshift, &gamma);
						Qinterp = fabs(InterpfromArray(knew, gamma, nx, differ));
						ipos = i-1+icor;
						corr = corr + Cweight[icor]*weights[ipos]*Qinterp;
						kshift = kshift + alpha;
					}
					Qplus[k] = Qplus[k] + corr;
				}
			}											/* end k loop */
		}												/* end kk loop */

	

		for (k=0; k<=nx; k++)							/* update Qarray */
			Qarr[k] = Qplus[k];
	}


	*answer = Qarr[n0];

	/* Memory freeage */
	free_dvector(Qplus,0,nx);
	free_dvector(Qarr,0,nx);
	free_dvector(differ,0,nx);
}


LGMErr ComputeBermudanAndEuropeanTimeSwap(
/* info about convolutions */
		long		nEx,			/* number of exercises */
		double		*zeta,			/* [0, 1, ..., nEx-1], values of zeta at the exercise dates */
		ConvParams	*parms,			/* convolution numerical constants */

/* info about today's discount curve and swaption/cap vols */ 
		Date		EvalDate,		/* tNow for deal evaluation */
		String		ycName,			/* yield curve name */
LGMErr	(*GetVol)(Date, Date, double, SRT_Boolean, double*),	/* swaption/cap vols */
LGMErr	(*GetBeta)(Date, Date, double*),		/* swaption/cap exponents (beta) */
		LGM_TS		*tsPtr,			/* calibrated LGM term structure */

/* information about the deal */
		LGMDealType	dealtype,
		void		*dealPtr,
		LGMErr		(*payofffunc)(),

		LGMCalParm	*CalReq,

/* output */
		double		*answer,		/* value of the deal */
		double		**xExBdry)		/* array of exercise points (in x) */
{
/* declarations */
	double *varr, *x, **payofftable;
	double *vmax, *reduction;
	double *x_exerarr;
	double *weights;
	double	***payofftable_payers = NULL;
	double	***payofftable_receivers = NULL;
	double	*receiver_european = NULL;
	double	*payer_european = NULL;
	double gridwidth, maxh, stencil, h, dx;
	double	stdev_european;
	double var, totalvar, totalstddev, previouszeta;
	long nx, n0, nz, m, j, k;
	int index_european;
	long killkinks;
	double price;
	LGMErr	error;
	int integrateagain;
	SrtCallTimeSwap *ptr = NULL;

	ptr = (SrtCallTimeSwap*) dealPtr; 

	integrateagain=0;
startagain:
/* eliminate trivial case */
	if (nEx<1)
	{	*answer = 0.0;
		return (NULL);
	}
/* unload the convolution parameters and set them to their default values if needed */
	gridwidth = parms->gridwidth;
	if (gridwidth<2.0 || gridwidth>8.0)
		gridwidth=6.0;

	
	
	nx = parms->nx;
	if (nx<32 || nx>576)
		nx=192;		
	n0 = ((long)(nx/2));  /* make sure nx is even */
	nx = 2*n0;
	
	maxh = parms->h;
	if (maxh<0.02 || maxh>0.5)
		maxh=0.0625;

	stencil = parms->stencil;
	if (stencil<2.0 || stencil>8.0)
		stencil=6.0;

	killkinks = parms->killkinks;
	if (killkinks!=0)
		killkinks=1;

/* Calculate the v[j] array (standard deviations) */
	varr = NULL; vmax = NULL; reduction = NULL;
	varr = (double*) srt_calloc(nEx, sizeof(double));
	vmax = (double*) srt_calloc(nEx, sizeof(double));
	reduction = (double*) srt_calloc(nEx, sizeof(double));

	if (varr==NULL || vmax==NULL || reduction==NULL)
	{	srt_free(varr); srt_free(vmax); srt_free(reduction);
		return ("first allocation failed in convolver");
	}

	totalvar=zeta[nEx-1];
	if (totalvar<1.e-12)
	{	totalvar=1.e-12;					/* Clearly an intrinsic value calculation the hard way */
		for (j=0; j<nEx; j++)
			zeta[j] = totalvar*((double)(j+1))/((double) nEx);
	}
	totalstddev = sqrt(totalvar);

/* find varr[j], the std dev from tEx[j-1] to tEx[j] */
/* if the change in variance is non-positive, signal with varr[j] = -1 */
	previouszeta = 0.;
	for (j=0; j<nEx; j++)
	{	var = zeta[j] - previouszeta;
		if (var<1.e-7*totalvar)
		{	varr[j] = -1.0;
			zeta[j] = previouszeta;
		}
		else
			varr[j] = sqrt(var);
		previouszeta = zeta[j];
	}
	
	vmax[0] = max(0,varr[0]);				/* vmax[j] is max of varr[k] for k<=j */
	for (j=1; j<nEx; j++)
		vmax[j] = max(vmax[j-1], varr[j]);
	
/* STEP 1: Set up x and z grids, and allocate space */
/* x[k] = (k-n0)*dx*reduction[j], k=0,1,..., nx, j=0,1,2,...,nEx-1 */
	dx = 2.0*totalstddev*gridwidth/((double)nx);

/* z[i] = (i-m)*h, i = 0, 1, ..., 2m=nz */
/* integration intervals from z[i] to z[i+1], for i = 0, 1, ..., 2m-1 = nz-1 */
	h = 0.9233333*dx/vmax[nEx-1];			/* 0.9233 to prevent staircase from roundoff */
	h = min(h, maxh);
  
	m = ((long) (1.0 + stencil/h));		
	nz = 2*m;
 
/* find out where grid spacing can be reduced */
/* can be reduced if reduced grid
	a. keeps totalstddev standard deviations (from zero to tEx[j]) on grid
	b. keeps delta_z less than delta_x
*/
	reduction[nEx-1] = 1.0;
	for (j=nEx-2; j>=0; j--)
	{	reduction[j] = reduction[j+1];
		if (varr[j]>0)
		{	while (sqrt(zeta[j])<0.5*totalstddev*reduction[j]
					&& h<0.5*dx*reduction[j]/vmax[j])
				reduction[j] = 0.5*reduction[j];
		}
	}

/* Allocate arrays */

	x=NULL; 
	payofftable=NULL; weights=NULL;
	if ((xExBdry) && (*xExBdry!=NULL))
		srt_free(*xExBdry);

	x = (double*) srt_calloc(nx+1, sizeof(double));
	payofftable = dmatrix(0, nEx-1, 0, nx);
	weights = (double*) srt_calloc(nz+1, sizeof(double));
	x_exerarr = (double*) srt_calloc(nEx, sizeof(double));

	/* New memory allocations */
	payofftable_payers = f3tensor(0,nEx-1,0,0,0,nx);
	payofftable_receivers = f3tensor(0,nEx-1,0,0,0,nx);
	receiver_european = dvector(0,nEx-1);
	payer_european = dvector(0,nEx-1);

	if (x==NULL ||  
		payofftable==NULL || weights==NULL || x_exerarr==NULL ||
		payofftable_payers == NULL || payofftable_receivers == NULL|| 
		receiver_european == NULL ||
		payer_european == NULL)
	{	srt_free(varr); srt_free(vmax); srt_free(reduction);
		srt_free(x);
		srt_free(weights); srt_free(x_exerarr);
		free_dmatrix(payofftable,0, nEx-1, 0, nx);
		free_f3tensor(payofftable_payers,0,nEx-1,0,0,0,nx);
		free_f3tensor(payofftable_receivers,0,nEx-1,0,0,0,nx);
		free_dvector(receiver_european,0,nEx-1);
		free_dvector(payer_european,0,nEx-1);
		return("allocation failed in convolver");
	}
	
	if(xExBdry) *xExBdry = x_exerarr;

/* STEP 2: Generate weights and payoff table */
	GenConvWeights(weights, h, m);
	for (k=0; k<=nx; k++)
		x[k] = dx*((double)(k-n0));

/* call payoff function to set up calculation */
	
	error = (*payofffunc)(payofftable, nx, x, reduction,			
							dealtype, dealPtr, CalReq,				
							EvalDate, ycName, tsPtr, GetVol, GetBeta);	
	/* gets the payoff of all european options */
	for (index_european=0;index_european<nEx;index_european++)
	{
		for (k=0;k<=nx;k++)
		{
			payofftable_receivers[index_european][0][k] = payofftable[index_european][k];
			payofftable_payers[index_european][0][k] = - payofftable[index_european][k];
		}
	}


	if (ptr->PayRec == SRT_PAYER)	
	{
		for (index_european=0;index_european<nEx;index_european++)
		{
			for (k=0;k<=nx;k++)
			{
				payofftable_receivers[index_european][0][k] = - payofftable[index_european][k];
				payofftable_payers[index_european][0][k] = payofftable[index_european][k];
			}
		}	
	}
	
	/* Computes the values of the receivers */
	for (index_european = 0;index_european<nEx;index_european++)
	{
		stdev_european = sqrt(zeta[index_european]);
		NewConvolverTimeSwap(1,nx,n0,nz,dx,h,m,killkinks,&(reduction[index_european]),weights,
						&stdev_european,payofftable_receivers[index_european],&price);
		receiver_european[index_european] = price;
	}

	/* Computes values for the payers */

	for (index_european = 0;index_european<nEx;index_european++)
	{
		stdev_european = sqrt(zeta[index_european]);
		NewConvolverTimeSwap(1,nx,n0,nz,dx,h,m,killkinks,&(reduction[index_european]),weights,
						&stdev_european,payofftable_payers[index_european],&price);	
		payer_european[index_european] = price;
	}
		
	/* Call put parity for the forward */
	for (index_european=0;index_european<nEx;index_european ++)
		(ptr->tForwardTS)[1][index_european] = receiver_european[index_european]-payer_european[index_european];



/* Adjusts the payoff of the bermudan with the forward correction */
	
	if (ptr->PayRec == SRT_RECEIVER)
	{
	for (j=0;j<nEx;j++)
		for (k=0;k<=nx;k++)
			payofftable[j][k] += ((ptr->tForwardTS)[0][j] - (ptr->tForwardTS)[1][j]);
	}
	else
	{
	for (j=0;j<nEx;j++)
		for (k=0;k<=nx;k++)
			payofftable[j][k] -= ((ptr->tForwardTS)[0][j] - (ptr->tForwardTS)[1][j]);
	}
	
	
	/* Computes the price of the bermudan */
	NewConvolverTimeSwap(nEx,nx,n0,nz,dx,h,m,killkinks,reduction,weights,
						varr,payofftable,answer);	
												
/* free up space allocated */
	srt_free(varr); srt_free(vmax);
	srt_free(reduction);
	srt_free(x); 
	free_dmatrix(payofftable,0, nEx-1, 0, nx);
	srt_free(weights);
	free_f3tensor(payofftable_payers,0,nEx-1,0,0,0,nx);
	free_f3tensor(payofftable_receivers,0,nEx-1,0,0,0,nx);
	free_dvector(receiver_european,0,nEx-1);
	free_dvector(payer_european,0,nEx-1);
	if(!xExBdry) srt_free(x_exerarr);

/* debug */
	if (integrateagain==1)
		goto startagain;

	return (NULL);
}






/*****************************************************************************************/
/* Integer truncator for convolver */
static long IntPart(double x, double* remainder)
{	long intpart;
	intpart = (long)x;
	*remainder = x - (double)intpart;
	if (*remainder < 0)
	{	*remainder = *remainder + 1.0;
		intpart = intpart - 1;
	}
	return (intpart);
}


/*****************************************************************************************/
/* Interpolator for convolver
Function: Interpolates value at k+dk from an array of values
TheArray[k], k=0, ..., nk. */
static double InterpfromArray(long k, double dk, long nk, double *TheArray)
{	double interpval, x;

	if (k>0 && k<(nk-1))
		interpval = TheArray[k] + dk*(TheArray[k+1] - TheArray[k])
					- 0.25*dk*(1.0-dk)*(TheArray[k-1] - TheArray[k] - TheArray[k+1] + TheArray[k+2]);

	else if (k<=0)
	{	x = dk + (double)k;
		interpval = TheArray[0] + x*(TheArray[1] - TheArray[0]) 
								- 0.5*x*(1.0-x)*(TheArray[0] - 2.0*TheArray[1] + TheArray[2]);
	}
	
	else
	{	x = dk + (double)(k+1-nk);
		interpval = TheArray[nk-1] + x*(TheArray[nk] - TheArray[nk-1])
							- 0.5*x*(1.0-x)*(TheArray[nk] - 2.0*TheArray[nk-1] + TheArray[nk-2]);
	}

	return (interpval);
}

/*****************************************************************************************/
/* Interpolator for convolver
Function: Interpolates value at k+dk from an array of values
TheArray[k], k=0, ..., nk. For time swap, extrapolates flat; this solves a problem for very
small maturities  */
static double InterpfromArrayTS(long k, double dk, long nk, double *TheArray)
{	double interpval, x;

	if (k>0 && k<(nk-1))
		interpval = TheArray[k] + dk*(TheArray[k+1] - TheArray[k])
					- 0.25*dk*(1.0-dk)*(TheArray[k-1] - TheArray[k] - TheArray[k+1] + TheArray[k+2]);

	else if (k<=0)
	{	x = dk + (double)k;
		interpval = TheArray[0] + x*(TheArray[1] - TheArray[0]) 
								- 0.5*x*(1.0-x)*(TheArray[0] - 2.0*TheArray[1] + TheArray[2]);
		/* TRY FOR EXTRAPOLATION FLAT */
		interpval = TheArray[0];
	}
	
	else
	{	x = dk + (double)(k+1-nk);
		interpval = TheArray[nk-1] + x*(TheArray[nk] - TheArray[nk-1])
							- 0.5*x*(1.0-x)*(TheArray[nk] - 2.0*TheArray[nk-1] + TheArray[nk-2]);
		/* TRY FOR EXTRAPOLATION FLAT */
		interpval = TheArray[nk];
	}

	return (interpval);
}




/**************************************************************************************/
/* Computes integration weights for convolver */
static void GenConvWeights(double *w, double h, long m)
{	long i, nz;
	double z;
	double  sum, sum2, sum4, a, b;

	nz = 2*m;
	for (i=0; i<=m; i++)
	{	z = h*((double)(i-m));
		w[i] = h*LGMsafeGauss(z);
	}

	for (i=0; i<m; i++)
		w[nz-i] = w[i];

/* End effects */
	z = -h*((double) m);
	w[0] = w[0] * 11.0/24.0 + norm(z);
	w[2] = w[2]*25.0/24.0;
 	w[nz] = w[0];
 	w[nz-2] = w[2];

/* Renormalize */
	sum = 0.0;
	sum2 = 0.0;
	sum4 = 0.0;
 	for (i=0; i<=nz; i++)
	{	z = h*((double)(i-m));
		sum = sum + w[i];
		sum2 = sum2 + w[i]*z*z;
		sum4 = sum4 + w[i]*z*z*z*z;
	}

	a = (sum4 - sum2)/(sum4*sum - sum2*sum2);
	b = (sum - sum2)/(sum4*sum - sum2*sum2);

 	for (i=0; i<=nz; i++)
	{	z = h*((double)(i-m));
		w[i] = w[i]*(a + b*z*z);
	}
	return;
}

/********************************************************************/
/* routine to calculate the correction weights for the kink removal */
static void GetCorrection(double* w, double x)
{	double xx, yy;
	if(x>1.0)
		x = 1.0;
	if (x<0.0)
		x = 0.0;
	xx = x*x;
	yy = (1.0 - x) * (1.0 - x);
	w[0] = yy*(2.0 - yy)/24.;
	w[1] = (1.0 - yy*(13.0 + 2.0*x - 3.0*xx))/24.0;
	w[2] = (1.0 - xx*(12.0 + 4.0*x - 3.0*xx))/24.0;
	w[3] = xx*(2.0 - xx)/24.0;

	return;
}



/**********************************************************************/

static double GetNormalVol(
double			mat,
Date			start, 
Date			end,
double			fwd,
double			strike, 
SRT_Boolean			isCap,
LGMErr			(*GetVol) (Date, Date, double, SRT_Boolean, double*),	
LGMErr			(*GetBeta) (Date, Date, double*))
{
	LGMErr		error1;
	Err			error2;
	double		beta;
	double		inp_vol;
	double		norm_vol;
	double		std, nstd;
	double		bs_price;

	if (mat < 1.0/365)
	{
		mat = 1.0/365;
	}

	if (fwd <= 0.0 || strike <= 0.0)
	{
		return 1.0e-08;
	} 

	error1 = GetVol (start, end, fwd, isCap, &inp_vol);
	error1 = GetBeta (start, end, &beta);

	std = pow (fwd, beta) * inp_vol * sqrt (mat);
	nstd = (strike - fwd) / std;

	error1 = GetVol (start, end, strike, isCap, &inp_vol);

	if (fabs (beta) < 1.0e-04)
	{
		return inp_vol;
	}

	if (fabs (nstd) > 5.0)
	{
		strike = fwd + 5.0 * std;
	}

	/*	Here we assume that input is lognormal */

	bs_price = srt_f_optblksch(
		fwd,
		strike,
		inp_vol,
		mat,
		1.0,
		SRT_CALL,
		SRT_PREMIUM);

	error2 = srt_f_optimpvol(
		bs_price,
		fwd,
		strike,
		mat,
		1.0,
		SRT_CALL,
		SRT_NORMAL,
		&norm_vol);

	return norm_vol;
}



/* ---------------------------------------------------------------------------------------------------------------------------------- */
/* Computes the forward LOGNORMAL vols of LIBOR using a fudge for the autocal timeswap */
LGMErr GetVol_FudgeTimeSwap(
						long	tNow,				/* Today (for the yield curve) */
						long	value_date,			/* Time in the future when we want the vol */
						long	start,				/* Start of fra */
						long	end,				/* End of FRA */
						double	strike,				/* Strike */
						double	fra,				/* Fwd */
						SrtCallTimeSwapPtr deal,	/* pointer to CTS */
						double	*result)			/* output:  volatility */
{
/* Variable declaration */
	LGMErr error = NULL;
	double maturity,sigma,alpha,beta,rho,volstrike/*, dStdDev*/;
	long startfudge,endfudge;
	SABR_VOL_TYPE atmtype;


/* Computes maturity */
	maturity = (start-value_date+0.0)/365.0;

/* Sliding or converging */
	if ( deal->timeswapvolmethod==2 ) /* Converging:  tkaes vol from today's market assuming that the market in the future is shifted. */ 
	{
		startfudge = start-(value_date-tNow);
		endfudge = end-(value_date-tNow);
	}
	else if (deal->timeswapvolmethod==3) /* Sliding: takes vol from today's market independent of where we are evaluating it. */
	{
		startfudge = start;
		endfudge = end;
	}
	else
		return("Unknown value of timeswapvolmethod parameter");

/* imposes maturity > 0 */
	if (startfudge<tNow)
	{
		startfudge = tNow+1;
		endfudge += (tNow+1) - startfudge;
	}


/* Check to see if there is a vol cutoff */
/*	if ( deal->dMaxStd > 0.0 )
	{
		error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,SABR_ATM_LOG,&sigma);
		if (error)
			return error;
		dStdDev = fra * exp( sigma*sqrt(maturity) );
		if ( strike > fra * dStdDev )
			strike = fra* dStdDev;
		else if ( strike < fra / dStdDev )
			strike = fra / dStdDev;
	}
*/

/* which ATM Vol ? */
	switch (deal->atmvolmethod)
	{
	case 1: 
		atmtype = SABR_ATM_LOG;
		break;
	case 2: 
		atmtype = SABR_ATM_NORM;
		break;
	case 3: 
		atmtype = SABR_ATM_BETA;
		break;
	}

/* Compute SABR parameters */
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,atmtype,&sigma);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,ALPHA,&alpha);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,BETA,&beta);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,RHO,&rho);
	if (error)
		return error;
	
/* Computes the vol woth the SABR formula */


	if (fra<=0.0 || strike <=0.0) /* Adds this condition to avoid the crash of MAD */
		volstrike = 0.0;
	else
	{
		if (deal->atmvolmethod==1)
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_LOGNORMAL,SRT_LOGNORMAL,&volstrike);
		else if (deal->atmvolmethod==2)
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_NORMAL,SRT_LOGNORMAL,&volstrike);
		else if (deal->atmvolmethod==3)
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_BETAVOL,SRT_LOGNORMAL,&volstrike);
		else
			return("Unknown value of atmvolmethod parameter");
	}

/* eliminates crazy results */
	if (volstrike < 0.0001)
		volstrike=0.0001;
	if (volstrike>1)
		volstrike=1;

/* fudge to match the forward value of the underlying by using an input volatility shift. */
	if ( deal->shiftvol )
		volstrike += deal->shiftvol * (1-exp(-(value_date-tNow+0.0)/(deal->tCpnPay[0]-deal->tCpnStart[0]+0.0)));

	*result = volstrike;
	return error;

}


/* Computes the forward vols of LIBOR using a fudge for the autocal timeswap */
LGMErr GetVol_FudgeTimeSwap_new(
						long	tNow,
						long	value_date,
						long	start,
						long	end,
						double	strike,
						double	fra,
						SrtCallTimeSwapPtr deal,
						double	*result,
						double  alpha,
						double  beta,
						double  sigma,
						double  rho)
{
	LGMErr error = NULL;
	double maturity,volstrike;


	/* for debug */
/*	counter = clock();*/
	
	/* Computes maturity */
	maturity = (start-value_date+0.0)/365.0;
	
	/* Computes the vol woth the SABR formula */
	if (fra <=0.0 || strike <=0.0)
	{
		volstrike=0.0;
	}
	else
	{
		if (deal->atmvolmethod==1)
		{
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_LOGNORMAL,SRT_LOGNORMAL,&volstrike);
		}
		else
		if (deal->atmvolmethod==2)
		{
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_NORMAL,SRT_LOGNORMAL,&volstrike);
		}
		else
		if (deal->atmvolmethod==3)
		{
			error = srt_f_optsarbvol(fra,strike,maturity,sigma,alpha,beta,rho,SRT_BETAVOL,SRT_LOGNORMAL,&volstrike);
		}
		else
		{
			return("Unknown value of atmvolmethod parameter");
		}

	}
	/* eliminates crazy results */
	if (volstrike < 0.0001)
		volstrike=0.0001;
	if (volstrike>1)
		volstrike=1;

	/* fudge to match the forward value of the underlying*/
	if (! deal->shiftvol ==0)
	{
		volstrike += deal->shiftvol * (1-exp(-(value_date-tNow+0.0)/(deal->tCpnPay[0]-deal->tCpnStart[0]+0.0)));
	}

	*result = volstrike;

	/* for debug */
/*	counter_getvolfudgetimeswap_new += clock() - counter; */

	return error;
}

LGMErr GetVol_SABRValues(
						long	tNow,
						long	value_date,
						long	start,
						long	end,
						SrtCallTimeSwapPtr deal,
						double	*alpha, double *beta, double *sigma, double* rho)
{
	LGMErr error = NULL;
	long startfudge,endfudge;
	SABR_VOL_TYPE atmtype;

	/* Sliding or converging */
	if (deal->timeswapvolmethod==2)
	{
		startfudge = start-(value_date-tNow);;
		endfudge = end-(value_date-tNow);
	}
	else 
	if (deal->timeswapvolmethod==3)
	{
		startfudge = start;
		endfudge = end;
	}
	else
	{
		return("Unknown value of timeswapvolmethod parameter");
	}

	/* imposes maturity > 0 */
	if (startfudge<tNow)
	{
		startfudge = tNow+1;
		endfudge += (tNow+1) - startfudge;
	}

	/* which ATM Vol ? */
	switch (deal->atmvolmethod)
	{
	case 1: atmtype = SABR_ATM_LOG;
		break;
	case 2: atmtype = SABR_ATM_NORM;
		break;
	case 3: atmtype = SABR_ATM_BETA;
		break;
	}

	/* Compute SABR parameters */
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,atmtype,sigma);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,ALPHA,alpha);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,BETA,beta);
	if (error)
		return error;
	error = (*(deal->GetSABRParam))(startfudge,endfudge,deal->vol_id,RHO,rho);

	return error;
}
	








/* Computes the payoff of the time swap for all (time,state variable) */
LGMErr BerTimeSwapPayoff(double **payoff, long nx, double *x, double *reduc, 
				LGMDealType dealtype, void *dealPtr, LGMCalParm	*CalReq,
				Date tNow, String ycName, LGM_TS *tsPtr,
				LGMErr	(*GetVol)(Date, Date, double,SRT_Boolean, double*),	/* swaption/cap vols */
				LGMErr	(*GetBeta)(Date, Date, double*))		/* swaption/cap exponents (beta) */
{
	LGMErr error=NULL;
	SrtCallTimeSwap *ptr=NULL;
	long i,j,k,l,m;
	long j0,nEx,i0,nCpn,ifirst;	/* number of exos */
	long lag=2;
	double *zEx=NULL;
	double strikelow,strikelowshift,strikeup,strikeupshift;
	double vol_accrued_low,vol_accrued_up;				/* vols for the floating digital adjustment */
	double putlowshift,putlow,putlowshift_adj,putlow_adj;
	double callup,callupshift,callupshift_adj,callup_adj;
	double G,zeta;
	SrtCrvPtr yldcrv;
	long spotlag;
	long	Start,End,tExp,start_accrued,end_accrued;
	double	df_start_accrued,df_end_accrued;
	double adjustmentup,adjustmentlow,correl,fra_accrued;
	double proba_between,proba_between_adj,proba_below,proba_below_adj,proba_above,proba_above_adj;
	double ***expected_ratio=NULL;
	double ***expected_ratio_adj = NULL;
	double ***expected_ratio_low=NULL;
	double ***expected_ratio_high=NULL;
	double	***tfra_accrued = NULL;
	double dfPay,expected_accrued_ratio,expected_accrued_ratio_adj;
	double expected_accrued_ratio_low, expected_accrued_ratio_high;
	double	tfras_discr[11];
	double	***tvols_discret = NULL;
	double	**tfras_subperiod = NULL;
	double	**tcms_subperiod = NULL;
	double	*tmaturities_subperiod = NULL;
	double ***tvols_subperiod = NULL;
	long	*tindex_firstpositive = NULL;
	long	num_subdiscret_obsfreq;
	long	num_subdiscret_statvar;
	long	**index_discret_statvar=NULL;
	long	*index_discret_obsfreq = NULL;
	long	current_index_obsfreq,current_index_statvar;
	long observation_freq;
	double xx,yy,w1,w2,w3,w4;   /* weights */
	long	index_for_adj_cms;
	double	fra_for_adj_cms,vol_for_adj_cms,level;
	/* for debug 
	clock_t counter1,counter2,counter3,counter4,counter_observationdf1,counter_observationdf2;
	*/
	double  alpha, beta, sigma, rho, dfstart_init, dfend_init, dfstart, dfend, mult;
	double df_start_accrued_init, df_end_accrued_init, level_init, cvg, dfPay_init;
	double G1, G2;
	static int firstTime = 1;
/*	static */ double*** observationdf=NULL;


	/* checks it is a timeswap deal */
	if (dealtype != CallTimeSwap)
		return ("wrong deal type in BerTimeSwappayoff");

	/* Starts the time counter */
/*	counter1 = clock();
	counter3 = 0;
	counter_getvolfudgetimeswap_new = 0;*/

	/* gets spot date and yield curve */
	yldcrv = lookup_curve(ycName) ;
	spotlag = get_spotlag_from_curve(yldcrv) ; 

	/* Recasts as correct deal type */
	ptr = (SrtCallTimeSwap*) dealPtr; 


	/* Number of atual exercises and coupons */
	j0 = ptr->FirstEx;
	nEx = ptr->nEx - j0;
	i0 = ptr->iSet[j0];			/* cpn dates i = i0, ..., ptr.nCpn */
	nCpn = ptr->nCpn - i0;		/* number of cpns: i = 0, ..., nCpn */ 
	observation_freq = ptr->observation_freq;
	index_for_adj_cms = observation_freq-1;
	num_subdiscret_obsfreq = ptr->num_subdiscret_obsfreq;	
	num_subdiscret_statvar = 11;

	/*Payoff initialisation  */
	for (j=0; j<nEx; j++)
	{	for (k=0; k<=nx; k++)
			payoff[j][k] = 0.0;  /*Payoff value at the exercise date j and node k */
	}     /* End */

	

	/* Memory allocation */
	zEx = (double*) srt_calloc(nEx, sizeof(double));	 	/* zeta[j0+j], j=0,...,nEx-1 */
	expected_ratio = f3tensor(0,nEx-1,0,nCpn-1,0,nx);
	expected_ratio_adj  = f3tensor(0,nEx-1,0,nCpn-1,0,nx);
	expected_ratio_low = f3tensor(0,nEx-1,0,nCpn-1,0,nx);
	expected_ratio_high  = f3tensor(0,nEx-1,0,nCpn-1,0,nx);
	tfra_accrued = f3tensor(0,nEx-1,0,nCpn-1,0,nx);
	tvols_discret = f3tensor(0,num_subdiscret_obsfreq-1,0,num_subdiscret_statvar-1,0,3);   /* low,up,lowshift,upshift */
	tfras_subperiod = dmatrix(0,observation_freq-1,0,nx);
	tcms_subperiod = dmatrix(0,observation_freq-1,0,nx);
	tmaturities_subperiod = dvector(0,observation_freq-1);
	tvols_subperiod = f3tensor(0,observation_freq-1,0,nx,0,3);  /*low,up,lowshift,upshift */
	tindex_firstpositive = lngvector(0,observation_freq-1);
	index_discret_statvar = lngmatrix(0,observation_freq-1,0,num_subdiscret_statvar-1);
	index_discret_obsfreq = lngvector(0,num_subdiscret_obsfreq-1);

	/* Computes indexes of discretisation for the observation days */
	for (current_index_obsfreq=0;current_index_obsfreq<num_subdiscret_obsfreq;
		current_index_obsfreq++)
		{
			index_discret_obsfreq[current_index_obsfreq] = (int) (current_index_obsfreq*
				(observation_freq-1.0)/(num_subdiscret_obsfreq-1.0)) ;
		}

	if (zEx==NULL /* || dPay==NULL  */   ) 
	{	
		error = "alloc. failed in CallTimeSwap payoff";
		goto cleanup;
	}

		observationdf = (double***) srt_calloc(ptr->nCpn, sizeof(double**));
		for (i=0;i<ptr->nCpn;i++)
		{
			(observationdf)[i] = (double**) srt_calloc(observation_freq,sizeof(double*));
			for (j=0;j<observation_freq;j++)
			{
				(observationdf)[i][j] = (double*) srt_calloc(2,sizeof(double));
			}
		}
		firstTime = 0;





		ifirst = 0;	

			for (i=ifirst; i< nCpn; i++)		/*Loop on coupons */		
			{
				for (l=0;l<=observation_freq-1;l++)
				{
					
					Start = ptr->observationdays[i+i0][l][0];
					End = ptr->observationdays[i+i0][l][1];
					observationdf[i+i0][l][0] = swp_f_df(tNow,Start,ycName);
					observationdf[i+i0][l][1] = swp_f_df(tNow,End,ycName);
				}
			}
	

	for (j=0; j<nEx; j++)	/* Loop on exo date */
	{
		zeta = LGMZetaFromTS(ptr->tSet[j+j0], tNow, tsPtr);		/* Get zeta at exercise date */

		ifirst = ptr->iSet[j+j0] - i0;				
		for (i=ifirst; i< nCpn; i++)		/*Loop on coupons */		
		{	
			/* computes the four strikes */
			strikelow = ptr->barriers[i+i0][0];
			strikelowshift = strikelow + ptr->call_spread_low;
			strikeup = ptr->barriers[i+i0][1];
			strikeupshift = strikeup + ptr->call_spread_up;
			
			/* computes all maturities and fras with the dynamic of Autocal */
			for (l=0;l<=observation_freq-1;l++)
			{
				/* Computes Start and End for the observation day */
				Start = ptr->observationdays[i+i0][l][0];
				End = ptr->observationdays[i+i0][l][1];
				/* Computes the maturity */
				tExp = add_unit(Start,-lag,SRT_BDAY,SUCCEEDING);
				tmaturities_subperiod[l] = (double) (tExp-ptr->tEx[j+j0]+0.0)/365.0; 
				tindex_firstpositive[l] = -1;
				dfstart_init = observationdf[i+i0][l][0]; //swp_f_df(tNow,Start,ycName);
				G1 = LGMGFromTS(Start,tsPtr);
				dfstart_init *= exp(-0.5*G1*G1*zeta);

				dfend_init = observationdf[i+i0][l][1]; //swp_f_df(tNow,End,ycName);
				G2 = LGMGFromTS(End,tsPtr);
				dfend_init *= exp(-0.5*G2*G2*zeta);

				/* computes the coverage and level */
				cvg = coverage(Start,End,BASIS_ACT_360);


				
				for (k=0;k<=nx;k++)	/*state variable */
				{
				
					dfstart = dfstart_init*exp(G1*x[k]*reduc[j]);
					dfend   = dfend_init*exp(G2*x[k]*reduc[j]);
					
					level = cvg * dfend;

					/* Computes the fra */
					tfras_subperiod[l][k] = (dfstart-dfend)/level;
					
					/****************************/

					/* tracks the first positive fra */
					if (0<tfras_subperiod[l][k] && -1==tindex_firstpositive[l])
					{
						tindex_firstpositive[l] = k;
					}
				}
			}

			/* computes the cms */
		
				error = GetVol_SABRValues(tNow,
									  ptr->tEx[j+j0],
									  ptr->observationdays[i+i0][index_for_adj_cms][0],
									  ptr->observationdays[i+i0][index_for_adj_cms][1],
									  ptr,
									  &alpha, &beta, &sigma, &rho);
			if (error) return error;

		

			for (k=0;k<=nx;k++)
			{
				fra_for_adj_cms = tfras_subperiod[index_for_adj_cms][k];
				
				error = GetVol_FudgeTimeSwap_new(tNow,ptr->tEx[j+j0],
						ptr->observationdays[i+i0][index_for_adj_cms][0],
						ptr->observationdays[i+i0][index_for_adj_cms][1],
						fra_for_adj_cms,fra_for_adj_cms,ptr,
						&vol_for_adj_cms, alpha, beta, sigma, rho);
				if (error) return error;
			

				for (l=0;l<=observation_freq-1;l++)
				{
					if (tfras_subperiod[l][k]>0)
					{
						tcms_subperiod[l][k] = tfras_subperiod[l][k]*exp(
							vol_for_adj_cms*vol_for_adj_cms*fra_for_adj_cms*
							tmaturities_subperiod[l]*
							(ptr->observationdays[i+i0][l][1]-
							ptr->tCpnPay[i+i0]+0.0)/365.0);
					}
				}
			}

			/* computes the indexes of discretisation for the state variable */
			/* prestep: computes boundaries strikes */
			for (l=0;l<observation_freq;l++)
			{
				tfras_discr[0] = tfras_subperiod[l][tindex_firstpositive[l]];
				tfras_discr[10] = tfras_subperiod[l][nx];
				if (0.5*strikelow>=1.5*strikelow-0.5*strikeup)
				{
					tfras_discr[1] = 0.5*strikelow;
				}
				else
				{
					tfras_discr[1] = 1.5*strikelow-0.5*strikeup;
				}
				if (0.5*(strikeup+tfras_discr[10])<1.5*strikeup-0.5*strikelow)
				{
					tfras_discr[9] = 0.5*(strikeup+tfras_discr[10]);
				}
				else
				{
					tfras_discr[9] = 1.5*strikeup-0.5*strikelow;
				}
				for (m=2;m<=8;m++)
				{
					tfras_discr[m] = strikelow + (m-2.0)/6.0*(strikeup-strikelow);
				}

				/* now, really computes the indexes */
				
				m=0;
				k=0;
				while (k<=nx && m<=10)
				{
					if (tfras_subperiod[l][k]>=tfras_discr[m])
					{
						index_discret_statvar[l][m] = k;
						m+=1;
					}
					k++;
				}
				while (m<=10)
				{
					index_discret_statvar[l][m] = nx;
					m++;
				}
			}




			/* computes the vols at the discretisation points */
			for (current_index_obsfreq=0;current_index_obsfreq<num_subdiscret_obsfreq;current_index_obsfreq++)
			{
				for (current_index_statvar=0;current_index_statvar<num_subdiscret_statvar;current_index_statvar++)
				{
					l = index_discret_obsfreq[current_index_obsfreq];
					k = index_discret_statvar[l][current_index_statvar];
					/* Computes Start and End for the observation day */
					Start = ptr->observationdays[i+i0][l][0];
					End = ptr->observationdays[i+i0][l][1];
					
	error = GetVol_SABRValues(tNow,
									  ptr->tEx[j+j0],
									  Start,
									  End,
									  ptr,
									  &alpha, &beta, &sigma, &rho);
					if (error) return error;
					

					
					
					/* gets volup and volupshift */
					
					error = GetVol_FudgeTimeSwap_new(tNow,ptr->tEx[j+j0],Start,End,strikeup,
							tfras_subperiod[l][k],ptr,&(tvols_discret[current_index_obsfreq][current_index_statvar][1]),
							alpha, beta, sigma, rho);
					if (error) return error;

					error = GetVol_FudgeTimeSwap_new(tNow,ptr->tEx[j+j0],Start,End,strikeupshift,
							tfras_subperiod[l][k],ptr,&(tvols_discret[current_index_obsfreq][current_index_statvar][3]),
							alpha, beta, sigma, rho);
					if (error) return error;

					  if (0<strikelow)
					{	/* gets vollow and vollowshift */
					

						  error = GetVol_FudgeTimeSwap_new(tNow,ptr->tEx[j+j0],Start,End,strikelow,
							tfras_subperiod[l][k],ptr,&(tvols_discret[current_index_obsfreq][current_index_statvar][0]),
							alpha, beta, sigma, rho);
						if (error) return error;

						error = GetVol_FudgeTimeSwap_new(tNow,ptr->tEx[j+j0],Start,End,strikelowshift,
							tfras_subperiod[l][k],ptr,&(tvols_discret[current_index_obsfreq][current_index_statvar][2]),
							alpha, beta, sigma, rho);
						if (error) return error;



					  }
				}
			}

			/* COMPUTES ALL VOLS FOR THE PERIOD BY INTERPOLATION OF VOLS AT DISCRETISATION POINTS */
			current_index_obsfreq = 0;
			for (l=0;l<=observation_freq-1;l++)	/* Loop on the observation days */
			{
				if (l>index_discret_obsfreq[current_index_obsfreq+1])
				{
					current_index_obsfreq += 1;
				}
				if (tmaturities_subperiod[index_discret_obsfreq[current_index_obsfreq+1]] == 
					tmaturities_subperiod[index_discret_obsfreq[current_index_obsfreq]])
				{
					xx = 0.5;
				}
				else
				{
					xx = (tmaturities_subperiod[l] - tmaturities_subperiod[index_discret_obsfreq[current_index_obsfreq]])/
						(tmaturities_subperiod[index_discret_obsfreq[current_index_obsfreq+1]]-
						tmaturities_subperiod[index_discret_obsfreq[current_index_obsfreq]]);
				}

				current_index_statvar = 0;

				for (k=tindex_firstpositive[l];k<=nx;k++)		/* Loop on the state variable */
				{
					if (tfras_subperiod[l][k]>tfras_subperiod[l][index_discret_statvar[l][current_index_statvar+1]])
					{
						current_index_statvar += 1;
					}
					
					yy = (tfras_subperiod[l][k]-tfras_subperiod[l][index_discret_statvar[l][current_index_statvar]])/
						(tfras_subperiod[l][index_discret_statvar[l][current_index_statvar+1]]-tfras_subperiod[l][index_discret_statvar[l][current_index_statvar]]);

					/*  vol is a weighted average of vols at discretisation points */
					w1 = (1-xx)*(1-yy);
					w2 = (1-xx)*yy;
					w3 = xx*(1-yy);
					w4 = xx*yy;
					for (m=0;m<=3;m++)
					{
						tvols_subperiod[l][k][m] =	w1*tvols_discret[current_index_obsfreq][current_index_statvar][m]+
													w2*tvols_discret[current_index_obsfreq+1][current_index_statvar][m]+
													w3*tvols_discret[current_index_obsfreq][current_index_statvar+1][m]+
													w4*tvols_discret[current_index_obsfreq+1][current_index_statvar+1][m];
					
		
					  if (tvols_subperiod[l][k][m]<0)
						{
							tvols_subperiod[l][k][m]=0;
						}
					}
				}		/* End of loop on state variable */
			}	/* End of loop on the observation days */



			if (! (0==ptr->gear[i+i0]))		/* Floating digital computation */
			{
					start_accrued = ptr->tCpnStart[i+i0];
					end_accrued = add_unit(start_accrued,ptr->undAccrued,SRT_MONTH,MODIFIED_SUCCEEDING);
					G1 = LGMGFromTS(start_accrued,tsPtr);
					df_start_accrued_init = swp_f_df(tNow,start_accrued,ycName);
					df_start_accrued_init *= exp(-0.5*G1*G1*zeta);
					G2 = LGMGFromTS(end_accrued,tsPtr);
					df_end_accrued_init = swp_f_df(tNow,end_accrued,ycName);
					df_end_accrued_init *= exp(-0.5*G2*G2*zeta);
					level_init = coverage(start_accrued,end_accrued,BASIS_ACT_360);

					error = GetVol_SABRValues(tNow,
									  ptr->tEx[j+j0],
									  start_accrued,
									  end_accrued,
									  ptr,
									  &alpha, &beta, &sigma, &rho);
					if (error) return error;
			}


					/* computes probabilities to be between barriers, and corresponding ratios */
			for (k=0; k<=nx; k++)	/* Loop on the state variable */	
			{		
				if (! (0==ptr->gear[i+i0]))		/* Floating digital computation */
				{
					df_start_accrued = df_start_accrued_init*exp(G1*x[k]*reduc[j]);
					df_end_accrued = df_end_accrued_init*exp(G2*x[k]*reduc[j]);
					level = level_init*df_end_accrued;
					fra_accrued = (df_start_accrued - df_end_accrued)/level;  
					tfra_accrued[j][i][k] = fra_accrued;
					error = GetVol_FudgeTimeSwap_new(
								tNow,
								ptr->tEx[j+j0],
								start_accrued,
								end_accrued,
								strikeup,
								fra_accrued,
								ptr,
								&vol_accrued_up,
								alpha, beta, sigma, rho);
					if (error) return error;

					error = GetVol_FudgeTimeSwap_new(
								tNow,
								ptr->tEx[j+j0],
								start_accrued, 
								end_accrued,
								strikelow,
								fra_accrued,
								ptr,
								&vol_accrued_low,
								alpha, beta, sigma, rho);
					if (error) return error; 
				}

				/* Initializes the ratio */
				expected_accrued_ratio = 0;
				expected_accrued_ratio_adj = 0;
				expected_accrued_ratio_low = 0;
				expected_accrued_ratio_high = 0;
				for (l=0;l<observation_freq;l++) 	/* Loop on the number of days */
				{
					if (k<tindex_firstpositive[l])	/* case of negative fra */
					{
						proba_above = 0;
						proba_above_adj = 0;
						if (strikelow==0)
						{
							proba_below = 0;
							proba_below_adj = 0;
						}
						else
						{
							proba_below = 1;
							proba_below_adj = 1;
						}
					}
					else
					{
						if (ptr->gear[i+i0]==0)
						{
							proba_below_adj = 0;
							proba_above_adj = 0;
						}
						else			/* Floating digital */
						{
						/* computes the correl (Libor,index) */
						correl = (tmaturities_subperiod[l] - tmaturities_subperiod[0])/
								(tmaturities_subperiod[observation_freq-1]-tmaturities_subperiod[0])*
								ptr->CorrelEnd + 
								(tmaturities_subperiod[observation_freq-1] - tmaturities_subperiod[l])/(
								tmaturities_subperiod[observation_freq-1]-tmaturities_subperiod[0])*
								ptr->CorrelStart;
						adjustmentup = exp(correl*vol_accrued_up*tvols_subperiod[l][k][1]*
										tmaturities_subperiod[0]); /* Changed the adjustment */
						callup_adj = srt_f_optblksch(tcms_subperiod[l][k]*adjustmentup,strikeup,
							tvols_subperiod[l][k][1],tmaturities_subperiod[l],1,SRT_CALL,PREMIUM);
						callupshift_adj = srt_f_optblksch(tcms_subperiod[l][k]*adjustmentup,strikeupshift,
							tvols_subperiod[l][k][3],tmaturities_subperiod[l],1,SRT_CALL,PREMIUM);
						proba_above_adj = (callup_adj - callupshift_adj)/ptr->call_spread_up;
						if (ptr->barriers[i+i0][0] == 0)
						{
							proba_below_adj = 0;
						}
						else
						{
							adjustmentlow = exp(correl*vol_accrued_low*tvols_subperiod[l][k][0]*
										tmaturities_subperiod[0]); /* Changed the adjustment */
							putlow_adj = srt_f_optblksch(tcms_subperiod[l][k]*adjustmentlow,strikelow,
							tvols_subperiod[l][k][0],tmaturities_subperiod[l],1,SRT_PUT,PREMIUM);
							putlowshift_adj = srt_f_optblksch(tcms_subperiod[l][k]*adjustmentlow,
								strikelowshift,tvols_subperiod[l][k][2],tmaturities_subperiod[l],1,
								SRT_PUT,PREMIUM);
							proba_below_adj = (putlowshift_adj - putlow_adj)/ptr->call_spread_low;
						}
						}
						/* UPPER BARRIER */

						/* call spread */
					/*	counter4 = clock();		*/	/* for debug */
						callup = srt_f_optblksch(tcms_subperiod[l][k],strikeup,tvols_subperiod[l][k][1],tmaturities_subperiod[l],
						1,SRT_CALL,PREMIUM);
						callupshift = srt_f_optblksch(tcms_subperiod[l][k],strikeupshift,tvols_subperiod[l][k][3],
							tmaturities_subperiod[l],1,SRT_CALL,PREMIUM);

					/*	counter4= clock() - counter4;
						counter3 += counter4;		*/

						/* probability above upper barrier */
						proba_above = (callup - callupshift)/ptr->call_spread_up;

						/* LOWER BARRIER */
						if (0==strikelow)
						{
							proba_below = 0;
						}
						else
						{
							/* put spread */
							putlow = srt_f_optblksch(tcms_subperiod[l][k],strikelow,tvols_subperiod[l][k][0],
								tmaturities_subperiod[l],1,SRT_PUT,PREMIUM);
							putlowshift = srt_f_optblksch(tcms_subperiod[l][k],strikelowshift,tvols_subperiod[l][k][2],
								tmaturities_subperiod[l],1,SRT_PUT,PREMIUM);

							/* probability below lower barrier */
							proba_below = (putlowshift - putlow)/ptr->call_spread_low;

						}
					}

					/* computes proba to be between the barriers */
					proba_between = 1 - proba_below - proba_above;
					expected_accrued_ratio += proba_between * (ptr->ratiodays_for_subperiod)[i+i0][l];
					proba_between_adj = 1 - proba_below_adj - proba_above_adj;
					expected_accrued_ratio_adj += proba_between_adj * (ptr->ratiodays_for_subperiod)[i+i0][l];
					expected_accrued_ratio_low += proba_below * (ptr->ratiodays_for_subperiod)[i+i0][l];
//					expected_accrued_ratio_high += (1.0 - proba_above) * (ptr->ratiodays_for_subperiod)[i+i0][l];

				}	/* End of loop on observation days */

				/* computes ratio   */
				expected_ratio[j][i][k] = expected_accrued_ratio;
				expected_ratio_adj[j][i][k] = expected_accrued_ratio_adj;
				expected_ratio_low[j][i][k] = expected_accrued_ratio_low;
				expected_ratio_high[j][i][k] = expected_accrued_ratio;				
			}	/* End of loop on state variable */
		}	/* End of loop on coupon */
	}		/* End of loop on exo date */




	

	

	/* Substract pv of exercise strike */
	mult = 1.0;
	if (ptr->PayRec == SRT_PAYER) mult = -1.0;

	/* Fills the payoff */
	for (j=0; j<nEx; j++)
	{
		zeta = LGMZetaFromTS(ptr->tEx[j+j0], tNow, tsPtr);		/* Get zeta at exercise date */

		ifirst = ptr->iSet[j+j0]-i0;
		for (i=ifirst; i<nCpn; i++)
		{
			Start = ptr->tCpnStart[i+i0];
			End = ptr->tCpnPay[i+i0];
			dfPay_init = swp_f_df(tNow,End,ycName);
			G = LGMGFromTS(End,tsPtr);
			dfPay_init *= exp(-0.5*G*G*zeta);

			for (k=0; k<=nx; k++)
				{
					/* Computes the discount factor until the end of the period */
					dfPay = dfPay_init*exp(G*x[k]*reduc[j]);
					/* Fills in the ratios */
					expected_accrued_ratio = expected_ratio[j][i][k];
					expected_accrued_ratio_adj = expected_ratio_adj[j][i][k];
					expected_accrued_ratio_low = expected_ratio_low[j][i][k];
					expected_accrued_ratio_high = expected_ratio_high[j][i][k];

					/* Formula for the payoff */
					/* Fix Coupon */
					payoff[j][k] += dfPay * ptr->tCvgCpn[i+i0] * ptr->tCpn[i+i0] 
							*expected_accrued_ratio;
			
					/* Floating Digital */
					if (! ptr->gear[i+i0]==0)
						payoff[j][k] += dfPay*ptr->tCvgCpn[i+i0]*ptr->gear[i+i0]*
							tfra_accrued[j][i][k]*expected_accrued_ratio_adj;
					/* Funding payment*/
					payoff[j][k] += dfPay*ptr->tFundingPayment[i+i0];
					if (nCpn-1==i && ptr->nofunding==0)
					{
						payoff[j][k] += dfPay;
					}
				}
		}

		dfPay_init = swp_f_df(tNow,ptr->tSet[j+j0],ycName);
		G = LGMGFromTS(ptr->tSet[j+j0],tsPtr);
		dfPay_init *= exp(-0.5*G*G*zeta);

		for (k=0; k<=nx; k++)
		{
			/* computes discount factor to settlement date */
			dfPay = dfPay_init*exp(G*x[k]*reduc[j]);
			payoff[j][k] -= dfPay*ptr->strike[j]; 
			payoff[j][k] *= mult;
		}
	}


	
	
	
	
	
	

	for (i=0;i<ptr->nCpn;i++)
	{
		for (j=0;j<observation_freq;j++)
		{
			free((observationdf)[i][j]);
		}
		free((observationdf)[i]);
	}
	free(observationdf);

	/* End time counter */
	/*counter2 = clock();
	ptr->reval_times[1] = counter2-counter1;
	ptr->reval_times[2] = counter3;
	ptr->reval_times[3] = counter_getvolfudgetimeswap_new;
	ptr->reval_times[4] =counter_observationdf1;*/
		



	
	/* Free and return */
	
 cleanup: 
    srt_free(zEx); 
	free_f3tensor(expected_ratio,0,nEx-1,0,nCpn-1,0,nx);
	free_f3tensor(expected_ratio_adj,0,nEx-1,0,nCpn-1,0,nx);
	free_f3tensor(expected_ratio_low,0,nEx-1,0,nCpn-1,0,nx);
	free_f3tensor(expected_ratio_high,0,nEx-1,0,nCpn-1,0,nx);
	free_f3tensor(tfra_accrued,0,nEx-1,0,nCpn-1,0,nx);
	free_f3tensor(tvols_discret,0,num_subdiscret_obsfreq-1,0,num_subdiscret_statvar-1,0,3);
	free_dmatrix(tfras_subperiod,0,observation_freq-1,0,nx);
	free_dmatrix(tcms_subperiod,0,observation_freq-1,0,nx);
	free_dvector(tmaturities_subperiod,0,observation_freq-1);
	free_f3tensor(tvols_subperiod,0,observation_freq-1,0,nx,0,3); 
	free_lngvector(tindex_firstpositive,0,observation_freq-1);
	free_lngmatrix(index_discret_statvar,0,observation_freq-1,0,num_subdiscret_statvar-1);
	free_lngvector(index_discret_obsfreq,0,num_subdiscret_obsfreq-1);

	return (error); 
}


