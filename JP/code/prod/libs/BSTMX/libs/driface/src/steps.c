/*===================================================================
* PROGRAM:		STEPS.C
**
** AUTHOR:		Joris Klaentschi.	22 September 1997.
**
** CONTAINS:	Tree step processing function for backward induction 
**				from end to start of a single step with constant
**				variance across step.
// 27 Jul 99	Moved to core library	JK
**=================================================================*/


#include "steps.h"

#include <math.h>

#define DBL_THRESHHOLD (double)1.0e-13	/* good to 12 decimal places */
#define DBL_EQUAL(a,b) (fabs((a) - (b)) < DBL_THRESHHOLD)
#define SQR(x)	((x)*(x))


GTO_EXPORT(int) DrTreesOptGen1DSingleStep(
		double*		vAssetEnd,			/* (I) underlying asset value at the end of the step */
		double*		vPayoffEnd,			/* (I) option payoff on underlying asset values at the end of the step */
		long		nEnd,				/* (I) number of nodes at end of step */
		double*		vAssetStart,		/* (I) underlying asset values at the start of the step */
		long		nStart,				/* (I) number of nodes at start of step */			
		double		stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
		double		hiLimit,			/* (I) if specified, value of high integration limit to use */	
		double		loLimit,			/* (I) if specified, value of low integration limit to use */	
		double*		vPayoffStart)		/* (O) calculated vector payoff values at the start of the step, not allocated here */
/*-----------------------------------------------------------------------------
** CONTAINS:	Implementation of continuous asset value / payoff value tree
**				valuation.  Backward induction from 'end' to 'start'. Parameters
**				and variables named in this convention.
**				Doesn't ZERO the RESULTS vector!!!
**---------------------------------------------------------------------------*/
{	
	/* Optimisation:
	The inner loop here accounts for about 70% of the time taken to value a Flexicap, so is worth optimising.
	It is essentially CPU bound in the two exp() functions - comment them out to see the speed increase.
	Inlining the call to N(x) saves about 50%.
	We get about another 25% by re-arranging the code until the optimiser kicks in - essentially trial and error.
	Don't "tidy up" this section unless you measure the impact on the speed.

	TODO:
	The right way to do this is
	for(i=...)
		for(j=0; j<nEnd;...)
		{
			Nz = exp() ...
			for(nLive=..)
				CORE LOOP HERE
		}

	which should give a big speed increase. Unfortunately, we currently have nEnd[nLive], so the whole tree
	needs re-structuring to acheive this.	JBL 6/Aug/98.

	// Original, unoptimised code follows for reference:
	//for each node in the start slice, add the value of this end slice node - dependent on start slice value 
		for (i=0; i<nStart; ++i)
		{
			// adjustments to interpolating polynomial coefficients 
			b = 2*ax*vAssetStart[i]*stdev + bx*stdev;
			c = cx + bx*vAssetStart[i] + ax*SQR(vAssetStart[i]);

			// integrate payoff function over end slice nodes j and j+1 and add value to summation on the start node i 
			// upper limit 
			z = (upper-vAssetStart[i]) * stddevInv;
			Nz = N(z);
			expBuffer = exp(-0.5*SQR(z)) * sqrt2piInv;
			upValue = c*Nz - b*expBuffer + a*(Nz-z*expBuffer);

			// lower limit 
			z = (lower-vAssetStart[i]) * stddevInv;
			Nz = N(z);
			expBuffer = exp(-0.5*SQR(z)) * sqrt2piInv;
			loValue = c*Nz - b*expBuffer + a*(Nz-z*expBuffer);

			// add to value of start cell 
			vPayoffStart[i] += (upValue-loValue);
		}
	=============================================
	*/

	double		x0, x1, x2;					/* polynomial interpolation - underlying asset values at set points */
	double		f0, f1, f2;					/* polynomial interpolation - corresponding option payoff values at set points  */
	double		ax, bx, cx;					/* 'pure' coefficients for polynomial interpolation function */
	double		a, b, c;					/* 'adjusted' coefficients used in polynomial interpolation function */
	long		i, j;						/* index variables */
	double		upper, lower;				/* upper and lower limits for each integration */
	double		z;							/* integration variable buffers */
	double		upValue, loValue;			/* value of integration at upper/lower limits */
	double		expBuffer, Nz,bxStdev,axStddev;				/* loop precalculations */
	double		k, zz, n,vas_i;
	
	/* for each block in the end slice */
	double sqrt2piInv = 1.0/sqrt(2*PI);
	double stddevInv = 1.0/stdev;

	/* Constants for normal expansion */
	#define N0		0.2316419 
	#define A1		1.330274429                 
	#define A2		1.821255978
	#define A3		1.781477937    
	#define A4		0.356563782          
	#define A5		0.31938153    

	for (j=0; j<=nEnd-3; ++j)
	{
		/* save a bit of calculation time */
		f0 = vPayoffEnd[j];	
		f1 = (nEnd > 1) ? vPayoffEnd[j+1] : f0;	
		f2 = (nEnd > 2) ? vPayoffEnd[j+2] : f1;
		x0 = vAssetEnd[j];	
		x1 = (nEnd > 1) ? vAssetEnd[j+1] : x0;	
		x2 = (nEnd > 2) ? vAssetEnd[j+2] : x1;

		/* integration limits */
		upper = (j >= nEnd-3) ? hiLimit : vAssetEnd[j+1]; 
		lower = (j == 0) ? loLimit : vAssetEnd[j]; 

		/* find the interpolating polynomial coefficients for the payoff function at the end slice - no #DIV0s! */
		ax = DBL_EQUAL(x2,x0) || DBL_EQUAL(x2,x1) || DBL_EQUAL(x1,x0) ? 0.0 : ((f2-f1)/(x2-x1) - (f1-f0)/(x1-x0)) / (x2-x0);
		bx = DBL_EQUAL(x1,x0) ? 0.0 : (f1-f0)/(x1-x0) - ax*(x1+x0);
		cx = f0 - ax*SQR(x0) - bx*x0;
		
		a = ax*SQR(stdev);
		bxStdev = bx*stdev; 
		axStddev = 2.0*ax*stdev;

		/* for each node in the start slice, add the value of this end slice node - dependent on start slice value */
		for (i=0; i<nStart; ++i)
		{
			/* adjustments to interpolating polynomial coefficients */
			vas_i = vAssetStart[i];
			b = axStddev*vas_i + bxStdev; 
			c = cx + bx*vas_i + ax*SQR(vas_i);

			/* integrate payoff function over end slice nodes j and j+1 and add value to summation on the start node i */ 
			/* upper limit: Compute expBuffer=dN(z) and Nz=N(z) in parallel. */
			z = (upper-vas_i) * stddevInv;
			if(fabs(z)<37)	/* Avoid underflow on Alpha */
				expBuffer = exp(-0.5*z*z) * sqrt2piInv;
			else
				expBuffer = 0.0;

			k = 1.0/(1.0 + N0*fabs(z));
			zz = ((((A1*k-A2)*k+A3)*k-A4)*k+A5)*k;
			n = expBuffer*zz;       
			if (z >= 0.0)
				Nz = (1.0-n);
			else
				Nz = n;

			upValue = - b*expBuffer + c*Nz + a*(Nz-z*expBuffer);

			/* lower limit: Compute dN(z) and N(z) in parallel. */
			z = (lower-vas_i) * stddevInv;
			if(fabs(z)<37)	/* Avoid underflow on Alpha */
				expBuffer = exp(-0.5*z*z) * sqrt2piInv;
			else
				expBuffer = 0.0;
			
			k = 1.0/(1.0 + N0*fabs(z));
			zz = ((((A1*k-A2)*k+A3)*k-A4)*k+A5)*k;
			n = expBuffer*zz;       
			if (z >= 0.0)
				Nz = (1.0-n);
			else
				Nz = n;

			loValue = - b*expBuffer + c*Nz  + a*(Nz-z*expBuffer);

			/* add to value of start cell */
			vPayoffStart[i] += (upValue-loValue);
		}


	}
	/* if we got this far, it's good */
	return SUCCESS;
}

/* Constants for normal expansion */
#undef N0	
#undef A1	              
#undef A2	
#undef A3	 
#undef A4	       
#undef A5	

GTO_EXPORT(int) DrTreesSplit1DSingleStep(
		double*		vAssetEnd,			/* (I) underlying asset value at the end of the step */
		double*		vPayoffEnd,			/* (I) option payoff on underlying asset values at the end of the step */
		long		nEnd,				/* (I) number of nodes at end of step */
		long		nEndOne,			/* (I) number of nodes at end in middle slice. */
		long		nEndTwo,			/* (I) number of nodes at end in bottom slice. */
		double*		vAssetStart,		/* (I) underlying asset values at the start of the step */
		long		nStart,				/* (I) number of nodes at start of step */			
		double		stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
		double		hiLimit,			/* (I) if specified, value of high integration limit to use */	
		double		loLimit,			/* (I) if specified, value of low integration limit to use */	
		double      oneLimit,			/* (I) first middle limit to break integral at */
		double      twoLimit,			/* (I) second middle limit to break integral at */
		double*		vPayoffStart)		/* (O) calculated vector payoff values at the start of the step, not allocated here */
/* TO DO INTEGRAL IN . . . ..
// 1 step  : nEnd = nEndOne = nEndTwo, integral [0...nEnd]
// 2 steps : nEndOne = nEndTwo, integral [0...nEndOne] [nEndOne...nEnd]
// 3 steps : integral [0...nEndOne] [nEndOne...nEndTwo] [nEndTwo...nEnd]
*/
{	
	int i;
	/* Calls 1DSingleStep twice on ecah part of the integral. */
	for(i=0;i<nStart;i++)
		vPayoffStart[i]=0.0;

	DrTreesOptGen1DSingleStep(
		vAssetEnd,			/* (I) underlying asset value at the end of the step */
		vPayoffEnd,			/* (I) option payoff on underlying asset values at the end of the step */
		nEndOne,			/* (I) number of nodes at end of step */
		vAssetStart,		/* (I) underlying asset values at the start of the step */
		nStart,				/* (I) number of nodes at start of step */			
		stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
		oneLimit,			/* (I) if specified, value of high integration limit to use */	
		loLimit,			/* (I) if specified, value of low integration limit to use */	
		vPayoffStart);		/* (O) calculated vector payoff values at the start of the step, not allocated here */

	if (nEndTwo-nEndOne > 0)
	{
		DrTreesOptGen1DSingleStep(
			vAssetEnd+nEndOne,	/* (I) underlying asset value at the end of the step */
			vPayoffEnd+nEndOne,	/* (I) option payoff on underlying asset values at the end of the step */
			nEndTwo-nEndOne,	/* (I) number of nodes at end of step */
			vAssetStart,		/* (I) underlying asset values at the start of the step */
			nStart,				/* (I) number of nodes at start of step */			
			stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
			twoLimit,			/* (I) if specified, value of high integration limit to use */	
			oneLimit,			/* (I) if specified, value of low integration limit to use */	
			vPayoffStart);		/* (O) calculated vector payoff values at the start of the step, not allocated here */
	}

	if (nEnd-nEndTwo > 0)
	{
		DrTreesOptGen1DSingleStep(
			vAssetEnd+nEndTwo,	/* (I) underlying asset value at the end of the step */
			vPayoffEnd+nEndTwo,	/* (I) option payoff on underlying asset values at the end of the step */
			nEnd-nEndTwo,		/* (I) number of nodes at end of step */
			vAssetStart,		/* (I) underlying asset values at the start of the step */
			nStart,				/* (I) number of nodes at start of step */			
			stdev,				/* (I) standard deviation across step - NOT same as ChaseTools volatility! */
			hiLimit,			/* (I) if specified, value of high integration limit to use */	
			twoLimit,			/* (I) if specified, value of low integration limit to use */	
			vPayoffStart);		/* (O) calculated vector payoff values at the start of the step, not allocated here */
	}

	return SUCCESS;
}



