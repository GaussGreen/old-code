/********************************************************************************
  File name: bnd_f_repocurve.c
  Repo Curve Stripper: for GC and Special
********************************************************************************/

#include "srt_h_all.h"
#include "utallhdr.h"
#include "swp_h_all.h"
#include "swp_h_vol.h" 
#include "math.h"


/********************************************************************************
  Input: TradeDate; N term dates and rates in increasing order of dates
  Output: Arrays of consecutive dates and corresponding forward and term rates
********************************************************************************/
Err srt_f_GCRepoCurve
(
	long		TradeDate,
	int			N,
	long*		TermDates,	/* [0 to N-1] */
	double*		TermRates,	/* [0 to N-1] */
	long*		FwdDates,	/* [0 to ...] */
	double*		DailyFwdRates,	/* [0 to ...] */
	double*		DailyTermRates	/* [0 to ...] */
)
{
	Err	err=NULL;
	int i;
	int NumDays;
	double *x, *G, G1i, G1f, *G2;
	double *Gfit;

	// Step 1: Transform variables; call cubic spline function.
	x = dvector(1, N);
	G = dvector(1, N);
	G2 = dvector(1, N);

	for (i=1; i<=N; i++)
	{
		x[i] = TermDates[i-1] - TradeDate;
		G[i] = log(1.0 + TermRates[i-1] * x[i] / 36000.0);	// G = -log(discount)
	}

	G1i = G[1]/x[1];
	G1f = 1.0e30;

	// Solve cubic spline for G2[] -- second derivative
	spline(x, G, N, G1i, G1f, G2);


	// Step 2: Get discount factor and forward overnight and term rates for every calendar day.
	NumDays = (int) x[N] ;
	Gfit = dvector(0, NumDays);

	Gfit[0] = 0;
	for ( i=(int)( x[1] ); i<=NumDays; i++ )
	{
		err = splint(x, G, G2, N, i, &Gfit[i]);
	}

	if (( int)( x[1] ) > 1 )		// Cases when O/N rate applies for more than 1 day -- weekend, etc
	{
		for (i=1; i<(int)( x[1] ); i++)
			Gfit[i] = G[1] * i / x[1];
	}

	for (i=1; i<=NumDays; i++)
	{
		FwdDates[i-1] = TradeDate + i ;
		DailyFwdRates[i-1] = ( exp( Gfit[i] - Gfit[i-1] ) - 1 ) * 36000 ;
		DailyTermRates[i-1] = ( exp( Gfit[i] ) - 1 ) * 36000 / i ;
	}


	// Free memory
	free_dvector(x, 1, N);
	free_dvector(G, 1, N);
	free_dvector(G2, 1, N);
	free_dvector(Gfit, 0, NumDays);
	
	return err;

}


/********************************************************************************
  Input: TradeDate; N term dates and rates in increasing order of dates
  Output: Arrays of consecutive dates and corresponding forward and term rates
********************************************************************************/
Err srt_f_SpecialRepoCurve
(
	long		TradeDate,
	int			N,
	long*		TermDates,	/* [0 to N-1] */
	double*		TermRates,	/* [0 to N-1] */
	long*		FwdDates,	/* [0 to ...] */
	double*		DailyFwdRates,	/* [0 to ...] */
	double*		DailyTermRates	/* [0 to ...] */
)
{
	Err err=NULL;
	int i, j;
	int NumDays;
	long *xn;
	double *disc, rateTemp, discTemp;
	xn = lngvector(1, N);
	disc = dvector(1, N);

	for (i=1; i<=N; i++)
	{
		xn[i] = TermDates[i-1] - TradeDate;
		disc[i] = 1.0 / (1.0 + TermRates[i-1]*xn[i]/36000);
	}

	NumDays = xn[N];

	// Get forward rates
	rateTemp = ( pow(1.0/disc[1], 1.0/xn[1]) - 1 ) * 36000;
	for (i=1; i<=xn[1]; i++)
	{
		FwdDates[i-1] = TradeDate + i ;
		DailyFwdRates[i-1] = rateTemp ;
	}

	for (i=2; i<=N; i++)
	{
		rateTemp = ( pow(disc[i-1]/disc[i], 1.0/(xn[i]-xn[i-1])) - 1 ) * 36000;
		for (j=xn[i-1]+1; j<=xn[i]; j++)
		{
			FwdDates[j-1] = TradeDate + j;
			DailyFwdRates[j-1] = rateTemp;
		}
	}

	// Get term rates
	discTemp = 1.0;
	for (i=1; i<=NumDays; i++)
	{
		discTemp /= 1.0 + DailyFwdRates[i-1]/36000 ;
		DailyTermRates[i-1] = ( 1.0 / discTemp - 1 ) * 36000 / i ;
	}

	// Free memory
	free_lngvector(xn, 1, N);
	free_dvector(disc, 1, N);

	return err;
}