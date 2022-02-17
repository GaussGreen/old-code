
//-------------------------------------------------------------------------------
//---- Calibration of LGM2F on two swaptions to price one coupon of a VolBond ---
//-------------------------------------------------------------------------------


#include "srt_h_all.h"
#include "opfnctns.h"
#include "math.h"
#include "swp_h_vol.h"

#define MAX_CPN			600
#define NUM_HERMITE		6
#define	CALPRES			1.0e-08
#define	NITER			10

/*	Static functions */

/*	Just a safer Gaussian */
static double static_lgmsafenorm(
	double		z)
{
	if (z < -10)
	{
		z = -10.0;
	}
	else if (z > 10)
	{
		z = 10.0;
	}

	return norm (z);	
}

/*	In order to help static_lgmystar: PV of coupons when reconstruction is
		df (t, Ti) = df (0, t, Ti) * exp (a[i] + b[i] * Gaussian)
		and Gaussian is y */
static double static_lgmiv(
	double		y,				/*	The y */
	int			ncpn,			/*	Num coupons */
	double		cpn[],			/*	Discounted coupons */
	double		a[],			/*	The ai */
	double		b[])			/*	The bi */
{
	int			i;
	double		val;

	val = 0.0;
	for (i=0; i<ncpn; i++)
	{
		val += cpn[i] * exp (a[i] + b[i] * y);
	}

	return val;
}

#define			YPREC			1.0e-08
/*	Solve y / static_lgmiv (y) = 0 with Newton */
static double static_solvelgmy(
	int			ncpn,			/*	Num coupons */
	double		cpn[],			/*	Discounted coupons */
	double		a[],			/*	The ai */
	double		b[],			/*	The bi */
	int			*dir)			/*	Output, whether iv is increasing in y
										1: decreasing, -1: increasing */
{
	double		y1, y2, iv1, iv2, temp;
	int			it;

	y1 = 0.0;

	iv1 = static_lgmiv(
		y1,
		ncpn, 
		cpn, 
		a,
		b);
	
	y2 = 1.0e-02;

	iv2 = static_lgmiv(
		y2, 
		ncpn, 
		cpn, 
		a,
		b);

	if (iv2 < iv1)
	{
		*dir = 1;
	}
	else
	{
		*dir = -1;
	}

	if (fabs (iv1) < YPREC || fabs (iv2) < YPREC) 
	{

		if (fabs (iv1) < fabs (iv2))
		{
			return y1;
		}
		else
		{
			return y2;
		}
	}

	for (it=1; it<25; it++)
	{
		if (fabs (y2 - y1) < YPREC || fabs (iv2 - iv1) < YPREC) break;
		
		temp = y2;
		y2 -= iv2 * (y2 - y1) / (iv2 - iv1);
		y1 = temp;
		iv1 = iv2;

		iv2 = static_lgmiv(
			y2,
			ncpn, 
			cpn, 
			a,
			b);
		
		if (fabs (iv2) < YPREC) break;
	}

	if (fabs (iv2) > YPREC)
	{
		/*	We should return an error here */
	}

	return y2;
}

/*	Solve y* in the model: 1F case */
static double static_lgmystar1F(
	int			ncpn,			/*	Num coupons */
	double		cpn[],			/*	Discounted coupons */
/*	G at coupon dates and zeta and G at exercise date, in order to define a and b */
	double		cpn_G[],			
	double		ex_zeta, 
	double		ex_sqzeta,
	double		ex_G,
	int			*dir)			/*	Output, whether iv is increasing in y
										1: decreasing, -1: increasing */
{
	int			i;
	double		a[MAX_CPN], b[MAX_CPN];

	/*	Set a, b */
	for (i=0; i<ncpn; i++)
	{
		a[i] = - 0.5 * ex_zeta * (cpn_G[i] - ex_G) * (cpn_G[i] - ex_G);
		b[i] = - (cpn_G[i] - ex_G) * ex_sqzeta;
	}

	/*	Solver */
	return static_solvelgmy (ncpn, cpn, a, b, dir);
}

/*	Precalculate exponential factors in the IV, and decide what variable to use
	in integration, and which one to use in y* solving  */
static void static_lgm2Fcalcexpfact(
	int			ncpn,
	double		cpn[],
	double		cpn_G1[],
	double		cpn_G2[],
	double		ex_zeta1,
	double		sqz1,
	double		ex_zeta2,
	double		sqz2,
	double		ex_zeta12,
	double		c1,
	double		c2,
	double		ex_G1,
	double		ex_G2,
	double		cstfact[],
	double		n1fact[],
	double		n2fact[],
	int			*highest)	/*	Index of the variable to be used for solving 
									- the one with the highest sensitivity */
{
	int			i;
	double		s1 = 0.0, 
				s2 = 0.0;
	double		half_z1 = 0.5 * ex_zeta1, half_z2 = 0.5 * ex_zeta2;
	double		c1sqz2 = c1 * sqz2, c2sqz2 = c2 * sqz2;

	for (i=0; i<ncpn; i++)
	{
		cstfact[i] = half_z1 * (cpn_G1[i] - ex_G1) * (cpn_G1[i] - ex_G1);
		cstfact[i] += half_z2 * (cpn_G2[i] - ex_G2) * (cpn_G2[i] - ex_G2);
		cstfact[i] += ex_zeta12 * (cpn_G1[i] - ex_G1) * (cpn_G2[i] - ex_G2);
		
		n1fact[i] = (cpn_G1[i] - ex_G1) * sqz1 + (cpn_G2[i] - ex_G2) * c1sqz2;
		
		n2fact[i] = (cpn_G2[i] - ex_G2) * c2sqz2; 

		s1 -= cpn[i] * n1fact[i] * exp (-cstfact[i]);
		s2 -= cpn[i] * n2fact[i] * exp (-cstfact[i]);
	}

	if (fabs (s1) > fabs (s2))
	{
		*highest = 1;
	}
	else
	{
		*highest = 2;
	}
}

/*	Solve y* in the model: 2F case */
static double static_lgmystar2F(
	int			ncpn,			/*	Num coupons */
	double		cpn[],			/*	Discounted coupons */
/*	As output from static_lgm2Fcalcexpfact */
	double		cstfact[],
	double		n1fact[],
	double		n2fact[],
	double		highest,
/*	Condition on Intergral variable = n */
	double		n,
	int			*dir)			/*	Output, whether iv is increasing in y
										1: decreasing, -1: increasing */
{
	int			i;

	double		a[MAX_CPN], b[MAX_CPN];

	for (i=0; i<ncpn; i++)
	{
		a[i] = -cstfact[i];

		if (highest == 1)
		{
			a[i] -= n2fact[i] * n;
			b[i] = -n1fact[i];
		}
		else
		{
			a[i] -= n1fact[i] * n;
			b[i] = -n2fact[i];
		}
	}

	/*	Solver */
	return static_solvelgmy (ncpn, cpn, a, b, dir);
}

/*	Calculate zeta1 from sigma, lambda in the 2F model */
static void static_lgmcalczeta1(
int				nsig,								/*	Number of Sigmas */
double			sig_t[],							/*	Sigma Times */
double			sig[],								/*	Sigmas */
/*	Lambda, Alpha, Beta, Rho */
double			lambda,
double			alpha,
double			gamma,
double			rho,
/*	Output */
int				nzeta,
double			zeta_t[],							/*	Zeta times */
double			zeta[])								/*	Output */
{
	int			i, j;
	double		t0;
	double		z, dz;
	double		zeta_at_sig_t[MAX_CPN];

	t0 = 0.0;
	z = 0.0;

	for (i=0; i<nsig; i++)
	{
		dz = sig[i] * sig[i] * (exp (2 * lambda * sig_t[i]) - exp (2 * lambda * t0)) / 2 / lambda;
		t0 = sig_t[i];
		z += dz;
		zeta_at_sig_t[i] = z;
	}

	i = 0;
	for (j=0; j<nzeta; j++)
	{
		while (i < nsig && sig_t[i] < zeta_t[j])
		{
			i++;
		}
		
		if (i == 0)
		{
			t0 = 0.0;
			z = 0.0;
		}
		else
		{
			t0 = sig_t[i-1];
			z = zeta_at_sig_t[i-1];
		}

		dz = sig[i] * sig[i] * (exp (2 * lambda * zeta_t[j]) - exp (2 * lambda * t0)) / 2 / lambda;
		z += dz;
		zeta_at_sig_t[i] = z;
	}
}

/*	Calculate zeta2 and zeta12 from zeta1, lambda, alpha, gamma and rho in the 2F model */
static void static_lgmcalczeta2zeta12(
int				n,									/*	Number of dates */
double			t[],								/*	Times */
double			zeta1[],							/*	Zeta1 */
/*	Lambda, Alpha, Beta, Rho */
double			lambda,
double			alpha,
double			gamma,
double			rho,
/*	Output */
double			zeta2[],
double			zeta12[])
{
	int			i;
	double		t0;
	double		z1, z2, z12, dz1, dz2, dz12;
	double		l1 = lambda, l2 = lambda + gamma;

	t0 = 0.0;
	z1 = z2 = z12 = 0.0;

	for (i=0; i<n; i++)
	{
		dz1 = zeta1[i] - z1;
		dz2 = dz1 * alpha * alpha
			* ((exp (2 * l2 * t[i]) - exp (2 * l2 * t0)) / 2 / l2)
			/ ((exp (2 * l1 * t[i]) - exp (2 * l1 * t0)) / 2 / l1);
		zeta2[i] = z2 + dz2;
		dz12 = dz1 * alpha * rho
			* ((exp ((l1 + l2) * t[i]) - exp ((l1 + l2) * t0)) / (l1 + l2))
			/ ((exp (2 * l1 * t[i]) - exp (2 * l1 * t0)) / 2 / l1);
		zeta12[i] = z12 + dz12;

		t0 = t[i];
		z1 = zeta1[i];
		z2 = zeta2[i];
		z12 = zeta12[i];
	}	
}

/*	Main functions */

/*	Value of European option on a general set of cash flows within LGM 2F */
double x[NUM_HERMITE+1], w[NUM_HERMITE+1];
double lgmopval2FVolBond(
	int			ncpn,								/*	Number of cash-flows */
	double		cpn[],								/*	Discounted Cash-Flows */
	double		cpn_G1[],							/*	G1 at cash-flow dates */
	double		cpn_G2[],							/*	G2 at cash-flow dates */
	double		ex_zeta1,							/*	Z1 at exercise date */
	double		ex_zeta2,							/*	Z2 at exercise date */
	double		ex_zeta12,							/*	Z12 at exercise date */
	double		ex_G1,								/*	G1 at exercise date */
	double		ex_G2)								/*	G2 at exercise date */
{
	double		cstfact[MAX_CPN], n1fact[MAX_CPN], n2fact[MAX_CPN];
	double		sqz1 = sqrt (ex_zeta1), sqz2 = sqrt (ex_zeta2);
	double		c1 = ex_zeta12 / sqz1 / sqz2, c1sqz2 = c1 * sqz2, 
				c2 = sqrt (1.0 - c1 * c1), c2sqz2 = c2 * sqz2;
	double		*X = &(x[1]), *W = &(w[1]);
	double		ystar;
	int			i, j;
	int			highest, dir;
	double		temp, val = 0.0;

	/*	Check for intrinsic */
	if (ex_zeta1 < 1.0e-10 && ex_zeta2 < 1.0e-10)
	{
		for (i=0; i<ncpn; i++)
		{
			val += cpn[i];
		}
		
		return val;
	}

	/*	Precalculate exponential factors in the IV */
	static_lgm2Fcalcexpfact(
		ncpn,
		cpn,
		cpn_G1,
		cpn_G2,
		ex_zeta1,
		sqz1,
		ex_zeta2,
		sqz2,
		ex_zeta12,
		c1,
		c2,
		ex_G1,
		ex_G2,
		cstfact,
		n1fact,
		n2fact,
		&highest);

	/*	Do Hermite integration */
	for (i=0; i<NUM_HERMITE; i++)
	{
		/*	Disregard points of little weight */
		if (fabs (X[i]) > 5.0) continue;

		temp = 0.0;
		for (j=0; j<ncpn; j++)
		{
			ystar = static_lgmystar2F(
				ncpn, 
				cpn,
				cstfact,
				n1fact,
				n2fact,
				highest,
				X[i] - (highest == 1 ? n2fact[j] : n1fact[j]),
				&dir);

			temp += cpn[j] * static_lgmsafenorm(
				dir * (ystar + (highest == 1 ? n1fact[j] : n2fact[j])));
		}

		val += W[i] * temp;
	}

	return val;
}





/*	Setup G function from lambda */
static void static_lgmsetupG(
	double		lambda,
	int			ncpn,								/*	Number of cash-flows */
	double		cpn_time[],							/*	Cash-Flow times */
	double		cpn_G[],							/*	Output: G at cash-flow dates
															G(T) = (1.0 - exp (- lambda * T )) */
	int			nex,								/*	Number of exercise dates */
	double		ex_time[],							/*	Exercise times */
	double		ex_G[])								/*	Output: G at exercise dates */
{
	int			i;

	for (i=0; i<ncpn; i++)
	{
		cpn_G[i] = (1.0 - exp (-lambda * cpn_time[i])) / lambda;
	}

	for (i=0; i<nex; i++)
	{
		ex_G[i] = (1.0 - exp (-lambda * ex_time[i])) / lambda;
	}

}

/*	Setup G2 function (2F) */
static void static_lgmsetupG2(
	double		lambda,
	double		gamma,
	int			ncpn,								/*	Number of cash-flows */
	double		cpn_time[],							/*	Cash-Flow times */
	double		cpn_G2[],							/*	Output: G2 at cash-flow dates
															G2(T) = (1.0 - exp (- lambda2 * T )) */
	int			nex,								/*	Number of exercise dates */
	double		ex_time[],							/*	Exercise times */
	double		ex_G2[])							/*	Output: G2 at exercise dates */
{
	int			i;
	double		lambda2 = lambda + gamma;

	for (i=0; i<ncpn; i++)
	{
		cpn_G2[i] = (1.0 - exp (-lambda2 * cpn_time[i])) / lambda2;
	}

	for (i=0; i<nex; i++)
	{
		ex_G2[i] = (1.0 - exp (-lambda2 * ex_time[i])) / lambda2;
	}
}

/*	Maximum decreasing factor allowed on variance
		ex: the following structure
				1jan2003	1.00%
				1jan2004	0.25%
		will be refused and replaced by the following
				1jan2003	1.00%
				1jan2004	0.50%
		even if the second option is not perfectly matched */
				
#define			MAX_FACT 0.25 /* 0.50 ^ 2 */


/*	Calibrate zeta to one swaption given G: 2F case */
Err lgm2FcalibzetaTo1SwaptionGivenG(
int				ncpn,								/*	Total number of cash-flow dates */
double			cpn_time[],							/*	Cash-Flow times */
double			cpn[],								/*	Discounted cash-flow */
double			cpn_G1[],							/*	G1 at cash-flow dates */
double			cpn_G2[],							/*	G2 at cash-flow dates */
double			ex_time[],							/*	Exercise times */
double			ex_G1[],							/*	G1 at exercise date */
double			ex_G2[],							/*	G2 at exercise date */
double			strike[],							/*	Strikes */
double			mkt_price[],						/*	Market prices */
double			ex_zeta[],							/*	Output: zetas (1) */
/*	Lambda, Alpha, gamma, rho */
double			lambda,
double			alpha,
double			gamma,
double			rho,
double			volInit)
{
	double		s1, s2, s_last, ds;
	int			it;
	double		l1 = lambda, l2 = lambda + gamma;
	double		z1, z2, z12;
	double		q1, q2, q12;
	double		t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
	double		pr1, pr2;
	double		exp_fact1, exp_fact2, exp_fact12;
	double		temp;
	double		minvar, maxvar, d;
	double		quad_var;
	int			niter;
	
	/*	Initial guess: total vol of 85 - 115 bps */
	/*	s = local variance of first factor */
	s1 = 0.0085 / sqrt (1.0 + alpha * alpha + 2.0 * rho * alpha);
	s1 *= s1;
	s2 = 0.0115 / sqrt (1.0 + alpha * alpha + 2.0 * rho * alpha);
	s2 *= s2;

	/*	Initialisation */
	quad_var = 0.0;
	s_last = volInit*volInit;
	ds = s2 / s1;
	zeta1 = zeta2 = zeta12 = 0.0;
	t = 0.0;

	/*	If only one option, no skipping last */

	exp_fact1 = ((exp (2 * l1 * ex_time[0]) - exp (2 * l1 * t)) / 2 / l1);
	exp_fact2 = ((exp (2 * l2 * ex_time[0]) - exp (2 * l2 * t)) / 2 / l2)
				/ ((exp (2 * l1 * ex_time[0]) - exp (2 * l1 * t)) / 2 / l1);
	exp_fact12 = ((exp ((l1 + l2) * ex_time[0]) - exp ((l1 + l2) * t)) / (l1 + l2))
		/ ((exp (2 * l1 * ex_time[0]) - exp (2 * l1 * t)) / 2 / l1);

	/*	Initial guess: last vol */
	s1 = s_last;
	/*	Initial guess 2: last vol * ds */
	s2 = ds * s1;
	/*	First price */
	dz1 = s1 * exp_fact1;
	z1 = zeta1 + dz1;

	dz2 = dz1 * alpha * alpha * exp_fact2;
	z2 = zeta2 + dz2;
		
	dz12 = dz1 * alpha * rho * exp_fact12;
	z12 = zeta12 + dz12;
	
	pr1 = lgmopval2FVolBond(ncpn, 
					cpn,
					cpn_G1,
					cpn_G2,
					z1, 
					z2,
					z12,
					ex_G1[0],
					ex_G2[0]);
			
	if (fabs (mkt_price[0] - pr1) < CALPRES) 
	{
		s_last = s1;
		quad_var += s1 * (ex_time[0] - t);
		ex_zeta[0] = zeta1 = z1;
		zeta2 = z2;
		zeta12 = z12;
		t = ex_time[0];
		goto ReturnLine;
	}

	/*	Second price */	
	dz1 = s2 * exp_fact1;
	q1 = zeta1 + dz1;
		
	dz2 = dz1 * alpha * alpha * exp_fact2;
	q2 = zeta2 + dz2;

	dz12 = dz1 * alpha * rho * exp_fact12;
	q12 = zeta12 + dz12;
	
	pr2 = lgmopval2FVolBond(ncpn, 
					cpn,
					cpn_G1,
					cpn_G2,
					q1, 
					q2,
					q12,
					ex_G1[0],
					ex_G2[0]);
			
	if (fabs (mkt_price[0] - pr2) < CALPRES) 
	{
		s_last = s2;
		quad_var += s2 * (ex_time[0] - t);
		ex_zeta[0] = zeta1 = q1;
		zeta2 = q2;
		zeta12 = q12;
		t = ex_time[0];
		goto ReturnLine;
	}

	niter = 3 * NITER;
	minvar = 1.0e-16;
	maxvar = 1.0;
			
	d = 0.0;

	for (it=1; it<niter; it++)
	{	
		temp = s2;

		/*	Newton iteration */
		s2 -= (pr2 - (mkt_price[0] + d)) * (s2 - s1) / (pr2 - pr1);

		/*	Out of lower bound */
		if (s2 < minvar)
		{
			/*	Calibrate to market price + 1bp */
			d = CALPRES;
			s2 = temp;
			s2 -= (pr2 - (mkt_price[0] + d)) * (s2 - s1) / (pr2 - pr1);
			/*	Still out of bounds */
			if (s2 < minvar)
			{
				s2 = minvar;
			}
		}

		/*	Out of upper bound */
		if (s2 > maxvar)
		{
			/*	Calibrate to market price - 1bp */
			d = -CALPRES;
			s2 = temp;
			s2 -= (pr2 - (mkt_price[0] + d)) * (s2 - s1) / (pr2 - pr1);
			/*	Still out of bounds */
			if (s2 > maxvar)
			{
				s2 = maxvar;
			}
		}

		if (fabs (s2 - temp) < 1.0e-16) break;

		s1 = temp;
		pr1 = pr2;
		
		/*	Reprice with new vol */
		dz1 = s2 * exp_fact1;
		q1 = zeta1 + dz1;
			
		dz2 = dz1 * alpha * alpha * exp_fact2;
		q2 = zeta2 + dz2;

		dz12 = dz1 * alpha * rho * exp_fact12;
		q12 = zeta12 + dz12;
		
		pr2 = lgmopval2FVolBond(ncpn, 
					cpn,
					cpn_G1,
					cpn_G2,
					q1, 
					q2,
					q12,
					ex_G1[0],
					ex_G2[0]);
			
		if (fabs ((mkt_price[0] + d) - pr2) < CALPRES) break;
	}

	/*	If failed, keep last vol */		
	if (fabs ((mkt_price[0] + d) - pr2) > CALPRES)
	{
		smessage ("Failed to calibrate for exercise date %d", 0);
		s2 = s_last;
			
		dz1 = s2 * exp_fact1;
		dz2 = dz1 * alpha * alpha * exp_fact2;
		dz12 = dz1 * alpha * rho * exp_fact12;

		q1 = zeta1 + dz1;
		q2 = zeta2 + dz2;
		q12 = zeta12 + dz12;
	}

	s_last = s2;
	quad_var += s2 * (ex_time[0] - t);
	ex_zeta[0] = zeta1 = q1;
	zeta2 = q2;
	zeta12 = q12;
	t = ex_time[0];

ReturnLine:
	return NULL;
}


//	Calibrate Long Swaption and price Short Swaption given lambda 
Err lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(
int				Longncpn,			//*	Total number of cash-flow dates (Large swaption)
int				Shortncpn,			//*	Total number of cash-flow dates (Short swaption)
double			cpn_time[],			//*	Cash-Flow times 
double			shortcpn[],			//*	short swaption coupons
double			longcpn[],			//*	long swaption coupons
double			ex_time[],			//*	Exercise times
double			LongStrike,			//*	Strike for Long Swaption 
double			LongPrice,			//*	Market price for Long Swaption 
double			ShortStrike,		//*	Strike for Short Swaption 
double			ex_zeta[],			//*	Output: zetas 
double			lambda,				//	Lambda 
//	Alpha, Gamma, Rho 
double			alpha,
double			gamma,
double			rho,
double			volInit,
double			*ShortPrice)					//	Short Swaption price as output 
{
	static double
		cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN], zeta2[MAX_CPN], zeta12[MAX_CPN];
	Err err;
	long nex = 1;
	double ex_lstrike[1];						//	Strike for long swaption
	double ex_lprice[1];					//	Market price for Long Swaption

	ex_lstrike[0] = LongStrike;
	ex_lprice[0] = LongPrice;

	//	Setup G function 
	static_lgmsetupG (lambda, Longncpn, cpn_time, cpn_G, nex, ex_time, ex_G);
	static_lgmsetupG2 (lambda, gamma, Longncpn, cpn_time, cpn_G2, nex, ex_time, ex_G2);

	err = lgm2FcalibzetaTo1SwaptionGivenG(
			Longncpn,
			cpn_time,
			longcpn,
			cpn_G,
			cpn_G2,
			ex_time,
			ex_G,
			ex_G2,
			ex_lstrike,
			ex_lprice,
			ex_zeta,
			lambda,
			alpha,
			gamma,
			rho,
			volInit);

		if (err)
		{
			return err;
		}

		static_lgmcalczeta2zeta12(
			nex,
			ex_time,
			ex_zeta,
			lambda,
			alpha,
			gamma,
			rho,
			zeta2,
			zeta12);

		*ShortPrice = lgmopval2FVolBond(Shortncpn, //	Number of cash-flow dates, including
														//	start and end date 
										shortcpn, //	
										cpn_G,
										cpn_G2,
										ex_zeta[0],
										zeta2[0],
										zeta12[0],
										ex_G[0],
										ex_G2[0]);

	return NULL;
}



Err VolBondCalibration(
	char *yc_name,						//	Name of the yield curve 
	char *vol_curve_name,				//	Name of the market vol curve 
	char *ref_rate_name,					//	Name of the reference rate 
	long startDate,
	long longEndDate,
	long shortEndDate,
	char *swaption_freq,				//	Frequency and basis of swaption 
	char *swaption_basis,
	double alpha,				//	Alpha, Gamma, Rho 
	double gamma,
	double rho,
	double *lambda,						//	Result : Calibrated Lambda and Sigma 
	double *sig)				
{
	int				i;
	long ex_date[1];
	double longstrike, shortstrike;
	double			ex_time[MAX_CPN],
					ex_lfwd[MAX_CPN],
					ex_llvl[MAX_CPN],
					ex_zeta[MAX_CPN];
	int				ex_cpn[MAX_CPN];
	long			cpn_date[MAX_CPN];
	double			cpn_time[MAX_CPN],
					shortcpn[MAX_CPN],
					longcpn[MAX_CPN];
	long			today; 
	double			dfi, dff;
	double			power;
	double			longMktPrice;
	double			shortMktPrice;
	double			vol;
	SrtCurvePtr		yc_ptr;
	Err				err = NULL;
	long nex;
	double ShortModelPrice;
	double shortfwd;
	long compt;
	double AbsoluteError;

	SwapDP shortSwapDP;
	SwapDP longSwapDP;
	long *shortfloat_pay_dates=NULL;
	long shortfloat_nb_pay_dates;
	long *shortfloat_fixing_dates=NULL;
	long *shortfloat_start_dates=NULL;
	long *shortfloat_end_dates=NULL;
	double *shortfloat_cvgs=NULL;
	double *shortfloat_spreads=NULL;
	long shortfloat_nb_dates;

	long *longfloat_pay_dates=NULL;
	long longfloat_nb_pay_dates;
	long *longfloat_fixing_dates=NULL;
	long *longfloat_start_dates=NULL;
	long *longfloat_end_dates=NULL;
	double *longfloat_cvgs=NULL;
	double *longfloat_spreads=NULL;
	long longfloat_nb_dates;

	long *shortfix_pay_dates=NULL;
	long shortfix_nb_pay_dates; 
	long *shortfix_start_dates=NULL;
	long *shortfix_end_dates=NULL;
	double *shortfix_cvgs=NULL;
	long shortfix_nb_dates;

	long *longfix_pay_dates=NULL;
	long longfix_nb_pay_dates; 
	long *longfix_start_dates=NULL;
	long *longfix_end_dates=NULL;
	double *longfix_cvgs=NULL;
	long longfix_nb_dates;

	double shortLVL, longLVL;
	double shortFwdSwapRate, longFwdSwapRate;
	double shortFloatPV, longFloatPV;
	double shortFloatPVWithoutSpread, longFloatPVWithoutSpread;

	double swp_rte, shortStrike, longStrike;

	long shortNcpn, longNcpn;

	double deltaLambda;
	double ShortModelPriceDeltaLambda;
	double ShortModelPriceDerivative;

	double volInit;

	volInit = *sig;
	
	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}
	today = get_today_from_curve (yc_ptr);

	HermiteStandard (x, w, NUM_HERMITE);

	err = swp_f_initSwapDP(startDate, shortEndDate, swaption_freq, swaption_basis, &shortSwapDP);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_initSwapDP(startDate, longEndDate, swaption_freq, swaption_basis, &longSwapDP);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&shortSwapDP, today, ref_rate_name, 
												&shortfloat_pay_dates,
												&shortfloat_nb_pay_dates,
												&shortfloat_fixing_dates,
												&shortfloat_start_dates,
												&shortfloat_end_dates,
												&shortfloat_cvgs,
												&shortfloat_spreads,
												&shortfloat_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&longSwapDP, today, ref_rate_name, 
												&longfloat_pay_dates,
												&longfloat_nb_pay_dates,
												&longfloat_fixing_dates,
												&longfloat_start_dates,
												&longfloat_end_dates,
												&longfloat_cvgs,
												&longfloat_spreads,
												&longfloat_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	
	err = swp_f_make_FixedLegDatesAndCoverages(&shortSwapDP, today, 
											&shortfix_pay_dates, 
											&shortfix_nb_pay_dates, 
											&shortfix_start_dates, 
											&shortfix_end_dates, 
											&shortfix_cvgs, 
											&shortfix_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FixedLegDatesAndCoverages(&longSwapDP, today, 
											&longfix_pay_dates, 
											&longfix_nb_pay_dates, 
											&longfix_start_dates, 
											&longfix_end_dates, 
											&longfix_cvgs, 
											&longfix_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	shortLVL = 0;
	for(i=1;i<shortfix_nb_pay_dates;++i)
	{
		shortLVL += shortfix_cvgs[i-1] * swp_f_df (today, shortfix_pay_dates[i], yc_name);
	}
	longLVL = 0;
	for(i=1;i<longfix_nb_pay_dates;++i)
	{
		longLVL += longfix_cvgs[i-1] * swp_f_df (today, longfix_pay_dates[i], yc_name);
	}

	shortFloatPV = 0;
	shortFloatPVWithoutSpread = 0;
	for(i=0;i<shortfloat_nb_dates;++i)
	{
		shortFloatPV += swp_f_df (today, shortfloat_pay_dates[i+1], yc_name)
			* ( shortfloat_cvgs[i] * shortfloat_spreads[i]
				+ swp_f_df (today, shortfloat_start_dates[i], yc_name) / swp_f_df (today, shortfloat_end_dates[i], yc_name) - 1
			  );
		shortFloatPVWithoutSpread += swp_f_df (today, shortfloat_pay_dates[i+1], yc_name)
			* (swp_f_df (today, shortfloat_start_dates[i], yc_name) / swp_f_df (today, shortfloat_end_dates[i], yc_name) - 1);
	}

	longFloatPV = 0;
	longFloatPVWithoutSpread = 0;
	for(i=0;i<longfloat_nb_dates;++i)
	{
		longFloatPV += swp_f_df (today, longfloat_pay_dates[i+1], yc_name)
			* ( longfloat_cvgs[i] * longfloat_spreads[i]
				+ swp_f_df (today, longfloat_start_dates[i], yc_name) / swp_f_df (today, longfloat_end_dates[i], yc_name) - 1
			  );
		longFloatPVWithoutSpread += swp_f_df (today, longfloat_pay_dates[i+1], yc_name)
			* (swp_f_df (today, longfloat_start_dates[i], yc_name) / swp_f_df (today, longfloat_end_dates[i], yc_name) - 1);
	}

	shortFwdSwapRate = shortFloatPV / shortLVL;
	longFwdSwapRate = longFloatPV / longLVL;

//	shortStrike = shortFwdSwapRate;
//	longStrike = longFwdSwapRate;
	shortStrike = shortFloatPVWithoutSpread / shortLVL;
	longStrike = longFloatPVWithoutSpread / longLVL;

//	1.)	Setup the bond schedule and its coupons 
	cpn_time[0] = (longfloat_pay_dates[0]-today) * YEARS_IN_DAY;
	longcpn[0] = -swp_f_df(today, longfloat_pay_dates[0], yc_name);
	for(i=1;i<longfix_nb_pay_dates;++i)
	{
		cpn_time[i] = (longfix_pay_dates[i] - today) * YEARS_IN_DAY;
		longcpn[i] = longStrike * longfix_cvgs[i-1] * swp_f_df(today, longfix_pay_dates[i], yc_name);
	}
	longNcpn = longfix_nb_pay_dates;
	longcpn[longfix_nb_pay_dates-1] = longcpn[longfix_nb_pay_dates-1] + swp_f_df(today, longfix_pay_dates[longfix_nb_pay_dates -1], yc_name);

	shortcpn[0] = -swp_f_df(today, shortfloat_pay_dates[0], yc_name);
	for(i=1;i<shortfix_nb_pay_dates;++i)
	{
		shortcpn[i] = shortStrike * shortfix_cvgs[i-1] * swp_f_df(today, shortfix_pay_dates[i], yc_name);

	}
	shortNcpn = shortfix_nb_pay_dates;
	shortcpn[shortNcpn-1] = shortcpn[shortNcpn-1] +  swp_f_df(today, shortfix_pay_dates[shortfix_nb_pay_dates-1], yc_name);

//	Exercise 
	ex_date[0] = add_unit (startDate, - 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	ex_time[0] = (ex_date[0] - today) * YEARS_IN_DAY;
	nex = 1;

	//	Underlyings 

	//	Long 

	dff = swp_f_df (today, cpn_date[longfix_nb_pay_dates-1], yc_name);
	dfi = swp_f_df (today, cpn_date[0], yc_name);

	ex_cpn[0] = 0;
	ex_llvl[0] = longLVL;
	ex_lfwd[0] = longFwdSwapRate;
	longstrike = longFloatPVWithoutSpread / longLVL;

	//	Market Price of Long Swaption 
	err = swp_f_vol(vol_curve_name, startDate, longEndDate, ex_lfwd[0], &vol, &power);
	if (err)
	{
		goto FREE_RETURN;
	}
	if (power > 0.5)
	{
		longMktPrice = srt_f_optblksch(
			longFwdSwapRate,
			longFwdSwapRate,
			vol,
			ex_time[0],
			ex_llvl[0],
			SRT_PUT,
			PREMIUM);
		}
	else
	{
		longMktPrice = srt_f_optblknrm(
			longFwdSwapRate,
			longFwdSwapRate,
			vol,
			ex_time[0],
			ex_llvl[0],
			SRT_PUT,
			PREMIUM);
	}


	//	Market Price of Short Swaption 
	dff = swp_f_df (today, cpn_date[shortfix_nb_pay_dates-1], yc_name);
	dfi = swp_f_df (today, cpn_date[0], yc_name);

	shortfwd = shortFwdSwapRate;
	shortstrike = shortFloatPVWithoutSpread / shortLVL;

	err = swp_f_vol(vol_curve_name, startDate, shortEndDate, shortfwd, &vol, &power);
	if (err)
	{
		goto FREE_RETURN;
	}
	if (power > 0.5)
	{
		shortMktPrice = srt_f_optblksch(
			shortfwd,
			shortfwd,
			vol,
			ex_time[0],
			shortLVL,
			SRT_PUT,
			PREMIUM);
		}
	else
	{
		shortMktPrice = srt_f_optblknrm(
			shortfwd,
			shortfwd,
			vol,
			ex_time[0],
			shortLVL,
			SRT_PUT,
			PREMIUM);
	}

	err = swp_f_ForwardRate(
			startDate,
			shortEndDate,
			swaption_freq,
			swaption_basis,
			yc_name,
			ref_rate_name,
			&swp_rte);
	if(err)
	{
		goto FREE_RETURN;
	}	
	

//	*lambda = 0.1;
	deltaLambda = 0.0001;
	AbsoluteError = 10;
	compt = 0;

	if(	(shortEndDate - startDate) / 365.0 < 0.2)
	{
		err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPrice);
	}
	else
	{
		while((AbsoluteError > CALPRES)&&(compt < 20))
		{
			err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPrice);
			if(fabs(ShortModelPrice - shortMktPrice) <= CALPRES)
			{
				break;
			}

			err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda + deltaLambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPriceDeltaLambda);
			if(fabs(ShortModelPriceDeltaLambda - shortMktPrice) <= CALPRES)
			{
				*lambda = *lambda + deltaLambda;
				break;
			}

			ShortModelPriceDerivative = (ShortModelPriceDeltaLambda - ShortModelPrice) / deltaLambda;

			*lambda = max(0.00001, *lambda + (shortMktPrice - ShortModelPrice) / ShortModelPriceDerivative );

			AbsoluteError = fabs(ShortModelPrice - shortMktPrice);

			compt++;
		}
	}

	//	3.)	Transform zeta into sigma 
	(*sig) = sqrt (ex_zeta[0] * 2 * (*lambda) / (exp (2 * (*lambda) * ex_time[0]) - 1.0));


FREE_RETURN:


	if (shortfloat_pay_dates)
	{
		free(shortfloat_pay_dates);
	}
	if (shortfloat_fixing_dates)
	{
		free(shortfloat_fixing_dates);
	}
	if (shortfloat_start_dates)
	{
		free(shortfloat_start_dates);
	}
	if (shortfloat_end_dates)
	{
		free(shortfloat_end_dates);
	}
	if (shortfloat_cvgs)
	{
		free(shortfloat_cvgs);
	}
	if (shortfloat_spreads)
	{
		free(shortfloat_spreads);
	}
	if (longfloat_pay_dates)
	{
		free(longfloat_pay_dates);
	}
	if (longfloat_fixing_dates)
	{
		free(longfloat_fixing_dates);
	}
	if (longfloat_start_dates)
	{
		free(longfloat_start_dates);
	}
	if (longfloat_end_dates)
	{
		free(longfloat_end_dates);
	}
	if (longfloat_cvgs)
	{
		free(longfloat_cvgs);
	}
	if (longfloat_spreads)
	{
		free(longfloat_spreads);
	}

	if (shortfix_pay_dates)
	{
		free(shortfix_pay_dates);
	}
	if (shortfix_start_dates)
	{
		free(shortfix_start_dates);
	}
	if (shortfix_end_dates)
	{
		free(shortfix_end_dates);
	}
	if (shortfix_cvgs)
	{
		free(shortfix_cvgs);
	}

	if (longfix_pay_dates)
	{
		free(longfix_pay_dates);
	}
	if (longfix_start_dates)
	{
		free(longfix_start_dates);
	}
	if (longfix_end_dates)
	{
		free(longfix_end_dates);
	}
	if (longfix_cvgs)
	{
		free(longfix_cvgs);
	}

	if (err)
	{
		return err;
	}

	return NULL;
}


Err VolBondCalibrationOld(
	char *yc_name,						//	Name of the yield curve 
	char *vol_curve_name,				//	Name of the market vol curve 
	char *ref_rate_name,					//	Name of the reference rate 
	long startDate,
	long longEndDate,
	long shortEndDate,
	char *swaption_freq,				//	Frequency and basis of swaption 
	char *swaption_basis,
	double alpha,				//	Alpha, Gamma, Rho 
	double gamma,
	double rho,
	double *lambda,						//	Result : Calibrated Lambda and Sigma 
	double *sig)				
{
	int				i, j, k;
	long ex_date[1];
	double longstrike, shortstrike;
	double			ex_time[MAX_CPN],
					ex_lfwd[MAX_CPN],
					ex_llvl[MAX_CPN],
					ex_zeta[MAX_CPN];
	int				ex_cpn[MAX_CPN];
	long			cpn_date[MAX_CPN];
	double			cpn_time[MAX_CPN],
					shortcpn[MAX_CPN],
					longcpn[MAX_CPN];
	long			today; 
	double			dfi, dff;
	double			power;
	double			longMktPrice;
	double			shortMktPrice;
	double			vol;
	SrtCurvePtr		yc_ptr;
	Err				err = NULL;
	long nex;
	double ShortModelPrice;
	double shortfwd;
	long compt;
	double AbsoluteError;

	SwapDP shortSwapDP;
	SwapDP longSwapDP;
	long *shortfloat_pay_dates=NULL;
	long shortfloat_nb_pay_dates;
	long *shortfloat_fixing_dates=NULL;
	long *shortfloat_start_dates=NULL;
	long *shortfloat_end_dates=NULL;
	double *shortfloat_cvgs=NULL;
	double *shortfloat_spreads=NULL;
	long shortfloat_nb_dates;

	long *longfloat_pay_dates=NULL;
	long longfloat_nb_pay_dates;
	long *longfloat_fixing_dates=NULL;
	long *longfloat_start_dates=NULL;
	long *longfloat_end_dates=NULL;
	double *longfloat_cvgs=NULL;
	double *longfloat_spreads=NULL;
	long longfloat_nb_dates;

	long *shortfix_pay_dates=NULL;
	long shortfix_nb_pay_dates; 
	long *shortfix_start_dates=NULL;
	long *shortfix_end_dates=NULL;
	double *shortfix_cvgs=NULL;
	long shortfix_nb_dates;

	long *longfix_pay_dates=NULL;
	long longfix_nb_pay_dates; 
	long *longfix_start_dates=NULL;
	long *longfix_end_dates=NULL;
	double *longfix_cvgs=NULL;
	long longfix_nb_dates;

	double shortLVL, longLVL;
	double shortFwdSwapRate, longFwdSwapRate;
	double shortFloatPV, longFloatPV;
	double shortFloatPVWithoutSpread, longFloatPVWithoutSpread;

	double swp_rte, shortStrike, longStrike;

	long shortNcpn, longNcpn;

	double doubletemp;

	double deltaLambda;
	double ShortModelPriceDeltaLambda;
	double ShortModelPriceDerivative;

	double volInit;

	volInit = *sig;
	
	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}
	today = get_today_from_curve (yc_ptr);

	HermiteStandard (x, w, NUM_HERMITE);

	err = swp_f_initSwapDP(startDate, shortEndDate, swaption_freq, swaption_basis, &shortSwapDP);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_initSwapDP(startDate, longEndDate, swaption_freq, swaption_basis, &longSwapDP);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&shortSwapDP, today, ref_rate_name, 
												&shortfloat_pay_dates,
												&shortfloat_nb_pay_dates,
												&shortfloat_fixing_dates,
												&shortfloat_start_dates,
												&shortfloat_end_dates,
												&shortfloat_cvgs,
												&shortfloat_spreads,
												&shortfloat_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&longSwapDP, today, ref_rate_name, 
												&longfloat_pay_dates,
												&longfloat_nb_pay_dates,
												&longfloat_fixing_dates,
												&longfloat_start_dates,
												&longfloat_end_dates,
												&longfloat_cvgs,
												&longfloat_spreads,
												&longfloat_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	
	err = swp_f_make_FixedLegDatesAndCoverages(&shortSwapDP, today, 
											&shortfix_pay_dates, 
											&shortfix_nb_pay_dates, 
											&shortfix_start_dates, 
											&shortfix_end_dates, 
											&shortfix_cvgs, 
											&shortfix_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	err = swp_f_make_FixedLegDatesAndCoverages(&longSwapDP, today, 
											&longfix_pay_dates, 
											&longfix_nb_pay_dates, 
											&longfix_start_dates, 
											&longfix_end_dates, 
											&longfix_cvgs, 
											&longfix_nb_dates);
	if(err)
	{
		goto FREE_RETURN;
	}

	shortLVL = 0;
	for(i=1;i<shortfix_nb_pay_dates;++i)
	{
		shortLVL += shortfix_cvgs[i-1] * swp_f_df (today, shortfix_pay_dates[i], yc_name);
	}
	longLVL = 0;
	for(i=1;i<longfix_nb_pay_dates;++i)
	{
		longLVL += longfix_cvgs[i-1] * swp_f_df (today, longfix_pay_dates[i], yc_name);
	}

	shortFloatPV = 0;
	shortFloatPVWithoutSpread = 0;
	for(i=0;i<shortfloat_nb_dates;++i)
	{
		shortFloatPV += swp_f_df (today, shortfloat_pay_dates[i+1], yc_name)
			* ( shortfloat_cvgs[i] * shortfloat_spreads[i]
				+ swp_f_df (today, shortfloat_start_dates[i], yc_name) / swp_f_df (today, shortfloat_end_dates[i], yc_name) - 1
			  );
		shortFloatPVWithoutSpread += swp_f_df (today, shortfloat_pay_dates[i+1], yc_name)
			* (swp_f_df (today, shortfloat_start_dates[i], yc_name) / swp_f_df (today, shortfloat_end_dates[i], yc_name) - 1);
	}

	longFloatPV = 0;
	longFloatPVWithoutSpread = 0;
	for(i=0;i<longfloat_nb_dates;++i)
	{
		longFloatPV += swp_f_df (today, longfloat_pay_dates[i+1], yc_name)
			* ( longfloat_cvgs[i] * longfloat_spreads[i]
				+ swp_f_df (today, longfloat_start_dates[i], yc_name) / swp_f_df (today, longfloat_end_dates[i], yc_name) - 1
			  );
		longFloatPVWithoutSpread += swp_f_df (today, longfloat_pay_dates[i+1], yc_name)
			* (swp_f_df (today, longfloat_start_dates[i], yc_name) / swp_f_df (today, longfloat_end_dates[i], yc_name) - 1);
	}

	shortFwdSwapRate = shortFloatPV / shortLVL;
	longFwdSwapRate = longFloatPV / longLVL;

	shortStrike = shortFwdSwapRate;
	longStrike = longFwdSwapRate;
//	shortStrike = shortFloatPVWithoutSpread / shortLVL;
//	longStrike = longFloatPVWithoutSpread / longLVL;


//	1.)	Setup the bond schedule and its coupons 
	i=1;
	j=1;
	k=1;
	cpn_time[0] = (longfix_pay_dates[0]-today) * YEARS_IN_DAY;
	longcpn[0] = -swp_f_df(today, longfix_pay_dates[0], yc_name);
	while(i + j < longfix_nb_pay_dates + longfloat_nb_pay_dates)
	{
		if(longfix_pay_dates[i] < longfloat_pay_dates[j])
		{
			cpn_time[k] = (longfix_pay_dates[i] - today) * YEARS_IN_DAY;
			longcpn[k] = longStrike * longfix_cvgs[i-1] * swp_f_df(today, longfix_pay_dates[i], yc_name);
			i++;
		}
		else
		{
			if(longfix_pay_dates[i] > longfloat_pay_dates[j])
			{
				cpn_time[k] = (longfloat_pay_dates[j] - today) * YEARS_IN_DAY;
				longcpn[k] = - longfloat_spreads[j-1] * longfloat_cvgs[j-1] * swp_f_df(today, longfloat_pay_dates[j], yc_name);
				j++;
			}
			else
			{
				cpn_time[k] = (longfix_pay_dates[i] - today) * YEARS_IN_DAY;
				longcpn[k] = (longStrike * longfix_cvgs[i-1] - longfloat_spreads[j-1] * longfloat_cvgs[j-1]) * swp_f_df(today, longfix_pay_dates[i], yc_name);
				i++;
				j++;
			}
		}
		k++;
	}
	longNcpn = k;
	longcpn[longNcpn-1] = longcpn[longNcpn-1] + swp_f_df(today, longfix_pay_dates[longfix_nb_pay_dates -1], yc_name);

	i=1;
	j=1;
	k=1;
	shortcpn[0] = -swp_f_df(today, shortfix_pay_dates[0], yc_name);
	while(i + j < shortfix_nb_pay_dates + shortfloat_nb_pay_dates)
	{
		if(shortfix_pay_dates[i] < shortfloat_pay_dates[j])
		{
			shortcpn[k] = shortStrike * shortfix_cvgs[i-1] * swp_f_df(today, shortfix_pay_dates[i], yc_name);
			i++;
		}
		else
		{
			if(shortfix_pay_dates[i] > shortfloat_pay_dates[j])
			{
				shortcpn[k] = - shortfloat_spreads[j-1] * shortfloat_cvgs[j-1] * swp_f_df(today, shortfloat_pay_dates[j], yc_name);
				j++;
			}
			else
			{
				shortcpn[k] = (shortStrike * shortfix_cvgs[i-1] - shortfloat_spreads[j-1] * shortfloat_cvgs[j-1]) * swp_f_df(today, longfix_pay_dates[i], yc_name);
				i++;
				j++;
			}
		}
		k++;
	}
	shortNcpn = k;
	shortcpn[shortNcpn-1] = shortcpn[shortNcpn-1] +  swp_f_df(today, shortfix_pay_dates[shortfix_nb_pay_dates-1], yc_name);

	for(k=0;k<longNcpn;++k)
	{
		doubletemp = longcpn[k];
		doubletemp = cpn_time[k];
	}

//	Exercise 
	ex_date[0] = add_unit (startDate, - 2, SRT_BDAY, MODIFIED_SUCCEEDING);
	ex_time[0] = (ex_date[0] - today) * YEARS_IN_DAY;
	nex = 1;

	//	Underlyings 

	//	Long 

	dff = swp_f_df (today, cpn_date[longfix_nb_pay_dates-1], yc_name);
	dfi = swp_f_df (today, cpn_date[0], yc_name);

	ex_cpn[0] = 0;
	ex_llvl[0] = longLVL;
	ex_lfwd[0] = longFwdSwapRate;
	longstrike = longFloatPVWithoutSpread / longLVL;

	//	Market Price of Long Swaption 
	err = swp_f_vol(vol_curve_name, startDate, longEndDate, ex_lfwd[0], &vol, &power);
	if (err)
	{
		goto FREE_RETURN;
	}
	if (power > 0.5)
	{
		longMktPrice = srt_f_optblksch(
			longFwdSwapRate,
			longFwdSwapRate,
			vol,
			ex_time[0],
			ex_llvl[0],
			SRT_PUT,
			PREMIUM);
		}
	else
	{
		longMktPrice = srt_f_optblknrm(
			longFwdSwapRate,
			longFwdSwapRate,
			vol,
			ex_time[0],
			ex_llvl[0],
			SRT_PUT,
			PREMIUM);
	}


	//	Market Price of Short Swaption 
	dff = swp_f_df (today, cpn_date[shortfix_nb_pay_dates-1], yc_name);
	dfi = swp_f_df (today, cpn_date[0], yc_name);

	shortfwd = shortFwdSwapRate;
	shortstrike = shortFloatPVWithoutSpread / shortLVL;

	err = swp_f_vol(vol_curve_name, startDate, shortEndDate, shortfwd, &vol, &power);
	if (err)
	{
		goto FREE_RETURN;
	}
	if (power > 0.5)
	{
		shortMktPrice = srt_f_optblksch(
			shortfwd,
			shortfwd,
			vol,
			ex_time[0],
			shortLVL,
			SRT_PUT,
			PREMIUM);
		}
	else
	{
		shortMktPrice = srt_f_optblknrm(
			shortfwd,
			shortfwd,
			vol,
			ex_time[0],
			shortLVL,
			SRT_PUT,
			PREMIUM);
	}

	err = swp_f_ForwardRate(
			startDate,
			shortEndDate,
			swaption_freq,
			swaption_basis,
			yc_name,
			ref_rate_name,
			&swp_rte);
	if(err)
	{
		goto FREE_RETURN;
	}	
	

//	*lambda = 0.1;
	deltaLambda = 0.0001;
	AbsoluteError = 10;
	compt = 0;


	if(	(shortEndDate - startDate) / 365.0 < 0.2)
	{
		err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPrice);
	}
	else
	{
		while((AbsoluteError > CALPRES)&&(compt < 20))
		{
			err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPrice);
			if(fabs(ShortModelPrice - shortMktPrice) <= CALPRES)
			{
				break;
			}

			err = lgm2FPriceShortSwaptionAndCalibrateLongSwaptionGivenLambda(longNcpn,shortNcpn,cpn_time,shortcpn,longcpn,ex_time,longstrike,longMktPrice,shortstrike,ex_zeta,
											*lambda + deltaLambda,
											alpha,
											gamma,
											rho,
											*sig,
											&ShortModelPriceDeltaLambda);
			if(fabs(ShortModelPriceDeltaLambda - shortMktPrice) <= CALPRES)
			{
				break;
			}

			ShortModelPriceDerivative = (ShortModelPriceDeltaLambda - ShortModelPrice) / deltaLambda;

			*lambda = max(0.001, *lambda + (shortMktPrice - ShortModelPrice) / ShortModelPriceDerivative );

			AbsoluteError = fabs(ShortModelPrice - shortMktPrice);

			compt++;
		}
	}


	//	3.)	Transform zeta into sigma 
	(*sig) = sqrt (ex_zeta[0] * 2 * (*lambda) / (exp (2 * (*lambda) * ex_time[0]) - 1.0));


FREE_RETURN:


	if (shortfloat_pay_dates)
	{
		free(shortfloat_pay_dates);
	}
	if (shortfloat_fixing_dates)
	{
		free(shortfloat_fixing_dates);
	}
	if (shortfloat_start_dates)
	{
		free(shortfloat_start_dates);
	}
	if (shortfloat_end_dates)
	{
		free(shortfloat_end_dates);
	}
	if (shortfloat_cvgs)
	{
		free(shortfloat_cvgs);
	}
	if (shortfloat_spreads)
	{
		free(shortfloat_spreads);
	}
	if (longfloat_pay_dates)
	{
		free(longfloat_pay_dates);
	}
	if (longfloat_fixing_dates)
	{
		free(longfloat_fixing_dates);
	}
	if (longfloat_start_dates)
	{
		free(longfloat_start_dates);
	}
	if (longfloat_end_dates)
	{
		free(longfloat_end_dates);
	}
	if (longfloat_cvgs)
	{
		free(longfloat_cvgs);
	}
	if (longfloat_spreads)
	{
		free(longfloat_spreads);
	}

	if (shortfix_pay_dates)
	{
		free(shortfix_pay_dates);
	}
	if (shortfix_start_dates)
	{
		free(shortfix_start_dates);
	}
	if (shortfix_end_dates)
	{
		free(shortfix_end_dates);
	}
	if (shortfix_cvgs)
	{
		free(shortfix_cvgs);
	}

	if (longfix_pay_dates)
	{
		free(longfix_pay_dates);
	}
	if (longfix_start_dates)
	{
		free(longfix_start_dates);
	}
	if (longfix_end_dates)
	{
		free(longfix_end_dates);
	}
	if (longfix_cvgs)
	{
		free(longfix_cvgs);
	}

	if (err)
	{
		return err;
	}

	return NULL;
}


