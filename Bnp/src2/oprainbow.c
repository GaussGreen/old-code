/*******************************************************************************
**
**	OPRAINBOW.C
**
**	pays the minimum of (K, nX * X, nY * Y)
**
**	X and Y are lognormal
**
*******************************************************************************/


/* ==========================================================================  
   include files
   ========================================================================== */


#include		<srt_h_all.h"
#include		<math.h"

double 	srt_f_optrainbow(	double	T,		/* maturity			*/
							double	df,		/* discount factor	*/
							double	K,
							double	nX,
							double	X0,
							double	sigX,
							double	nY,
							double	Y0,
							double	sigY,
							double	rho)
{
double	d1, d2, d3, d4;
double	rho1, rho2;
double	sigXY;
double	res;

	sigXY = sqrt(sigX * sigX + sigY * sigY - 2 * rho * sigX * sigY);
	d1 = (log(nX * X0 / K) - 0.5 * sigX * sigX * T) / (sigX * sqrt(T));
	d2 = (log(nY * Y0 / K) - 0.5 * sigY * sigY * T) / (sigY * sqrt(T));
	d3 = (log(nY * Y0 / nX / X0) - 0.5 * sigXY * sigXY * T) / (sigXY * sqrt(T));
	d4 = (log(nX * X0 / nY / Y0) - 0.5 * sigXY * sigXY * T) / (sigXY * sqrt(T));
	rho1 = (sigX - rho * sigY) / sigXY;
	rho2 = (sigY - rho * sigX) / sigXY;


	res = K * bivar(d1, d2, rho);
	res += nX * X0 * bivar(-d1 - sigX * sqrt(T), d3, rho1);
	res += nY * Y0 * bivar(-d2 - sigY * sqrt(T), d4, rho2);
	res *= df;

	return res;
}

