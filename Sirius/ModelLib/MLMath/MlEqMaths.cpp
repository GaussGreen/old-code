#include "stdafx.h"

#include "MlEqMaths.h"
#include "cMatrix.h"
#include <math.h>

#include "nrutil.h"
#include "ran2.h"
#include "MlEqNormal.h"
#include "utility.h"
#include "Rnd.h"
#include "hpsort.h"
#include <set>


static const double EPS = 1.0e-07;
static const double BRENTEPS = 3.0e-12;

const int MlEqMaths::MATH_ROUND_NONE =0;      
const int MlEqMaths::MATH_ROUND_DOWN =1;      
const int MlEqMaths::MATH_ROUND_UP = 2;       
const int MlEqMaths::MATH_ROUND_TO_ZERO = 3;  
const int MlEqMaths::MATH_ROUND_FROM_ZERO = 4;
const int MlEqMaths::MATH_ROUND_NEAREST = 5;




double MlEqMaths::dsgn(double x) { if (x < 0.0) return -1.0; if (x > 0.0) return 1.0; return 0.0; }
double MlEqMaths::dsqr(double x) { return x * x; }
double MlEqMaths::dcube(double x) { return x * x * x; }
double MlEqMaths::dmax(double x, double y) { if (x > y) return x; return y; }
double MlEqMaths::dmin(double x, double y) { if (x < y) return x; return y; }
int    MlEqMaths::isgn(double x) { if (x < 0.0) return -1; if (x > 0.0) return 1; return 0; }
int    MlEqMaths::imax(int x, int y) { if (x > y) return x; return y; }
int    MlEqMaths::imin(int x, int y) { if (x < y) return x; return y; }

double MlEqMaths::MSIGN(double a,double b)
{
	if( b >= 0.0 ){return fabs(a);}
	else{
		return -1.0*fabs(a);
	}
}

/*#ifndef __MSIGN
#define __MSIGN
#define MSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
*/


double MlEqMaths::Min(double x,double y)
{
	if ( x < y )
		return x;
	else
		return y;
}

double MlEqMaths::Max(double x,double y)
{
	if ( x > y )
		return x;
	else
		return y;
}


void MlEqMaths::choldc(CMatrix& a,CVector& p)
{
	int i,j,k;
	double sum;

	int n = a.rows();
	for ( i = 0 ; i < n; i++ )
	{
		for ( j = i ; j < n; j++ )
		{
			for ( sum = a[i][j], k = i-1; k >= 0; k-- )
				sum -= a[i][k]*a[j][k];

			if ( i == j )
			{
				if ( sum <= 0.0-1e-30 )
				{
					throw("Cholesky decomposition failed");
				}
				if ( sum > -1e-30 && sum <= 0 ){
					sum = 1e-30;
				}

				p[i] = sqrt(sum);
			}
			else
			{
				a[j][i] = sum/p[i];
			}
		}
	}
}


double MlEqMaths::binomial(int n,int k)
{
	double res =  floor(0.5+exp(MlEqMaths::factln(n)-MlEqMaths::factln(k)-MlEqMaths::factln(n-k)));
	return res;
}

double MlEqMaths::factln(int n)
{
	static double a[101];
	if (n < 0){
		throw("Negative factorial in routine factln");
	}
	if (n <= 1){
		return 0.0;
	}
	if (n <= 100) {
		return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	}
	else
	{
		double res =  gammln(n+1.0);
		return res;
	}
}

double MlEqMaths::gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}




double MlEqMaths::Hermite(int n,double x,double y)
{
	if ( n < 0 ){
		throw("hermite index must be greater or equal to zero");
	}
	if ( n == 0 ){
		return 1.0;
	}
	else if ( n == 1 ){
		return x;
	}

	double prevVal		= Hermite(n-1,x,y);
	double prevPrevVal	= Hermite(n-2,x,y);
	double z = Hermite(n,x,y,prevVal,prevPrevVal);
	return z;
}


double MlEqMaths::Hermite(int n,double x,double y,double prevVal,double prevPrevVal)
{
	if ( n < 0 ){
		throw("hermite index must be greater or equal to zero");
	}

	if ( n == 0 ){
		return 1.0;
	}
	else if ( n == 1 ){
		return x;
	}

	double z = x*prevVal-(double)(n-1)*y*prevPrevVal;
	return z;
}


////////////////////////////////////////////////


double MlEqMaths::linearInterp(const CVector& x, const CVector& y,double xx,int method,bool& isEdge)
{
	isEdge = false;

	if ( method )
		throw("only method = 0 implemented at present" );

	int i;
	double res;
	int n = x.getsize();

	if ( n == 1 )
		return y[0];

	locate(x,xx,i);

	if ( i == -1 )
	{
			isEdge = true;
			res = y[0];
			return res;
	}
	else if ( i == n-1 )
	{
			isEdge = true;
			res = y[n-1];
			return res;
	}
	else
	{
		isEdge = false;
		res = y[i] + ( xx - x[i] )/( x[i+1] - x[i] )* ( y[i+1]-y[i] );
		return res;
	}


}

double MlEqMaths::linearInterp(const CVector& x, const CVector& y,double xx,int method)
{
	bool isEdge = false;
	return linearInterp(x, y,xx,method,isEdge);
};


double MlEqMaths::constantInterp(const CVector& x, const CVector& y,double xx)
{
	bool isEdge=0;
	return MlEqMaths::constantInterp(x, y,xx,isEdge);
};

double MlEqMaths::constantInterp(const CVector& x, const CVector& y,double xx,bool & isEdge)
{

	int i;
	double res;
	int n = x.getsize();

	if ( n == 1 )
		return y[0];

	locate(x,xx,i);

	if ( i == -1 )
	{
			res = y[0];
			isEdge = true;
			return res;
	}
	else if ( i == n-1 )
	{
			isEdge = true;
			res = y[n-1];
			return res;
	}
	else
	{
		isEdge = false;
		res = y[i];
		return res;
	}
};

/*void MlEqMaths::Hpsort(CVector& x,GVector< long >& map)
{
	map.resize(x.getsize());
	hpsort(x-1,map-1);
	for ( int i = 0 ; i < map.getsize(); i++ ){
		map[i]--;
	}
}*/

void MlEqMaths::Hpsort(CVector& x)
{
	hpsort(x.getsize(),x-1);
}

void MlEqMaths::locate(const CVector& xx,double x,int& i) // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
{
	locate(xx,x,xx.getsize(),i);
}

void MlEqMaths::locate(const CVector& xx,double x,int n,int& i) // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
{
//	take only fist n out of array

	int ju,jm,jl;
	int ascnd;

	if ( n <= 0 )
		throw("vector size smaller equal zero in interpolation" );


	jl=0;
	ju=n+1;
	ascnd=(xx[n-1]>xx[0]);
	while (ju-jl > 1)
	{
		jm=(ju+jl)>>1;
		if(x>xx[jm-1]==ascnd)
			jl=jm;
		else
			ju=jm;
	}

	i=jl-1;

}


void MlEqMaths::locate(const GVector<long>& xx,long x,int& i) // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
{

	int ju,jm,jl;
	int ascnd;
    int n=xx.getsize();

	if ( n <= 0 )
		throw("vector size smaller equal zero in interpolation" );


	jl=0;
	ju=n+1;
	ascnd=(xx[n-1]>xx[0]);
	while (ju-jl > 1)
	{
		jm=(ju+jl)>>1;
		if(x>xx[jm-1]==ascnd)
			jl=jm;
		else
			ju=jm;
	}

	i=jl-1;

}

void MlEqMaths::Locate(const GVector<long>& xx,long x,int& i) // searches array xx to find i such that xx[ii] <= x <= xx[i+1]	
{
	MlEqMaths::locate(xx,x,i);
	int n = xx.getsize();

/*	if ( i > n ){
		throw("error in locate index");
	}
	else if ( i == n-1 ){
		return;
	}
	else
	{
		i++;
	}
*/
	if ( i == n ){
		i--;
	}else if ( i == -1 ){
		i++;}

	if ( i < n-1 && xx[i+1] == x ){
		i++;}


}


double MlEqMaths::linearInterp(int ndim,CVector& x,CVector& y,double xx,int method)
{

	assert(false);
	return 0.0;

}



double MlEqMaths::dround(double x, int how)
{
	/* Generic double-to-double rounding */
	switch (how)
	{
		case MlEqMaths::MATH_ROUND_DOWN:
			x = floor(x);
			break;
		case MlEqMaths::MATH_ROUND_UP:
			x = ceil(x);
			break;
		case MlEqMaths::MATH_ROUND_TO_ZERO:
			if (x >= 0)
				x = floor(x);
			else
				x = ceil(x);
			break;
		case MlEqMaths::MATH_ROUND_FROM_ZERO:
			if (x >= 0)
				x = ceil(x);
			else
				x = floor(x);
			break;
		case MlEqMaths::MATH_ROUND_NEAREST:
			if (x >= 0)
				x = floor(x + 0.5);
			else
				x = ceil(x - 0.5);
			break;
	}
	return x;
}

int MlEqMaths::iround(double x, int how)
{
	/* Generic double-to-integer rounding */
	x = MlEqMaths::dround(x, how);
	if (x < (double)INT_MIN)
		return INT_MIN;
	if (x > (double)INT_MAX)
		return INT_MAX;
	return (int)x;
}


int MlEqMaths::is_integer(double x)
{
	/* To be integer has to be representable. */
	if ((x < (double)INT_MIN) || (x > (double)INT_MAX))
		return 0;
	if (x == MlEqMaths::dround(x, MlEqMaths::MATH_ROUND_NEAREST))
		return 1;
	return 0;
} /* is_integer */


double MlEqMaths::Min(const double& a, const double& b)
{
	if ( a <= b )
		return a;
	else 
		return b;
}

double MlEqMaths::Max(const double& a, const double& b)
{
	if ( a >= b )
		return a; 
	else
		return b;
}
int MlEqMaths::Min(const int& a, const int& b)
{
	if ( a <= b )
		return a;
	else 
		return b;
}

int MlEqMaths::Max(const int& a, const int& b)
{
	if ( a >= b )
		return a; 
	else
		return b;
}

bool MlEqMaths::equal(const double& a, const double& b)
{
	if ( fabs( a - b ) < EPS )
		return true; 
	else 
		return false; 
}

bool MlEqMaths::deqz(double x)
{
//    double DBL_EPSILON   =  2.2204460492503131e-016; /* smallest such that 1.0+DBL_EPSILON != 1.0 */

	return	(fabs(x) < 2.2204460492503131e-016);
}

double MlEqMaths::math_factorial(int x, int y)
{
	double factorial;

//	assert(x >= 0 && y >= 0);

	if (x == y)
		return 1.0;

	factorial = 1;
	if (y < x)
	{
		int index;

		for (index = y + 1; index <= x; index++)
			factorial *= (double)index;
	}
	else
	{
		int index;

		for (index = x + 1; index <= y; index++)
			factorial /= (double)index;
	}

	return factorial;
}



int MlEqMaths::math_cubic_roots(double a, double b, double c, double d, double* roots)
{
	double Q_term;
	double R_term;
	double R2_Q3_term;

//	assert(fabs(a) > 2.2204460492503131e-016);
	b /= a;
	c /= a;
	d /= a;

	Q_term = (b * b - 3.0 * c) / 9.0;
	R_term = (2.0 * b * b * b - 9.0 * b * c + 27.0 * d) / 54.0;
	R2_Q3_term = R_term * R_term - Q_term * Q_term * Q_term;

	if (R2_Q3_term < 0.0)
	{
		const double tmp = -2.0 * sqrt(Q_term);
		double       theta_term;
		double       x1;
		double       x2;
		double       x3;

		/* There are three real roots (but two maybe equal) */
		theta_term = R_term / pow(Q_term, 3.0 / 2.0);
//		assert((theta_term >= -1.0) && (theta_term <= 1.0));
		theta_term = acos(theta_term);

		x1 = tmp * cos(theta_term / 3.0) - b / 3.0;
		x2 = tmp * cos((theta_term + 2 * 3.141592653589793238462643) / 3.0) - b / 3.0;
		x3 = tmp * cos((theta_term - 2 * 3.141592653589793238462643) / 3.0) - b / 3.0;
		if ((x1 <= x2) && (x1 <= x3))
		{
			roots[0] = x1;
			if (x2 <= x3)
			{
				roots[1] = x2;
				roots[2] = x3;
			}
			else
			{
				roots[1] = x3;
				roots[2] = x2;
			}
		}
		else if (x2 <= x3)
		{
			roots[0] = x2;
			if (x1 <= x3)
			{
				roots[1] = x1;
				roots[2] = x3;
			}
			else
			{
				roots[1] = x3;
				roots[2] = x1;
			}
		}
		else
		{
			roots[0] = x3;
			if (x1 <= x2)
			{
				roots[1] = x1;
				roots[2] = x2;
			}
			else
			{
				roots[1] = x2;
				roots[2] = x1;
			}
		}
		// return number of roots found
		return 3;
	}
	else
	{
		double A_term;
		double B_term;
		double x1;

		if (R_term < 0.0)
			A_term = pow(fabs(R_term) + sqrt(R2_Q3_term), 1.0 / 3.0);
		else
			A_term = - pow(fabs(R_term) + sqrt(R2_Q3_term), 1.0 / 3.0);

		if (fabs(A_term) <= 2.2204460492503131e-016)
			B_term = 0.0;
		else
			B_term = Q_term / A_term;

		x1 = (A_term + B_term) - b / 3.0;

		roots[0] = x1;
		// return number of roots found
		return 1;
	}
	// return number of roots found
	return 0;
}


// root finding stuff
////////////////////////////////////////////////////////////////////////////////////////////



typedef struct _rootfind_nonzero_t
{
	bool (*fn)(double x,void* vp,double* f);
	int (*dfndx)(double x,void* vp,double* f);
	void* vp;
	double diff_y;
}
rootfind_nonzero_t;


inline int deqz(double x) {if (fabs(x) < EPS) return 1;return 0 ;};

int rootfind_bracket(bool (*fn)(double x,void* vp,double* f), double* lower_x, double* upper_x, double factor,
				int max_tries, double find_y, void* vp)
{
	double f1, f2;
	bool retval;
	int tries;

	/* Bracket a root between two x-ordinates */
	if (*lower_x >= *upper_x)
	{
		throw(  "Input error: Lower limits must be smaller than upper limits" );
	}
	retval = (*fn)(*lower_x, vp, &f1);
	if (retval != true)  
		return retval;
	f1 -= find_y;
	retval = (*fn)(*upper_x, vp, &f2);
	if (retval != true)
		return retval;
	f2 -= find_y;
	for (tries = 2; tries < max_tries; tries++)
	{
		if ((f1 * f2) <= 0.0)
			return true;
		if (fabs(f1) < fabs(f2))
		{
			/* Move the lower guess */
			*lower_x -= factor * (*upper_x - *lower_x);
			retval = (*fn)(*lower_x, vp, &f1);
			if (retval != true)
				return retval;
			f1 -= find_y;
		}
		else
		{
			/* Move the upper guess */
			*upper_x += factor * (*upper_x - *lower_x);
			retval = (*fn)(*upper_x, vp, &f2);
			if (retval != true)
				return retval;
			f2 -= find_y;
		}
	}

	throw( "Root not found" );
}



/****************************************************************************/

static bool rootfind_bisection(bool (*fn)(double x,void* vp,double* f), double x1, const double* y1, double x2, const double* y2, 
					  double xacc, int max_tries, void* vp, double* found_x)
{
	double dx, f, fmid, xmid, rtb;
	bool retval;
	int tries = 0;

	/* Use the bisection method to find the root */
	if (y1 != NULL)
		f = *y1;
	else
	{
		tries++;
		retval = (*fn)(x1, vp, &f);
		if (retval != true)
			return retval;
	}
	if (y2 != NULL)
		fmid = *y2;
	else
	{
		tries++;
		retval = (*fn)(x2, vp, &fmid);
		if (retval != true)
			return retval;
	}
	//assert((f * fmid) <= 0);
	if (f < 0.0)
	{
		/* Slopes up */
		dx = x2 - x1;
		rtb = x1;
	}
	else
	{
		/* Slopes down */
		dx = x1 - x2;
		rtb = x2;
	}
	for (; tries <= max_tries; tries++)
	{
		/* Bisect the interval and replace interval limit which has same sign as */
		/* function evaluated at the bisection with the bisected point.          */
		dx *= 0.5;
		xmid = rtb + dx;
		retval = (*fn)(xmid, vp, &fmid);
		if (retval != true)
			return retval;
		if (fmid <= 0.0)
			rtb = xmid;
		if ((fabs(dx) < xacc) || (fmid == 0.0))
		{
			/* Solution found */
			*found_x = rtb;
			return true;
		}
	}

	throw( "Root not found" );
}

/****************************************************************************/


static bool rootfind_secant(bool (*fn)(double x,void* vp,double* f), double x1, const double* y1, double x2, const double* y2, 
                           double xacc, int max_tries, void* vp, double* found_x)
{
	double fl, f, dx, xl, rts;
	bool retval;
	int tries = 0;

	/* Ensure that |fl| < |f| by swapping the ends of the interval. */
	if (y1 != NULL)
		fl = *y1;
	else
	{
		tries++;
		retval = (*fn)(x1, vp, &fl);
		if (retval != true)
			return retval;
	}
	if (y2 != NULL)
		f = *y2;
	else
	{
		tries++;
		retval = (*fn)(x2, vp, &f);
		if (retval != true)
			return retval;
	}
	if (fabs(fl) < fabs(f))
	{
		/* Need to swap to get the smaller value */
		double swap = fl;
		fl = f;
		f = swap;
		xl = x2;
		rts = x1;
	}
	else
	{
		xl = x1;
		rts = x2;
	}
	for (; tries <= max_tries; tries++)
	{
		/* Using (f - fl)/(rts - xl) as an estimate of the derivative of fn, perform an */
		/* approximate iteration of the Newton-Rhapson method.  On the next iteration   */
		/* (xl,fl) is replaced by (rts,f) which is replaced by the updated solution.    */
		dx = f - fl;
		if (deqz(dx))
			break;
		dx = (xl - rts) * f / (f - fl);
		xl = rts;
		fl = f;
		rts += dx;
		retval = (*fn)(rts, vp, &f);
		if(retval != true)
			return retval;
		if ((fabs(dx) < xacc) || (f == 0.0))
		{
			*found_x = rts;
			return true;
		}
	}
	*found_x = rts;
	throw( "Root not found" );
}

/****************************************************************************/



static bool rootfind_false_position(bool (*fn)(double x,void* vp,double* f), 
					       double x1, const double* y1, double x2, const double* y2, 
					       double xacc, int max_tries, void* vp, double* found_x)
{
	double fl, fh, xl, xh, dx, del, f, rtf;
	bool retval;
	int tries = 0;

	/* Swap the ends of the interval so that the solution slopes up from (xl,fl) to  */
	/* (xh, fh).                                                                     */
	if (y1 != NULL)
		fl = *y1;
	else
	{
		tries++;
		retval = (*fn)(x1, vp, &fl);
		if (retval != true)
			return retval;
	}
	if (y2 != NULL)
		fh = *y2;
	else
	{
		tries++;
		retval = (*fn)(x2, vp, &fh);
		if (retval != true)
			return retval;
	}
	//assert((fl * fh) <= 0.0);
	if (fl < fh)
	{
		xl = x1;
		xh = x2;
	}
	else
	{
		double swap = fl;
		fl = fh;
		fh = swap;
		xl = x2;
		xh = x1;
	}
	rtf = 0.5 * (xl + xh);
	for (; tries <= max_tries; tries++)
	{
		/* Use secant to bisect the interval and then replace the interval limit that */
		/* has the same sign as the function evaluated at the bisection.              */
		dx = xh - xl;
		rtf = xl + dx * fl / (fl - fh);
		retval = (*fn)(rtf, vp, &f);
		if (retval != true)
			return retval;
		if (f < 0.0)
		{
			/* Move the lower limit */
			del = xl - rtf;
			xl = rtf;
			fl = f;
		}
		else
		{
			/* Move the upper limit */
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		if ((fabs(del) < xacc) || (f == 0.0))
		{
			/* Found root */
			*found_x = rtf;
			return true;
		}
	}
	*found_x = rtf;
	throw(  "Root not found" );
}


/****************************************************************************/

inline double rootfind_sign(double a, double b)
{
	/* Return 'a' but with the sign of 'b' */
	if (b >= 0)
		return fabs(a);
	return -fabs(a);
}

/****************************************************************************/

static bool rootfind_ridders(bool (*fn)(double x,void* vp,double* f),
					double x1, const double* y1, double x2, const double* y2, 
					double xacc, int max_tries, void* vp, double* found_x)
{
	double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;
	bool retval;
	int tries = 0;

	/* Use Ridders' method to find the root */
	if (y1 != NULL)
		fl = *y1;
	else
	{
		tries++;
		retval = (*fn)(x1, vp, &fl);
		if (retval != true)
			return retval;
	}
	if (fl == 0.0)
	{
		/* Root is on the lower boundary */
		*found_x = x1;
		return true;
	}
	if (y2 != NULL)
		fh = *y2;
	else
	{
		tries++;
		retval = (*fn)(x2, vp, &fh);
		if (retval != true)
			return retval;
	}
	if (fh == 0.0)
	{
		/* Root is on the upper boundary */
		*found_x = x2;
		return true;
	}
	//assert((fl * fh) < 0.0);

	/* Main interval reduction loop */
	xl = x1;
	xh = x2;
	ans = 0.5 * (xl + xh);
	for (; tries <= max_tries; tries++)
	{
		xm = 0.5 * (xl + xh);
		retval = (*fn)(xm, vp, &fm);
		if (retval != true)
			return retval;
		if (fm == 0.0)
		{
			/* Lucky strike */
			*found_x = xm;
			return true;
		}
		s = fm * fm - fl * fh;
		if (s <= 0.0)
		{
			/* Function is not well-behaved */
			*found_x = ans;
			throw(  "Root not found" );
		}
		if (fl >= fh)
			xnew = xm + (xm - xl) * fm / sqrt(s);
		else
			xnew = xm - (xm - xl) * fm / sqrt(s);
		if ((tries > 1) && (fabs(xnew - ans) < xacc))
		{
			/* Didn't get very far during this attempt */
			*found_x = xnew;
			return true;
		}
		ans = xnew;
		retval = (*fn)(ans, vp, &fnew);
		if (retval != true)
			return retval;
		if (fnew == 0.0)
		{
			/* Lucky strike again */
			*found_x = ans;
			return true;
		}
		if (rootfind_sign(fm, fnew) != fm)
		{
			/* Keep the root bracketed */
			xl = xm;
			fl = fm;
			xh = ans;
			fh = fnew;
		}
		else if (rootfind_sign(fl, fnew) != fl)
		{
			/* Keep the root bracketed */
			xh = ans;
			fh = fnew;
		}
		else if (rootfind_sign(fh, fnew) != fh)
		{
			/* Keep the root bracketed */
			xl = ans;
			fl = fnew;
		}
		else
			throw(  "unreachable code reached");
		if (fabs(xh - xl) < xacc)
		{
			/* That's close enough */
			*found_x = ans;
			return true;
		}
	}
	*found_x = ans;
	throw(  "root not found");
}





/****************************************************************************/

static bool rootfind_brent(bool (*fn)(double x,void* vp,double* f), double x1, const double* y1, double x2, const double* y2, 
                          double xacc, int max_tries, void* vp, double* found_x)
{
	double a, b, c, d, e, min1, min2, fa, fb, fc, p, q, r, s, tol1, xm;
	bool retval;
	int tries = 0;

	/* Use Brent's method to find the root */
	a = x1;
	b = x2;
	c = x2;
	d = b - a;
	e = d;
	if (y1 != NULL)
		fa = *y1;
	else
	{
		tries++;
		retval = (*fn)(a, vp, &fa);
		if (retval != true)
			return retval;
	}
	if (y2 != NULL)
		fb = *y2;
	else
	{
		tries++;
		retval = (*fn)(b, vp, &fb);
		if (retval != true)
			return retval;
	}
	fc = fb;
	//assert((fa * fb) <= 0);

	/* Main interval reduction loop */
	for (; tries <= max_tries; tries++)
	{
		/* Adjust the bounding interval */
		if ((fb * fc) > 0.0)
		{
			c = a;
			fc = fa;
			d = b - a;
			e = d;
		}
		if (fabs(fc) < fabs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}

		/* Convergence check */
		tol1 = 2.0 * BRENTEPS * fabs(b) + 0.5 * xacc;
		xm = 0.5 * (c - b);
		if ((fabs(xm) <= tol1) || (fb == 0.0))
		{
			/* Found root */
			*found_x = b;
			return true;
		}
		if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb)))
		{
			/* Attempt inverse quadratic interpolation */
			s = fb / fa;
			if (a == c)
			{
				/* Optimised for repeated function values */
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}
			else
			{
				/* Full interpolation */
				q = fa / fc;
				r = fb / fc;
				p = s * (2.0 * xm * q * (q-r) - (b-a) * (r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)
				q = -q;
			else
				p = -p;
			min1 = 1.5 * xm * q - fabs(tol1 * q);
			min2 = 0.5 * fabs(e * q);
			if ((p < min1) && (p < min2))
			{
				/* Accept the interpolation */
				e = d;
				d = p / q;
			}
			else
			{
				/* Interpolation failed, so use bisection */
				d = xm;
				e = d;
			}
		}
		else
		{
			/* Bounds decreasing too slowly, so use bisection */
			d = xm;
			e = d;
		}
		a = b;
		fa = fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += rootfind_sign(tol1, xm);
		retval = (*fn)(b, vp, &fb);
		if ( fabs(fb) < xacc )
		{
			*found_x = b;
			return true;
			// found case where algorithm had to be forced to stop here
		}
		if (retval != true)
			return retval;
	}
	*found_x = b;
	throw(  "root not found");

}

/****************************************************************************/


#define NTRY	50


bool zbrac(bool (*func)(double,void* vp,double *f), void* vp,
			  double * x1, 
			  double * x2,
			  double * pf1,
			  double * pf2,
			  double   prec=0,
			  double   factor=1.6)
{
	/*
		Bracketed root function [Num. Rec. C 9.1].
		Input:	double (*)(double)
				double *
				double *
				double *
				double *
				double 
				double
		Output:	int
		Description: Given a function func and an initial guessed range x1 to x2, the routine
			expands the range geometrically until a root is bracketed by the returned values
			x1 and x2 (in which case zbrac returns true) or until the range becomes unacceptably
			large (in which case zbrac returns false).
			If the function evaluation is within prec, set *x1 and *x2 to the value tested and
			*f1 and *f2 to the resulting function values.
	*/

	double f1, f2;

	if ( *x1 == *x2 )
	{
		throw( "Bad initial range in zbrac");

		return true;
	}

	bool retval;

	retval = (*func)(*x1,vp,&f1);
	if (retval != true)
		return retval;

	retval = (*func)(*x2,vp,&f2);
	if (retval != true)
		return retval;

	for (int j=1; j<=NTRY; j++)
	{
		if ( fabs(f1) < prec )
		{
			*x2 = *x1;
			if (pf1 != 0) *pf1 = f1;
			if (pf2 != 0) *pf2 = f1;
			return true;
		}

		if ( fabs(f2) < prec )
		{
			*x1 = *x2;
			if (pf1 != 0) *pf1 = f2;
			if (pf2 != 0) *pf2 = f2;
			return true;
		}

		if ( f1 * f2 < 0.0 )
		{
			if (pf1 != 0) *pf1 = f1;
			if (pf2 != 0) *pf2 = f2;
			return true;
		}

		if ( fabs(f1) < fabs(f2) )
			(*func)(*x1 += factor * (*x1 - *x2), vp, &f1);
		else
			(*func)(*x2 += factor * (*x2 - *x1), vp, &f2);
	}

	if (pf1 != 0) *pf1 = f1;
	if (pf2 != 0) *pf2 = f2;

	return false;
}

#undef NTRY

void zbrak(bool (*fx)(double,void* vp,double *f),void *vp,double x1,double x2,int n,double xb1[],double xb2[],int *nb )
{

	int nbb,i;
	double x,fp,fc,dx;

	nbb = 0;
	dx = (x2-x1)/n;
	
	(*fx)(x=x1,vp,&fp);
	for ( i = 1;i<=n;i++ )
	{
		(*fx)(x+= dx,vp,&fc);
		if ( fc*fp < 0.0 ){
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb == nbb)return;
		}

		fp = fc;
	}
	*nb=nbb;
}

bool rootfind_brent_grow_range(	bool (*fn)(double,void* vp, double* f),
								double x1, 
								const double* y1, 
								double x2, 
								const double* y2,
								double xacc, 
								int max_tries, 
								void* vp, 
								double* found_x, 
								int look_outside_user_bounds)
{
	double fx1, fx2;
	bool retval = true;

	if ( y1 == NULL ) {
		retval = (*fn)(x1, vp, &fx1);
		if (retval != true)
			return retval;
	}
	else
		fx1 = *y1;

	if ( y2 == NULL ) {
		retval = (*fn)(x2, vp, &fx2);
		if (retval != true)
			return retval;
	}
	else
		fx2 = *y2;

	int nb;

	if ( fx1*fx2 > 0.0 )
	{

	//  reset limits start searching for better limits inside

//		CVector xb1(1);
//		CVector xb2(1);

		int n = 64;//32;//8;//4;sos
		nb = 1;

		CVector xb1(n+1);
		CVector xb2(n+1);

		zbrak(fn,vp,x1,x2,n,xb1.getPtr(),xb2.getPtr(),&nb);

		if ( nb > 0 )// new bounds found
		{
			x1 = xb1[1];
			x2 = xb2[1];

		}
		else if ( look_outside_user_bounds )
		{

		// search outside

			double prec = 0.0;
			double factor = 1.2;//sos1.6;


			retval = zbrac(fn, vp,
							&x1, 
							&x2,
							&fx1,
							&fx2,
							prec,
							factor);

			if (retval!=true)
				return retval;

		}
	}

	if ( (fx1*fx2 > 0) && ( nb <= 0 ) )
	{

		throw(  "Bad initial range in zbrac");
	}


	retval = (*fn)(x1, vp, &fx1);
	if (retval != true)
		return retval;
	
	(*fn)(x2, vp, &fx2);
	if (retval != true)
		return retval;

	return rootfind_brent(fn, x1, &fx1, x2, &fx2, xacc, max_tries, vp, found_x);
	
}





static double rootfind_cblfr_guess(double x1, double x2, double x3, double d1, double d2, double d3)
{
	static const double lower_interp = 0.1;
	static const double upper_interp = 0.9;
	double x;

	/* We have a zero denominator for 'CBLFR' so use the secant method in three points */
	if ((d1 * d2) < 0.0)
	{
		/* Guesses 1 and 2 bracket the solution */
		if ((d1 * d3) > 0.0)
		{
			/* So do guesses 2 and 3 */
			if ((d1 * d2) < (d2 * d3))
			{
				/* Interpolate between x2 and x3 */
				x = d3 / (d3 - d2);
				if (x < lower_interp)
					x = lower_interp;
				else if (x > upper_interp)
					x = upper_interp;
				return x3 + x * (x2 - x3);
			}
			else
			{
				/* Interpolate between x1 and x2 */
				x = d1 / (d1 - d2);
				if (x < lower_interp)
					x = lower_interp;
				else if (x > upper_interp)
					x = upper_interp;
				return x1 + x * (x2 - x1);
			}
		}
		else
		{
			/* So do guesses 1 and 3 */
			if ((d1 * d3) < (d1 * d2))
			{
				/* Interpolate between x1 and x2 */
				x = d1 / (d1 - d2);
				if (x < lower_interp)
					x = lower_interp;
				else if (x > upper_interp)
					x = upper_interp;
				return x1 + x * (x2 - x1);
			}
			else
			{
				/* Interpolate between x1 and x3 */
				x = d1 / (d1 - d3);
				if (x < lower_interp)
					x = lower_interp;
				else if (x > upper_interp)
					x = upper_interp;
				return x1 + x * (x3 - x1);
			}
		}
	}
	else if ((d1 * d3) < 0.0)
	{
		/* Guesses 1 and 3 bracket the solution */
		if ((d1 * d3) < (d2 * d3))
		{
			/* Interpolate between x2 and x3 */
			x = d3 / (d3 - d2);
			if (x < lower_interp)
				x = lower_interp;
			else if (x > upper_interp)
				x = upper_interp;
			return x3 + x * (x2 - x3);
		}
		else
		{
			/* Interpolate between x3 and x2 */
			x = d1 / (d1 - d2);
			if (x < lower_interp)
				x = lower_interp;
			else if (x > upper_interp)
				x = upper_interp;
			return x1 + x * (x1 - x3);
		}
	}

	/* All three delta have the same sign, so use a very crude extrapolation */
	if (fabs(d1) < fabs(d2))
		x = x1;
	else
		x = x2;
	if (fabs(d1) < fabs(d3))
		x += x1;
	else
		x += x3;
	if (fabs(d2) < fabs(d3))
		x += x2;
	else
		x += x3;
	x -= x1 + x2 + x3;
	return x / 1.5;
}

/****************************************************************************/

static bool rootfind_cblfr(bool (*fn)(double x,void* vp,double* f), double x1, const double* y1, double x2, const double* y2, 
                          const double* p_x3, const double* y3, double xacc, int max_tries, void* vp, double* found_x)
{
	double x3, d1, d2, d3, guess;
	bool retval;
	int tries = 0;

	/* If the user gives all three then we don't need to estimate the third one. */
	if (p_x3 != NULL)
	{
		x3 = *p_x3;

		if (y3 != NULL)
			d3 = *y3;
		else
		{
			tries++;
			retval = fn(x3, vp, &d3);
			if (retval != true)
				return retval;
		}
		if (fabs(d3) < xacc)
		{
			/* That's a spot of luck */
			*found_x = x3;
			return true;
		}

		if (y2 != NULL)
			d2 = *y2;
		else
		{
			tries++;
			retval = (*fn)(x2, vp, &d2);
			if (retval != true)
				return retval;
		}
		if (fabs(d2) < xacc)
		{
			/* That's almost as lucky */
			*found_x = x2;
			return true;
		}

		if (y1 != NULL)
			d1 = *y1;
		else
		{
			tries++;
			retval = (*fn)(x1, vp, &d1);
			if (retval != true)
				return retval;
		}
		if (fabs(d1) < xacc)
		{
			/* That's taking the biscuit */
			*found_x = x1;
			return true;
		}
	}
	else
	{
		/* Evaluate at x1 and x2. */
		if (y1 != NULL)
			d1 = *y1;
		else
		{
			tries++;
			retval = (*fn)(x1, vp, &d1);
			if (retval != true)
				return retval;
		}
		if (fabs(d1) < xacc)
		{
			/* That's a spot of luck */
			*found_x = x1;
			return true;
		}

		if (y2 != NULL)
			d2 = *y2;
		else
		{
			tries++;
			retval = (*fn)(x2, vp, &d2);
			if (retval != true)
				return retval;
		}
		if (fabs(d2) < xacc)
		{
			/* That's a spot of luck */
			*found_x = x2;
			return true;
		}

		/* Get x3 by interpolation or bisection. */
		if ( deqz(MlEqMaths::dsqr(d1 - d2)))
			x3 = (x1 + x2) / 2.0;
		else
			x3 = x1 - (x2 - x1) * MlEqMaths::Min(5.0, (double) (fabs(d1 / (d2 - d1))) ) * MlEqMaths::dsgn(d1 / (d2 - d1));
		retval = (*fn)(x3, vp, &d3);
		if (retval != true)
			return retval;
		if (fabs(d3) < xacc)
		{
			/* That's taking the biscuit */
			*found_x = x3;
			return 0;
		}
	}

	/* Use the 'CBLFR' algorithm to find the root */
	//assert(x2 != x3);
	for (; tries <= max_tries; tries++)
	{
		/* Fetch the next best guess given three previous guesses and deltas */
		guess = d2 * d3 * (x2 - x3) + d1 * d2 * (x1 - x2) + d3 * d1 * (x3 - x1);
		if (deqz(guess))
			guess = rootfind_cblfr_guess(x1, x2, x3, d1, d2, d3);
		else
			guess = x3 + (d3 * (d1 - d2) * (x3 - x1) * (x2 - x3)) / guess;

		/* Shuffle and test */
		x1 = x2;
		x2 = x3;
		x3 = guess;
		d1 = d2;
		d2 = d3;
		retval = (*fn)(x3, vp, &d3);
		if (retval != true)
			return retval;
		if (fabs(d3) < xacc)
		{
			/* That's my boy */
			*found_x = x3;
			return true;
		}
	}
	*found_x = x3;
	throw(  "root not found");
}

/****************************************************************************/

static bool rootfind_newton(bool (*fn)(double x,void* vp,double* f), int (*dfndx)(double x,void* vp,double* f) , double x1, const double* y1, 
                           const double* dydx1, double xacc, int max_tries, void* vp, double* found_x)
{
	double deriv;
	double function;
	double x = x1;
	double oldx;
	double h;
	double oldh;
	int    tries;
	int    count;
	bool    retval;
	int    n_zero_deriv;

    h = 0.0;
	oldx = x1;
	oldh = 0.0;
	n_zero_deriv = 0;
	tries = 0;
	count = 0;
	while (count <= max_tries)
	{
		tries++;

		/* Evaluate the derivative. */
		if ((tries == 1) && (dydx1 != NULL))
			deriv = *dydx1;
		else
		{
			int retval = (*dfndx)(x, vp, &deriv);
			if (retval != 1)
				return false;
		}
		if (deqz(deriv))
		{
			n_zero_deriv++;
			if (n_zero_deriv > 5)
			{
				*found_x = 0.0;
				throw( "root not found");
			}

			/* Derivative is zero. Half the latest adjustment and try the iteration again. */
			if (deqz(h))
			{
				/* Last adjustment was zero, i.e. first iteration, bump the initial guess. */
				x += fabs(x1) / 1000.0;
			}
			else
			{
				/* Halve the adjustment and try again. */
				h /= 2.0;
				oldx += h;
			}

			continue;
		}
		n_zero_deriv = 0;

		/* Evaluate the function. */
		if ((tries == 1) && (y1 != NULL))
			function = *y1;
		else
		{
			count++;
			retval = (*fn)(x, vp, &function);
			if (retval != true)
				return retval;
		}

		/* Calculate the adjustment. */
		h = - function / deriv;

		/* Test for convergence. */
		if (fabs(h) < xacc)
		{
			*found_x = x + h;
			return true;
		}

		/* Are we oscillating around the root? */
		if (fabs(h + oldh) < 1.0e-8)
		{
			/* Root lies in the middle somewhere. */
			h /= 2.0;
		}

		/* Store the old iteration and update. */
		oldx = x;
		oldh = h;
		x += h;
	}

	*found_x = x1;
	throw(  "root not found");

} /* rootfind_newton */


/****************************************************************************/

static bool rootfind_nonzero (double x, void* vp, double* f)
{
	rootfind_nonzero_t* nonzero = (rootfind_nonzero_t*)vp;
	bool retval;

	/* Callback to use when the required function is not to be solved for zero */
	retval = nonzero->fn(x, nonzero->vp, f);
	if (retval != true)
		return retval;
	*f -= nonzero->diff_y;
	return true;
}

/****************************************************************************/

static int rootfind_nonzero_deriv (double x, void* vp, double* dfdx)
{
	rootfind_nonzero_t* nonzero = (rootfind_nonzero_t*)vp;

	/* Callback to use when the required function is not to be solved for zero */
	return nonzero->dfndx(x, nonzero->vp, dfdx);
}


int rootfind_solve(int rootfind_flag, bool (*fn)(double x,void* vp,double* f), double lower_x, double upper_x,
				   double accuracy, int max_tries, double find_y, void* vp, double* found_x,
				   int (*dfndx)(double x,void* vp,double* f) )
{
	static const double bracket_factor = 1.618;
	static const int bracket_tries = 50;
	rootfind_nonzero_t nonzero;
	int tries;


	/* Make sure the root is bracketed, if necessary */
	if ((rootfind_flag != eROOTFIND_SECANT) && (rootfind_flag != eROOTFIND_CBLFR)&&
	    (rootfind_flag != eROOTFIND_NEWTON) && (rootfind_flag != eROOTFIND_BRENT_GROWRANGE) &&
		(rootfind_flag != eROOTFIND_BRENT_NOGROW))
	{
		/* Only the secant and CBLFR methods can get away without bracketing the root */
		tries = rootfind_bracket(fn, &lower_x, &upper_x, bracket_factor, bracket_tries, find_y, vp);
		if (tries < 0)
			return tries;
	}

	/* If we're not solving for zero, fix-up the callback */
	if (find_y != 0.0)
	{
		/* Use rootfind_nonzero as a proxy to offset the whole curve */
		nonzero.fn = fn;
		nonzero.dfndx = dfndx;
		nonzero.vp = vp;
		nonzero.diff_y = find_y;
		fn = rootfind_nonzero;
		dfndx = rootfind_nonzero_deriv;
		vp = &nonzero;
	}

	/* Call the appropriate solver */
	switch (rootfind_flag)
	{
		case eROOTFIND_BISECTION:
			/* Bisection method */
			return rootfind_bisection(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_SECANT:
			/* Secant method */
			return rootfind_secant(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_FALSE_POSITION:
			/* False position method */
			return rootfind_false_position(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_RIDDERS:
			/* Ridders' method */
			return rootfind_ridders(fn, lower_x, NULL, upper_x, NULL,  accuracy, max_tries, vp, found_x);

		case eROOTFIND_BRENT:
			/* Brent's method */
			return rootfind_brent(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_CBLFR:
			/* CBLFR method */
			return rootfind_cblfr(fn, lower_x, NULL, upper_x, NULL, NULL, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_NEWTON:
			/* Newtons method. */
			//assert(dfndx != NULL);
			return rootfind_newton(fn, dfndx, (lower_x + upper_x) / 2.0, NULL, NULL, accuracy, max_tries, vp, found_x);

		case eROOTFIND_BRENT_GROWRANGE:
			return rootfind_brent_grow_range(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x, 1);	
			;

		case eROOTFIND_BRENT_NOGROW:
			return rootfind_brent_grow_range(fn, lower_x, NULL, upper_x, NULL, accuracy, max_tries, vp, found_x, 0);
			;
	}
	throw(  "input error");
}



//1.0e-20;

void MlEqMaths::ludcmp(CMatrix& a, GVector<int>& indx, double& d)
{
	const double TINY = 1e-30;
	int n = a.rows();
	int i,imax,j,k;
	double big,dum,sum,temp;

	CVector vv(n);

	indx.resize(n);

	d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i-1][j-1])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i-1]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i-1][j-1];
			for (k=1;k<i;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i-1][j-1];
			for (k=1;k<j;k++)
				sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
			if ( (dum=vv[i-1]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax-1][k-1];
				a[imax-1][k-1]=a[j-1][k-1];
				a[j-1][k-1]=dum;
			}
			d = -(d);
			vv[imax-1]=vv[j-1];
		}
		indx[j-1]=imax;
		if (a[j-1][j-1] == 0.0) a[j-1][j-1]=TINY;
		if (j != n) {
			dum=1.0/(a[j-1][j-1]);
			for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
		}
	}
//	free_vector(vv,1,n);
}

void MlEqMaths::lubksb(CMatrix& a, GVector<int>& indx, CVector& b)
{
	int n = a.rows();
	if ( indx.getsize() == 0 ){
		indx.resize(n);
	}else if ( indx.getsize() != n ){
		throw("incorrect size of indx array encountered");
	}
				
	int i,ii=0,ip,j;
	double sum;
	for (i=1;i<=n;i++) {
		ip=indx[i-1];
		sum=b[ip-1];
		b[ip-1]=b[i-1];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
		else if (sum) ii=i;
		b[i-1]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i-1];
		for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1]=sum/a[i-1][i-1];
	}
}



/****************************************************************************/


int rootfind_solve_ex(int rootfind_flag, bool (*fn)(double x,void* vp,double* f), double accuracy, int max_tries, double find_y, 
                      void* vp, double* found_x, int (*dfndx)(double x,void* vp,double* f) , double x1, const double* y1,
                      const double* dydx1, const double* x2, const double* y2, const double* dydx2, const double* x3,
                      const double* y3, const double* dydx3)
{
	rootfind_nonzero_t nonzero;


	/* If we're not solving for zero, fix-up the callback */
	if (find_y != 0.0)
	{
		/* Use rootfind_nonzero as a proxy to offset the whole curve */
		nonzero.fn = fn;
		nonzero.dfndx = dfndx;
		nonzero.vp = vp;
		nonzero.diff_y = find_y;
		fn = rootfind_nonzero;
		dfndx = rootfind_nonzero_deriv;
		vp = &nonzero;
	}

	/* Call the appropriate solver */

	switch (rootfind_flag)
	{
		case eROOTFIND_BISECTION:
			/* Bisection method */
			//assert(x2 != NULL);
			return rootfind_bisection(fn, x1, y1, *x2, y2, accuracy, max_tries, vp, found_x);

		case eROOTFIND_SECANT:
			/* Secant method */
			//assert(x2 != NULL);
			return rootfind_secant(fn, x1, y1, *x2, y2, accuracy, max_tries, vp, found_x);

		case eROOTFIND_FALSE_POSITION:
			/* False position method */
			//assert(x2 != NULL);
			return rootfind_false_position(fn, x1, y1, *x2, y2, accuracy, max_tries, vp, found_x);

		case eROOTFIND_RIDDERS:
			/* Ridders' method */
			//assert(x2 != NULL);
			return rootfind_ridders(fn, x1, y1, *x2, y2, accuracy, max_tries, vp, found_x);

		case eROOTFIND_BRENT:
			/* Brent's method */
			//assert(x2 != NULL);
			return rootfind_brent(fn, x1, y1, *x2, y2, accuracy, max_tries, vp, found_x);

		case eROOTFIND_CBLFR:
			/* CBLFR method */
			//assert(x2 != NULL);
			return rootfind_cblfr(fn, x1, y1, *x2, y2, x3, y3, accuracy, max_tries, vp, found_x);

		case eROOTFIND_NEWTON:
			/* Newtons method. */
			//assert(dfndx != NULL);
			return rootfind_newton(fn, dfndx, x1, y1, dydx1, accuracy, max_tries, vp, found_x);
	}
	throw(  "input error");

}



#undef NR_ITMAX
#undef NR_EPS
#undef NR_FPMIN

#define EPS 3.0e-11

#undef EPS


void MlEqMaths::rkck(double y[], double dydx[], int n, double x, double h, double yout[],
	double yerr[], void (*derivs)(double, double [], double[], void *),
	void *vp)
{
	int i;
	static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;
	double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

	ak2=dvector(1,n);
	ak3=dvector(1,n);
	ak4=dvector(1,n);
	ak5=dvector(1,n);
	ak6=dvector(1,n);
	ytemp=dvector(1,n);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(*derivs)(x+a2*h,ytemp,ak2,vp);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(*derivs)(x+a3*h,ytemp,ak3,vp);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(*derivs)(x+a4*h,ytemp,ak4,vp);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(*derivs)(x+a5*h,ytemp,ak5,vp);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(*derivs)(x+a6*h,ytemp,ak6,vp);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
	free_dvector(ytemp,1,n);
	free_dvector(ak6,1,n);
	free_dvector(ak5,1,n);
	free_dvector(ak4,1,n);
	free_dvector(ak3,1,n);
	free_dvector(ak2,1,n);
}



#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void MlEqMaths::rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
		  double yscal[], double *hdid, double *hnext,
		  void (*derivs)(double, double [], double[], void *),
		  void *vp)
{

	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;
	for (;;) {
		MlEqMaths::rkck(y,dydx,n,*x,h,ytemp,yerr,derivs,vp);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? DMAX(htemp,0.1*h) : DMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) {
			nrerror("stepsize underflow in rkqs");
		}
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON

// in gasdev, "if (*idum < 0) iset=0" is an added statement
// also uses ran2 (NOT ran1 as in NR)

double MlEqMaths::gasdev(long *idum)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}	

	

double MlEqMaths::gasdev1(long *idum,double lo,double lup)
{	
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	double rand1,rand2;
	
	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {

			rand1 = (lup-lo)*ran2(idum)+lo;
//			rand2 = (lup-lo)*(1.0-MlEqMaths::ran2(&idum)+0.5)+lo;
			rand2 = ran2(idum);

			v1=2.0*rand1-1.0;
			v2=2.0*rand2-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);


		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}	


double MlEqMaths::RoundDouble(double doValue, int nPrecision)
{
	static const double doBase = 10.0;
	double doComplete5, doComplete5i;
	
	doComplete5 = doValue * pow(doBase, (double) (nPrecision + 1));
	
	if(doValue < 0.0)
		doComplete5 -= 5.0;
	else
		doComplete5 += 5.0;
	
	doComplete5 /= doBase;
	modf(doComplete5, &doComplete5i);
	
	return doComplete5i / pow(doBase, (double) nPrecision);
}



void fitPolynomial::initialize(CVector& initialGuess,CMatrix& data,
					   CMatrix& xValBounds,CMatrix& ObjectiveBounds,
					   double InitialTolerance,double FinalTolerance,
					   double StoppingTolerance,
					   CMatrix& NonZeroPartialDerivativesSpecification,int outputFlag)
{

	if ( ObjectiveBounds.rows() == 0 )
	{
		ObjectiveBounds.resize(1,2);
		ObjectiveBounds[0][0] = -1e99;
		ObjectiveBounds[0][1] = 1e99;
	}

	m_degree = initialGuess.getsize()-1;
	m_data = data;

	iVector linearVariableIndex;
/*	iVector linearVariableIndex(initialGuess.getsize());
	for ( int i = 0 ; i < linearVariableIndex.getsize(); i++ )
	{
		linearVariableIndex[i] = 1;
	}

*/
	LSGRGSolver::initialize(  initialGuess,linearVariableIndex,
							  xValBounds,ObjectiveBounds,
							  InitialTolerance,FinalTolerance,
							  StoppingTolerance,
							  NonZeroPartialDerivativesSpecification,outputFlag);
}


void  fitPolynomial::ObjectiveFcn(double* gVals,double* xVals)
{

		double objVal = 0.0,z,t;
		int idata,i;

		for (  idata = 0; idata < m_data.rows(); idata++ )
		{
			z = 0.0;
			t = m_data[idata][0];
			for ( i = 0; i <= m_degree; i++ )
			{
				z += xVals[i]*pow(t,(double)i);
			}

			objVal += pow(z-m_data[idata][1],2.0);
		}

		*gVals = objVal;
}




void randomGenerator::throwNext(CVector& randoms,int fillflag,int useRandoms)
{
	
	if ( useRandoms )
	{

		for ( int i = 0 ; i < randoms.getsize(); i++ )
		{
			if ( m_icounter == randoms.getsize() )
				m_icounter = 0;

			randoms[i] = m_filledRandoms[m_icounter][i];

			m_icounter++;
		}
	}
	else
	{
		generateRandoms(randoms);

		if ( fillflag )
		{
			for ( int i = 0 ; i < m_filledRandoms.cols(); i++ )
			{
				m_filledRandoms[m_icounter][i] = randoms[i];
			}
		}
	}

	return;
}


double randomGenerator::throwNext(int fillflag,int useRandoms)
{
	
	double res=0;

	if ( useRandoms )
	{
		if ( m_icounter == 1)
				m_icounter = 0;

		res = m_filledRandoms[m_icounter][0];
		m_icounter++;
		
	}
	else
	{
		res = generateRandom();
		if ( fillflag ){
			m_filledRandoms[m_icounter][0] = res;
		}
	}

	return res;
}


void randomGenerator::generateRandoms(CVector& randoms,int ipath)
{
	for ( int i = 0 ; i < randoms.getsize(); i++ ){
		randoms[i] = MlEqMaths::gasdev(&m_seed);
	}
}


void randomGenerator::generateRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp,bool generateGauss)
{

//  randoms[iasset][istep*numberOfFactors[iasset]+jfactor]

	int istep,iasset,ifactor,n,m;

	if ( ipath == 0 )
	{
//		check for consistency of matrix set up once :

		if ( numberOfFactors.getsize() != randoms.getsize() ){
			throw("incorrect dimensioning of arrays encountered");
		}

		for ( iasset = 0 ; iasset < randoms.getsize(); iasset++ )
		{
			n = randoms[iasset].getsize();
			m = numberOfFactors[iasset];
			if ( n%m != 0 ){
				throw("incorrect dimensioning encountered in generateRandoms");
			}
			if ( randoms[0].getsize()/numberOfFactors[0] != randoms[iasset].getsize()/numberOfFactors[iasset] ){
				throw("incorrect dimensioning encountered in generateRandoms");
			}
		}
	}

	int nsteps = randoms[0].getsize()/numberOfFactors[0];

	int ndim = rndTemp.getsize();
	if ( rndTemp.getsize() == 0 )
	{
		for ( iasset = 0 ; iasset < randoms.getsize(); iasset++ ){
			ndim += numberOfFactors[iasset];
		}
		ndim *= nsteps;
		rndTemp.resize(ndim);
	}

	if ( generateGauss ){
		generateRandoms(rndTemp,ipath);// this is virtual function call
	}
	else
	{
		generateUniformRandoms(rndTemp,ipath);// this is virtual function call
	}

	int factor = 0;
	for ( iasset = 0 ; iasset < randoms.getsize(); iasset++ ){
		factor += numberOfFactors[iasset];
	}


	int k,l;
	for ( istep = 0 ; istep < nsteps; istep++ )
	{
		m = istep*factor; 
		k = -1;
		for ( iasset = 0 ; iasset < randoms.getsize(); iasset++ )
		{
			for (ifactor = 0; ifactor < numberOfFactors[iasset]; ifactor++ )
			{	
				k++;
				n = istep*numberOfFactors[iasset] + ifactor;
				randoms[iasset][n] = 0.0;
				for ( l = 0; l <= k; l++ ){	
					randoms[iasset][n] += cholesky.getCholesky(istep,k,l)*rndTemp[l+m];
				}	
			}		
		}			
	}				
}





void randomGenerator::generateUniformRandoms(CVector& randoms,int ipath)
{
	double x;
	for ( int i = 0 ; i < randoms.getsize(); i++ )
	{
		x = MlEqMaths::gasdev(&m_seed);
		randoms[i] = normal(x);
	}
}


void randomGenerator::generateUniformRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp)
{
	generateRandoms(randoms,numberOfFactors, cholesky,ipath,rndTemp,false);
}








randomGenerator::randomGenerator(long idum,int numberScenariosToBeStored,int numberRandomsPerScenario)
{
	initialize(idum,numberScenariosToBeStored,numberRandomsPerScenario);
}


double randomGenerator::generateRandom()
{
	return MlEqMaths::gasdev(&m_seed);
}

double randomGenerator::generateUniformRandom()
{
	double x = MlEqMaths::gasdev(&m_seed);
	x = normal(x);
	return x;
}



void randomGenerator::initialize(long idum,int numberScenariosToBeStored,int numberRandomsPerScenario)
{
	m_seed		= idum;
	m_icounter	= 0;

	if ( numberScenariosToBeStored )
	{
		if ( numberRandomsPerScenario == 0 ){
			throw("error in setting up random numbers");
		}

		m_filledRandoms.resize(numberScenariosToBeStored,numberRandomsPerScenario);
	}

}

randomGenerator::~randomGenerator()
{
	m_filledRandoms.resize(0,0);
}

void randomGenerator::average(CVector& results,int numberPaths)
{
	for ( int i = 0 ; i < results.getsize(); i++ ){
		results[i] /= numberPaths;
	}
}




void Diophantine::generateRandoms(CVector& randoms,int ipath)
{
	ipath += offset();
	m_pIntegrator->getIntegrationPoint( ipath, m_line,*m_pDomain,*m_pParameters);

	for ( int i = 0 ; i < randoms.getsize();i++)
	{
		randoms[i] = m_line[i];
		randoms[i] = normal_inv(randoms[i]);
	}
}

void Diophantine::generateUniformRandoms(CVector& randoms,int ipath)
{
	ipath += offset();
	m_pIntegrator->getIntegrationPoint( ipath, m_line,*m_pDomain,*m_pParameters);

	for ( int i = 0 ; i < randoms.getsize();i++)
	{
		randoms[i] = m_line[i];
	}
}


void createCholesky(CMatrix &cholesky,CMatrix& correl)
{
	// remember lower triangle of correl has been changed !

	int m_numberOfFactors = correl.rows();
	CVector diag(m_numberOfFactors);
	MlEqMaths::choldc(correl,diag);
	int ifactor,jfactor;

	cholesky.resize(m_numberOfFactors,m_numberOfFactors);
	for ( ifactor = 0 ; ifactor < m_numberOfFactors; ifactor++ )
		{
		  for ( jfactor = 0 ; jfactor <= ifactor; jfactor++ )
			{
				if ( ifactor == jfactor ){
					cholesky[ifactor][jfactor]  += diag[ifactor];
				}
				else{
					cholesky[ifactor][jfactor]	+= correl[ifactor][jfactor];
				}
			}
		}
}


void Cholesky::initialize(CMatrix& correl,	int numberOfPeriods)
{
	m_numberOfPeriods	= numberOfPeriods;
	m_numberOfFactors	= correl.rows();

	createCholesky(m_cholesky,correl);
}

Cholesky::Cholesky(CMatrix& correl,int numberOfPeriods)
{
	initialize(correl,numberOfPeriods);
}

double Cholesky::getCholesky(int isclice,int ifactor,int jfactor)
{
	return m_cholesky[ifactor][jfactor];
}



void Diophantine::average(CVector& results,int numberPaths)
{
	double C = m_pWeights->C(m_Neff);
	double factor = C;

	if ( m_AmIsConst )
	{
      double Am = m_pWeights->A(m_Neff, 0);
	  factor *= Am;
	}
	if ( m_JacobiIsConst ){
		factor *= m_jacobi;
	}
	  
	for ( int i = 0 ; i < results.getsize(); i++ ){
		results[i] *= factor;
	}
}


void Diophantine::completeValue(double& payoff,int ipath)
{
	ipath += offset();
	double fac=1.0;
	if ( m_AmIsConst == 0 ){
		fac = m_pWeights->A(m_Neff, ipath);
	}

	if ( m_JacobiIsConst == 0 )
	{
		m_pIntegrator->getIntegrationPoint( ipath,m_line,*m_pDomain,*m_pParameters);
		fac *= m_pPeriodization->product(m_line);
	}

	if ( !m_JacobiIsConst || !m_AmIsConst ){
		payoff *= fac;
	}

}

int Diophantine::offset()
{
	return m_Nrange.first;
}

void Diophantine::initialize(	
								int ndimension,
								int npaths,
								CAlpha * pAlpha ,	
								CWeights * pWeights,
								CDomain * pDomain,
								CPeriodization * pPeriodization,
								CNTintegrator* pIntegrator,
								CNTparameters* pParameters)
{


	m_pAlpha			= pAlpha;
	m_pWeights			= pWeights;
	m_pDomain			= pDomain;
	m_pPeriodization	= pPeriodization;
	m_pIntegrator		= pIntegrator;
	m_pParameters		= pParameters;
	m_numberDimensions  = ndimension;
	m_nPaths			= npaths;

	m_line.resize(ndimension);
	m_Neff	= m_pWeights->Neff(npaths);
	m_pWeights->Nrange(m_Nrange, m_Neff);



	// check for constancy of jacobioan
	m_JacobiIsConst = m_pPeriodization->jacobianIsConstant();

	if ( m_JacobiIsConst )
	{
		int ipath = m_Nrange.first ;
		m_pIntegrator->getIntegrationPoint( ipath, m_line, *m_pDomain, *m_pParameters);	
		m_jacobi = m_pPeriodization->product(m_line);
	}

	// check for constancy of Am
	m_AmIsConst = m_pWeights->AisConstant();
}

Diophantine::Diophantine(int ndimension,
				int npaths,
				CAlpha * pAlpha ,	
				CWeights * pWeights,
				CDomain * pDomain,
				CPeriodization * pPeriodization,
				CNTintegrator* pIntegrator,
				CNTparameters* pParameters)
{

	initialize(ndimension,npaths,pAlpha ,pWeights,pDomain,pPeriodization,
				pIntegrator,pParameters);
}




void MlEqMaths::odeint(CVector& ystart,int nvar,double x1,double x2,double eps,double h1,double hmin,int*nok,int* nbad,
									void (*derivs)(double,double[],double[],void *),
									void (*rkqs)(double[],double[],int,double*,double,double,double[],double*,double*,void (*)(double,double[],double[],void *), void *),
									double* xp,double** yp,int& kount,void *vp,const int kmax,double dxsav,double scale )
{

	int nstp,i;
	double xsav,x,hnext,hdid,h;
	CVector yscal(nvar+1);
	CVector y(nvar+1);
	CVector dydx(nvar+1);
	
	const int MAXSTP = 10000;
	const double TINY = 1e-30+eps/scale;
	x = x1;
	h = MSIGN(h1,x2-x1);//sos
	*nok = (*nbad) =kount = 0;
	
	for (i = 1; i <= nvar; i++) y[i]=ystart[i];
	if(kmax > 0) xsav = x-dxsav*2.0;
	for (nstp = 1; nstp <= MAXSTP; nstp++)
	{
		(*derivs)(x,y.getPtr(),dydx.getPtr(),vp); 
		for (i = 1 ; i<= nvar; i++)
			yscal[i] = MlEqMaths::Max(fabs(y[i])+fabs(dydx[i]*h),TINY);
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav))
		{
			xp[++kount] = x;
			for (i = 1; i<=nvar; i++) yp[i][kount] = y[i];
			xsav = x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h = x2-x;
		(*rkqs)(y.getPtr(),dydx.getPtr(),nvar,&x,h,eps,yscal.getPtr(),&hdid,&hnext,derivs,vp);
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0)
		{
			for (i = 1; i<= nvar; i++) ystart[i] = y[i];
			if (kmax)
			{
				xp[++kount] = x;
				for (i = 1; i<= nvar; i++) yp[i][kount] = y[i];
			}
			return;
		}
		if (fabs(hnext) <= hmin) 
		{
			hnext = hmin;
		}
		h = hnext;
	}
	
	throw("Too many steps in routine odeint" );
}



void MlEqMaths::svdcmp(CMatrix& a, CVector& w, CMatrix& v)
{
	int n = a.cols();
	int m = a.rows();

	if ( w.getsize() == 0 ){
		w.resize(n);
	}else if ( w.getsize() != n ){
		throw("w array has wrong size");
	}
	if ( v.rows() == 0 ){
		v.resize(n,n);
	}
	if ( v.cols() != v.rows() ){
		throw("input matrix v must be a square matrix");
	}
	else if ( v.rows() != n ){
		throw(" v matrix has wrong size");
	}

	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;

	CVector rv1(n);
//	rv1=dvector(0,n-1);

	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i-1]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k-1][i-1]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k-1][i-1] /= scale;
					s += a[k-1][i-1]*a[k-1][i-1];
				}
				f=a[i-1][i-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i-1][i-1]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
					f=s/h;
					for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
				}
				for (k=i;k<=m;k++) a[k-1][i-1] *= scale;
			}
		}
		w[i-1]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i-1][k-1]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i-1][k-1] /= scale;
					s += a[i-1][k-1]*a[i-1][k-1];
				}
				f=a[i-1][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i-1][l-1]=f-g;
				for (k=l;k<=n;k++) rv1[k-1]=a[i-1][k-1]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j-1][k-1]*a[i-1][k-1];
					for (k=l;k<=n;k++) a[j-1][k-1] += s*rv1[k-1];
				}
				for (k=l;k<=n;k++) a[i-1][k-1] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i-1])+fabs(rv1[i-1])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i-1][k-1]*v[k-1][j-1];
					for (k=l;k<=n;k++) v[k-1][j-1] += s*v[k-1][i-1];
				}
			}
			for (j=l;j<=n;j++) v[i-1][j-1]=v[j-1][i-1]=0.0;
		}
		v[i-1][i-1]=1.0;
		g=rv1[i-1];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i-1];
		for (j=l;j<=n;j++) a[i-1][j-1]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
				f=(s/a[i-1][i-1])*g;
				for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
			}
			for (j=i;j<=m;j++) a[j-1][i-1] *= g;
		} else for (j=i;j<=m;j++) a[j-1][i-1]=0.0;
		++a[i-1][i-1];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l-1])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm-1])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i-1];
					rv1[i-1]=c*rv1[i-1];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i-1];
					h=pythag(f,g);
					w[i-1]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j-1][nm-1];
						z=a[j-1][i-1];
						a[j-1][nm-1]=y*c+z*s;
						a[j-1][i-1]=z*c-y*s;
					}
				}
			}
			z=w[k-1];
			if (l == k) {
				if (z < 0.0) {
					w[k-1] = -z;
					for (j=1;j<=n;j++) v[j-1][k-1] = -v[j-1][k-1];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l-1];
			nm=k-1;
			y=w[nm-1];
			g=rv1[nm-1];
			h=rv1[k-1];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i-1];
				y=w[i-1];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j-1]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj-1][j-1];
					z=v[jj-1][i-1];
					v[jj-1][j-1]=x*c+z*s;
					v[jj-1][i-1]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j-1]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj-1][j-1];
					z=a[jj-1][i-1];
					a[jj-1][j-1]=y*c+z*s;
					a[jj-1][i-1]=z*c-y*s;
				}
			}
			rv1[l-1]=0.0;
			rv1[k-1]=f;
			w[k-1]=x;
		}
	}

}



double ARTrackingPursue::scalarProduct(CVector&x, CVector& y)
{
	double res = 0;
	if ( x.getsize() != y.getsize() ){
		throw("scalarproduct of vectors with different sizes attemped");
	}
	for ( int i = 0 ; i < x.getsize(); i++ ){
		res += x[i]*y[i];
	}
	res /= x.getsize();
	return res;
}

double ARTrackingPursue::scalarProduct(int j, CVector& x)
{
	double res = 0;
	for ( int i = 0 ; i < x.getsize(); i++ ){
		res += x[i]*m_projections[i][j];
	}
	res /= x.getsize();
	return res;
}

void ARTrackingPursue::multiplyVec(CVector& x,double f)
{
	for (int i = 0; i < x.getsize(); i++ ){
		x[i] *= f;
	}
}

double ARTrackingPursue::variance(CVector& x)
{
	int ndim = x.getsize();
	double var = 0.0;
	for (int i = 0; i < ndim; i++ ){
		var  += x[i]*x[i];
	}
	var  /= ndim;
	var   = sqrt(var);
	return var;
}





void ARTrackingPursue::project()
{
	int i,j,iter,ivec;
	int ndim = m_weights.getsize();
	int nvec = m_projections.cols();
	double proj,maxproj,var;

	CVector x = m_targetVec;
	double initialVar = variance(x);	
	for ( iter = 0 ; iter < m_maxiter; iter++ )
	{	
//		determine maximal projection
		j = 0;			
		maxproj =  scalarProduct(0, x)/(m_stdev[0]*m_stdevx);
		for ( ivec = 1; ivec < nvec; ivec++ )
		{	
			proj =  scalarProduct(ivec, x)/(m_stdev[ivec]*m_stdevx);
			if ( proj > maxproj ){
				maxproj = proj;
				j		= ivec;
			}
		}
		
//		determine upgrated vector
		maxproj *= m_stdev[j]*m_stdevx;
		double z		=  m_epsilon;
		m_weights[j]	+= z;
		for ( i = 0 ; i < ndim; i++ ){	
			x[i] -= z*m_projections[i][j];	
		}	
//		check stopping criteria
		var = variance(x);
//		if ( var > initialVar ){
//			throw("weird result encountered in ARtracking pursue algorithm");
//		}
		if ( var <= initialVar*m_stoppingTolerance ){
			break;
		}
						
	}

	if ( iter == m_maxiter ){
		throw("ARtracking algorithm has not converged");
	}

}


void ARTrackingPursue::initialize(CVector& targetVec,CMatrix& projections,double epsilon,double varReduction,int maxiter)
{
	m_targetVec			= targetVec;
	m_projections		= projections;
	m_epsilon			= epsilon;
	m_stoppingTolerance = (1.0-varReduction);
	m_maxiter			= maxiter;
	m_weights.resize(m_projections.rows());

	if ( m_projections.rows() != targetVec.getsize() ){
		throw("projection vectors have inconsitent arraysize");
	}

	if ( m_weights.getsize() == 0 ){
		m_weights.resize(m_projections.cols());
	}

	m_stdev.resize(m_projections.cols());

	int j;
	double var;
	for ( int i = 0; i < m_stdev.getsize(); i++ )
	{
		var = 0.0;
		for ( j = 0; j < m_projections.rows(); j++ ){
			var += pow(m_projections[j][i],2.0);
		}
		var /= m_projections.rows();
		m_stdev[i] = sqrt(var);
	}
	
	m_stdevx = 0.0;
	for ( int i = 0; i < m_targetVec.getsize(); i++ ){
		m_stdevx += pow(m_targetVec[i],2.0);
	}
	m_stdevx /= m_targetVec.getsize();
	m_stdevx = sqrt(m_stdevx);
	
}
	
	
bool invertMatrix(CMatrix& inverse,const CMatrix& in)
{	
	bool isInvertible = true;
	int n;	
	if ( in.rows() != in.cols() ){
		throw("only square matrices can be inverted, Franzl");
	}
	int ndim = in.rows();

	inverse.resize(ndim,ndim);
	CVector w;
	CMatrix v,u;

	u = in;
	MlEqMaths::svdcmp(u,w,v);

//  construct the inverse now

	double eps = 1e-30;

	double wmax = w[0];
	for ( n = 1 ; n < ndim; n++ )
	{
		if ( fabs(w[n]) > wmax ) {
			wmax = fabs(w[n]);
		}
	}
	for ( n = 0 ; n < ndim; n++ )
	{
		if ( fabs(w[n]) < eps*wmax ) {
			w[n] = 0.0;
			isInvertible = false;
		}
	}

	for ( int i = 0 ; i < ndim; i++ )
	{
		for ( int j = 0 ; j < ndim; j++ )
		{
			for ( int k = 0 ; k < ndim; k++ )
			{
				if ( w[k] > eps){
					inverse[i][j] += v[i][k]/w[k]*u[j][k];
				}
			}
		}
	}

	return isInvertible;
}

static void dgauleg(double x1, double x2,double* x, double*  w, int n)
{
//  arrays start from 1
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
}


void MlEqMaths::dGauleg( double x1, double x2, CVector& x, CVector& w, int n ,bool normaliseWeights /*=false*/)
{
	dgauleg(x1,x2,x.getPtr()-1, w.getPtr()-1, n);

	if (normaliseWeights)
	{
		double sq = 1.0/sqrt(2.0*3.141592654);	
		double fac = 0.0;	
		for (int i = 0 ; i < n; i++ ){
			fac += w[i]*sq*exp(-0.5*pow(x[i],2.0));
		}
		fac = ( normal(x2)-normal(x1) )/fac;
		for ( int i = 0 ; i < n; i++ ){
			w[i] *= fac;
		}
	}
}



void merge( CVector& outPutSet, CVector& inPutArray_a,CVector& inPutArray_b ,double tolerance)
{

	throw( "the tolerance function needs to be implemented");

	std::vector<double> stl_inPutArray_a;
	std::vector<double> stl_inPutArray_b;

	STLVectorFromCVector(stl_inPutArray_a,inPutArray_a);
	STLVectorFromCVector(stl_inPutArray_b,inPutArray_b);
	std::vector<double> stl_outPutSet;

	merge(stl_outPutSet, stl_inPutArray_a,stl_inPutArray_b ,tolerance);

	CVectorFromSTLVector(outPutSet,stl_outPutSet);
}


void merge( std::vector<double>& outPutSet, const std::vector<double>& inPutArray_a,const std::vector<double>& inPutArray_b ,double tolerance)
{
	std::map<double, bool> mapOut;
	
	for (std::vector<double>::const_iterator it = inPutArray_a.begin(); it != inPutArray_a.end(); it++){
		mapOut[*it] = true;
	}

	for (std::vector<double>::const_iterator it = inPutArray_b.begin(); it != inPutArray_b.end(); it++){
		mapOut[*it] = true;
	}

	outPutSet.clear();
	for (std::map<double, bool>::const_iterator it = mapOut.begin(); it != mapOut.end(); it++){
		outPutSet.push_back(it->first);
	}
	
}



void merge( std::vector<long>& outPutSet, const std::vector<long>& inPutArray_a,const std::vector<long>& inPutArray_b )
{
	std::map<long, bool> mapOut;
	
	for (std::vector<long>::const_iterator it = inPutArray_a.begin(); it != inPutArray_a.end(); it++){
		mapOut[*it] = true;
	}

	for (std::vector<long>::const_iterator it = inPutArray_b.begin(); it != inPutArray_b.end(); it++){
		mapOut[*it] = true;
	}

	outPutSet.clear();
	for (std::map<long, bool>::const_iterator it = mapOut.begin(); it != mapOut.end(); it++){
		outPutSet.push_back(it->first);
	}

	
}


void merge( GVector<long>& outPutSet, GVector<long>& inPutArray_a,GVector<long>& inPutArray_b )
{

	std::vector<long> stl_inPutArray_a(inPutArray_a.getsize());
	std::vector<long> stl_inPutArray_b(inPutArray_b.getsize());

	for ( int i = 0 ; i < inPutArray_a.getsize(); i++ ){
		stl_inPutArray_a[i] = inPutArray_a[i];
	}

	for ( int i = 0 ; i < inPutArray_b.getsize(); i++ ){
		stl_inPutArray_b[i] = inPutArray_b[i];
	}

	std::vector<long> stl_outPutSet;
	merge( stl_outPutSet,stl_inPutArray_a,stl_inPutArray_b );

	outPutSet.resize(stl_outPutSet.size());
	for ( int i = 0 ; i < outPutSet.getsize(); i++ ){
		outPutSet[i] = stl_outPutSet[i];
	}

}

void merge( std::vector<long>& outPutSet,std::vector<int>& mapping, const std::vector<long>& inPutArray_a,const std::vector<long>& inPutArray_b )
{
	merge( outPutSet,inPutArray_a,inPutArray_b );
	
	mapping.resize(inPutArray_a.size());
	for ( int i = 0 ; i < inPutArray_a.size(); i++ )
	{
		for ( int k = 0 ; k < outPutSet.size(); k++ )
		{
			if ( inPutArray_a[i] == outPutSet[k] )
			{
				mapping[i] = k;
				break;
			}

			if ( k == outPutSet.size() ){
				throw("mapping failed");
			}

		}
	}


}


void RKI_func_wrapper(double x, double y[], double dydx[], void *vp) 
{
	RKI_class *p = (RKI_class *) vp;
	p->func(x, dydx, p->data);
}

void Runge_Kutta_Integrate_new_vector(CVector& integral, void (*fcn)(double, double[], void *p), const CVector& limits, double h_min, double eps, double bump, void *vp, double h_init, double scale)
{ 		
  int nok,nbad,kount,n,k;
  double dxsav = 0.0;
 
  RKI_class rki;

  rki.func = fcn;
  rki.data = vp;
   	
  int nvar = integral.getsize();
  	
  CVector ystart(nvar+1);
  for ( n = 0 ; n <= nvar; n++ )
  	ystart[n] = 0.0;
  	
  for ( n = 0 ; n < nvar; n++ )
	integral[n] = 0.0;
	  
  int nlimits = limits.getsize()-1;
  if ( limits[0] >= limits[nlimits] )
  {			
  	  return;
  }		
  for ( n = 0 ; n < nlimits; n++ )
  {		
	  for ( k = 0 ; k <= nvar; k++ )
		ystart[k] = 0.0;
		
	  MlEqMaths::odeint(ystart,nvar,limits[n],limits[n]*(1.0+bump+EPS),1e-8,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];


	  for ( k = 1 ; k <= nvar; k++ )
		ystart[k] = 0.0;

	  MlEqMaths::odeint(ystart,nvar,limits[n+1]*(1.0-bump+EPS),limits[n+1],1e-8,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];

	  for ( k = 1 ; k <= nvar; k++ )
		ystart[k] = 0.0;

	  MlEqMaths::odeint(ystart,nvar,limits[n]*(1.0+bump+EPS),limits[n+1]*(1.0-bump+EPS),eps,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);
	  for ( k = 0 ; k < nvar; k++ )
		integral[k] += ystart[k+1];
  }		
} 		


void Runge_Kutta_Integrate_new(double& integral,void (*fcn)(double, double[], void *p),const CVector& limits,double h_min,double eps,double bump, void *vp, double h_init ,double scale)
{ 		

  CVector ystart(2);
  ystart[0] = ystart[1] = 0.0;
  int nok,nbad,kount;
  double dxsav = 0.0;

  RKI_class rki;

  rki.func = fcn;
  rki.data = vp;
   		
  int n;
  int nlimits = limits.getsize()-1;
  		
  if (limits[0] >= limits[nlimits])
  {		
  	  integral = 0.0;
  	  return;
  }
  			
  integral = 0.0;
  for (n = 0 ; n < nlimits; n++)
  {	  	
  	      	
  	  ystart[0]= 0.0;
  	  ystart[1] = 0.0;

	  MlEqMaths::odeint(ystart,1,limits[n],limits[n]*(1.0+bump+EPS),1e-8,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);
  	  integral += ystart[1];
  	  ystart[0]= 0.0;
  	  ystart[1] = 0.0; 
	  MlEqMaths::odeint(ystart,1,limits[n+1]*(1.0-bump+EPS),limits[n+1],1e-8,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);

	  integral += ystart[1];
  	  ystart[0]= 0.0;
  	  ystart[1] = 0.0; 

	  MlEqMaths::odeint(ystart,1,limits[n]*(1.0+bump+EPS),limits[n+1]*(1.0-bump+EPS),eps,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);
  	  integral += ystart[1];  		
  }		
} 		


void Runge_Kutta_Integrate_buckets(CVector& integralBuckets,void (*fcn)(double, double[], void *p),const CVector& limits,double h_min,double eps,double bump, void *vp, double h_init ,double scale)
{ 		

  CVector ystart(2);
  ystart[0] = ystart[1] = 0.0;
  int nok,nbad,kount;
  double dxsav = 0.0;

  RKI_class rki;

  rki.func = fcn;
  rki.data = vp;
   		
  int n;
  int nlimits = limits.getsize()-1;  		
  integralBuckets.resize(nlimits);

  if (limits[0] >= limits[nlimits]){		
  	  return;
  }
  			

  for (n = 0 ; n < nlimits; n++)
  {	  	
  	  ystart[0]= 0.0;
  	  ystart[1] = 0.0;

	  MlEqMaths::odeint(ystart,1,limits[n],limits[n+1],1e-8,h_init,h_min,&nok,&nbad,RKI_func_wrapper,MlEqMaths::rkqs,NULL,NULL,kount,(void *) &rki,0,dxsav,scale);

  	  integralBuckets[n] = ystart[1];  		
  }		

} 		




#undef NRANSI
#undef MIN_INT  
#undef INT_MAX  

