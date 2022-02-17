/* ==========================================================================
   FILE_NAME:	num_f_bessel.c

   PURPOSE:     compute modified Bessel function I and K 
                (cf NRC pages 237ff)    
   ========================================================================== */

#include   "utallhdr.h"
#include	<math.h"
#include   "num_h_gamma.h"

/* -------------------------------------------------------------------------
   Modified Bessel functions of integer order 
   ------------------------------------------------------------------------- */

/* Modified Bessel function I0(x) for any real x */

double bessi0(double x)
{
	double  ax,ans;
	double  y;

	if ((ax=fabs(x)) < 3.75) 
	{
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} 
	else 
	{
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

/* Modified Bessel function K0(x) for any real x */

double bessk0(double x)
{
	double bessi0(double x);
	double y,ans;

	if (x <= 2.0) 
	{
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	}
	else 
	{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans;
}

/* Modified Bessel function I1(x) for any real x */

double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) 
	{
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} 
	else 
	{
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

/* Modified Bessel function K1(x) for any real x */

double bessk1(double x)
{
	double bessi1(double x);
	double y,ans;

	if (x <= 2.0) 
	{
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} 
	else 
	{
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return ans;
}

/* Modified Bessel function In(x) for any real x and n > 2 */

double bessk(int n, double x)
{
	double bessk0(double x);
	double bessk1(double x);
	int j;
	double bk,bkm,bkp,tox;

	bkm = bessk0(x);
	bk  = bessk1(x);
	
	if (n == 0) 
	{
		return bkm;
	}
	else
	if (n == 1)
	{
		return bk;
	}
	
	tox=2.0/x;
	for (j=1;j<n;j++) 
	{
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}

	return bk;
}


#define   ACC     40.0
#define   BIGNO   1.0e10
#define   BIGNI   1.0e-10

/* Modified Bessel function In(x) for any real x and n >=0 */

double bessi(int n, double x)
{
	double  bessi0(double x);
	double  bessi1(double x);
	int     j;
	double  bi,bim,bip,tox,ans;

	if (n == 0) 
	{
		return bessi0(x);
	}
	else
	if (n == 1)
	{
		return bessi1(x);
	}
	if (x == 0.0)
	{
		return 0.0;
	}
	else 
	{
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) 
		{
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > BIGNO) 
			{
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if (j == n) 
				ans=bip;
		}
		ans *= bessi0(x)/bi;
		return x < 0.0 && (n & 1) ? -ans : ans;
	}
}


#undef   ACC
#undef   BIGNO
#undef   BIGNI

/* (C) Copr. 1986-92 Numerical Recipes Software. */


/* ========================================================================= */

/* -------------------------------------------------------------------------
   Modified Bessel functions of fractionnal order 
   ------------------------------------------------------------------------- */

#define   NUSE1          7
#define   NUSE2          8

#define   BESSEL_EPS     1.0e-16
#define   FPMIN          1.0e-60
#define   MAXIT          10000
#define   XMIN           2.0


/* -------------------------------------------------------------------------
   Evaluates Gamma1 and Gamma2 by Chebyshev expansion for fabs(x) < 0.5.
   Also returns 1/Gamma(1+x) and 1/Gamma(1-x).
   -------------------------------------------------------------------------- */
 
Err beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	double chebev(double a, double b, double c[], int m, double x);
	double xx;
	static double c1[] = {
		-1.142022680371172e0,
		 6.516511267076e-3,
		 3.08709017308e-4,
		-3.470626964e-6,
		 6.943764e-9,
		 3.6780e-11,
		-1.36e-13};
	static double c2[] = {
		 1.843740587300906e0,
		-0.076852840844786e0,
		 1.271927136655e-3,
		-4.971736704e-6,
		-3.3126120e-8,
		 2.42310e-10,
		-1.70e-13,
		-1.0e-15};

	xx = 8.0*x*x-1.0;
	*gam1 = chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2 = chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);

	return NULL;
}

/* -------------------------------------------------------------------------
   Returns the scaled modified Bessel function ri=e^{-x}*I_{nu}(x) , 
   rk=e^{x}*K_{nu}(x) and their derivatives rip and rkp for positive x and for xnu >0.
   The relative accuracy is within one or two significant digits of BESSEL_EPS.
   FPMIN is a number close to the machine smallest floating point number. 
   -------------------------------------------------------------------------- */

Err bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp)
{
	Err beschb(double x, double *gam1, double *gam2, double *gampl,
		double *gammi);
	int i,l,nl;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;

	if (x <= 0.0 || xnu < 0.0) 
	{
		return serror("Bad arguments in bessik");
	}
	nl=(int)(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) 
		h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) 
	{
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < BESSEL_EPS) 
			break;
	}
	if (i > MAXIT)
	{
		return serror("x too large in bessik; try asymptotic expansion");
	}
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--)
	{
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < XMIN)
	{
		x2=0.5*x;
		pimu=SRT_PI*xmu;
		fact = (fabs(pimu) < BESSEL_EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < BESSEL_EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*BESSEL_EPS) break;
		}
		if (i > MAXIT) 
		{
			return serror("bessk series failed to converge");
		}
		/* NAB: Changed to scale by e^{-x}
		rkmu=sum;
		rk1=sum1*xi2;
		*/
		rkmu=sum*exp(x);
        rk1=sum1*xi2*exp(x);
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < BESSEL_EPS) break;
		}
		if (i > MAXIT)
		{
			return serror("bessik: failure to converge in cf2");
		}
		h=a1*h;
		/* NAB changed to scale properly by e^x
		rkmu=sqrt(SRT_PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi; 
		*/
		rkmu=sqrt(SRT_PI/(2.0*x))/s;
        rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	*ri=(rimu*ril1)/ril;
	*rip=(rimu*rip1)/ril;
	for (i=1;i<=nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	*rk=rkmu;
	*rkp=xnu*xi*rkmu-rk1;

	return 0; 
}

#undef NUSE1
#undef NUSE2

#undef   BESSEL_EPS     
#undef   FPMIN          
#undef   MAXIT          
#undef   XMIN 


/* (C) Copr. 1986-92 Numerical Recipes Software . */

/* ------------------------------------------------------------------------ */
/* 
	 Wrapper function for bessik(NumRec) so we 
	 return only I_nu(x).
*/
Err I_nu(double xnu, double x, double *i_nu)
{

	double rk,rip,rkp;
	Err err = NULL;

	if (xnu >= 0) 
	{
		err = bessik(x,xnu, i_nu,&rk,&rip,&rkp);
		if (err)
			return err;
	}
	else
	{
		err = bessik(x,-xnu, i_nu,&rk,&rip,&rkp);
		if (err)
			return err;
		*i_nu = 2*exp(-2*x)*rk*sin(-xnu*SRT_PI)/SRT_PI + *i_nu;
	}
	
	return NULL;
}

/* ------------------------------------------------------------------------ */
/* 
	 wrapper function for bessik(NumRec) so we 
	 return only K_nu(x).
*/
Err K_nu(double xnu, double x, double *k_nu)
{

	double ri,rip,rkp;
	Err err = NULL;

	err = bessik(x,xnu,&ri,k_nu,&rip,&rkp);
	if (err)
		return err;

	return NULL;

}




