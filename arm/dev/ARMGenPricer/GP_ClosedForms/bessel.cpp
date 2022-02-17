/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *
 *	\file bessel.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>
#include <complex>
#include <limits>

#include "gpclosedforms/long_double.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/bessel.h"

#include "gpbase/numericconstant.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;


/// Implementation of the subroutine given in Numerical Recipes  in C++ for the computation ofthe Bessel Function J and Y and I and K 
/// for integer dimension and non integer (fractional) dimensions
/// the structure of the subroutine have been made a little bit more transparent and natural


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-12
#define ARM_CF_MAXIT 1000

inline double SIGN(const double &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

double chebev(const double a, const double b, vector<double> &c, const int m, const double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"x not in range in routine chebev");
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>0;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}


void beschb(const double x, double &gam1, double &gam2, double &gampl, double &gammi)
{
	const int NUSE1=7, NUSE2=8;
	int i;
	static const double c1_d[7] = {
		-1.142022680371168e0,6.5165112670737e-3,
		3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
		3.67795e-11,-1.356e-13};
	static const double c2_d[8] = {
		1.843740587300905e0,-7.68528408447867e-2,
		1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
		2.423096e-10,-1.702e-13,-1.49e-15};
	double xx;
	static vector<double> c1(7),c2(8);
	for(i=0;i<7;i++) c1[i]=c1_d[i];
	for(i=0;i<8;i++) c2[i]=c2_d[i];

	xx=8.0*x*x-1.0;
	gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	gampl= gam2-x*gam1;
	gammi= gam2+x*gam1;
}


double bessi0(const double x)
{
	double ax,ans,y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}


double bessi1(const double x)
{
	double ax,ans,y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y=3.75/ax;
		ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2));
		ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
		ans *= (exp(ax)/sqrt(ax));
	}
	return x < 0.0 ? -ans : ans;
}

double bessel_integer_I(const int n, const double x)
{
	const double ACC=200.0;
	const int IEXP=numeric_limits<double>::max_exponent/2;
	int j,k;
	double bi,bim,bip,dum,tox,ans;

	if (n ==0) return bessi0(x); 
	if (n ==1) return bessi1(x); 
	
	if (x*x <= 8.0*numeric_limits<double>::min()) return 0.0;
	else {
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for (j=2*(n+int(sqrt(ACC*n)));j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			dum=frexp(bi,&k);
			if (k > IEXP) {
				ans=ldexp(ans,-IEXP);
				bi=ldexp(bi,-IEXP);
				bip=ldexp(bip,-IEXP);
			}
			if (j == n) ans=bip;
		}
		ans *= bessi0(x)/bi;
		return x < 0.0 && (n & 1) ? -ans : ans;
	}
}

void bessik(const double x, const double xnu, double &ri, double &rk, double &rip, double &rkp)
{
	const int MAXIT=10000;
	const double EPS=numeric_limits<double>::epsilon();
	const double FPMIN=numeric_limits<double>::min()/EPS;
	const double CF_XMIN=2.0;
	const double CF_PI=3.141592653589793;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
	int i,l,nl;

	if (x <= 0.0 || xnu < 0.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bad arguments in bessik");
	nl=int(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=0;i<MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i >= MAXIT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"x too large in bessik; try asymptotic expansion");
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl-1;l >= 0;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < CF_XMIN) {
		x2=0.5*x;
		pimu=CF_PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,gam1,gam2,gampl,gammi);
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
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bessk series failed to converge");
		rkmu=sum;
		rk1=sum1*xi2;
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
		for (i=1;i<MAXIT;i++) {
			a -= 2*i;
			c = -a*c/(i+1.0);
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
			if (fabs(dels/s) <= EPS) break;
		}
		if (i >= MAXIT) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bessik: failure to converge in cf2");
		h=a1*h;
		rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	ri=(rimu*ril1)/ril;
	rip=(rimu*rip1)/ril;
	for (i=1;i <= nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	rk=rkmu;
	rkp=xnu*xi*rkmu-rk1;
}

double bessj0(const double x)
{
	double ax,z,xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934945152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

double bessj1(const double x)
{
	double ax,z,xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

double bessel_integer_J(const int n, const double x)
{
	const double ACC=160.0;
	const int IEXP=numeric_limits<double>::max_exponent/2;
	bool jsum;
	int j,k,m;
	double ax,bj,bjm,bjp,dum,sum,tox,ans;

	if (n ==0) return bessj0(x) ;
	if (n ==1) return bessj1(x) ;
	ax=fabs(x);
	if (ax*ax <= 8.0*numeric_limits<double>::min()) return 0.0;
	else if (ax > double(n)) {
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+int(sqrt(ACC*n)))/2);
		jsum=false;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			dum=frexp(bj,&k);
			if (k > IEXP) {
				bj=ldexp(bj,-IEXP);
				bjp=ldexp(bjp,-IEXP);
				ans=ldexp(ans,-IEXP);
				sum=ldexp(sum,-IEXP);
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans : ans;
}

void bessjy(const double x, const double xnu, double &rj, double &ry, double &rjp, double &ryp)
{
	const int MAXIT=10000;
	const double EPS=numeric_limits<double>::epsilon();
	const double FPMIN=numeric_limits<double>::min()/EPS;
	const double CF_XMIN=2.0;
	const double CF_PI=3.141592653589793;
	double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
		fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl,
		rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
		temp,w,x2,xi,xi2,xmu,xmu2;
	int i,isign,l,nl;

	if (x <= 0.0 || xnu < 0.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bad arguments in bessjy");
	nl=(x < CF_XMIN ? int(xnu+0.5) : MAX(0,int(xnu-x+1.5)));
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	w=xi2/CF_PI;
	isign=1;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=0;i<MAXIT;i++) {
		b += xi2;
		d=b-d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b-1.0/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=c*d;
		h=del*h;
		if (d < 0.0) isign = -isign;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i >= MAXIT)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"x too large in bessjy; try asymptotic expansion");
	rjl=isign*FPMIN;
	rjpl=h*rjl;
	rjl1=rjl;
	rjp1=rjpl;
	fact=xnu*xi;
	for (l=nl-1;l>=0;l--) {
		rjtemp=fact*rjl+rjpl;
		fact -= xi;
		rjpl=fact*rjtemp-rjl;
		rjl=rjtemp;
	}
	if (rjl == 0.0) rjl=EPS;
	f=rjpl/rjl;
	if (x < CF_XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,gam1,gam2,gampl,gammi);
		ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
		e=exp(e);
		p=e/(gampl*PI);
		q=1.0/(e*PI*gammi);
		pimu2=0.5*pimu;
		fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
		r=PI*pimu2*fact3*fact3;
		c=1.0;
		d = -x2*x2;
		sum=ff+r*q;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*(ff+r*q);
			sum += del;
			del1=c*p-i*del;
			sum1 += del1;
			if (fabs(del) < (1.0+fabs(sum))*EPS) break;
		}
		if (i > MAXIT)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bessy series failed to converge");
		rymu = -sum;
		ry1 = -sum1*xi2;
		rymup=xmu*xi*rymu-ry1;
		rjmu=w/(rymup-f*rymu);
	} else {
		a=0.25-xmu2;
		p = -0.5*xi;
		q=1.0;
		br=2.0*x;
		bi=2.0;
		fact=a*xi/(p*p+q*q);
		cr=br+q*fact;
		ci=bi+p*fact;
		den=br*br+bi*bi;
		dr=br/den;
		di = -bi/den;
		dlr=cr*dr-ci*di;
		dli=cr*di+ci*dr;
		temp=p*dlr-q*dli;
		q=p*dli+q*dlr;
		p=temp;
		for (i=1;i<MAXIT;i++) {
			a += 2*i;
			bi += 2.0;
			dr=a*dr+br;
			di=a*di+bi;
			if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
			fact=a/(cr*cr+ci*ci);
			cr=br+cr*fact;
			ci=bi-ci*fact;
			if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
			den=dr*dr+di*di;
			dr /= den;
			di /= -den;
			dlr=cr*dr-ci*di;
			dli=cr*di+ci*dr;
			temp=p*dlr-q*dli;
			q=p*dli+q*dlr;
			p=temp;
			if (fabs(dlr-1.0)+fabs(dli) <= EPS) break;
		}
		if (i >= MAXIT) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"cf2 failed in bessjy");
		gam=(p-f)/q;
		rjmu=sqrt(w/((p-f)*gam+q));
		rjmu=SIGN(rjmu,rjl);
		rymu=rjmu*gam;
		rymup=rymu*(p+q/gam);
		ry1=xmu*xi*rymu-rymup;
	}
	fact=rjmu/rjl;
	rj=rjl1*fact;
	rjp=rjp1*fact;
	for (i=1;i<=nl;i++) {
		rytemp=(xmu+i)*xi2*ry1-rymu;
		rymu=ry1;
		ry1=rytemp;
	}
	ry=rymu;
	ryp=xnu*xi*rymu-ry1;
}

double bessk0(const double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return ans;
}


double bessk1(const double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return ans;
}

double bessel_integer_K(const int n, const double x)
{
	int j;
	double bk,bkm,bkp,tox;

	if (n ==0) return bessk0(x) ;
	if (n ==1) return bessk1(x) ;
	tox=2.0/x;
	bkm=bessk0(x);
	bk=bessk1(x);
	for (j=1;j<n;j++) {
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}


double bessy0(const double x)
{
	double z,xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
			+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438
			+y*(47447.26470+y*(226.1030244+y*1.0))));
		ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}


double bessy1(const double x)
{
	double z,xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13
			+y*(-0.5153438139e11+y*(0.7349264551e9
			+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.4244419664e12
			+y*(0.3733650367e10+y*(0.2245904002e8
			+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

double bessel_integer_Y(const int n, const double x)
{
	int j;
	double by,bym,byp,tox;

	if (n ==0) return bessy0(x);
	if (n ==1) return bessy1(x);
	tox=2.0/x;
	by=bessy1(x);
	bym=bessy0(x);
	for (j=1;j<n;j++) {
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}


double bessel_fractional_I(double nu, double x)
{
	double I,Ip,K,Kp;
	bessik(x,nu,I,K,Ip,Kp);
	return I;
}

double bessel_fractional_K(double nu, double x)
{
	double I,Ip,K,Kp;
	bessik(x,nu,I,K,Ip,Kp);
	return K;
}

double bessel_fractional_I_Derx(double nu, double x)
{
	double I,Ip,K,Kp;
	bessik(x,nu,I,K,Ip,Kp);
	return Ip;
}

double bessel_fractional_K_Derx(double nu, double x)
{
	double I,Ip,K,Kp;
	bessik(x,nu,I,K,Ip,Kp);
	return Kp;
}



double bessel_fractional_J(double nu, double x)
{
	double J,Jp,Y,Yp;
	bessjy(x,nu,J,Y,Jp,Yp);
	return J;
}

double bessel_fractional_Y(double nu, double x)
{
	double J,Jp,Y,Yp;
	bessjy(x,nu,J,Y,Jp,Yp);
	return Y;
}

double bessel_fractional_J_Derx(double nu, double x)
{
	double J,Jp,Y,Yp;
	bessjy(x,nu,J,Y,Jp,Yp);
	return Jp;
}

double bessel_fractional_Y_Derx(double nu, double x)
{
	double J,Jp,Y,Yp;
	bessjy(x,nu,J,Y,Jp,Yp);
	return Yp;
}

double Hankel_Symbol(double alpha,double n)
{
	if(n>alpha) return 0;		/// pour l'instant afin de preserver le calcul de bessel_fractional_I_LargeZ_Log
	long_double r=gammaLD(0.5+alpha+n)/gammaLD(0.5+alpha-n);
	return r.todouble()/factorial(n);
}


double bessel_fractional_I_LargeZ_Log(double nu,double z, int n)
{
	double sum=1.;
	int signi=1;
	for(int i=1;i<n;i++)
	{
		signi=-signi;
		sum+=Hankel_Symbol(nu,i)*signi/pow(2.*z,i);
	}
	return z-log(z)/2.-ARM_NumericConstants::ARM_HALF_LOG_2PI+log(sum);
}

CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/