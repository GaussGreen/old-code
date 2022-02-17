#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

#include <math.h>
#include <stdio.h>
#include "beta.h"
#include "error2.h"


/* ---------------------------------------------------------
// beta
// from NR
*/
double beta(double a, double b)
{
    return exp(gammln(a)+gammln(b)-gammln(a+b));
}

/* ---------------------------------------------------------
// betacf
// from NR
*/ 
double betacf(double a, double b, double x)
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) DR_Error("a or b too big, or MAXIT too small in betacf");
	return h;
}

/* ----------------------------------------------------
// gammln
// from NR
*/
double gammln(double xx)
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

/* --------------------------------------------------------
// betai
// from NR
*/
double betai(double a, double b, double x)
{
	double betacf(double a, double b, double x);
	double gammln(double xx);
	double bt;

	if (x < 0.0 || x > 1.0) DR_Error("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

/* --------------------------------------------------------
// gser
// from NR
*/
void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) DR_Error("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=MAXIT;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		DR_Error("a too large, MAXIT too small in routine gser");
		return;
	}
}

/* --------------------------------------------------------
// gcf
// from NR
*/
void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=MAXIT;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) DR_Error("a too large, MAXIT too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

/* --------------------------------------------------------
// gammp
// from NR
*/
double gammp(double a, double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) DR_Error("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

/* ---------------------------------------------------------
// Beta_Distribution
*/
double Beta_Distribution(double x, double a, double b)
{
    if(x==0||x==1) return 0.0; 
    return pow(1.- x,b-1.)*pow(x,a-1.)/beta(a,b);
}

/* ---------------------------------------------------------
// Beta_Mean
*/
double Beta_Mean(double alpha, double beta)
{
    return alpha / (alpha+beta);
}

/* ---------------------------------------------------------
// Beta_Variance
*/
double Beta_Variance(double alpha, double beta)
{
    double sum = alpha + beta;
    return alpha*beta / (sum*sum*(sum+1));
}

/* ---------------------------------------------------------
// Beta_Skew
*/
double Beta_Skew(double alpha, double beta)
{
    double sum = alpha + beta;
    return 2*(beta-alpha)*sqrt((1+sum)/(alpha*beta)) / (2+sum);
}

/* ---------------------------------------------------------
// Beta_Kurtosis
*/
double Beta_Kurtosis(double alpha, double beta)
{
    double sum = alpha + beta;
    double prod = alpha*beta;
    double res = alpha*alpha*alpha;
    res += alpha*alpha*(1.-2*beta);
    res += beta*beta*(1+beta);
    res -= 2*prod*(2+beta);
    res /= prod*(sum+2)*(sum+3);
    return res;
}

/* ---------------------------------------------------------
// Beta_Mode
*/
double Beta_Mode(double alpha, double beta)
{
    return (alpha-1.)/(alpha+beta-2.);
}

#undef MAXIT
#undef EPS
#undef FPMIN

