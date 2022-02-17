/****************************************************************************************
  Computes the price of a Time Swap in the context of a Resetable Time Swap by optimizing 
the range using a Newton algorithm 
******************************************************************************************/







#include   "utallhdr.h"
#include   <OPFNCTNS.H>
#include   <num_h_gamma.h"
#include   <srt_h_resetable.h"
#include	<math.h"


#define pi 3.14159265358979
#define ITMAX 10
#define ZEPS 1.0e-10
/* Here ITMAX is the maximum allowed number of iterations; ZEPS is
a small number that protects against trying to achieve fractional accuracy for a minimum that
happens to be excatly zero */
#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f);

/* Function which obtains the minimum value of a function of 1 variable when the derivative 
is given. (taken from Numerical recipes in C)*/

double Newdbrent(double ax,
			   double bx,
			   double cx,
			   double (*f)(double,double*,int,double*,double*,double **,double),
			   double (*df)(double,double*,int,double*,double*,double **,double),
			   double tol,
			   double *xmin,
			   double *DRS,
			   int nmat,
			   double *MatShiftes,
			   double  *isFriday,
			   double **param,
			   double width)
{
int iter,ok1,ok2;
double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

a=(ax<cx ? ax:cx);
b=(ax>cx ? ax:cx);
x=w=v=bx;
fw=fv=fx=(*f)(x,DRS,nmat,MatShiftes,isFriday,param,width);
dw=dv=dx=(*df)(x,DRS,nmat,MatShiftes,isFriday,param,width);
for (iter=1;iter<=ITMAX;iter++) {
	xm=0.5*(a+b);
	tol1=tol*fabs(x)+ZEPS;
	tol2=2.0*tol1;
	if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
		*xmin=x;
		return fx;
	}
	if (fabs(e)>tol1) {
		d1=2.0*(b-a);
		d2=d1;
		if (dw != dx) d1=(w-x)*dx/(dx-dw);
		if (dv != dx) d2 = (v-x)*dx/(dx-dv);
		u1=x+d1;
		u2=x+d2;
		ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
		ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
		olde=e;
		e=d;
		if (ok1 || ok2) {
			if (ok1 && ok2) 
				d=(fabs(d1) < fabs(d2) ? d1 : d2);
			else if (ok1)
				d=d1;
			else
				d=d2;
			if (fabs(d) <= fabs(0.5*olde)) {
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			} else{
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			
		}
	} else {
		d=0.5*(e=(dx>=0.0 ? a-x:b-x));
	}
	if (fabs(d) >= tol1) {
		u=x+d;
		fu=(*f)(u,DRS,nmat,MatShiftes,isFriday,param,width);
	} else {
		u=x+SIGN(tol1,d);
		fu=(*f)(u,DRS,nmat,MatShiftes,isFriday,param,width);
		if (fu>fx) {
			*xmin=x;
			return fx;
		}
	}
	du=(*df)(u,DRS,nmat,MatShiftes,isFriday,param,width);
	if (fu <= fx) {
		if (u>=x) a=x; else b=x;
		MOV3(v,fv,dv,w,fw,dw)
		MOV3(w,fw,dw,x,fx,dx)
		MOV3(x,fx,dx,u,fu,du)
	} else {
		if (u<x) a=u; else b=u;
		if (fu <= fw || w==x) {
			MOV3(v,fv,dv,w,fw,dw)
			MOV3(w,fw,dw,u,fu,du)
		} else if (fu<fv || v==x || x==w) {
			MOV3(v,fv,dv,u,fu,du)
		}
		}
}
return fx;
}


/* Price of a call spread option which approximates the price of an option whose payoff is
1 if at maturity, the forward is in the range [K-width+eps,K+eps]*/


double PrixDoubleCallSpreadMerton(double forward,
				 double *param, /* 5 parameters for the dynamics of the forward: sig,U1,lambda1,U2,lambda2  */ 
				 double mat,
				 double K,
				 double width, 
				 double callspread,
				 double eps)
	 		
{
double strike1;
double strike2;
double spread;	
double lambda1p,lambda1;
double lambda2p,lambda2;
double N1max;
double N2max;
double d1up,d1upsp,d1dn,d1dnsp;
double d2up,d2upsp,d2dn,d2dnsp;
int n1,n2;
double premium1up=0;
double  premium2up=0;
double premium1upsp=0;
double  premium2upsp=0;
double premium1dn=0;
double  premium2dn=0;
double premium1dnsp=0;
double  premium2dnsp=0;
double premium;
double proba1, proba1p,proba2,proba2p;
double probainf;
double U1,U2,sigma,T,dFwd;

sigma = param[1];
U1 = param[2];
lambda1 = param[3];
U2 = param[4];
lambda2 = param[5];
T=mat;
dFwd = forward;

spread = callspread;
	
strike2=K-width+eps;
strike1=K+eps;

probainf = 0.000001;

lambda1p = lambda1*(1+U1);
lambda2p = lambda2*(1+U2);


N1max = lambda1p*T+10*sqrt(lambda1p*T);
N2max = lambda2*T+10*sqrt(lambda2*T);

for (n1 = 0;n1<N1max;n1++)
{
	for (n2 = 0;n2<N2max;n2++)
	{


		
		proba1p= exp(-lambda1p*T)*pow(lambda1p*T,n1)/fact(n1);
		proba2p = exp(-lambda2p*T)*pow(lambda2p*T,n2)/fact(n2);
		
		proba1= exp(-lambda1*T)*pow(lambda1*T,n1)/fact(n1);
		proba2 = exp(-lambda2*T)*pow(lambda2*T,n2)/fact(n2);
 
		if((proba1p*proba2p>probainf)||(proba1*proba2>probainf))
		{

			d1up = (log(dFwd/strike1)-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1upsp = (log(dFwd/(strike1+spread))-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1dn = (log(dFwd/strike2)-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1dnsp = (log(dFwd/(strike2-spread))-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));

			d2up = d1up-sigma*sqrt(T);
			d2upsp = d1upsp-sigma*sqrt(T);
			d2dn = d1dn-sigma*sqrt(T);
			d2dnsp = d1dnsp-sigma*sqrt(T);

			premium1up +=proba1p*proba2p*norm(d1up);
			premium2up +=proba1*proba2*norm(d2up);

			premium1upsp +=proba1p*proba2p*norm(d1upsp);
			premium2upsp +=proba1*proba2*norm(d2upsp);


			premium1dn +=proba1p*proba2p*norm(d1dn);
			premium2dn +=proba1*proba2*norm(d2dn);

			premium1dnsp +=proba1p*proba2p*norm(d1dnsp);
			premium2dnsp +=proba1*proba2*norm(d2dnsp);


			if (strike1 <= 0 || strike1 <= -spread || strike2 <= 0|| strike2 <= spread ||
				dFwd <= 0|| U1 <= -1||U2 <= -1 || sigma == 0)
			{
				smessage("divide by 0");
			}
		}

	}
}

premium = (dFwd*premium1dnsp-(strike2-spread)*premium2dnsp)-(dFwd*premium1dn-strike2*premium2dn)
+(dFwd*premium1upsp-(strike1+spread)*premium2upsp)-(dFwd*premium1up-strike1*premium2up);
premium = premium/spread;

return(premium);
}


/* Exact price of the option which is only approximated above, basically a digital option on 
a range */

/* We use the Call spread as it provides more smoothness in the derivative and thus ensures 
faster convergence of the Newton Algorithm.
*/

double PrixDoubleDigitMerton(double forward,
				 double *param, 
				 double mat,
				 double K,
				 double width, 
				 double eps)
	 		
{
int N1max,N2max; /* Maximum number of jumps*/
int N1,N2;
double strike1,strike2;
double d1,d2;
double probasauts;
double probainf;  /* Threshold according to  which we consider that more jumps are very unlikely */
double prix;

N1max = (int) (param[3]*mat+10*sqrt(param[3]*mat));
N2max = (int) (param[5]*mat+10*sqrt(param[5]*mat));
probainf=0.000001;
prix=0;
for (N1=0;N1<=N1max;N1++)
for (N2=0;N2<=N2max;N2++)
{
	probasauts = exp(-param[3]*mat-param[5]*mat+N1*log(param[3]*mat)+
		N2*log(param[5]*mat))/(fact(N1)*fact(N2));

	if (probasauts>probainf)
	{
	strike1=K-width+eps;
	strike2=K+eps;
	d1=(param[3]*param[2]*mat+param[5]*param[4]*mat+0.5*param[1]
		*param[1]*mat+log(strike1/(forward*exp(N1*log(1+param[2])+N2
		*log(1+param[4])))))/(param[1]*sqrt(mat));

	d2=(param[3]*param[2]*mat+param[5]*param[4]*mat+0.5*param[1]
		*param[1]*mat+log(strike2/(forward*exp(N1*log(1+param[2])+N2
		*log(1+param[4])))))/(param[1]*sqrt(mat));
		prix+=probasauts*(norm(d2)-norm(d1));
	}
}
return(prix);
}



/*Computes the derivative of the price coming from the Call spreads*/

double DerivPrixDoubleCallSpreadMerton(double forward,
				 double *param,  
				 double mat,
				 double K,
				 double width, 
				 double eps)
	 		
{
double strike1;
double strike2;
double spread;	
double lambda1p,lambda1;
double lambda2p,lambda2;
double N1max;
double N2max;
double d1up,d1upsp,d1dn,d1dnsp;
double d2up,d2upsp,d2dn,d2dnsp;
int n1,n2;
double premium1up=0;
double premium2up=0;
double premium3up=0;
double premium1upsp=0;
double premium2upsp=0;
double premium3upsp=0;
double premium1dn=0;
double premium2dn=0;
double premium3dn =0; 
double premium1dnsp=0;
double premium2dnsp=0;
double premium3dnsp =0;
double deriv;
double proba1, proba1p,proba2,proba2p;
double probainf;
double U1,U2,sigma,T,dFwd;

sigma = param[1];
U1 = param[2];
lambda1 = param[3];
U2 = param[4];
lambda2 = param[5];
T = mat;
dFwd = forward;

spread = 0.00125;
	
strike2=K-width+eps;
strike1=K+eps;

probainf = 0.000001;

lambda1p = lambda1*(1+U1);
lambda2p = lambda2*(1+U2);


N1max = lambda1p*T+10*sqrt(lambda1p*T);
N2max = lambda2*T+10*sqrt(lambda2*T);

for (n1 = 0;n1<N1max;n1++)
{
	for (n2 = 0;n2<N2max;n2++)
	{


		
		proba1p= exp(-lambda1p*T)*pow(lambda1p*T,n1)/fact(n1);
		proba2p = exp(-lambda2p*T)*pow(lambda2p*T,n2)/fact(n2);
		
		proba1= exp(-lambda1*T)*pow(lambda1*T,n1)/fact(n1);
		proba2 = exp(-lambda2*T)*pow(lambda2*T,n2)/fact(n2);

		if((proba1p*proba2p>probainf)||(proba1*proba2>probainf))
		{

			d1up = (log(dFwd/strike1)-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1upsp = (log(dFwd/(strike1+spread))-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1dn = (log(dFwd/strike2)-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));
			d1dnsp = (log(dFwd/(strike2-spread))-U1*lambda1*T-U2*lambda2*T+n1*log(1+U1)+
					n2*log(1+U2)+sigma*sigma*T/2)/(sigma*sqrt(T));

			d2up = d1up-sigma*sqrt(T);
			d2upsp = d1upsp-sigma*sqrt(T);
			d2dn = d1dn-sigma*sqrt(T);
			d2dnsp = d1dnsp-sigma*sqrt(T);

			premium1up +=proba1p*proba2p*(-1/(strike1*sigma*sqrt(T)))*exp(-d1up*d1up/2)/sqrt(2*SRT_PI);
			premium2up +=proba1*proba2*(-1/(strike1*sigma*sqrt(T)))*exp(-d2up*d2up/2)/sqrt(2*SRT_PI);
			premium3up +=proba1*proba2*norm(d2up);



			premium1upsp +=proba1p*proba2p*(-1/((strike1+spread)*sigma*sqrt(T)))*exp(-d1upsp*d1upsp/2)/sqrt(2*SRT_PI);
			premium2upsp +=proba1*proba2*(-1/((strike1+spread)*sigma*sqrt(T)))*exp(-d2upsp*d2upsp/2)/sqrt(2*SRT_PI);
			premium3upsp +=proba1*proba2*norm(d2upsp);	

			premium1dn +=proba1p*proba2p*(-1/(strike2*sigma*sqrt(T)))*exp(-d1dn*d1dn/2)/sqrt(2*SRT_PI);
			premium2dn +=proba1*proba2*(-1/(strike2*sigma*sqrt(T)))*exp(-d2dn*d2dn/2)/sqrt(2*SRT_PI);
			premium3dn +=proba1*proba2*norm(d2dn);

			premium1dnsp +=proba1p*proba2p*(-1/((strike2-spread)*sigma*sqrt(T)))*exp(-d1dnsp*d1dnsp/2)/sqrt(2*SRT_PI);
			premium2dnsp +=proba1*proba2*(-1/((strike2-spread)*sigma*sqrt(T)))*exp(-d2dnsp*d2dnsp/2)/sqrt(2*SRT_PI);
			premium3dnsp +=proba1*proba2*norm(d2dnsp);
		
		if (strike1 <= 0 || strike1 <= -spread || strike2 <= 0|| strike2 <= spread ||
				dFwd <= 0|| U1 <= -1||U2 <= -1 || sigma == 0)
			{
				smessage("divide by 0");
			}
		
		}

	}
}

deriv = (dFwd*premium1dnsp-(strike2-spread)*premium2dnsp-premium3dnsp)-(dFwd*premium1dn-strike2*premium2dn-premium3dn)
+(dFwd*premium1upsp-(strike1+spread)*premium2upsp-premium3upsp)-(dFwd*premium1up-strike1*premium2up-premium3up);

deriv = deriv/spread;


return(deriv);
}


/* Derivative of the digital price */

double DerivPrixDoubleDigitMerton(double forward,
				 double *param, 
				 double mat,
				 double K,
				 double width, 
				 double eps)
				
{
int N1max,N2max; 
int N1,N2;
double strike1,strike2;
double d1,d2;
double probasauts;
double probainf;  
double deriv;

N1max = (int) (param[3]*mat+10*sqrt(param[3]*mat));
N2max = (int) (param[5]*mat+10*sqrt(param[5]*mat));
probainf=0.000001;
deriv=0;

for (N1=0;N1<=N1max;N1++)
for (N2=0;N2<=N2max;N2++)
{
	probasauts = exp(-param[3]*mat-param[5]*mat+N1*log(param[3]*mat)+
		N2*log(param[5]*mat))/(fact(N1)*fact(N2));
	if (probasauts>probainf)
	{
	strike1=K-width+eps;
	strike2=K+eps;
	d1=(param[3]*param[2]*mat+param[5]*param[4]*mat+0.5*param[1]
			*param[1]*mat+log(strike1/(forward*exp(N1*log(1+param[2])+N2
			*log(1+param[4])))))/(param[1]*sqrt(mat));
	d2=(param[3]*param[2]*mat+param[5]*param[4]*mat+0.5*param[1]
			*param[1]*mat+log(strike2/(forward*exp(N1*log(1+param[2])+N2	
			*log(1+param[4])))))/(param[1]*sqrt(mat));
	
	deriv+=probasauts/(param[1]*sqrt(mat))*(exp(-0.5*d2*d2)/strike2-
			exp(-0.5*d1*d1)/strike1)/sqrt(2*pi);
	}
}
return(deriv);
}

/* Computes the price of a Time Swap with a fixed range*/
double PrixTimeSwap(
			double eps,
			double *DRS,    /* Term structure of the DRS. tab[1..nmat] */
			int nmat,    /* Number of dates +1 (the first date of the period)*/
			double *mat,   /* tab[1..nmat] */
			double *IsFriday,  /*tab[1..nmat] */
			double **param,  /*tab[1..nmat][1..5] *. term structure of the jump parameters */
			double width,     /* width of the range */
			double callspread			)
{
double prix;
int i;

prix=0;
for (i=2;i<=nmat;i++)  /* Computes the approximate price of the range digital everyday, 
					   making sure that friday counts as 3 days*/
	prix+=PrixDoubleCallSpreadMerton(
		DRS[i],
		param[i],
		mat[i],
		DRS[1],
		width,
		callspread,
		eps)*(1+2*IsFriday[i]);
return prix;
}

/* Computes the derivative of the price of the Time Swap computed above */	
double DerivPrixTimeSwap(
			double eps,
			double *DRS,
			int nmat,
			double *mat,
			double *IsFriday,
			double **param,
			double width)
{
double deriv;
int i;

deriv=0;
for (i=2;i<=nmat;i++)
	deriv+=DerivPrixDoubleCallSpreadMerton(
		DRS[i],
		param[i],
		mat[i],
		DRS[1],
		width,
		eps)*(1+2*IsFriday[i]);
return(-deriv);
}			


double Optimise_Time_Swap(int nmat,
						 double *DRS,
						 double width,
						 double lowermin,
						 double lowermax,
					     double *eps,
						 double callspread,
               			 double *mat,
					     double *IsFriday,
					     double **param) 
{
int iter = 0;
double pleft, pright, pmid;
double xleft, xright, xmid;
double xmidleft,xmidright;
double pmidleft,pmidright;
double precision;
int maxiter = 15;

xleft = lowermin;
xright = lowermax;
xmid = (lowermin+lowermax)/2;
precision = 1.0;
pleft = PrixTimeSwap(xleft,DRS, nmat, mat,IsFriday, param, width, callspread);
pright = PrixTimeSwap(xright,DRS, nmat, mat,IsFriday, param, width, callspread);



	while ((iter < maxiter)&&(precision>=0.0001))
	{

		pmid = PrixTimeSwap(xmid,DRS, nmat, mat,IsFriday, param, width, callspread);

		xmidleft = (xleft+xmid)/2;
		xmidright = (xright+xmid)/2;

		pmidleft = PrixTimeSwap(xmidleft,DRS, nmat, mat,IsFriday, param, width, callspread);
		pmidright = PrixTimeSwap(xmidright,DRS, nmat, mat,IsFriday, param, width, callspread);
		
		if ((pmid>pmidleft)&&(pmid>pmidright))
		{
			xleft = xmidleft;
			xright = xmidright;
		}
		else if ((pmidleft>pmid)&&(pmidleft>pmidright))
		{
			xright = xmid;
		}

		else
		{
			xleft = xmid;
		}

		precision = fabs(pmidleft-pmidright);
		xmid = (xleft+xright)/2;
		iter++;

		
	}

return pmid;

}

/* Function which returns the optimal range (eps) and the corresponding Time Swap price*/

Err NewtonJumps(int nmat, /* = 1+number of range digitals within the period */
				double *DRS, /* term structure of the  DRS at the  reset date
							 : [1..nmat] */
				double *maturites, /* [1..nmat] */
				double  *isFriday,  /* [1..nmat] vector assigning 1 to fridays and 0 to other days */
				double **param, /*parameters of the dynamics of each  DRS[1..nmat][1..5] */
				double width,
				double callspread,
				double *eps,   
				double *prix)
{
double *MatShiftes;
int i;
double lowermin,lowermax,prec;

lowermin= -width;    
lowermax= 2*width;
prec= 0.0001;       /* precision for the optimization of the range */

MatShiftes=dvector(1,nmat);
for (i=1;i<=nmat;i++) MatShiftes[i]=(maturites[i]-maturites[1])/365.0;
*prix=Optimise_Time_Swap(nmat,
						 DRS,
						 width,
						 lowermin,
						 lowermax,
					     eps,
						 callspread,
               			 MatShiftes,
					     isFriday,
					     param);  /*  Calls the  Newton algorithm*/
free_dvector(MatShiftes,1,nmat);

return NULL;

}
