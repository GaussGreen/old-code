#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/SCIDconvolution.hpp"

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////// a few SCIDconvolution algorithms and approximations //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void SCIDconvolution::setParameters(int nbNames, array<double> kmin, array<double> kmax, double paraSize)
{
	QLIB_VERIFY(nbNames>0,"no names");
	QLIB_VERIFY(kmin.size()==kmax.size(),"not the same number of lower and upper attachment points");
	for (int i=0; i<kmin.size(); i++) 
		QLIB_VERIFY(kmin[i]<kmax[i], "lower attachment point greater than upper attachment points");

	m_nbNames = nbNames;
	m_kmin = kmin;
	m_kmax = kmax;
	m_probas.resize(nbNames+1);
	m_lgd.resize(nbNames);
	m_notional.resize(nbNames);
	m_recovery.resize(nbNames);
	m_properConvolutionSize = double(nbNames) + floor(paraSize*nbNames);
	setPast(0,0);
}


double SCIDconvolution::average(double *x)
{
	double mean = 0;
	for (int m=0; m<m_nbNames; m++) mean += x[m];
	return mean / double(m_nbNames);
}

void SCIDconvolution::setPast(double previousRealizedLoss, double previousNotionalLoss)
{
	m_previousRealizedLoss = previousRealizedLoss;
	m_previousNotionalLoss = previousNotionalLoss;
}

void SCIDconvolution::setLGD(double *recovery, double *notional, double rho)
{
	int m;
	m_meanLGD = 0;

	for (m=0; m<m_nbNames; m++) 
	{
		m_notional[m] = notional[m];
		m_recovery[m] = recovery[m];
		m_lgd[m] = notional[m]*(1-recovery[m]*rho);
		m_meanLGD += m_lgd[m];
	}
	m_meanLGD /= double(m_nbNames);
}

/*double SCIDconvolution::PriceWithFullConvolution(double *survival, double *ETL)
{
	double lossUnit;
	array<double> out(m_properConvolutionSize+1);
	ProperConvolution(
		m_nbNames,
		&m_notional[0],
		survival,
		&m_recovery[0],
		m_properConvolutionSize,
		0,
		&out[0],
		&lossUnit);
	double cdf=0;
	double loss = 0, mean = 0;
	for (int m=0; (m<=m_properConvolutionSize)&&(cdf<1-1e-9); m++)
	{
		mean += out[m]*loss;
		cdf+=out[m];
		for (int k=0; k<m_kmin.size(); k++)
			ETL[k] += out[m]*PayoffTranche(m_previousRealizedLoss + loss, m_kmin[k], m_kmax[k]);
		loss += lossUnit;
	}
	return mean;
}
*/

double SCIDconvolution::PriceWithConvolutionDefaultedNames(double *survival, double averageLoss, double *ETL,
														   double weight)
{
	int m,n;
	m_probas[0]=1.0;
	for (m=0; m<m_nbNames; m++) m_probas[m+1]=0;

	double mean = 0;
	for (m=0; m<m_nbNames; m++)
	{
		for (n=m+1; n>=1; n--) m_probas[n]=m_probas[n]*survival[m]+m_probas[n-1]*(1-survival[m]);
		m_probas[0]*=survival[m];
		mean += (1-survival[m])*averageLoss;
	}
	double cdf = 0, loss = 0;
	for (m=0; (m<m_nbNames)&&(cdf-m_probas[m]<1-1e-12); m++)
	{
		cdf+=m_probas[m];
		for (int k=0; k<m_kmin.size(); k++)
			ETL[k] += weight*m_probas[m]*PayoffTranche(m_previousRealizedLoss + loss, m_kmin[k], m_kmax[k]);
		loss += averageLoss;
	}
	return mean;
}

double SCIDconvolution::PriceWithBinomial(double survival, double averageLoss, double *ETL,
										  double weight)
{
	double pdf = pow(survival, m_nbNames);
	double cdf = pdf;
	double loss = 0.0;
	for (int m=0; (m<m_nbNames)&&(cdf-pdf<1-1e-12); m++)
	{
		for (int k=0; k<m_kmin.size(); k++)
			ETL[k] += weight*pdf*PayoffTranche(m_previousRealizedLoss+ (1-m_previousNotionalLoss)*m*averageLoss, m_kmin[k], m_kmax[k]);
		pdf *= (1.0 - survival)*(m_nbNames - m) / (survival*(m+1));
		cdf+= pdf;
		loss += averageLoss;
	}
	return (1-survival)*m_nbNames;
}


double SCIDconvolution::PriceWithGaussian(double mean, double variance, double *ETL,
										  double weight)
{
	if ((variance>0)&&(mean>0))
	{
		double inv = 1/(1 - m_previousNotionalLoss);
		mean *= inv;
		variance *= inv*inv;
		for (int k=0; k<m_kmin.size(); k++) 
			ETL[k] += weight*(GaussianCall(inv*(m_kmin[k]-m_previousRealizedLoss),mean,variance)
						- GaussianCall(inv*(m_kmax[k]-m_previousRealizedLoss),mean,variance))
						* (1 - m_previousNotionalLoss);
	}
	else
		for (int k=0; k<m_kmin.size(); k++)
			ETL[k] += weight*PayoffTranche(m_previousRealizedLoss+ (1-m_previousNotionalLoss)*mean, m_kmin[k], m_kmax[k]);

	return mean;
//	static const sqrt2pi=2.5066282746310004990396590352614;
//	if ((variance>0)&&(mean>0))
//	{
//		double stdev = sqrt(variance);
//		for (int k=0; k<m_kmin.size(); k++) 
//			ETL[k] += m_kmax[k]-m_kmin[k]-(m_kmax[k]-mean)*GtoNormalCum((m_kmax[k]-mean)/stdev)
//				          +(m_kmin[k]-mean)*GtoNormalCum((m_kmin[k]-mean)/stdev)
//						  + stdev /sqrt2pi*(exp(-0.5*(m_kmin[k]-mean)*(m_kmin[k]-mean)/variance) 
//						  - exp(-0.5*(m_kmax[k]-mean)*(m_kmax[k]-mean)/variance));
//	}
//	return mean;
}

double SCIDconvolution::PriceWithBeta(double mean, double variance, double *ETL,
									  double weight)
{
	if ((variance>0)&&(mean>0))
	{
		double inv = 1/(1 - m_previousNotionalLoss);
		mean *= inv;
		variance *= inv*inv;
		for (int k=0; k<m_kmin.size(); k++) 
			ETL[k] += weight*(BetaCall(inv*(m_kmin[k]-m_previousRealizedLoss),mean,variance)
						- BetaCall(inv*(m_kmax[k]-m_previousRealizedLoss),mean,variance))
						* (1 - m_previousNotionalLoss);
	}
	else
		for (int k=0; k<m_kmin.size(); k++)
			ETL[k] += weight*PayoffTranche(m_previousRealizedLoss+ (1-m_previousNotionalLoss)*mean, m_kmin[k], m_kmax[k]);

	return mean;
}


double SCIDconvolution::computeTEL(int method, double *survival, double *ETL, double weight)
{
	QLIB_VERIFY(method>0 && method<5, "unknown convolution method");
	double averageSP, mean, variance;
	int m;
	switch (method)
	{
	case 1: 
		return PriceWithConvolutionDefaultedNames(survival, m_meanLGD, ETL, weight);
		break;
	case 2:
		averageSP = average(survival);
		return PriceWithBinomial(averageSP, m_meanLGD, ETL, weight);
		break;
	case 3:
		mean = 0;
		variance = 0;
		for (m=0; m<m_nbNames; m++) 
		{
			mean += (1-survival[m])*m_lgd[m];
			variance += survival[m]*(1-survival[m])*m_lgd[m]*m_lgd[m];
		}
		return PriceWithGaussian(mean,variance,ETL, weight);
		break;
	case 4:
		mean = 0;
		variance = 0;
		for (m=0; m<m_nbNames; m++) 
		{
			mean += (1-survival[m])*m_lgd[m];
			variance += survival[m]*(1-survival[m])*m_lgd[m]*m_lgd[m];
		}
		return PriceWithBeta(mean, variance, ETL, weight);
		break;
	default: 
		return 0;
		break;
	}
}

 double SCIDconvolution::computeTELfuture(
							int method,
							double *survival, 
							double *defaultTime, double time,  // name is alive if defaultTime>time...
							double *ETL,
							double weight)  // use the Beta Approximation...
{
	if (method==4)
	{
		double mean = 0;
		double variance = 0;
		for (int m=0; m<m_nbNames; m++) 
			if (defaultTime[m]>time)
			{
				mean += (1-survival[m])*m_lgd[m];
				variance += survival[m]*(1-survival[m])*m_lgd[m]*m_lgd[m];
			}
		return PriceWithBeta(mean,variance,&ETL[0]);
	}
	else 
	{
		for (int m=0; m<m_nbNames; m++)
			if (defaultTime[m]<=time) survival[m]=1;
//		if (method==1)
			return PriceWithConvolutionDefaultedNames(survival, m_meanLGD, ETL);
	//	else 
	//		return PriceWithFullConvolution(survival, ETL);
	}
}


/* incomplete beta function continued fraction expansion */
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
static double betacf(double a, double b, double x)
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
    if (m > MAXIT) return h;  // should throw exception
//        throw ModelException(routine, "a or b too big, or MAXIT too small");

    return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

/* incomplete beta function */
static double betai(double a, double b, double x)
{
    double bt,y;
    double out;

	x = Maths::max(0.0,x);
	x = Maths::min(x,1.0);

    if (x == 0.0 || x == 1.0) 
        bt=0.0;
    else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));

    if (x < (a+1.0)/(a+b+2.0)) 
    {
        y = betacf(a,b,x);
        out = bt * y / a;
    }
    else
    {
        y = betacf(b,a,1.0-x);
        out = 1.0 - bt * y / b;
    }
    return out;
}


/**
 * Calculate E((L-K)^+)
 * Assuming that L follows a beta distribution
 * with a given expected loss and variance.
 *
 * we use an "incomplete beta" integ close formula 
 */
double SCIDconvolution::BetaCall(
    double K,      /* (I) lower strike   */
    double mean,   /* (I) expected loss  */
    double var)    /* (I) loss variance  */
{
    double a,b,x, y1,y2;                          
    double e;
    /* strike is >= 1 */
    if (1-K<1e-12)
        return 0.0;
    /* strike is <= 0 */
    if (K<1e-12)
        return mean - K;
    /* variance is 0 */
    if (var<1e-12)
        return Maths::max(mean-K,0.);
    x = mean*(1.-mean)/var - 1.; 
    if ((mean<0) || (mean>1)) return 0;
    if (x<1e-12)
        return (1-K)*mean;
    /* general case */
    a = mean*x;
    b = (1.-mean)*x;
    y1 = betai(b,a+1.,1.-K);
    y2 = betai(b,a   ,1.-K);

    e = a/(a+b) * y1 - K * y2;
    e = Maths::max(0.,e);
    e = Maths::min(1.-K, e);
    return e;
}

double SCIDconvolution::GaussianCall(
    double K,      /* (I) lower strike   */
    double mean,   /* (I) expected loss  */
    double var)    /* (I) loss variance  */
{
	static double const sqrt2pi=2.5066282746310004990396590352614;
	if (var>0)
	{
		double stdDev = sqrt(var);
		double newStrike = (K-mean)/stdDev;
		return stdDev*(exp(-newStrike*newStrike*0.5)/sqrt2pi - newStrike*(1-N1(newStrike)) );
	}
	return Maths::max(mean-K,0.0);

}

DRLIB_END_NAMESPACE
