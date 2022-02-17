#ifndef __SCIDCONVOLUTION_HPP
#define __SCIDCONVOLUTION_HPP


#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

// THERE UNTIL ANTOINE AND SEB'S SCIDconvolution IS READY...

inline double PayoffTranche(double x, double L, double U)
{
	if (x<L) return 0;
	if (x<U) return (x-L);
	return U-L;
};

class MARKET_DLL SCIDconvolution
{
public:
	SCIDconvolution(){};
	SCIDconvolution(int nbNames,
				array<double> kmin,
				array<double> kmax,
				double paraSize = 1.0) { setParameters(nbNames, kmin, kmax, paraSize); };
	void setParameters(int nbNames,
				array<double> kmin,
				array<double> kmax,
				double paraSize = 1.0);

	void setPast(double previousRealizedLoss, double previousNotionalLoss);
	void setLGD(double *recovery, double *notional, double rho = 1);
	int getNbTranches() { return m_kmin.size(); }
	double computeTEL(int method, double *survival, double *ETL, double weight = 1.0);  // anything possible 
	double computeTELfuture(int method, 
							double *survival, 
							double *defaultTime, 
							double time,  // name is alive if defaultTime>time...
							double *ETL,
							double weight = 1.0);  

	double GaussianCall(double K, double mean, double var);
	inline double BetaCall(double K,      /* (I) lower strike   */
					  double mean,   /* (I) expected loss  */
					  double var);   /* (I) loss variance  */

// SCIDconvolution algorithm
private:
	array<double> m_probas, m_kmin, m_kmax;
	array<double> m_lgd, m_notional, m_recovery;
	double m_meanLGD;
	int m_nbNames, m_properConvolutionSize;
	double m_previousRealizedLoss, m_previousNotionalLoss;

	double PriceWithConvolutionDefaultedNames(double *survival, double averageLoss, double *ETL, double weight=1.0);
	double PriceWithBinomial(double survival, double averageLoss, double *ETL, double weight=1.0);     
	double PriceWithGaussian(double mean, double variance, double *ETL, double weight=1.0);
	double PriceWithBeta(double mean, double variance, double *ETL, double weight=1.0);
	
	double average(double *x);
};


DRLIB_END_NAMESPACE

#endif
