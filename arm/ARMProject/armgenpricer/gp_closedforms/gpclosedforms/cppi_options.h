/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file cppi_options.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_CPPI_OPTIONS_H
#define _GP_CF_CPPI_OPTIONS_H


#include "firsttoinc.h"
#include "gpbase/port.h"
#include <vector>
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpmodels/sabr_eq.h"

#include "expt.h"
using std::vector;

class ARM_ZeroCurve;
//class ARM_SABR_Eq;

CC_BEGIN_NAMESPACE(ARM)

/// density of the integrated volatility squared in a black and sholes setting

double IntegratedVolSquarredDensity(double x,double mub,double sigmab);

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Classes for building an equity binomial tree
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
inline double max(double u,double v)
{
	return ((u>=v)? u : v);
}

inline double min(double u,double v)
{
	return ((u>=v)? v : u);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTree
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

struct BinomialTree 
{
	int itsSize;
	double UpProbability;
	double u;
	double d;
	/// CapitalizationFactor is like exp(r dt)
	double CapitalizationFactor;
	double BorrowingSpreadingFactor;
	double T;
	int path_state_size;
	
private: 
	vector<double>* itsUnderlyingContainer;
	vector<double>* itsPathStateContainer;
	vector<double>* itsPathOptionContainer;
	
public:
	BinomialTree(int n, int size,double T0);
	
	~BinomialTree();
	
	inline double get_underlyingstate(int i,int j)
	{
		return (*itsUnderlyingContainer)[i*itsSize+j];
	}
	
	inline double get_pathstate(int i,int j,int k)
	{
		return (*itsPathStateContainer)[(i*itsSize+j)*path_state_size+k];
	}
	
	inline double get_pathoption(int i,int j,int k)
	{
		return (*itsPathOptionContainer)[(i*itsSize+j)*path_state_size+k];
	}
	
	
	
	inline void set_underlyingstate(int i,int j,double v)
	{
		///		if((i<0)||(i>=itsSize)||(j<0)||(j>=itsSize)) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"set_underlyingstate : bad indexes i or j");
		
		(*itsUnderlyingContainer)[i*itsSize+j]=v;
		return;
	}
	
	inline void set_pathstate(int i,int j,int k,double v)
	{
	/*		if((i<0)||(i>=itsSize)||(j<0)||(j>=itsSize)||(k<0)||(k>=path_state_size)) 	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"set_pathstate : bad indexes i or j");
	if((k>0)&&(get_pathstate(i,j,k-1)>v))
	{
	///__asm int 3;
	double previous=get_pathstate(i,j,k-1);
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"set_pathstate : non monotonicity from below");
	}
	if((k<path_state_size-1)&&(get_pathstate(i,j,k+1)<v))
	{
	///__asm int 3;
	double further=get_pathstate(i,j,k+1);
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"set_pathstate : non monotonicity from upper");
	}
		*/		
		(*itsPathStateContainer)[(i*itsSize+j)*path_state_size+k]=v;
		return;
	}
	
	inline void unsafe_set_pathstate(int i,int j,int k,double v)
	{
		(*itsPathStateContainer)[(i*itsSize+j)*path_state_size+k]=v;
		return;
	}
	inline void set_pathoption(int i,int j,int k,double v)
	{
		///		if((i<0)||(i>=itsSize)||(j<0)||(j>=itsSize)||(k<0)||(k>=path_state_size)) throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"set_pathoption : bad indexes i or j");
		(*itsPathOptionContainer)[(i*itsSize+j)*path_state_size+k]=v;
		return;
	}
	
	
	double Price();
	
	/// interpole_pathoption computes the interpolated value of the option, knowing the option values (i,j) with 
	/// a pathstate equal to the porfolio value after a move equal to v
	double interpole_pathoption(int i,int j,double v);
	
	virtual double PortfolioValueNextStep(double currenttime,double PortValue0, double S0, double S1,double D0,double D1)=0;
	
		
};

////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTreeMepiVanillaOption
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

struct BinomialTreeMepiVanillaOption : public BinomialTree
{
	double PATHSTATE_LIMITDOWN;
	double PATHSTATE_LIMITUP;
	double S;
	double K;
	double Vol;
	double r;
	double BorrowingSpread;
	double YearlyFees;
	double InitialPortfolio;
	double InitialMinimumGaranty;
	double MinimumGaranty;
	double MinimumExposure;
	double MaximumExposure;
	double RiskFactor;
	double VolRatio;
	double FeesPerPeriod;
	int AveragingPeriodNb;
	int CallOrPut;
	
	BinomialTreeMepiVanillaOption(int n0, int nz0,double p0,double K0,double Vol0,double r0,double BorrowingSpread0,double YearlyFees0,double T0,
			double minExp,double MaximumExposure0,double riskFactor,double g0,double g,double minport,double maxport,int AveragingPeriodNb0,int CallOrPut0);

	void VanillaOptionInit(double Vol0, double initialPortfolio);

///  the result of add_pathstate is that v should be between get_pathstate(i,j,0) and get_pathstate(i,j,path_state_size-1)
	void add_pathstate(int i,int j,double v);
	
	double PortfolioValueNextStep(double currenttime,double PortfolioPrecedingValue, double S0, double S1,double D0,double D1);
	
	double InitialExposition();
};


/// we derive a new classe that adds the stochasticity of the volatility

////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTreeMepiVanillaOptionSV
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

struct BinomialTreeMepiVanillaOptionSV : public BinomialTreeMepiVanillaOption
{

	double VolOfVol;
	double VolDrift;
	int	LegendreNb;

	BinomialTreeMepiVanillaOptionSV(int n0, int nz0,int nblg, double p0,double K0,double Vol0,double VolDrift0,
		double VolOfVol0,double r0,double b0,double YearlyFees0,double T0,double minExp,double maxExp,double riskFactor,double g0,double g,double minport,double maxport,int AveragingPeriodNb,int CallOrPut):
		BinomialTreeMepiVanillaOption(n0,nz0,p0,K0,Vol0,r0,b0,YearlyFees0,T0,minExp,maxExp,riskFactor,g0,g,minport,maxport,AveragingPeriodNb,CallOrPut),
		VolDrift(VolDrift0),VolOfVol(VolOfVol0),LegendreNb(nblg)
	{}
	


	double PriceSV();

/// to compute the delta, we add 1% of the exposed value to the portfolio, it simulates an increase of 1% ofthe underlying
	double PriceSV_delta();
};

////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												Global Functions 
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

double mepi_VanillaOption_STOBS(double date,double p0,double K, double T, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_SABR_Eq* SABR_Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,
						  int CallOrPut,int n,int nz,int nbLegendre);

double mepi_VanillaOption_STOBS_delta(double date,double p0,double K, double T, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_SABR_Eq* SABR_Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,
						  int CallOrPut,int n,int nz,int nbLegendre);

double mepi_VanillaOption_SABR(double date,double p0,double K, double maturity_date, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_PricingModel* Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,int resetfreq);

double mepi_VanillaOption_SABR_delta(double date,double p0,double K, double T, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_PricingModel* Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,double n,double shift);



double Export_Fund_VanillaOption(double f,double K,double T,double drift,double sig,
								 double VolDrift,double VolVol,double callput,int nbsteps);


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

