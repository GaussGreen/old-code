
#ifndef _ICM_BASEPRICER_H
#define _ICM_BASEPRICER_H

#include "ICMKernel/pricer/icm_pricer.h"

class ICM_Credit_Index; 
class ICM_ModelMultiCurves; 

class ICM_BasePricer : public ICM_Pricer
{
public:
	ICM_BasePricer();
	virtual ~ICM_BasePricer();
protected :
	void Init() {itspPricer = NULL;
				 itsPVoptimum = 0.;}
public:
	void Set(ICM_Credit_Index* sec,ICM_ModelMultiCurves*mod,const ICM_Parameters&params,const ARM_Date&asof) ;
	virtual void DoPriceVector(qVECTMETH measure) ;
private: 
	void doPriceIndexBase(); // compute base from fair spread
	void doPriceMarketBase(); // compute base from = of MTM using implied as coupon spread 
	void doPriceTrueMarketBase(); // compute base from = of MTM using coupon spread 
private:
	virtual double ComputeSpread(const double& MtM = 0.) ; 
 	virtual	double ComputeImpliedVol(const double& Price)  ;
	virtual double Accrued() ; 
	virtual double FeeLegPV ()  ;
	virtual double DefLegPV () ; 
	virtual double ComputeSensitivity(const qSENSITIVITY_TYPE& typesensi, 
		const std::string& plot,
		const std::string& label,
		double epsilon, double epsvalueGamma ) ; 
	virtual double ComputeDuration() ;	
private:
	ICM_BasePricer(const ICM_BasePricer &); // NA
	ICM_BasePricer & operator=(const ICM_BasePricer &); // NA 
	double Evaluate(const double& spread);
// FIXMEFRED: mig.vc8 (28/05/2007 10:39:15):missing return type
	void ResetSpread(const double& spread);
	ICM_Pricer* itspPricer; // for BasePricerMethod : TrueMarketPrice
	double itsPVoptimum;
}; 


#endif // _ICM_BASEPRICER_H 
