/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file AnalyticIRModel.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_ANALYTICIRMODEL_H
#define _INGPMODELS_ANALYTICIRMODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/pricingmodelir.h"
#include <string>
CC_USING_NS(std,string)

/// forward declaration
class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

class ARM_AnalyticIRModel : public ARM_PricingModelIR
{
private:
	void ThrowError( const string& mssg );
	void ValidatePricingStates( const ARM_PricingStatesPtr& states ) const;

protected:
	enum CalibSecurityType
	{
		ARM_SWAPTION_TYPE = 0, 
		ARM_CAP_TYPE
	};
	static string CalibSecurityTypeString[];
	static void ComputeTenorAndExpiry( ARM_Security* security, double& expiry, double& tenor, CalibSecurityType& type, double asOfDate );
	static CalibSecurityType GetInstType( ARM_Security* security );

public:
	ARM_AnalyticIRModel( const ARM_ZeroCurvePtr& zc=ARM_ZeroCurvePtr(NULL), const ARM_ModelParams* params=NULL, ARM_DensityFunctor* densityFct = NULL)
	: ARM_PricingModelIR(zc,params,densityFct) {}
	ARM_AnalyticIRModel( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, ARM_DensityFunctor* densityFct = NULL);
	ARM_AnalyticIRModel( const ARM_AnalyticIRModel& rhs);
    ARM_AnalyticIRModel& operator=(const ARM_AnalyticIRModel& rhs);
	virtual ~ARM_AnalyticIRModel() {};

	static double GetMatchingTenor( double tenor );
    /// Reconstruction formula
	virtual ARM_VectorPtr DiscountFactor( const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

    /// Calibration purpose
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

    /// numerical part of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
    virtual void PostInit(){};
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const {return ARM_PricingStatesPtr(NULL);}
    virtual std::vector<double>* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos ) {return NULL;}
// FIXMEFRED: mig.vc8 (22/05/2007 18:11:55):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}

    virtual bool SupportBackwardInduction() const {	return false;}
    virtual bool SupportForwardInduction()  const {	return false;}
	virtual bool SupportAnalyticMarginal()  const {	return false;}

    // Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const {};

	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const {};
	
	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const {};
    
    virtual bool NeedsToCholeskyDecomposeFactors( ) const {return false;}
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0;}

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
			const string& curveName,
			double evalTime,
			int callPut,
			double startTime,
			double endTime,
			double resetTime,
			double payTime,
			double payPeriod,
			double notional,
			double coeffLong,
			double coeffShort,
			const std::vector<double>& strikes,
			double swapLongFloatStartTime,
			double swapLongFloatEndTime,
			const std::vector<double>& swapLongFixPayTimes,
			const std::vector<double>& swapLongFixPayPeriods,
			double swapShortFloatStartTime,
			double swapShortFloatEndTime,
			const std::vector<double>& swapShortFixPayTimes,
			const std::vector<double>& swapShortFixPayPeriods,
			const ARM_PricingStatesPtr& states) const;

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const {}
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const {}
	double cDuration(double Rate, double Term, int F, double t) const;
	double cConvexity(double Rate, double Term, int F, double t) const;
	double DecompRate(double Rate, int Freq) const;
	double InverseDecompRate(double Rate, int Freq) const;
	double ForwardDecap( double ForwardF, double Spread, double VolF, double Freq, double Time ) const;
	double VolDecap( double ForwardF, double VolF, double Freq ) const;	
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
