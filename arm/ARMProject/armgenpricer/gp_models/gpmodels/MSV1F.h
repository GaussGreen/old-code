/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV1F.h
 *
 *  \brief 
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_MSV1F_H
#define _INGPMODELS_MSV1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"

#include "gpbase/port.h"
//#include "gpbase/karmglob.h"
#include "MSV.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class ARM_MarkovSV1F
// \brief
//  1 factor Markov Stochastic  pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_MarkovSV1F : public ARM_MarkovSV
{
private :
	double itsForwardTerm;
	ARM_DateStrip* itsFixDateStrip;
	ARM_Curve* itsSwapRateCurve;
	ARM_Curve* itsSwapRateDerivativeCurve;
	bool itsIsSwapRate;
public:

	ARM_MarkovSV1F( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, double forwardTerm=1 , bool isSwapRate = true);
	ARM_MarkovSV1F(const ARM_MarkovSV1F& rhs);
	virtual ~ARM_MarkovSV1F();

    ARM_MarkovSV1F& operator = (const ARM_MarkovSV1F& rhs);

	/// DF
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	/// only for variable notional swaptions (numerical integration)
	/// if std swaption, this method will call the ARM_HullWhite method
	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;


	virtual ARM_VectorPtr VanillaSmiledSwaption(
		const string& curveName,
        double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
        double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixTimes,
        const std::vector<double>& fixPayPeriods,
        const std::vector<double>& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		const std::vector<double>& data,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	virtual ARM_VectorPtr VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
		double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
		double fwdEndTime,
		double fwdPeriod,
		const std::vector<double>& strikesPerState,
		int capFloor,
		const ARM_PricingStatesPtr& states) const;

	/// Variable Notional swaptions can be priced in HW1F
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

    // Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const std::vector<double>& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual void ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalCorrels( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map );

	virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	/// Accessors
	virtual ARM_DateStrip* GetDateStrip () const {return itsFixDateStrip;}
	virtual ARM_Curve* GetSwapRateCurve () const {return itsSwapRateCurve;}
	virtual ARM_Curve* GetSwapRateDerivativeCurve () const {return itsSwapRateDerivativeCurve;}
	

//	ARM_DateStrip* GetFloatDateStrip( const ARM_Date& startDate, const ARM_Date& endDate ) const;
	double ComputeFwdTerm (double evalTerm) const;

	std::vector<double>& ComputeFwdTerm_t(double evalTime,double F0,const ARM_PricingStatesPtr& states) const;

	double ComputeEta( double evalTime, double x) const;

	ARM_VectorPtr ComputeEtaStates( double evalTime, const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr LocalDiscounts(size_t timeIdx, double dt, const ARM_PricingStatesPtr& states) const;

	double ComputeDerivativeSt( const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		double* annuityValue,			// To store the Annuity Value
		double* swapRate,				// To store the SwapRate Value
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_PricingStatesPtr& states,
		double x = 0) const;

	double ComputeEquivalentVolatilityt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods) const;

	double ComputeEquivalentShiftt( 
		const string& curveName,
		double evalTime,  
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		double lambdaValue) const;
	
	/// function for the generic tree (2nd generation)
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// function for the generic tree (2nd generation)
    virtual void EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;


     /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );

	////////////////////////////////////////////////////////////////////////////////////////////////
	/// Datas that any irmodel must provide to work in a multiasset environment when used with infaltion
	////////////////////////////////////////////////////////////////////////////////////////////////

	/// This is supposed to provide Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
	/// Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
	virtual double IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const;

	/// Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
	virtual double IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const;

	/// Scalar product: returns Int_startTime^endTime Gamma(s,bondMaturity) * otherModelVolatility(s) ds
	virtual double VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const;

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;
	//////////////////////////////////////////////////////////////////////////////////////////////
    /// Standard ARM object support
	//////////////////////////////////////////////////////////////////////////////////////////////
	virtual ARM_Object* Clone() const;
	virtual string ExportShortName() const { return "LSV1F";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
