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
 *	\date November 2005
 */


#ifndef _INGPMODELS_FRMSV1F_H
#define _INGPMODELS_FRMSV1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

#include "gpbase/port.h"
//#include "gpbase/karmglob.h"
#include "gpmodels/SVModels.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsSFRM;
struct ARM_VanillaSwaptionArgSFRM;
//-----------------------------------------------------------------------------
// \class ARM_MarkovSV1F
// \brief
//  1 factor Markov Stochastic  pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_FRMSV1F : public ARM_SVModels
{
public:
	/// --------------- fast calibration caching arguments -----------------
    CC_IS_MUTABLE ARM_VanillaSwaptionArgSFRM* itsCurrentArg;
public:
	ARM_FRMSV1F( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL);
	ARM_FRMSV1F(const ARM_FRMSV1F& rhs);
	virtual ~ARM_FRMSV1F();

    ARM_FRMSV1F& operator = (const ARM_FRMSV1F& rhs);

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
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;

	/// Variable Notional swaptions can be priced in HW1F
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

    // Give local drifts and variances w.r.t. a given schedule
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalCorrels( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map );

	virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

	ARM_VectorPtr LocalDiscounts(size_t timeIdx, double dt, const ARM_PricingStatesPtr& states) const;

	double ComputeEquivalentVolatilityt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods) const;

	double ComputeEquivalentShiftt( 
		const string& curveName,
		double evalTime,  
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double lambdaValue) const;


	
	ARM_VanillaSwaptionArgSFRM* GetVanillaSwaptionArg( 
		const string& curveName,
		double evalTime,
		double swapResetTime,
		double swapNotional,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
        const ARM_GP_Matrix& strikesPerState,
		int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true		) const;

	
	/// function for the generic tree (2nd generation)
	virtual void VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol = true) const;

	/// function for the generic tree (2nd generation)
    virtual void EulerLocalDrifts(const ARM_GP_Vector& timeSteps,
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
	virtual string ExportShortName() const { return "LFMSV";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
