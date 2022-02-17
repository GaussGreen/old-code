/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW1F.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date September 2003
 */


#ifndef _INGPMODELS_HW1F_H
#define _INGPMODELS_HW1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

#include "gpbase/port.h"
//#include "gpbase/karmglob.h"
#include "HW.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class HW1FVarianceFunction
// \brief
//  Variance function to find the time to reach a given variance
//-----------------------------------------------------------------------------

class HW1FVarianceFunction
{
private:
    const ARM_ModelParams* itsModelParams;

public:
    HW1FVarianceFunction(const ARM_ModelParams* modelParams);
    virtual ~HW1FVarianceFunction();
    double operator () ( double x ) const;
};


//-----------------------------------------------------------------------------
// \class ARM_HullWhite1F
// \brief
//  1 factor Hull & White pricing model for closed form,
//  backward and forward diffusion abilities
//-----------------------------------------------------------------------------

class ARM_HullWhite1F : public ARM_HullWhite
{
public:
	ARM_HullWhite1F( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, const ARM_BoolVector& soFormulaFlags=ARM_BoolVector(2,true)  );
	ARM_HullWhite1F(const ARM_HullWhite1F& rhs);
	virtual ~ARM_HullWhite1F();

    ARM_HullWhite1F& operator = (const ARM_HullWhite1F& rhs);

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

		ARM_VectorPtr VanillaSpreadOptionLet(
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
			const ARM_GP_Vector& strikes,
			double swapLongFloatStartTime,
			double swapLongFloatEndTime,
			const ARM_GP_Vector& swapLongFixPayTimes,
			const ARM_GP_Vector& swapLongFixPayPeriods,
			double swapShortFloatStartTime,
			double swapShortFloatEndTime,
			const ARM_GP_Vector& swapShortFixPayTimes,
			const ARM_GP_Vector& swapShortFixPayPeriods,
			const ARM_PricingStatesPtr& states) const;

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

	virtual double VarianceToTime(double var,double minTime,double maxTime) const;

	virtual bool NeedsToCholeskyDecomposeFactors( ) const { return false; }

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(1, true); };

	ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

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
	
	//Calculation of the short rate needed for PDE computing when the numeraire is cash
	virtual ARM_VectorPtr RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded=false ) const;
	//Calculation of the coefficients of the PDE of a HW1f model at time of index timeIdx	
	virtual void UpdatePDE3DCoeffs(
		size_t timeIdx,
		const ARM_GP_MatrixPtr& numMethodStates,
		ARM_GP_VectorPtr& qxx,
		ARM_GP_VectorPtr& qyy,
		ARM_GP_VectorPtr& qzz,
		ARM_GP_VectorPtr& qxy,
		ARM_GP_VectorPtr& qyz,
		ARM_GP_VectorPtr& qzx,
		ARM_GP_VectorPtr& px,
		ARM_GP_VectorPtr& py,
		ARM_GP_VectorPtr& pz,
		ARM_GP_VectorPtr& o
		) const;


	/// for cash numeraire
    virtual ARM_VectorPtr LocalDiscounts(size_t timeIdx, double dt, const ARM_PricingStatesPtr& states) const;


	//////////////////////////////////////////////////////////////////////////////////////////////
    /// Standard ARM object support
	//////////////////////////////////////////////////////////////////////////////////////////////
	virtual ARM_Object* Clone() const;
	virtual string ExportShortName() const { return "LHW1M";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
