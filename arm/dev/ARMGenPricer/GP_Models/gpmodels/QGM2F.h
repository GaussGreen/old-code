/*!
 *
 * Copyright (c) CDC IXIS CM May 2006 Paris
 *
 *	\file QGM2F.h
 *
 *  \brief 
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date May 2006
 */


#ifndef _INGPMODELS_QGM2F_H
#define _INGPMODELS_QGM2F_H


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "gpbase/gpmatrix.h"

#include <vector>
CC_USING_NS(std,vector)

#include <map>
CC_USING_NS(std,map)

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsQGM2F;

//-----------------------------------------------------------------------------
// \class ARM_QGM2F
// \brief
//  Class for 1F Quadratic Gaussian Model
//-----------------------------------------------------------------------------

class ARM_QGM2F : public ARM_PricingModelIR
{
private :
    void CopyNoCleanUp(const ARM_QGM2F& rhs);

    /// Precomputed surface of the reconstruction formula
    ARM_GP_Vector      itsTime;        // t
    ARM_GP_Vector      itsMaturity;    // T
    ARM_GP_Matrix       itsA;           // a(t,T)
    ARM_GP_Matrix       itsB;           // b(t,T)
    ARM_GP_Matrix       itsC;           // c(t,T)

    
    
public:
	ARM_QGM2F(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQGM2F& params);
	ARM_QGM2F(const ARM_QGM2F& rhs);
	virtual ~ARM_QGM2F();

    ARM_QGM2F& operator = (const ARM_QGM2F& rhs);

	void InitABC(const ARM_GP_Vector& times, const vector< ARM_GP_Vector >& maturities);
	
    /// Computation of A(t,T), B(t,T) and C(t,T)
    double A(double t, double T) const;
    double B(double t, double T) const;
    double C(double t, double T) const;

    /// Reconstruction formula
	virtual ARM_VectorPtr DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const;

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
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

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
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;
    
    virtual ARM_VectorPtr VanillaDigital(
		const string& curveName, 
		double evalTime,
		double payTime,			
		double period,
        double payNotional,
		double fwdResetTime,	
		double fwdStartTime,
        double fwdEndTime,
 		double fwdPeriod,
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;  

public:
    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter);
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);
    virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);
    virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const;

    virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos );
// FIXMEFRED: mig.vc8 (25/05/2007 15:28:21):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return true;}

    
    virtual void IntegratedLocalDrifts(
		const ARM_GP_Vector& timeSteps,
        ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const;

	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const;

    virtual void NumMethodStateGlobalVariances(
        const ARM_GP_Vector& timeSteps,
        ARM_MatrixVector& variances) const;

    virtual bool NeedsToCholeskyDecomposeFactors( ) const {return false;}

	virtual ARM_BoolVector NeedMCIntegProcess() const { return ARM_BoolVector(1, true); };

    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0;}

	virtual ARM_VectorPtr  VanillaSpreadOptionLet(const string& curveName,
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

	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LQM1F";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
