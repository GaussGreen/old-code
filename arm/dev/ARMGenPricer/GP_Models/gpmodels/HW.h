/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW.h
 *
 *  \brief 
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_HW_H
#define _INGPMODELS_HW_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

///-----------------------------------------------------------------------------
/// \class HWBoundaryD1Function
/// \brief
///  Derivative of the boundary exercice function
///  for vanilla swaption closed form formula
///-----------------------------------------------------------------------------
class HWBoundaryD1Function
{
private:
    vector< double > itsPolyCoef;
    vector< double > itsExpCoef;

public:
    HWBoundaryD1Function(int size=0);
    virtual ~HWBoundaryD1Function();

    void SetPolyCoef(const vector< double >& polyCoef);
    void SetExpCoef(const vector< double >& expCoef);
    double operator () ( double x ) const;
};


//-----------------------------------------------------------------------------
// \class HWBoundaryFunction
// \brief
//  Boundary exercice function
//  for vanilla swaption closed form formula
//-----------------------------------------------------------------------------

class HWBoundaryFunction
{
private:
    vector< double > itsPolyCoef;
    vector< double > itsExpCoef;
    HWBoundaryD1Function* itsBoundaryD1;

public:
    HWBoundaryFunction(int size=0);
    virtual ~HWBoundaryFunction();

    void SetPolyCoef(const vector< double >& polyCoef);
    void SetExpCoef(const vector< double >& expCoef);

    double operator () ( double x ) const;
    HWBoundaryD1Function* Derivative() const;
};


//-----------------------------------------------------------------------------
// \class ARM_HullWhite
// \brief
//  Interface class for Hull & White pricing models
//-----------------------------------------------------------------------------

class ARM_HullWhite : public ARM_PricingModelIR
{
private :
    /// Precomputed data
	// !!! Just used for the equity model !!!
	ARM_MatrixVector itsLocalStdDevs;
	ARM_GP_MatrixPtr itsRelativeDrifts;

	// flags to manage SO pricing : normal approx or exact computation
	/// [0] : approx/exact
	/// [1] : normal/deepITM SO (used if exact computation)
	ARM_BoolVector itsSOFormulaFlags;
	
protected :
	bool isFwdStart(const ARM_GP_Vector& fixNotional, int& idx) const;

	void ComputeFlowZc(
		const string& curveName, 
		double evalTime,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
        const ARM_PricingStatesPtr& states,
		ARM_VectorPtr& zcFloatStart,
		vector< ARM_VectorPtr >& zcFlowCouponPay,
        vector< double >& flowPayTimes,
        vector< double >& couponPayPeriod ) const;


    ARM_VectorPtr VanillaSwaptionPrice(
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
        const ARM_VectorPtr& zcFloatStart,
        const vector< ARM_VectorPtr >& zcFlowPay,
        const vector< double >& zcCoef,
        const vector< double >& stdDev,
        const vector< double >& drift,
	    const ARM_GP_Matrix& strikesPerState,
        int payRec,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike,
		int idx = -1) const;

public:
	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

    static const double VOL_NORMAL_MAX;

	ARM_HullWhite( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL, const ARM_BoolVector& soFormulaFlags=ARM_BoolVector(2,true) );
	ARM_HullWhite(const ARM_HullWhite& rhs);
	virtual ~ARM_HullWhite();

    ARM_HullWhite& operator = (const ARM_HullWhite& rhs);

	bool IsApproxSOFormula() const {return itsSOFormulaFlags[0];}
	void SetIsApproxSOFormula(bool isApproxSOFormula) {itsSOFormulaFlags[0]=isApproxSOFormula;}
	bool IsDeepITMSOFormula() const {return !(itsSOFormulaFlags[1]);}
	void SetDeepITMSOFormula(bool isDeepITMSOFormula) {itsSOFormulaFlags[1] = !isDeepITMSOFormula;}
	const ARM_BoolVector& GetSOFormulaFlags() const { return itsSOFormulaFlags; }
	void SetSOFormulaFlags(const ARM_BoolVector& soFormulaFlags) {itsSOFormulaFlags = soFormulaFlags;}

    /// Closed form formulas
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
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;
    
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


	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;


	// Take care this function is just used to make the equity model work !!!!
	void MCModelStatesFromToNextTimeSpecialEquity(ARM_PricingStatesPtr& states,int timeIndex) const;

	virtual double UnderlyingCorrelation(  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const;

	virtual double UnderlyingCovariance (  string	underlyingType,
										   double	fromTime,
										   double	toTime,
										   double	startTime1,
										   double   endTime1,
										   double	startTime2,
										   double   endTime2,
										   double	startTime3,
										   double   endTime3,
										   double	startTime4,
										   double   endTime4) const;

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
	
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    /// Default initialisation of the model
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return true;}
	virtual bool NeedArrowDebreuPrices() const { return true; }

	/// no post init hence nothing
	virtual void PostInit(){};

	virtual ARM_GP_Vector* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );

// FIXMEFRED: mig.vc8 (23/05/2007 11:39:09):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	/// Volatilities and correlation time steps for PDEs
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;

	/// HW is MEanReverting
	virtual bool IsMeanRevertingCompatible() const { return true;}

    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_HullWhite : abstract class !"); }

	/// factor in front of sigma(t)*exp(-lambda*(Te-t)) in swap rate normal dynamics
	static double SwapRateVolFactor   ( double lambda,
										double expiryTime,
										double swapFloatStartTime,	// adjusted ...
										double swapFloatEndTime,	// adjusted ...
										double dfStart,
										double dfEnd,
										const ARM_GP_Vector& swapFixPayTimes,
										const ARM_GP_Vector& swapFixPayPeriods,
										const ARM_GP_Vector& swapFixPayDfs,
										// optional 
										double* swapRate = NULL,	
										double payTime = 0.0,
										double* numeraireVolFactor = NULL) ;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
