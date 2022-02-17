/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SVModels.h
 *	\ Class for general SVModels (MSV & FRMSV) that use the Pieterbarg Approach
 *  \brief 
 *
 *	\author  A  Triki
 *	\version 1.0
 *	\date November 2005
 */ 


#ifndef _INGPMODELS_SVMODELS_H
#define _INGPMODELS_SVMODELS_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpinfra/pricingmodelir.h"
#include "gpnumlib/odefunctions.h"


#include <vector>
CC_USING_NS(std,vector)

class ARM_IRIndex;

CC_BEGIN_NAMESPACE( ARM )
//-----------------------------------------------------------------------------
// \class ARM_MarkovSV
// \brief
//  Interface class for Hull & White pricing models
//-----------------------------------------------------------------------------

/// Internal structures To Store Swaption Pricing Data
struct ARM_GeneralSwaptionData
{	
	// Store Some Swaption elts for the quick computation of Lambdat and Shiftt
	ARM_GeneralSwaptionData(){};
	ARM_GeneralSwaptionData (int sizeVector);
	
	std::vector<double> itsNewSchedule;
	// Datas for the calculation of the equivalent constant VolOfVol
	std::vector<double> itsEpsPart_1;	 // Vector to Store Int (Eps² exp (2 * theta *t ) dt , ti, ti+1 )
	std::vector<double> itsEpsPart_2;  // Vector to Store Int ( exp (2 * theta *t ) dt , ti, ti+1 )
	std::vector<double> itsEpsPart_3;  // Vector to Store Int (sigma² exp (- theta *t ) dt , ti, ti+1 )
	
	
	// Datas for the calculation of the equivalent Shift
	std::vector<double> itsBetaPart_1; // to store Integral(sigma²,0,ti)
	std::vector<double> itsBetaPart_2; // to store Integral(sigma²*sinh(ks)/k,0,ti)
}; 


class ARM_SVModels : public ARM_PricingModelIR
{
protected :
    /// Precomputed data
	double itsIntegratedVol;	
	ARM_GeneralSwaptionData itsStoredData;
	double ComputeInitialPoint( double cstShift);

public:

	virtual	double ComputePhi(const ARM_ODEFunc& test,
							double cstShiftValue,
							double resetTime) const;

	virtual double ComputePhi_Zero( double cstShift,
							double resetTime) const;

	virtual	double ComputeMu(double cstShiftValue,
							double integratedVol) const;



	virtual ARM_GP_MatrixPtr MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const;

    static const double VOL_NORMAL_MAX;

	ARM_SVModels( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params=NULL );
	ARM_SVModels(const ARM_SVModels& rhs);
	virtual ~ARM_SVModels();

    ARM_SVModels& operator = (const ARM_SVModels& rhs);

	/// Functions to update At, Bt and Ct in the Riccati Equation f' = At f² + Bt f + Ct
	virtual double UpdateRiccatiCt(const string& curveName,
									double t,
									double floatStartTime,
									double floatEndTime,
									const std::vector<double>& fixPayTimes,
									const std::vector<double>& fixPayPeriods);
	/// For Tests
	virtual double UpdateRiccatiAt(double t);
	virtual double UpdateRiccatiBt(double t);

	virtual void PrecomputeDatas(
				const string& curveName,
				double resetTime, 
				double floatStartTime,
				double floatEndTime,
				const std::vector<double>& fixPayTimes,
				const std::vector<double>& fixPayPeriods)	;

	virtual ARM_DateStrip* GetFloatDateStrip( const ARM_Date& startDate, const ARM_Date& endDate,const ARM_IRIndex* pIndex ) const;

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
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

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
		const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const;  


	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const;



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
	
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;

    virtual bool SupportBackwardInduction() const {	return true;}
    virtual bool SupportForwardInduction()  const {	return true;}
	virtual bool SupportAnalyticMarginal()  const {	return false;}

	virtual double ComputeEquivalentCstVolatility( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_PricingStatesPtr& states) const;

	virtual double ComputeEquivalentCstShift( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		double effetiveCstVolOfVol,
		double* globalVol) const;

	virtual double ComputeEquivalentCstVolOfVol( 
		const string& curveName,
		double resetTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods) const;

	virtual double ComputeEquivalentVolatilityt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods) const = 0;

	virtual double ComputeEquivalentShiftt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		double lambdaValue) const = 0;

	/// no post init hence nothing
	virtual void PostInit(){};

	virtual std::vector<double>* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos );

// FIXMEFRED: mig.vc8 (28/05/2007 12:03:41):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	/// Volatilities and correlation time steps for PDEs
	virtual ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;

	/// MSV is MEanReverting
	virtual bool IsMeanRevertingCompatible() const { return true;}

	virtual bool ValidateModelParams(const ARM_ModelParams& params) const = 0;

	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_SVModels : abstract class !"); }

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
