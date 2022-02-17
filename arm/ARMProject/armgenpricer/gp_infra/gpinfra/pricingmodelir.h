/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodelir.h
 *
 *  \brief interest rate version of the pricing model
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_PRICINGMODELIR_H
#define _INGPINFRA_PRICINGMODELIR_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"

#include "typedef.h"
#include "pricingmodel.h"
#include "pricingfunctionir.h"

CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported


///////////////////////////////////////////////////////
/// \class ARM_PricingModelIR
/// \brief
/// This abstract class is the standard
/// interface for interest rate pricing models
/// Derived from the ARM_PricingModel
///////////////////////////////////////////////////////
class ARM_PricingModelIR :  public ARM_PricingModel,
                            public ARM_PricingFunctionIR
{
public:
	ARM_PricingModelIR( const ARM_ZeroCurvePtr& zc=ARM_ZeroCurvePtr(NULL), const ARM_ModelParams* params=NULL, ARM_DensityFunctor* densityFct = NULL); 
	ARM_PricingModelIR(const ARM_PricingModelIR& rhs);
	virtual ~ARM_PricingModelIR();
    ARM_PricingModelIR& operator = (const ARM_PricingModelIR& rhs);

    /// ================== elementary pricing functions =============
	/// Defauft Libor provided but may be redefined by
    /// a Libor Market Model
	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double fwdResetTime,    // for convexity adjustment...
        double payTime,         //... in derived classes
        const ARM_PricingStatesPtr& states) const;

	/// Defauft annuity provided but may be redefined
	virtual ARM_VectorPtr Annuity(
		const string& curveName, 
        double evalTime,
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const;

	/// Defauft annuity provided but may be redefined
	/// The calculation integrate the nominal
	virtual ARM_VectorPtr AnnuityWithNominal(
		const string& curveName, 
        double evalTime,
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fixNominal,
        const ARM_PricingStatesPtr& states) const;

	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods,
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
        const std::vector<double>& margin,
        bool isDbleNotional,
		const ARM_VectorPtr& PreComputedAnnuity,
		const ARM_PricingStatesPtr& states) const;

	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuityAndNominal(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods,
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
        const std::vector<double>& margin,
		const ARM_VectorPtr& PreComputedAnnuity,
		const std::vector<double>& floatNotional,
		const ARM_PricingStatesPtr& states) const;
	
	/// Defauft Swap provided but may be redefined by
    /// a Swap Market Model
	virtual ARM_VectorPtr NPVSwap(
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& floatPayTimes, 
		const std::vector<double>& floatPayPeriods, 
		const std::vector<double>& margin,
        bool isDbleNotional,
		const std::vector<double>& FixNotional,
		const std::vector<double>& FloatNotional,
		const ARM_GP_Matrix& strikesPerState,
        int payRec,
        const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& PayTimes, 
		const std::vector<double>& PayPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const;
		
	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const;

	/// conrtol variate on libor
	virtual ARM_VectorPtr LiborWithControlVariate( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double fwdResetTime,    // for convexity adjustment...
        double payTime,         //... in derived classes
        const ARM_PricingStatesPtr& states) const;

     virtual double ImpliedVol(const ARM_VanillaArg& arg) const ;
	
	virtual int GetType() const;
    
    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_PricingModelIR : abstract class !"); }
};



CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
