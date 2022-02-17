/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodelinflation.h
 *
 *  \brief inflation version of the pricing model
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */


#ifndef _INGPINFRA_PRICINGMODELINFLATION_H
#define _INGPINFRA_PRICINGMODELINFLATION_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

#include "typedef.h"
#include "pricingmodelir.h"
#include "pricingfunctioninflation.h"
#include "gpinflation/infcurv.h"

CC_BEGIN_NAMESPACE( ARM ) /// macro for namespace ... define namespace only if supported


///////////////////////////////////////////////////////
/// \class ARM_PricingModelIR
/// \brief
/// This abstract class is the standard
/// interface for interest rate pricing models
/// Derived from the ARM_PricingModel
///////////////////////////////////////////////////////
class ARM_PricingModelInflation :  public ARM_PricingModel,
								   public ARM_PricingFuncInflation
{
private: 
	ARM_InfCurvPtr itsInfCurve;
	ARM_PricingModelIR* itsIRModel;
	ARM_PricingModelPtr itsIRModelPtr; /// This is a kludge until I find something better. 
	double itsCorrelWithBonds;

	/// ============================================================= ///
    /// ================== elementary pricing functions ============= ///
	/// ============================================================= ///

	ARM_GP_VectorPtr DefaultYoYSwapLeg(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const;

	ARM_GP_VectorPtr DefaultOATSwapLeg(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const;

	ARM_VectorPtr DefaultAnnuity(
		const string& curveName,
		double evalTime,
		const ARM_DateStripPtr& PayDateStrip, 
		const ARM_PricingStatesPtr& states) const;



public:
	/// ============================================================= ///
    /// ================== constructors/desctructors ================ ///
	/// ============================================================= ///

	ARM_PricingModelInflation( const ARM_InfCurvPtr& infc=ARM_InfCurvPtr(NULL), 
		const ARM_ModelParams* params=NULL  ); 
	ARM_PricingModelInflation( const ARM_PricingModelInflation& rhs);
	virtual ~ARM_PricingModelInflation();
    ARM_PricingModelInflation& operator = (const ARM_PricingModelInflation& rhs);

	/// ============================================================= ///
    /// ====================== accessors  =========================== ///
	/// ============================================================= ///

	inline const ARM_PricingModelIR* getIRModel() const { return itsIRModel;}
	inline void setIRModel( ARM_PricingModelIR* IRModel ) { itsIRModel = IRModel;}
	inline void setIRModel( const ARM_PricingModelPtr& IRModel )
		{ itsIRModelPtr = IRModel; itsIRModel = dynamic_cast<ARM_PricingModelIR*> (&*IRModel); if (!itsIRModel) ARM_THROW( ERR_INVALID_ARGUMENT, " Inflation Models are supposed to work with IRModels. " ); }
	inline ARM_InfCurvPtr getInfCurve() const { return itsInfCurve;}
	inline void setInfCurve( const ARM_InfCurvPtr& InfCurve ) { itsInfCurve = InfCurve;}
	inline double getCorrelWithBonds() const { return itsCorrelWithBonds; }
	inline void setCorrelWithBonds( double CorrelWithBonds ) { itsCorrelWithBonds = CorrelWithBonds; }
	virtual ARM_Date GetAsOfDate() const { return itsInfCurve->GetAsOf(); }


	/// ============================================================= ///
    /// ========== DiscountFactor(from pricingmodel ================= ///
	/// ============================================================= ///

    virtual ARM_VectorPtr DiscountFactor( const string& curveName, double evalTime, double maturityTime, 
        const ARM_PricingStatesPtr& states) const;

	/// ============================================================= ///
    /// ================== elementary pricing functions ============= ///
	/// ============================================================= ///

	/// CPI Spot
	virtual ARM_GP_VectorPtr CPISpot( 
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, string DCFLag, long DailyInterp,
		string ResetLag,
		const ARM_PricingStatesPtr& states) const { return ARM_GP_VectorPtr(NULL); }

	/// CPI Forward
	virtual ARM_GP_VectorPtr CPIForward(
		const string& InfcurveName, 
		double evalTime, 
		double CPITime, 
		double FixingTime, 
		const ARM_PricingStatesPtr& states) const { return ARM_GP_VectorPtr(NULL); }

	/// Convexity Adjustment
	virtual ARM_GP_VectorPtr ConvexityAdjustment(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const { return ARM_GP_VectorPtr(NULL); }

	/// forward Ratio
	virtual ARM_GP_VectorPtr ForwardCPIRatio(
		const string& InfcurveName,
        double evalTime,
		double tenor,
		double CPITime,
		double maturityTime,
        const ARM_PricingStatesPtr& states) const { return ARM_GP_VectorPtr(NULL); }

	/// Year on year swap rate
	virtual ARM_GP_VectorPtr YoYSwapRate(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const;

 	virtual ARM_GP_VectorPtr YoYSwap(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr OATSwapRate(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr OATSwap(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		const ARM_DateStripPtr& fixedDateStrip,
		double itsCoupon,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_GP_VectorPtr YoYCapFloor(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		int CapFloor,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const = 0;

	virtual ARM_GP_VectorPtr OATCapFloor(
		const string& irCurveName,
		const string& infCurveName, 
		double evalTime,
		double Strike,
		double FloatMargin, 
		int CapFloor,
		const ARM_DateStripPtr& numDateStrip,
		const ARM_DateStripPtr& denomDateStrip,
		double itsSpread,
		const ARM_PricingStatesPtr& states) const {ARM_THROW( ERR_INVALID_ARGUMENT, "OATCapFloor Not implemented in this model!" ); }

	virtual ARM_GP_VectorPtr ZCCap( ) const { return ARM_GP_VectorPtr(NULL); }
	
	virtual int GetType() const = 0;
    
    /// Standard ARM object support
	virtual string toString(const string& indent="",const string& nextIndent="") const { return indent + string("ARM_PricingModelInflation : abstract class !"); }
};



CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
