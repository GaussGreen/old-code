/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file cfmethod.h
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date May 2005
 */


#ifndef _INGPNUMMETHODS_CFMETHOD_H
#define _INGPNUMMETHODS_CFMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/typedef.h"
#include "gpnumlib/typedef.h"
#include "argconvdefault.h"


CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////
/// \class ARM_CFMethod
/// \brief class for Monte Carlo
/// This numerical method valuate all the keyword with a closed form
/// at the asof date of the model
//////////////////////////////////////////////
class ARM_CFMethod : public ARM_NumMethod
{
public:
	//enum
	enum CFMethodType
	{
		Analytic =0,
		Integral,
		Unknown ,
	};
	typedef ARM_CFMethod::CFMethodType  ARM_CFMethodType;
	// Constructor & Destructor
	explicit ARM_CFMethod(CFMethodType cf_method = Integral, 
			const ARM_GP_Matrix& Parameters=ARM_GP_Matrix(0));
	ARM_CFMethod(const ARM_CFMethod & rhs);
	ASSIGN_OPERATOR(ARM_CFMethod)
	virtual ~ARM_CFMethod() {};
private:
	ARM_CFMethodType itsCFmethod;
	ARM_GP_Matrix itsCFparameters;
public:
	/// forward pricing
    virtual GP_PricingDirection GetPricingDirection() const {return ARM_NumMethod::GP_FWDLOOKING;}
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return ARM_NumMethod::GP_FWDLOOKING;}
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return 0; }

    /// Initialisation of the numerical method
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
	virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model );

	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {}
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const {};

    /// Numerical induct for a method
	virtual ARM_PricingStatesPtr Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, double toTime);

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Does we do not compute otherpayoffs for too large simluations
	virtual bool GetOtherPayoffsFlag() const { return true; }

	/// accessors to the buckets
	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const;

	//Accessor to the type of CF
	ARM_CFMethodType GetCFMethodType()  {return itsCFmethod;}
	ARM_GP_Matrix GetCFParameters() {return itsCFparameters;}
	void SetCFParameters( ARM_GP_Matrix cFparameters) {itsCFparameters=cFparameters;}


	/// ---------- Control Variate
	/// ability to have control variate on various instruments
	virtual void ControlVariateOnInstrument( ARM_VectorPtr& numericInstrument, double correctValue,
	const string& curveName, double evalTime, const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const {};

	virtual double ConvertEvalDate(double evalDate, double asOfDate) const { return asOfDate; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

