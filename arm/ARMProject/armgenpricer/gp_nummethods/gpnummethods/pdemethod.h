/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pdemethod.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPINFRA_PDEMETHOD_H
#define _INGPINFRA_PDEMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/typedef.h"
#include "gpnumlib/typedef.h"
#include "pdenumericalschemes.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_DiscretisationScheme;


class ARM_PDEMethod : public ARM_NumMethod
{
private:
	size_t itsAdditionalTimeStepsNb;  // additional TimeSteps per Year specified by user. 
	double itsStdDevNb;				  // How many StdDev LastEvent

	double itsFirstInductTime;

	size_t itsTimeItersNb;
	size_t itsSpaceItersNb;
	
    void CopyNoCleanUp(const ARM_PDEMethod& rhs);
    void CleanUp();

	bool itsOtherPayoffsFlag;

	string itsMethodName;

	ARM_PDENumericalScheme * itsPDENumericalScheme;

	/// To handle case where we call induct several times
	/// on same time interval (HK bootstrap)
	int itsPrevLastTimeIdx;
	double itsPrevToTime;

public:
	ARM_PDEMethod( ARM_PDENumericalScheme * numericalScheme, size_t timeItersNb, size_t spaceItersNb, double stdDevNb );
	ARM_PDEMethod(const ARM_PDEMethod& rhs);
	ARM_PDEMethod& operator=(const ARM_PDEMethod& rhs);

	virtual ~ARM_PDEMethod();

	/// accessors
	const double& getStdDevNb () {return itsStdDevNb;}

	void setStdDevNb(double stdDevNb) 
	{	itsStdDevNb = stdDevNb;
		if( itsPDENumericalScheme)
			itsPDENumericalScheme->setStdDevNb(stdDevNb);
	}

	virtual ARM_PDENumericalScheme* GetPDENumericalScheme() const {return itsPDENumericalScheme;}
	virtual bool SupportFrontierComputation() const { return GetPDENumericalScheme()->IsOneDim(); }

	/// backward pricing
    virtual GP_PricingDirection GetPricingDirection() const {return ARM_NumMethod::GP_BCKWDLOOKING;}
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return ARM_NumMethod::GP_BCKWDLOOKING;}
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[GetTimeSteps()->size()-1]; }
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const;

    /// Initialisation of the numerical method
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
	virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model );

	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {}

    /// Numerical induct for a method
	virtual ARM_PricingStatesPtr Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, double toTime);

	// Return the global variance of the numerical method
	virtual const ARM_MatrixVector& GetNumMethodStateGlobalVars() const {return itsPDENumericalScheme->GetGlobalVars(); };

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Does we do not compute otherpayoffs for too large simluations
	virtual bool GetOtherPayoffsFlag() const { return itsOtherPayoffsFlag; }

	/// returns the vector of buckets
// FIXMEFRED: mig.vc8 (22/05/2007 18:10:55): cast
	virtual ARM_VectorPtr GetBuckets() const { return static_cast<ARM_VectorPtr>(new std::vector<double>(1,1));} 
	virtual size_t GetBucketIndex() const {return 0; };
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const std::vector<double>& eventTimes) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const;

	/// ---------- Control Variate
	/// ability to have control variate on various instruments
	virtual void ControlVariateOnInstrument( ARM_VectorPtr& numericInstrument, double correctValue,
		const string& curveName, double evalTime, const ARM_PricingStatesPtr& states, const ARM_PricingModel& model ) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LPDEM";}
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

