/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file numerical.h
 *  \brief Model fitter that calibrates a whole smile numerically
 * 
 *	\author  A. Schauly
 *	\version 1.0
 *	\date August 2005
 */

#ifndef _INGPCALIB_NUMCALIB_H
#define _INGPCALIB_NUMCALIB_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpcalib/modelfitter.h"
#include "gpnumlib/typedef.h"
#include "gpinfra/pricingstates.h"
#include "modelfitterdes.h"
#include "gpnumlib/numfunction.h"


CC_BEGIN_NAMESPACE( ARM )

/// ----------------------------------------------------
/// **** Function to solve for bootstrap calibration ***
/// ----------------------------------------------------
class ARM_ModelParamsMF;
class ARM_PDENumericalScheme;
class ARM_PricingStates;

struct BootstrapPriceFunc : public CC_NS( ARM_GP, UnaryFunc)<double,double>
{
private:
	ARM_PricingModel*		itsModel;
	ARM_ModelParamsMF*		itsModelParams;
	ARM_VanillaArg*			itsArgVanilla;
	ARM_PDENumericalScheme* itsPdeScheme;
	ARM_GP_VectorPtr		itsTimeSteps;
	mutable ARM_PricingStatesPtr	itsPricingStates;
	ARM_GP_Vector*			itsResetDates;
	int						itsCurIndex;

public:
	BootstrapPriceFunc  () {};
	~BootstrapPriceFunc () {};

	/// to be called before before each new step of the bootstrap
	void Init (	ARM_PricingModel*			model,
				ARM_ModelParamsMF*			modelParams,
				ARM_VanillaArg*				argVanilla,
				ARM_PDENumericalScheme*		pdeScheme,
				const ARM_GP_VectorPtr&		timeSteps,
				const ARM_PricingStatesPtr&	pricingStates,
				ARM_GP_Vector*				resetDates,
				int							curIndex);

	/// What it is for
	double operator () (double x) const;

	/// Handle variance squeezes
	void setVarianceSqueeze(double lowerBound) const;
};


/// ----------------------------------------
/// **** class ARM_NumericalModelFitter ***
/// ----------------------------------------
class ARM_DateStrip;

class ARM_NumericalModelFitter : public ARM_ModelFitter
{
private: 
	ARM_DateStripPtr					itsCalibSchedule; /// ASSOCIATION (owned by calib method)
	ARM_VanillaSecDensityPtrVector		itsCalibSecDensities;
	ARM_PricingStatesPtr				itsCalibStates;

	double itsStartTime, itsEndTime;

public:
	
	/// constructor / destructor ) operator =
	ARM_NumericalModelFitter() ;
    ARM_NumericalModelFitter(const ARM_NumericalModelFitter& rhs);

	ARM_NumericalModelFitter(	ARM_PricingModel* model,
								ARM_DateStrip* calibSchedule, 
								const ARM_VanillaSecDensityPtrVector& calibSecDensities,
								const ARM_StdPortfolioPtr portfolio		= ARM_StdPortfolioPtr(NULL), 
								
								/// -> Bootstrap specific
								const size_t max_iter					= SolverConstant::DefaultMax_Iter,
								bool getDetails							= false, 
								double fxTolerance						= SolverConstant::DefaultFxTolerance,
								double xTolerance						= SolverConstant::DefaultXTolerance,
								const size_t dicho_max_iter				= SolverConstant::DefaultDichoMax_Iter,
								double dichoXTolerance					= SolverConstant::DefaultDichoXTolerance,
								ARM_SolverType type						= ARM_ModelFitterSolverType::NewtonRaphson,
								size_t factorNb							= 0,
								size_t NbCalibIter						= 1,
								ARM_MktTargetType = ARM_CalibrationTarget::PriceTarget);


	virtual ~ARM_NumericalModelFitter();

	ARM_NumericalModelFitter& operator = (const ARM_NumericalModelFitter& rhs);
    
	/// accessors
	ARM_DateStripPtr						getCalibSchedule ()		const {return itsCalibSchedule;}
	const ARM_VanillaSecDensityPtrVector&	getCalibSecDensities()	const {return itsCalibSecDensities;}
	
	/// tenor of security in terms of number of periods in underlying calib schedule
	/// int ComputeCalibSecPeriodNumber (size_t resetIndex) const;
	const ARM_IntVector& getCalibSecPayDatesRelIndexes(size_t resetIndex) const;// {return itsCalibSecDensities[resetIndex]->getPayDatesRelIndexes();}
	const ARM_IntVector& getCalibSecPayDatesAbsIndexes(size_t resetIndex) const;// {return itsCalibSecDensities[resetIndex]->getPayDatesRelIndexes();}
	const ARM_GP_Vector& getCalibSecInterestTerms(size_t resetIndex)	  const;// {return itsCalibSecDensities[resetIndex]->getInterestTerms();}


	/// to set/get calibration dates
	virtual void setStartAndEndCalibrationTimes( double startTime, double endTime) { itsStartTime = startTime; itsEndTime =endTime;}

	/// calibration
    virtual void Calibrate();

	/// Update States During Calibration
	void UpdateStates( const ARM_PricingStatesPtr& states, ARM_GP_VectorPtr& probas, size_t toResetIdx ) const;

	bool BoostrapMode () const {return !GetPortfolio().IsNull();}

	/// standard ARM Object support
    virtual ARM_Object* Clone() const { return new ARM_NumericalModelFitter(*this); }
	virtual string ModelFitterName() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	/// --- All following methods are
	/// --- specific to HK bootstrap calibration
	virtual void BootstrapCalibration();
	virtual void Initialise() const ;
	virtual ARM_IntVectorPtr ComputePortfolioSecuritiesIdx() const ;
	void Validate ();
	BootstrapPriceFunc	itsFunction;
	UnaryFuncWithNumDerivative<double> * itsFunctionWithNumDerivative;
    T_SolverWithInitialGuess< UnaryFuncWithNumDerivative<double> > * itsSolver;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
