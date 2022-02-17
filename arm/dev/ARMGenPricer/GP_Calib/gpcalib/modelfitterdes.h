/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *! \file modelfitterdes.h
 *
 *  \brief file to describe the model fitters in the interface
 *
 *	\author   A. Schauly & A. Triki
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPCALIB_MODELFITTERDES_H
#define _INGPCALIB_MODELFITTERDES_H

///gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"

///gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/optimizer.h"

#include "typedef.h"
#include "calibmethod.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

/// base class
class ARM_ModelFitterDes : public ARM_RootObject
{
private:
	
	///Parameters exclusily used by Solver
	double itsXTolerance;
	double itsGradTolerance;

	///Parameters exclusily used by Optimizer
	double itsStepMax;
	bool itsLocalSearch;

	/// Parameters used by both
	size_t itsMax_iter;
	bool itsGetDetails;
	double itsFxTolerance;

	// Parameter for dichotomy
	double itsDichoXTol;
	double itsDichoMaxIter;

	ARM_SolverType itsSolverType; 
	ARM_OptimizerType itsOptimizerType;

	/// for iterative Calibration
	size_t itsNbCalibIter;
	
	void CopyNoCleanUp(const ARM_ModelFitterDes& rhs );

public:
	
	/// Constructor to build Nag optimizer
	ARM_ModelFitterDes( 
		ARM_OptimizerType algoType ,
		size_t max_iter			= OptimizerConstant::DefaultMax_Iter, 
		double Precision		= OptimizerConstant::DefaultFxTolerance, 
		double stepMax			= OptimizerConstant::DefaultStepMax,
		bool localSearch		= OptimizerConstant::DefaultLocalSearch,
		bool GetDetails			= OptimizerConstant::DefaultPrint);

	/// Constructor to build solver
	ARM_ModelFitterDes( 
		ARM_SolverType algoType	,
		size_t max_iter			= SolverConstant::DefaultMax_Iter , 
		double fxTol			= SolverConstant::DefaultFxTolerance, 
		double xTol				= SolverConstant::DefaultXTolerance,  
		double gradTol			= SolverConstant::DefaultGradTolerance,
		size_t dichoMaxIter		= SolverConstant::DefaultDichoMax_Iter,
		double dichoXTol		= SolverConstant::DefaultDichoXTolerance,
		bool printLevel			= SolverConstant::DefaultPrint);


	~ARM_ModelFitterDes() {};

	ARM_ModelFitterPtr CreateModelFitter(
		ARM_PricingModel* model,
		const ARM_StdPortfolioPtr portfolio, 
		const ARM_ModelParamVector&  CalibParams,
		ARM_MethodType methodType, 
		ARM_ModelFitterPtr& linkedmodelfitter,
        ARM_ModelFitterPtr& previousmodelfitter,
		ARM_MktTargetType targetType = ARM_CalibrationTarget::PriceTarget,
		size_t FactorNb = 0,
		size_t nbIter = 1, 
		ARM_DateStrip* numSchedule = NULL,
		const ARM_VanillaSecDensityPtrVector& numSecDensities = ARM_VanillaSecDensityPtrVector(0)	) const;

	/// accessors
	inline size_t GetMax_iter() const					{ return itsMax_iter;		}
	inline bool GetGetDetails() const					{ return itsGetDetails;		}
	inline double GetPrecision() const					{ return itsFxTolerance;	}
	inline double GetStepMax() const					{ return itsStepMax;		}
	inline double GetLocalSearch() const				{ return itsLocalSearch;	}
	inline ARM_SolverType GetSolverType() const			{ return itsSolverType;		}
	inline ARM_OptimizerType GetOptimizerType() const	{ return itsOptimizerType;	}
	inline size_t GetNbCalibIter() const				{ return itsNbCalibIter;	}
	inline void SetNbCalibIter(size_t NbCalibIter) 		{ itsNbCalibIter = NbCalibIter; }

	///setters
	inline void SetGetDetails(bool value)			{ itsGetDetails=value;		}
	inline void SetOptimizerType(ARM_OptimizerType optimizerType) 	{  itsOptimizerType = optimizerType;	}


	/// standard ARM support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
