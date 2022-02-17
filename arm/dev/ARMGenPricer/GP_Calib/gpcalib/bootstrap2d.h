/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file bootstrap2d.h
 *
 *  \brief 
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPCALIB_BOOTSTRAP2D_H
#define _INGPCALIB_BOOTSTRAP2D_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpnumlib/typedef.h"
#include "modelfitter.h"

#include "modelfitterdes.h"

/// forward declaration
class ARM_Security;

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class FunctionToSolve;

class ARM_Bootstrap2D : public ARM_ModelFitter
{
private:
	
	FunctionToSolveWithDerivativePtr itsFunction;
    ModifiedNRSolverPtr itsSolver;

	/// functions
    void CopyNoCleanUp(const ARM_Bootstrap2D& rhs );
    void CleanUp();
    double RootFinding(double mktPrice, 
		ARM_Security* sec,
		double guess_initial,
		double lowerbound,
		double upperbound,
		int index, 
		double time, 
		double expiry = 0.0 );
    
	/// standard methods
	void Validate();
    void Initialise();

public:
    ARM_Bootstrap2D( ARM_PricingModel* model, 
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& calibParam,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter	= ARM_ModelFitterPtr(NULL ),
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
    ARM_Bootstrap2D(const ARM_Bootstrap2D& rhs);
	ARM_Bootstrap2D& operator = (const ARM_Bootstrap2D& rhs);
    virtual ~ARM_Bootstrap2D();


    /// the core of the calibration
	virtual void Calibrate();

	///Accessors
	inline ModifiedNRSolverPtr GetSolver() const {return itsSolver;}

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	string  detailedString( const string& indent ) const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
