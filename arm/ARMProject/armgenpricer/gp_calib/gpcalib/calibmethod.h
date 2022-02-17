/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file calibmethod.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#ifndef _INGPCALIB_CALIBMETHOD_H
#define _INGPCALIB_CALIBMETHOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/numericconstant.h"


/// gpcalib
#include "typedef.h"
#include "gpinfra/calibdirection.h"
#include "gpcalib/vanillasecuritydensity.h"

/// ARM Kernel
//#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelFitter;
class ARM_PricingModel;
class ARM_Warning;


///////////////////////////////////////////////////////////////
/// \class ARM_
/// \brief
///  Interface class for calibration of pricing model
///////////////////////////////////////////////////////////////

class ARM_CalibMethod : public ARM_RootObject
{
public:
    /// ugly public private mixed up because of the conventional
    /// way to show private stuff before public
    

private:
    ARM_StdPortfolioPtr  itsPortfolio;
    ARM_ModelParamVector itsCalibParams;

	/// Type of calibration startegy
	ARM_MethodType itsMethodType;

    /// can be only of the enum type above
    ARM_MktTargetType itsTargetFuncType;
    ARM_CalibDirection itsCalibDirection;
	
	// Calib Method Descriptor
	ARM_ModelFitterDes* itsModelFitterDes;

    ///linked method that we use to calibrate other parameter(s)
    ARM_CalibMethod* itsLinkedMethod;

    ///Previous method that we use to calibrate previous parameter(s)
    ARM_CalibMethod* itsPreviousMethod;

	///Next method 
    ARM_CalibMethod* itsNextMethod;

    /// model fitter implements all the logic of calibration
	ARM_ModelFitterPtr itsModelFitter;

	/// for multi-factor model
	size_t itsFactorNb; 

	///To reiterate modelfitter routine
	size_t itsNbIteration;

	/// For Numerical Calibration Only
	ARM_DateStripPtr					itsNumSchedule;
	ARM_VanillaSecDensityPtrVector		itsNumSecDensities;
	
	/// ...
	void CleanUp();
    void CopyNoCleanUp(const ARM_CalibMethod& rhs );
    
    /// before setting the linked method... does some validations
	void Validate() const;
	vector<ARM_StdPortfolioPtr> DividePortfolio(ARM_StdPortfolioPtr);

    ARM_CalibDirection DefaultCalibDirection() const;

    bool itsIsCalibMethodShared;
    bool itsDoesValidation;

	string toString( const string& indent, bool printPreviousMethods ) const;

public:
	/// constructor with method type
	ARM_CalibMethod(ARM_StdPortfolioPtr Pf,
            const ARM_ModelParamVector&  CalibParams,
            ARM_MethodType methodType,
            const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER ,
            ARM_MktTargetType targetFuncType		= ARM_CalibrationTarget::PriceTarget,
            ARM_CalibMethod* linkedMethod			= NULL,
            ARM_CalibMethod* PreviousMethod			= NULL,
            bool isCalibMethodShared				= false,
			size_t FactorNb							= 0,
			size_t nbIter							= 1,
			bool validate							= true,
			const ARM_DateStripPtr&	numSchedule		= ARM_DateStripPtr(NULL),
			const ARM_VanillaSecDensityPtrVector& numSecDensities = ARM_VanillaSecDensityPtrVector(0) );
	
	/// constructor with model fitter descriptor
	ARM_CalibMethod(ARM_StdPortfolioPtr Pf,
            const ARM_ModelParamVector&  CalibParams,
			ARM_MethodType methodType,
            ARM_ModelFitterDes* modelFitterDes,
            ARM_MktTargetType targetFuncType		= ARM_CalibrationTarget::PriceTarget,
            ARM_CalibMethod* linkedMethod			= NULL,
            ARM_CalibMethod* PreviousMethod			= NULL,
            bool isCalibMethodShared				= false,
			size_t FactorNb							= 0,
			size_t nbIter							= 1,
			bool validate							= true,
			const ARM_DateStripPtr&	NumSchedule		= ARM_DateStripPtr(NULL),
			const ARM_VanillaSecDensityPtrVector& NumSecDensities = ARM_VanillaSecDensityPtrVector(0) );

	ARM_CalibMethod(const ARM_CalibMethod& rhs);
	virtual ~ARM_CalibMethod();
	ARM_CalibMethod& operator = (const ARM_CalibMethod& rhs);

    /// How many modelparmas
	int size() const { return itsCalibParams.size(); }

	/// does some validations with the model
	void DefaultValidateWithModel(const ARM_PricingModel& model) const;

    /// method to create the core of the calibration
    void InitialiseModelFitter( ARM_PricingModel* model );
    void Initialise(ARM_PricingModel* model);
    void Calibrate( ARM_PricingModel* pricingModel );
	
    /// accessors
    inline size_t CalibParamSize() const								{ return itsCalibParams.size();					}
    inline ARM_StdPortfolioPtr GetPortfolio() const						{ return itsPortfolio;							}
    inline void SetPortfolio(const ARM_StdPortfolioPtr& portfolio)		{itsPortfolio = portfolio;						}
    inline ARM_ModelParamVector& GetCalibParams()						{ return itsCalibParams;						}
    inline const ARM_ModelParamVector& GetCalibParams() const			{ return itsCalibParams;						}
    inline const ARM_ModelParam* GetCalibParam(int index = 0) const		{ return itsCalibParams[index];					}
    inline ARM_ModelParam* GetCalibParam(int index = 0)					{ return itsCalibParams[index];					}
    inline ARM_CalibDirection GetCalibDirection() const					{ return itsCalibDirection;						}
    inline ARM_ModelFitterPtr& GetModelFitter()							{ return itsModelFitter;						}
    inline bool GetIsCalibMethodShared() const							{ return itsIsCalibMethodShared;				}
    inline void SetIsCalibMethodShared(bool isCalibMethodShared )		{itsIsCalibMethodShared = isCalibMethodShared;	}
	inline const ARM_ModelFitterDes* const GetModelFitterDes() const	{ return itsModelFitterDes;						}
	inline ARM_ModelFitterDes* GetModelFitterDes()						{ return itsModelFitterDes;						}
	inline void SetCalibParam(ARM_ModelParam* calibParam )				{itsCalibParams[0] = calibParam;				}
	inline size_t GetFactorNb() const									{ return itsFactorNb;							}
	inline size_t GetNbIteration() const                                { return itsNbIteration;						} 
	inline void SetNbIteration(int nbIter)                              { itsNbIteration = nbIter;     					} 
	inline ARM_MethodType GetMethodType() const							{ return itsMethodType;							}
	inline void SetMethodType(ARM_MethodType methodType) 				{ itsMethodType = methodType;					}

	/// Get how long time the calibration needed
	double GetDuration() const;
   
	inline ARM_CalibMethod* GetlinkedMethod() const     { return itsLinkedMethod;   }
    void SetlinkedMethod(ARM_CalibMethod* method)       {itsLinkedMethod = method;  }
    inline ARM_CalibMethod* GetPreviousMethod() const   { return itsPreviousMethod; }  
    void SetPreviousMethod(ARM_CalibMethod* method)     {itsPreviousMethod = method;}
	inline ARM_CalibMethod* GetNextMethod() const		{ return itsNextMethod; }  
    void SetNextMethod(ARM_CalibMethod* method)			{itsNextMethod = method;}
    
    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual void View(char* id = NULL, FILE* ficOut = NULL) const;

	ARM_Warning* GetWarning() const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/