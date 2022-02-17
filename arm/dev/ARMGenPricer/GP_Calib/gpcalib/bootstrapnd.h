/*§
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file bootstrapnd.h
 *
 *  \brief 
 *
 *	\author  E.M Ezzine
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_BOOTSTRAPND_H
#define _INGPCALIB_BOOTSTRAPND_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "modelfitter.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_BootstrapND : public ARM_ModelFitter
{
private:

//	string* itsFileName;	/// corresponding fileName
    void CleanUp();
    void CopyNoCleanUp(const ARM_BootstrapND& rhs ){};

	/// standard methods
	void Validate();
    void Initialise();

public:
    ARM_BootstrapND( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& modelParams ,
        ARM_ModelFitterPtr& linkedModelFitter			= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter			= ARM_ModelFitterPtr(NULL ),
        const size_t max_iter							= ARM_NumericConstants::ARM_GP_MAX_ITER,
		bool getDetails									= false,
		size_t factorNb									= 0,
		size_t NbCalibIter								= 1,
        ARM_MktTargetType	                            = ARM_CalibrationTarget::PriceTarget );

    ARM_BootstrapND(const ARM_BootstrapND& rhs);
	ARM_BootstrapND& operator = (const ARM_BootstrapND& rhs);
    virtual ~ARM_BootstrapND();

    virtual void Calibrate();

    /// standard ARM Object support
    virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
