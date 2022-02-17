/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file bootstrapnd.cpp
 *  \brief bootstrap n dimensional model fitter
 *
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpcalib/bootstrapnd.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"

/// kernel
//#include <inst/portfolio.h>

/// gpcalib
#include "gpcalib/optimise.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/targetfunc.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: copy constructor
///	Returns: a built object
///	Action : builds the objet
////////////////////////////////////////////////////
ARM_BootstrapND::ARM_BootstrapND(const ARM_BootstrapND& rhs)
:   ARM_ModelFitter(rhs)
{
   CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_BootstrapND:: ~ARM_BootstrapND()
{
    CleanUp();
}

////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Assignment operator
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_BootstrapND& ARM_BootstrapND::operator = (const ARM_BootstrapND& rhs)
{
    if( this != &rhs )
	{
		ARM_Object::operator=( rhs );
		/// Delete the old ones.
		CleanUp();
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: CleanUp
///	Returns: 
///	Action : deletes any dynamic object defined at the object level
////////////////////////////////////////////////////
void ARM_BootstrapND::CleanUp( )
{
}


///////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Clone
///	Returns: 
///	Action : clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_BootstrapND::Clone() const
{
	return new ARM_BootstrapND(*this); 
}

////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: constructor
///	Returns: 
///	Action : constructor
////////////////////////////////////////////////////
ARM_BootstrapND::ARM_BootstrapND( ARM_PricingModel* model, 
    const ARM_StdPortfolioPtr portfolio, 
    const ARM_ModelParamVector& modelParam,
    ARM_ModelFitterPtr& linkedModelFitter,
    ARM_ModelFitterPtr& previousModelFitter,
    size_t max_iter,
	bool getDetails,
	size_t factorNb,
	size_t NbCalibIter,
    ARM_MktTargetType  targetType)
:   
	ARM_ModelFitter(model,
		portfolio,
		modelParam,
		linkedModelFitter, 
		previousModelFitter, 
		max_iter,
		getDetails,
		factorNb,
		NbCalibIter, 
		targetType)
{
    Validate();
    Initialise();
};


////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Validate
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_BootstrapND::Validate()
{
    if(GetPortfolio() == ARM_StdPortfolioPtr(NULL))
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"optimise calibration needs non null portfolio!" );

    if( GetCalibParams().empty() )
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"optimise calibration needs at least one calib param!" ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Initialise
///	Returns: 
///	Action : ....
////////////////////////////////////////////////////
void ARM_BootstrapND::Initialise()
{
	/*int size = GetCalibParams().size();
	GetPortfolio()->sort();
	for(int i= 0; i<size; ++i)
	{
		GetPricingModel()->AdviseBreakPointTimes( GetPortfolio(), GetCalibParam(i), GetCalibParam(i)->GetAdviseBreakPointTimes() );
		GetPricingModel()->GetModelParams()->MergeModelParam(GetCalibParam(i));
	}*/
}


////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: Calibrate
///	Returns: 
///	Action :  call nag optimiser and calibrate
////////////////////////////////////////////////////
void ARM_BootstrapND::Calibrate()
{
    /*int i,k;
    int cpsize = GetCalibParams().size();    
    /// Uses the local object to use optimise routine
    
    vector<ARM_StdPortfolio*>* portfolios = GetPortfolio()->BuiltPortfolioByMaturity();
    int size = portfolios->size();
    for(k = 1; k<size; ++k)
    { 
        ARM_ModelParamVector calibparams(cpsize);
        for(i=0; i< cpsize; ++i)
		{
			ARM_CurveModelParam* curveCalibParam = dynamic_cast<ARM_CurveModelParam*>(GetCalibParam(i));
			if( curveCalibParam )
				calibparams[i] = curveCalibParam->GetCalibParam(k,k+1);
			else
				throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "type not supported!" );
		}

        ARM_Optimise optimize(GetPricingModel(),ARM_StdPortfolioPtr((*portfolios)[k]),calibparams,
                        ARM_ModelFitterPtr(NULL), ARM_ModelFitterPtr(NULL),GetMax_Iter());
		optimize.DoCalibrationProcess();
    } */
}

////////////////////////////////////////////////////
///	Class  : ARM_BootstrapND
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_BootstrapND::toString( const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << indent << "=======> BOOTSTRAP_ND FITTER <====== \n\n";
    os << indent << "Corresponding errors\n";

    /// for unix compatibility
    /// we use printf
    /// more traditional C++ should use std::iostream manipulation
    char temp1[5];
    char temp2[30];
    SetUpError();
    for(int j=0; j<GetError()->GetRowsNb(); ++j)
    {
        os << indent << "[";
        sprintf( temp1, "%2d", j );
        os << temp1 << "]=  ";
		sprintf( temp2, "%+6.6E (=%+4.4f %s)", (*GetError())(j,0), (*GetError())(j,1)*100.0, "%" );
        os << temp2 << "\n";
    }

	os << "\n";
	os << indent << "===========> End Of Infos    <=========\n";
    return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

