/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file nummethod.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#include "gpbase/ostringstream.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/gpmatrixtriangular.h"

#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"

#include <iomanip>
CC_USING_NS(std,setw);
CC_USING_NS(std,fixed);
CC_USING_NS(std,setprecision);

/// numerical method should evaluate in reverse order
/// in Monte Carlo should go last row to first row
/// in backward trees/PDES should go first row to last row
/// recursive call will make sure to call appropriate methods

/// order of backwardisation/MC is the following
/// initialisation
/// browse ParseTree and execute accordingly
/// tree should have only one cashflow!


#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

const bool isTree2G = true;

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Const  : GP_PricingDirectionTxt
/// Meaning: text corresponding to the enum
////////////////////////////////////////////////////

const string ARM_NumMethod::GP_PricingDirectionTxt[] =
{
	"Ambiguous",
	"Forward Looking",
	"Backward Looking",
	"Forward Backward Looking"
};


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_NumMethod::CopyNoCleanUp(const ARM_NumMethod& rhs)
{
    itsNbSteps					= rhs.itsNbSteps;
    itsLastTimeIdx				= rhs.itsLastTimeIdx;
    itsTimeSteps				= rhs.itsTimeSteps ? (ARM_GP_Vector*) rhs.itsTimeSteps->Clone() : NULL;
	itsSampler					= rhs.itsSampler ? static_cast<ARM_SamplerBase*>(rhs.itsSampler->Clone()) : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_NumMethod::CleanUp()
{
    delete itsTimeSteps;
	delete itsSampler;
}



////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_NumMethod::ARM_NumMethod(const ARM_SamplerBase* sampler)
:	
	itsNbSteps(0), 
	itsLastTimeIdx(0), 
	itsTimeSteps(NULL),
	itsSampler(sampler?static_cast<ARM_SamplerBase*>(sampler->Clone()):NULL)
{}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_NumMethod::ARM_NumMethod( const ARM_NumMethod& rhs )
:	ARM_RootObject( rhs )	
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////

ARM_NumMethod::~ARM_NumMethod()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: operator =
///	Returns: itself
///	Action : Affectation
////////////////////////////////////////////////////

ARM_NumMethod& ARM_NumMethod::operator=( const ARM_NumMethod& rhs)
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );
		CleanUp();
		CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: SetSampler
///	Returns:
///	Action : Set the sampler
////////////////////////////////////////////////////
void ARM_NumMethod::SetSampler(const ARM_SamplerBase* const sampler)
{
    itsSampler = sampler ? static_cast< ARM_SamplerBase* >(sampler->Clone()) : NULL;
}


////////////////////////////////////////////////////
///	Class   : ARM_NumMethod
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_NumMethod::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

    os << CC_NS(std,setprecision)(4) << "\n\n";
    
    os << indent << "Nb Steps : " << itsNbSteps << "\n";
    os << indent << "Diffusion Schedule : \n";
    if(itsTimeSteps)
    {
        os << indent << "Last Time Idx : " << itsLastTimeIdx;
        os << "   Time lags : \n";

        for(int i=0;i<itsTimeSteps->size();++i)
            os << indent << "t=" << fixed << setprecision(2) << setw(8) << (*itsTimeSteps)[i] << "\n";
    }
    else
        os << indent << "not initialised !\n";

	if( itsSampler)
        os << "Sampler    :  " << itsSampler->toString( indent, nextIndent ) << "\n";

    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: SetTimeSteps 
///	Returns: 
///	Action : Set the diffusion schedule
////////////////////////////////////////////////////
void ARM_NumMethod::SetTimeSteps(const ARM_GP_Vector& timeSteps)
{
	if (&timeSteps != itsTimeSteps)
	{
		delete itsTimeSteps;
		itsTimeSteps = (ARM_GP_Vector*) (const_cast<ARM_GP_Vector&>(timeSteps)).Clone();
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: ConvertEvalDate
///	Returns: void
///	Action : Return the evaluation date
////////////////////////////////////////////////////
double ARM_NumMethod::ConvertEvalDate(double evalDate, double asOfDate) const
{
	return evalDate;
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: GetNumMethodStateGlobalVars 
///	Returns: ARM_MatrixVector
///	Action : returns for each time step the global VCV matrix from
///          asOfDate to current time step
////////////////////////////////////////////////////
const ARM_MatrixVector& ARM_NumMethod::GetNumMethodStateGlobalVars() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented method GetNumMethodStateGlobalVars" );
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: GetSpotProbabilities 
///	Returns: ARM_GP_VectorPtr
///	Action : returns spot probabilities at the slice timeIdx
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_NumMethod::GetSpotProbabilities(size_t timeIdx) const
{
	return ARM_GP_VectorPtr(NULL);
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: GetArrowDebreuPrices 
///	Returns: ARM_GP_VectorPtr
///	Action : returns the Arrow Debreu prices at the slice timeIdx
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_NumMethod::GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented method GetArrowDebreuPrices" );
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: ProcessPaidPayoffs 
///	Returns: 
///	Action : Those functions are used to apply 
/// Importance Sampling
////////////////////////////////////////////////////

void ARM_NumMethod::ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{
	// By Default it does nothing
}

////////////////////////////////////////////////////
///	Class  : ARM_NumMethod
///	Routine: ProcessUnPaidPayoffs
///	Returns: 
///	Action : Those functions are used to apply 
/// Importance Sampling
////////////////////////////////////////////////////

void ARM_NumMethod::ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const
{
	// By Default it does nothing
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

