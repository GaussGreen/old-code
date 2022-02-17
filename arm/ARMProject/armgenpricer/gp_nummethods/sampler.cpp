/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file sampler.cpp
 *
 *  \brief sampler object are used to discretise in time according to 
 *      various hypotheses
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date November 2004
 */

#include "gpnummethods/sampler.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_SamplerBase
///	Routine:
///	Returns:
///	Action : downcast to a 1D or ND sampler Base
////////////////////////////////////////////////////
ARM_Sampler1DBase* ARM_SamplerBase::ToSampler1DBase()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Sampler1DBase" );
}

const ARM_Sampler1DBase* ARM_SamplerBase::ToSampler1DBase() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Sampler1DBase" );
}

ARM_SamplerNDBase* ARM_SamplerBase::ToSamplerNDBase()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SamplerNDBase" );
}

const ARM_SamplerNDBase* ARM_SamplerBase::ToSamplerNDBase() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to SamplerNDBase" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SamplerBase
///	Routine: GetGlobalVCV
///	Returns: ARM_MatrixVector
///	Action : returns for each time step the global VCV matrix from
///          asOfDate to current time step
////////////////////////////////////////////////////
const ARM_MatrixVector& ARM_SamplerBase::GetGlobalVCV() const
{
    const const ARM_GP_Vector& timeSteps = *(itsModel->GetNumMethod()->GetTimeSteps());
    if(itsGlobalVCV.size() != timeSteps.size())
    {
        /// Free memory then rebuild global VCV matrixes
        if(itsGlobalVCV.size() > 0)
	        DeletePointorVector<ARM_GP_Matrix>(itsGlobalVCV);

        ARM_MatrixVector localVCV;
        itsModel->NumMethodStateLocalGlobalVariances(timeSteps,localVCV,itsGlobalVCV);
	    DeletePointorVector<ARM_GP_Matrix>(localVCV);
    }
    return itsGlobalVCV;        
}


CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

