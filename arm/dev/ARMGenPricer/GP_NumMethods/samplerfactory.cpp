/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file samplerfactory.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnummethods/samplerfactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"


/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/markoviandriftsampler.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SamplerFactoryData
///	Routine: toString
///	Returns: 
///	Action : Object dump
////////////////////////////////////////////////////
string ARM_SamplerFactoryData::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "Sampler Factory Data\n";
	os << indent << "--------------------\n";
	os << indent << "\nSampler=" << ARM_ArgConvReverse_SamplerType.GetString(itsSamplerType) << "\n";
    if(itsSamplerDatas.size()>1)
        os << indent << "   MinStdDev="<<itsSamplerDatas[0] << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_SamplerFactoryImp
///	Routine: CreateSampler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SamplerBase* ARM_SamplerFactoryImp::CreateSampler(const ARM_SamplerFactoryData& samplerFactoryData)
{
    return CreateSampler(
		samplerFactoryData.GetNbDims(),
		samplerFactoryData.GetSamplerType(),
		samplerFactoryData.GetSamplerDatas(),
		samplerFactoryData.GetScheduler());
}

////////////////////////////////////////////////////
///	Class  : ARM_SamplerFactoryImp
///	Routine: CreateSampler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SamplerBase* ARM_SamplerFactoryImp::CreateSampler(
		int nbDims,
        int samplerType, 
		const ARM_GP_Vector& samplerDatas,
		ARM_SchedulerBase* scheduler)
{
    ARM_SamplerBase* sampler            = NULL;

    /// Sampler creation
    size_t nbDatas = samplerDatas.size();
    switch(samplerType)
    {
    case ARM_SamplerBase::NormalCentred:
		if(nbDims == 1)
            sampler = new ARM_NormalCentredSampler1D(scheduler);
        else
		{
			if(nbDatas>0)
				sampler = new ARM_NormalCentredSamplerND(scheduler,samplerDatas[0]);
			else
				sampler = new ARM_NormalCentredSamplerND(scheduler);
		}
        break;

    case ARM_SamplerBase::MeanReverting:
        if(nbDims == 1)
        {
            if(nbDatas>0)
                sampler = new ARM_MeanRevertingSampler1D(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MeanRevertingSampler1D(scheduler);
        }
        else
        {
            if(nbDatas>1)
                sampler = new ARM_MeanRevertingSamplerND(scheduler,samplerDatas[0],samplerDatas[1]);
			else if(nbDatas>0)
				sampler = new ARM_MeanRevertingSamplerND(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MeanRevertingSamplerND(scheduler);
        }
        break;

    case ARM_SamplerBase::DriftedMeanReverting:
        if(nbDims == 1)
        {
            if(nbDatas>0)
                sampler = new ARM_DriftedMeanRevertingSampler1D(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_DriftedMeanRevertingSampler1D(scheduler);
        }
        else
        {
            if(nbDatas>0)
                sampler = new ARM_DriftedMeanRevertingSamplerND(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_DriftedMeanRevertingSamplerND(scheduler);
        }
        break;

    case ARM_SamplerBase::MarkovianDrift:
        if(nbDims==1)
        {
            if(nbDatas>0)
                sampler = new ARM_MarkovianDriftSampler1D(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MarkovianDriftSampler1D(scheduler);
        }
        
		else if(nbDims==2)
        {
            if(nbDatas>0)
                sampler = new ARM_MarkovianDriftSampler2D(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MarkovianDriftSampler2D(scheduler);
        }

		else if(nbDims==3)
        {
            if(nbDatas>0)
                sampler = new ARM_MarkovianDriftSampler3D(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MarkovianDriftSampler3D(scheduler);
        }

        else
        {
            if(nbDatas>0)
                sampler = new ARM_MarkovianDriftSamplerND(scheduler,samplerDatas[0]);
            else
                sampler = new ARM_MarkovianDriftSamplerND(scheduler);
        }
        break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown sampler type");
    }

	/// return the result
    return sampler;
}


ARM_SingletonHolder<ARM_SamplerFactoryImp> ARM_SamplerFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/