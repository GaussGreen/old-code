/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file impsamplerfactory.cpp
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnummethods/impsamplerfactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"

/// gpnummethods
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/impsampler.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ImpSamplerFactoryData
///	Routine: toString
///	Returns: 
///	Action : Object dump
////////////////////////////////////////////////////
string ARM_ImpSamplerFactoryData::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "Importance Sampler Factory Data\n";
	os << indent << "--------------------\n";
	os << indent << "\nSampler=" << ARM_ArgConvReverse_ImpSamplerType.GetString(itsImpSamplerType) << "\n";
	os << indent << "\nalphar=" << itsAlpha->toString() << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_ImpSamplerImp
///	Routine: CreateImpSampler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ImpSampler* ARM_ImpSamplerFactoryImp::CreateImpSampler(const ARM_ImpSamplerFactoryData& impSamplerFactoryData)
{
    return CreateImpSampler(
		impSamplerFactoryData.GetImpSamplerType(),
		impSamplerFactoryData.GetAlpha());
}

////////////////////////////////////////////////////
///	Class  : ARM_SamplerFactoryImp
///	Routine: CreateImpSampler
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_ImpSampler* ARM_ImpSamplerFactoryImp::CreateImpSampler(
        int impSamplerType, 
		const ARM_MultiCurve* alpha)
{
    ARM_ImpSampler* impSampler            = NULL;

    /// Sampler creation
    switch(impSamplerType)
    {
	case ARM_ImpSampler::DummyImpSampler:
		impSampler = new ARM_DummyImpSampler();
		break;
    case ARM_ImpSampler::PropImpSampler:
		impSampler = new ARM_PropImpSampler(alpha);
		break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown importance sampler type");
    }

	/// return the result
    return impSampler;
}


ARM_SingletonHolder<ARM_ImpSamplerFactoryImp> ARM_ImpSamplerFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/