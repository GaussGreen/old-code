/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treefactory.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/treefactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"

/// gpinfra
#include "gpinfra/pricingstates.h"


/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/tree1d.h"
#include "gpnummethods/tree2d.h"
#include "gpnummethods/tree3d.h"
#include "gpnummethods/treend.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/samplerfactory.h"
#include "gpnummethods/slice.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/schedulerfactory.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/truncatorfactory.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"
#include "gpnummethods/argconvdefault.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_TreeFactoryData
///	Routine: toString
///	Returns: 
///	Action : Object dump
////////////////////////////////////////////////////
string ARM_TreeFactoryData::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

	os << indent << "Tree ND Factory Data\n";
	os << indent << "--------------------\n";
	os << indent << "\nScheduler=" << ARM_ArgConvReverse_SchedulerType.GetString(itsSchedulerType) << "\n";
    if(itsSchedulerDatas.size()>1)
        os << indent << "   NbSteps="<<itsSchedulerDatas[0] << "\n";
	os << indent << "\nSampler=" << ARM_ArgConvReverse_SamplerType.GetString(itsSamplerType) << "\n";
    if(itsSamplerDatas.size()>1)
        os << indent << "   MinStdDev="<<itsSamplerDatas[0] << "\n";
	os << indent << "\nTruncator=" << ARM_ArgConvReverse_TruncatorType.GetString(itsTruncatorType) << "\n";
    if(itsTruncatorDatas.size()>1)
        os << indent << "   NbStdDev="<<itsTruncatorDatas[0] << "\n";
	os << indent << "\nReconnector=" << ARM_ArgConvReverse_ReconnectorType.GetString(itsReconnectorType) << "\n";
    os << indent << "\nProbas Computation = " << (itsProbasFlag ? "On" : "Off") << "\n";

	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeFactoryImp
///	Routine: CreateTreeND
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeBase* ARM_TreeFactoryImp::CreateTreeND(size_t nbDims,const ARM_TreeFactoryData& treeFactoryData)
{
    return CreateTreeND(nbDims,
        treeFactoryData.GetSchedulerType(), treeFactoryData.GetSchedulerDatas(),
        treeFactoryData.GetSamplerType(), treeFactoryData.GetSamplerDatas(),
        treeFactoryData.GetTruncatorType(),treeFactoryData.GetTruncatorDatas(),
        treeFactoryData.GetProbasFlag(),
        treeFactoryData.GetReconnectorType(),treeFactoryData.GetSmootherType());
}

////////////////////////////////////////////////////
///	Class  : ARM_TreeFactoryImp
///	Routine: CreateTreeND
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_TreeBase* ARM_TreeFactoryImp::CreateTreeND(size_t nbDims,
        int schedulerType, const ARM_GP_Vector& schedulerDatas,
        int samplerType, const ARM_GP_Vector& samplerDatas,
        int truncatorType, const ARM_GP_Vector& truncatorDatas,
        bool probasFlag,
        int reconnectorType, int smootherType)
{
    /// Scheduler & sampler consistency
    if( (samplerType == ARM_SamplerBase::MeanReverting  || samplerType == ARM_SamplerBase::MarkovianDrift)
        && schedulerType == ARM_SchedulerBase::ConstantVariance)
    {
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": sampler & scheduler types are inconsistent" );
    }

    ARM_SamplerBase* sampler            = NULL;
    ARM_SchedulerBase* scheduler        = NULL;
    ARM_TruncatorBase* truncator        = NULL;
    ARM_ReconnectorBase* reconnector    = NULL;
    ARM_SmootherBase* smoother          = NULL;
	ARM_TreeBase* tree                  = NULL;

	int nbSteps;

	// Scheduler creation
	scheduler = ARM_SchedulerFactory.Instance()->CreateScheduler(
		schedulerType,
		schedulerDatas,
		nbSteps);

	// Sampler creation
	sampler = ARM_SamplerFactory.Instance()->CreateSampler(
		nbDims,
		samplerType,
		samplerDatas,
		scheduler);

    /// Truncator creation
	truncator = ARM_TruncatorFactory.Instance()->CreateTruncator(
		nbDims,
		truncatorType,
		truncatorDatas);


    /// Reconnector creation
    switch(reconnectorType)
    {
    case ARM_ReconnectorBase::DoNothing:
        reconnector = new ARM_ReconnectorDoNothing;
        break;

    case ARM_ReconnectorBase::Mean:
        reconnector = new ARM_ReconnectorMean;
        break;

    case ARM_ReconnectorBase::Variance:
        reconnector = new ARM_ReconnectorVar;
        break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown reconnector type");
    }
    truncator->SetReconnector(reconnector);

    /// Smoother creation
    switch(smootherType)
    {
    case ARM_SmootherBase::DoNothing:
        smoother = new ARM_SmootherDoNothing;
        break;

    case ARM_SmootherBase::Linear:
        smoother = new ARM_SmootherLinear;
        break;

    case ARM_SmootherBase::Quadratic:
        smoother = new ARM_SmootherQuadratic;
        break;

    case ARM_SmootherBase::Cubic:
        smoother = new ARM_SmootherCubic;
        break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown smoother type");
    }

    /// Tree creation
    switch(nbDims)
    {
    case 1:
        tree = new ARM_Tree1D(sampler,truncator,reconnector,smoother,probasFlag);
        break;

    case 2:
        tree = new ARM_Tree2D(sampler,truncator,reconnector,smoother,probasFlag);
        break;

    case 3:
        tree = new ARM_Tree3D(sampler,truncator,reconnector,smoother,probasFlag);
        break;

    default:
        tree = new ARM_TreeND(sampler,truncator,reconnector,nbDims,smoother,probasFlag);
        break;
    }

    tree->SetNbSteps(nbSteps);

	/// to avoid memory leak, delete objects
	delete sampler;
	delete scheduler;
	delete truncator;
	delete reconnector;
	delete smoother;

	/// return the result
    return tree;
}


ARM_SingletonHolder<ARM_TreeFactoryImp> ARM_TreeFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/