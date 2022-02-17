/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file samplerbase.h
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPNUMMETHODS_SAMPLERBASE_H
#define _INGPNUMMETHODS_SAMPLERBASE_H

#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "typedef.h"
#include "schedulerbase.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;

//// structure returned by a process Sampler
struct ARM_TimeStepsAndSlices;

/// forward declaration
struct ARM_Sampler1DBase;
struct ARM_SamplerNDBase;

/////////////////////////////////////////////
/// abstract class for all the process sampler
/// with the defaul method for ComputeTimeSteps
/// put a minimum and maximum nb of steps before 1st event data
/////////////////////////////////////////////
struct ARM_SamplerBase : public ARM_RootObject
{
	/// constructor, copy constructor, assignment operator, destructor
	ARM_SamplerBase( const ARM_SchedulerBasePtr& scheduler = ARM_SchedulerBasePtr(NULL) )
	:	itsModel(NULL), itsScheduler(scheduler), itsGlobalVCV(0) {}
	ARM_SamplerBase( const ARM_SamplerBase& rhs )
	:	ARM_RootObject(rhs), itsModel(rhs.itsModel), itsScheduler(rhs.itsScheduler), itsGlobalVCV(0)
    { DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsGlobalVCV, itsGlobalVCV );}
	ARM_SamplerBase& operator=(const ARM_SamplerBase& rhs )
	{
		if( this != &rhs )
		{
			ARM_RootObject::operator=(rhs);
			itsModel	 = rhs.itsModel;
			itsScheduler = rhs.itsScheduler;
			DeletePointorVector<ARM_GP_Matrix>(itsGlobalVCV);
            DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsGlobalVCV, itsGlobalVCV );
		}
		return *this;
	}
    virtual ~ARM_SamplerBase() { DeletePointorVector<ARM_GP_Matrix>(itsGlobalVCV); }

    enum SamplerType
	{
		NormalCentred = 0,
		MeanReverting,
		DriftedMeanReverting,
        MarkovianDrift
    };

    /// Downcasts
    virtual ARM_Sampler1DBase* ToSampler1DBase();
    virtual const ARM_Sampler1DBase* ToSampler1DBase() const;
    virtual ARM_SamplerNDBase* ToSamplerNDBase();
    virtual const ARM_SamplerNDBase* ToSamplerNDBase() const;

	/// 1D equivalent for hybrid model
	virtual ARM_SamplerBase* CorrespondingSampler1D(size_t dim=0) const = 0;

    /// Accessor
    void SetModel(const ARM_PricingModel* mod ) { itsModel=mod;}
    const ARM_PricingModel* const GetModel() const { return itsModel;}

	ARM_SchedulerBasePtr GetScheduler() const{ return itsScheduler; }
	void SetScheduler( const ARM_SchedulerBasePtr& scheduler ) { itsScheduler = scheduler; }

    virtual bool IsIntegratedSampling() const { return true; };

    virtual ARM_TimeStepsAndSlices* Init( const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices=true) = 0;

    const ARM_MatrixVector& GetGlobalVCV() const;
    void SetGlobalVCV(const ARM_MatrixVector& globalVCV) { DeletePointorVector<ARM_GP_Matrix>(itsGlobalVCV);itsGlobalVCV=globalVCV; }

    /// X to Z space global converter
    virtual void ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const = 0;

private:
    const ARM_PricingModel* itsModel; /// set but not cloned hence never deleted
	ARM_SchedulerBasePtr itsScheduler;

    mutable ARM_MatrixVector itsGlobalVCV;
};

CC_END_NAMESPACE()

#endif