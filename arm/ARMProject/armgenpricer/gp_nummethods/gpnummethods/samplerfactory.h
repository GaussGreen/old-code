/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file samplerfactory.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPNUMMETHODS_SAMPLERFACTORY_H
#define _INGPNUMMETHODS_SAMPLERFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
struct ARM_SamplerBase;
struct ARM_SchedulerBase;

class ARM_SamplerFactoryData : public ARM_RootObject
{
private:
	int nbDims;
    int itsSamplerType;
    std::vector<double> itsSamplerDatas;
	ARM_SchedulerBase* itsScheduler;

public:
	enum {OneDim, MultiDim} DimType;


    ARM_SamplerFactoryData(int samplerType, const std::vector<double>& samplerDatas)
        : itsSamplerType(samplerType), itsSamplerDatas(samplerDatas)
    {}

    virtual ~ARM_SamplerFactoryData() {}
	virtual ARM_Object* Clone() const { return new ARM_SamplerFactoryData(*this); }

    virtual string toString(const string& indent="",const string& nextIndent="") const;

	int GetNbDims() const { return nbDims; };
    int GetSamplerType() const { return itsSamplerType; }
    const std::vector<double>& GetSamplerDatas() const { return itsSamplerDatas; }
	ARM_SchedulerBase* GetScheduler() const { return itsScheduler; };
};


struct ARM_SamplerFactoryImp
{
    ARM_SamplerBase* CreateSampler(
		int nbDims,
		int samplerType,
		const std::vector<double>& samplerDatas,
		ARM_SchedulerBase* scheduler);

    ARM_SamplerBase* CreateSampler(const ARM_SamplerFactoryData& samplerData);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_SamplerFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_SamplerFactoryImp>;
};

extern ARM_SingletonHolder<ARM_SamplerFactoryImp> ARM_SamplerFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/