/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file impsamplerfactory.h
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPNUMMETHODS_IMPSAMPLERFACTORY_H
#define _INGPNUMMETHODS_IMPSAMPLERFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_ImpSampler;

class ARM_ImpSamplerFactoryData : public ARM_RootObject
{
private:
    int itsImpSamplerType;
    const ARM_MultiCurve* itsAlpha;

public:
    ARM_ImpSamplerFactoryData(int samplerType, const ARM_MultiCurve* alpha)
        : itsImpSamplerType(samplerType), itsAlpha(alpha)
    {}

    virtual ~ARM_ImpSamplerFactoryData() {}
	virtual ARM_Object* Clone() const { return new ARM_ImpSamplerFactoryData(*this); }

    virtual string toString(const string& indent="",const string& nextIndent="") const;

    int GetImpSamplerType() const { return itsImpSamplerType; }
    const const ARM_MultiCurve* GetAlpha() const { return itsAlpha; }
};


struct ARM_ImpSamplerFactoryImp
{
    ARM_ImpSampler* CreateImpSampler(
		int impSamplerType,
		const ARM_MultiCurve* alpha);

    ARM_ImpSampler* CreateImpSampler(const ARM_ImpSamplerFactoryData& impSamplerData);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_ImpSamplerFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_ImpSamplerFactoryImp>;
};

extern ARM_SingletonHolder<ARM_ImpSamplerFactoryImp> ARM_ImpSamplerFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

