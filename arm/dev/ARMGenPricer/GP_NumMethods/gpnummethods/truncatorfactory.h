/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file trucatorfactory.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPNUMMETHODS_TRUNCATORFACTORY_H
#define _INGPNUMMETHODS_TRUNCATORFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
struct ARM_TruncatorBase;

class ARM_TruncatorFactoryData : public ARM_RootObject
{
private:
	int itsNbDims;
	int itsTruncatorType;
    ARM_GP_Vector itsTruncatorDatas;

public:
    ARM_TruncatorFactoryData(int truncatorType, const ARM_GP_Vector& truncatorDatas)
        : itsTruncatorType(truncatorType), itsTruncatorDatas(truncatorDatas)
    {}

    virtual ~ARM_TruncatorFactoryData() {}
	virtual ARM_Object* Clone() const { return new ARM_TruncatorFactoryData(*this); }

    virtual string toString(const string& indent="",const string& nextIndent="") const;

	int GetNbDims() const { return itsNbDims; }
    int GetTruncatorType() const { return itsTruncatorType; }
    const ARM_GP_Vector& GetTruncatorDatas() const { return itsTruncatorDatas; }
};


struct ARM_TruncatorFactoryImp
{
    ARM_TruncatorBase* CreateTruncator(
		int nbDims,
		int truncatorType,
		const ARM_GP_Vector& truncatorDatas);

    ARM_TruncatorBase* CreateTruncator(const ARM_TruncatorFactoryData& truncatorData);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_TruncatorFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_TruncatorFactoryImp>;
};

extern ARM_SingletonHolder<ARM_TruncatorFactoryImp> ARM_TruncatorFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/