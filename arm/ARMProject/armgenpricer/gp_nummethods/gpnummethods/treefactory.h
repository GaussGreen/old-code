/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treefactory.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREEFACTORY_H
#define _INGPNUMMETHODS_TREEFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_TreeBase;

class ARM_TreeFactoryData : public ARM_RootObject
{
private:
    int itsSchedulerType;
    std::vector<double> itsSchedulerDatas;
    int itsSamplerType;
    std::vector<double> itsSamplerDatas;
    int itsTruncatorType;
    std::vector<double> itsTruncatorDatas;
    bool itsProbasFlag;
    int itsReconnectorType;
    int itsSmootherType;

public:
    ARM_TreeFactoryData(int schedulerType, const std::vector<double>& schedulerDatas,
        int samplerType, const std::vector<double>& samplerDatas,
        int truncatorType, const std::vector<double>& truncatorDatas, bool probasFlag,
        int reconnectorType, int smootherType)
        : itsSchedulerType(schedulerType), itsSchedulerDatas(schedulerDatas),
          itsSamplerType(samplerType), itsSamplerDatas(samplerDatas),
          itsTruncatorType(truncatorType), itsTruncatorDatas(truncatorDatas),
          itsProbasFlag(probasFlag),
          itsReconnectorType(reconnectorType), itsSmootherType(smootherType)
    {}

    virtual ~ARM_TreeFactoryData() {}
	virtual ARM_Object* Clone() const { return new ARM_TreeFactoryData(*this); }

    virtual string toString(const string& indent="",const string& nextIndent="") const;

    int GetSchedulerType() const { return itsSchedulerType; }
    const std::vector<double>& GetSchedulerDatas() const { return itsSchedulerDatas; }
    int GetSamplerType() const { return itsSamplerType; }
    const std::vector<double>& GetSamplerDatas() const { return itsSamplerDatas; }
    int GetTruncatorType() const { return itsTruncatorType; }
    const std::vector<double>& GetTruncatorDatas() const { return itsTruncatorDatas; }
    bool GetProbasFlag() const { return itsProbasFlag; }
    void SetProbasFlag(bool probaFlag) { itsProbasFlag=probaFlag; }
    int GetReconnectorType() const { return itsReconnectorType; }
    int GetSmootherType() const { return itsSmootherType; }
};



struct ARM_TreeFactoryImp
{
    ARM_TreeBase* CreateTreeND(size_t nbDims,
        int schedulerType, const std::vector<double>& schedulerDatas,
        int samplerType, const std::vector<double>& samplerDatas,
        int truncatorType, const std::vector<double>& truncatorDatas,
        bool probasFlag, int reconnectorType, int smootherType);

    ARM_TreeBase* CreateTreeND(size_t nbDims,const ARM_TreeFactoryData& treeFactoryData);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_TreeFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_TreeFactoryImp>;
};

extern ARM_SingletonHolder<ARM_TreeFactoryImp> ARM_TreeFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

