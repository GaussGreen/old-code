//----------------------------------------------------------------------------
//   Group       : Credit QR
//
//   Filename    : VolRequestMQ.cpp
//
//   Description : CMS Vol request for MQ Vol
//
//   Author      : Keith Jia
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VolRequestMQCMS.hpp"
#include "edginc/TimeMetric.hpp"


DRLIB_BEGIN_NAMESPACE

VolRequestMQCMS::~VolRequestMQCMS() 
{}

VolRequestMQCMS::VolRequestMQCMS(ExpiryConstSP tenor,
                                 TimeMetricConstSP metric,
                                 DateTime valueDate,
                                 DateTime resetDate,
                                 DateTime swapStartDate,
                                 DateTime swapMatDate,
                                 DateTime payDate,
                                 double swapMatZeroRate,
                                 double payDateZeroRate,
                                 int freq,
                                 double parSwapRate,
                                 double fwdAnnuity,
                                 bool isCashSettled,
                                 double beta,
                                 double RichardsonOrder)
    :
    VolRequestMQ(TYPE,
                 tenor,
                 metric,
                 valueDate,
                 resetDate,
                 swapStartDate,
                 swapMatDate,
                 payDate,
                 swapMatZeroRate,
                 payDateZeroRate,
                 freq,
                 parSwapRate,
                 fwdAnnuity),
    isCashSettled(isCashSettled),
    beta(beta),
    RichardsonOrder(RichardsonOrder)
{}

VolRequestMQCMS::VolRequestMQCMS(const CClassConstSP& clazz)
    : VolRequestMQ(clazz)
{}


bool VolRequestMQCMS::getIsCashSettled() const
{
    return isCashSettled;
}

double VolRequestMQCMS::getBeta() const
{
    return beta;
}

double VolRequestMQCMS::getRichardsonOrder() const
{
    return RichardsonOrder;
}


// static methods
void VolRequestMQCMS::load(CClassSP& clazz)
{
    REGISTER(VolRequestMQCMS, clazz);
    SUPERCLASS(VolRequestMQ);
    EMPTY_SHELL_METHOD(defaultVolR);
    FIELD(isCashSettled, "");
    FIELD(beta, "mean reversion");
    FIELD(RichardsonOrder, "integration order");
}

IObject* VolRequestMQCMS::defaultVolR()
{
    return new VolRequestMQCMS();
}

CClassConstSP const VolRequestMQCMS::TYPE = CClass::registerClassLoadMethod(
    "VolRequestMQCMS", typeid(VolRequestMQCMS), load);


DRLIB_END_NAMESPACE
