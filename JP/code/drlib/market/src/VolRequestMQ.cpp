//----------------------------------------------------------------------------
//   Group       : Credit QR
//
//   Filename    : VolRequestMQ.cpp
//
//   Description : Vol request for MQ Vol
//
//   Author      : Keith Jia
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Class.hpp"
#include "edginc/VolRequestMQ.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/TimeMetric.hpp"


DRLIB_BEGIN_NAMESPACE

VolRequestMQ::~VolRequestMQ() 
{}

VolRequestMQ::VolRequestMQ(const CClassConstSP& clazz,
                           ExpiryConstSP tenor,
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
                           double fwdAnnuity)
    :
    CVolRequest(clazz), tenor(tenor), metric(metric), resetDate(resetDate), 
    swapStartDate(swapStartDate),
    payDate(payDate), swapMatZeroRate(swapMatZeroRate), 
    payDateZeroRate(payDateZeroRate), freq(freq), parSwapRate(parSwapRate), 
    fwdAnnuity(fwdAnnuity)
{
    dVolStart = 0;
    dReset = metric->yearFrac(valueDate, resetDate);
    dSwapStart = metric->yearFrac(valueDate, swapStartDate);
    dSwapMat = metric->yearFrac(swapStartDate, swapMatDate);
    dPayDate = metric->yearFrac(swapStartDate, payDate);
}


DateTime VolRequestMQ::getResetDate() const
{
    return resetDate;
}

double VolRequestMQ::getVolStart() const
{
    return dVolStart;
}

double VolRequestMQ::getReset() const
{
    return dReset;
}

double VolRequestMQ::getSwapStart() const
{
    return dSwapStart;
}

double VolRequestMQ::getSwapMat() const
{
    return dSwapMat;
}

double VolRequestMQ::getPayDate() const
{
    return dPayDate;
}

int    VolRequestMQ::getFreq() const
{
    return freq;
}

double VolRequestMQ::getParFwd() const
{
    return parSwapRate;
}

double VolRequestMQ::getFwdAnnuity() const
{
    return fwdAnnuity;
}


double VolRequestMQ::getSwapMatZR() const
{
    return swapMatZeroRate;
}

double VolRequestMQ::getPayZR() const
{
    return payDateZeroRate;
}

ExpiryConstSP VolRequestMQ::getTenor() const
{
    return tenor;
}


VolRequestMQ::VolRequestMQ(const CClassConstSP& clazz)
    : CVolRequest(clazz)
{}



void VolRequestMQ::load(CClassSP& clazz)
{
    REGISTER(VolRequestMQ, clazz);
    SUPERCLASS(CVolRequest);

    FIELD(metric, "time metric for year fraction calculation");
    FIELD(resetDate, "option expiry or reset date");
    FIELD(swapStartDate, "");
    FIELD(payDate, "");
    FIELD(swapMatZeroRate, "zero rate at swap maturity");
    FIELD(payDateZeroRate, "zero rate at payment date");
    FIELD(freq, "swap frequency");
    FIELD(parSwapRate, "fwd par swap rate");
    FIELD(fwdAnnuity, "fwd annuity");

    //transient
    FIELD(dVolStart, "");
    FIELD_MAKE_TRANSIENT(dVolStart);
    FIELD(dReset, "");
    FIELD_MAKE_TRANSIENT(dReset);
    FIELD(dSwapStart, "");
    FIELD_MAKE_TRANSIENT(dSwapStart);
    FIELD(dSwapMat, "");
    FIELD_MAKE_TRANSIENT(dSwapMat);
    FIELD(dPayDate, "");
    FIELD_MAKE_TRANSIENT(dPayDate);
}

CClassConstSP const VolRequestMQ::TYPE = CClass::registerClassLoadMethod(
    "VolRequestMQ", typeid(VolRequestMQ), load);


DRLIB_END_NAMESPACE
