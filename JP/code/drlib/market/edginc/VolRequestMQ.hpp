//----------------------------------------------------------------------------
//   Group       : Credit QR
//
//   Filename    : VolRequestMQ.hpp
//
//   Description : Vol request MultiQ related vol
//
//   Author      : Keith Jia
//
//----------------------------------------------------------------------------

#ifndef _VOLREQUESTMQ_HPP
#define _VOLREQUESTMQ_HPP

#include "edginc/VolRequest.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(TimeMetric);
//FORWARD_DECLARE(Expiry);


/** Captures how to interpolate a volatility, i.e., it captures both
    instrument specific data needed to calculate volatilities/variance
    together with the methodology to be used for the interpolation.
    The methodology is captured by the type of the VolRequest used */
class MARKET_DLL VolRequestMQ: public CVolRequest {
public:
    static CClassConstSP const TYPE;
    virtual ~VolRequestMQ();
    VolRequestMQ(const CClassConstSP& clazz,
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
                 double fwdAnnuity);
    
    DateTime getResetDate() const;
    double getVolStart() const;
    double getReset() const;
    double getSwapStart() const;
    double getSwapMat() const;
    double getPayDate() const;
    int    getFreq() const;
    double getParFwd() const;
    double getFwdAnnuity() const;
    double getSwapMatZR() const;
    double getPayZR() const;
    ExpiryConstSP getTenor() const;

protected:

    VolRequestMQ(const CClassConstSP& clazz = TYPE);

private:
    //=========================================================================
    // DATA FIELDS
    //=========================================================================
    ExpiryConstSP tenor;
    TimeMetricConstSP metric;
    DateTime resetDate;
    DateTime swapStartDate;
    DateTime swapMatDate;
    DateTime payDate;
    double swapMatZeroRate;
    double payDateZeroRate;
    int freq;
    double parSwapRate;
    double fwdAnnuity;

    //transient
    double dVolStart;
    double dReset;
    double dSwapStart;
    double dSwapMat;
    double dPayDate;

private:
    static void load(CClassSP& clazz);

};

DECLARE(VolRequestMQ);

DRLIB_END_NAMESPACE

#endif
