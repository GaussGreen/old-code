//----------------------------------------------------------------------------
//   Group       : Credit QR
//
//   Filename    : VolRequestMQCMS.hpp
//
//   Description : CMS Vol request for MultiQ related vol
//
//   Author      : Keith Jia
//
//----------------------------------------------------------------------------

#ifndef _VOLREQUESTMQCMS_HPP
#define _VOLREQUESTMQCMS_HPP

#include "edginc/VolRequestMQ.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** Captures how to interpolate a volatility, i.e., it captures both
    instrument specific data needed to calculate volatilities/variance
    together with the methodology to be used for the interpolation.
    The methodology is captured by the type of the VolRequest used */
class MARKET_DLL VolRequestMQCMS: public VolRequestMQ {
public:
    static CClassConstSP const TYPE;

    virtual ~VolRequestMQCMS();
    
    VolRequestMQCMS(ExpiryConstSP tenor,
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
                    double RichardsonOrder);

    bool getIsCashSettled() const;
    double getBeta() const;
    double getRichardsonOrder() const;

protected:

    VolRequestMQCMS(const CClassConstSP& clazz = TYPE);

private:
    //=========================================================================
    // DATA FIELDS
    //=========================================================================
    bool   isCashSettled;
    double beta;
    double RichardsonOrder;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultVolR();
};

DECLARE(VolRequestMQCMS);

DRLIB_END_NAMESPACE

#endif
