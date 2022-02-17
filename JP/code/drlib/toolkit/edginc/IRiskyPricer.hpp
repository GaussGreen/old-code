//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : IRiskyPricer.hpp
//
//   Description : Interface that must be implemented if an Instrument needs to use
//                 a 'risky price' (i.e., a price assuming that equity is subject to 
//                 the risk of default) for a given sensitivity.
//                 Adapted from IRiskyPricer.cpp
//
//   Author      : Milan Kovacevic
//
//   Date        : 19th November 2002
//
//
//----------------------------------------------------------------------------
#ifndef I_RISKY_PRICER_HPP
#define I_RISKY_PRICER_HPP

#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL IRiskyPricer {
public:
    static CClassConstSP const TYPE;

    virtual ~IRiskyPricer();

    /** If true, set instument to price with risky equity, otherwise set it to price with risk free equity */
    virtual void setRisky(bool flag = true) = 0;

    /** Return true if the instrument is currently pricing assuming risky equity, otherwise false */
    virtual bool isRisky() const = 0;

    /** Returns the effective maturity to be used for building a risky curve */
    virtual DateTime getEffMaturityDate() const = 0;

    /** Set the effective maturity to be used for building a risky curve */
    virtual void setEffMaturityDate(const DateTime & maturityDate) = 0;

protected:
    IRiskyPricer();

};

DRLIB_END_NAMESPACE
#endif
