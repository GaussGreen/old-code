//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MRSpotVolProcessedBS.hpp
//
//   Description : What a mean reverting 'Spot' volatility can do.
//                 Here 'spot' means in the traditional IR sense
//
//   Author      : Mark A Robson
//
//   Date        : 14 Dec 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SPOT_VOL_PROCESSED_HPP
#define EDR_SPOT_VOL_PROCESSED_HPP
#include "edginc/VolProcessed.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/** What a 'Spot' volatility can do.
    Here 'spot' means in the traditional IR sense */
class MARKET_DLL MRSpotVolProcessed: public CObject,
                          public virtual IVolProcessed{
public:
    static CClassConstSP const TYPE; // in MRSpotVolRequest.cpp

    ~MRSpotVolProcessed(); // in MRSpotVolRequest.cpp
    /** What this vol mean reverts to */
    virtual double meanReversion() const = 0;

    /** Returns spot vols between pairs of dates ie between
     [initialStartDate, subsequentDates[0]], 
     [subsequentDates[0], subsequentDates[1], ... */
    virtual void spotVol(const DateTime&      initialStartDate,
                         const DateTimeArray& subsequentDates,
                         DoubleArray&         vol) const = 0;

    /** Returns the dates (possibly none!) that the spot vol is defined on */
    virtual DateTimeArraySP getSpotVolDates() const = 0;
protected:
    MRSpotVolProcessed(CClassConstSP clazz); // in MRSpotVolRequest.cpp
private:
    static void load(CClassSP& clazz);// in MRSpotVolRequest.cpp
    MRSpotVolProcessed(const MRSpotVolProcessed &rhs);
    MRSpotVolProcessed& operator=(const MRSpotVolProcessed& rhs);
    
};

DECLARE(MRSpotVolProcessed);

DRLIB_END_NAMESPACE
#endif
