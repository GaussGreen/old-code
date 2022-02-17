//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IQMCStateVariableBase.hpp
//
//   Description : An interface to all the state variables produced via IQMCDiffusibleAsset
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IQMCSTATEVARIABLEBASE_HPP
#define EDR_IQMCSTATEVARIABLEBASE_HPP

#include "edginc/SVPath.hpp"
#include "edginc/ISVBase.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/TemplateIdx.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE


class IQMCMomentMatched
{
public:
    virtual ~IQMCMomentMatched(){}
    // if you need Moment matching - implement these three functions in the concrete SV classes
    virtual void resetMMCorrection() {}
    virtual void accumulateMMCorrection() {}
    virtual void setMMCorrection() {}
};


class MCARLO_DLL IQMCStateVariableBase : public virtual ISVBase, public IQMCMomentMatched, public SVPath
{
public:
    IQMCStateVariableBase() {}
    virtual const SVPath& path() const { return *this; }
    virtual void prepare(bool mm) = 0;
};
DECLARE(IQMCStateVariableBase);


class MCARLO_DLL IQMCStateVariableSpot : public IQMCStateVariableBase {

public:

    IQMCStateVariableSpot(const DateTimeArray& _measurementDates) :
        measurementDates(_measurementDates),
        measurementDatesIdx() 
        {
            // Because of the  pastPathGen we cannot assert that
            // ASSERT(!measurementDates.empty()); // no fake SVs
        }

    const DateTimeArray& getMeasurementDates() const {return measurementDates;}
protected:

    const DateTimeArray&   measurementDates;
    vector<SpotIdx> measurementDatesIdx;
};
DECLARE(IQMCStateVariableSpot);



class MCARLO_DLL IQMCStateVariableExpected : public IQMCStateVariableBase {

public:

    IQMCStateVariableExpected(const DateTime &             _measurementDate,
                       const DateTimeArray &        _futureDates,
                       bool                         _doLog) :
        measurementDate(_measurementDate),
        futureDates(_futureDates),
        doLog(_doLog) 
        {
            // Because of the  pastPathGen we cannot assert that
            //ASSERT(!futureDates.empty());   // no fake SVs
        }

    const DateTime& getMeasurementDate() const {return measurementDate;}
    const DateTimeArray& getFutureDates() const {return futureDates;}

protected:

    DateTime measurementDate;
    FwdIdx   measurementDateIdx;
    const DateTimeArray&  futureDates;   // FIXME: use reference (or better an SP) to the array
    vector<FwdIdx> futureDatesIdx;

    bool     doLog;
    // something for caching?
};
DECLARE(IQMCStateVariableExpected);

DRLIB_END_NAMESPACE
#endif // EDR_IQMCSTATEVARIABLEBASE_HPP

