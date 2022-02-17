//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CliquetVolRequest.hpp
//
//   Description : Cliquet vol request 
//
//   Author      : Andre Segger
//
//   Date        : 02 Nov 2001
//
//
//----------------------------------------------------------------------------

#ifndef CLIQUET_VOLREQUEST_HPP
#define CLIQUET_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequestLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** The cliquet is defined for the purposes of volatility
    interpolation by a list of date-times, the live cliquet start
    dates, and interpolation levels. The term live indicates that all
    provided cliquet start dates must be in the future, with the
    possible exception of the first.  In this special case, the
    corresponding interpolation level must be supplied as an absolute
    level, otherwise levels are relative (%).

    For a detailed description of the methodology used, see the
    document "Methodology for Cliquet Volatility Interpolation" in the
    EDR Research database */
class MARKET_DLL CliquetVolRequest: public CVolRequestLN {
public:
    static CClassConstSP const TYPE;
    friend class CliquetVolRequestHelper;

    /** firstCliquetInFuture means if futStrikeDates[0] > value date.
        interpLevels are floored at 0 (ie negative strikes set to 0 - makes
        vega matrix nicer) */
    CliquetVolRequest(bool                 firstCliquetInFuture,
                      const DateTimeArray& futStrikeDates,
                      const DateTime&      endDate,
                      const DoubleArray&   interpLevels);

    /** Returns the end date for vol interpolation */
    const DateTime& getEndDate() const;

    /** Scale the interpolation object */
    virtual void scale(double scaleFactor); // multiplicative scaling

    /** Returns the start date for forward starting volatility 
        interpolation. */
    virtual DateTime getStartDate() const;

    const DateTimeArray& getCliqStartDates() const;

    /** The expectPercFor1stCliq must match the firstCliquetInFuture flag */
    const DoubleArray& getInterpLevels(bool expectPercFor1stCliq) const;

    /** Returns the sensitive strike */
    virtual void getSensitiveStrike(double spot,
                                    CDoubleArraySP sensitiveStrikes) const;

    /** Returns an array of CVolRequestLN, one for each cliquet */
    virtual CVolRequestLNArray getRequestsArray() const;

    /** Validates that this vol request is configured correctly for forward
        starting/spot starting. If fwdStartExpected is true and the request
        corresponds to an absolute strike, for example, an exception will
        be thrown. Similarly if fwdStartExpected is false and the request
        holds a % strike. */
    void validateFwdStart(bool fwdStartExpected) const;

    virtual void validatePop2Object();

protected:
    bool          firstCliquetInFuture;
    DateTimeArray liveCliqStartDates;
    DateTime      endDate;
    DoubleArray   interpLevels;

private:
    CliquetVolRequest();
    CliquetVolRequest(const CliquetVolRequest &rhs);
    CliquetVolRequest& operator=(const CliquetVolRequest& rhs);
};

typedef smartConstPtr<CliquetVolRequest> CliquetVolRequestConstSP;
typedef smartPtr<CliquetVolRequest> CliquetVolRequestSP;

DRLIB_END_NAMESPACE

#endif
