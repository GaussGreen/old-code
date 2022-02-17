//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearStrikeTSVolReq.hpp
//
//   Description : Linear Strike vol request which handles term structure
//                 for forward starting options 
//
//   Author      : Stephen Hope
//
//   Date        : 23 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef LINEARSTRIKE_TS_VOLREQUEST_HPP
#define LINEARSTRIKE_TS_VOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequestLNStrike.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class represents a volatility interpolation algorithm which
    handles term structure of the spread for forward starting options */

class MARKET_DLL LinearStrikeTSVolRequest: public VolRequestLNStrike {
public:
    static CClassConstSP const TYPE;
    friend class LinearStrikeTSVolRequestHelper;

    /** the interp level is floored at 0.0 (make vega matrix nicer) */
    LinearStrikeTSVolRequest(double          interpLevel,
                             const DateTime& startDate,
                             const DateTime& endDate,
                             bool            isFwdStarting);

    /** Returns the interpolation level. The expectPerc must match up with
        the isFwdStarting flag */
    double getInterpLevel(bool expectPerc) const;

    /** Returns the end date for vol interpolation */
    const DateTime& getEndDate() const;

    /** Scale the interpolation object */
    virtual void scale(double scaleFactor); // multiplicative scaling

    /** Returns the start date for forward starting volatility 
        interpolation. */
    virtual DateTime getStartDate() const;

    /** sets the strike of this request using the supplied strike */
    virtual void setStrike(double strike);

    /** Returns the strike as a percentage. If the strike is held internally
        as a level, the supplied spot is used to turn into a percentage */
    virtual double getPercStrike(double spot) const;

    /** Validates that this vol request is configured correctly for forward
        starting/spot starting. If fwdStartExpected is true and the request
        corresponds to an absolute strike, for example, an exception will
        be thrown. Similarly if fwdStartExpected is false and the request
        holds a % strike. */
    void validateFwdStart(bool fwdStartExpected) const;
protected:
    virtual void sensitiveStrikes(double               spot, 
                                  const DoubleArraySP& strikes) const;

    double       interpLevel;
    DateTime     startDate;
    DateTime     endDate;
    bool         isFwdStarting;
private:
    LinearStrikeTSVolRequest();
    LinearStrikeTSVolRequest(const LinearStrikeTSVolRequest &rhs);
    LinearStrikeTSVolRequest& operator=(const LinearStrikeTSVolRequest& rhs);
};

typedef smartConstPtr<LinearStrikeTSVolRequest> LinearStrikeTSVolRequestConstSP;
typedef smartPtr<LinearStrikeTSVolRequest> LinearStrikeTSVolRequestSP;

DRLIB_END_NAMESPACE

#endif
