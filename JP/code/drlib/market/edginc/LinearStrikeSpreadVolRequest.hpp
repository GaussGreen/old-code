//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LinearStrikeSpreadVolRequest.hpp
//
//   Description : Linear Strike with spread vol request 
//
//   Author      : Andrew J Swain
//
//   Date        : 14 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef LINEARSTRIKESPREADVOLREQUEST_HPP
#define LINEARSTRIKESPREADVOLREQUEST_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VolRequestLNStrike.hpp"

DRLIB_BEGIN_NAMESPACE

/** This class represents a volatility interpolation algorithm which
    does the traditional straight forward vol interpolation except
    the fwd starting spread is supplied, not derived */
class MARKET_DLL LinearStrikeSpreadVolRequest: public VolRequestLNStrike {
public:
    static CClassConstSP const TYPE;
    friend class LinearStrikeSpreadVolRequestHelper;

    /** the interp level is floored at 0.0 (makes vega matrix nicer) */
    LinearStrikeSpreadVolRequest(double          interpLevel,
                                 const DateTime& startDate,
                                 const DateTime& endDate,
                                 double          spread);

    /** Returns the interpolation level. The expectPerc must be true
        for this vol request */
    double getInterpLevel(bool expectPerc) const;
    /** Returns the end date for vol interpolation */
    const DateTime& getEndDate() const;

    /** return the spread */
    double getSpread() const;

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
    double       spread;
private:
    LinearStrikeSpreadVolRequest();
    LinearStrikeSpreadVolRequest(const LinearStrikeSpreadVolRequest &rhs);
    LinearStrikeSpreadVolRequest& operator=(const LinearStrikeSpreadVolRequest& rhs);
};

typedef smartConstPtr<LinearStrikeSpreadVolRequest> LinearStrikeSpreadVolRequestConstSP;
typedef smartPtr<LinearStrikeSpreadVolRequest> LinearStrikeSpreadVolRequestSP;

DRLIB_END_NAMESPACE

#endif
