//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolRequestLNStrike.hpp
//
//   Description : Abstract LN vol request for a strike-centric world
//
//   Author      : Andrew J Swain
//
//   Date        : 21 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef VOLREQUESTLNSTRIKE_HPP
#define VOLREQUESTLNSTRIKE_HPP

#include "edginc/VolRequestLN.hpp"

DRLIB_BEGIN_NAMESPACE

/** Abstract LN vol request for a strike-centric world */
class MARKET_DLL VolRequestLNStrike : public CVolRequestLN {
public:
    static CClassConstSP const TYPE;
    /** Scale the interpolation object */
    virtual void scale(double scaleFactor) = 0; // multiplicative scaling

    /** Returns the start date for forward starting volatility interpolation */
    virtual DateTime getStartDate() const = 0;

    /** Returns the sensitive strike */
    virtual void getSensitiveStrike(double spot, DoubleArraySP strikes) const;

    /** sets the strike of this request using the supplied strike. This method
        is somewhat dubious since different vol requests treat the strike
        differently.  */
    virtual void setStrike(double strike) = 0;

    /** Returns the strike as a percentage. If the strike is held internally
        as a level, the supplied spot is used to turn into a percentage */
    virtual double getPercStrike(double spot) const = 0;

    virtual ~VolRequestLNStrike();
protected:
    virtual void sensitiveStrikes(double               spot, 
                                  const DoubleArraySP& strikes) const = 0;

    VolRequestLNStrike(const CClassConstSP& clazz);
private:
    VolRequestLNStrike(const VolRequestLNStrike& rhs);
    VolRequestLNStrike& operator=(const VolRequestLNStrike& rhs);
};

// smart pointers for VolRequestLNStrike
typedef smartConstPtr<VolRequestLNStrike> VolRequestLNStrikeConstSP;
typedef smartPtr<VolRequestLNStrike> VolRequestLNStrikeSP;

// array of vol request
typedef array<VolRequestLNStrikeSP, VolRequestLNStrike> VolRequestLNStrikeArray;

// smart pointers for VolRequestLNStrikeArray
typedef smartConstPtr<VolRequestLNStrikeArray> VolRequestLNStrikeArrayConstSP;
typedef smartPtr<VolRequestLNStrikeArray> VolRequestLNStrikeArraySP;

DRLIB_END_NAMESPACE

#endif
