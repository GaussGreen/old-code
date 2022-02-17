//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedSpot.hpp
//
//   Description : A Generator of MC EXPECTED Asset Spot State Variables
//                 aka Forward Prices
//
//   Date        : 19 Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef SVGenExpectedSpot_HPP
#define SVGenExpectedSpot_HPP

#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE
/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVExpectedFX : public virtual ISVBase
{
public:
    virtual ~SVExpectedFX() {}
    virtual double getFwdFX(int idx) const = 0;
    virtual double getFwd(int idx) const = 0; // legacy name
};
DECLARE(SVExpectedFX);

class MCARLO_DLL SVExpectedEQ : public virtual ISVBase
{
public:
    virtual ~SVExpectedEQ() {}
    virtual double getFwdPrice(int idx) const = 0;
    virtual double getFwd(int idx) const = 0; // legacy name
};
DECLARE(SVExpectedEQ);

/** A Generator of MC Expected Spot State Variables.
    It's important to understand the difference between regular spot prices
    (as per SVGenSpot) and expected spot prices.
    Here we are saying, at a given date in the simulation, what is the
    expected spot at that date between two future dates.
    Note that here for simplicity we don't support an array of expected spots.
    Might want to review this at some point.
 */
class MCARLO_DLL SVGenExpectedSpot:    virtual public IElemStateVariableGen,
                                    public virtual VirtualDestructorBase
{
public:

    /** Interface for the state variable that SVGenExpectedSpot produces is
        either SVExpectedFXSP or SVExpectedEQSP and is in the file:
        IQMCStateVariables.hpp */

    /** Constructor - for computing on calcDate the expected spot values on
        dates */
    SVGenExpectedSpot(int                   assetIdx, // into MultiMarketFactors
                   const DateTime&       calcDate, // When to compute the
                   const DateTimeArray&  dates); // expected spot on these dates

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP              oldStateVar,
                                 IStateVariableGen::IStateGen*  pathGen) const;

    /** Returns a MC Expected Spot state variable which then
        provides access to the values etc. This is the method that
        products should call to get an "SVGenExpectedSpot::IStateVar". */
    //SVExpectedFXSP getExpSpotSV(IStateVariableGen::IStateGen* pathGen) const;
    SVPathSP getExpSpotSV(IStateVariableGen::IStateGen* pathGen) const;
    SVExpectedFXSP getExpFXSV(IStateVariableGen::IStateGen* pathGen) const;
    SVExpectedEQSP getExpEQSV(IStateVariableGen::IStateGen* pathGen) const;

    /** Returns the date on which the expected value should be computed */
    const DateTime&   getCalcDate() const;

    /** Returns the dates for when the expected spot should be computed */
    const DateTimeArray& getDates() const;

    /** Create an SVExpectedFX which has zero for all past values */
    SVExpectedFXSP createHistoricZeroSV(const DateTime& today) const;

    /** Returns the index of the associated asset (into MultiMarketFactors) */
    int getAssetIdx() const;

    /** Returns the number of dates */
    int numDates() const;

    /** If log of forward rate is required instead of the rate itself */
    bool logRequired() const { return false; }

    void attachSVGen(IElemStateVariableGenVisitor*) const;
private:
    class ZeroSV;
    ///// fields (excluding refCount) ////////
    int               assetIdx;  // identifies asset (in MultiMarketFactors)
    DateTime          calcDate; // when to compute expected value
    DateTimeArray     dates;   // dates to compute expected value to
};

typedef smartPtr<SVGenExpectedSpot> SVGenExpectedSpotSP;

DRLIB_END_NAMESPACE

#endif // SVGenExpectedSpot_HPP
