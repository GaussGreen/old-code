//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedEnergyFuture.hpp
//
//   Description : A Generator of MC EXPECTED Survival Discount Factor State Variables
//
//   Author      : Lawrence Siu
//
//
//----------------------------------------------------------------------------

#ifndef SVGenExpectedEnergyFuture_HPP
#define SVGenExpectedEnergyFuture_HPP

#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/EnergyFuturesCurve.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE
/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVExpEnergyFuture : public virtual ISVBase
{
public:
    virtual ~SVExpEnergyFuture() {}
    virtual const DateTimeArray & getFutureDates() const = 0;
    virtual double getExpFuturePrice( int idx ) const = 0;
};
DECLARE(SVExpEnergyFuture);

class MCARLO_DLL SVGenExpectedEnergyFuture:
    virtual public IElemStateVariableGen,
    public virtual VirtualDestructorBase
{
public:

    /** Note: The SVGenExpectedSurvivalDiscFactor::IStateVar has been obsoleted by the
        class IQSVGenExpectedSurvivalDiscFactorSV.   */

    /** Constructor - for computing expected survival discount values on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
    SVGenExpectedEnergyFuture(	const DateTime & _calcDate,				// when to compute
								EnergyFuturesCurveConstSP _futureCurve,
								const DateTimeArray & _dates,			// list of future maturities
								bool _computeLog						// true: do log
								);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.
        The return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP                oldStateVar,
                                    IStateVariableGen::IStateGen*   pathGen) const;

	/** for 'visitor' model of sorting SVGens,
	many thanks to Vladimir for suggesting it */
    virtual void attachSVGen(IElemStateVariableGenVisitor *) const;

    /** Returns a Expected Energy future state variable which then
        provides access to the values etc. This is the method that
        products should call to get an SVExpSurvDiscFactor. */
    SVExpEnergyFutureSP getSVExpEnergyFuture( IStateVariableGen::IStateGen* pathGen) const;

    /** For use by Path Generators (past or future) that want to use
        determinstic energy future prices */
    SVExpEnergyFuture* getPastSV(const DateTime& today, bool doingPast) const;

    /** Retrieve the Energy future curve associated with this SVGenExpectedEnergyFuture */
    EnergyFuturesCurveConstSP getEnergyFutureCurve() const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getDates() const;

//    /** Returns the date to which we pv to */
//    const DateTime& getPVDate() const;

    /** Returns the date on which the expected value should be computed */
    const DateTime& getCalcDate() const;

    /** is the log of the expected discount factor wanted */
    bool logRequired() const;

    /** name of underlying credit */
    const string getName() const;

private:
    class PastSV;
    friend class PastSV;

    DateTime					calcDate;           // when to compute expected value
    DateTimeArray				dates;
    EnergyFuturesCurveConstSP   futureCurve;
    bool						computeLog;
};

//typedef smartPtr<SVGenExpectedEnergyFuture> SVGenExpectedEnergyFutureSP;
DECLARE(SVGenExpectedEnergyFuture);

DRLIB_END_NAMESPACE

#endif
