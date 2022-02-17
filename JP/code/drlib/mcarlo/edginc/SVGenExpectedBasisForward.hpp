//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenExpectedBasisFwdSpread.hpp
//
//   Description : A Generator of MC EXPECTED Basis Forward Spread State Variables
//
//   Author      : Anatoly V Morosov
//
//
//----------------------------------------------------------------------------

#ifndef SVGenExpectedBasisFwdSpread_HPP
#define SVGenExpectedBasisFwdSpread_HPP

#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE
/** A type-specific interface to the MC State Variable */
class MCARLO_DLL SVExpectedBasisFwdSpread : public virtual ISVBase
{
public:
    virtual ~SVExpectedBasisFwdSpread() {}
    virtual double getFwdSpread( int idx ) const = 0;
    virtual double firstSpread() const = 0;
};
DECLARE(SVExpectedBasisFwdSpread);

/** A Generator of MC Expected Basis Forward Spread State Variables.
    
 */
class MCARLO_DLL SVGenExpectedBasisFwdSpread: 
    virtual public IElemStateVariableGen,
    virtual public VirtualDestructorBase {
public:
    /** Constructor - for computing expected forward spread on calcDate between
        pvDate and each date in dates. Note that currently past values for this
        are not supported. In particular a value of 0 is returned when calcDate
        is in the past (even if computeLog is true) */
    SVGenExpectedBasisFwdSpread( IBasisIndexCurveConstSP  basisCurve,
                              const DateTime&        calcDate,   // when to compute
                              const DateTimeArray&   resetDates);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  
        The return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP                oldStateVar,
                                    IStateVariableGen::IStateGen*   pathGen) const;

    /** Returns a Expected Survival Discount Factor state variable which then
        provides access to the values etc. This is the method that
        products should call to get an SVExpSurvDiscFactor. */
    
    SVExpectedBasisFwdSpreadSP getSVExpectedBasisFwdSpread(
        IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
    IBasisIndexCurveConstSP getBasisCurve() const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getDates() const;
    
    /** Returns the date to which we pv to */
    const DateTime& getPVDate() const;
    
    /** Returns the date on which the expected value should be computed */
    const DateTime& getCalcDate() const;
    
    /** name of underlying credit */ 
    const string getName() const;
    
    /** For use by Path Generators (past or future) that want to use
        determinstic rates */
    SVExpectedBasisFwdSpread* determinsticSV(
        const DateTime& today,
        bool            doingPast) const;

    virtual void attachSVGen(IElemStateVariableGenVisitor*) const;    // support visitor pattern

private:
    class DeterminsticSV;
    friend class DeterminsticSV;
    DateTime                calcDate;           // when to compute expected value = pvDate
    DateTimeArray           dates;
    IBasisIndexCurveConstSP basisCurve;
};

DECLARE(SVGenExpectedBasisFwdSpread);

DRLIB_END_NAMESPACE

#endif
