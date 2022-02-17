//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenDateOfDefault.hpp
//
//   Description : A Generator of MC Portfolio Loss State Variables
//
//   Author      : Lawrence Siu
//
//
//----------------------------------------------------------------------------

#ifndef SVGenPortfolioLoss_HPP
#define SVGenPortfolioLoss_HPP

#include "edginc/StateVariableClient.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE


/** A type-specific interface to the MC Portfolio Loss State Variable */
class MCARLO_DLL SVPortfolioLoss : public virtual ISVBase
{
public:
    virtual ~SVPortfolioLoss() {}
    // get a vector of cumulative portfolio losses, with each loss corresponds to a measure date
    virtual DoubleArrayConstSP getPortfolioLoss() = 0; 
};
DECLARE(SVPortfolioLoss);


/** A Generator of MC Date Of Default State Variables. Can return simulated date of default. */

class MCARLO_DLL SVGenPortfolioLoss: public virtual IStateVariableGen,
                                     public virtual IStateVariableClient,
                                     public virtual VirtualDestructorBase
{
public:

    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector. Implementations typically call
    IStateVariableCollector::append */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Constructor. */
    SVGenPortfolioLoss(
        const DateTime& _today,
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DoubleArraySP _notionals,
        const DateTimeArraySP _portfolioLossDates);      

    /** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP oldStateVar,
        IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a smart pointer to a portfolio loss SV */
    SVPortfolioLossSP getPortfolioLossSV(IStateVariableGen::IStateGen* pathGen) const;

    /** For use by Path Generators (past or future) that want to use deterministic rates */
    SVPortfolioLoss* determinsticSV(bool doingPast) const;

private:
    /** Basic validation. Should be called by all constructors after
    population of all fields */
    void validate(    
        const DateTime& _today,
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DoubleArraySP _notionals,
        const DateTimeArraySP _portfolioLossDates);

    // create an internal vector of date of default SV gens
    void createDDSvGens(
        const DateTime & _today, 
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DateTimeArraySP _portfolioLossDates);

    class DeterminsticSV;
    DateTime today;
    DoubleArraySP notionals;
    DateTimeArraySP portfolioLossDates;
    vector<SVGenDateOfDefaultSP> ddSvGens; // vector of date of default SV gens

};

DECLARE(SVGenPortfolioLoss);

//typedef refCountPtr<SVGenPortfolioLoss> SVGenPortfolioLossSP;
//typedef vector<SVGenPortfolioLossSP> SVGenPortfolioLossArray;
//typedef refCountPtr<SVGenPortfolioLossArray> SVGenPortfolioLossArraySP;


DRLIB_END_NAMESPACE


#endif
