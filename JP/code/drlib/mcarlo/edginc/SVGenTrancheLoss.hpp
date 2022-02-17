//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenTrancheLoss.hpp
//
//   Description : A Generator of MC Tranche Loss State Variables
//
//   Author      : Lawrence Siu
//
//
//----------------------------------------------------------------------------

#ifndef SVGenTrancheLoss_HPP
#define SVGenTrancheLoss_HPP

#include "edginc/StateVariableClient.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/IConvolutor.hpp"

DRLIB_BEGIN_NAMESPACE


/** A type-specific interface to the MC Portfolio Loss State Variable */
class MCARLO_DLL SVTrancheLoss : public virtual ISVBase
{
public:
    virtual ~SVTrancheLoss() {}
    // get a vector of cumulative portfolio losses, with each loss corresponds to a measure date
    virtual CDoubleMatrixConstSP getTrancheLosses() = 0; 
};
DECLARE(SVTrancheLoss);


/** A Generator of MC Date Of Default State Variables. Can return simulated date of default. */

class MCARLO_DLL SVGenTrancheLoss  : public virtual IStateVariableGen,
                                     public virtual IStateVariableClient,
                                     public virtual VirtualDestructorBase
{
public:

    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector. Implementations typically call
    IStateVariableCollector::append */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const;

    /** Constructor. */
    SVGenTrancheLoss(
        const DateTime& _today,
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DoubleArraySP _notionals,
        const DateTimeArraySP _portfolioLossDates,
        const DoubleArraySP _lowerPct,
        const DoubleArraySP _upperPct,
		const IConvolutorSP _convolutor);      

    /** Create the corresponding State Variable for this State
    Variable Generator (from IStateVariableGen interface). The
    previous IStateVariableSP (may be null) should be passed in.  The
    return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP oldStateVar,
        IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a smart pointer to a portfolio loss SV */
    SVTrancheLossSP getTrancheLossSV(IStateVariableGen::IStateGen* pathGen) const;

    /** For use by Path Generators (past or future) that want to use deterministic rates */
    SVTrancheLoss* determinsticSV(bool doingPast) const;

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

    void createSDFSvGens(
        const DateTime & _today, 
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DateTimeArraySP _portfolioLossDates);

    class DeterminsticSV;
    DateTime today;
    DoubleArraySP notionals;
    DateTimeArraySP portfolioLossDates;
    DoubleArraySP lowerPct; // tranches
    DoubleArraySP upperPct; // tranches
	IConvolutorSP convolutor;
    vector<SVGenDateOfDefaultSP> ddSvGens; // vector of date of default SV gens
    vector<SVGenSurvivalDiscFactorSP> sdfSvGens; // vector of conditional SV gens

};

DECLARE(SVGenTrancheLoss);

DRLIB_END_NAMESPACE


#endif
