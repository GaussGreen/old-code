//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenAggregatedSurvivalDiscFactor.hpp
//
//   Description : A Generator that captures dates for both MC Survival Discount 
//                 Factor and MC EXPECTED Survival Discount Factor State Variables
//
//   Author      : Eva X Strasser
//
//
//----------------------------------------------------------------------------

#ifndef SVGenAggregatedSurvivalDiscFactor_HPP
#define SVGenAggregatedSurvivalDiscFactor_HPP

#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/ISVBase.hpp"

DRLIB_BEGIN_NAMESPACE

class IQMCDiffusibleCreditSpreadBase;

/** A type-specific interface to the MC Discount Factor State Variable */
class MCARLO_DLL SVAggregatedSurvDiscFactor : public virtual ISVBase
{
public:
    virtual ~SVAggregatedSurvDiscFactor() {}
    virtual IQMCDiffusibleCreditSpreadBase* getAsset(void) const = 0;
};
DECLARE(SVAggregatedSurvDiscFactor);

class IElemStateVariableGenVisitor;
class QMCGenDiffusibleCredit;
/**
    In some cases we need to create SV for different assets share the same dates. (ex: SimpathI can have 10000 assets and we need DFs for the same diffusion and forward dates. In this case, it is a serious overkill to create DateTimeArrays again and again. Instead we can create an AggregatedSV, that will share the set of dates with similar SVs (i.e. there will be one AggregatedSV per asset). Moreover, since timeline is the same, the corresponding Idx arrays are the same! so we can share them as well.
    To accommodate the fact that different assets may have different maturity dates, we pass that information as well, so it is used when maxDiffusion/maxCurveDate are calculated.
    Note that currently we join all forward dates, even though it would not be much more difficult to keep them separately.
 */
class MCARLO_DLL SVGenAggregatedSurvivalDiscFactor: 
    virtual public IElemStateVariableGen,
    virtual public VirtualDestructorBase {
public:
    SVGenAggregatedSurvivalDiscFactor(
            DateTimeArraySP   _sdfDates,            // set of dates for discount factors
            SpotIdxArraySP    _sdfIdxSP,
            DateTimeArraySP   _esdfRequestedDates,  //  union of {t_i} for all expected factors
            FwdIdxArraySP     _esdfReqIdxSP,
            DateTimeArraySP   _esdfForwardDates,    // union of all {T_j} for all expected factors
            FwdIdxArraySP     _esdfForIdxSP,
            const DateTime&        _maxDiffDate,        // asset specific max diffusion
            const DateTime&        _maxCurveDate,       // asset specific max curve
            ICDSParSpreadsConstSP  _cdsParSpreadCurve,
            bool                   _computeLog);
    
    virtual IStateVariableSP create(IStateVariableSP oldStateVar,
                                    IStateGen*     stateGen) const;

    /** Returns a Expected Survival Discount Factor state variable which then
    provides access to the values etc. This is the method that
    products should call to get an SVExpSurvDiscFactor. */

    SVAggregatedSurvDiscFactorSP getIQSVGenAggregatedSurvivalDiscFactorSV(
            IStateVariableGen::IStateGen* pathGen) const;

    /** Retrieve the CDS par spread curve associated with this SVGenExpectedSurvivalDiscFactor */
    ICDSParSpreadsConstSP getCDSParSpreadCurve() const;

    /** Retrieve the dates for which a discount factor is required */
    const DateTimeArray& getSDFDates() const;
    /** Retrieve measurement the dates for which a expected discount factor is required */
    const DateTimeArray& getESDFRequestedDates() const;
    
    /** Retrieve forward dates for which a expected discount factor is required */
    const DateTimeArray& getESDFForwardDates() const;
   
    /** is the log of the expected discount factor wanted */
    bool logRequired() const;

    /** name of underlying credit */
    const string getName() const;

    virtual void attachSVGen(IElemStateVariableGenVisitor*) const;    // support visitor pattern
    DateTime getMaxMatDate() const; 
    DateTime getMaxCurveDate() const; 
private:
    void                    check(DateTimeArraySP dates, const DateTime& bound);
    DateTime                maxMaturityDate;    // asset specific maxMaturity date
    DateTime                maxCurveDate;       // asset specific maxCurveDate
    DateTimeArraySP         sdfDates;           // when to compute expected value
    DateTimeArraySP         esdfRequestedDates;             // when to discount to
    DateTimeArraySP         esdfForwardDates;
    ICDSParSpreadsConstSP   cdsParSpreadCurve;
    bool                    computeLog;
public: // add accessors later...    
    SpotIdxArraySP          sdfIdxSP;
    FwdIdxArraySP           esdfReqIdxSP;
    FwdIdxArraySP           esdfForIdxSP;
    friend class QMCGenDiffusibleCredit; // FIXME: write proper accessors
};

DECLARE(SVGenAggregatedSurvivalDiscFactor);

DRLIB_END_NAMESPACE

#endif
