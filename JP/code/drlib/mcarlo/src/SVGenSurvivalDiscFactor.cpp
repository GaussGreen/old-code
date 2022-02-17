//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenSurvivalDiscFactor.cpp
//
//   Description : A Generator of MC Survival Discount Factor State Variables
//
//   Author      : Eva X Strasser
//
//   Date        : 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"
//#include "edginc/MCPathConfigSRMGenSV.hpp"
DRLIB_BEGIN_NAMESPACE

/** Constructor - from  an array of dates. 
    For computing discount factors between today and each date in dates */
SVGenSurvivalDiscFactor::SVGenSurvivalDiscFactor(const DateTime&          today,
                                           ICDSParSpreadsConstSP    cdsParSpreadCurve,
                                           const DateTimeArray&     maturityDates):
today(today), cdsParSpreadCurve(cdsParSpreadCurve), maturityDates(maturityDates)
{
    static const string routine = "SVGenSurvivalDiscFactor::SVGenSurvivalDiscFactor";
    ASSERT(!maturityDates.empty());
    try {
        validate();
    } catch(exception& e){
        throw ModelException(e, routine);
    }
}

/** Constructor - from a single date */
SVGenSurvivalDiscFactor::SVGenSurvivalDiscFactor(const DateTime&        today,
                                           ICDSParSpreadsConstSP  cdsParSpreadCurve, 
                                           const DateTime&        maturityDate):
today(today), cdsParSpreadCurve(cdsParSpreadCurve), maturityDates(1, maturityDate){
    static const string routine = "SVGenSurvivalDiscFactor::SVGenSurvivalDiscFactor";
    try {
        validate();
    } catch(exception& e){
        throw ModelException(e, routine);
    }
}

/** Retrieve the CDS par spread curve associated with this 
    SVGenSurvivalDiscFactor */
ICDSParSpreadsConstSP SVGenSurvivalDiscFactor::getCDSParSpreadCurve() const {
    return cdsParSpreadCurve;
}

/** Retrieve the dates for which a survival discount factor is required */
const DateTimeArray& SVGenSurvivalDiscFactor::getDates() const {
    return maturityDates;
}

/** Create the corresponding State Variable for this State 
    Variable Generator (from IStateVariableGen interface). 
    The previous IStateVariableSP (may be null) should be passed in.  
    The return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenSurvivalDiscFactor::create(
    IStateVariableSP              oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const {
    return getSVSurvivalDiscFactor(pathGen);
}

SVSurvivalDiscFactorSP SVGenSurvivalDiscFactor::getSVSurvivalDiscFactor(
    IStateVariableGen::IStateGen* pathGen) const{
    SVSurvivalDiscFactorSP survivalDiscFactorSV(&dynamic_cast<SVSurvivalDiscFactor&>(*pathGen->create(this)));
    return survivalDiscFactorSV;
}

/** Basic validation. Should be called by all constructors after 
    population of all fields */
void SVGenSurvivalDiscFactor::validate() {
    static const string routine = "SVGenSurvivalDiscFactor::validate";
    try {
        DateTime::ensureIncreasing(maturityDates, "Maturity dates", 0);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


class SVGenSurvivalDiscFactor::DeterminsticSV: public SVQmcSurvivalDiscFactor {
public:
    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    
    /** All elements are inside array */
    virtual double element(int idx) const {
        return (*this)[idx];
    }

    DeterminsticSV(ICDSParSpreadsConstSP    cdsParSpreadCurve, 
                   const DateTime&          today,
                   const DateTimeArray&     maturityDates,
                   bool                     doingThePast):
    SVQmcSurvivalDiscFactor(NULL, maturityDates),
    doingThePast(doingThePast)
    {
        static const string routine = "SVGenSurvivalDiscFactor::DeterminsticSV::DeterminsticSV";
        try {
            // Figure out the start and end of loop
            int numMatDatesInPast = today.numPastDates(maturityDates);
            int end = doingThePast? numMatDatesInPast: maturityDates.size();
        
            // Only compute survival discount factors for future maturity dates
            survivalDiscountFactors = vector<double>(end, 0.0);
            DefaultRatesSP defaultRates = cdsParSpreadCurve->defaultRates(); 
            for (int i = numMatDatesInPast; i < end; i++){
                survivalDiscountFactors[i] = defaultRates->calcDefaultPV(today, maturityDates[i]);
            }
            
            // Create path and first discount factor
            if (survivalDiscountFactors.empty())
                 survivalDiscountFactors.push_back(0); // so firstDF returns 0
            
            SVPath::initialize(&survivalDiscountFactors[0],
                             doingThePast? 0 : numMatDatesInPast,
                             end);
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
private:
    vector<double> survivalDiscountFactors;
    bool           doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVSurvivalDiscFactor* SVGenSurvivalDiscFactor::determinsticSV(bool doingPast) const{
    return new DeterminsticSV(cdsParSpreadCurve, today, maturityDates, doingPast);
}

/** implementing 'visitor' model */
void SVGenSurvivalDiscFactor::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}


DRLIB_END_NAMESPACE
