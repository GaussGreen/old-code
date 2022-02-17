//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenPortfolioLoss.cpp
//
//   Description : A Generator of MC Portfolio Loss State Variables
//
//   Author      : Lawrence Siu
//
//   Date        : 
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"

#if 0

#include "edginc/SVGenPortfolioLoss.hpp"
#include "edginc/SVQmcPortfolioLoss.hpp" // the implementation of SV

DRLIB_BEGIN_NAMESPACE

/** Constructor - from  an array of dates. 
    For computing default times between today and each date in dates */
SVGenPortfolioLoss::SVGenPortfolioLoss(
    const DateTime& _today,
    const vector<ICDSParSpreadsConstSP>& _cdsCurves,
    const DoubleArraySP _notionals,
    const DateTimeArraySP _portfolioLossDates)
    :
    today(_today),
    notionals(_notionals),
    portfolioLossDates(_portfolioLossDates),
    ddSvGens()
{
    static const string routine = "SVGenPortfolioLoss::SVGenPortfolioLoss";
    try {
        validate(_today, _cdsCurves, _notionals, _portfolioLossDates);
        createDDSvGens(_today, _cdsCurves, _portfolioLossDates);
    } 
    catch(exception& e){
        throw ModelException(e, routine);
    }
}

void SVGenPortfolioLoss::validate(    
    const DateTime& _today,
    const vector<ICDSParSpreadsConstSP>& _cdsCurves,
    const DoubleArraySP _notionals,
    const DateTimeArraySP _portfolioLossDates)
{
    const string & method = "SVGenPortfolioLoss::validate";
    if (_cdsCurves.size() != _notionals->size())
        throw ModelException(method, "Number of cds curves not equal the number of notionals.");
    if (_cdsCurves.empty())
        throw ModelException(method, "Empty cds curve vector.");
    if (_notionals->empty())
        throw ModelException(method, "Empty notional vector.");
    
    DateTime::ensureStrictlyIncreasing(*_portfolioLossDates, "Portfolio loss observation dates", true);

    DateTimeArray futureLossDates(_today.getFutureDates(*_portfolioLossDates));
    if (futureLossDates.empty())
        throw ModelException(method, "All portfolio loss dates are before or equal to today.");

    // check non-positive notionals?
}

// create an internal vector of date of default SV gens
void SVGenPortfolioLoss::createDDSvGens(
    const DateTime & _today, 
    const vector<ICDSParSpreadsConstSP>& _cdsCurves,
    const DateTimeArraySP _portfolioLossDates)
{
    // create date of default SV Gens:
    for (size_t i = 0; i < _cdsCurves.size(); ++i)
        ddSvGens.push_back(SVGenDateOfDefaultSP(new SVGenDateOfDefault(_today, _cdsCurves[i], _portfolioLossDates->back())));
}

/** Appends 'true' (ie non derived) state variable generators
required to the supplied collector. Implementations typically call
IStateVariableCollector::append */
void SVGenPortfolioLoss::collectStateVars(IStateVariableCollectorSP svCollector) const
{
    for (size_t i = 0; i < ddSvGens.size(); ++i)
        svCollector->append(ddSvGens[i].get());
}

/** Create the corresponding State Variable for this State
Variable Generator (from IStateVariableGen interface). The
previous IStateVariableSP (may be null) should be passed in.  The
return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenPortfolioLoss::create(
    IStateVariableSP oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const 
{ 
    return getPortfolioLossSV(pathGen); 
}

SVPortfolioLossSP SVGenPortfolioLoss::getPortfolioLossSV(
    IStateVariableGen::IStateGen* pathGen) const 
{
    // create a vector of date of default SVSPs locally
    vector<SVDateOfDefaultSP> ddSvs;
    for (size_t i = 0; i < ddSvGens.size(); ++i)
    {
        SVDateOfDefaultSP ddSv = ddSvGens[i]->getSVDateOfDefault(pathGen);
        ddSvs.push_back(ddSv);
    }

    SVPortfolioLossSP portfolioLossSV(new SVQmcPortfolioLossFullMC(ddSvs, notionals, portfolioLossDates));
    return portfolioLossSV; // portfolio loss SVSP
}


class SVGenPortfolioLoss::DeterminsticSV: public SVQmcPortfolioLossFullMC
{
public:
    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    
    /** All elements are inside array */
    virtual double element(int idx) const {
        return 0.0; //(*this)[idx]; // should throw an exception here?
    }

    DeterminsticSV(    
        const DateTime& _today,
        const vector<ICDSParSpreadsConstSP>& _cdsCurves,
        const DoubleArraySP _notionals,
        const DateTimeArraySP _portfolioLossDates,
        bool _doingThePast)
        :
        SVQmcPortfolioLossFullMC(vector<SVDateOfDefaultSP>(), DoubleArraySP(NULL), _portfolioLossDates),
        doingThePast(_doingThePast),
        today(today)
    {
        static const string routine = "SVGenPortfolioLoss::DeterminsticSV::DeterminsticSV";
        try {
            SVPath::initialize(doingThePast ? 0 : 0, 0); // TODO: Might fail here, fix it!

        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    bool doingThePast;
    DateTime today;
};

/** For use by Path Generators that want to use determinstic rates. */
SVPortfolioLoss* SVGenPortfolioLoss::determinsticSV(bool doingPast) const{
    return new DeterminsticSV(
        today, 
        vector<ICDSParSpreadsConstSP>(), 
        DoubleArraySP(NULL), 
        DateTimeArraySP(NULL), 
        doingPast);
}

/** implementing 'visitor' model */
/*
void SVGenDateOfDefault::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}
*/

DRLIB_END_NAMESPACE
#endif
