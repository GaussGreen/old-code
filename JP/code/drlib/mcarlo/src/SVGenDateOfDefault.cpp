//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenDateOfDefault.cpp
//
//   Description : A Generator of MC Default Date State Variables
//
//   Author      : Lawrence Siu
//
//   Date        : 
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE

/** Constructor - from  an array of dates. 
    For computing default times between today and each date in dates */
SVGenDateOfDefault::SVGenDateOfDefault(
    const DateTime&         _today,
    ICDSParSpreadsConstSP   _cdsParSpreadCurve,
    const DateTime&         _maxDate)
    :
    today(_today), cdsParSpreadCurve(_cdsParSpreadCurve), maxDate(_maxDate) {
    static const string routine = "SVGenDateOfDefault::SVGenDateOfDefault";
    try {
        validate();
    } catch(exception& e){
        throw ModelException(e, routine);
    }
}

/** Create the corresponding State Variable for this State 
    Variable Generator (from IStateVariableGen interface). 
    The previous IStateVariableSP (may be null) should be passed in.  
    The return object may or may not be the same as oldStateVar. */
IStateVariableSP SVGenDateOfDefault::create(
    IStateVariableSP              oldStateVar,
    IStateVariableGen::IStateGen* pathGen) const 
{
    return getSVDateOfDefault(pathGen);
}

SVDateOfDefaultSP SVGenDateOfDefault::getSVDateOfDefault(
    IStateVariableGen::IStateGen* pathGen) const {
    SVDateOfDefaultSP ddSV(&dynamic_cast<SVDateOfDefault&>(*pathGen->create(this)));
    return ddSV; // default date SVSP
}


class SVGenDateOfDefault::DeterminsticSV: public SVQmcDateOfDefault {
public:
    /** Returns true if this state variable is being used to 'simulate'
        the past. This is a useful method for users of state variables - 
        need to see how hard it is to implement */
    virtual bool doingPast() const{
        return doingThePast;
    }
    
    /** All elements are inside array */
    virtual double element(int idx) const {
        return 0.0; //(*this)[idx];
    }

    /** asset specific name of one of the generic functions */
    inline DateTime getDateOfDefault() const { return maxDate.rollDate(1); }     // guaranteed to be alive

    DeterminsticSV(ICDSParSpreadsConstSP    cdsParSpreadCurve, 
                   const DateTime&          _today,
                   const DateTime&          _maxDate,
                   bool                     _doingThePast):
    SVQmcDateOfDefault(NULL),
    doingThePast(_doingThePast),
    today(_today),
    maxDate(_maxDate)
    {
        static const string routine = "SVGenDateOfDefault::DeterminsticSV::DeterminsticSV";
        try {
            /*
            // Figure out the start and end of loop
            int numMatDatesInPast = today.numPastDates(maturityDates);
            int end = doingThePast? numMatDatesInPast: maturityDates.size();
        
            // Only compute default dates for future maturity dates
            defaultDates = DateTimeArray(end, today); // already defaulted

            for (int i = numMatDatesInPast; i < end; i++){
                defaultDates[i] = maturityDates[i]; // no default
            
            // Create path and first default date
            if (defaultDates.empty())
                 defaultDates.push_back(today); // so firstDT returns 0
            
            Path::initialize(&defaultDates[0],
                             doingThePast? 0 : numMatDatesInPast,
                             end);
                             */
            SVPath::initialize(doingThePast ? 0 : 0, 0); // TODO: Might fail here, fix it!

        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }

private:
    //DateTimeArray   defaultDates;
    DateTime        today;
    DateTime        maxDate;
    bool            doingThePast;
};

/** For use by Path Generators that want to use determinstic rates. */
SVDateOfDefault* SVGenDateOfDefault::determinsticSV(bool doingPast) const{
    return new DeterminsticSV(cdsParSpreadCurve, today, maxDate, doingPast);
}

/** implementing 'visitor' model */
void SVGenDateOfDefault::attachSVGen(IElemStateVariableGenVisitor* sv) const
{
    sv->processSVGen(this);
}

DRLIB_END_NAMESPACE
