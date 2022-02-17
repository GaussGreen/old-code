//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DayCountConventionFactory.cpp
//
//   Description : Factory class for building DayCountConventions
//
//   Author      : Andrew J Swain
//
//   Date        : 5 November 2000
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/Actual365.hpp"
#include "edginc/Actual365FJ.hpp"
#include "edginc/ActualActual.hpp"
#include "edginc/B30360.hpp"
#include "edginc/B30360F.hpp"
#include "edginc/B30E360.hpp"
#include "edginc/B30EP360.hpp"
#include "edginc/B30E360I.hpp"
#include "edginc/Business252.hpp"
#include "edginc/ActualActualAFB.hpp"
#include "edginc/ActualActualISMA.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE
DayCountConvention* DayCountConventionFactory::make(const string& dcc) {
    static string routine = "DayCountConventionFactory::make";

    // Convert dcc string to upper case
    //string ucDcc = dcc;
    //transform(ucDcc.begin(), ucDcc.end(), ucDcc.begin(), toupper);
    string ucDcc;
    for (string::const_iterator i = dcc.begin(); i != dcc.end(); ++i)
    {
        ucDcc.push_back(toupper(*i));
    }

    // Storage for Bus/252 data
    char holName[512];
    int denom;
    char stop[2];

    if      (ucDcc == "ACT/360" || ucDcc == "ACTUAL/360") {
        return new Actual360();
    }
    else if (ucDcc == "ACT/365F" || ucDcc == "ACTUAL/365F") {
        return new Actual365F();
    }
    else if (ucDcc == "ACT/365" || ucDcc == "ACTUAL/365") {
        return new Actual365();
    }
    else if (ucDcc == "ACT/365FJ" || ucDcc == "ACTUAL/365FJ") {
        return new Actual365FJ();
    }
    else if (ucDcc == "ACT/ACT" || ucDcc == "ACTUAL/ACTUAL") {
        return new ActualActual();
    }
    else if (ucDcc == "B30/360" || ucDcc == "30/360") {
        return new B30360();
    }
    else if (ucDcc == "B30/360F" || ucDcc == "30/360F") {
        return new B30360F();
    }
    else if (ucDcc == "B30E/360" || ucDcc == "30E/360") {
        return new B30E360();
    }
    else if (ucDcc == "B30E/360I" || ucDcc == "30E/360I") {
        return new B30E360I();
    }
    else if (ucDcc == "B30E+/360" || ucDcc == "30E+/360") {
        return new B30EP360();
    }
    else if (ucDcc == "B30EP/360" || ucDcc == "30EP/360") {
        return new B30EP360();
    }
    else if (ucDcc == "ACT/ACT AFB" || ucDcc == "ACTUAL/ACTUAL AFB") {
        return new ActualActualAFB();
    }
    else if (ucDcc == "ACT/ACT ISMA 1") {
        return ActualActualISMA::make(1);
    }
    else if (ucDcc == "ACT/ACT ISMA 2") {
        return ActualActualISMA::make(2);
    }    
    else if (ucDcc == "ACT/ACT ISMA 4") {
        return ActualActualISMA::make(4);
    }
    else if (ucDcc == "ACT/ACT ISMA 12") {
        return ActualActualISMA::make(12);
    }
    else if ((sscanf(ucDcc.c_str(), "BUS/%d%1s"     , &denom, stop) == 1 ||
              sscanf(ucDcc.c_str(), "BUSINESS/%d%1s", &denom, stop) == 1) &&
             denom == 252) {
        return new Business252();
    }
    else if ((sscanf(ucDcc.c_str(), "BUS/%d(%511[^)])%1s"     , &denom, holName, stop) == 2 ||
              sscanf(ucDcc.c_str(), "BUSINESS/%d(%511[^)])%1s", &denom, holName, stop) == 2) &&
             denom == 252) {

        // Rescan the holiday name to preserve the case
        sscanf(dcc.c_str(), "%*[^(](%511[^)]", holName);

        return new Business252(holName);
    }
    else {
        ModelException e(routine + ": Unknown day count " + dcc);
        e.addMsg(routine + ": Failed");
        throw e;
    }
}

DayCountConvention* DayCountConventionFactory::clone(const DayCountConvention* dcc) {
    return make(dcc->toString());
}

DayCountConventionFactory::DayCountConventionFactory() {
    // empty
}

class DayCountConventionAddin: public CObject{
    static CClassConstSP const TYPE;
    
/** The string representing the day count convention i.e. 30/360, Act/ACT, ... */
    string dcc;

/** Optional market data */
    MarketDataConstSP marketData;

    static IObjectSP createDCC(DayCountConventionAddin *params){
        static const string routine("DayCountConventionAddin::createDCC");
        try{
            DayCountConventionSP dcc(DayCountConventionFactory::make(params->dcc));

            if (params->marketData.get())
            {
                // May need to fetch holidays for Business/252 DCC
                params->marketData->fetchForNonMarketObjects(dcc, 
                                                             IModelConstSP(new NonPricingModel()),
                                                             "Business252");
            }                        

            return dcc;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    DayCountConventionAddin():  CObject(TYPE){}
    
    static void load(CClassSP& clazz){
        REGISTER(DayCountConventionAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDCC);
        FIELD(dcc, "Day Count Convention");
        FIELD(marketData, "market data");
        FIELD_MAKE_OPTIONAL(marketData);
        Addin::registerClassObjectMethod(
            "DAY_COUNT_CONVENTION",
            Addin::MARKET,
            "Create a Day Count Convention object from a string",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)createDCC);
    }
    
    static IObject* defaultDCC(){
        return new DayCountConventionAddin();
    }
};


CClassConstSP const DayCountConventionAddin::TYPE = CClass::registerClassLoadMethod(
    "DayCountConventionAddin", typeid(DayCountConventionAddin), load);


DRLIB_END_NAMESPACE

