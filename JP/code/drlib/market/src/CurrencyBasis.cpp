//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CurrencyBasis.cpp
//
//   Description : Captures ccy basis curve and conventions
//
//   Author      : Andrew J Swain
//
//   Date        : 6 July 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_CURRENCYBASIS_CPP
#include "edginc/CurrencyBasis.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/B30360.hpp"

using namespace std;  // string

DRLIB_BEGIN_NAMESPACE

const string CurrencyBasis::SWAP_DAYS_ACTUAL = "ACT";
const string CurrencyBasis::SWAP_DAYS_BOND   = "30";

/** Optimized hashCode for performance : use for caching only */
int CurrencyBasis::hashCodeOpt() const {
    int hCode = (size_t) TYPE;
    hCode ^= expiries->hashCode();
    hCode ^= spreads.hashCode();
    return hCode;
}

// bread and butter validation
void CurrencyBasis::validatePop2Object() {
    static const string method("CurrencyBasis::validatePop2Object");
    try {
        if (expiries->size() != spreads.size()) {
            string msg = Format::toString(
                "expiries and spread arrays are of different lengths (%d vs %d)",
                expiries->size(), spreads.size());
            throw ModelException(method, msg);
        }
        if (expiries->empty()) {
            throw ModelException(method, "expiries and spreads are empty");
        }
        
        for (int i = 1; i < expiries->size(); ++i) {
            if ((*expiries)[i]->toDate(today) <= (*expiries)[i-1]->toDate(today)) {
                string msg = Format::toString(
                    "expiries are out of order - expiry %d [%s] is not before expiry %d [%s]",
                    i-1, (*expiries)[i-1]->toString().c_str(),
                    i, (*expiries)[i]->toString().c_str());
                throw ModelException(method, msg);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method, 
                             "Validation failure for " + name + " currency basis");
    }
}

// what's my name?
string CurrencyBasis::getName() const {
    return name;
}

// validate today
void CurrencyBasis::getMarket(const IModel* model, const MarketData* market) {
    static const string method = "CurrencyBasis::getMarket";

    try
    {
        // validate today
        market->GetReferenceDate(today);

        // get holidays
        hols.getData(model, market);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

// return basis spread to add to a given cash/swap rate and benchmark
double CurrencyBasis::basis(double        rate,
                            const Expiry* expiry,
                            const string& rateType,
                            int           mmrtDenom,
                            int           mmrtFreq,
                            const string& swapDays,
                            int           swapDenom,
                            int           swapFreq) const {
    static const string method("CurrencyBasis::basis");
    try {
        // check we understand the rate type
        if (rateType != YieldCurve::MMRT_RATE 
			&& rateType != YieldCurve::FUTURE_RATE
			&& rateType != YieldCurve::SWAP_RATE) {
            throw ModelException(method, "unsupported rate type: " + rateType);
        }

        // compute adjustment
        // several of these quantities could be transients as they are
        // immutable post-construction
        double adjust;

        if (rateType == YieldCurve::MMRT_RATE || rateType == YieldCurve::FUTURE_RATE) {
            adjust = 1.0;
        }
        else if (adjusted) {
            adjust = 1.0;
        }
        else {
            double floatAdj = 365.0/mmrtDenom/mmrtFreq;
            double fixAdj   = (swapDays==SWAP_DAYS_ACTUAL?365.0:360.0)/swapDenom;
            
            double adjDenom = pow(1+(rate*fixAdj)/swapFreq, (double)swapFreq/mmrtFreq) - 1.0;
            adjust   = rate*floatAdj/adjDenom;
        }

        // now get the basis
        DateTime target = expiry->toDate(today);
        double   spread = CashFlow::interpolate(*expiries.get(),
                                                spreads,
                                                today,
                                                target,
                                                true);  // extendFlat
                                
        return (spread * adjust);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// return basis spread to add to a given benchmark
double CurrencyBasis::cashBasis(const Expiry& expiry) const
{
    return basis(0.0, &expiry, YieldCurve::MMRT_RATE, 0, 0, "", 0, 0);
}

double CurrencyBasis::futuresBasis(const Expiry& expiry) const
{
    return basis(0.0, &expiry, YieldCurve::FUTURE_RATE, 0, 0, "", 0, 0);
}

double CurrencyBasis::swapBasis(double rate,
             const Expiry&             expiry,
             const DayCountConvention& mmrtDcc,
             int                       mmrtFreq,
             const DayCountConvention& swapDcc,
             int                       swapFreq) const
{
    static string method ="CurrencyBasis::swapBasis";

    int mmrtDenom;

    if (Actual360::TYPE->isInstance(mmrtDcc))
    {
        mmrtDenom = 360;
    }
    else if (Actual365::TYPE->isInstance(mmrtDcc))
    {
        mmrtDenom = 365;
    }
    else if (Actual365F::TYPE->isInstance(mmrtDcc))
    {
        mmrtDenom = 365;
    }
    else
    {
        string msg = Format::toString(
            "%s is an unsupported money market day count convention for calculating currency basis adjustments",
            mmrtDcc.toString().c_str());
        throw ModelException(method, msg);
    }

    string swapDays;
    if (Actual360::TYPE->isInstance(swapDcc))
    {
        swapDays = CurrencyBasis::SWAP_DAYS_ACTUAL;
    }
    else if (B30360::TYPE->isInstance(swapDcc))
    {
        swapDays = CurrencyBasis::SWAP_DAYS_BOND;
    }
    else if (Actual365::TYPE->isInstance(swapDcc))
    {
        swapDays = CurrencyBasis::SWAP_DAYS_ACTUAL;
    }
    else if (Actual365F::TYPE->isInstance(swapDcc))
    {
        swapDays = CurrencyBasis::SWAP_DAYS_ACTUAL;
    }
    else
    {
        string msg = Format::toString(
            "%s is an unsupported swap day count convention for calculating currency basis adjustments",
            swapDcc.toString().c_str());
        throw ModelException(method, msg);
    }

    return basis(
        rate, 
        &expiry, 
        YieldCurve::SWAP_RATE, 
        mmrtDenom, 
        mmrtFreq, 
        swapDays, 
        swapDcc.daysPerYear(), 
        swapFreq);
}


// is this the same as another CurrencyBasis
bool CurrencyBasis::equals(const CurrencyBasis* cb2) const {
    if (!cb2) {
        return false;
    }

    if (name != cb2->name) {
        return false;
    }
    if (today != cb2->today) {
         return false;
    }       

    if (expiries->size() != cb2->expiries->size() ||
        spreads.size()  != cb2->spreads.size()) {
        return false;
    }

    ExpiryArray exp1 = *(expiries.get());
    ExpiryArray exp2 = *(cb2->expiries.get());
    for (int i = 0; i < expiries->size(); i++) {
        if (!exp1[i]->equals(exp2[i].get())){
            return false;
        }        
        if (!Maths::equals(spreads[i], cb2->spreads[i])){
            return false;
        }
    }

    if (hols.get() && !hols->equals(cb2->hols.get()))
        return false;

    if (adjusted != cb2->adjusted) {
        return false;
    }

    return true;
}

const HolidayWrapper& CurrencyBasis::getHolidays() const {
	return hols;
}

/** for reflection */
CurrencyBasis::CurrencyBasis():MarketObject(TYPE), adjusted(false) {}



/*********************
 * Rho Parallel tweak
 *********************/

/** Returns name identifying currency basis for a Rho Parallel Tweak */
string CurrencyBasis::sensName (CurrencyBasisRhoParallelTweak* shift) const{
    return name;
}


/** Shifts the object using given shift */
bool CurrencyBasis::sensShift (CurrencyBasisRhoParallelTweak* shift){
    static const string method = "CurrencyBasis::sensShift";
	
    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            for (int i = 0; i < spreads.getLength(); ++i) {
                spreads[i] += shiftSize;
            }
        }
    } 
    catch (exception& e) {
       throw ModelException(e, method,
                            "CurrencyBasisRhoParallelTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}



/**********************
 * Rho Pointwise tweak
 **********************/

/** Returns name identifying currency basis for a Pointwise Tweak */
string CurrencyBasis::sensName (CurrencyBasisRhoPointwiseTweak* shift) const {
    return name;
}

/** Shifts the object using given TweakSizes */
bool CurrencyBasis::sensShift (CurrencyBasisRhoPointwiseTweak* shift) {
    static const string method = "CurrencyBasis::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            int i = shift->getExpiry()->search(expiries.get());
            spreads[i] += shiftSize;
        }
    }
    catch (exception& e) {
       throw ModelException(e, method,
                           "CurrencyBasisRhoPointwiseTweak failed for " + getName());
    }
    return false; // No sub-components need tweaking
}


ExpiryArrayConstSP CurrencyBasis::sensExpiries (CurrencyBasisRhoPointwiseTweak* shift) const {
    return expiries;
}


/*********************
 * Spread Level tweak
 *********************/

/** Returns name identifying currency basis for a Rho Parallel Tweak */
string CurrencyBasis::sensName (CurrencyBasisSpreadLevelTweak* shift) const{
    return name;
}


/** 
 * SETS the object using given shift - note that this is not a shift in the
 * traditional way but rather a scenario to set level 
 */
bool CurrencyBasis::sensShift (CurrencyBasisSpreadLevelTweak* shift){
    static const string method = "CurrencyBasis::sensShift";
	
    try {
        double shiftSize = shift->getShiftSize();

        for (int i = 0; i < spreads.getLength(); ++i) {
            spreads[i] = shiftSize;
        }
    } 
    catch (exception& e) {
       throw ModelException(e, method,
                            "CurrencyBasisSpreadLevelTweak failed for " + getName());
   }
   return false; // No sub-components need tweaking
}


/***************
 * Helper class
 ***************/
class CurrencyBasisHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CurrencyBasis, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(TweakableWith<CurrencyBasisRhoParallelTweak>);
        IMPLEMENTS(PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>);
        IMPLEMENTS(TweakableWith<CurrencyBasisSpreadLevelTweak>);
        EMPTY_SHELL_METHOD(defaultCurrencyBasis);
        FIELD(today, "today's date");
        FIELD_MAKE_OPTIONAL(today);
        FIELD(name, "identifier for market cache");
        FIELD(expiries, "expiries");
        FIELD(spreads, "ccy basis spreads");
        FIELD(hols, "basis holidays");
        FIELD(adjusted, "basis spreads input for swaps are already adjusted");
        FIELD_MAKE_OPTIONAL(adjusted);
        // by default don't get CurrencyBasis
        MarketDataFetcher::setDefaultRetrievalMode(CurrencyBasis::TYPE,
                                                   false, NULL);
    }
    

    static IObject* defaultCurrencyBasis(){
        return new CurrencyBasis();
    }
};



CClassConstSP const CurrencyBasis::TYPE = CClass::registerClassLoadMethod(
    "CurrencyBasis", typeid(CurrencyBasis), CurrencyBasisHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(CurrencyBasisWrapper);



/*************************
 * Ccy Basis Interpolator,
 * addin to help testing
 *************************/
class CcyBasisInterpolator: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    double          rate;
    ExpirySP        expiry;
    string          rateType;
    int             mmrtDenom;
    int             mmrtFreq;
    string          swapDays;
    int             swapDenom;
    int             swapFreq;
    CurrencyBasisSP ccyBasis;
    MarketDataSP    market;


    static double interp(CcyBasisInterpolator* params) {
        static const string routine = "CcyBasisInterpolator::interp";
        try {
            ClosedForm model;
            params->ccyBasis->getMarket(&model, params->market.get());

            return (params->ccyBasis->basis( params->rate,
                                             params->expiry.get(),
                                             params->rateType,
                                             params->mmrtDenom,
                                             params->mmrtFreq,
                                             params->swapDays,
                                             params->swapDenom,
                                             params->swapFreq));
        } 
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    CcyBasisInterpolator():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(CcyBasisInterpolator, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCcyBasisInterpolator);
        FIELD(rate, "rate");
        FIELD(expiry, "expiry");
        FIELD(rateType, "rate type");
        FIELD(mmrtDenom, "money market denominator");
        FIELD(mmrtFreq, "money market frequency");
        FIELD(swapDays, "swap day count days");
        FIELD(swapDenom, "swap day count denominator");
        FIELD(swapFreq, "swap frequency");
        FIELD(ccyBasis, "ccy basis curve");
        FIELD(market, "market");

       Addin::registerInstanceDoubleMethod("CCY_BASIS_SPREAD",
                                           Addin::RISK,
                                           "Compute ccy basis spread",
                                           TYPE,
                                           (Addin::DoubleMethod*)interp);
    }

    static IObject* defaultCcyBasisInterpolator(){
        return new CcyBasisInterpolator();
    }
};

   
CClassConstSP const CcyBasisInterpolator::TYPE = CClass::registerClassLoadMethod(
    "CcyBasisInterpolator", typeid(CcyBasisInterpolator), load);


DRLIB_END_NAMESPACE

