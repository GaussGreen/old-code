//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : BootstrappedYieldCurve.cpp
//
//   Description : A classical yield curve of cash, futures & swap rates
//
//   Author      : Richard Appleton
//
//   Date        : 25th November 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/B30360.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ThetaFwdRate.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/RiskyCDSCurve.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/VolProcessedBS.hpp"

#include "edginc/BadDayFollowing.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/SimpleCashFlowStream.hpp"
#include "edginc/StubPlacement.hpp"
#include "edginc/Equals.hpp"
#include <vector>
#include <algorithm>
#include ext_hash_map



DRLIB_BEGIN_NAMESPACE

const string BootstrappedYieldCurve::DEFAULT_FUTURES_MATURITY = "3M";

static const double MIN_RATE = -1.0;
static const double MAX_RATE = +10.0;

#define CS_CURVE_NAME "CREDIT_SPREADS"

/** Hash function (needed by hash_map template) */
struct HashUtil{
    size_t operator()(const BootstrappedYieldCurveConstSP c) const {
        return c->zeroCurveHash();
    }
};

/** Equal function (needed by hash_map template) */
struct EqualUtil {
    bool operator()(const BootstrappedYieldCurveConstSP c1, const BootstrappedYieldCurveConstSP c2) const {
        return c1->zeroCurveEquals(c2.get());
    }
};


/** The hash table used for caching (BootstrappedYieldCurve, Zero Curves) pairs */
typedef hash_map<BootstrappedYieldCurveConstSP, ZeroPairSP, HashUtil, EqualUtil> CashSwapCurveGlobalCacheHashtable;


/**
 * Global cache used for all BootstrappedYieldCurve objects.
 * Essentially a hash table with :
 * - Key  = the BootstrappedYieldCurve itself
 * - Data = the associated Zero Curves
 * A "tweaked" BootstrappedYieldCurve is treated as any other BootstrappedYieldCurve (i.e. we
 * don't have to know the tweak size)
 * */
class CashSwapCurveGlobalCache {
public:
    /**
     * Gets the Zero Curves associated to that BootstrappedYieldCurve :
     * - computes and stores the Zero Curves if not already in the hashtable
     * - just returns the Zero Curves if already computed
     * */
    static const ZeroPairConstSP get(const BootstrappedYieldCurve* yc) {

        ASSERT(yc == NULL || yc->getRefCount() > 0); // if yc is valid but yc->refCount is 0,
        BootstrappedYieldCurveConstSP key(yc); // then yc will be deleted on return when temp shared ptr is destroyed

        // Calls "find(key)" only once to improve performance
        CashSwapCurveGlobalCacheHashtable::iterator iter = hashtable.find(key);
        if (iter == hashtable.end()) {
            // Key not found in the hash table : need to compute the
            // Zero Curves and insert the (BootstrappedYieldCurve, Zero Curves) pair
            // in the hash table

            // First test if the cache is not full
            if (hashtable.size() >= MAX_SIZE) {
                // Cache is full, so we simply empty it !
                // (might do something more clever here, but not sure...)
                hashtable.clear();
            }

            const BootstrappedYieldCurveConstSP newKey(dynamic_cast<BootstrappedYieldCurve*>(yc->clone()));
            // Here is the time consuming operation !
            ZeroPairSP newData = yc->zeroCurve();

            // Insert the new (key, data) pair into hashtable
            hashtable[newKey] = newData;
            return newData;
        } else {
            return iter->second;
        }
    }
private:

    // The cache itself ! (i.e. a (BootstrappedYieldCurve, Zero Curves) hash table)
    static CashSwapCurveGlobalCacheHashtable hashtable;

    // Maximum number of entries in the cache
    static const unsigned int MAX_SIZE;
};

/** The cache itself ! (i.e. a (BootstrappedYieldCurve, Zero Curves) hash table) */
CashSwapCurveGlobalCacheHashtable CashSwapCurveGlobalCache::hashtable(MAX_SIZE);

/** Maximum number of entries in the cache */
const unsigned int CashSwapCurveGlobalCache::MAX_SIZE = 1000;


const ZeroCurve& BootstrappedYieldCurve::ZCAccessor::get(const BootstrappedYieldCurve& yc) const
{
    // "local" caching
    if (zPair.get() == 0)
    {
        // "global" caching
        zPair = CashSwapCurveGlobalCache::get(&yc);
    }

    return yc.useProjectionCurve ? *zPair.get()->growZC.get() : *zPair.get()->discZC.get();
}

/** Clears locally cached zero curve */
void BootstrappedYieldCurve::ZCAccessor::reset() const
{
    zPair.reset();
}

BootstrappedYieldCurve::ZCAccessor::ZCAccessor(): zPair() {}


/** Clears locally cached zero curve */
void BootstrappedYieldCurve::fieldsUpdated(const CFieldArray& fields)
{
    // this is called by the infrastructure when a component field has been
    // changed
    zc.reset();
}


BootstrappedYieldCurve::BootstrappedYieldCurve(CClassConstSP clazz)
: YieldCurve(clazz),
  spotOffset(0),
  floatIvl(NULL),
  allExpiries(NULL),
  nosort(false),
  useProjectionCurve(false),
  isIndexCurve_(false)
{
}


BootstrappedYieldCurve::BootstrappedYieldCurve(
    const string&                  pCcy,
    const string&                  pName,
    const DateTime&                pToday,
    int                            pSpotOffset,
    const HolidayWrapper&          pHolidays,
    const ZeroCurveBenchmarkArray& pBenchmarks,
    const IZeroCurveFactory&       pFactory,
    const DayCountConvention&      pMoneyMarketDcc,
    const MaturityPeriod*          pFuturesMaturity,
    const IRVolBaseWrapper&        pIrVol,
    const BadDayConvention*        pBadDayConvention,
    const DayCountConvention&      pFixedDcc,
    const MaturityPeriod&          pFixedIvl,
    const DayCountConvention*      pFloatDcc,
    const MaturityPeriod*          pFloatIvl,
    const ExpiryArray*             pFixDates,  // if not NULL float rate has fixings (floatRateFixed = true)
    const DoubleArray*             pFixRates,  // if passed must be same size as fixDates
    const DayCountConvention*      pBasisDcc,
    const MaturityPeriod*          pBasisIvl,
    const CurrencyBasisWrapper&    pCcyBasis,
    bool                           isIndexCurve)
    :
    YieldCurve(TYPE),
    ccy(pCcy),
    name(pName),
    today(pToday),
    spotOffset(pSpotOffset),
    hols(pHolidays),
    zcMethod(const_cast<IZeroCurveFactory*>(&pFactory)),
    benchmarks(pBenchmarks),
    moneyMarketDayCount(const_cast<DayCountConvention*>(&pMoneyMarketDcc)),
    irVol(pIrVol),
    badDayConvention(const_cast<BadDayConvention*>(pBadDayConvention)),
    fixedDcc(&pFixedDcc),
    fixedIvl(&pFixedIvl),
    floatDcc(pFloatDcc),
    floatIvl(pFloatIvl),
    fixDates(pFixDates),
    fixRates(pFixRates),
    basisDcc(pBasisDcc),
    basisIvl(pBasisIvl),
    ccyBasis(pCcyBasis),
    futMaturity(pFuturesMaturity ? const_cast<MaturityPeriod*>(pFuturesMaturity)
        : new MaturityPeriod(DEFAULT_FUTURES_MATURITY)),
    allExpiries(NULL),
    nosort(false),
    useProjectionCurve(false),
    isIndexCurve_(isIndexCurve)
{
    static const string method = "BootstrappedYieldCurve::BootstrappedYieldCurve";

    try
    {
        if (spotOffset == 0)
        {
            // doesn't really matter what the holidays are in this case
            valueDate = today;
        }
        else if (hols.get())
        {
            // if we have holidays (otherwise done in getMarket)
            // calculate transient fields
            valueDate = hols->addBusinessDays(today, spotOffset);
        }
        else
        {
            throw ModelException(method,
                "Cannot set base date as holidays must be specified if spotOffset > 0");
        }

        zc.reset();

        // (validation for ccy basis usage done by curve factory)
        zcMethod->validate(
            today,
            valueDate,
            benchmarks,
            moneyMarketDayCount.get(),
            fixedDcc.get(),
            floatDcc.get(),
            basisDcc.get(),
            fixedIvl.get(),
            floatIvl.get(),
            basisIvl.get(),
            badDayConvention.get(),
            fixDates.get() != NULL,           // float rate fixed
            fixDates.get(),
            fixRates.get(),
            IRVolBaseConstSP(irVol.get()),
            HolidayConstSP(hols.get()),
            ccyBasis,
            futMaturity.get());
    }
    catch (exception& e)
    {
        throw ModelException(e, method, e.what());
    }
}


// not very robust at the moment.Need to check expiries too.
IYieldCurveSP BootstrappedYieldCurve::makeRiskyCurve(
    const CreditSpreadCurve& spreadCurve,
    const DateTime*          maturityDate) const
{
    static const string method = "BootstrappedYieldCurve::makeRiskyCurve";

    try
    {
        int i;
        BootstrappedYieldCurveSP risky(cloneYieldCurve(true));

        // create adjusted credit spread curve which matches the maturities of the yield curve

        ExpiryArrayConstSP creditSpreadExpiries = spreadCurve.getExpiries();

        DateTime benchmarkDate;
        for (i=0; i<benchmarks.size(); i++)
        {
            benchmarkDate = benchmarks[i]->getBenchmarkDate(today);
            risky->benchmarks[i]->setRate(risky->benchmarks[i]->getRate() +
                spreadCurve.getCurrentSpread(today,benchmarkDate));
        }

        risky->name = name + spreadCurve.getName();
        risky->ccy = ccy;

        // make sure that we have a valid risky curve [pv() causes curve to be built]
        try
        {
            risky->pv(benchmarks.back()->getBenchmarkDate(today));
            return IYieldCurveSP(risky.release());
        }
        catch (exception&)
        {
        }

        // if we got here the risky curve is not valid.
        if (maturityDate != NULL)
        {
            // make proper adjustment to risky rates if this does not
            // affect bond prices
            for (i=benchmarks.size()-1; i>0; i--)
            {
                if (benchmarks[i]->isTurn())
                {
                    throw ModelException(method, "Cannot create risky curve when using turns");
                }

                // we should ignore futures
                int j = i - 1;
                if (benchmarks[i]->isFuture())
                    continue;
                else
                {
                    while (j > -1 && benchmarks[j]->isFuture())
                        j--;

                    if (j < 0)  // there are only futures below i in the array
                        break;
                }

                if (benchmarks[j]->getEnd()->toDate(today) >= (*maturityDate))
                {
                    for (int k=i; k < benchmarks.size(); k++)
                        risky->benchmarks[k]->setRate(risky->benchmarks[j]->getRate());

                    // NB. rate change will cause hash code to differ from cached curve
                    risky->zc.reset();

                    // make sure that we have a valid risky curve [pv() causes curve to be built]
                    try
                    {
                        risky->pv(benchmarks.back()->getEnd()->toDate(today));
                        return IYieldCurveSP(risky.release());
                    }
                    catch (exception&)
                    {
                    }
                }
                else
                    break;
            }
        }

        // i == 0
        throw ModelException(method, "Cannot generate a valid risky curve");
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

/** make a credit spread curve from default rates */
CreditSpreadCurveSP BootstrappedYieldCurve::makeCreditSpreadCurve(
	const string&        credName,
	const CashFlowArray& defaultRates,
	double               recovery) const
{
    static const string method = "BootstrappedYieldCurve::makeCreditSpreadCurve";

    try
    {
        // get 'risky' curve
        IYieldCurveSP riskyCurve(zcMethod->makeRiskyCurve(*this, defaultRates, recovery));

        // backout risky benchmarks from risky yield curve
        ZeroCurveBenchmarkArraySP riskyBenchmarks(dynamic_cast<ZeroCurveBenchmarkArray*>(benchmarks.clone()));
        adjust(*riskyBenchmarks, valueDate, riskyCurve.get());

        /* all OK so determine shifts to the yield curve */
        DoubleArray credRates(0);
        ExpiryArraySP credExpiries(new ExpiryArray(0));

        for (int idx = 0; idx < benchmarks.size(); idx++)
        {
            // original ignored futures, and adjust() rejects turns
            if (benchmarks[idx]->isCash() || benchmarks[idx]->isSwap())
            {
                credRates.push_back((*riskyBenchmarks)[idx]->getRate() - benchmarks[idx]->getRate());
                credExpiries->push_back(ExpirySP(dynamic_cast<Expiry*>(benchmarks[idx]->getEnd()->clone())));
            }
        }

        return CreditSpreadCurveSP(new CreditSpreadCurve(credName, credExpiries.get(), credRates, recovery));
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


IYieldCurveSP BootstrappedYieldCurve::createForwardCurve(
    const DateTime& forwardDate) const
{
    static const string method = "BootstrappedYieldCurve::createForwardCurve";

    try
    {
        BootstrappedYieldCurve* fwdCurve = cloneYieldCurve(true);
        adjust(fwdCurve->benchmarks, forwardDate, NULL);
        return IYieldCurveSP(fwdCurve);
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


void BootstrappedYieldCurve::adjust(
    ZeroCurveBenchmarkArray& adjusted,
    const DateTime&          date,
    const IYieldCurve*       yc) const
{
    static string method = "BootstrappedYieldCurve::adjust";
    bool isForForward = (yc == NULL);
    const char* what = isForForward ? "forward" : "credit spread";

    for (int i=0; i<benchmarks.size(); i++)
    {
        // calculate the forward rate
        if (benchmarks[i]->isCash())
        {
            MoneyMarketBenchmark& benchmark = dynamic_cast<MoneyMarketBenchmark&>(*benchmarks[i]);
            const DayCountConvention* dcc = benchmark.getDcc();

            double rate;
            DateTime expiry = benchmark.getEnd()->toDate(date);

            if (isForForward)
            {
                // calculate forward money market rate
                rate = fwd(date,
                           expiry,
                           dcc ? dcc : moneyMarketDayCount.get(),
                           CompoundBasis::SIMPLE);
            }
            else
            {
                // calculate money market rate from the risky zero curve passed to this func
                double discount = yc->pv(date, expiry);
                rate = RateConversion::discountToRate(discount,
                                                      date,
                                                      expiry,
                                                      dcc ? dcc : moneyMarketDayCount.get(),
                                                      CompoundBasis::SIMPLE);
            }

            adjusted[i]->setRate(rate);
        }
        else if (benchmarks[i]->isFuture())
        {
            string msg = Format::toString("Cannot create %s curve when using futures", what);
            throw ModelException(method, msg);
        }
        else if (benchmarks[i]->isTurn())
        {
            string msg = Format::toString("Cannot create %s curve when using turns", what);
            throw ModelException(method, msg);
        }
        else if (benchmarks[i]->isSwap())
        {
            SwapBenchmark& benchmark = dynamic_cast<SwapBenchmark&>(*benchmarks[i]);

            // validate is at par with no instrument specific settings
            if (!Maths::equals(1.0,benchmark.getPrice()))
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using turns", what);
                throw ModelException(method, msg);
            }

            if (benchmark.hasAdjustment())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using swap adjustments", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getFixedIvl())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap fixed leg interval", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getFloatIvl())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap floating leg interval", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getBasisIvl())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap basis leg interval", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getFixedDcc())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap fixed leg day count convention", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getFloatDcc())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap floating leg day count convention", what);
                throw ModelException(method, msg);
            }

            if (benchmark.getBasisDcc())
            {
                string msg = Format::toString("Cannot create %s curve "
                    "when using custom swap basis leg day count convention", what);
                throw ModelException(method, msg);
            }

            double rate;
            if (isForForward)
            {
                // calculate forward swap rate
                rate = couponRate(date,
                                  benchmark.getEnd()->toDate(date),
                                  *fixedIvl,
                                  true,
                                  fixedDcc.get());
            }
            else
            {
                // calculate swap rate from the risky curve passed to this func
                rate = yc->couponRate(date,
                                  benchmark.getEnd()->toDate(date),
                                  *fixedIvl,
                                  true,
                                  fixedDcc.get());
            }

            adjusted[i]->setRate(rate);
        }
    }
}


BootstrappedYieldCurve* BootstrappedYieldCurve::getScaledCurve( const double& scalingFactor) const
{
    static const string method = "BootstrappedYieldCurve::getScaledCurve";

    try {
        BootstrappedYieldCurve* scaledCurve = cloneYieldCurve(true);

        for (int i=0; i<benchmarks.size(); i++) {
            // calculate forward money market rate
            scaledCurve->benchmarks[i]->setRate( scaledCurve->benchmarks[i]->getRate() * scalingFactor);
        }

        return scaledCurve;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}



BootstrappedYieldCurve::~BootstrappedYieldCurve() {
    // empty
}

/** @return Yield curve's currency */
string BootstrappedYieldCurve::getCcy() const {
    return ccy;
}

/** @return Yield curve's name - used to identify YC_RHO sensitivities */
string BootstrappedYieldCurve::getName() const {
    return name;
}

/** Returns tradeDate + n business days where n = spot offset */
DateTime BootstrappedYieldCurve::settles(const DateTime& tradeDate) const{
    return hols->addBusinessDays(tradeDate, spotOffset);
}

/** @return Yield curve's spot date */
DateTime BootstrappedYieldCurve::getSpotDate() const {
    return valueDate;
}

/** @return Yield curve's value date */
DateTime BootstrappedYieldCurve::getToday() const {
    return today;
}

/** Useful accessor methods */
ExpiryArrayConstSP BootstrappedYieldCurve::getExpiries() const
{
    int numBmarks = benchmarks.size();
    ExpiryArraySP expiries = ExpiryArraySP(new ExpiryArray(numBmarks));

    for (int i=0; i<numBmarks; i++)
    {
        const Expiry* end = benchmarks[i]->getEnd();
        (*expiries)[i] = ExpirySP::attachToRef(const_cast<Expiry*>(end));
    }

    return expiries;
}

StringArrayConstSP BootstrappedYieldCurve::getInstruments() const
{
    int numBmarks = benchmarks.size();
    StringArraySP instruments = StringArraySP(new StringArray(numBmarks));

    for (int i=0; i<numBmarks; i++)
    {
        //turns and futures are also cash, so process those types first
        if (benchmarks[i]->isTurn())
        {
            (*instruments)[i] = YieldCurve::TURN_RATE;
        }
        else if (benchmarks[i]->isFuture())
        {
            (*instruments)[i] = YieldCurve::FUTURE_RATE;
        }
        else if (benchmarks[i]->isCash())
        {
            (*instruments)[i] = YieldCurve::MMRT_RATE;
        }
        else if (benchmarks[i]->isSwap())
        {
            (*instruments)[i] = YieldCurve::SWAP_RATE;
        }
        else
        {
            throw ModelException("BootstrappedYieldCurve::getInstruments",
                                 "Unknown instrument type");
        }
    }

    return instruments;
}

/** @return spot offset */
int BootstrappedYieldCurve::getSpotOffset() const
{
    return spotOffset;
}


/** @return bad day convention */
BadDayConventionConstSP BootstrappedYieldCurve::getBadDayConvention() const
{
    return badDayConvention;
}


/** @return holidays */
const HolidayWrapper& BootstrappedYieldCurve::getHolidays() const
{
    return hols;
}


/** @return curve factory */
IZeroCurveFactorySP BootstrappedYieldCurve::getCurveFactory() const
{
    return zcMethod;
}


/** @return curve benchmarks */
const ZeroCurveBenchmarkArray& BootstrappedYieldCurve::getBenchmarks() const
{
    return benchmarks;
}


/** @return money market day count convention */
DayCountConventionConstSP BootstrappedYieldCurve::getMoneyMarketDcc() const
{
    return moneyMarketDayCount;
}


/** @return fixed leg day count convention */
DayCountConventionConstSP BootstrappedYieldCurve::getFixedDcc() const
{
    return fixedDcc;
}


/** @return floating leg day count convention */
DayCountConventionConstSP BootstrappedYieldCurve::getFloatDcc() const
{
    return floatDcc.get() ? floatDcc : fixedDcc;
}

/** @return basis leg day count convention */
DayCountConventionConstSP BootstrappedYieldCurve::getBasisDcc() const
{
    return basisDcc.get() ? basisDcc : getFloatDcc();
}


/** @return fixed leg day period */
MaturityPeriodConstSP BootstrappedYieldCurve::getFixedIvl() const
{
    return fixedIvl;
}


/** @return float leg day period */
MaturityPeriodConstSP BootstrappedYieldCurve::getFloatIvl() const
{
    return floatIvl.get() ? floatIvl : fixedIvl;
}


/** @return basis leg day period */
MaturityPeriodConstSP BootstrappedYieldCurve::getBasisIvl() const
{
    return basisIvl.get() ? basisIvl : getFloatIvl();
}


/** @return fixing dates */
ExpiryArrayConstSP BootstrappedYieldCurve::getFixDates() const
{
    return fixDates;
}


/** @return IR vol */
const IRVolBaseWrapper& BootstrappedYieldCurve::getIrVol() const
{
    return irVol;
}


/** Optimized hashCode for performance : use for caching only */
int BootstrappedYieldCurve::zeroCurveHash() const {
    int hCode = (size_t) getClass();
    hCode ^= hash_string(ccy);
    hCode ^= today.hashCode();
    hCode ^= valueDate.hashCode();
    hCode ^= ZeroCurveBenchmark::hashCode(benchmarks);
    hCode ^= CBool::hashCode(!ccyBasis.isEmpty());

    if (!ccyBasis.isEmpty() && ccyBasis.get())
    {
        hCode ^= ccyBasis->hashCodeOpt();
    }

    if (moneyMarketDayCount.get())
    {
        hCode ^= moneyMarketDayCount->hashCode();
    }

    if (badDayConvention.get())
    {
        hCode ^= badDayConvention->hashCode();
    }

    if (zcMethod.get())
    {
        hCode ^= zcMethod->hashCode();
    }

    if (futMaturity.get())
    {
        hCode ^= futMaturity->hashCode();
    }

    if (fixedDcc.get()) hCode ^= fixedDcc->hashCode();
    if (floatDcc.get()) hCode ^= floatDcc->hashCode();
    if (basisDcc.get()) hCode ^= basisDcc->hashCode();
    if (fixedIvl.get()) hCode ^= fixedIvl->hashCode();
    if (floatIvl.get()) hCode ^= floatIvl->hashCode();
    if (basisIvl.get()) hCode ^= basisIvl->hashCode();
    if (fixDates.get()) hCode ^= fixDates->hashCode();
    if (fixRates.get()) hCode ^= fixRates->hashCode();

    return hCode;
}

/**
 * Optimized equalTo for performance : use for caching only
 * No need to compare IR vols for caching !
 * No need to compare useProjectionCurve because we store pairs in Cache
 */
bool BootstrappedYieldCurve::zeroCurveEquals(const IYieldCurve* curve) const
{
    try
    {
        if (this == curve)    // obvious first test
            return true;

        if (!curve || getClass() != curve->getClass())
            return false;

        const BootstrappedYieldCurve* yc2 = STATIC_CAST(BootstrappedYieldCurve, curve);

        if (!today.equals(yc2->today))
            return false;

        // avoid comparing holidays (and spot offset) by comparing value date
        if (!valueDate.equals(yc2->valueDate))
            return false;

        if (!equalsClass<BadDayConvention>(badDayConvention, yc2->badDayConvention))
            return false;

        // assume day count conventions are unique by class - should really
        // put equals method on them
        if (!equalsClass<DayCountConvention>(moneyMarketDayCount, yc2->moneyMarketDayCount))
            return false;

        if (!equalsClass<DayCountConvention>(fixedDcc, yc2->fixedDcc))
            return false;

        if (!equalsValue<MaturityPeriod>(fixedIvl, yc2->fixedIvl))
            return false;

        if (!equalsClass<DayCountConvention>(floatDcc, yc2->floatDcc))
            return false;

        if (!equalsValue<MaturityPeriod>(floatIvl, yc2->floatIvl))
            return false;

        if (!equalsClass<DayCountConvention>(basisDcc, yc2->basisDcc))
            return false;

        if (!equalsValue<MaturityPeriod>(basisIvl, yc2->basisIvl))
            return false;

        if (!ZeroCurveBenchmark::equals(benchmarks,yc2->benchmarks))
            return false;

        if (!equalsValue<IZeroCurveFactory>(zcMethod, yc2->zcMethod))
            return false;

        if (!ccyBasis.isEmpty() != !yc2->ccyBasis.isEmpty())
            return false;

        if (!ccyBasis.isEmpty())
        {
            const CurrencyBasis* cb = ccyBasis.get();
            if (cb)
            {
                return (cb->equals(yc2->ccyBasis.get()));
            }

            return false;
        }

        // this code also gets triggered on construction, so may not have data
        // from market cache yet
        if (badDayConvention.get())
        {
            if (hols.get() && !hols->equals(yc2->hols.get()))
            {
                return false;
            }
        }

        return true;
    }
    catch (exception& e)
    {
        throw ModelException(e, "BootstrappedYieldCurve::equalTo");
    }
}


/** Returns a key used to optimise repeated calculations of
    discount factors */
YieldCurve::IKey* BootstrappedYieldCurve::logOfDiscFactorKey() const
{
    return zc.get(*this).logOfDiscFactorKey();
}


/** Compute discount factor between two dates
 * @param lodate Lower date
 * @param hidate Upper date
 * @return Discount factor between lodate & hidate
 */
double BootstrappedYieldCurve::pv(const DateTime& lodate, const DateTime& hidate) const
{
    static const string method = "BootstrappedYieldCurve::pv";
    try
    {
        if (lodate.equals(hidate))
        {
            return 1.0;
        }
        else
        {
            const ZeroCurve& zcurve = zc.get(*this);
            return zcurve.pv(lodate, hidate);
        }
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

/** Compute discount factor between value date and a date
 * @param date To get discount factor for
 * @return Discount factor between value date & given date
 */
double BootstrappedYieldCurve::pv(const DateTime& date) const
{
    return pv(today, date);
}

/** Interpolate zero coupon rate at a date
 * @param date Interpolate zero coupon rate here
 * @return Zero coupon rate at given date
 */
double BootstrappedYieldCurve::zero(const DateTime& date) const
{
    const ZeroCurve& zcurve = zc.get(*this);
    return zcurve.zeroCouponRate(date);
}

/** Returns rateMaturity->toDate(rateStart) bad day adjusted
    either using bad day convention supplied (if non null)
    together with the yield curve holidays or the yield curve's
    bad day convention and holidays */
DateTime BootstrappedYieldCurve::rateMaturity(
    const DateTime&         rateStart,
    const Expiry*           rateMaturity,
    const BadDayConvention* bdc) const
{ // optional
    DateTime mat(rateMaturity->toDate(rateStart));
    if (bdc)
    {
        return bdc->adjust(mat, hols.get());
    }
    // as we don't adjust the dates when we bootstrap the zero curve we should
    // not adjust here either
    return mat;
}


/** Interpolate a forward rate between two dates
 * see CompoundBasis for basis values
 */
double BootstrappedYieldCurve::fwd(
    const DateTime&           lodate,
    const DateTime&           hidate,
    const DayCountConvention* dcc,
    int                       basis) const
{
    static const string method = "BootstrappedYieldCurve::fwd";
    try
    {
        const ZeroCurve& zcurve = zc.get(*this);
        return zcurve.fwd(lodate, hidate, dcc, basis);
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


/** Interpolate a forward rate between two dates
 * if CMS, uses swapFrequency for the compounding basis
 */
double BootstrappedYieldCurve::fwd(
    const DateTime&           refixDate,
    const Expiry*             rateMaturity,
    const BadDayConvention*   bdc, // optional
    const DayCountConvention* dcc,
    const bool                isCMS) const
{
    static const string method = "BootstrappedYieldCurve::fwd";
    try
    {
        DateTime lodate = bdc->adjust(refixDate, hols.get());
        DateTime hidate = bdc->adjust(rateMaturity->toDate(lodate), hols.get());

        if (isCMS)
        {
            return couponRate(lodate, hidate, *fixedIvl, false /* stub at end */, dcc);
        }
        else
        {
            return fwd(lodate, hidate, dcc, CompoundBasis::SIMPLE);
        }
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


double BootstrappedYieldCurve::fwd(const DateTime&    payDate,
                       const DateTime&           refixDate,
                       const Expiry*             rateMaturity,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const
{
    static const string method = "BootstrappedYieldCurve::fwd";

    try
    {
        if (!irVol)
        {
            throw ModelException(method, "Volatility must be supplied for adjusted forwards from " + name);
        }

        //for convenience
        Actual365F          act365F;
        DateTime            rateEnd  = bdc->adjust(rateMaturity->toDate(refixDate),hols.get());
        double              yF_delay = dcc->years(refixDate, payDate);

        CVolRequestSP fwdVolRequest(new SwapMaturityVolRequest(rateMaturity));
        CVolProcessedSP fwdVolCurve(irVol->getProcessedVol(fwdVolRequest.get(), 0));
        CVolProcessedBS& fwdVolBS = dynamic_cast<CVolProcessedBS&>(*fwdVolCurve);

        BenchmarkDate bMark(refixDate);
        CVolRequestSP delayVolRequest(new SwapMaturityVolRequest(&bMark));
        CVolProcessedSP delayVolCurve(irVol->getProcessedVol(delayVolRequest.get(), 0));
        CVolProcessedBS& delayVolBS = dynamic_cast<CVolProcessedBS&>(*delayVolCurve);

        //get the unadjusted forward
        double forwardRate = fwd(refixDate,rateMaturity,bdc,dcc,isCMS);

        //now calculate convexity and delay adjustments
        //taken from GtoSwapConvexityAdj2 and GtoSwapDelayAdj2
        //but here we apply both adjustments since they are always used in conjunction

        //get cxv/delay inputs, derived from this CashSwapCurve
        double fwdZeroRate = fwd(refixDate,rateEnd,&act365F,CompoundBasis::ANNUAL);
        double yrs2reset   = act365F.years(valueDate, refixDate);
        double maturity    = act365F.years(refixDate,rateEnd);
        int    frequency   = isCMS ? fixedIvl->approxAnnualFrequency() : 1;
        double fwdVol      = fwdVolBS.CalcVol(valueDate, refixDate);
        double delayRate   = fwd(refixDate,payDate,dcc,CompoundBasis::SIMPLE);
        double delayVol    = delayVolBS.CalcVol(valueDate, payDate);
        double delayModDur = yF_delay / (1 + (delayRate * yF_delay));
        double correlation = 1.0;

        //validate derived inputs
        if (forwardRate < 0.0)
        {
            throw ModelException(method, "Negative forward rate.");
        }
        if (fwdZeroRate < 0.0)
        {
            throw ModelException(method, "Negative forward zero rate.");
        }
        if (yrs2reset < 0.0)
        {
            throw ModelException(method, "Reset in the past.");
        }
        if (maturity < 0.0)
        {
            throw ModelException(method, "Maturity in the past.");
        }
        if (fwdVol < 0.0)
        {
            throw ModelException(method, "Negative forward volatility.");
        }
        if (delayRate < 0.0)
        {
            throw ModelException(method, "Invalid delay rate.");
        }
        if (delayVol  < 0.0)
        {
            throw ModelException(method, "Negative delay volatility.");
        }
        if (delayModDur < 0.0)
        {
            throw ModelException(method, "Invalid modified duration.");
        }
        if (fabs(correlation) > 1.0)
        {
            throw ModelException(method, "Invalid correlation.");
        }

        //common variables
        double paymentsNo;  /* N : total no. of coupon pmnts for the bond */
        double fwDiscount;  /* Z{^[N]} : yield corresponding to forwardDiscount */
        double zcBondPrice; /* Z{[N]}(t;Y{^}) : price of zero coupon bond maturing in <maturity> years */
        double denom;       /* 1 - fwDiscount */
        double myRatio;     /* numerator/denom */
        double cFactor;     /* common factor used by both convexity & delay*/

        paymentsNo  = maturity * frequency;
        fwDiscount  = 1.0/pow((1.0 + fwdZeroRate),maturity);
        zcBondPrice = 1.0 / pow((1.0+forwardRate/frequency), paymentsNo);
        denom       = 1.0 - fwDiscount;              /* 1-Z{^,[N]}(t) */
        myRatio     = (1.0 - zcBondPrice)/denom;
        cFactor     = forwardRate * fwdVol * yrs2reset * myRatio;

        //convexity
        /* Convert fwdZCYield from annual compounding to compounding=frequency.
         */
        double fwdCompRate = RateConversion::discountToRate(fwDiscount,
                                                     refixDate,
                                                     rateEnd,
                                                     &act365F, //dcc,
                                                     frequency);

        double f17SecondFactor = 1.0 - paymentsNo * fwdCompRate * fwDiscount /
                                 (frequency + fwdCompRate) / (1 - zcBondPrice);
        double convexAdj = cFactor * fwdVol * myRatio * f17SecondFactor;

        //delay
        double delayAdj = - cFactor * correlation * delayRate *
                           delayModDur * delayVol;

        //now apply adjustments
        double rate = forwardRate + convexAdj + delayAdj;
        return rate;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


/** risky discount factor - the bond spreads are implicit in the yield, so this
    will do the same as pv */
double BootstrappedYieldCurve::riskyPV(const DateTime& lodate,
                              const DateTime& hidate,
                              double          cashFlow,
                              double          recoveryNotional) const {
    return pv(lodate, hidate) * cashFlow;
}


double BootstrappedYieldCurve::firstFixingRate() const
{
    if (fixDates.get() != NULL)
    {
        DateTime fixDate = DateFwdThenAdjust(
            getSpotDate(), *floatIvl, 1, *badDayConvention, *hols.get());

        for (int i = 0 ; i < fixDates->size() ; i++)
        {
            DateTime date = (*fixDates)[i]->toDate(getSpotDate());
            date = badDayConvention->adjust(date, hols.get());

            if (date.equals(fixDate))
            {
                return (*fixRates)[i];
            }
        }
    }

    return 0.0;
}


double BootstrappedYieldCurve::parSwapRate(const DateTime& maturity) const
{
    ZeroPairConstSP zp = CashSwapCurveGlobalCache::get(this);
    const ZeroCurve* discounting = zp->get(false);
    const ZeroCurve* estimating = zp->get(true);

    bool valueFloating = (discounting != estimating)
        || zcMethod->isFloatLegValued();

    return SwapTool::swapRate(
        *discounting,
        getSpotDate(),
        maturity,
        *fixedIvl,
        *fixedDcc,
        valueFloating,
        1.0,                                                // ... ie. at par
        estimating,
        floatIvl.get(),
        floatDcc.get(),
        (fixDates.get() != NULL && fixDates->size() > 0),   // fixed?
        firstFixingRate(),
        false,                                              // convexityDelayAdj
        irVol.get(),
        StubPlacement("AUTO"),
        *badDayConvention,
        *badDayConvention,
        *badDayConvention,
        *hols.get());
}


/** grab the dates used in the zero curve */
DateTimeArray BootstrappedYieldCurve::zeroDates() const
{
    return zc.get(*this).getDates();
}


/** drive which style of zero curve is used. Default is discounting */
void BootstrappedYieldCurve::setProjectionCurve(bool pUseProjectionCurve) const
{
    // projection curve only differs from discounting when using ccy basis
    useProjectionCurve = pUseProjectionCurve;
}


/** overrides CObject version to allow for easy default */
bool BootstrappedYieldCurve::accept(ICollector* collector) const
{
    if (!CClass::invokeAcceptMethod(this, collector))
    {
        // if no method registered try  vol
        return irVol->accept(collector);
    }

    return false;
}


/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
IVolProcessed* BootstrappedYieldCurve::getProcessedVol(
    const CVolRequest* volRequest) const
{
    if (irVol.isEmpty())
    {
        throw ModelException("BootstrappedYieldCurve::getProcessedVol",
                             "Pricing model requires IR Vol be supplied for "
                             "BootstrappedYieldCurve with name "+getName());
    }
    return irVol->getProcessedVol(volRequest, this);
}


// build a zero curve
ZeroPairSP BootstrappedYieldCurve::zeroCurve() const
{
    static const string method = "BootstrappedYieldCurve::zeroCurve";

    try
    {
        // curve type built is determined by factory type supplied
        return zcMethod->bootstrap(
            today,
            valueDate,
            benchmarks,
            moneyMarketDayCount.get(),
            fixedDcc.get(),
            floatDcc.get(),
            basisDcc.get(),
            fixedIvl.get(),
            floatIvl.get(),
            basisIvl.get(),
            badDayConvention.get(),
            fixDates.get() != NULL,     // float rate fixed?
            fixDates.get(),
            fixRates.get(),
            IRVolBaseConstSP(irVol.get()),
            HolidayConstSP(hols.get()),
            ccyBasis,
            futMaturity.get());
    }
    catch (exception &e)
    {
        string msg = Format::toString("zero curve failed for %s (%s)",
            getCcy().c_str(), getName().c_str());
        throw ModelException(e, method, msg);
    }
}


IPublicObject* BootstrappedYieldCurve::toPublicObject() const
{
    static const string method = "BootstrappedYieldCurve::toPublicObject";

    try
    {
        BootstrappedCurveAddin* iface = NULL;

        iface = new BootstrappedCurveAddin(
            ccy,
            name,
            today,
            spotOffset,
            hols,
            benchmarks,
            *zcMethod,
            moneyMarketDayCount->toString(),
            futMaturity.get(),
            irVol,
            badDayConvention.get() ? badDayConvention->toString() : "",
            fixedDcc->toString(),
            fixedIvl.get(),
            floatDcc.get() ? floatDcc->toString() : "",
            floatIvl.get(),
            fixDates.get(),
            fixRates.get(),
            basisDcc.get() ? basisDcc->toString() : "",
            basisIvl.get(),
            ccyBasis,
            isIndexCurve_);

        return iface;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

/** Override clone method to copy our extra data over */
IObject* BootstrappedYieldCurve::clone() const
{
    return cloneYieldCurve(false);
}

/** Same as usual clone() but optionally doesn't copy the local ZC cache */
BootstrappedYieldCurve* BootstrappedYieldCurve::cloneYieldCurve(bool withoutCache) const
{
    // first clone all the registered fields
    IObject* copy = CObject::clone();
    BootstrappedYieldCurve* yc = dynamic_cast<BootstrappedYieldCurve*>(copy);
    if (!yc)
    {
        delete copy;
        throw ModelException("BootstrappedYieldCurve::cloneYieldCurve"); // shouldn't happen
    }

    yc->zc = zc; // 'structure type' copy
    yc->valueDate = valueDate;

    if (withoutCache)
        yc->zc.reset();

    return yc;
}


/** Returns name identifying yield curve for IR delta pointwise */
string BootstrappedYieldCurve::sensName(const IRRatePointwise* shift) const
{
    return name;
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryWindowArrayConstSP
BootstrappedYieldCurve::sensQualifiers(const IRRatePointwise* shift) const
{
    return sensQualifiers((RatePointwise*) NULL);
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome BootstrappedYieldCurve::sensShift(const PropertyTweak<IRRatePointwise>& shift)
{
    static const string method = "BootstrappedYieldCurve::sensShift";

    try
    {
        checkCache();
        double shiftSize = shift.coefficient;
        if (!Maths::isZero(shiftSize))
        {
            int i = ZeroCurveBenchmark::search(benchmarks, *shift.qualifier->expiry);

            benchmarks[i]->setRate(benchmarks[i]->getRate() + shiftSize);

            // NB. rate change will cause hash code to differ from cached curve
            zc.reset();
        }
    }
    catch (exception& e)
    {
        throw ModelException(e,method,"Failed for "+getCcy()+ " ("+getName()+")");
    }

    // none of our components has a rho type sensitivity
    return TweakOutcome(shift.coefficient, false);
}

void BootstrappedYieldCurve::sensRestore(const PropertyTweak<IRRatePointwise>& shift)
{
    checkCache();
    double shiftSize = shift.coefficient;
    if (!Maths::isZero(shiftSize))
    {
        int i = ZeroCurveBenchmark::search(benchmarks, *shift.qualifier->expiry);

        benchmarks[i]->setRate(benchmarks[i]->getRate() - shiftSize);

        // NB. rate change will cause hash code to differ from cached curve
        zc.reset();
    }
}


/** Returns name identifying yield curve for rho parallel */
string BootstrappedYieldCurve::sensName(RateShift* shift) const
{
    return getCcy();
}

/** Shifts the object using given shift */
bool BootstrappedYieldCurve::sensShift(RateShift* shift)
{
    for (int i = 0; i < benchmarks.getLength(); i++)
    {
        double shiftSize = shift->shiftSize(today,benchmarks[i]->getBenchmarkDate(today));
        benchmarks[i]->setRate(Maths::max(benchmarks[i]->getRate() + shiftSize, 0.0)); // floor to zero
    }

    // NB. rate change will cause hash code to differ from cached curve
    zc.reset();

    return false; // none of our components has a rho type sensitivity
}


/** Returns name identifying yield curve for additive or multiplicative
 * weighted shift */
string BootstrappedYieldCurve::sensName(YCWeightedShift* shift) const
{
    return getCcy();
}

/** Shifts the object using given shift */
bool BootstrappedYieldCurve::sensShift(YCWeightedShift* shift)
{
    DoubleArray rates(benchmarks.size());
    checkCache();

    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        rates[i] = benchmarks[i]->getRate();
    }

    shift->shiftArray(allExpiries, &rates, today);

    for (int j = 0 ; j < benchmarks.size() ; j++)
    {
        benchmarks[j]->setRate(rates[j]);
    }

    // NB. rate change will cause hash code to differ from cached curve
    zc.reset();

    return false; // do not continue tweaking
}


/** Returns name identifying yield curve for rho parallel */
string BootstrappedYieldCurve::sensName(const RateParallel* shift) const
{
    return getCcy();
}

/** Shifts the object using given shift */
TweakOutcome BootstrappedYieldCurve::sensShift(const PropertyTweak<RateParallel>& shift)
{
    static const string method = "BootstrappedYieldCurve::sensShift";

    try
    {
        double shiftSize = shift.coefficient;
        if (!Maths::isZero(shiftSize))
        {
            for (int i = 0; i < benchmarks.getLength(); i++)
            {
                benchmarks[i]->setRate(benchmarks[i]->getRate() + shiftSize);
            }

            // NB. rate change will cause hash code to differ from cached curve
            zc.reset();
        }

        // none of our components has a vega type sensitivity
        return TweakOutcome(shift.coefficient, false);
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

/** Restores the object to its original form */
void BootstrappedYieldCurve::sensRestore(const PropertyTweak<RateParallel>& shift)
{
    double shiftSize = shift.coefficient;
    if (!Maths::isZero(shiftSize))
    {
        for (int i = 0; i < benchmarks.getLength(); i++)
        {
            benchmarks[i]->setRate(benchmarks[i]->getRate() - shiftSize);
        }

        // NB. rate change will cause hash code to differ from cached curve
        zc.reset();
    }
}

/** Returns the name of the yield curve - used to determine whether
    to tweak the object */
string BootstrappedYieldCurve::sensName(const RatePointwise* shift) const
{
    return getCcy();
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryWindowArrayConstSP BootstrappedYieldCurve::sensQualifiers(const RatePointwise* shift) const
{
    if (!allExpiries.get())
    { // this should never happen
        throw ModelException("BootstrappedYieldCurve::sensExpiries",
                             "expiry cache has not been set");
    }

    return ExpiryWindow::series(allExpiries);
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome BootstrappedYieldCurve::sensShift(const PropertyTweak<RatePointwise>& shift)
{
    static const string method = "BootstrappedYieldCurve::sensShift";

    try
    {
        checkCache();
        double shiftSize = shift.coefficient;
        if (!Maths::isZero(shiftSize))
        {
            // minor hackette - rho works at a ccy (not curve) level
            // so there may be many (say) EUR curves with different names
            // and potentially different expiries. It's pot luck which one
            // gets picked to provide the expiry list, so avoid failing
            // it we're tweaking a curve with a subset of the expiries
            int i;
            try
            {
                i = ZeroCurveBenchmark::search(benchmarks, *shift.qualifier->expiry);
            }
            catch (ModelException)
            {
                // exit gracefully if DR code fails (presumably on looking
                // for the expiry), but dump out on anything more sinister
                return TweakOutcome(shift.coefficient, false);
            }

            benchmarks[i]->setRate(benchmarks[i]->getRate() + shiftSize);

            // NB. rate change will cause hash code to differ from cached curve
            zc.reset();
        }
    }
    catch (exception& e)
    {
        throw ModelException(e,method,"Failed for "+getCcy()+
                             " ("+getName()+")");
    }

    // none of our components has a rho type sensitivity
    return TweakOutcome(shift.coefficient, false);
}

/** Restores the object to its original form */
void BootstrappedYieldCurve::sensRestore(const PropertyTweak<RatePointwise>& shift)
{
    checkCache();
    double shiftSize = shift.coefficient;
    if (!Maths::isZero(shiftSize))
    {
        // minor hackette - rho works at a ccy (not curve) level
        // so there may be many (say) EUR curves with different names
        // and potentially different expiries. It's pot luck which one
        // gets picked to provide the expiry list, so avoid failing
        // it we're tweaking a curve with a subset of the expiries
        int i;
        try
        {
            i = ZeroCurveBenchmark::search(benchmarks, *shift.qualifier->expiry);
        }
        catch (ModelException)
        {
            // exit gracefully if DR code fails (presumably on looking
            // for the expiry), but dump out on anything more sinister
            return;
        }

        benchmarks[i]->setRate(benchmarks[i]->getRate() - shiftSize);

        // NB. rate change will cause hash code to differ from cached curve
        zc.reset();
    }
}

bool BootstrappedYieldCurve::sensShift(Theta* shift)
{
    try
    {
        DateTime rolledDate = shift->rollDate(today);

        if (!(rolledDate.equals(today)))
        {
            today = rolledDate;

            bool isThetaFwdRate = !!(dynamic_cast<ThetaFwdRate*>(shift));

            if (!isThetaFwdRate)
            {
                // update value date
                valueDate = hols->addBusinessDays(today, spotOffset);
                expiryCache();
                zc.reset();
            }
        }
    }
    catch (exception& e){
        throw ModelException(e, "BootstrappedYieldCurve::sensShift (theta)");
    }
    return true; // need to shift the IR vol possibly
}

void BootstrappedYieldCurve::checkCache()
{
    if (!allExpiries.get())
    {
        expiryCache();
    }
}



// if we have futures in the curve, we need to rank the expiries in
// date order (for last sens date to work) and then flag which ones
// are futures (for the tweaking to work)
void BootstrappedYieldCurve::expiryCache()
{
    static const string method = "BootstrappedYieldCurve::expiryCache";

    try
    {
        if (valueDate.empty())
        {
            throw ModelException(method, "Curve value date has not been set for " + getName());
        }

        /*
         * For consistency with the ALIB ZC3 approach futures are sorted
         * according to value date, not 'today'.
         *
         * Also original QLib curve implementation did not sort inputs unless
         * some futures were used : for performance reasons this is retained
         * if using the original curve building interface (CashSwapCurve).
         */
        if (nosort)
        {
            for (int j = 0 ; j < benchmarks.size() ; j++)
            {
                if (benchmarks[j]->isFuture())
                {
                    ZeroCurveBenchmark::sort(benchmarks, valueDate);
                    break;
                }
            }
        }
        else
        {
            ZeroCurveBenchmark::sort(benchmarks, valueDate);
        }

        allExpiries = ExpiryArraySP(new ExpiryArray(benchmarks.size()));
        for (int i = 0 ; i < benchmarks.size() ; i++)
        {
            const Expiry* expiry = benchmarks[i]->isFuture() ? benchmarks[i]->getStart() : benchmarks[i]->getEnd();
            (*allExpiries)[i] = ExpirySP(const_cast<Expiry*>(expiry));
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

void BootstrappedYieldCurve::acceptValueDateCollector(
    BootstrappedYieldCurve*       yieldCurve,
    CValueDateCollector* collector)
{
    collector->valueDateValidate(yieldCurve->today,
                                 yieldCurve->getName());
}

void BootstrappedYieldCurve::acceptWrapperNameCollector(BootstrappedYieldCurve* yc,
                                               WrapperNameCollector* collector)
{
    collector->addName(yc->getName());
}

void BootstrappedYieldCurve::acceptYieldNameCollector(BootstrappedYieldCurve* yc,
                                             YieldNameCollector* collector)
{
    collector->addName(yc->getName());
}

/** Records name of iso code against yc name in market data object.
    It is invoked ONCE only
    - immediately after this object is placed in the cache. */
void BootstrappedYieldCurve::initialise(MarketData* market)
{
    market->setYieldCurveISOCode(name, ccy);
}

/** populate from market cache */
void BootstrappedYieldCurve::getMarket(const IModel* model, const MarketData* market)
{
    if (today.empty())
    {
        today = market->GetReferenceDate();
    }

    hols.getData(model, market);
    if (!ccyBasis.isEmpty())
    {
        ccyBasis.getData(model, market);
    }

    // if specified, ask for the vol
    if (!irVol.isEmpty())
    {
        // NB This might result in a null vol if the model reckons we don't
        // need it
        irVol.getData(model, market);
    }

    // calculate transient fields
    valueDate = hols->addBusinessDays(today, spotOffset);
    expiryCache();
}

/** return either the growth or discount curve  */
ZeroCurveSP BootstrappedYieldCurve::get(bool growthCurve) const
{
    static const string method = "BootstrappedYieldCurve::get";

    try
    {
        //TODO : could use ZCAccessor here to have "local" caching
        ZeroPairConstSP zp = CashSwapCurveGlobalCache::get(this);

        // clone the ZeroCurve to avoid modifying data in cache !
        const ZeroCurve* zc = zp->get(growthCurve);
        return ZeroCurveSP(dynamic_cast<ZeroCurve*>(zc->clone()));
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


CashFlowArraySP BootstrappedYieldCurve::getRatesAndDates() const
{
    static const string method = "BootstrappedYieldCurve::getRatesAndDates";

    try
    {
        const ZeroCurve& zcurve = zc.get(*this);
        return zcurve.getRatesAndDates();
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


//-----------------------------
//Duration::IParhandler methods
//-----------------------------

//return benchmarks for par instruments
//these will determine the durations to calculate unless specified
ExpiryArrayConstSP BootstrappedYieldCurve::getParBenchmarks() const
{
    ExpiryArraySP result(new ExpiryArray());
    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        if (benchmarks[i]->isCash() || benchmarks[i]->isSwap())
        {
            result->push_back(ExpirySP(const_cast<Expiry*>(benchmarks[i]->getEnd())));
        }
    }

    return ExpiryArrayConstSP(result.release());
}

//return closed form solution
//must be implemented if supportsDurationClosedForm = true
ExpiryResultArraySP BootstrappedYieldCurve::getDuration(const Duration* durationObj) const
{
    static const string method = "BootstrappedYieldCurve::getDuration";

    try
    {
        ExpiryArrayConstSP benchmarks = durationObj->getBenchmarks();
        if (!benchmarks) {
            benchmarks = getParBenchmarks();
        }

        //build storage for the duration results
        ExpiryResultArraySP durCalcs(new ExpiryResultArray(0));

        //for each tweakPoint
        for(int j=0;j<benchmarks->size();j++) {
            //get duration for this benchmark
            ExpirySP maturity((*benchmarks)[j]);

            //find maturity in BootstrappedYieldCurve expiries
            //-throws exception if not found
            int i = maturity->search(allExpiries.get());

            //store maturity date of this par instrument;
            DateTime matDate = this->benchmarks[i]->getBenchmarkDate(valueDate);
            //holiday adjust TODO

            double duration;

            if (this->benchmarks[i]->isCash())
            {
                //equivalent to ALIB_MM_SENS
                double yearFraction = floatDcc->years(valueDate, matDate);
                duration = yearFraction / (1 + (this->benchmarks[i]->getRate() * yearFraction));
            }
            else if (this->benchmarks[i]->isFuture())
            {
                continue;
            }
            else if (this->benchmarks[i]->isTurn())
            {
                throw ModelException(method, "Curve cannot have adjustments");
            }
            else if (this->benchmarks[i]->isSwap())
            {
                //clone this untweaked curve
                BootstrappedYieldCurveSP aCloneCSC(cloneYieldCurve(true));
                aCloneCSC->setProjectionCurve(true);

                //shift the appropriate rate by 1 bp
                aCloneCSC->benchmarks[i]->setRate(aCloneCSC->benchmarks[i]->getRate() + 0.0001);

                //equivalent to ALIB_BOND_DUR_EFF
                //with annual frequency and tweaked zero rate
                double zeroRate = aCloneCSC->zero(matDate);
                //always act/365F year fraction
                double yearFraction = matDate.daysDiff(valueDate)/365.0;
                duration = (1-pow(1+zeroRate,-yearFraction)) / aCloneCSC->benchmarks[i]->getRate();
            }

            //create result element
            ExpiryResult bmarkDuration(maturity, duration);
            //and store
            durCalcs->push_back(bmarkDuration);
        }

        return durCalcs;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

//return a means of tweaking the MarketObject in a pointwise manner
//assumption is that the results of this twaeker is an ExpiryResultArray.....
VectorRiskPropertySensitivityConstSP BootstrappedYieldCurve::getPointwiseTweaker() const
{
    //1 bp shift explictly required
    return VectorRiskPropertySensitivityConstSP(new RhoPointwise(0.0001));
}

//return a par instrument for the specified maturity
InstrumentSP BootstrappedYieldCurve::getParInstrument(const ExpirySP maturity) const
{
    static const string method = "BootstrappedYieldCurve::getParInstrument";

    try
    {
        //find maturity in BootstrappedYieldCurve expiries
        //-throws exception if not found
        int i = maturity->search(allExpiries.get());

        //store maturity date of this par instrument;
        DateTime matDate = benchmarks[i]->getBenchmarkDate(valueDate);

        CashFlowArray cfl(0);
        if (benchmarks[i]->isCash())
        {
            if (!floatDcc)
            {
                throw ModelException(method, "swap floating leg day count convention not supplied");
            }

            //initial notional exchange at valueDate
            CashFlow cfs(valueDate, -1);
            cfl.push_back(cfs);

            //final notional exchange 1 + coupon
            CashFlow cfm(matDate, 1 + (benchmarks[i]->getRate() * floatDcc->years(valueDate, matDate)));
            cfl.push_back(cfm);
        }
        else if (benchmarks[i]->isSwap())
        {
            //multiple cashflows occurring at swapFrequency intervals
            CashFlow  cf;
            DateTime  cfDate;
            DateTime  prevCfDate;    //store previous date for accrual period calculations
            double    cfAmount;

            if (!fixedDcc)
            {
                throw ModelException(method, "swap fixed leg day count convention not supplied");
            }

            //initial notional exchange at valueDate
            cf.date = valueDate;
            cf.amount = -1;
            cfl.push_back(cf);
            prevCfDate = valueDate;

            //regular cashflows
            bool done = false;
            int count = 1;
            while (!done)
            {
                //cashflow date is regular from valueDate
                cfDate = MaturityPeriod::toDate(count * 12 / fixedIvl->approxAnnualFrequency(), "M", valueDate);

                //up to but not including the final cashflow date
                if (cfDate < matDate)
                {
                    //rate * accrual period year fraction * notional
                    cfAmount = benchmarks[i]->getRate() * fixedDcc->years(prevCfDate, cfDate); //notional of 1
                    cf.date = cfDate;
                    cf.amount = cfAmount;

                    //add to the list & continue loop
                    cfl.push_back(cf);
                    prevCfDate = cfDate;
                    count++;
                }
                else
                {
                    done = true;
                }
            }

            //final notional exchange 1 + coupon
            cfAmount = 1 + (benchmarks[i]->getRate() * fixedDcc->years(prevCfDate, matDate)); //notional of 1
            cf.date = matDate;
            cf.amount = cfAmount;
            cfl.push_back(cf);
        }

        SimpleCashFlowStreamSP cfs(new SimpleCashFlowStream(&cfl, name));
        return cfs;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

bool
BootstrappedYieldCurve::isIndexCurve() const
{
    return isIndexCurve_;
}

/*
 * Reflection and addin support.
 */

/** Invoked when Class is 'loaded' */
/* static */
void BootstrappedYieldCurve::load(CClassSP& clazz)
{
    REGISTER(BootstrappedYieldCurve, clazz);
    SUPERCLASS(YieldCurve);
    IMPLEMENTS(IDeterministicYieldCurve);
    IMPLEMENTS(IPrivateObject);
    IMPLEMENTS(IGetMarket);
    IMPLEMENTS(IRestorableWithRespectTo<RateParallel>);
    IMPLEMENTS(IRestorableWithRespectTo<RatePointwise>);
    IMPLEMENTS(IRestorableWithRespectTo<IRRatePointwise>);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(IRiskyCurve);
    IMPLEMENTS(Duration::IParHandlerWithClosedForm);
    IMPLEMENTS(Duration::IParHandlerWithoutClosedForm);
    IMPLEMENTS(YCWeightedShift::IShift);
    EMPTY_SHELL_METHOD(defaultBootstrappedYieldCurve);
    FIELD(ccy,                  "currency name");
    FIELD(name,                 "rate index");
    FIELD(today,                "today");
    FIELD(spotOffset,           "spot offset");
    FIELD(hols,                 "holidays");
    FIELD       (moneyMarketDayCount,  "money market day count");
    FIELD(benchmarks,           "benchmarks");
    FIELD       (zcMethod,             "zero curve methodology");
    FIELD(ccyBasis,             "currency basis");
    FIELD(irVol,                "Interest rate volatility");
    FIELD       (futMaturity,          "length of future");
    FIELD       (badDayConvention,     "swap date adjustment convention");
    FIELD       (fixedDcc,             "swap fixed leg day count convention");
    FIELD       (fixedIvl,             "swap fixed leg period");
    FIELD       (floatDcc,             "swap floating leg day count convention");
    FIELD       (floatIvl,             "swap floating leg period");
    FIELD       (fixDates,             "fixing dates");
    FIELD       (fixRates,             "fixing rates");
    FIELD       (basisDcc,             "swap basis leg day count convention");
    FIELD       (basisIvl,             "swap basis leg period");
    FIELD       (isIndexCurve_,         "act as index curve");

    FIELD_MAKE_OPTIONAL(today);
    FIELD_MAKE_OPTIONAL(ccyBasis);
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_MAKE_OPTIONAL(futMaturity);
    FIELD_MAKE_OPTIONAL(badDayConvention);
    FIELD_MAKE_OPTIONAL(floatDcc);
    FIELD_MAKE_OPTIONAL(floatIvl);
    FIELD_MAKE_OPTIONAL(basisDcc);
    FIELD_MAKE_OPTIONAL(basisIvl);
    FIELD_MAKE_OPTIONAL(isIndexCurve_);

    // transient fields
    FIELD_NO_DESC       (allExpiries);
    FIELD_NO_DESC(useProjectionCurve);
    FIELD_NO_DESC(nosort);

    FIELD_MAKE_TRANSIENT(allExpiries);
    FIELD_MAKE_TRANSIENT(useProjectionCurve);
    FIELD_MAKE_TRANSIENT(nosort);

    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptValueDateCollector);
    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptWrapperNameCollector);
    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptYieldNameCollector);
    clazz->setPrivate(); // hide this class
    // record how to build it (since implements IPrivateObject)
    clazz->addConstructorClass(BootstrappedCurveAddin::TYPE);
}

/* static */
IObject* BootstrappedYieldCurve::defaultBootstrappedYieldCurve()
{
    return new BootstrappedYieldCurve();
}

CClassConstSP const BootstrappedYieldCurve::TYPE = CClass::registerClassLoadMethod(
    "BootstrappedYieldCurve", typeid(BootstrappedYieldCurve), BootstrappedYieldCurve::load);


// Addin for bootstrapping yield curves

/* static */
IObject* BootstrappedCurveAddin::defaultBootstrappedCurveAddin()
{
    return new BootstrappedCurveAddin();
}


/** Invoked when Class is 'loaded' */
/* static */
void BootstrappedCurveAddin::load(CClassSP& clazz)
{
    clazz->setDRIProxyType(BootstrappedYieldCurve::TYPE); //use BootstrappedYieldCurve for dri
    REGISTER(BootstrappedCurveAddin, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IPublicObject);
    EMPTY_SHELL_METHOD(defaultBootstrappedCurveAddin);
    FIELD(ccy,              "currency name");
    FIELD(name,             "yield curve name");
    FIELD(today,            "today");
    FIELD(spotOffset,       "spot offset");
    FIELD(hols,             "holidays");
    FIELD(benchmarks,       "benchmarks");
    FIELD       (factory,          "bootstrapping methodology");
    FIELD(moneyMarketDcc,   "money market day count");
    FIELD       (futuresMaturity,  "maturity of futures rates, eg. 3M or 1I");
    FIELD(irVol,            "Interest rate volatility");
    FIELD(badDayConvention, "swap bad day convention");
    FIELD(fixedDcc,         "swap fixed leg day count convention");
    FIELD       (fixedIvl,         "swap fixed leg period");
    FIELD(floatDcc,         "swap floating leg day count convention");
    FIELD       (floatIvl,         "swap floating leg period");
    FIELD       (fixDates,         "fixing dates");
    FIELD       (fixRates,         "fixing rates");
    FIELD(basisDcc,         "swap basis leg day count convention");
    FIELD       (basisIvl,         "swap basis leg period");
    FIELD(ccyBasis,         "currency basis");
    FIELD(isIndexCurve,     "whether to act as index curve");

    FIELD_MAKE_OPTIONAL(today);
    FIELD_MAKE_OPTIONAL(futuresMaturity);
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_MAKE_OPTIONAL(badDayConvention);
    FIELD_MAKE_OPTIONAL(floatDcc);
    FIELD_MAKE_OPTIONAL(floatIvl);
    FIELD_MAKE_OPTIONAL(fixDates);
    FIELD_MAKE_OPTIONAL(fixRates);
    FIELD_MAKE_OPTIONAL(basisDcc);
    FIELD_MAKE_OPTIONAL(basisIvl);
    FIELD_MAKE_OPTIONAL(ccyBasis);
    FIELD_MAKE_OPTIONAL(isIndexCurve);

    Addin::registerConstructor("BOOTSTRAP_YIELD_CURVE",
                               Addin::MARKET,
                               "Bootstrap a zero curve",
                               BootstrappedCurveAddin::TYPE);
}


BootstrappedCurveAddin::BootstrappedCurveAddin() : CObject(TYPE), isIndexCurve(false)
{
}


BootstrappedCurveAddin::BootstrappedCurveAddin(
    const string&                  ccy,
    const string&                  name,
    const DateTime&                today,
    int                            spotOffset,
    const HolidayWrapper&          hols,
    const ZeroCurveBenchmarkArray& benchmarks,
    const IZeroCurveFactory&       factory,
    const string&                  moneyMarketDcc,   // may be ""
    const MaturityPeriod*          futuresMaturity,
    const IRVolBaseWrapper&        irVol,
    const string&                  badDayConvention, // may be ""
    const string&                  fixedDcc,         // may be ""
    const MaturityPeriod*          fixedIvl,
    const string&                  floatDcc,         // may be ""
    const MaturityPeriod*          floatIvl,
    const ExpiryArray*             fixDates,
    const DoubleArray*             fixRates,
    const string&                  basisDcc,         // may be ""
    const MaturityPeriod*          basisIvl,
    const CurrencyBasisWrapper&    ccyBasis,
    bool                           isIndexCurve)
    : CObject(TYPE),
    ccy(ccy), name(name), today(today), spotOffset(spotOffset), hols(hols),
    benchmarks(benchmarks),
    factory(const_cast<IZeroCurveFactory*>(&factory)),
    moneyMarketDcc(moneyMarketDcc),
    futuresMaturity(futuresMaturity ? dynamic_cast<MaturityPeriod*>(futuresMaturity->clone()) : NULL),
    irVol(irVol),
    badDayConvention(badDayConvention),
    fixedDcc(fixedDcc),
    fixedIvl(fixedIvl ? dynamic_cast<MaturityPeriod*>(fixedIvl->clone()) : NULL),
    floatDcc(floatDcc),
    floatIvl(floatIvl ? dynamic_cast<MaturityPeriod*>(floatIvl->clone()) : NULL),
    fixDates(fixDates ? dynamic_cast<ExpiryArray*>(fixDates->clone()) : NULL),
    fixRates(fixRates ? dynamic_cast<DoubleArray*>(fixRates->clone()) : NULL),
    basisDcc(basisDcc),
    basisIvl(basisIvl ? dynamic_cast<MaturityPeriod*>(basisIvl->clone()) : NULL),
    ccyBasis(ccyBasis),
    isIndexCurve(isIndexCurve)
{
}


IPrivateObject* BootstrappedCurveAddin::toPrivateObject() const
{
    static const string method = "BootstrappedCurveAddin::toPrivateObject";

    try
    {
        if (!factory)
        {
            throw ModelException(method, "bootstrapping methodology not defined");
        }

        if (moneyMarketDcc.empty())
        {
            throw ModelException(method, "Money market day count convention must be defined");
        }

        if (fixedDcc.empty())
        {
            throw ModelException(method, "Swap fixed day count convention must be defined");
        }

        if (!fixedIvl.get())
        {
            throw ModelException(method, "Swap fixed interval must be defined");
        }

        return new BootstrappedYieldCurve(ccy,
                                 name,
                                 today,
                                 spotOffset,
                                 hols,
                                 benchmarks,
                                 *factory,
                                 *getDcc(moneyMarketDcc),
                                 futuresMaturity.get(),
                                 irVol,
                                 getBdc(badDayConvention).get(),
                                 *getDcc(fixedDcc),
                                 *fixedIvl,
                                 getDcc(floatDcc).get(),
                                 floatIvl.get(),
                                 fixDates.get(),
                                 fixRates.get(),
                                 getDcc(basisDcc).get(),
                                 basisIvl.get(),
                                 ccyBasis,
                                 isIndexCurve);
    }
    catch (exception &e)
    {
        throw ModelException(e, method, e.what());
    }
}


BadDayConventionSP  BootstrappedCurveAddin::getBdc(const string& bdc) const
{
    return BadDayConventionSP(bdc.empty() ? NULL : BadDayConventionFactory::make(bdc));
}


DayCountConventionSP  BootstrappedCurveAddin::getDcc(const string& dcc) const
{
    return DayCountConventionSP(dcc.empty() ? NULL : DayCountConventionFactory::make(dcc));
}


CClassConstSP const BootstrappedCurveAddin::TYPE =
CClass::registerClassLoadMethod("BootstrappedCurveAddin", typeid(BootstrappedCurveAddin),
                                BootstrappedCurveAddin::load);


// other addins

class ForwardCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    DateTime            forwardDate;
    BootstrappedYieldCurveSP yieldCurve;

    static IObjectSP getForwardSwapCurve(ForwardCurveAddin* params) {
        static const string routine = "ForwardCurveAddin::getForwardSwapCurve";
        try {
            IYieldCurveSP fwdCurve = IYieldCurveSP(params->yieldCurve->createForwardCurve(params->forwardDate));
            return fwdCurve;

        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ForwardCurveAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ForwardCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultForwardCurveAddin);
        FIELD(forwardDate, "forward spot date of the yield curve");
        FIELD(yieldCurve,         "spot yield curve");
        Addin::registerClassObjectMethod("GET_FORWARD_YIELD_CURVE",
                                         Addin::MARKET,
                                         "Returns the forward yield curve for a given date",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getForwardSwapCurve);

    }

    static IObject* defaultForwardCurveAddin(){
        return new ForwardCurveAddin();
    }

};

CClassConstSP const ForwardCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "ForwardCurveAddin", typeid(ForwardCurveAddin), load);

class ZeroCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    BootstrappedYieldCurveSP     yieldCurve;
    bool                extractGrowth;

    static IObjectSP getZeroCurve(ZeroCurveAddin* params) {
        static const string routine = "ZeroCurveAddin::getZeroCurve";
        try {
            ZeroCurveSP     zc        = params->yieldCurve->get(params->extractGrowth);
            CashFlowArraySP cashFlows = zc->getRatesAndDates();


            DateTime::DateArraySP dateArray(new DateTime::DateArray(0));
            DoubleArraySP         doubleArray(new DoubleArray(0));

            for (int i = 0; i < cashFlows->size(); i++) {
                dateArray->push_back(DateTime::DateSP(new DateTime::Date((*cashFlows)[i].date.getDate())));
                doubleArray->push_back((*cashFlows)[i].amount);
            }

            ObjectArraySP objArray(new ObjectArray(0));

            objArray->push_back(dateArray);
            objArray->push_back(doubleArray);

            return objArray;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ZeroCurveAddin():  CObject(TYPE), extractGrowth(true) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ZeroCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZeroCurveAddin);
        FIELD(yieldCurve,         "spot yield curve");
        FIELD(extractGrowth, "return growth curve or discount");
        FIELD_MAKE_OPTIONAL(extractGrowth);
        Addin::registerClassObjectMethod("ZERO_CURVE_EXTRACT",
                                         Addin::MARKET,
                                         "Returns the zero curve for a yield curve",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getZeroCurve);

    }

    static IObject* defaultZeroCurveAddin(){
        return new ZeroCurveAddin();
    }

};

CClassConstSP const ZeroCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "ZeroCurveAddin", typeid(ZeroCurveAddin), load);


DRLIB_END_NAMESPACE
