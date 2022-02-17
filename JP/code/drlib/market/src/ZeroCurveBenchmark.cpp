//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZeroCurveBenchmark.cpp
//
//   Description : Input benchmarks for bootstrapping curve.
//
//   Author      : Richard Appleton
//
//   Date        : 16th November 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ZeroCurveBenchmark.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Maths.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Equals.hpp"
#include "edginc/Addin.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Hashtable.hpp" // for hash_string
#include "edginc/BenchmarkDate.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include <algorithm>



DRLIB_BEGIN_NAMESPACE
using namespace std;


/* static */
ZeroCurveBenchmarkArraySP ZeroCurveBenchmark::toArray(
    const ExpiryArray* starts,
    const ExpiryArray& expiries,
    const StringArray& instruments,
    const IntArray*    includeFlags,
    const DoubleArray& rates,
    const DoubleArray* prices,
    const DoubleArray* adjustments,
    const DateTime*    refDate,
    const Holiday*     holidays)
{
    static const string method = "ZeroCurveBenchmark::toArray";

    // verify input array sizes

    if (starts && starts->size() > 0 && instruments.size() != starts->size())
    {
        string msg = Format::toString(
            "Inconsistent start date array size [%d, expected %d]", 
            starts->size(), instruments.size());
        throw ModelException(method, msg);
    }

    if (instruments.size() != expiries.size())
    {
        string msg = Format::toString( 
            "Inconsistent end date array size [%d, expected %d]",
            expiries.size(), instruments.size());
        throw ModelException(method, msg);
    }

    if (includeFlags && includeFlags->size() > 0 && instruments.size() != includeFlags->size())
    {
        string msg = Format::toString( 
            "Inconsistent include flags array size [%d, expected %d]",
            includeFlags->size(), instruments.size());
        throw ModelException(method, msg);
    }

    if (instruments.size() != rates.size())
    {
        string msg = Format::toString(
            "Inconsistent rate array size [%d, expected %d]",
            rates.size(), instruments.size());
        throw ModelException(method, msg);
    }

    if (prices && prices->size() > 0 && instruments.size() != prices->size())
    {
        string msg = Format::toString( 
            "Inconsistent prices array size [%d, expected %d]",
            prices->size(), instruments.size());
        throw ModelException(method, msg);
    }

    if (adjustments && adjustments->size() > 0 && instruments.size() != adjustments->size())
    {
        string msg = Format::toString(
            "Inconsistent adjustments array size [%d, expected %d]",
            adjustments->size(), instruments.size());
        throw ModelException(method, msg);
    }

    ZeroCurveBenchmarkArraySP benchmarks(new ZeroCurveBenchmarkArray());
    benchmarks->reserve(instruments.size());

    for (int i = 0 ; i < instruments.size() ; i++)
    {
        append(
            *benchmarks, 
            i, 
            starts, 
            expiries, 
            instruments, 
            includeFlags, 
            rates, 
            prices, 
            adjustments,
            refDate,
            holidays);
    }

    return benchmarks;
}


/* static */
void ZeroCurveBenchmark::append(
    ZeroCurveBenchmarkArray& benchmarks, 
    int                      i, 
    const ExpiryArray*       starts,
    const ExpiryArray&       expiries,
    const StringArray&       instruments,
    const IntArray*          includeFlags,
    const DoubleArray&       rates,
    const DoubleArray*       prices,
    const DoubleArray*       adjustments,
    const DateTime*          refDate,
    const Holiday*           holidays)
{
    static const string method = "ZeroCurveBenchmark::toArray";

    int includeFlag = includeFlags && includeFlags->size() > 0 ? (*includeFlags)[i] : 1;

    // NULL instruments are allowed if we are skipping this particular instrument
    if (includeFlag == 0)
    {
        return;
    }
    else if (instruments[i].empty())
    {
        string msg = Format::toString("NULL instrument for element %d", i);
        throw ModelException(method, msg);
    }
    else
    {
        ZeroCurveBenchmark* benchmark;
        const Expiry* start = starts && starts->size() > 0 ? (*starts)[i].get() : NULL;
        double price = prices && prices->size() > 0 ? (*prices)[i] : 0.0;

        if (CString::equalsIgnoreCase(instruments[i], YieldCurve::MMRT_RATE,1)
            || CString::equalsIgnoreCase(instruments[i], YieldCurve::MMRT_RATE2,1))
        {
            benchmark = new MoneyMarketBenchmark(
                start,
                *expiries[i],
                rates[i],
                includeFlag,
                NULL);
        }
        else if (CString::equalsIgnoreCase(instruments[i], YieldCurve::FUTURE_RATE,1))
        {
            /*
                * For backwards compatibility if no start dates are provided
                * the future is assumed to start on the 'expiry' date and be 
                * futuresMaturity long
                */
            benchmark = new FuturesBenchmark(
                start ? *start : *expiries[i],
                start ? expiries[i].get() : NULL,
                !Maths::isZero(price) ? price : rates[i],
                !Maths::isZero(price),
                includeFlag,
                NULL,
                adjustments && adjustments->size() > 0 ? &(*adjustments)[i] : NULL,
                "");
        }
        else if (CString::equalsIgnoreCase(instruments[i], YieldCurve::BRAZIL_FUTURE_PRICE,1))
        {
            if( ! refDate || ! holidays )
            {
                throw ModelException( method,
                    "Reference date and holidays must be provided for brazil futures benchmark" );
            }

            benchmark = new BrazilFuturesBenchmark(
                *expiries[ i ],
                rates[ i ],
                true,
                includeFlag,
                *refDate,
                *holidays );
        }
        else if (CString::equalsIgnoreCase(instruments[i], YieldCurve::TURN_RATE,1))
        {
            benchmark = new TurnBenchmark(
                instruments[i].substr(1), 
                start, *expiries[i], rates[i], includeFlag, NULL);
        }
        else if (CString::equalsIgnoreCase(instruments[i], YieldCurve::SWAP_RATE,1)
                || CString::equalsIgnoreCase(instruments[i], YieldCurve::SWAP_RATE2,1))
        {
            benchmark = new SwapBenchmark(
                start,
                *expiries[i],
                rates[i],
                adjustments && adjustments->size() > 0 ? &(*adjustments)[i] : NULL,
                price,
                includeFlag,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL,
                NULL);
        }
        else
        {
            string msg = Format::toString("Invalid instrument type: %s", instruments[i].c_str());
            throw ModelException(method, msg);
        }

        for (int j = 0 ; j < benchmarks.size() ; j++)
        {
            if (benchmark->equals(*benchmarks[j]))
            {
                const Expiry* expiry = benchmark->getEnd();
                if (!expiry)
                {
                    expiry = benchmark->getStart();
                }

                string msg = Format::toString(
                    "Duplicate benchmarks exist for %s [benchmarks %d and %d]", 
                    expiry->toString().c_str(), j, i);
                throw ModelException(method, msg);
            }
        }

        benchmarks.push_back(ZeroCurveBenchmarkSP(benchmark));
    }
}


/**
 * Returns index of benchmark with specified expiry
 *
 * Uses end expiry, except for futures where a match is against the start expiry.
 */
/* static */
int ZeroCurveBenchmark::search(const ZeroCurveBenchmarkArray& benchmarks, const Expiry& expiry)
{
    static const string method = "ZeroCurveBenchmark::search";

    for (int i = 0; i < benchmarks.size(); i++)
    {
        const Expiry* end = benchmarks[i]->isFuture() ? 
            benchmarks[i]->getStart() : benchmarks[i]->getEnd();

        if (expiry.equals(end))
        {
            return i;
        }
    }

    string msg = Format::toString("Benchmark with expiry (%s) not found", expiry.toString().c_str());
    throw ModelException(method, msg);
}


struct ResolvedZeroCurveBenchmark
{
    ZeroCurveBenchmarkSP benchmark;
    DateTime             date;
};


// comparator for STL sort algorithm
class ZeroCurveBenchmarkComparator
{
public:
    bool operator()(const ResolvedZeroCurveBenchmark& lhs, const ResolvedZeroCurveBenchmark& rhs) const
    {
        if (lhs.date == rhs.date)
        {
            return !lhs.benchmark->isFuture();
        }
        else
        {
            return lhs.date < rhs.date;
        }
    }
};


/* static */
void ZeroCurveBenchmark::sort(ZeroCurveBenchmarkArray& benchmarks, const DateTime& valueDate)
{
    // toDate() is quite slow, so do it just once before sorting
    int length = benchmarks.size();
    vector<ResolvedZeroCurveBenchmark> data(length);
    for (int i = 0 ; i < length ; i++)
    {
        data[i].benchmark = benchmarks[i];
        data[i].date = benchmarks[i]->getBenchmarkDate(valueDate);
    }

    std::sort(data.begin(), data.end(), ZeroCurveBenchmarkComparator());

    // copy sorted data back to input array
    for (int j = 0 ; j < length ; j++)
    {
        benchmarks[j] = ZeroCurveBenchmarkSP(data[j].benchmark);
    }
}


/* static */
bool ZeroCurveBenchmark::equals(const ZeroCurveBenchmarkArray& lhs, const ZeroCurveBenchmarkArray& rhs)
{
    if (&lhs != &rhs)
    {
        int len = lhs.size();

        if (len != rhs.size())
            return false;

        for (int i = 0; i < len; i++)
        {
            if (!lhs[i]->equals(*rhs[i]))
                return false;
        }
    }

    return true;
}


/* static */
int ZeroCurveBenchmark::hashCode(const ZeroCurveBenchmarkArray& benchmarks)
{
    int hCode = 0;

    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        hCode ^= benchmarks[i]->hashCode();
    }

    return hCode;
}


ZeroCurveBenchmark::ZeroCurveBenchmark(const CClassConstSP& pClazz)
: CObject(pClazz), start(NULL), end(NULL), rate(0.0), includeFlag(1)
{
}


ZeroCurveBenchmark::ZeroCurveBenchmark(
    const Expiry*        pStart, 
    const Expiry*        pEnd,
    double               pRate, 
    int                  pIncludeFlag,
    const CClassConstSP& pClazz)
    : CObject(pClazz),
      start(pStart),
      end(pEnd), 
      rate(pRate),
      includeFlag(pIncludeFlag)
{
}


ZeroCurveBenchmark::~ZeroCurveBenchmark()
{
}


int ZeroCurveBenchmark::hashCode() const
{
    int hCode = 0;
    hCode ^= CDouble::hashCode(rate);
    return hCode;
}


bool ZeroCurveBenchmark::equals(const ZeroCurveBenchmark& other) const
{
    if (typeid(*this) != typeid(other))
        return false;

    if (!equalsValue<Expiry>(start, other.start))
        return false;

    if (!equalsValue<Expiry>(end, other.end))
        return false;

    if (!Maths::equals(rate,other.rate))
        return false;

    if (includeFlag != other.includeFlag)
        return false;

    return true;
}


bool ZeroCurveBenchmark::isCash() const
{
    return false;
}


bool ZeroCurveBenchmark::isFuture() const
{
    return false;
}


bool ZeroCurveBenchmark::isTurn() const
{
    return false;
}


bool ZeroCurveBenchmark::isSwap() const
{
    return false;
}


const Expiry* ZeroCurveBenchmark::getStart() const
{
    return start.get();
}


const Expiry* ZeroCurveBenchmark::getEnd() const
{
    return end.get();
}


DateTime ZeroCurveBenchmark::getBenchmarkDate(const DateTime& date) const
{
    return end->toDate(date);
}


double ZeroCurveBenchmark::getRate() const
{
    return rate;
}


int ZeroCurveBenchmark::getIncludeFlag() const
{
    return includeFlag;
}


bool ZeroCurveBenchmark::hasAdjustment() const
{
    return false;
}


double ZeroCurveBenchmark::getAdjustment() const
{
    return 0.0;
}


void ZeroCurveBenchmark::setRate(double pRate)
{
    rate = pRate;
}


MoneyMarketBenchmark::MoneyMarketBenchmark()
: ZeroCurveBenchmark(TYPE), dcc(NULL)
{
}


MoneyMarketBenchmark::MoneyMarketBenchmark(
    const Expiry*              pStart, 
    const Expiry&              pEnd,
    double                     pRate, 
    int                        pIncludeFlag,
    const DayCountConvention*  pDcc)
    : ZeroCurveBenchmark(pStart, &pEnd, pRate, pIncludeFlag, TYPE), dcc(pDcc)
{
    if (includeFlag != 1)
    {
        string msg = Format::toString(
            "Invalid benchmark for cash instrument [%s to %s]",
            start.get() ? start->toString().c_str() : "spot", 
            end->toString().c_str());
        throw ModelException("MoneyMarketBenchmark::MoneyMarketBenchmark", msg);
    }
}


bool MoneyMarketBenchmark::isCash() const
{
    return true;
}


IObject* MoneyMarketBenchmark::clone() const
{
    MoneyMarketBenchmark* cloned = new MoneyMarketBenchmark();
    cloned->start = start;
    cloned->end = end;
    cloned->rate = rate;
    cloned->includeFlag = includeFlag;
    cloned->dcc = dcc;
    return cloned;
}


int MoneyMarketBenchmark::hashCode() const
{
    int hCode = ZeroCurveBenchmark::hashCode();

    if (dcc.get())
    {
        hCode ^= dcc->hashCode();
    }

    return hCode;
}


bool MoneyMarketBenchmark::equals(const ZeroCurveBenchmark& other) const
{
    if (!ZeroCurveBenchmark::equals(other))
        return false;

    const MoneyMarketBenchmark& rhs = static_cast<const MoneyMarketBenchmark&>(other);

    if (!equalsClass<DayCountConvention>(dcc, rhs.dcc))
        return false;

    return true;
}


const DayCountConvention* MoneyMarketBenchmark::getDcc() const
{
    return dcc.get();
}


const double FuturesBenchmark::ZERO_PERCENT = 10000.0;


FuturesBenchmark::FuturesBenchmark()
: ZeroCurveBenchmark(TYPE), type(""), adjustment(0.0), dcc(NULL)
{
}


FuturesBenchmark::FuturesBenchmark(
    const Expiry&             pStart, 
    const Expiry*             pEnd,
    double                    pQuote, 
    bool                      pQuotedAsPrice,
    int                       pIncludeFlag,
    const DayCountConvention* pDcc,
    const double*             pAdjustment,
    const string&             pType)
    : ZeroCurveBenchmark(&pStart, pEnd, 
        (pQuotedAsPrice ? (ZERO_PERCENT - pQuote) / ZERO_PERCENT : pQuote), 
        pIncludeFlag, TYPE),
    type(pType), adjustment(pAdjustment ? *pAdjustment : 0.0), dcc(pDcc)
{
    static const string method = "FuturesBenchmark::FuturesBenchmark";

    if (includeFlag != 1 && includeFlag != 2)
    {
        string msg = Format::toString(
            "Invalid benchmark for futures instrument [%s to %s]",
            start->toString().c_str(),
            end.get() ? end->toString().c_str() : "-");
        throw ModelException(method, msg);
    }
}


DateTime FuturesBenchmark::getBenchmarkDate(const DateTime& date) const
{
    return start->toDate(date);
}


bool FuturesBenchmark::isCash() const
{
    return false;
}


bool FuturesBenchmark::isFuture() const
{
    return true;
}


void FuturesBenchmark::copyTo(FuturesBenchmark* copy) const
{
    copy->start = start;
    copy->end = end;
    copy->rate = rate;
    copy->includeFlag = includeFlag;
    copy->type = type;
    copy->adjustment = adjustment;
    copy->dcc = dcc;
}


IObject* FuturesBenchmark::clone() const
{
    FuturesBenchmark* clone = new FuturesBenchmark;
    copyTo(clone);
    return clone;
}


int FuturesBenchmark::hashCode() const
{
    int hCode = ZeroCurveBenchmark::hashCode();
    hCode ^= hash_string(type);     // avoid creating CString object
    hCode ^= CDouble::hashCode(adjustment);

    if (dcc.get())
    {
        hCode ^= dcc->hashCode();
    }

    return hCode;
}


bool FuturesBenchmark::equals(const ZeroCurveBenchmark& other) const
{
    if (!ZeroCurveBenchmark::equals(other))
        return false;

    const FuturesBenchmark& rhs = static_cast<const FuturesBenchmark&>(other);

    if (!equalsClass<DayCountConvention>(dcc, rhs.dcc))
        return false;

    if (adjustment != rhs.adjustment)
        return false;

    if (type != rhs.type)
        return false;

    return true;
}


const string& FuturesBenchmark::getType() const
{
    return type;
}


bool FuturesBenchmark::hasAdjustment() const
{
    return !Maths::isZero(adjustment);
}


double FuturesBenchmark::getAdjustment() const
{
    return adjustment;
}


const DayCountConvention* FuturesBenchmark::getDcc() const
{
    return dcc.get();
}


BrazilFuturesBenchmark::BrazilFuturesBenchmark() :
    FuturesBenchmark()
{}

BrazilFuturesBenchmark::BrazilFuturesBenchmark(
    const Expiry&   expiry,
    double          quote,
    bool            quotedAsPrice,
    int             includeFlag,
    const DateTime& refDate,
    const Holiday&  holidays)
    :
    FuturesBenchmark(
        expiry,
        NULL,
        quote,
        quotedAsPrice,
        includeFlag,
        NULL,
        NULL,
        "")
{
    if( quotedAsPrice )
    {
        const DateTime & date = adjustExpiry( *getStart(), refDate, holidays )->toDate( refDate );

        // convert price to rate
        setRate(
            pow(
                100000. / quote,
                252. / holidays.businessDaysDiff( refDate, date ) )
            - 1. );
    }
}

double BrazilFuturesBenchmark::getPrice(
    const DateTime& refDate,
    const Holiday&  holidays) const
{
    const DateTime & date = adjustExpiry( *getStart(), refDate, holidays )->toDate( refDate );

    // convert rate to price
    return
        100000. /
        pow(
            1. + getRate(),
            holidays.businessDaysDiff( refDate, date ) / 252. );
}

ExpiryConstSP BrazilFuturesBenchmark::adjustExpiry(
    const Expiry&   expiry,
    const DateTime& refDate,
    const Holiday&  holidays)
{
    static const string method = "BrazilFuturesBenchmark::adjustExpiry";

    const BenchmarkDate * benchmarkDate =
        dynamic_cast< const BenchmarkDate * >( &expiry );

    if( benchmarkDate )
    {
        DateTime date = benchmarkDate->toDate();
        DateTime prevDate = holidays.addBusinessDays( date, -1 );
        if( prevDate.toMDY().month == date.toMDY().month )
        {
            throw ModelException(method,
                Format::toString("Expiry date for futures benchmark should be the first business day of the month, '%s' specified", date.toString().c_str()));
        }

        return ExpiryConstSP( benchmarkDate );
    }
    else
    {
        const MaturityTimePeriod * maturityTimePeriod =
            dynamic_cast< const MaturityTimePeriod * >( &expiry );
        const MaturityPeriod * maturityPeriod = maturityTimePeriod ?
            maturityTimePeriod->getMaturityPeriod().get() :
            dynamic_cast< const MaturityPeriod * >( &expiry );

        if( ! maturityPeriod )
        {
            throw ModelException(method,
                Format::toString("Unsupported expiry type '%s' for futures benchmark", expiry.getClass()->getName().c_str()));
        }

        int count;
        string interval;
        maturityPeriod->decompose( count, interval );
        if( interval[0] != 'J' )
        {
            throw ModelException(method,
                Format::toString("Expiry interval for futures benchmark should be 'J', '%c' specified", interval[0]));
        }

        return ExpirySP( new BenchmarkDate(
            holidays.addBusinessDays(
                refDate.rollDateInMonths( count - 1 ).returnEndOfMonth( false ), 1 ) ) );
    }
}

IObject* BrazilFuturesBenchmark::clone() const
{
    BrazilFuturesBenchmark* clone = new BrazilFuturesBenchmark;
    copyTo(clone);
    return clone;
}


const string TurnBenchmark::ABSOLUTE = "A";
const string TurnBenchmark::IMPLIED  = "I";
const string TurnBenchmark::RELATIVE = "R";
const string TurnBenchmark::UTURN    = "U";


TurnBenchmark::TurnBenchmark()
: ZeroCurveBenchmark(TYPE), type(RELATIVE), dcc(NULL)
{
}


TurnBenchmark::TurnBenchmark(
    const string&             pType,
    const Expiry*             pStart, 
    const Expiry&             pEnd,
    double                    pRate,
    int                       pIncludeFlag,
    const DayCountConvention* pDcc)
  : ZeroCurveBenchmark(pStart, &pEnd, pRate, pIncludeFlag, TYPE), type(pType), dcc(pDcc)
{
    static const string method = "TurnBenchmark::TurnBenchmark";

    if (includeFlag != 1)
    {
        string msg = Format::toString(
            "Invalid benchmark for turn [%s to %s]",
            start.get() ? start->toString().c_str() : "spot", 
            end->toString().c_str());
        throw ModelException(method, msg);
    }

    // turn type defaults to 'relative'
    if (type.empty())
    {
        type = RELATIVE;
    }

    if (CString::equalsIgnoreCase(type,ABSOLUTE,1)
     || CString::equalsIgnoreCase(type,IMPLIED,1)
     || CString::equalsIgnoreCase(type,RELATIVE,1))
    {
        // valid types
    }
    else if (CString::equalsIgnoreCase(type,UTURN,1))
    {
        throw ModelException(method, "U_turns are not yet implemented");
    }
    else
    {
        // end-user interface prefixes turn type with 'A' in instrument type
        string msg = Format::toString("Invalid turn type: A%s", type.c_str());
        throw ModelException(method, msg);
    }
}


bool TurnBenchmark::isCash() const
{
    return false;
}


bool TurnBenchmark::isTurn() const
{
    return true;
}


IObject* TurnBenchmark::clone() const
{
    TurnBenchmark* cloned = new TurnBenchmark();
    cloned->start = start;
    cloned->end = end;
    cloned->rate = rate;
    cloned->includeFlag = includeFlag;
    cloned->type = type;
    cloned->dcc = dcc;
    return cloned;
}


int TurnBenchmark::hashCode() const
{
    int hCode = ZeroCurveBenchmark::hashCode();
    hCode ^= hash_string(type);     // avoid creating CString object

    if (dcc.get())
    {
        hCode ^= dcc->hashCode();
    }

    return hCode;
}


bool TurnBenchmark::equals(const ZeroCurveBenchmark& other) const
{
    if (!ZeroCurveBenchmark::equals(other))
        return false;

    const TurnBenchmark& rhs = static_cast<const TurnBenchmark&>(other);

    if (!equalsClass<DayCountConvention>(dcc, rhs.dcc))
        return false;

    if (type != rhs.type)
        return false;

    return true;
}


bool TurnBenchmark::isAbsolute() const
{
    return CString::equalsIgnoreCase(type, ABSOLUTE, 1);
}


bool TurnBenchmark::isRelative() const
{
    return CString::equalsIgnoreCase(type, RELATIVE, 1);
}


bool TurnBenchmark::isImplied() const
{
    return CString::equalsIgnoreCase(type, IMPLIED, 1);
}


const DayCountConvention* TurnBenchmark::getDcc() const
{
    return dcc.get();
}


SwapBenchmark::SwapBenchmark() 
: ZeroCurveBenchmark(TYPE), 
  price(1.0), 
  fixedIvl(NULL), fixedDcc(NULL), 
  floatIvl(NULL), floatDcc(NULL), 
  basisIvl(NULL), basisDcc(NULL), 
  adjustment(0.0)
{
}


SwapBenchmark::SwapBenchmark(
    const Expiry*             pStart, 
    const Expiry&             pEnd,
    double                    pRate, 
    const double*             pAdjustment,
    double                    pPrice,
    int                       pIncludeFlag,
    const MaturityPeriod*     pFixedIvl,
    const DayCountConvention* pFixedDcc,
    const MaturityPeriod*     pFloatIvl,
    const DayCountConvention* pFloatDcc,
    const MaturityPeriod*     pBasisIvl,
    const DayCountConvention* pBasisDcc)
    : ZeroCurveBenchmark(pStart, &pEnd, pRate, pIncludeFlag, TYPE), 
    price(Maths::isZero(pPrice) ? 1.0 : pPrice), 
    fixedIvl(pFixedIvl), 
    fixedDcc(pFixedDcc), 
    floatIvl(pFloatIvl), 
    floatDcc(pFloatDcc),
    basisIvl(pBasisIvl), 
    basisDcc(pBasisDcc),
    adjustment(pAdjustment ? *pAdjustment : 0.0)
{
    static const string method = "SwapBenchmark::SwapBenchmark";

    if (includeFlag < 1 || includeFlag > 2)
    {
        string msg = Format::toString(
            "[%s] Invalid swap benchmark", pEnd.toString().c_str());
        throw ModelException(method, msg);
    }

    if (includeFlag == 2 && (pAdjustment == NULL))
    {
        string msg = Format::toString(
            "[%s] Must provide adjustments for benchmark flag 2",
            pEnd.toString().c_str());
        throw ModelException(method, msg);
    }
}


bool SwapBenchmark::isSwap() const
{
    return true;
}


IObject* SwapBenchmark::clone() const
{
    SwapBenchmark* cloned = new SwapBenchmark();
    cloned->start = start;
    cloned->end = end;
    cloned->rate = rate;
    cloned->includeFlag = includeFlag;
    cloned->price = price;
    cloned->fixedIvl = fixedIvl;
    cloned->fixedDcc = fixedDcc;
    cloned->floatIvl = floatIvl;
    cloned->floatDcc = floatDcc;
    cloned->basisIvl = basisIvl;
    cloned->basisDcc = basisDcc;
    cloned->adjustment = adjustment;
    return cloned;
}


int SwapBenchmark::hashCode() const
{
    int hCode = ZeroCurveBenchmark::hashCode();
    hCode ^= CDouble::hashCode(price);
    hCode ^= CDouble::hashCode(adjustment);

    if (fixedIvl.get())
        hCode ^= fixedIvl->hashCode();

    if (fixedDcc.get())
        hCode ^= fixedDcc->hashCode();

    if (floatIvl.get())
        hCode ^= floatIvl->hashCode();

    if (floatDcc.get())
        hCode ^= floatDcc->hashCode();

    if (basisIvl.get())
        hCode ^= basisIvl->hashCode();

    if (basisDcc.get())
        hCode ^= basisDcc->hashCode();

    return hCode;
}


bool SwapBenchmark::equals(const ZeroCurveBenchmark& other) const
{
    if (!ZeroCurveBenchmark::equals(other))
        return false;

    const SwapBenchmark& rhs = static_cast<const SwapBenchmark&>(other);

    if (adjustment != rhs.adjustment)
        return false;

    if (price != rhs.price)
        return false;

    if (!equalsClass<DayCountConvention>(fixedDcc, rhs.fixedDcc))
        return false;

    if (!equalsClass<DayCountConvention>(floatDcc, rhs.floatDcc))
        return false;

    if (!equalsClass<DayCountConvention>(basisDcc, rhs.basisDcc))
        return false;

    if (!equalsValue<MaturityPeriod>(fixedIvl, rhs.fixedIvl))
        return false;

    if (!equalsValue<MaturityPeriod>(floatIvl, rhs.floatIvl))
        return false;

    if (!equalsValue<MaturityPeriod>(basisIvl, rhs.basisIvl))
        return false;

    return true;
}


bool SwapBenchmark::hasAdjustment() const
{
    return !Maths::isZero(adjustment);
}


double SwapBenchmark::getAdjustment() const
{
    return adjustment;
}


double SwapBenchmark::getPrice() const
{
    return price;
}


const MaturityPeriod* SwapBenchmark::getFixedIvl() const
{
    return fixedIvl.get();
}


const MaturityPeriod* SwapBenchmark::getFloatIvl() const
{
    return floatIvl.get();
}


const MaturityPeriod* SwapBenchmark::getBasisIvl() const
{
    return basisIvl.get();
}


const DayCountConvention* SwapBenchmark::getFixedDcc() const
{
    return fixedDcc.get();
}


const DayCountConvention* SwapBenchmark::getFloatDcc() const
{
    return floatDcc.get();
}


const DayCountConvention* SwapBenchmark::getBasisDcc() const
{
    return basisDcc.get();
}


/*
 * Reflection support
 */

/** Invoked when Class is 'loaded' */
void ZeroCurveBenchmark::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ZeroCurveBenchmark, clazz);
    SUPERCLASS(CObject);
    FIELD       (start,         "start");
    FIELD       (end,           "end");
    FIELD(rate,          "rate");
    FIELD(includeFlag,   "include flag");
    FIELD_MAKE_OPTIONAL(start);
    FIELD_MAKE_OPTIONAL(includeFlag);
}


CClassConstSP const ZeroCurveBenchmark::TYPE = 
    CClass::registerClassLoadMethod("ZeroCurveBenchmark", typeid(ZeroCurveBenchmark), load);


DEFINE_TEMPLATE_TYPE(ZeroCurveBenchmarkArray);


/** Invoked when Class is 'loaded' */
static IObject* defaultMoneyMarketBenchmark()
{
    return new MoneyMarketBenchmark();
}


void MoneyMarketBenchmark::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MoneyMarketBenchmark, clazz);
    SUPERCLASS(ZeroCurveBenchmark);
    EMPTY_SHELL_METHOD(defaultMoneyMarketBenchmark);
    FIELD(dcc, "Day count convention");
    FIELD_MAKE_OPTIONAL(dcc);
}


CClassConstSP const MoneyMarketBenchmark::TYPE = 
    CClass::registerClassLoadMethod("MoneyMarketBenchmark", typeid(MoneyMarketBenchmark), load);


/** Invoked when Class is 'loaded' */
static IObject* defaultFuturesBenchmark()
{
    return new FuturesBenchmark();
}


void FuturesBenchmark::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FuturesBenchmark, clazz);
    SUPERCLASS(ZeroCurveBenchmark);
    EMPTY_SHELL_METHOD(defaultFuturesBenchmark);
    FIELD(type,          "futures type");
    FIELD(adjustment,    "adjustment");
    FIELD       (dcc,           "Day count convention");
    FIELD_MAKE_OPTIONAL(type);
    FIELD_MAKE_OPTIONAL(adjustment);
    FIELD_MAKE_OPTIONAL(dcc);
}


CClassConstSP const FuturesBenchmark::TYPE = 
    CClass::registerClassLoadMethod("FuturesBenchmark", typeid(FuturesBenchmark), load);


/** Invoked when Class is 'loaded' */
static IObject* defaultTurnBenchmark()
{
    return new TurnBenchmark();
}


void TurnBenchmark::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(TurnBenchmark, clazz);
    SUPERCLASS(ZeroCurveBenchmark);
    EMPTY_SHELL_METHOD(defaultTurnBenchmark);
    FIELD(type, "adjustment type - A(bsolute), R(elative) or I(mplied)");
    FIELD       (dcc,  "Day count convention");
    FIELD_MAKE_OPTIONAL(type);
    FIELD_MAKE_OPTIONAL(dcc);
}


CClassConstSP const TurnBenchmark::TYPE = 
    CClass::registerClassLoadMethod("TurnBenchmark", typeid(TurnBenchmark), load);


/** Invoked when Class is 'loaded' */
static IObject* defaultSwapBenchmark()
{
    return new SwapBenchmark();
}


void SwapBenchmark::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SwapBenchmark, clazz);
    SUPERCLASS(ZeroCurveBenchmark);
    EMPTY_SHELL_METHOD(defaultSwapBenchmark);
    FIELD(price,      "price");
    FIELD       (fixedIvl,   "fixed leg interval");
    FIELD       (fixedDcc,   "fixed leg day count convention");
    FIELD       (floatIvl,   "floating leg interval");
    FIELD       (floatDcc,   "floating leg day count convention");
    FIELD       (basisIvl,   "basis leg interval");
    FIELD       (basisDcc,   "basis leg day count convention");
    FIELD(adjustment, "adjustment");
    FIELD_MAKE_OPTIONAL(price);
    FIELD_MAKE_OPTIONAL(fixedIvl);
    FIELD_MAKE_OPTIONAL(fixedDcc);
    FIELD_MAKE_OPTIONAL(floatIvl);
    FIELD_MAKE_OPTIONAL(floatDcc);
    FIELD_MAKE_OPTIONAL(basisIvl);
    FIELD_MAKE_OPTIONAL(basisDcc);
    FIELD_MAKE_OPTIONAL(adjustment);
}


CClassConstSP const SwapBenchmark::TYPE = 
    CClass::registerClassLoadMethod("SwapBenchmark", typeid(SwapBenchmark), load);


/** Addin for building handles to an array of benchmarks */
class BenchmarkArrayAddin: public CObject
{
    static CClassConstSP const TYPE;

    static IObject* defaultBenchmarkArrayAddin()
    {
        return new BenchmarkArrayAddin();
    }

    ExpiryArraySP startDates;
    ExpiryArray   expiries;
    StringArray   instruments;
    IntArraySP    includeFlags;
    DoubleArray   rates;
    DoubleArraySP prices;
    DoubleArraySP adjustments;

    static IObjectSP createBenchmarkArray(BenchmarkArrayAddin* params)
    {
        return ZeroCurveBenchmark::toArray(
            params->startDates.get(),
            params->expiries,
            params->instruments,
            params->includeFlags.get(),
            params->rates,
            params->prices.get(),
            params->adjustments.get());
    }

    /** for reflection */
    BenchmarkArrayAddin() : CObject(TYPE)
    {
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        REGISTER(BenchmarkArrayAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBenchmarkArrayAddin);

        FIELD       (startDates,   "start dates");
        FIELD(expiries,     "expiries");
        FIELD(instruments,  "instruments");
        FIELD       (includeFlags, "include flags");
        FIELD(rates,        "rates");
        FIELD       (prices,       "prices");
        FIELD       (adjustments,  "adjustments");

        FIELD_MAKE_OPTIONAL(startDates);
        FIELD_MAKE_OPTIONAL(includeFlags);
        FIELD_MAKE_OPTIONAL(prices);
        FIELD_MAKE_OPTIONAL(adjustments);

        Addin::registerClassObjectMethod(
            "BENCHMARKS",
            Addin::MARKET,
            "Constructs an array of benchmarks",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*) createBenchmarkArray);
    }
};


CClassConstSP const BenchmarkArrayAddin::TYPE = CClass::registerClassLoadMethod(
    "BenchmarkArrayAddin", typeid(BenchmarkArrayAddin), load);


DRLIB_END_NAMESPACE
