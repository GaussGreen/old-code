//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : ZeroCurveBenchmark.hpp
//
//   Description : Input benchmarks for curve bootstrapping
//
//   Author      : Richard Appleton
//
//   Date        : 15th November 2005
//
//----------------------------------------------------------------------------

#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/Holiday.hpp"


DRLIB_BEGIN_NAMESPACE
class BadDayConvention;
class DayCountConvention;


class ZeroCurveBenchmark;
typedef smartConstPtr<ZeroCurveBenchmark>               ZeroCurveBenchmarkConstSP;
typedef smartPtr<ZeroCurveBenchmark>                    ZeroCurveBenchmarkSP;
typedef array<ZeroCurveBenchmarkSP, ZeroCurveBenchmark> ZeroCurveBenchmarkArray;
typedef smartPtr<ZeroCurveBenchmarkArray>               ZeroCurveBenchmarkArraySP;


class MARKET_DLL ZeroCurveBenchmark : public CObject
{
public:
    /**
     * Converts arrays of values into an array of benchmarks.
     */
    static ZeroCurveBenchmarkArraySP toArray(
        const ExpiryArray* starts,
        const ExpiryArray& expiries,
        const StringArray& instruments,
        const IntArray*    includeFlags,
        const DoubleArray& rates,
        const DoubleArray* prices,
        const DoubleArray* adjustments,
        const DateTime*    refDate = NULL,
        const Holiday*     holidays = NULL);

    /**
     * Add benchmark to benchmark array.
     */
    static void append(
        ZeroCurveBenchmarkArray& benchmarks, 
        int                      i,             // index of arrays to use for new benchmark
        const ExpiryArray*       starts,
        const ExpiryArray&       expiries,
        const StringArray&       instruments,
        const IntArray*          includeFlags,
        const DoubleArray&       rates,
        const DoubleArray*       prices,
        const DoubleArray*       adjustments,
        const DateTime*          refDate = NULL,
        const Holiday*           holidays = NULL);

    static bool equals(const ZeroCurveBenchmarkArray& lhs, const ZeroCurveBenchmarkArray& rhs);

    /**
     * Returns index of benchmark in array with specified expiry.
     */
    static int search(const ZeroCurveBenchmarkArray& benchmarks, const Expiry& expiry);

    /**
     * Sort array into expiry order with tenors resolved from value date.
     */
    static void sort(ZeroCurveBenchmarkArray& benchmarks, const DateTime& valueDate);

    /**
     * Calculate hash code used in zero curve caching.
     */
    static int hashCode(const ZeroCurveBenchmarkArray& benchmarks);

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    virtual ~ZeroCurveBenchmark();

    virtual IObject* clone() const = 0;
    virtual int      hashCode() const;
    virtual bool     equals(const ZeroCurveBenchmark& other) const;

    const Expiry*    getStart() const;  // never NULL for futures
    const Expiry*    getEnd() const;    // may be NULL for futures
    virtual DateTime getBenchmarkDate(const DateTime& date) const;
    double           getRate() const;
    int              getIncludeFlag() const;
    virtual bool     hasAdjustment() const;
    virtual double   getAdjustment() const;

    void setRate(double rate);

    // convenience type tests
    virtual bool isCash() const;
    virtual bool isFuture() const;
    virtual bool isTurn() const;
    virtual bool isSwap() const;

protected:
    ZeroCurveBenchmark(const CClassConstSP& klass = TYPE);

    ZeroCurveBenchmark(
        const Expiry*        start, 
        const Expiry*        end,
        double               rate, 
        int                  includeFlag,
        const CClassConstSP& klass = TYPE);

    ExpiryConstSP start;
    ExpiryConstSP end;
    double        rate;
    int           includeFlag;
};


class MARKET_DLL MoneyMarketBenchmark : public ZeroCurveBenchmark
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    MoneyMarketBenchmark();

    MoneyMarketBenchmark(
        const Expiry*             start, 
        const Expiry&             end,
        double                    rate, 
        int                       includeFlag,
        const DayCountConvention* dcc);

    bool isCash() const;

    IObject*                  clone() const;
    int                       hashCode() const;
    bool                      equals(const ZeroCurveBenchmark& other) const;
    const DayCountConvention* getDcc() const;

private:
    DayCountConventionConstSP dcc;
};


class MARKET_DLL FuturesBenchmark : public ZeroCurveBenchmark
{
public:
    static const double ZERO_PERCENT;

    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    FuturesBenchmark();

    FuturesBenchmark(
        const Expiry&             start, 
        const Expiry*             end,
        double                    quote, 
        bool                      quotedAsPrice,
        int                       includeFlag,
        const DayCountConvention* dcc,
        const double*             adjustment,
        const string&             type);

    bool isCash() const;
    bool isFuture() const;

    IObject*                  clone() const;
    int                       hashCode() const;
    bool                      equals(const ZeroCurveBenchmark& other) const;
    bool                      hasAdjustment() const;
    double                    getAdjustment() const;
    const string&             getType() const;
    const DayCountConvention* getDcc() const;
    DateTime                  getBenchmarkDate(const DateTime& date) const;

protected:
    void copyTo(FuturesBenchmark* copy) const;

private:
    //    int    basis;     // eg. CompoundBasis::SIMPLE
    string                    type; // futures type "1M", "3M", "AUD", "USD_TBILL"
    double                    adjustment;
    DayCountConventionConstSP dcc;
};

class MARKET_DLL BrazilFuturesBenchmark : public FuturesBenchmark
{
public:
    BrazilFuturesBenchmark();

    BrazilFuturesBenchmark(
        const Expiry&   expiry,
        double          quote,
        bool            quotedAsPrice,
        int             includeFlag,
        const DateTime& refDate,
        const Holiday&  holidays);

    static ExpiryConstSP adjustExpiry(
        const Expiry&   expiry,
        const DateTime& refDate,
        const Holiday&  holidays);

    double getPrice(
        const DateTime& refDate,
        const Holiday&  holidays) const;

    IObject* clone() const;
};

class MARKET_DLL TurnBenchmark : public ZeroCurveBenchmark
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    // codes for turn type
    static const string ABSOLUTE;
    static const string IMPLIED;
    static const string RELATIVE;
    static const string UTURN;

    TurnBenchmark();

    TurnBenchmark(
        const string&             type,
        const Expiry*             start, 
        const Expiry&             end,
        double                    rate,
        int                       includeFlag,
        const DayCountConvention* dcc);

    bool isCash() const;
    bool isTurn() const;

    IObject*                  clone() const;
    int                       hashCode() const;
    bool                      equals(const ZeroCurveBenchmark& other) const;
    bool                      isAbsolute() const;
    bool                      isRelative() const;
    bool                      isImplied() const;
    const DayCountConvention* getDcc() const;

private:
    string                    type;
    DayCountConventionConstSP dcc;
};


class MARKET_DLL SwapBenchmark : public ZeroCurveBenchmark
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    SwapBenchmark();
    
    SwapBenchmark(
        const Expiry*             start, 
        const Expiry&             end,
        double                    rate, 
        const double*             adjustment,
        double                    price,
        int                       includeFlag,
        const MaturityPeriod*     fixedIvl,
        const DayCountConvention* fixedDcc,
        const MaturityPeriod*     floatIvl,
        const DayCountConvention* floatDcc,
        const MaturityPeriod*     basisIvl,
        const DayCountConvention* basisDcc);

    bool isSwap() const;

    IObject*                  clone() const;
    int                       hashCode() const;
    bool                      equals(const ZeroCurveBenchmark& other) const;
    bool                      hasAdjustment() const;
    double                    getAdjustment() const;
    double                    getPrice() const;
    const MaturityPeriod*     getFixedIvl() const;
    const MaturityPeriod*     getFloatIvl() const;
    const MaturityPeriod*     getBasisIvl() const;
    const DayCountConvention* getFixedDcc() const;
    const DayCountConvention* getFloatDcc() const;
    const DayCountConvention* getBasisDcc() const;

private:
    double                    price;
    MaturityPeriodConstSP     fixedIvl;
    DayCountConventionConstSP fixedDcc;
    MaturityPeriodConstSP     floatIvl;
    DayCountConventionConstSP floatDcc;
    MaturityPeriodConstSP     basisIvl;
    DayCountConventionConstSP basisDcc;
    double                    adjustment;
};



DRLIB_END_NAMESPACE

#endif // BENCHMARK_HPP
