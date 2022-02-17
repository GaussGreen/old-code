//----------------------------------------------------------------------------
//
//   Group       : QR Core Analytics Team
//
//   Filename    : LMParabCurve.hpp
//
//   Description : ALIB style local market parabolic curve
//
//   Author      : Andrew Greene
//
//   Date        : 07 April 2006
//
//----------------------------------------------------------------------------

#ifndef LM_PARAB_CURVE_HPP
#define LM_PARAB_CURVE_HPP

#include "edginc/CashFlow.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/ZC3ZeroCurve.hpp"

DRLIB_BEGIN_NAMESPACE

class DayCountConvention;
class LMParabCurve;

typedef refCountPtr<LMParabCurve> LMParabCurveSP;
typedef refCountPtr<const LMParabCurve> LMParabCurveConstSP;

// Holds a zero curve which uses parabolic coefficients to
// allow different interpolation methods on each benchmark interval.
class MARKET_DLL LMParabCurve
{
public:

    // Represents a Forward Rate Agreement as
    // start date, maturity date and negative
    // logarithm of the forward discount factor.
    struct MARKET_DLL FwdRateIvl
    {
        DateTime startDate;
        DateTime matDate;
        double   rate;
    };

    // Defines a year end effect interval.
    // This type of information is often used
    // to adjust a zero curve.
    struct MARKET_DLL YearEndEffect
    {
        FwdRateIvl ivl;
        bool       isRateSpread;
    };

// -----------------------------------------------------------------------------
// From yldgen.h

    // Constructs a curve given a FwdRateIvl list
    // and a YearEndEffect list.
    // (Based on Alib function LmGenTParabCurveFromFras)
    LMParabCurve
    (
        const DateTime&              baseDate,         // (I) base date for curve
        const vector<FwdRateIvl>*    frl,              // (I)
        const vector<YearEndEffect>* yel,              // (I)
        HolidayConstSP               noAccrue,         // (I) dates with no accrual
        const vector<DateTime>*      flatRateDates,    // (I) the rate will be constant in between flatRateDates,
                                                       //     last date must be <= flatParabBdryDate
        const DateTime&              flatParabBdryDate // (I) date at which interpolation transitions
                                                       //     from flat to parabolic
    );

// -----------------------------------------------------------------------------
// From lmzcutil.h

    // Constructs a curve associated to Money Market Instruments
    //   (1) Money Market Points
    //   (2) Fra Points
    //   (3) Year End Effects
    //   (4) Relative Money Market Points, ie, priced relative to the curve
    //       built out of the previous instruments
    // It has a two stage generation to allow for relative benchmarks.
    // (Based on Alib function LmMMParabCurve)
    LMParabCurve
    (
        const DateTime&           baseDate,          // (I) base date for zero curve.  Discounts from this date
        const vector<DateTime>*   mmMatDates,        // (I) vector of MM maturity dates
        const vector<double>*     mmRates,           // (I) vector of MM rates
        const vector<DateTime>*   fraStartDates,     // (I) vector of Fra start dates
        const vector<DateTime>*   fraMatDates,       // (I) vector of Fra maturity dates
        const vector<double>*     fraRates,          // (I) vector of Fra rates
        const vector<DateTime>*   yeStartDates,      // (I) vector of Year End adjustment start dates
        const vector<DateTime>*   yeEndDates,        // (I) vector of Year End adjustment end dates
        const vector<double>*     yeRates,           // (I) vector of Year End adjustment rates
        const vector<bool>*       yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                     //     spread, otherwise in absolute level
        const vector<DateTime>*   relMMDates,        // (I) vector of relative MM maturity dates
        const vector<double>*     relMMSpreads,      // (I) vector of spreads for relative MM points
        const DayCountConvention& mmDayCount,        // (I) Common Day Count Convention for MM, Fras and YE
        int                       rateBasis,         // (I) Common Rate Basis for MM, Fras and YE
        const vector<DateTime>*   flatRateDates,     // (I) the rate will be constant in between flatRateDates,
                                                     //     last date must be <= flatParabBdryDate
        const DateTime&           flatParabBdryDate, // (I) date at which interpolation transitions
                                                     //     from flat to parabolic
        HolidayConstSP            noAccrue,          // (I) no accrual dates 
        const string&             interpType         // (I) smooth forwards or linear
    );

// -----------------------------------------------------------------------------
// From lmzcgen.h

    // Constructs a Zero Curve associated to Money Market Instruments and CashFlows
    //
    //     (1) Money Market Points
    //     (2) Fra Points
    //     (3) Year End Effects
    //     (4) Relative Money Market Points, ie, priced relative to the curve
    //         built out of the previous instruments
    //     (5) CFL Points
    //
    // It has a two stage generation to allow for relative MM benchmarks.
    // (Based on Alib function GtoLmCFLCurve)
    LMParabCurve
    (
        const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
        const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
        const vector<double>*        mmRates,           // (I) vector of MM rates
        const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
        const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
        const vector<double>*        fraRates,          // (I) vector of Fra rates
        const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
        const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
        const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
        const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                        //     spread, otherwise in absolute level
        const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
        const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
        const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
        int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
        const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
        const vector<CashFlowArray>* cfl,               // (I) vector of cashflow arrays
        const vector<double>*        cflPrice,          // (I) vector of dirty prices for cashflow arrays
        const DateTime&              flatParabBdryDate, // (I) date at which interpolation transitions
                                                        //     from flat to parabolic
        HolidayConstSP               noAccrue,          // (I) no accrual dates 
        const string&                interpType         // (I) smooth forwards or linear
    );

    // Creates a Zero Curve associated to Money Market Instruments and Bonds
    //     (1) Money Market Points
    //     (2) Fra Points
    //     (3) Year End Effects
    //     (4) Relative Money Market Points, ie, priced relative to the curve
    //         built out of the previous instruments
    //     (5) Bond Points
    // It has a two stage generation to allow for relative MM benchmarks.
    // (Based on Alib function GtoLmBondCurve)
    static LMParabCurveSP genFromBonds
    (
        const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
        const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
        const vector<double>*        mmRates,           // (I) vector of MM rates
        const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
        const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
        const vector<double>*        fraRates,          // (I) vector of Fra rates
        const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
        const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
        const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
        const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                        //     spread, otherwise in absolute level
        const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
        const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
        const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
        int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
        const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
        const vector<DateTime>*      bondMat,           // (I) bond maturity dates
        const vector<double>*        bondRate,          // (I) bond rates
        int                          bondFreq,          // (I) bond number of compounding periods per year
        const DayCountConvention&    bondDayCount,      // (I) bond day count convention
        const BadDayConvention&      accBadDayConv,     // (I) bond bad day convention for interest rate accrual
        const BadDayConvention&      payBadDayConv,     // (I) bond bad day convention for payment dates
        const vector<double>*        bondPrice,         // (I) bond dirty prices
        const DateTime&              flatParabBdryDate, // (I) date at which interpolation transitions
                                                        //     from flat to parabolic
        HolidayConstSP               noAccrue,          // (I) no accrual dates 
        const string&                interpType         // (I) smooth forwards or linear
    );

    // Creates a Zero Curve associated to Money Market Instruments and Swaps
    //     (1) Money Market Points
    //     (2) Fra Points
    //     (3) Year End Effects
    //     (4) Relative Money Market Points, ie, priced relative to the curve
    //         built out of the previous instruments
    //     (5) Swap Points
    // It has a two stage generation to allow for relative MM benchmarks.
    // (Based on Alib function GtoLmSwapCurve)
    static LMParabCurveSP genFromSwaps
    (
        const DateTime&              baseDate,          // (I) base date for zero curve.  Discounts from this date
        const vector<DateTime>*      mmMatDates,        // (I) vector of MM maturity dates
        const vector<double>*        mmRates,           // (I) vector of MM rates
        const vector<DateTime>*      fraStartDates,     // (I) vector of Fra start dates
        const vector<DateTime>*      fraMatDates,       // (I) vector of Fra maturity dates
        const vector<double>*        fraRates,          // (I) vector of Fra rates
        const vector<DateTime>*      yeStartDates,      // (I) vector of Year End adjustment start dates
        const vector<DateTime>*      yeEndDates,        // (I) vector of Year End adjustment end dates
        const vector<double>*        yeRates,           // (I) vector of Year End adjustment rates
        const vector<bool>*          yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                        //     spread, otherwise in absolute level
        const vector<DateTime>*      relMMDates,        // (I) vector of relative MM maturity dates
        const vector<double>*        relMMSpreads,      // (I) vector of spreads for relative MM points
        const DayCountConvention&    mmDayCount,        // (I) common day count convention for MM, Fras and YE
        int                          rateBasis,         // (I) common rate basis for MM, Fras and YE
        const vector<DateTime>*      flatRateDates,     // (I) the rate will be constant in between flatRateDates,
        const vector<DateTime>*      swapMat,           // (I) bond maturity dates
        const vector<double>*        swapRate,          // (I) bond rates
        int                          swapFreq,          // (I) bond number of compounding periods per year
        const DayCountConvention&    swapDayCount,      // (I) bond day count convention
        const BadDayConvention&      accBadDayConv,     // (I) bond bad day convention for interest rate accrual
        const BadDayConvention&      payBadDayConv,     // (I) bond bad day convention for payment dates
        const DateTime&              flatParabBdryDate, // (I) date at which interpolation transitions
                                                        //     from flat to parabolic
        HolidayConstSP               noAccrue,          // (I) no accrual dates 
        const string&                interpType         // (I) smooth forwards or linear
    );

// -----------------------------------------------------------------------------
// From lmzcutil.h

    // Given definite dates, discs, parabolic coefficients, calculates the
    // disc from start date to mat date.
    // (Based on Alib function LmParabCurveDisc)
    double disc
    (
        const DateTime& startDate, // (I) start date for discount factor
        const DateTime& endDate    // (I) end date for discount factor
    )
    const;

    // Calculates the fwd rate with the given day count and basis
    // (Based on Alib function LmMMRate)
    double mmRate
    (
        const DateTime&           startDate,    // (I) start date
        const DateTime&           endDate,      // (I) end date
        const DayCountConvention& dayCountConv, // (I) day count convention
        int                       basis         // (I) rate basis
    )
    const;

// -----------------------------------------------------------------------------
// From tpcurve.h

    // Returns an equivalent ZC3ZeroCurve
    // (based on Alib function LmConvTPCToSZC)
    ZC3ZeroCurveSP toZeroCurve();

// -----------------------------------------------------------------------------

private:

    // Holds a date, rate and parabolic coefficients
    // for interpolation on the period from previous
    // benchmark date to fDate
    struct MARKET_DLL RatePt
    {
        DateTime date;
        double disc;

        // coeff[0] is constant parabolic coefficient
        // coeff[1] is linear   parabolic coefficient
        // coeff[2] is leading  parabolic coefficient
        double coeff[3];
    };

// -----------------------------------------------------------------------------
// From bondgen.h

    // Extends an LMParabCurve constructed as an MM stub curve,
    // using a list of cashflows with their corresponding price
    // (Based on Alib function LmGenTParabCurveFromCF)
    void addCF
    (
        const vector<CashFlowArray>& cf,               // (I) vector of cashflow arrays to be added
        const vector<double>&        cfPrice,          // (I) array of prices
        const DateTime&              flatParabBdryDate // (I) date at which interpolation transitions
                                                       //     from flat to parabolic
    );

    // Adds the next fra period to the current built curve
    // using flat forwards.
    // (Based on Alib function LmAddNextCFFlat)
    void addNextCFFlat
    (
        const CashFlowArray& cashFlow, // (I) cashflow to be added
        double               cfPrice   // (I) PV of the cashflow to curve base date
    );

    // Introduces the next period of the curve as
    // a parabolic segment.
    // (Based on Alib function LmAddNextCFParab)
    void addNextCFParab
    (
        const CashFlowArray& currCF,      // (I) cashflows to be added
        double               currPrice,   // (I) price of the current cashflows
        const CashFlowArray& nextCF,      // (I) next cashflows to be added
        double               nextPrice,   // (I) price of the next cashflows
        double&              lastInstRate // (I/O) instantaneous forward rate at end of period
    );

    // (Based on Alib function LmAddFirstCFLinear)
    void addFirstCFLinear
    (
        const CashFlowArray& currCF,      // (I) cashflows to be added
        double               currPrice,   // (I) price of the current cashflows
        const CashFlowArray& nextCF,      // (I) next cashflows to be added
        double               nextPrice,   // (I) price of the next cashflows
        double&              lastInstRate // (O) instantaneous forward rate at end of period
    );

    // Introduces the last period of the curve as
    // a parabolic segment.
    // (Based on Alib function LmAddLastCFParab)
    void addLastCFParab
    (
         const CashFlowArray& currCF,      // (I) cashflows to be added
         double               currPrice,   // (I) price of the current cashflows
         double               lastInstRate // (I) instantaneous forward rate at end of period
    );

    // Takes the last two segments of a curve and computes
    // the inverse duration weighted average for parabolic
    // interpolation.
    // (Based on Alib function LmParabInstBdryRate)
    double instBdryRate
    (
        const CashFlowArray& currCF, // (I) cashflows to be added
        const CashFlowArray& nextCF  // (I) next cashflow to be added
    )
    const;

    // Takes the last two segments of a curve and computes
    // the inverse duration weighted average for parabolic
    // interpolation.
    // (Based on Alib function LmParabInstBdryRateLast)
    double instBdryRateLast
    (
        const CashFlowArray& currCF,      // (I) cashflows to be added
        double               lastInstRate // (I) instantaneous forward rate at start of period
    )
    const;

// -----------------------------------------------------------------------------
// From linzcgen.h

    // Given a start Date, and discount to start of interval,
    // compute its linear zeros approximation coefficients.
    // Uses definition of Zero Rate as (1+rate)^(accrueDays/timeDenom).
    // (Based on Alib function LmLinZerosCoeffsTParabRatePt)
    static void linZerosCoeffsRatePt
    (
        RatePt&         rp,        // (I/O) parabolic rate interval
        const DateTime& baseDate,  // (I) needed to compute the zero rate
        const DateTime& startDate, // (I) start date of the interval
        double          startDisc, // (I) discount factor to start of period
        double          timeDenom, // (I) time denominator for zero coupon rates
        HolidayConstSP  noAccrue   // (I) dates with no accrual
    );

    // Constructs a curve given a FwdRateIvl list
    // and a YearEndEffect list.
    // (Based on Alib function LmLinZerosTParabCurveFromFras)
    void linZerosGenFromFras
    (
        const DateTime&              bd,  // (I) base date for curve
        const vector<FwdRateIvl>&    frl, // (I)
        const vector<YearEndEffect>& yel, // (I)
        HolidayConstSP               na   // (I) dates with no accrual
    );

    // Adds the next fra period to the current built curve
    // using approximation to linear zeros.
    void addNextFraLinearZeros
    (
        DateTime&                 lastDate, // (I/O) last date entered into curve
        double&                   lastDisc, // (I/O) discount up to lastDate
        const vector<FwdRateIvl>& fraList,  // (I) fras to be included
        size_t&                   idxFra,   // (I/O) index of Fra being used currently
        double                    timeDenom // (I) time denominator for zero coupon rates
    );

    // Adds the next CF to the current built curve
    // using linear zeros approximation
    // (Based on Alib function LmAddNextCFLinearZeros)
    void addNextCFLinearZeros
    (
        const CashFlowArray& cashFlow, // (I) cashflow to be added
        double               cfPrice,  // (I) PV of the cashflow to curve base date        
        double               timeDenom // (I) time denominator for zero coupon rates
    );

    // Extends an LMParabCurve constructed as an MM stub curve,
    // using a list of cashflows with their corresponding price
    // (Based on Alib function LmLinZerosTParabCurveFromCF)
    void linZerosAddCF
    (
        const vector<CashFlowArray>& cf,     // (I) vector of cashflow arrays to be added
        const vector<double>&        cfPrice // (I) array of prices
    );

// -----------------------------------------------------------------------------
// From lmzcutil.h

    // Constructs a curve associated to Money Market Instruments
    //   (1) Money Market Points
    //   (2) Fra Points
    //   (3) Year End Effects
    //   (4) Relative Money Market Points, ie, priced relative to the curve
    //       built out of the previous instruments
    // It has a two stage generation to allow for relative benchmarks.
    // (Based on Alib function LmMMParabCurve)
    void genFromMM
    (
        const DateTime&           bd,                // (I) base date for zero curve.  Discounts from this date
        const vector<DateTime>&   mmMatDates,        // (I) vector of MM maturity dates
        const vector<double>&     mmRates,           // (I) vector of MM rates
        const vector<DateTime>&   fraStartDates,     // (I) vector of Fra start dates
        const vector<DateTime>&   fraMatDates,       // (I) vector of Fra maturity dates
        const vector<double>&     fraRates,          // (I) vector of Fra rates
        const vector<DateTime>&   yeStartDates,      // (I) vector of Year End adjustment start dates
        const vector<DateTime>&   yeEndDates,        // (I) vector of Year End adjustment end dates
        const vector<double>&     yeRates,           // (I) vector of Year End adjustment rates
        const vector<bool>&       yeIsSpread,        // (I) if TRUE the Year End adjustment is used as a
                                                     //     spread, otherwise in absolute level
        const vector<DateTime>&   relMMDates,        // (I) vector of relative MM maturity dates
        const vector<double>&     relMMSpreads,      // (I) vector of spreads for relative MM points
        const DayCountConvention& mmDayCount,        // (I) Common Day Count Convention for MM, Fras and YE
        int                       rateBasis,         // (I) Common Rate Basis for MM, Fras and YE
        const vector<DateTime>&   flatRateDates,     // (I) the rate will be constant in between flatRateDates,
                                                     //     last date must be <= flatParabBdryDate
        const DateTime&           flatParabBdryDate, // (I) date at which interpolation transitions
                                                     //     from flat to parabolic
        HolidayConstSP            na,                // (I) no accrual dates 
        const string&             interpType         // (I) smooth forwards or linear
    );

    // Error checking on inputs to mm
    // before allowing analytics processing.
    //
    // Checks for:
    //   (1) MMDates should be in increasing order
    //   (2) FraStartDates[i] < FraMatDates[i]
    //   (3) YEStartDates[i] < YEMatDates[i]
    //   (4) relMMDates should be in increasing order,
    //       and smaller than the biggest of MMDates or FraMatDates
    // (Based on Alib function LmValidateMMData)
    static void validateMMData
    (
        const vector<DateTime>& mmMatDates,    // (I) vector of MM maturity dates
        const vector<double>&   mmRates,       // (I) vector of MM rates
        const vector<DateTime>& fraStartDates, // (I) vector of Fra start dates
        const vector<DateTime>& fraMatDates,   // (I) vector of Fra maturity dates
        const vector<double>&   fraRates,      // (I) vector of Fra rates
        const vector<DateTime>& yeStartDates,  // (I) vector of Year End adjustment start dates
        const vector<DateTime>& yeEndDates,    // (I) vector of Year End adjustment end dates
        const vector<double>&   yeRates,       // (I) vector of Year End adjustment rates
        const vector<bool>&     yeIsSpread,    // (I) if TRUE the Year End adjustment is used as a
                                               //     spread, otherwise in absolute level
        const vector<DateTime>& relMMDates,    // (I) vector of relative MM maturity dates
        const vector<double>&   relMMSpreads,  // (I) vector of spreads for relative MM points
        const string&           interpType     // (I) Smooth forwards or linear
    );

// -----------------------------------------------------------------------------
// From tfraivl.h

    // Builds a vector of FwdRateIvls from MM rates
    static void appendFwdRateIvlListMM
    (
        const vector<DateTime>&   startDates, // (I)
        const vector<DateTime>&   matDates,   // (I)
        const vector<double>&     rates,      // (I)
        const DayCountConvention& dayCount,   // (I)
        int                       rateBasis,  // (I)
        vector<FwdRateIvl>&       list        // (I/O) List to be added to
    );

    static void appendYearEndEffectListMM
    (
        const vector<DateTime>&   startDates,   // (I)
        const vector<DateTime>&   matDates,     // (I)
        const vector<double>&     rates,        // (I)
        const vector<bool>&       isRateSpread, // (I)
        const DayCountConvention& dayCount,     // (I)
        int                       rateBasis,    // (I)
        vector<YearEndEffect>&    list          // (I/O) List to be added to
    );

    // True if within each interval startDate[i] < matDate[i]
    // (Based on Alib function LmIsValidTFwdRateIvlList)
    static bool isValidFwdRateIvlList
    (
        const vector<FwdRateIvl>& list // (I) List to be checked
    );

    // True if the FRA list is stripped.
    // That is, matDate[i] <= startDate[i+1]
    // (Based on Alib function LmIsStripTFwdRateIvlList)
    static bool isStripFwdRateIvlList
    (
        const vector<FwdRateIvl>& list // (I) FwdRateIvl list to be checked
    );

    // True if YearEndEffect dates are in
    // ascending order.  matDate[i] <=startDate[i+1]
    // and within each YearEndEffect startDate precedes
    // matDate.
    // (Based on Alib function LmIsOrderedTYearEndEffectList)
    static bool isOrderedYearEndEffectList
    (
        const vector<YearEndEffect>& list // (I) YearEndEffect list to be checked
    );

    // Sorts all FRAs.  Returns false
    // if two FRAs share same startDate and matdate.
    // ie if the TFwdRateIvlList is not monotonic.
    // (Based on Alib function LmSortTFwdRateIvlList)
    static bool sortFwdRateIvlList
    (
        vector<FwdRateIvl>& listFras // (I/O) FRA dates to sort
    );

    // Predicate for FRA date comparison
    static bool compareFra
    (
        const LMParabCurve::FwdRateIvl& fra1, // (I) fra1
        const LMParabCurve::FwdRateIvl& fra2  // (I) fra2
    );

    // Creates a new list of Fras which is sorted by increasing
    // maturity date first and 
    // by increasing start date for those with same matDate
    // Converts input clusters of FRAs which have
    // same matDate and into an equal number of 
    // consecutive FRAs over the same period spanned
    // by the cluster.
    // (based on Alib function LmSortAndNormalizeFRAList)
    static bool sortAndNormalizeFRAList
    (
        const vector<FwdRateIvl>& list,   // (I)   list of Fras
        vector<FwdRateIvl>&       newList // (I/O) result
    );

    // Locates all FRAs which are clustered at
    // the endDate of indexFRA FRA and normalizes
    // the input list of Fras in place
    // (based on Alib function lmNormalize1FRADate)
    static bool normalize1FRADate
    (
        size_t indexFRA,         // (I) index of FRA Date to check
        vector<FwdRateIvl>& list // (I/O) FRA dates to normalize and sort
    );

    // Combines a no accrue file and the dates covered
    // by an array of year end effects into a
    // merged holiday file which is stored in the cache
    // called by name in holidayYE.
    // (based on Alib function LmYearEndAndNoAccrueCombine)
    static void yearEndAndNoAccrueCombine
    (
        const vector<YearEndEffect>& listYE, // (I)
        HolidayConstSP noAccrue,             // (I)   current no accrue
        HolidaySP& noAccrueYE                // (I/O) new no accrue
    );

    // Creates a new list of Fras which differs from
    // the input list of Fras by:
    // (1) being sorted and normalized
    // (2) having all year end effect rates subtracted from the
    // fra rate if the year end effects fall within the Fra
    // period.  The input list of fras does not have to be
    // sorted and normalized.
    // (3) year end effects before first fra date or after
    // last fra maturity are just discarded, and no warning is issued.
    // User must readjust accordingly.
    // (Based on Alib function LmFRAListAdjustYearEnd)
    static void fraListAdjustYearEnd
    (
        const vector<YearEndEffect>& yel,    // (I)
        const vector<FwdRateIvl>&    frl,    // (I) can be unsorted/unnormalised
        vector<FwdRateIvl>&          newList // (I/O) result
    );

// -----------------------------------------------------------------------------
// From tpcurve.h

    // Given an LMParabCurve whose dates, discs, parabolic
    // coefficients are defined, calculates the
    // disc at a desired date by using the parabolic coefficients
    // of the immediate interval.
    // (Based on Alib function LmTParabCurveInterpDisc)
    double interpDisc
    (
        DateTime desiredDate // (I) where interpolation is done
    )
    const;

    // Given an LMParabCurve whose dates, discs, parabolic
    // coefficients are defined, calculates either the
    // instantaneous or the overnight rate at a desired date
    // by using the parabolic coefficients of the immediate interval.
    double instRate
    (
        DateTime desiredDate, // (I) where interpolation is done
        int rateType          // (I) LM_INST_RATE, LM_INST_RATE_SLOPE or LM_OVERNIGHT_RATE
    )
    const;

    // Given an LMParabCurve whose dates, discs, parabolic
    // coefficients are defined, calculates the parameters for
    // instantaneous rate, forward rate or disc at a desired date.
    // (Based on Alib function LmTParabCurveInterpParam)
    void interpParam
    (
        DateTime desiredDate, // (I) where interpolation is done.
        double& prevDisc,     // (O) Discount at start of appropriate interval
        int&    daysIvl,      // (O) days from start of appropriate interval
        double  coeff[3]      // (O) curve coefficients of appropriate interval
    ) 
    const;

    // Finds the curve indices bracketing a desired date 
    // in a TParabCurve for interpolation.
    // (Based on Alib function LmTParabCurveIndexFind)
    void indexFind
    (
        DateTime date, // (I) Date to intepolate at
        int& low,      // (O) lower bracketing point
        int& high      // (O) upper bracketing point
    )
    const;

    // Predicates for rate point date comparison
    static bool compareRatePt
    (
        const RatePt& rp1, // (I) rate point 1
        const RatePt& rp2  // (I) rate point 2
    );

// -----------------------------------------------------------------------------
// From yldgen.h

    // Given a start Date, the instantaneous start rate
    // and discount to start of interval, 
    // end date of next interval, and the discount factor to
    // the maturity of the next interval compute its parabolic
    // coefficients.
    // (Based on Alib function LmCalcParabolicCoeffsTParabRatePt)
    static void calcParabolicCoeffsRatePt
    (
        RatePt&         rp,            // (I/O) parabolic rate interval
        const DateTime& startDate,     // (I) start date of the interval
        double          startDisc,     // (I) discount factor to start of period
        double          startInstRate, // (I) instantaneous forward rate at start of period
        const DateTime& nextDate,      // (I) date for following interval
        double          nextDisc,      // (I) discount factor to the end of next period
        double&         endInstRate,   // (O) instantaneous forward rate at end of period
        HolidayConstSP  noAccrue       // (I) dates with no accrual
    );

    // Given a start Date, the instantaneous start rate
    // and discount to start of interval, compute its linear coefficients.
    // (Based on Alib function LmCalcLinearCoeffsTParabRatePt)
    static void calcLinearCoeffsRatePt
    (
        RatePt&         rp,            // (I/O) parabolic rate interval
        const DateTime& startDate,     // (I) start date of the interval
        double          startDisc,     // (I) discount factor to start of period
        double          startInstRate, // (I) instantaneous forward rate at start of period
        double&         endInstRate,   // (O) instantaneous forwad rate at end of period
        HolidayConstSP  noAccrue       // (I) dates with no accrual
    );

    // Given a start Date, and discount to start of interval,
    // compute its flat coefficients.
    // (Based on Alib function LmCalcFlatCoeffsTParabRatePt)
    static void calcFlatCoeffsRatePt
    (
        RatePt&         rp,        // (I/O) parabolic rate interval
        const DateTime& startDate, // (I) start date of the interval
        double          startDisc, // (I) discoiunt factor to start of period
        HolidayConstSP  noAccrue   // (I) dates with no accrual
    );

    // Given a Fra maturing after the curve end
    // it gives back the date and discount which follows.
    // (Based on Alib function LmNextFraDisc)
    void nextFraDisc
    (
        const FwdRateIvl& fra,  // (I) fra rate interval 
        DateTime&         date, // (O) date for next discount point
        double&           disc  // (O) next disc factor
    )
    const;

    // Reintroduce the year end effects into the curve
    // (Based on Alib function LmAddYEToCurve)
    void addYE
    (
        const vector<YearEndEffect>& yel,      // (I)
        HolidayConstSP               noAccrue  // (I) dates with no accrual (Not 
                                               //     necessarily equal to this->noAccrue
                                               //     which includes ye dates 
    );

    // Constructs a curve given a FwdRateIvl list
    // and a YearEndEffect list.
    // (Based on Alib function LmGenTParabCurveFromFras)
    void genFromFras
    (
        const DateTime&              bd,               // (I) base date for curve
        const vector<FwdRateIvl>&    frl,              // (I)
        const vector<YearEndEffect>& yel,              // (I)
        HolidayConstSP               na,               // (I) dates with no accrual
        const vector<DateTime>&      flatDates,        // (I) the rate will be constant in between flatRateDates,
                                                       //     last date must be <= flatParabBdryDate
        const DateTime&              flatParabBdryDate // (I) date at which interpolations transitions
                                                       //     from flat to parabolic
    );

    // Adds the next flat period to the current built curve
    // (Based on Alib function LmAddNextFlatPeriod)
    void addNextFlatPeriod
    (
        const vector<DateTime>&   flatDates,    // (I) list of flat dates still to cover
        size_t&                   idxFlatDate,  // (I/O) current flat date index
        size_t                    numFlatDates, // (I) number of flat dates
        DateTime&                 lastFlatDate, // (I/O) last flat date covered
        double&                   lastDisc,     // (I/O) discount up to lastFlatDate
        const vector<FwdRateIvl>& fraList,      // (I) fras to be included
        size_t&                   idxFra        // (I/O) index of Fra being used currently
    );

    // Adds the next fra period to the current built curve
    // using flat forwards.
    // (Based on Alib function LmAddNextFraFlat)
    void addNextFraFlat
    (
        DateTime&                 lastDate, // (I/O) last date entered into curve
        double&                   lastDisc, // (I/O) discount up to lastDate
        const vector<FwdRateIvl>& fraList,  // (I) fras to be included
        size_t&                   idxFra    // (I/O) index of Fra being used currently
    );

    // Introduces the first period of the curve as
    // a linear segment.  Done for fully Parabolic curves.
    // (Based on Alib function LmAddFirstFraLinear)
    void addFirstFraLinear
    (
        const vector<FwdRateIvl>& fraList,      // (I) fras to be included
        double&                   currInstRate  // (O) instantaneous rate at the end of first period
    );

    // Introduces the next period of the curve as
    // a parabolic segment.
    // (Based on Alib function LmAddNextFraParab)
    void addNextFraParab
    (
        DateTime&                 lastDate,     // (I/O) last date entered into curve
        double&                   lastDisc,     // (I/O) discount up to lastDate
        double&                   lastInstRate, // (I/O) instantaneous forward rate at lastDate
        const vector<FwdRateIvl>& fraList,      // (I) fras to be included
        size_t&                   idxFra        // (I/O) index of Fra being used currently
    );

    // Introduces the next period of the curve as
    // a parabolic segment.
    // (Based on Alib function LmAddBdryFraFlatParab)
    void addBdryFraFlatParab
    (
        const DateTime&           prevStartDate, // (I) date at which the flat period started
        DateTime&                 lastDate,      // (I/O) last date entered into curve
        double&                   lastDisc,      // (I/O) discount up to lastDate 
        double&                   lastInstRate,  // (I/O) instantaneous forward rate at lastDate
        const vector<FwdRateIvl>& fraList,       // (I) fras to be included
        size_t&                   idxFra         // (I/O) index of Fra being used currently
    );

    // Computes flat forward rate from discount factors
    // with holidays for time intervals.  Flat forward
    // rate satisfies:
    //
    // exp(-FFRate*(business days(t1,t2))=
    // Discount(t2)/Discount(t1)
    // Returns SUCCESS or FAILURE if the time interval
    // is zero or negative or rates are wacky.
    // (Based on Alib function LmGetFlatFwdFromDiscountHoliday)
    static double getFlatFwdFromDiscountHoliday
    (
        const DateTime& leftDate,  // (I)
        double          leftRate,  // (I)
        const DateTime& rightDate, // (I)
        double          rightRate, // (I)
        HolidayConstSP  holiday    // (I) used for time between 1st, 2nd dates
    );

// -----------------------------------------------------------------------------

    vector<RatePt> ratePts;  // Dates & rates
    DateTime baseDate;       // Discount date
    bool areCoeffsValid;     // True if parabolic coefficients have
                             //   all been calculated or defined via
                             //   a constructor
    HolidayConstSP noAccrue; // All dates for which there is no accrual
};

DRLIB_END_NAMESPACE

#endif
