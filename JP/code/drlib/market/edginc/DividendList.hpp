//   Filename    : DividendList.hpp
//
//   Description : DividendList representation
//
//   Author      : Stephen Hope
//
//   Date        : 08 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef DIVIDENDLIST_HPP
#define DIVIDENDLIST_HPP

#include "edginc/Object.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Dividend.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/MuSpecial.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/OutputRequestUtil.hpp"


using namespace std;

DRLIB_BEGIN_NAMESPACE

class IYieldCurve;
class CAsset;

class MARKET_DLL DividendList : public CObject {
public:
    static CClassConstSP const TYPE;
    friend class DividendListHelper;
    friend class DivListCreateAddin;
    friend class Equity;    // I only need to make Equity::convertDollarDivs(...)
                            // friend of DividendList, but I'm confused about the syntax

    /** Constructor */
    DividendList(const DividendArray&  divArray);

    /** Create a DividendList full of all divs due to be paid > startDate, 
        and <= endDate. Returns an empty dividend list if no dividends match 
        this condition. It is the callers responsibility to free the memory 
        allocated by this function */
    DividendList* getAllDivsBetweenDates(
        const DateTime&   startDate,         // (I) No divs before this date 
        const DateTime&   endDate)const;     // (I) No divs after this date */

    /** same as getAllDivsBetweenDates, but ignores continuous dividends */
    DividendList* getAllDiscreteDivsBetweenDates(
        const DateTime&   startDate,         // (I) No divs before this date 
        const DateTime&   endDate)const;     // (I) No divs after this date */

    /** return a DateTimeArray made up of all the pay dates in the dividend list */
    DateTimeArrayConstSP getPayDates()const;

    /** return a DateTimeArray made up of all the ex-div dates in the dividend list */
    DateTimeArrayConstSP getExDivDates()const;

    /** return a doubleList made up of all the amounts in a div list */
    DoubleArrayConstSP getDivAmounts()const;

    /** return a dividend list which is the result of merging the passed
        dividend list with 'this' dividend list. It is the callers 
        responsibility to free the memory allocated by the function. */
    DividendList* merge(const DividendList& otherList)const;

    /** scales all dividends by multFact ie amount -> amount*multFact */
    void scale(double multFact);

    /** MU_PARALLEL shift method. Note that the Equity class inherits
        from MuParallel and this method is wrapped in the equity
        shift method */
    bool sensShift(MuParallel* shift, const DateTime& valueDate); 
  
    /** MU_S shift method. Note that the Equity class inherits
        from MuSpecial and this method is wrapped in the equity
        shift method */
    bool sensShift(MuSpecial* shift, 
                   double spotPrice,
                   const DateTime& expiry); 
    
    /** MU_POINTWISE shift method. Note that the Equity class inherits
        from MuPointwise and this method is wrapped in the equity
        shift method */
    bool sensShift(MuPointwise* shift, 
                                 const DateTime& valueDate,
                                 DateTime & bucketStartDate,
                                 DateTime & bucketEndDate);

    /* Access the dividend array */

    /** return the dividend array  */
    const DividendArray& getArray()const;

    /** finds any continuous dividend stricly != 0 paying at theDate -
        if there is one returns true and the amount else returns false
        and 0.0 (A cts div of 0.0 will cause false to be returned) */
    bool ctsDivPayingAtDate(
        const DateTime&      theDate,
        double&              ctsDiv) const; /* (O) ctsDiv (cts basis) 
                                               paying at theDate */

    /** Returns true if this DividendList contains any CONTINUOUS dividends,
	along with the date and index of the first continuous dividend */
    bool hasContinuousDividend(DateTime& firstContDivDate)const; // (O) ex-div of first cont div 

    /** Returns true if this DividendList contains any AMOUNT dividends */
    bool hasDollarDividend() const;

    /** Returns true if this DividendList contains any PERCENT dividends */
    bool hasYieldDividend() const;

    /** Given a dividend list, returns a dividend list where all dollar / yield 
        dividends falling on a same date / time are aggregated. If more than 1
        continuous dividends fall on the same date / time, only the last dividend
        is kept. */
    static DividendList* aggregate(const DividendList& input,
                                   const IYieldCurve*  yc);

    /** Given a date, find the corresponding dividend with exDate
        immediately after that date.
        Returns 0 for the following cases:
        DividendList is empty
        theDate >= last exDate */  
    DividendConstSP getNextDivFromDate(const DateTime& date) const;

    void setElement(int idx, const Dividend& Dividend);

    /** Converts any yield dividends to absolute amounts using fwd price of
        asset. Validates against continuous dividends. Returns true if
        any yield divs were found */
    bool convertYieldDivs(const CAsset* asset);

    /** Converts any dollar dividends to yields using fwd price of
        asset. Validates against continuous dividends. Returns true if
        any dollar divs were found */
    bool convertDollarDivs(const CAsset* asset);

    /** Use to adjust divs (primarily pay dates) in adjustDivs(). In
        other words a callback mechanism. It is used by the 
        DividendCollector. */
    class MARKET_DLL IDivAdjuster{
    public:
        /** Adjusts the dividend. Useful for adjusting pay dates and/or
            the size. */
        virtual void adjustDividend(Dividend& div) = 0;
        virtual ~IDivAdjuster();
    };

    /** Adjusts each of the dividends using the adjustDividend() method
        on adjuster. */
    void adjust(IDivAdjuster* adjuster);

    /** Scales any dollar dividends in the future (based upon pay
        date) by fwd rate whilst those in the past are scaled by the
        spot. Yield divs are scaled by fwd rate(pay date)/fwd rate(ex
        div date). Continuous divs are left alone */
    void makeStruck(const CAsset*    fxAsset,
                    const DateTime&  valDate);

    /** Update known cashflows object with details of dividends (dollar divs
        are known, whilst yields aren't */
    void updateKnownCashFlows(OutputRequestUtil::KnownCashFlows& knownCFs,
                              const string&                      isoCode,
                              const DateTime&                    valueDate);
    
///////////////  NOT YET IMPLEMENTED ///////////////////////////////

#if 0
    /** tweak all dividends between two dates by a percentage amount */
    double tweakDivsBetweenDates(
        const DateTime&  startDate,  // (I) shift after this date 
        const DateTime&  endDate,    // (I) shift up to this date 
        double           shift);     // (I) shift (d => d*(1+shift)) 
   
    /** From a div list get all the dividends whose pay date satisfies
        startDate < pay Date <= endDate. */
    DividendList* dividendFuturePayDates(
        const DateTime&       startDate)const;   // (I) first div date

    /** Modifiy the pay dates on a given dividend list so as to make them
        represent a cumulative mode of payment. */
    void dividendShiftForCumulative(
        const DateTimeArray&   cumDates,     // (I) List of cumulative dates 
        const DateTimeArray&   payDates);    // (I) List of payment dates 

#endif

private:
    DividendList();
    DividendList& operator=(const DividendList& rhs);
    DividendList(const DividendList& rhs);
    
    virtual void validatePop2Object();

    /** Given a date, find the corresponding dividend idx with exDate
        immediately after that date , exDivDate[idx] STRICT > date.
        Returns -1 for the following cases:
        DividendList is empty
        theDate >= last exDate */  
    int findDivArrayIndexFromDate(
        const DateTime&   theDate)const;       // (I) the date

    /** Same as findDivArrayIndexFromDate except returns an iterator. 
        Returns divArray.end() for the following cases:
        DividendList is empty
        theDate < 1st exDate 
        theDate > last exDate */
    vector<Dividend>::iterator findDivArrayIteratorFromDate(
        const DateTime&   theDate);       // (I) the date
    
    /** insert a dividend at a given index, ex-div date = pay date. 
        Returns an iterator to the newly inserted dividend */
     vector<Dividend>::iterator dividendInsert(
         vector<Dividend>::iterator  idx,  // (I) array index to place div 
         Dividend::TDividendType  divType,       // (I) 
         const  DateTime&   exDivDate,           // (I) 
         double             divAmount);          // (I)

    /** tweak all divs for idx1 <= div < idx2 - note array indexes. 
        No validation */
    void dividendTweakBetweenIndexes(vector<Dividend>::iterator  idx1, 
                                     vector<Dividend>::iterator  idx2, 
                                     double  shiftSize);   
    
    
    /** given start and end date, report cts divs in effect at those time 
        plus index of appropriate dividend which marks the first div after 
        start date and first div after end date. Note index is array index, 
        also indexes can be equal to divList->fNumItems (ie one past last div).
        Validation only for start date and end date*/
    void findCtsDivsForBucket(
        const DateTime&  startDate,      // (I) start date for bucket 
        const DateTime&  endDate,        // (I) end date for bucket (>= start) 
        bool&            ctsDivAtStart,  /* (O) true - there is a cts div at 
                                            start */
        double&          ctsDivStart,    // (O) what value it is 
        vector<Dividend>::iterator&   divStartIdx, /* (O) index for first
                                                      div after start */
        bool&            ctsDivAtEnd,   /* (O) true - there is a cts div at
                                           end date */
        double&          ctsDivEnd,     // (O) what value it is 
        vector<Dividend>::iterator&   divEndIdx); /* (O) index for
                                                     first div after
                                                     end */

    /** Sorting method. 
        Re-order order of dividends falling on a same day in the following order
            1. PERCENT
            2. AMOUNT
            3. CONTINUOUS */
    void sort();

    /** Sorting method. 
        Transported from Numerical recipes */
    static void piksr2(int iDivBegin, 
                       int iDivEnd, 
                       CIntArray& divTypes, 
                       CDividendArray& divs);

    /* Data Members */
    DividendArraySP  divArray;
};

typedef smartPtr<DividendList>               DividendListSP;
typedef smartConstPtr<DividendList>          DividendListConstSP;
#ifndef QLIB_DIVIDENDLIST_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<DividendList>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<DividendList>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<DividendList>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<DividendList>);
#endif
typedef DividendList                         CDividendList;
typedef array<DividendListSP, DividendList>  DividendListArray;
#ifndef QLIB_DIVIDENDLIST_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<DividendListSP _COMMA_ DividendList>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<DividendListSP _COMMA_ DividendList>);
#endif

DRLIB_END_NAMESPACE

#endif
