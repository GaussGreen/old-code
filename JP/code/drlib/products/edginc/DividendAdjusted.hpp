//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DividendAdjBase.hpp
//
//   Description : Dividend adjusted option instrument
//
//   Date        : 11/20/2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VANILLADIVADJ_HPP
#define EDR_VANILLADIVADJ_HPP
#include "edginc/Class.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Asset.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SensitiveStrikes.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/InstrumentSettlement.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Theta.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Average.hpp"

DRLIB_BEGIN_NAMESPACE

/** Shell structure for spot/hybrid/ratio payoffs */
class PRODUCTS_DLL DividendAdjBase: public Generic1Factor, 
                       public LastSensDate {
public:
    static CClassConstSP const TYPE;

    /** instrument validation */
    virtual void Validate();

    /** retrieve market data needed by Vanilla - just valueDate, asset and
        discount yield curve */
    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market);


    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which 
        this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*      model);

    virtual CSensControl* AlterControl(const IModel*       modelParams,
                                       const CSensControl* sensControl) const;

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    // if useFwd is true, record forward to fwdDate in setting fx and df */
    void thetaHelper(Theta*          shift, 
                     const DateTime& newDate, 
                     bool            useFwd = false, 
                     const DateTime& fwdDate = DateTime(0,0));
  
    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    // override base
    virtual void validatePop2Object();

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
 //   virtual bool priceDeadInstrument(CControl* control, CResults* results) const;


    ////////////////////////// some helper methods ////////////////////
    // get zero rate from valueDate to date
    static double getZeroRate(const DateTime& date, const DividendAdjBase* inst);
    // get fwd zero rate from date1 to date2
    static double getFwdZeroRate(const DateTime& date1, 
                                 const DateTime& date2, 
                                 const DividendAdjBase* inst);

    // get spot fx rate for the date
    static double getFXValue(const DateTime& date, const CAsset* struckAsset);
    //get historic data
    static double getHist(const CashFlowArray& histValues, DateTime date);
    // get growth factor from date1 to date2
    static double growthFactor(const DateTime& date1, const DateTime& date2, const DividendAdjBase* inst);
    // returns dollar div (convert yield if needed)
    static double div2dollar(const Dividend& div, const CAsset* asset);
	// calculate sub total using data from subDates & subValues, summed according to periodDates.
	// last date of periodDates is inputted separately to allow inclusive/exclusive counting
	static void subTotal(const DateTimeArray& periodDates, const DateTime& lastDate, bool lastDateInclusive, 
                     const DateTimeArray& subDates, const DoubleArray& subValues, DoubleArray& total);
	/** finds corresponding pay dates for the supplied ex dates, using the supplied divArray.  
		If the pay date cannot be found, uses ex date. */
	static void lookUpPayDates(const DateTimeArray& exDates, const DividendArray& divArray, DateTimeArray& foundPayDates);
	/** returns true if ex date is in the div array, whence foundPayDate is assigned to exDate's pay date. 
				false otherwise, with foundPayDate set to 0 */
	static bool findDivPayDate(const DateTime& exDate, const DividendArray& divArray, DateTime& foundPayDate);

    // util function to calc proportional strike weight adjustments.
    // first sum to each sub date and the 1-sum.
    static void propAdjust(const DateTimeArray& periodDates, const DateTimeArray& subDates,
                    const DoubleArray& subValues, DoubleArray& prop);

    static void processOutputRequest(const string&   outputRequest,
                                     const double value,
                                     Control*        control, 
                                     CResults*       results);

    // returns the difference between a div and an assumed div, according to our convention of multiple divs / qtr.
    static double computeDivDiff(int    numDvPd,   // number of divs in the period
                                 double amtDvPd,   // total dollar amount of divs in the period
                                 double amtDv,     // dollar amount of the div
                                 double amtAssm);   // dollar amount of the assumed divs in the period


    double computeDivDiff(const Dividend& exDiv);

protected:

    // calculate FX spot and growth factor for eq div dates
    void calcFXandDFdiv(DoubleArray& fx, DoubleArray& gf) const;

    // calculate FX spot and growth factor for assumed dates
    void calcFXandDFassumed(const DoubleArray& periodSum,
                                   DoubleArray& fx,
                                   DoubleArray& gf) const;

    // build a average spot instrument
    AverageSP buildAvgSpot(SampleListConstSP outSample, double strike) const;

    /** get equity dividends */
    DividendArray getEqDivs() const;

    /** if eq divs contain yield - will throw for continous yield */
    bool hasYield() const;

    /* true = this is a ex-div date or last date of assumed period with no div */
    bool recordHist() const;

	/** for assumedDivsInStruckCcy case, we convert assumed divs into base ccy before pricing */
	CashFlowArraySP assumedDivs2BaseCcy() const;

    /** gets num div dates on the fly.  This allows proper treatment for MuSpecial */
    CashFlowArraySP getNumDivDates(const DividendArray& currentEqDivs) const;

private:
    friend class DividendAdjBaseHelper;
    friend class DividendAdjustedClosedForm;
    friend class DividendAdjustedFLClosedForm;

    friend class VanillaTSO;

    DividendAdjBase(const DividendAdjBase& rhs);
    DividendAdjBase& operator=(const DividendAdjBase& rhs);

protected:
	DividendAdjBase(CClassConstSP const type);  // for derived class

    // interface 
    double          strike;
    DateTime        matDate;
    CashFlowArray   assumedDivs;
    CashFlowArray   histFXRates;
    CashFlowArray   histDiscRates;
    CashFlowArray   numExDivDates;

    // optional
	bool			avgWeightAdjusted;
    bool            usePayDate; // to do : what's the exact definition ?

	// new optional inputs for scaleTypeA 
	string			scaleType;
	CashFlowArray   spotBeforeExDates;
	CashFlowArray   spotBeforeAssumedDates;

    // transient in base, but exposed for new Flow thru child class
    bool            isCall;
    string          divadjType;
    bool            fwdAdj;
    SampleListSP    avgOut;
	bool			assumedDivsInStruckCcy;

    // transient
    // base case avg spot instrument, used e.g. to call getSensitiveStrikes
    AverageSP       avgInst;
    DividendArray   eqDivArr;
};

typedef smartPtr<DividendAdjBase> DividendAdjBaseSP;

DRLIB_END_NAMESPACE
#endif
