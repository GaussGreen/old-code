//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : KnockInFwd.hpp
//
//   Description : Knock-in forward
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : October 25, 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_KIFWD_HPP
#define EDR_KIFWD_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/ClosedFormLN.hpp"
// #include "edginc/SampleList.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/PhysicalDelivery.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL KnockInFwd: public Generic1Factor, public CClosedFormLN::IIntoProduct, 
                  public LastSensDate {
public:
    static CClassConstSP const TYPE;

    static const string SPREAD_NONE;
    static const string SPREAD_INCREASING;
    static const string SPREAD_DECREASING;
    static const string SPREAD_SYMMETRIC;
    static const string SPREAD_SYMMETRICPCT;
    static const double MINIMUM_ALPHA;

    /** instrument validation */
    virtual void Validate();

    /** retrieve market data needed by VolVarSwap - just valueDate, asset and
        discount yield curve */
//    void GetMarket(const IModel*          model, 
//                   const CMarketDataSP    market);

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

    /** indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*       model);

    virtual void validatePop2Object();

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);
  
//    virtual DateTime getValueDate() const;
 
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** price a dead instrument until settlement - exercised, expired, knocked out etc.
    returns true if it is dead (and priced), false if it is not dead */
    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    void getDigStrikes(int index, double *lowStrike, double *highStrike) const;
    
    bool getUpDownBool() const; // up = true, down = false
    double getUpDownDouble() const; // up = true, down = false

	// useful utility functions
	bool isEndMonitoringPeriod(const DateTime& testDate) const;
	const DateTime& getPayDate(const DateTime& today) const;
	DateTime getNextPaymentDate(const DateTime& today) const;
	DateTime getFirstMonitoringDate(const DateTime& settleDate) const;

    // gets relevant info for monitoring period containing evalDate
// gets relevant info for monitoring period containing evalDate
// gets relevant info for monitoring period containing evalDate
    void getPeriodHistContractsAndStrikes(const DateTime& evalDate, 
                                      DoubleArray&                pdContracts, 
                                      DoubleArray&                pdStrikes, 
                                      double&                     putValue) const;
    
    // gets num contracts for period.  Required for computing num shares for divs.
    double getPeriodHistContracts(const DateTime& evalDate) const;

    // returns value of knocked in forward contracts for all known periods
    double productValue(double spot,
                                const DoubleArray& pdContracts, 
                                const DoubleArray& pdStrikes) const;

    // returns value of knocked in forward contracts for one known period
    double productValue(double spot,
                                double nbContracts,
                                double strike) const;

    // sets aggregateNumContracts to the sum of pdContracts, aggregateStrike to the center of mass of strikes, weighted
    // by pdContract.  
    void getDeliveryInfo(const DoubleArray& pdContracts, 
                                const DoubleArray& pdStrikes, 
                                double &aggregateNumContracts, 
                                double &aggregateStrike) const;
    
    
    // true if div treatment is "PassThroughOnPayDate" or "PassThroughOnPayDateDollar"
    bool includeDivs() const;

	// MATRIX functions
	void getMATRIXinfo(CashFlowArray &cfl, PhysicalDeliveryArray &pda) const;

	DateTimeArraySP getPaymentDates() const;

private:
    friend class KnockInFwdHelper;
    friend class KnockInFwdClosedForm;

    KnockInFwd();
    KnockInFwd(const KnockInFwd& rhs);
    KnockInFwd& operator=(const KnockInFwd& rhs);
    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

protected:
    KnockInFwd(CClassConstSP clazz);


    DateTimeArraySP      monitorDates;           // One date for eack ki fwd
    DoubleArraySP        barrierLevels;          // ki Barrier
    DoubleArraySP        numberOfShares;         // number of shares underlying ki per unit of contract
    DateTimeArraySP      forwardMatDates;        // forward trade dates
    DoubleArraySP        forwardStrikes;         // strike of forward
    DoubleArraySP        histMonSamples;         // level of u/l on hist sample dates
    DoubleArraySP        histTradeSamples;       // level of u/l on hist trade dates
    string               spreadType;             // increasing, decreasing, or none
    double               strikeSpread;           // spread for barrier call spread approx
    string               tradeDateStyle;         // OneForOne, AllOnOrBefore
    string               upDown;                 // u/d up and in or down and out
    bool                 usePutSchedule;         // include the puts
    DoubleArraySP        putStrikes;             // 
    string               divTreatment;           // none, passThroughEx, passThroughPay, passThroughAtSettle

    // Transient
    DateTimeArraySP      fullForwardMatDates;    // full forward trade date list
    DoubleArraySP        fullHistTradeSamples;   // full level of u/l on hist trade dates
};

typedef smartPtr<KnockInFwd> KnockInFwdSP;

DRLIB_END_NAMESPACE
#endif
