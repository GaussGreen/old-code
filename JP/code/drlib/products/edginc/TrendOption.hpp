//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TrendOption.hpp
//
//   Description : Option on the number of times the return on an asset changes direction
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : November 25, 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_TREND_HPP
#define EDR_TREND_HPP

#include "edginc/Generic1Factor.hpp"
#include "edginc/ClosedFormLN.hpp"
// #include "edginc/SampleList.hpp"
#include "edginc/LastSensDate.hpp"


DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL TrendOption: public Generic1Factor, public CClosedFormLN::IIntoProduct, 
                   public LastSensDate {
public:
    static CClassConstSP const TYPE;

    static const string SPREAD_NONE;
    static const string SPREAD_UP;
    static const string SPREAD_DOWN;
    static const string SPREAD_SYMMETRIC;
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

    void getDigStrikes(double strike, double *lowStrike, double *highStrike) const;

private:
    friend class TrendOptionHelper;
    friend class TrendOptionClosedForm;

    TrendOption();
    TrendOption(const TrendOption& rhs);
    TrendOption& operator=(const TrendOption& rhs);
    static void load(CClassSP& clazz);

    void addOutputRequests(Control* control,
                           Results* results,
                           const double& fairValue,
                           const double& indVol) const;

protected:
    TrendOption(CClassConstSP clazz);

    bool                 isCall;
    double               strike;
    double               trendMult;              // multiplies the number of trend switches in the payoff. redundant
    DateTime             maturityDate; // $unregistered
    DateTimeArraySP      monitorDates;           // One date for sample
    DoubleArraySP        histMonSamples;         // level of u/l on hist sample dates
    string               spreadType;             // increasing, decreasing, or none
    double               strikeSpread;           // spread for barrier call spread approx

};

typedef smartPtr<TrendOption> TrendOptionSP;

DRLIB_END_NAMESPACE
#endif
