#include "edginc/config.hpp"
#include "edginc/EffectiveCurve.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"

DRLIB_BEGIN_NAMESPACE

/** An aberration -- to be removed. Too much to get rid of in first phase
    of refactoring CDO code. Please do not use this class unless writing
    a payoff that uses ConvolutionProduct */
class FlatFwdZeroCurve: public CObject,
                        public virtual IDiscountCurve //minimal implementation
{
    typedef DateTimeArray::const_iterator DateIter;
    typedef DoubleArray::const_iterator   DoubleIter;
public:

    //-------------------------
    // FlatFwdZeroCurve methods
    //-------------------------
    static CClassConstSP const TYPE;

    FlatFwdZeroCurve(const DateTime& baseDate,
                     const DateTime& spotDate,
                     DateIter        dateBegin,
                     DateIter        dateEnd, 
                     DoubleIter      discFactorBegin,
                     DoubleIter      discFactorEnd)
    : CObject(TYPE),
      today(baseDate),
      spotDate(spotDate),
      date(dateBegin, dateEnd), 
      zero(discFactorBegin, discFactorEnd)
    {
        // convert from discount factors to zero rates
        Actual365F act365F;
        for (int i = 0; i < zero.size(); i++)
        {
            zero[i] = RateConversion::discountToRate(zero[i],
                                                    today, date[i],
                                                    &act365F,
                                                    CompoundBasis::CONTINUOUS);
        }
    }

    FlatFwdZeroCurve(const DateTime&   today,
                     YieldCurveConstSP yc)
    : CObject(TYPE),
      today(today),
      spotDate(yc->getSpotDate()),
      date(yc->zeroDates())
    {
        // remove 'spotDate' from zero dates (zero rate is undefined for spot date)
        for (vector<DateTime>::iterator iter(date.begin()); iter != date.end();)
        {
            if (!iter->equals(spotDate))
            {
                ++iter;
            }
            else
            {
                iter = date.erase(iter);
            }
        }

        // resize zero rates array
        zero.resize(date.size());

        Actual365F act365F;
        double annualRate;
        for (int i = 0; i < date.size(); i++)
        {
            // convert from annual to continuous
            annualRate = yc->zero(date[i]);
            zero[i] = RateConversion::rateConvert(
                spotDate,
                date[i],
                annualRate,
                &act365F,
                CompoundBasis::ANNUAL,
                &act365F,
                CompoundBasis::CONTINUOUS);
        }
    }


    void toCcmCurve(DoubleArray &tR,       // (O) time points
                    DoubleArray &ZR) const // (O) zero prices
    {
        tR.resize(date.size());
        ZR.resize(date.size());
        for (int i = 0; i < date.size(); i++){
            tR[i] = Actual365F::yearFraction(today, date[i]);
            ZR[i] = zeroPrice(date[i]);
        }
    }

    static FlatFwdZeroCurveSP Create(const DateTime&           today,
                                     const YieldCurveWrapper&  yc)
    {
        return FlatFwdZeroCurveSP(new FlatFwdZeroCurve(today, yc.getSP()));
    }

    static FlatFwdZeroCurveSP Create(
        const DateTime&              today,
        const ICDSParSpreadsWrapper& cdsParSpreads,
        const YieldCurveWrapper&     discount)
    {
        DefaultRatesSP defaultRates(cdsParSpreads->defaultRates());
        CashFlowArraySP cleanFwdSpreadCurve(defaultRates->getCleanSpreadCurve());
        return FlatFwdZeroCurveSP(new FlatFwdZeroCurve(today,*cleanFwdSpreadCurve));
    }
    
    //// empty constructor
    FlatFwdZeroCurve() : CObject(TYPE){}

    virtual ~FlatFwdZeroCurve(){}

    //// usual clone method - hopefully not called now
    virtual IObject* clone() const
    {
        //static const string method("FlatFwdZeroCurve::clone");
        //throw ModelException(method, "Not yet implemented.");

        FlatFwdZeroCurve* copy = new FlatFwdZeroCurve();
        copy->date = date; // stl copy of a Vector
        copy->zero = zero; // stl copy of a Vector
        copy->today = today;
        copy->spotDate = spotDate;
        return copy;
    }

    //// returns pv to supplied date
    virtual double zeroPrice(const DateTime &d) const
    {
        double r1 = interpolate(today);
        double z1 = computeDiscount(today, r1);
        double r2 = interpolate(d);
        double z2 = computeDiscount(d, r2);
        return z2 / z1;
    }

    //// constructor from fwd rates
    FlatFwdZeroCurve(const DateTime&        today,
                     const CashFlowArray&   fwdRates)
    : CObject(TYPE),
      today(today),
      spotDate(today),
      date(fwdRates.size()),
      zero(fwdRates.size())
    {
        double integratedRate;
        double timeSum;
        for (int i = 0; i < date.size(); i++){
            date[i] = fwdRates[i].date;
            // convert from forward to spot rates
            if (i == 0){
                zero.front() = fwdRates.front().amount;
                timeSum = Actual365F::yearFraction(today, date.front());
                integratedRate = zero.front() * timeSum;
            } else {
                double yearFrac = Actual365F::yearFraction(date[i-1], date[i]);
                integratedRate +=  yearFrac * fwdRates[i].amount;
                timeSum += yearFrac;
                zero[i] = integratedRate/timeSum;
            }
        }
    }

    //-----------------------
    // IDiscountCurve methods
    //-----------------------
    virtual string getCcy() const
    {
        static const string method("FlatFwdZeroCurve::getCcy");
        throw ModelException(method, "Not yet implemented.");
    }

    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const
    {
        static const string method("FlatFwdZeroCurve::pv(date1,date2)");
        throw ModelException(method, "Not yet implemented.");
    }

    virtual double pv(const DateTime& date) const
    {
        static const string method("FlatFwdZeroCurve::pv(date)");
        throw ModelException(method, "Not yet implemented.");
    }

    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const
    {
        static const string method("FlatFwdZeroCurve::pv(cashflows,date)");
        throw ModelException(method, "Not yet implemented.");
    }

    class FFdiscountKey : public IDiscountCurve::IKey
    {
    public:
        FFdiscountKey(const FlatFwdZeroCurve* ffzc)
        {
            curve = ffzc;
        }

        //-----------------------------
        // IDiscountCurve::IKey methods
        //-----------------------------
        /** Returns the log of the risky discount factor between the two dates */
        virtual double calc(const DateTime&  loDate,
                            const DateTime&  hiDate)
        {
            //unoptimised
            double loPv = curve->zeroPrice(loDate);
            double hiPv = curve->zeroPrice(hiDate);
            return log(hiPv/loPv);
        }
    private:
        const FlatFwdZeroCurve* curve;
    };

    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const
    {
        return new FFdiscountKey(this);
    }

    virtual DateTimeArray zeroDates() const
    {
        static const string method("FlatFwdZeroCurve::zeroDates");
        throw ModelException(method, "Not yet implemented.");
    }

    //-----------------------------------------------
    // IGetMarket methods - exposed by IDiscountCurve
    //-----------------------------------------------
    virtual void getMarket(const IModel* model, const MarketData* market)
    {
        static const string method("FlatFwdZeroCurve::getMarket");
        throw ModelException(method, "Not yet implemented.");
    }

    static void load(CClassSP& clazz)
    {
        REGISTER(FlatFwdZeroCurve, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDiscountCurve);
        //IMPLEMENTS(IDiscountCurve::IKey);

        FIELD(today,    "today");
        FIELD(spotDate, "spot date");
        FIELD(date,     "term structure dates");
        FIELD(zero,     "term structure rates");
    }


protected:        
    //// constructor from 2 FlatFwdZeroCurves
    FlatFwdZeroCurve(FlatFwdZeroCurve& c1,
                     FlatFwdZeroCurve& c2)
    : CObject(TYPE),
      today(c1.today),
      spotDate(c1.spotDate)
    {
        if (!today.equals(c2.today)){
            throw ModelException("FlatFwdZeroCurve", "Curves have different "
                                "value dates");
        }
        date = DateTime::merge(c1.date, c2.date);
    }

private:
    double computeDiscount(const DateTime& date,
                           double          rate) const
    {
        Actual365F act365F;
        return RateConversion::rateToDiscount(
            rate,
            spotDate,
            date,
            &act365F,
            CompoundBasis::CONTINUOUS);
    }

    double interpolate(const DateTime &d) const
    {
        int dt = d.getDate();
        if (dt <= date.front().getDate()){
            return zero.front(); // extrapolate flat
        }
        // indentify where to interpolate from
        int lo;
        int hi;
        if (dt >= date.back().getDate()){
            // extrapolate flat forward off the end
            lo = date.size()-2; // let's hope there are >= 2 points
            hi = date.size()-1;
        } else {
            // Do a binary search to find the lower and upper bounds
            lo = 0;
            hi = zero.size() -1;
            while ((hi - lo) > 1) {
                int mid = (hi + lo) >> 1;  // compute a mid point
                if (dt >= date[mid].getDate()) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }
        if (date[lo].getDate() == dt) {
            return zero[lo];
        }
        if (date[hi].getDate() == dt) {
            return zero[hi];
        }
        // flat fwd interpolation (ie linear in r*t)
        double yf1 = Actual365F::yearFraction(spotDate, date[lo]);
        double yf2 = Actual365F::yearFraction(spotDate, date[hi]);
        double yf = Actual365F::yearFraction(spotDate, d);
        double rate = zero[lo] * yf1 + ((zero[hi] * yf2 - zero[lo] * yf1)/(yf2-yf1)) * (yf - yf1);
        rate /= yf;
        return rate;
    }

    ///// fields ////
    DateTime today;
    DateTime spotDate;
    DateTimeArray date;
    //// NB Rates are continuous quoted spot rates
    DoubleArray   zero;
};

/**
 * RiskyZCurve is just a combination of 2 FlatFwdZeroCurve
 * (zero rates and clean spreads).
 * It derives from FlatFwdZeroCurve to keep CDO code unchanged.
 * */
class RiskyZCurve: public FlatFwdZeroCurve {
public:    
    // Constructor
    RiskyZCurve(FlatFwdZeroCurveSP c1, FlatFwdZeroCurveSP c2)
    : FlatFwdZeroCurve(*c1, *c2),
      curve1(c1),
      curve2(c2)
    {}

    // Discount method
    virtual double zeroPrice(const DateTime &d) const
    {
        return curve1->zeroPrice(d) * curve2->zeroPrice(d);
    }

    // Clone - hopefully not used
    virtual IObject* clone() const
    {
        static const string method("RiskyZCurve::clone");
        throw ModelException(method, "Not yet implemented.");

        FlatFwdZeroCurve* c1clone = dynamic_cast<FlatFwdZeroCurve*>(curve1->clone());
        FlatFwdZeroCurve* c2clone = dynamic_cast<FlatFwdZeroCurve*>(curve2->clone());

        return new RiskyZCurve(
            FlatFwdZeroCurveSP(c1clone),
            FlatFwdZeroCurveSP(c2clone));
    }

private:
    // FIELDS
    FlatFwdZeroCurveSP curve1;
    FlatFwdZeroCurveSP curve2;
};

//-----------------------
// EffectiveCurve methods
//-----------------------

const string EffectiveCurve::FLAT_FORWARD = "FLAT_FORWARD";
const string EffectiveCurve::LINEAR = "LINEAR";

//default destructor
EffectiveCurve::~EffectiveCurve()
{}

//default constructor
EffectiveCurve::EffectiveCurve() : CObject(TYPE)
{}

////Constructor with a risk free discount curve
EffectiveCurve::EffectiveCurve(const DateTime&         valueDate,
                               const YieldCurveConstSP dsc,
                               const DateTimeArray&    expectedLossDates,
                               const DoubleArray&      expectedLosses,
                               const string&           interpStyle)
    : CObject(TYPE),
      valueDate(valueDate),
      interpolationStyle(interpStyle)
{
    //build the discount curve
    discount = IDiscountCurveSP::constCast(dsc); //used as is
    //and the flat fwd version
    ffDiscount.reset(new FlatFwdZeroCurve(valueDate, dsc));

    //copy the losses and timeline
    expLoss.reset(new DoubleArray(expectedLosses));
    timeline.reset(new DateTimeArray(expectedLossDates));

    //build the default rates
    buildDefRatesFromLossesAndTimeline();

    validate();
}

////Constructor with a risky discount curve
EffectiveCurve::EffectiveCurve(const DateTime&                valueDate,
                               const YieldCurveConstSP        dsc,
                               const SingleCreditAssetConstSP cpty, //combines with dsc to become risky
                               const DateTimeArray&           expectedLossDates,
                               const DoubleArray&             expectedLosses,
                               const string&                  interpStyle)
    : CObject(TYPE),
      valueDate(valueDate),
      interpolationStyle(interpStyle)
{
    //build the risky discount curve
    //...get the counterparty default rates
    DefaultRatesSP cptyDefaultRates(cpty->getParSpreadCurve()->defaultRates());
    //...and build an EffectiveCurve to use for risky discounting
    discount.reset(
        new EffectiveCurve(
            valueDate,
            dsc,
            cptyDefaultRates));

    //...and the flatFwd version
    CashFlowArraySP cleanFwdSpreadCurve(cptyDefaultRates->getCleanSpreadCurve());
    FlatFwdZeroCurveSP d = FlatFwdZeroCurveSP(
        new FlatFwdZeroCurve(valueDate, dsc));
    FlatFwdZeroCurveSP c = FlatFwdZeroCurveSP(
        new FlatFwdZeroCurve(valueDate, *cleanFwdSpreadCurve));
    ffDiscount.reset(new RiskyZCurve(d,c));

    //copy the losses and timeline
    expLoss.reset(new DoubleArray(expectedLosses));
    timeline.reset(new DateTimeArray(expectedLossDates));

    //build the default rates
    buildDefRatesFromLossesAndTimeline();

    validate();
}

////Cpty clean spreads are the losses
EffectiveCurve::EffectiveCurve(const DateTime&                valueDate,
                               const YieldCurveConstSP        dsc,
                               const SingleCreditAssetConstSP cpty,
                               const string&                  interpStyle)
    : CObject(TYPE),
      valueDate(valueDate),
      interpolationStyle(interpStyle)

{
    //build the discount curve
    discount = IDiscountCurveSP::constCast(dsc); //used as is
    //and the flat fwd version
    ffDiscount.reset(new FlatFwdZeroCurve(valueDate, dsc));

    //get the counterparty default rates
    defaultRates = cpty->getParSpreadCurve()->defaultRates();

    //build the losses
    buildLossesAndTimelineFromDefRates();

    validate();
}

////Constructor with a pre-built DefaultRates object
EffectiveCurve::EffectiveCurve(const DateTime&           valueDate,
                               const YieldCurveConstSP   dsc,
                               const DefaultRatesConstSP defRates)
    : CObject(TYPE),
      valueDate(valueDate),
      interpolationStyle(FLAT_FORWARD)
{
    discount = IDiscountCurveSP::constCast(dsc); //used as is
    ffDiscount.reset(new FlatFwdZeroCurve(valueDate, dsc));

    defaultRates = CONST_POINTER_CAST<DefaultRates>(defRates); //used as is

    //build the losses
    buildLossesAndTimelineFromDefRates();

    validate();
}

//convert the spot losses (discount factors) into forwards
void EffectiveCurve::buildDefRatesFromLossesAndTimeline()
{
    //default rates are continuous forwards; expLoss are discount factors
    Actual365F act365f;

    //first date is today - need to lose that
    int numSpots = timeline->size();
    DateTimeArraySP times = DateTimeArraySP(new DateTimeArray(numSpots-1));
    DoubleArraySP   spots = DoubleArraySP(new DoubleArray(numSpots-1));

    //convert from discount factor to spot rate
    for (int i=1; i<numSpots; i++)
    {
        DateTime thisDate = (*timeline)[i];
        double   thisDf   = (*expLoss)[i];
        double   spotRate = 0.0;
        //"unconstructed" losses can have 0's...so handle that
        if (thisDf != 0.0)
        {
            spotRate = RateConversion::discountToRate(
                thisDf, valueDate, thisDate, &act365f, CompoundBasis::CONTINUOUS);
        }

        //store
        (*times)[i-1] = thisDate;
        (*spots)[i-1] = spotRate;
    }

    //convert from spots to forwards
    DoubleArraySP fwds = RateConversion::spotsToForwards(
        valueDate, *(times.get()), *(spots.get()), &act365f);


    //build the default rates
    defaultRates.reset(
        new CDSHelper::CParSpreadDefaultRates(
            valueDate,
            *(times.get()), 
            *(fwds.get()),
            interpolationStyle));
}

//build up the timeline & losses from the default rates
void EffectiveCurve::buildLossesAndTimelineFromDefRates()
{
    //default rates are continuous forwards; expLoss are discount factors
    Actual365F act365f;
    //CashFlowArraySP fwdSpds = defaultRates->getCleanSpreadCurve();
    DateTimeArray fwdDates = defaultRates->getDates();
    DoubleArray fwdRates = defaultRates->getRates();

    int numFwds = fwdDates.size();

    timeline.reset(new DateTimeArray(numFwds+1));
    expLoss.reset(new DoubleArray(numFwds+1));

    //insert valueDate
    (*timeline)[0] = valueDate;
    (*expLoss)[0] = 1.0;

    //convert forwards into spots
    DoubleArraySP spots = RateConversion::forwardsToSpots(
        valueDate,
        fwdDates,
        fwdRates,
        &act365f);

    //convert spots into discount factors
    for (int i = 0; i < numFwds; i++)
    {
        (*expLoss)[i+1] = RateConversion::rateToDiscount(
            (*spots)[i], valueDate, fwdDates[i], &act365f, CompoundBasis::CONTINUOUS);
    }

}

void EffectiveCurve::validate()
{
    static const string method("EffectiveCurve::validate");

    //basic validation
    if (timeline->empty())
    {
        throw ModelException(method, "Timeline must be non empty");
    }
    if ((*timeline)[0] != valueDate)
    {
        throw ModelException(method, "First point of the timeline must be today");
    }

    //validate the interp type
    if (interpolationStyle == FLAT_FORWARD)
    {
        interpStyle = CCMPriceUtil::FLATFWD; //see comment against field in header
    }
    else if (interpolationStyle == LINEAR)
    {
        interpStyle = CCMPriceUtil::LINEAR; //see comment against field in header
    }
    else
    {
        throw ModelException(method, "Invalid interpolation style. Expected " +
                                     FLAT_FORWARD + " or " + LINEAR);
    }
}


////Utility method used in constructing the integrator
DateTimeArrayConstSP EffectiveCurve::buildIntegrationDates(
    const DateTime&      startDate,     //start of integration
    const DateTime&      endDate,       //end of integration
    const double         recoveryDelay,
    IDiscountCurveSP     discount,
    DateTimeArrayConstSP defaultDates) 
{
    static const string method = "EffectiveCurve::buildIntegrationDates";

    try {
        // get merged discount and default dates
        DateTimeArraySP mergedDates = 
            mergeDates(recoveryDelay, discount, defaultDates);

        //trim for our integration period
        DateTime::removeOutliers(*(mergedDates.get()), startDate, endDate);

        // then put startDate back in
        if (mergedDates->front().getDate() != startDate.getDate()) {
            mergedDates->insert(mergedDates->begin(), startDate);
        }

        //lastly put endDate back in
        if (mergedDates->back().getDate() != endDate.getDate()) {
            mergedDates->push_back(endDate);
        }

        //form new SP to return
        DateTimeArrayConstSP intDates = 
            DateTimeArrayConstSP(new DateTimeArray(*(mergedDates.get())));

        return intDates;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


////Utility method for merging dates
DateTimeArraySP EffectiveCurve::mergeDates(double               riskfreeDelay,
                                           IDiscountCurveSP     discount,
                                           DateTimeArrayConstSP defaultDates) 
{
    //combination of default dates and discount dates

    vector<const DateTimeArray*> dateArraysToMerge;
    dateArraysToMerge.reserve(2);
    
    //firstly discount dates
    DateTimeArray zcDates = discount->zeroDates();
    dateArraysToMerge.push_back(&zcDates);

    //adjust for recovery delay
    for (int i=0; i<zcDates.size();i++) {
        zcDates[i] = zcDates[i].rollDate(-riskfreeDelay);
    }

    // Default dates
    dateArraysToMerge.push_back(defaultDates.get());
    
    // now merge
    DateTimeArray mergedDates;
    DateTime::merge(dateArraysToMerge, mergedDates);

    //and return
    return DateTimeArraySP(new DateTimeArray(mergedDates));
}


//// Convert losses into Act/365F, Annual rates
CashFlowArraySP EffectiveCurve::asZeroRates() const
{
    Actual365F act365f;

    CashFlowArraySP cleanSpreads = defaultRates->convertToSpot(valueDate, act365f);

    //insert a 0 cashflow at valueDate to retain results
    CashFlow cf(valueDate,0);
    cleanSpreads->insert(cleanSpreads->begin(),cf);

    return cleanSpreads;
}

//----------------------------
// IDiscountCurveRisky methods
//----------------------------
/**Returns the probability that there are no default events from d1 to d2,
    conditional that there are no default events up to d1*/
double EffectiveCurve::survivalProb(const DateTime& d1,
                                    const DateTime& d2) const
{
    return defaultRates->calcDefaultPV(d1, d2);
}

/**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
double EffectiveCurve::survivalProb(const DateTime& dt) const
{
    static const string method("EffectiveCurve::survivalProb");
/*
    if (interpolationStyle == FLAT_FORWARD)
    {
        
        // build apropriate data types and use ZeroPrice
        //FlatFwdZeroCurve z(valueDate,
        //                   valueDate,
        //                   timeline->begin()+1, timeline->end(),
        //                   expLoss->begin()+1,  expLoss->end());
        //return z.zeroPrice(dt);
        
        return defaultRates->calcDefaultPV(valueDate,dt);

    }
    else if (interpolationStyle == LINEAR)
    {
        return CCMPriceUtil::linearLossInterpolate(dt,
                                                   *(timeline.get()),
                                                   *(expLoss.get()));
    }
    else
    {
        //shouldn't happen since validated on construction
        throw ModelException(method, "Unsupported loss interpolation style");
    }
*/
    //default rates now handles linear loss interpolation
    return defaultRates->calcDefaultPV(valueDate,dt);
}

/**Returns the market recovery rate on defaulted assets given default at
    time defaultDate. This allows the possibility of a term-structure of recovery rates. */
double EffectiveCurve::getRecovery(const DateTime& defaultDate) const
{
    static const string method("EffectiveCurve::getRecovery");
    throw ModelException(method, "Not yet implemented.");
}

/**Recovery-rate for immediate default.*/
double EffectiveCurve::getRecovery() const
{
    static const string method("EffectiveCurve::getRecovery");
    throw ModelException(method, "Not yet implemented.");
}

/**Returns the value at paymentDate (and conditional on no default before then)
    * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
    * in case of default between startDate and endDate, and zero otherwise.
    * The default payment is made with a delay of recoveryDelay calendar days.
    * Whether integration is done continuously or discretely and what,
    * if any approximations,
    * are used is up to the curve implementation. */
double EffectiveCurve::protectionPV(const DateTime&     paymentDate,
                                    const DateTime&     startDt,
                                    const DateTime&     endDt,
                                    RecoveryType        recTyp,
                                    double              recoveryDelay) const //number of delay days
{
    static const string method("EffectiveCurve::protectionPV");

    double pv = 0.;

    if (interpolationStyle == LINEAR)
    {
        double Ts = Actual365F::yearFraction(valueDate, startDt);
        double Te = Actual365F::yearFraction(valueDate, endDt); 
        double d  = recoveryDelay/365.;
        int numDates = timeline->size(); // for ease
        DoubleArray tL(numDates);
        int idx;
        for (idx = 0; idx < numDates; ++idx){
            tL[idx] = Actual365F::yearFraction(valueDate, (*timeline)[idx]);
        }

        //here the ZR array are just the values from the FlatFwdZeroCurve
        //...ie continuously compounded rates, converted to discount factors
        //and the tR are year fractions from the value date
        DoubleArray tR;
        DoubleArray ZR;
        ffDiscount->toCcmCurve(tR, ZR);

        // now integrate
        pv = CCMPriceUtil::contingentPrice(Ts, Te, d,
                                           tR, ZR,
                                           tL, *(expLoss.get()),
                                           interpStyle);
    }

    if (interpolationStyle == FLAT_FORWARD)
    {
        // use the CDSPricer Integration via the default rates
        DayCountConventionConstSP dcc = DayCountConventionConstSP(new Actual365F());
        IDecretionCurveConstSP prepaySP = 
            IDecretionCurveConstSP(
                new NullDecretionCurve("NullPrepayCurve"));

        // get the integration dates
        DateTimeArrayConstSP integrationDates = 
            buildIntegrationDates(startDt,
                                  endDt,
                                  recoveryDelay,
                                  discount,
                                  defaultRates->getDefaultDates());

        //now build the integrator
        CDSPricer::Integrator integrator(
            valueDate,
            integrationDates,
            ffDiscount, //discount,
            defaultRates,
            prepaySP,
            dcc,
            CDSPricer::DELAYED_DISCOUNTING,
            recoveryDelay);

        pv = 0.0;

        if (!(integrator.empty()))
        {
            while (integrator.nextStep())
            {
                pv += integrator.integralOfPVTimesProbDensity();
            }
        }
    }

    return pv;
}

/**Returns the value at paymentDate (and conditional on no default before then)
    * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
    * in case of default between startDate and endDate, and zero otherwise.
    * The default payment is made on recoveryDate, no matter when the default happens.
    * Whether integration is done continuously or discretely and what,
    * if any approximations,
    * are used is up to the curve implementation. */
double EffectiveCurve::protectionPV(const DateTime&     paymentDate,
                                    const DateTime&     startDt,
                                    const DateTime&     endDt,
                                    RecoveryType        recTyp,
                                    const DateTime&     recoveryDate) const
{
    static const string method("EffectiveCurve::protectionPV");
    throw ModelException(method, "Not yet implemented.");
}

/**Returns the value at paymentDate (and conditional on no default before then)
    * of a sequence of payments, with simple linear accrued-interest
    * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
    * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for
    * protectionPV. This allows the curve to compute default-accrual PVs itself
    * in a way which reflects the type of curve (e.g. flat forwards, etc.)
    * The accrual periods are the intervals between consecutive payment dates. To have a first
    * accrual period, you should set the first payment to zero, and then the first date is
    * the beginning of the first accrual period.
    * Payments before paymentDate are ignored, except to the extent that they affect accrued
    * interest due.
    * Any default payments are made with a delay of recoveryDelay calendar days after default.
    * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
    * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
    * this method should return the same value as pv(payments, paymentDate), inherited
    * from IDiscountCurve. */
double EffectiveCurve::annuityPV(const CashFlowArray&    payments,
                                 const DateTime&         paymentDate,
                                 RecoveryType            accruedRecTyp,
                                 double                  recoveryDelay,
                                 DateTime                accrueStartDate) const
{
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual365F());
    CashFlowArray*              paymentsCopy = new CashFlowArray(payments);
    CashFlowArraySP             paymentsSP = CashFlowArraySP(paymentsCopy);
    double                      protectionValue = .0;
    double                      accrualValue = .0;
    double                      couponValue = .0;

    if (payments.empty())
    {
        return .0;
    }
    
    if (accrueStartDate.empty())
      {
	accrueStartDate = payments[0].date;
      }
    
    if (accruedRecTyp != IDiscountCurveRisky::RECOVER_0)
    {
		YieldCurveSP myYC = YieldCurveSP::dynamicCast(discount);
        CDSPricer   cdsPricer 
            = CDSPricer(paymentsSP,                // paymentDates,
                        defaultRates,			  // defaultRates,
                        1.0,                    // notional,
                        0.0,                    // recovery,
                        true,                   // payAccruedFeed,
                        dcc,                    // swpAccrualDCC
                        paymentDate,            // valueDate,
                        accrueStartDate,            // protectionStartDate,
                        payments[payments.size()-1].date,            // protectionEndDate,
                        accrueStartDate,            // accruedStartDate,
						myYC);
        
        if (recoveryDelay == .0)
        {
            cdsPricer.calculateDefaultPayment(  true, 
                                                protectionValue, 
                                                accrualValue,
                                                CDSPricer::NO_DELAY_DISCOUNTING);
        }
        else
        {
            cdsPricer.calculateDefaultPayment(  true, 
                                                protectionValue, 
                                                accrualValue,
                                                CDSPricer::DELAYED_DISCOUNTING,
                                                recoveryDelay);
        }

        if (accruedRecTyp == IDiscountCurveRisky::RECOVER_R)
        {
            accrualValue *= getRecovery();
        }
        else if (accruedRecTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
        {
            accrualValue *= (1-getRecovery());
        }
    }
    
    couponValue = pv(payments, paymentDate);

    return couponValue + accrualValue;
}

/**Returns the value at paymentDate (and conditional on no default before then)
    * of a sequence of payments, with simple linear accrued-interest
    * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
    * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for
    * protectionPV. This allows the curve to compute default-accrual PVs itself
    * in a way which reflects the type of curve (e.g. flat forwards, etc.)
    * The accrual periods are the intervals between consecutive payment dates. To have a first
    * accrual period, you should set the first payment to zero, and then the first date is
    * the beginning of the first accrual period.
    * Payments before paymentDate are ignored, except to the extent that they affect accrued
    * interest due.
    * Any default payments are made on recoveryDate, no matter when the default happens.
    * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
    * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
    * this method should return the same value as pv(payments, paymentDate), inherited
    * from IDiscountCurve. */
double EffectiveCurve::annuityPV(
                        const CashFlowArray&    payments,
                        const DateTime&         paymentDate,
                        RecoveryType            accruedRecTyp,
                        const DateTime&         recoveryDate,
			DateTime                accrueStartDate) const 
{
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual365F());
    CashFlowArray*              paymentsCopy = new CashFlowArray(payments);
    CashFlowArraySP             paymentsSP = CashFlowArraySP(paymentsCopy);
    double                      protectionValue = .0;
    double                      accrualValue = .0;
    double                      couponValue = .0;

    if (payments.empty())
    {
        return .0;
    }
    
    if (accrueStartDate.empty())
    {
        accrueStartDate = payments[0].date;
    }
    if (accruedRecTyp != IDiscountCurveRisky::RECOVER_0)
    {
		YieldCurveSP myYC = YieldCurveSP::dynamicCast(discount);
        CDSPricer cdsPricer = CDSPricer(paymentsSP,                // paymentDates,
								        defaultRates,			 // defaultRates,
				                        1.0,                    // notional,
								        0.0,                    // recovery,
				                        true,                   // payAccruedFeed,
								        dcc,                    // swpAccrualDCC
				                        paymentDate,            // valueDate,
							accrueStartDate,            // protectionStartDate,
				                        payments[payments.size()-1].date,            // protectionEndDate,
							accrueStartDate,            // accruedStartDate,
				                        myYC);
        
        cdsPricer.calculateDefaultPayment(  true, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::NO_TIME_DISCOUNTING);
        
        
        if (accruedRecTyp == IDiscountCurveRisky::RECOVER_R)
        {
            accrualValue *= getRecovery();
        }
        else if (accruedRecTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
        {
            accrualValue *= (1-getRecovery());
        }
    }
    
    couponValue = pv(payments, paymentDate);

    return couponValue + accrualValue * discount->pv(paymentDate, recoveryDate);
}


//------------------------------------------------------
// Additional Riskless PV methods for flexibility
//------------------------------------------------------

/** Compute discount factor between two dates. This is the price
    * at time date1 of a zero-coupon (discount) bond that pays 1
    * at time date2 if it still exists then.
    * Note that, if this bond is risky, this method
    * in IDiscountCurveRisky means the PV of such a bond which knocks
    * out on default with no recovery at all, and the price is
    * contingent on no default before date1.
    * @param date1 payment settlement date (conditional on existence at date1)
    * @param date2 payment/zero-coupon bond maturity date
    * @return Discount factor between date1 & date2
    */
double EffectiveCurve::risklessPV(const DateTime& date1,
                                  const DateTime& date2) const
{
    return discount->pv(date1, date2);
}

/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double EffectiveCurve::risklessPV(const DateTime& date) const
{
    //return discount->pv(date);
    return ffDiscount->zeroPrice(date);
}

/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double EffectiveCurve::risklessPV(const CashFlowArray& cashFlows,
                                  const DateTime&      baseDate) const
{
    return discount->pv(cashFlows, baseDate);
}

//-----------------------
// IDiscountCurve methods
//-----------------------

/** @return Discounting curve's currency - typically this is the
    ISO code eg "USD" (although it is up to the client - you may
    assume however that yield curves in the same currency return the
    same value in getCcy(). Note it is NOT the name of the yield
    curve (which might be eg "GBP-LIBOR"). */
string EffectiveCurve::getCcy() const
{
    return discount->getCcy();
}

/** Compute discount factor between two dates. This is the price
    * at time date1 of a zero-coupon (discount) bond that pays 1
    * at time date2 if it still exists then.
    * Note that, if this bond is risky, this method
    * in IDiscountCurveRisky means the PV of such a bond which knocks
    * out on default with no recovery at all, and the price is
    * contingent on no default before date1.
    * @param date1 payment settlement date (conditional on existence at date1)
    * @param date2 payment/zero-coupon bond maturity date
    * @return Discount factor between date1 & date2
    */
double EffectiveCurve::pv(const DateTime& date1,
                          const DateTime& date2) const
{
    return defaultRates->calcTotalDefaultPV(date1, date2) * discount->pv(date1, date2);
}

/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double EffectiveCurve::pv(const DateTime& date) const
{
    static const string method("EffectiveCurve::pv");

    double z1 = ffDiscount->zeroPrice(date);
    double z2 = defaultRates->calcDefaultPV(valueDate,date);
    return z1 * z2;
}

/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double EffectiveCurve::pv(const CashFlowArray& cashFlows,
                          const DateTime&      baseDate) const
{
    double value = .0;
    for (int i=0; i < cashFlows.size(); ++i)
    {
       const DateTime cfDate = cashFlows[i].date;
		if (cfDate.isGreater(baseDate))
                    value += cashFlows[i].amount * pv(baseDate, cfDate);
    }
    return value;
}

/** return the bootstrapped dates */
DateTimeArray EffectiveCurve::zeroDates() const
{
    //combine discount and default dates
    DateTimeArraySP mergedDates = 
        mergeDates(0.0, discount, defaultRates->getDefaultDates());

    return *(mergedDates.get());
}

// accessor methods for logOfDiscFactorKey
IDiscountCurve::IKey* EffectiveCurve::getDiscountKey() const
{
    return discount->logOfDiscFactorKey();
}

DefaultRates::IKey* EffectiveCurve::getRiskyKey() const
{
    return defaultRates->logOfDefaultPVKey();
}

/** Returns a key used to optimise repeated calculations of discount
    factors/forward rate. The calc method for this key returns the 
    natural logarithm of the discount factor (or equivalently the
    product of the forward rate (continuous, Act/365F) and the negative
    year fraction (Act/365F) betweeen the two dates.
    The default implementation has no performance improvements. */
//this is questionable !
IDiscountCurve::IKey* EffectiveCurve::logOfDiscFactorKey() const
{
    return new RiskyLogOfDiscFactorKey(this);
}

/** Returns a DefaultRates object which gives access to 
    useful functionality including "default rates", aka clean default
    spreads. The information needed to bootstrap the clean spreads is
    obtained from this object.
    The DefaultRate returned is immutable and therefore it will not be
    changed - which means that there is no need to clone it */
DefaultRatesSP EffectiveCurve::getDefaultRates() const
{
    return defaultRates;
}


//-----------------------------------------------
// IGetMarket methods - exposed by IDiscountCurve
//-----------------------------------------------

//// populate from market cache
void EffectiveCurve::getMarket(const IModel* model, const MarketData* market)
{
    //nothing to do;
}

//-------------------------------
// EffectiveCurve private methods
//-------------------------------

IObject* EffectiveCurve::defaultConstructor()
{
    return new EffectiveCurve();
}


void EffectiveCurve::load(CClassSP& clazz)
{
    REGISTER(EffectiveCurve, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IDiscountCurveRisky);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(valueDate,          "Reference date");
    FIELD(interpolationStyle, "How the curve should be interpolated");
    FIELD       (timeline,           "The curve dates");
    FIELD       (expLoss,            "The curve data");

    //FIELD(interpStyle, "a slightly different form of interpolationStyle");
    //FIELD       (defaultRates,"re-use existing default rates class [essentially just dates & rates]");
    FIELD       (discount,    "may be risk-free or risky");
    FIELD       (ffDiscount,  "for backwards compatability");

    //FIELD_MAKE_TRANSIENT(interpStyle);
    //FIELD_MAKE_TRANSIENT(defaultRates);
    FIELD_MAKE_TRANSIENT(discount);
    FIELD_MAKE_TRANSIENT(ffDiscount);

}

CClassConstSP const FlatFwdZeroCurve::TYPE = CClass::registerClassLoadMethod(
    "FlatFwdZeroCurve", typeid(FlatFwdZeroCurve), load);

CClassConstSP const EffectiveCurve::TYPE = CClass::registerClassLoadMethod(
    "EffectiveCurve", typeid(EffectiveCurve), load);

DRLIB_END_NAMESPACE
