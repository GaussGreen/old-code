//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : UntweakableCDSParSpread.hpp
//
//   Description : A CDSParSpread that is built from clean default spreads
//                 Motivation: for porting SRM3 tests
//                 DO NOT USE THIS FOR ANYTHING ELSE
//
//   Author      : Mark A Robson
//
//   Date        : 25 August 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/CDSParSpreadsBase.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/ICDSParSpreadsCache.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"


DRLIB_BEGIN_NAMESPACE

class CParSpreadDefaultRates2: public DefaultRates {
public:
    // constructor
    CParSpreadDefaultRates2(const DateTime&       valueDate,
                            const DateTimeArray&  benchmarkDates, 
                            const CDoubleArray&   annualizedCleanSpreads):
        valueDate(valueDate), benchmarkDates(benchmarkDates),
        annualizedCleanSpreads(annualizedCleanSpreads)
    {}

    double calcDefaultPV(const DateTime& fromDate,
                         const DateTime& toDate) const 
    {
        double pvFactor = 0.0;
        int idx = 1;
        double nbDates = benchmarkDates.size();
        while(idx < nbDates - 1 && toDate > benchmarkDates[idx]) {
            idx++;
        }
        idx -= 1;

        double T = (double)(toDate.getDate() - fromDate.getDate())/365.0;
        if (fromDate.getDate() == benchmarkDates[idx].getDate()) {
            pvFactor = 1.0;
            return pvFactor;
        } else 
        if (toDate.getDate() == benchmarkDates[idx].getDate()) {
            double r = annualizedCleanSpreads[idx];
            pvFactor = pow(1.0 + r, -T);
            return pvFactor;
        } else {
           
            double denominator = (double)(benchmarkDates[idx+1].getDate() - benchmarkDates[idx].getDate()); 
            double weight1 = (double)(benchmarkDates[idx+1].getDate() - toDate.getDate()) / denominator;
            double weight2 = (double)(toDate.getDate() - benchmarkDates[idx].getDate()) / denominator;
            double r = weight1 * annualizedCleanSpreads[idx] + weight2 * annualizedCleanSpreads[idx+1];
            pvFactor = pow(1.0 + r, -T);
            return pvFactor; 
        }    
    }
    CashFlowArraySP getCleanSpreadCurve() const {
        int nb = annualizedCleanSpreads.size();
        CashFlowArraySP cleanSpreadsWithDates(new CashFlowArray(nb));
        for (int i = 0; i < nb; i++) {
            (*cleanSpreadsWithDates)[i].date = benchmarkDates[i];
            (*cleanSpreadsWithDates)[i].amount = annualizedCleanSpreads[i];
        }
        return cleanSpreadsWithDates;
    }
    DateTimeArraySP getDefaultDates() const {
        return DateTimeArraySP(new DateTimeArray(benchmarkDates));
    }
    const CDoubleArray& getRates() const {
        return annualizedCleanSpreads;
    }
    const DateTimeArray& getDates() const {
        return benchmarkDates;
    }
    const DateTime& getValueDate() const {
        return valueDate;
    }

private: 
    DateTime        valueDate;
    DateTimeArray   benchmarkDates;
    DoubleArray     annualizedCleanSpreads;

};

class UntweakableCDSParSpread: public MarketObject,
                               public virtual ICDSParSpreads,
                               public virtual Duration::IParHandlerWithClosedForm,
                               public virtual ITweakableWithRespectTo<ParSpreadParallel>,
                               public virtual ITweakableWithRespectTo<ParSpreadPointwise>,
                               public virtual ParSpreadLevel::IShift,
                               public virtual ParSpreadParallelShift::IShift,
                               public virtual ParSpreadPropShift::IShift,
                               public virtual Theta::Shift
{
    //  fields ///
    string                 name;
    DateTime               valueDate;
    DateTimeArray          dates;
    DoubleArray            cleanSpreads; // continuous spot rates
    int                    parSwapFreq;
    DayCountConventionSP   parDCC;
    double                 parRecovery;
    bool                   parAccrueFee; /* recovery rate convention: true means
                                            you get accrual + notional */
    int                    spotOffset;
    YieldCurveWrapper      discount; // corresponding discount curve
    bool                   hasDefaulted; // has a default happened
    DateTime               defaultDate;  // if so, when
    ICDSVolWrapper         cdsVol;   /* holds the cds vol, optional */
    bool                   useDiscreteCompounding;
    
public:
    static CClassConstSP const TYPE;
    
    virtual ~UntweakableCDSParSpread(){}

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const{
        return name;
    }

    virtual string getParSpreadsName() const{
        throwError(NULL);
        return "";
    }

    void validatePop2Object() {
        // Initialise the cache here - could be done in the constructor
        // but it would be created for clones also (just to be released
        // almost inmediateley in the clone() function), so do it here 
        cache = ICDSParSpreadsCacheSP(new ICDSParSpreadsCache);
    }

    IObject* clone() const{
        UntweakableCDSParSpread* copy = DYNAMIC_CAST(
               UntweakableCDSParSpread, MarketObject::clone());
        copy->cache = cache;
        return copy;
    }

    virtual void getMarket(const IModel* model, const MarketData* market){
        market->GetReferenceDate(valueDate);
        if (!discount.isEmpty()){
            // this should be mandatory
            discount.getData(model, market);
        }
        if (!cdsVol.isEmpty()){
            // NB Model might decline to get a vol
            cdsVol.getData(model, market);
        }
    }

    /** @return Yield curve's currency [isocode] */
    virtual string getCcy() const{
        if (discount.isEmpty()){
            throw ModelException("UntweakableCDSParSpread::getCcy", 
                                 "No discount curve specified for "+name);
        }
        return discount->getCcy();
    }

    /** return CDS Par Spreads discount yield curve name (not the isocode) */
    virtual string getYCName() const{
        return discount->getName();
    }

    /** Throws an exception*/
    DoubleArrayConstSP getParSpreads() const {
        throwError(NULL);
        return *(new DoubleArrayConstSP());
    }

    void bucketShift(const DoubleArray& shifts)
    {
        throwError(NULL);
    }

    /** Returns the expiries (relative: eg. 1M) */
    ExpiryArrayConstSP getParSpreadsExpiries() const {
        throwError(NULL);
        return ExpiryArrayConstSP();
    }

    bool isFeeAccrued() const {
        throwError(NULL);
        return false;
    }

    int getSwapFrequency() const {
        throwError(NULL);
        return 0;
    }

    /** Returns a DefaultRates object which gives access to useful
        functionality including "default rates", aka clean default
        spreads */
    DefaultRatesSP defaultRates(YieldCurveConstSP  discount,
                                BadDayConventionSP bdc) const
    {
        if(useDiscreteCompounding) {
            return DefaultRatesSP(new CParSpreadDefaultRates2(valueDate, 
                                                              dates, 
                                                              cleanSpreads));
        } else {
            // think we need to take cleanSpreads and compute forward spreads
            Actual365F act365F;
            DoubleArray fwdSpread(dates.size());
            for (int i = 0; i < dates.size(); i++){
                try{
                    if (i == 0){
                        fwdSpread[0] = 
                            RateConversion::rateConvert(valueDate, dates[0],
                                                        cleanSpreads[0],
                                                        &act365F,
                                                        CompoundBasis::ANNUAL,
                                                        &act365F, 
                                                        CompoundBasis::CONTINUOUS);
                    } else {
                        // convert from annual to continuous
                        double ctsRate1 =
                            RateConversion::rateConvert(valueDate, dates[i-1],
                                                        cleanSpreads[i-1], &act365F,
                                                        CompoundBasis::ANNUAL,
                                                        &act365F, 
                                                        CompoundBasis::CONTINUOUS);
                        double ctsRate2 =
                            RateConversion::rateConvert(valueDate, dates[i],
                                                        cleanSpreads[i], &act365F,
                                                        CompoundBasis::ANNUAL,
                                                        &act365F, 
                                                        CompoundBasis::CONTINUOUS);
                        double yf1 = Actual365F::yearFraction(valueDate,dates[i-1]);
                        double yf2 = Actual365F::yearFraction(valueDate, dates[i]);
                        double fwdRate = ctsRate2*yf2 - ctsRate1*yf1;
                        fwdRate /= (yf2-yf1);
                        fwdSpread[i] = fwdRate;
                    }
                } catch (exception& e) {
                    throw ModelException(e, 
                                         "UntweakableCDSParSpread::defaultRates",
                                         "For curve with name "+name+" at "+
                                         dates[i].toString());
                }
            }   
            // now create the DefaultRatesSP object
            return DefaultRatesSP(new CDSHelper::CParSpreadDefaultRates(valueDate,
                                                                        dates,
                                                                        fwdSpread));
        }
    }

    /** Same as above but uses day count convention/bad
        day convention from this object (or defaults it as
        appropriate) similarly for yield curve */
    virtual DefaultRatesSP defaultRates() const {
        return getCachedDefaultRates(this, UNTWEAKABLE_CDS_PAR_SPREAD);
    }

        /** return the bootstrapped dates */
    virtual DateTimeArray zeroDates() const
    {
        DefaultRatesSP dR = defaultRates();
        DateTimeArraySP dates = dR->getDefaultDates();
        return *(dates.get());
    }

    // accessor methods for logOfDiscFactorKey
    virtual IDiscountCurve::IKey* getDiscountKey() const
    {
        return discount->logOfDiscFactorKey();
    }

    virtual DefaultRates::IKey* getRiskyKey() const
    {
        DefaultRatesSP dR = defaultRates();
        return dR->logOfDefaultPVKey();
    }

    /** Returns a key used to optimise repeated calculations of discount
        factors/forward rate. The calc method for this key returns the 
        natural logarithm of the discount factor (or equivalently the
        product of the forward rate (continuous, Act/365F) and the negative
        year fraction (Act/365F) betweeen the two dates.
        The default implementation has no performance improvements. */
    virtual IDiscountCurve::IKey* logOfDiscFactorKey() const
    {
        return new RiskyLogOfDiscFactorKey(this);
    }

    /** Returns the recovery rate */
    virtual double getRecovery() const{
        return parRecovery;
    }

    /** Returns true if the name has defaulted */
    virtual bool defaulted() const{
        return hasDefaulted;
    }
    
    /** Returns the date that this name defaulted. If a default has not
        happened an 'empty' DateTime is returned */
    virtual const DateTime& getDefaultDate() const{
        return defaultDate;
    }

    /** Returns a processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const
    {
        if (cdsVol.isEmpty()){
            throw ModelException("UntweakableCDSParSpread::getProcessedVol",
                                 "Pricing model requires CDS Vol be supplied for "
                                 "UntweakableCDSParSpread curve with name " + getName());
        }
        return cdsVol->getProcessedVol(volRequest, this);
    }

    /** Returns the same as getName() */
    string getCorrelationName() const{
        return getName();
    }

    virtual bool sensShift(Theta* shift){
        DateTime newDate(shift->rollDate(valueDate));
        if (!newDate.equals(valueDate)){
            throwError(shift);
        }
        return true;
    }

    virtual string sensName(const ParSpreadParallel*) const{
        return name;
    }

    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadParallel>&){
        throwError(NULL);
        return TweakOutcome(0, false);
    }

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const ParSpreadPointwise*) const{
        return name;
    }
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadPointwise>&){
        throwError(NULL);
        return TweakOutcome(0, false);
    }
    virtual ExpiryWindowArrayConstSP sensQualifiers(const ParSpreadPointwise* shift)
            const{
        throwError(NULL);
        return ExpiryWindowArrayConstSP();
    }

    // Proportional shift support
    virtual string sensName(ParSpreadPropShift* shift) const{
        return name;
    }
    virtual bool sensShift(ParSpreadPropShift* shift){
        throwError(shift);
        return false;
    }

    // Parallel shift (scenario) support
    virtual string sensName(ParSpreadParallelShift* shift) const{
        return name;
    }
    virtual bool sensShift(ParSpreadParallelShift* shift){
        throwError(shift);
        return false;
    }

    // Flat override support
    virtual string sensName(ParSpreadLevel* shift) const{
        return name;
    }
    virtual bool sensShift(ParSpreadLevel* shift){
        throwError(shift);
        return false;
    }
    
    /** Returns the date when to stop tweaking after lastDate */
    virtual const DateTime stopTweaking(const DateTime& lastDate) const {
        throw ModelException("UntweakableCDSParSpread::stopTweaking", 
                             "You really can't shift these CDSParSpreads!");
    }

    /** Calculates the maximum tweak sizes for this ICDSParSpreads */
    virtual DoubleArrayConstSP calculateTweakSizes() const
    {
        throw ModelException("UntweakableCDSParSpread::calculateTweakSizes", 
                             "You really can't shift these CDSParSpreads!");
    }
    
    /** Returns the day count convention used for par CDS */
    virtual DayCountConvention* dayCountConv() const {
        return parDCC.get();
    }

    HolidayConstSP getHolidays() const {
        throwError(0);
        return HolidayConstSP();
    }

    BadDayConventionConstSP getBadDayConvention() const {
        throwError(0);
        return BadDayConventionSP();
    }
    
    virtual DateTime spotDate(
        const DateTime& valueDate) const
    {
        return valueDate.rollDate(spotOffset);  
    }
    
    
#ifdef CDS_BACKWARD_COMPATIBILITY
    /** Set the bad day convention */
    virtual void setBadDayConvention(BadDayConventionSP bdc) {
        // nothing to do
    }

    /** Set the holidays */
    virtual void setHolidays(HolidaySP holidays) {
        // nothing to do
    }
#endif

    /** ICreditCurve implementation */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate,
                                    const BadDayConvention* bdc,
                                    const Holiday* hols) const {
        throw ModelException(
            "UntweakableCDSParSpread::getCurrentSpread", 
            "Not implemented");
    }

    /** ICreditCurve implementation */
    virtual double getCurrentSpread(const DateTime& valueDate,
                                    const DateTime& maturityDate) const {
        throw ModelException(
            "UntweakableCDSParSpread::getCurrentSpread", 
            "Not implemented");
    }

    /** IDiscountCurveRisky partial implementation */
    virtual double pv(const DateTime& date1, 
                      const DateTime& date2) const { throw ModelException("not implemented yet"); }
    
    virtual double pv(const DateTime& date) const { return -9999; }
    
    virtual double pv(const CashFlowArray& cashFlows,
                      const DateTime&      baseDate) const { throw ModelException("not implemented yet"); }



    // -----------------------------------------
    //  IBadDayAdjuster method
    //------------------------------------------
    /** Returns date bad day adjusted using the bad day convention
     * and holidays in this object */
    virtual DateTime badDayAdjust(const DateTime& date) const {
        // Return the original date with no adjustment
        return date;
    }

    /** Add a number of business days to a date */
    virtual DateTime addBusinessDays(const DateTime& from, int busDays) const {
        // Roll the "from" date, with no business days adjustments
        return from.rollDate(busDays);
    }


    //------------------------------
    // Cache-related methods
    //------------------------------

    /** Returns a DefaultRates object from the cache.
     * This method collects from the cache the default rates associated to the
     * ICDSParSpreads passed as parameter (calculating them if not there in
     * the first place). 
     * It is typically invoked from the "defaultRates" method if required (there
     * can be another layer of caching there) and "external classes" are 
     * strongly DISCOURAGED from using this method directly: the "defaultRates" 
     * method should be used instead */
    const DefaultRatesSP getCachedDefaultRates(
        const ICDSParSpreads* cds,
        const TypeOfEntry     entryType) const 
   {
        static const string method = 
            "UntweakableCDSParSpread::getCachedDefaultRates";

        if (!cache.get()) {
            throw ModelException(method, 
                                 "Cache is NULL trying to obtain default rates");
        }
        return cache->getEntry(cds, entryType);
    }
    
    /** Provides a facility for previously constructed DefaultRates
     *  objects to be added to the cache */
    const bool cacheDefaultRates(const ICDSParSpreads* cds,
                                 const TypeOfEntry     entryType,
                                 const DefaultRatesSP  entry)
    {
        static const string method = 
            "UntweakableCDSParSpread::getCachedDefaultRates";

        if (!cache.get()) {
            throw ModelException(method, 
                                 "Cache is NULL trying to obtain default rates");
        }
        return cache->addEntry(cds, entryType, entry);
    }

    /** Returns a DefaultRates object.
     * This method does NOT use the (potentially availabe) caching mechanism
     * which avoids the slow process of calculating default rates everytime.
     * It is invoked from inside the cache to calculate the default rates
     * for the first time and its direct use from any external classes is
     * strongly DISCOURAGED: the "defaultRates" method should be used instead.
     * Here we are in UntweakableCDSParSpread so the way to compute the clean 
     * spreads is bootstrapping... so do that */
    const DefaultRatesSP computeDefaultRatesForCache() const {
        // Note this is being done the other way around in this class... 
        // See the "Quanto caching" document in the QLib DB
        if (!discount){
            throw ModelException("UntweakableCDSParSpread::computeDefaultRatesForCache",
                                 "Discount curve not specified for "+name);
        }
        return defaultRates(discount.getSP(),
                            BadDayConventionSP(new BadDayNone()));
    }
    
    /** Resets the pointer to the cache - The cache will not be used
     * in this object from this point onwards:
     * This is invoked for the CDSParSpreadsBase clone stored inside the cache, 
     * and is required so that the cache is only pointed at from external 
     * CDSParSpreadsBase. Otherwise the cache object would never get freed */
    void resetCache() const {
        cache.reset();
    }

    /** Hash code function used for caching - the cache needs improved
     * performance compared with the default "hashCode" function in CObjet: only
     * the required components are hashed (see comment in equalToOpt below) */
    int hashCodeOpt() const {
        // It is easier to compare for equality, so return always the same
        // value (arbitrarily, 1) and do the actual comparison in equaltoOpt
        return 1;
    }


    /** Comparison function used for caching - the cache needs improved
     * performance compared with the default "equalTo" function in CObjet: only
     * the required components are compared. This means that the objects being 
     * compared are not required to be identical: only equal as far as the 
     * default rates calculation is concerned (e.g, if the volType is not used 
     * in the calculation, the equalToOpt method should not fail if the volTypes
     * of both objects are different) */
    bool equalToOpt(const ICDSParSpreads* cds) const {
        // Since we cannot tweak this curve, comparing if this object equals 
        // cds is enough.
        return (this == cds);
    }

    
    //------------------------------
    // Duration::IParHandler methods
    //------------------------------
    //return benchmarks for par instruments
    //these will determine the durations to calculate unless specified
    ExpiryArrayConstSP getParBenchmarks() const {
        return getParSpreadsExpiries();
    }

    //return closed form solution
    ExpiryResultArraySP getDuration(const Duration* duration) const {
        return CDSParSpreadsBase::CDSGetDuration(duration, 
                                                 getParBenchmarks(),
                                                 getYCName(),
                                                 getName());
    }

private:
    static IObject* defaultConstructor(){
        return new UntweakableCDSParSpread();
    }

    static void throwError(IObject* shift){
        throw ModelException("UntweakableCDSParSpread",
                             "You really can't shift these CDSParSpreads!");
    }

    UntweakableCDSParSpread(): MarketObject(TYPE),
        useDiscreteCompounding(false) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(UntweakableCDSParSpread, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(ICDSParSpreads);
        IMPLEMENTS(ITweakableWithRespectTo<ParSpreadParallel>);
        IMPLEMENTS(Duration::IParHandlerWithClosedForm);
        IMPLEMENTS(ITweakableWithRespectTo<ParSpreadPointwise>);
        IMPLEMENTS(ParSpreadLevel::IShift);
        IMPLEMENTS(ParSpreadParallelShift::IShift);
        IMPLEMENTS(ParSpreadPropShift::IShift);
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "Identifies this object");
        FIELD(valueDate, "today");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(dates, "Dates for clean spreads");
        FIELD(cleanSpreads, "The [spot] clean spreads");
        FIELD(parSwapFreq, "swap freq for par spreads (1, 2, 4)");
        FIELD(parDCC, "day-count-convention");
        FIELD(parRecovery, "recovery rate for bonds underlying"
                     " par swaps");
        FIELD(parAccrueFee, "make accrued payment on default?");
        FIELD(spotOffset, "gives effective date of par curve");
        FIELD(discount, "Discount yield curve");
        FIELD(hasDefaulted, "has a default happened");
        FIELD(defaultDate, "If a default has happened, when");
        FIELD_MAKE_OPTIONAL(defaultDate);
        FIELD(cdsVol, "Volatility of associated CDS Swaptions");
        FIELD_MAKE_OPTIONAL(cdsVol);
        FIELD(useDiscreteCompounding, "whether to use discrete or continuous compounding");
        FIELD_MAKE_OPTIONAL(useDiscreteCompounding);
    }

    mutable ICDSParSpreadsCacheSP cache;     // <ICDSParSpreads, cleanSpreads> cache $unregistered
};

CClassConstSP const UntweakableCDSParSpread::TYPE = 
CClass::registerClassLoadMethod(
    "UntweakableCDSParSpread", typeid(UntweakableCDSParSpread), load);


/* external symbol to allow class to be forced to be linked in */
bool UntweakableCDSParSpreadLinkIn(){
    return true;
}

DRLIB_END_NAMESPACE

