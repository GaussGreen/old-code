//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : CDSParSpreadsBase.hpp
//
//   Description : Base class for CDS par spreads
//
//   Date        : August 30, 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/CDSPricer.hpp"
#include "edginc/ParCDS.hpp"
#include "edginc/ParSpreadRhoPointwise.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/BadDayNone.hpp"
#include "edginc/ICDSParSpreadsCache.hpp"
#include "edginc/NullDecretionCurve.hpp"
#include "edginc/DecretionCurve.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/CreditIndexBasis.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"
#include "edginc/CDSHelper.hpp"

DRLIB_BEGIN_NAMESPACE


/*******************************************************************************
 *                                CDSParSpreadsBase
 *******************************************************************************/


IObject* CDSParSpreadsBase::clone() const{
    CDSParSpreadsBase* copy = DYNAMIC_CAST(CDSParSpreadsBase, MarketObject::clone());
    copy->cleanSpreads = cleanSpreads;
    copy->cache = cache;
    return copy;
}

CDSParSpreadsBase::~CDSParSpreadsBase(){}

CDSParSpreadsBase::CDSParSpreadsBase(CClassConstSP clazz): 
    BootstrappableCDSParSpreads(clazz), 
    parRecovery(0.0), parSwapFreq(1), 
    parBDC("No value")/*parBDC("N")*/, // should be mandatory
    parAccrueFee(true), spotOffset(1), 
    priceViaCreditSpreadCurve(false) 
{}

CDSParSpreadsBase::CDSParSpreadsBase
(string                   name,
 DateTime                 valueDate,
 int                      spotOffset,
 double                   parRecovery,
 int                      parSwapFreq,
 ExpiryArraySP            expiries,
 CDoubleArraySP           spreads,
 CDoubleArraySP           upfronts,
 bool                     parAccrueFee,
 string                   parDCC,       // day count convention for par CDS
 string                   parBDC,       // bad day convention for par CDS
 HolidayConstSP           parHolSP,     // holidays for par CDS
 YieldCurveConstSP        discountSP,   // corresponding discount curve
 DecretionCurveConstSP    prepaySP      // prepay curve
 ):
    BootstrappableCDSParSpreads(TYPE),
    name(name), valueDate(valueDate), parRecovery(parRecovery),
    parSwapFreq(parSwapFreq), parDCC(parDCC), parBDC(parBDC),
    parAccrueFee(parAccrueFee), spotOffset(spotOffset),
    priceViaCreditSpreadCurve(false), creditSpreadCurve(0),csRestore(0)
{
    static const string method = "CDSParSpreadsBase::CDSParSpreadsBase";

    try
    {
        //YieldCurve
        if(!discountSP.get())
        {
            throw ModelException(method, "missing discount curve!");
        }else{
            YieldCurveSP ysp = YieldCurveSP(dynamic_cast<YieldCurve*>(discountSP->clone()));
            discount.setObject(MarketObjectSP(ysp));
        }
        
        //Holiday
        if(parHolSP.get()) 
        {
            HolidaySP holsp = HolidaySP(dynamic_cast<Holiday*>(parHolSP->clone()));
            parHols.setObject(MarketObjectSP(holsp));
        } else {
            parHols.setObject(MarketObjectSP(Holiday::weekendsOnly()));
        }
        
        //prepay
        if(!prepaySP.get())
        {
            prepay.setObject(MarketObjectSP(NullDecretionCurveSP(new NullDecretionCurve(name))));
        }else{
            DecretionCurveSP p = DecretionCurveSP(dynamic_cast<DecretionCurve*>(prepaySP->clone()));
            prepay.setObject(MarketObjectSP(p));
        }
        
        //parSpread
        if(!spreads.get())
        {
            throw ModelException(method, "missing par spreads curve!");
        } else {
            ExpiryArraySP expi = ExpiryArraySP(new ExpiryArray(*expiries));
            DoubleArray sprds = DoubleArray(*spreads);
            DoubleArraySP upf = DoubleArraySP(0);
            if(upfronts.get())
            {
                upf = DoubleArraySP(new DoubleArray(*upfronts));
            }
            //ParSpreadCurveSP parSpreadSP = 
            //  ParSpreadCurveSP(new ParSpreadCurve(name, expi, sprds, upf));
            
            ParSpreadCurveSP parSpreadSP = 
                ParSpreadCurveSP(new ParSpreadCurve(name, expiries, *spreads, upfronts));
            
            parSpreads.setObject(MarketObjectSP(parSpreadSP));
        }
        
        parDayCountConv = DayCountConventionSP((DayCountConventionFactory::make(parDCC)));
        parBadDayConv = BadDayConventionSP((BadDayConventionFactory::make(parBDC)));

        validatePop2Object();
    }
    catch(exception& e)
    {
        throw ModelException(e, method, "failed to construct a CDSParSpreadsBase!" + name);
    }
}

void CDSParSpreadsBase::validatePop2Object(){
    static const string method = "CDSParSpreadsBase::validatePop2Object";
    try {
        if (!(parSwapFreq == 1 || parSwapFreq == 2 || parSwapFreq == 4 
              || parSwapFreq == 12)) {
            throw ModelException(method,
                                 "Par swap frequency (" + 
                                 Format::toString(parSwapFreq) + 
                                 ") must be 1, 2, 4, or 12");   
        }
        
        if (!isRecoveryValid()) {
            throw ModelException(method, 
                                 "Recovery out of [0,1] (Recovery=" +
                                 Format::toString(parRecovery) + ")");   
        }
        
        parDayCountConv.reset(DayCountConventionFactory::make(parDCC));
        if (!(parBDC == "No value")) {
            parBadDayConv.reset(BadDayConventionFactory::make(parBDC));
        }
#ifdef CDS_BACKWARD_COMPATIBILITY
        else {
            //default to None
            parBadDayConv.reset(new BadDayNone());
        }
#endif

        
        // Initialise the cache here - could be done in the constructor
        // but it would be created for clones also (just to be released
        // almost inmediateley in the clone() function), so do it here 
        cache = ICDSParSpreadsCacheSP(new ICDSParSpreadsCache);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Clears cached clean spreads */
void CDSParSpreadsBase::fieldsUpdated(const CFieldArray& fields){
    // this is called by the infrastructure when a component field has
    // been changed - reset the "local" cleanSpreads cache
    cleanSpreads.reset();
}


string CDSParSpreadsBase::getName() const {
    return name;
}

string CDSParSpreadsBase::getParSpreadsName() const {
    return parSpreads->getName();
}

DayCountConvention* CDSParSpreadsBase::dayCountConv() const {
    return parDayCountConv.get();
}

DateTime CDSParSpreadsBase::spotDate(const DateTime& valueDate) const {
    return valueDate.rollDate(spotOffset);
}

int CDSParSpreadsBase::getSwapFrequency() const {
    return parSwapFreq;
}

/** Allow setting the par spreads. Used, e.g., to allow interpolating
 * the original par spreads - note this method is protected */
void CDSParSpreadsBase::setParSpreadCurveWrapper(ParSpreadCurveWrapper curveWrapper) {
    parSpreads = curveWrapper;
    // Reset the "local" cleanSpreads cache
    cleanSpreads.reset();    
}

int CDSParSpreadsBase::getSpotOffset() const
{
    return spotOffset;
}

DoubleArrayConstSP CDSParSpreadsBase::getUpfrontFees() const
{
    return parSpreads->getUpfronts();
}

// void CDSParSpreadsBase::acceptWrapperNameCollector(CDSParSpreadsBase* parSpreads, 
//                                                    WrapperNameCollector* collector)
// {
//     collector->addName(parSpreads->getName());
// }

/** Get today and the par spreads curve from the market data cache */
void CDSParSpreadsBase::getMarket(const IModel* model, const MarketData *market){
    static const string method = "CDSParSpreadsBase::getMarket";
    market->GetReferenceDate(valueDate);
    parSpreads.getData(model, market);
    if (!discount.isEmpty()){
        // this should be mandatory
        discount.getData(model, market);
    }
    if (!parHols.isEmpty()){
        // ideally this should be mandatory
        parHols.getData(model, market);
    } else {
        // default it to weekends only
        parHols.setObject(MarketObjectSP(Holiday::weekendsOnly()));
    }

    if(!prepay.isEmpty()){
        prepay.getData(model, market);
    } else {
        //default to NullDecretionCurve
        
        prepay.setObject(MarketObjectSP
                         (NullDecretionCurveSP
                          (new NullDecretionCurve( getName() ))));
        
    }

    // if specified, ask for the vol
    if (!cdsVol.isEmpty()){
        // NB This might result in a null vol if the model reckons we don't 
        // need it
        cdsVol.getData(model, market);
    }
    
#if 0
    /** Removed this code as it's not clear why it's here and it doesn't
        work for fixed benchmarks and we reckon there's no requirement
        that the dates are "on cycle" */

    // Validates the expiries. This can't be done in validatePop2Object
    // since if the parSpreads are in the cache and not the object it will
    // fail - this routine is called by the CDS product validate
    try {
        ExpiryArrayConstSP expiries = parSpreads->getExpiries();

        // validate that all of the maturities are on cycle with the
        // swap frequency.  Use a not-so-arbitrary date as a base date
        // to test from. The on-cycle routine is pretty crappy and
        // automatically fails for dates after the 28th of the month.
        DateTime pseudoValueDate("01-JAN-2001", DateTime::START_OF_DAY);
        for (int i = 0; i< expiries->size(); i++) {
            DateTimeArray tmpSwapDates(
                SwapTool::simpleSwapDates(
                    pseudoValueDate,
                    (*expiries)[i]->toDate(pseudoValueDate),
                    parSwapFreq));
            
            if (SwapTool::onCycle(pseudoValueDate,     
                                  tmpSwapDates[tmpSwapDates.size()-1],   
                                  (tmpSwapDates.size()-1)*(4/parSwapFreq),
                                  "Q") == false) {
                throw ModelException(method, "Par swap maturities must be on"
                                     " cycle with the par swap frequency");
            }      
        }    
    } catch (exception& e) {
        throw ModelException(e, method, "For " + name + " CDSParSpreadsBase");
    }
#endif 
}

/** @return [discount] yield curve's currency [isocode] */
string CDSParSpreadsBase::getCcy() const{
    if (discount.isEmpty()){
        throw ModelException("CDSParSpreadsBase::getCcy", "No discount curve "
                             "specified");
    }
    return discount->getCcy();
}

/** return CDS Par Spreads discount yield curve name (not the isocode) */
string CDSParSpreadsBase::getYCName() const{
    return discount->getName();
}

/** bdc and hols to be removed, until then they are optional and default to
    values from this object */
DateTimeArray CDSParSpreadsBase::getExpiryDates(DateTime                today, 
                                                const BadDayConvention* bdc,
                                                const Holiday*          hols) const{
    static const string method = "CDSParSpreadsBase::getExpiryDates";
    try {
        if (!bdc){
            if (!parBadDayConv.get()) {
                throw ModelException(
                    method,
                    "Attempting to set an undefined bad day convention !");
            }
            bdc = parBadDayConv.get();
        }
        if (!hols){
            if (!parHols.get()) {
                throw ModelException(
                    method,
                    "Attempting to set undefined holidays !");
            }
            hols = parHols.get();
        }

        ExpiryArrayConstSP expiries = parSpreads->getExpiries();

        // first turn expiries into absolute dates
        int numExpiries = expiries->size();
        DateTimeArray benchmarks(numExpiries);

        for (int i = 0; i < numExpiries; i++) {
            benchmarks[i] = bdc->adjust((*expiries)[i]->toDate(today), hols);
        }
        return benchmarks;
    } catch (exception &e) {
        throw ModelException(e, method, "Failed to convert offsets to dates");
    }
}

ExpiryArrayConstSP CDSParSpreadsBase::getExpiries() const {
    return parSpreads->getExpiries();
}

int CDSParSpreadsBase::numSpreads() {
    return parSpreads->getExpiries()->size();
}


/** Allows a shift to multiple points on the curve simultaneously */
// tenorShifts is expected to be the same size as the expiries/spreads
// shifts are additive
// modifies the object
// This method is outside the usual sensitivity framework
void CDSParSpreadsBase::bucketShift(const DoubleArray& tenorShifts)
{
    parSpreads->bucketShift(tenorShifts);
}

/** Shifts the object using given shift (see Theta::Shift)*/
bool CDSParSpreadsBase::sensShift(Theta* shift){
    valueDate = shift->rollDate(valueDate);
    // Reset the "local" cleanSpreads cache
    cleanSpreads.reset();
    return true; // deal with any components
}

/** Returns name identifying CDS Par Curve for RECOVERY_TWEAK */
string CDSParSpreadsBase::sensName(IRecoveryPerturb* shift) const{
    return name; 
}

/** Shifts the object using given shift */
bool CDSParSpreadsBase::sensShift(IRecoveryPerturb* shift){
    static const string method = 
        "CDSParSpreadsBase::sensShift(IRecoveryPerturb* shift)";

    double existingRecovery = parRecovery;
    parRecovery = shift->applyShift(parRecovery);
        
    // Checks if the tweaked recovery is valid !
    if (!isRecoveryValid()) {
        // must undo our change
        double shiftRecovery = parRecovery;
        parRecovery = existingRecovery;
        throw ModelException(method,
                             "Shifted recovery out of [0,1] "
                             "(Shifted recovery="
                             + Format::toString(shiftRecovery) + 
                             ")");
    }
    // Reset the "local" cleanSpreads cache
    cleanSpreads.reset();

    return false; /* none of our components has a RecoveryTweak
                   * type sensitivity */
}

/** Restores the object to its original form */
void CDSParSpreadsBase::sensRestore(IRecoveryPerturb* shift){
    parRecovery = shift->undoShift(parRecovery);
    // Reset the "local" cleanSpreads cache
    cleanSpreads.reset();
}


/** Returns the date when to stop tweaking after lastDate */
const DateTime CDSParSpreadsBase::stopTweaking(const DateTime& lastDate) const {
    DateTimeArray expiryDates = getExpiryDates();
    
    // assumes expiryDates array is sorted
    for(int i=0;i<expiryDates.size();i++) {
        if (expiryDates[i].isGreater(lastDate)) {
            return expiryDates[i];
        }
    } 
    return lastDate;
}

/**Returns the market recovery rate on defaulted assets given default at 
   time defaultDate. This allows the possibility of a term-structure of recovery rates. */
double CDSParSpreadsBase::getRecovery(const DateTime& defaultDate) const 
{
    return getRecovery();
}


double CDSParSpreadsBase::getRecovery() const {
    return parRecovery;
}

bool CDSParSpreadsBase::getAccrueFee() const {
    return parAccrueFee;
}

/** to be removed - should use bdc and hols from this object. Until
    then bdc and hols are optional, if null then values from this are
    used */
double CDSParSpreadsBase::getCurrentSpread(const DateTime& valueDate,
                                       const DateTime& maturityDate,
                                       const BadDayConvention* bdc,
                                       const Holiday* hols) const {
    if (priceViaCreditSpreadCurve)
    {
        return creditSpreadCurve->getCurrentSpread(valueDate, maturityDate, 
                                                   bdc? bdc: parBadDayConv.get(),
                                                   hols? hols: parHols.get());
    }
    else
    {
        return parSpreads->getCurrentSpread(valueDate, maturityDate, 
                                            bdc? bdc: parBadDayConv.get(),
                                            hols? hols: parHols.get());
    }
}

double CDSParSpreadsBase::getCurrentSpread(const DateTime& valueDate,
                                       const DateTime& maturityDate) const
{
    if (priceViaCreditSpreadCurve)
    {
        return creditSpreadCurve->getCurrentSpread(valueDate, maturityDate, 
                                                   parBadDayConv.get(), parHols.get());
    }
    else
    {
        return parSpreads->getCurrentSpread(valueDate, maturityDate, 
                                            parBadDayConv.get(), parHols.get());
    }
}

/** Returns FALSE, indicating the name has NOT defaulted.
 * Derived classes can override this method to keep track of defaults */
bool CDSParSpreadsBase::defaulted() const{
    return false;
}
    
/** Returns the date that this name defaulted. Since the name has
 * not happened, an 'empty' DateTime is returned.
 * Derived classes can override this method to keep track of defaults */
const DateTime& CDSParSpreadsBase::getDefaultDate() const{
    static DateTime emptyDate;
    return emptyDate;
}


/** Returns a DefaultRates object which gives access to 
    useful functionality including "default rates", aka clean default
    spreads. The information needed to bootstrap the clean spreads is
    obtained from this object 
    The DefaultRate returned is immutable and therefore will not be
    changed - which means that there is no need to clone it */
DefaultRatesSP CDSParSpreadsBase::defaultRates() const{
    if (!cleanSpreads){
        // Obtain the clean spread associated to this CDSParSpread 
        // from the cache
        cleanSpreads = getCachedDefaultRates(this, BASE_CDS_PAR_SPREADS);
    }
    return cleanSpreads;
}

/** return the bootstrapped dates */
DateTimeArray CDSParSpreadsBase::zeroDates() const
{
    DefaultRatesSP dR = defaultRates();
    DateTimeArraySP dates = dR->getDefaultDates();
    return *(dates.get());
}

// accessor methods for logOfDiscFactorKey
IDiscountCurve::IKey* CDSParSpreadsBase::getDiscountKey() const
{
    return discount->logOfDiscFactorKey();
}

DefaultRates::IKey* CDSParSpreadsBase::getRiskyKey() const
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
IDiscountCurve::IKey* CDSParSpreadsBase::logOfDiscFactorKey() const
{
    return new RiskyLogOfDiscFactorKey(this);
}

/** Checks that recovery is in [0,1] */
bool CDSParSpreadsBase::isRecoveryValid() const {
    return (parRecovery >= 0.0 && parRecovery <= 1.0);
}                                            


/** Returns an processed vol - which combines the vol market data with the
    instrument data in the volRequest */
CVolProcessed* CDSParSpreadsBase::getProcessedVol(
    const CVolRequest* volRequest) const{
    if (cdsVol.isEmpty()){
        throw ModelException("CDSParSpreadsBase::getProcessedVol",
                             "Pricing model requires CDS Vol be supplied for "
                             "CDSParSpreadsBase curve with name "+getName());
    }
    return cdsVol->getProcessedVol(volRequest, this);
}

/** Returns the same as getName() */
string CDSParSpreadsBase::getCorrelationName() const{
    return getName();
}


/** Invoked when this class is 'loaded' */
class CDSParSpreadsBaseHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CDSParSpreadsBase, clazz);
        SUPERCLASS(BootstrappableCDSParSpreads);
        IMPLEMENTS(ICreditCurve);
        IMPLEMENTS(Theta::IShift);
        IMPLEMENTS(CreditSpreadRhoParallel::RestorableShift);
        IMPLEMENTS(CreditSpreadRhoPointwise::IRestorableShift);
        IMPLEMENTS(RecoveryShiftBase::IRestorableShift);
        IMPLEMENTS(ParSpreadMaxTweakSize::IShift);
        IMPLEMENTS(Duration::IParHandlerWithClosedForm);
        IMPLEMENTS(Duration::IParHandlerWithoutClosedForm);
        EMPTY_SHELL_METHOD(defaultCDSParSpreadsBase);
        FIELD(name, "name");
        FIELD(valueDate, "today");
        FIELD_MAKE_OPTIONAL(valueDate);// retrieve from market data cache
        FIELD(parSpreads, "par Spread Curve");
        FIELD(parRecovery,
                     "recovery rate for bonds underlying par swaps");
        FIELD(parSwapFreq, "swap freq for par spreads (1, 2, 4)");
        FIELD(parDCC, "day-count-convention string for par CDS");
        FIELD(parBDC, "bad day convention for par CDS");
        FIELD_MAKE_OPTIONAL(parBDC); // to be made mandatory ideally
        FIELD(parHols, "holidays for par CDS");
        FIELD_MAKE_OPTIONAL(parHols); // to be made mandatory ideally
        FIELD(parAccrueFee, "make accrued payment on default?");
        FIELD_MAKE_OPTIONAL(parAccrueFee);
        FIELD(spotOffset, "gives effective date of par curve");
        FIELD_MAKE_OPTIONAL(spotOffset);
        FIELD(discount, "Discount yield curve");
        FIELD(prepay, "ABS prepay curve");
        FIELD_MAKE_OPTIONAL(prepay);
        FIELD(cdsVol, "Volatility of associated CDS Swaptions");
        FIELD_MAKE_OPTIONAL(cdsVol);
        FIELD_NO_DESC(parDayCountConv);
        FIELD_MAKE_TRANSIENT(parDayCountConv);
        FIELD_NO_DESC(parBadDayConv);
        FIELD_MAKE_TRANSIENT(parBadDayConv);
        FIELD(priceViaCreditSpreadCurve, "");
        FIELD_MAKE_TRANSIENT(priceViaCreditSpreadCurve);
        FIELD(creditSpreadCurve,"");
        FIELD_MAKE_TRANSIENT(creditSpreadCurve);
        FIELD(csRestore,"");
        FIELD_MAKE_TRANSIENT(csRestore);

        //        ClassSetAcceptMethod(CDSParSpreadsBase::acceptWrapperNameCollector);
    }

    static IObject* defaultCDSParSpreadsBase(){
        return new CDSParSpreadsBase();
    }
};

// ---------------------------------
// ICDSBootstrappable implementation
// ---------------------------------

/** Returns the corresponding discount curve */
YieldCurveConstSP CDSParSpreadsBase::getYieldCurve() const {
    return discount.getSP();
}

/** Returns the value date */
const DateTime& CDSParSpreadsBase::getValueDate() const {
    return valueDate;
}

/** Returns the 'effective date' = date used as reference to
 * compute cash flow dates */
const DateTime CDSParSpreadsBase::getEffDate() const {
    return spotDate(getValueDate());
}

/** Returns bad day convention used to adjust cash flow dates */
BadDayConventionConstSP CDSParSpreadsBase::getBadDayConvention() const {
    return parBadDayConv;
}

/** Returns holidays used to adjust cash flow dates */
HolidayConstSP CDSParSpreadsBase::getHolidays() const {
    return parHols.getSP();
}

/** Returns the expiries dates (absolute: eg. 27/02/2007) */
const DateTimeArray CDSParSpreadsBase::getExpiryDates() const {
    return getExpiryDates(
        getEffDate(),
        getBadDayConvention().get(),
        getHolidays().get());    
}

/** Returns the protection end date corresponding to expiry = index */
const DateTime CDSParSpreadsBase::getProtectionEndDate(int index) const{
    //unadjusted
    ExpiryArrayConstSP expiries = parSpreads->getExpiries();
    DateTime protEndDate = (*expiries)[index]->toDate(getEffDate());
    return protEndDate;
}

/** Returns TRUE if paying fee accrued */
bool CDSParSpreadsBase::isFeeAccrued() const {
    return getAccrueFee();
}

/** Returns CDS par spreads */
DoubleArrayConstSP CDSParSpreadsBase::getParSpreads() const {
    return parSpreads->getParSpreads();
}

/** Returns the expiries (relative: eg. 1M) */
ExpiryArrayConstSP CDSParSpreadsBase::getParSpreadsExpiries() const {
    return getExpiries();
}

/** Return prepay curve */

IDecretionCurveConstSP CDSParSpreadsBase::getPrepayCurve() const {
    return prepay.getSP();
}

DecretionCurveConstSP CDSParSpreadsBase::getPrepayCurveObj() const {
    return prepay.getSP();
}

CashFlowArraySP CDSParSpreadsBase::getCashFlowArray(int index) const
{
    CashFlowArraySP cfsp = BootstrappableCDSParSpreads::getCashFlowArray(index);
    DoubleArrayConstSP upf = parSpreads->getUpfronts();
    if(upf.get()) {
        if(cfsp->front().date.equals(getEffDate())){
           cfsp->front().amount += (*upf)[index];
        }
    }
    return cfsp;
}

CashFlowArraySP CDSParSpreadsBase::getCashFlowArray(CreditIndexBasisConstSP indexBasis,
                                                    int index) const
{
    CashFlowArraySP cfsp = BootstrappableCDSParSpreads::getCashFlowArray(indexBasis, index);
    DoubleArrayConstSP upf = parSpreads->getUpfronts();
    if(upf.get()) {
        if(cfsp->front().date.equals(getEffDate())){
           cfsp->front().amount += (*upf)[index];
        }
    }
    return cfsp;
}


//------------------------------
// Duration::IParHandler methods
//------------------------------

//return a means of tweaking the MarketObject in a pointwise manner
//assumption is that the results of this is a VectorResult.....
VectorRiskPropertySensitivityConstSP CDSParSpreadsBase::getPointwiseTweaker() const {
    //1 bp shift explictly required
    return VectorRiskPropertySensitivityConstSP(newParSpreadRhoPointwise(0.0001));
}

//return a par instrument for the specified maturity
InstrumentSP CDSParSpreadsBase::getParInstrument(const ExpirySP maturity) const {
    static const string method = "CDSParSpreadsBase::getParInstrument";

    try {
        ExpiryArrayConstSP expiries = parSpreads->getExpiries();
        DoubleArrayConstSP spreads = parSpreads->getParSpreads();

        //find the specified maturity
        int i = maturity->search(expiries.get());
        //now construct the par instrument
        ParCDSSP cds(new ParCDS(discount.getName(), getName(),
                                spotDate(valueDate), maturity, 
                                (*spreads)[i], parSwapFreq));
        return cds;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

//return benchmarks for par instruments
//these will determine the durations to calculate unless specified
ExpiryArrayConstSP CDSParSpreadsBase::getParBenchmarks() const {
    return getParSpreadsExpiries();
}

//return closed form solution
ExpiryResultArraySP CDSParSpreadsBase::getDuration(const Duration* duration) const {
    return CDSGetDuration(duration, 
                          getParBenchmarks(),
                          getYCName(),
                          getName());
}

//static function to return closed form solution
ExpiryResultArraySP CDSParSpreadsBase::CDSGetDuration(const Duration* duration,
                                                      ExpiryArrayConstSP altBenchmarks,
                                                      string altYCName,
                                                      string altName){
    MarketDataSP market = duration->getMarket();
    IModelSP model = duration->getModel();

    ExpiryArrayConstSP benchmarks = duration->getBenchmarks();
    if (!benchmarks) {
        benchmarks = altBenchmarks;
    }
    
    // Set the discount yield curve name set as a Duration parameter if any,
    // or otherwise the one in the cdsParSpreads
    string discountYCName;
    YieldCurveWrapper durDiscount = duration->getDiscount();
    if (durDiscount.isEmpty()) {
        // use the yield curve of cdsParSpreads
        discountYCName = altYCName;
        durDiscount.setName(altYCName);
    } else {
        // use the discount curve we are given
        discountYCName = durDiscount->getName();
    }
    model->setDomesticYCName(discountYCName);

    //fetch the par spreads
    ICDSParSpreadsWrapper cdsParSpreads;
    cdsParSpreads.setName(altName);
    ICDSParSpreads::getMarketData(model.get(),
                                  market.get(),
                                  discountYCName,
                                  cdsParSpreads);

    //fetch the discount (may not have been done....)
    durDiscount.getData(model.get(), market.get());

    return CDSPricer::riskyDurations(ICDSParSpreadsConstSP(cdsParSpreads.get()),
                                     durDiscount.getSP(),
                                     benchmarks);

}


#ifdef CDS_BACKWARD_COMPATIBILITY

/** Set the bad day convention */
void CDSParSpreadsBase::setBadDayConvention(BadDayConventionSP bdc) {
    // Test if bad day convention is already here
    if (!parBadDayConv.get() || parBDC == "No value") {
        if (!bdc.get()) {
            throw ModelException(
                "CDSParSpreadsBase::setBadDayConvention",
                "Attempting to set an undefined bad day convention !"); 
        } else {
            parBadDayConv = bdc;
        }
    }
}

/** Set the holidays */
void CDSParSpreadsBase::setHolidays(HolidaySP holidays) {
    // Test if holidays are already here
    if (!this->parHols.getMO().get()) {
        if (!holidays.get()) {
            throw ModelException(
                "CDSParSpreadsBase::setHolidays",
                "Attempting to set undefined holidays !"); 
        } else {
            parHols.setObject(holidays);
        }
    }
}

#endif    

/** Build a credit spread curve from the cds curve */
/** Taken from ConvBond implementation */
void CDSParSpreadsBase::buildCreditSpreadCurve() const {
    static const string method = "CDSParSpreadsBase::buildCreditSpreadCurve";
    try {
        //-since we are building the credit spread curve
        //-we need to (potentially) override the priceViaCreditSpreadCurve flag
        //-to ensure the spreads come from the par spreads
        bool priceFlag = priceViaCreditSpreadCurve;
        priceViaCreditSpreadCurve = false;

        // convert the default rates to ANNUAL 365F basis before building the risky zero curve
        DefaultRatesSP psDefRates = defaultRates();
        CashFlowArraySP defaultRates = psDefRates->convertToSpot(valueDate, Actual365F());

        CreditSpreadCurveSP csCurve = discount->makeCreditSpreadCurve
            (name, *defaultRates, parRecovery);

        // copy it to the transient field
        creditSpreadCurve = csCurve;
        // make backup for restore too
        csRestore = CreditSpreadCurveSP(copy(csCurve.get()));

        //reset the state of the pricing flag
        priceViaCreditSpreadCurve = priceFlag;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** Returns name identifyer for CREDIT_SPREAD_RHO_PARALLEL */
string CDSParSpreadsBase::sensName(CreditSpreadRhoParallel* shift) const {
    return name;
}

/** Shifts the object using given shift */
bool CDSParSpreadsBase::sensShift(CreditSpreadRhoParallel* shift) {
    static const string method = 
        "CDSParSpreadsBase::sensShift(CreditSpreadRhoParallel)";
    try {
        //build credit spreads if necessary
        if (!(creditSpreadCurve.get())) {
            // May have been built earlier for POINTWISE TWEAK
            buildCreditSpreadCurve();
        }

        // Do the shift
        creditSpreadCurve->sensShift(shift);

        // Reset the "local" cleanSpreads cache
        cleanSpreads.reset();

        // Set to price via credit spreads - do this last (in case of anything
        // failing above)
        priceViaCreditSpreadCurve = true;

        //no sub elements support this tweak
        return false;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
    return false;
}

/** Restores the object to its original form */
void CDSParSpreadsBase::sensRestore(CreditSpreadRhoParallel* shift) {
    static const string method = 
        "CDSParSpreadsBase::sensRestore(CreditSpreadRhoParallel)";
    try {
        creditSpreadCurve = CreditSpreadCurveSP(copy(csRestore.get()));
        priceViaCreditSpreadCurve = false;
        // Reset the "local" cleanSpreads cache
        cleanSpreads.reset();
    } catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** Returns name identifying yield curve for CREDIT_SPREAD_RHO_POINTWISE */
string CDSParSpreadsBase::sensName(CreditSpreadRhoPointwise* shift) const {
    return name;
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryArrayConstSP CDSParSpreadsBase::sensExpiries(CreditSpreadRhoPointwise* shift) const {
    //build credit spreads if necessary
    if (!(creditSpreadCurve.get())) {
        // May have been built earlier for POINTWISE TWEAK
        buildCreditSpreadCurve();
    }

    // get the market data for the parSpreads
    return creditSpreadCurve->getExpiries();
}


/** Shifts the object using given shift */
bool CDSParSpreadsBase::sensShift(CreditSpreadRhoPointwise* shift) {
    static const string method = "CDSParSpreadsBase::sensShift(CreditSpreadRhoPointwise)";
    try {
        if (!(creditSpreadCurve.get())) {
            // May have been built earlier for POINTWISE TWEAK
            buildCreditSpreadCurve();
        }

        // Do the shift
        creditSpreadCurve->sensShift(shift);

        // Reset the "local" cleanSpreads cache
        cleanSpreads.reset();

        // Set to price via credit spreads - do this last (in case of anything
        // failing above)
        priceViaCreditSpreadCurve = true;

        return false; // dont want to go on a tweak the par curve rho which inherits from CreditSpreadCurve
    } catch (exception &e) {
        throw ModelException(e, method);
    }
    return false;
}

/** Restores the object to its original form */
void CDSParSpreadsBase::sensRestore(CreditSpreadRhoPointwise* shift) {
    static const string method = "CredDefSwap::sensRestore";
    try {
        priceViaCreditSpreadCurve = false;
        // Reset the "local" cleanSpreads cache
        cleanSpreads.reset();
        creditSpreadCurve = CreditSpreadCurveSP(copy(csRestore.get()));
    } catch (exception &e) {
        throw ModelException(e, method);
    }
}

string CDSParSpreadsBase::sensName(ParSpreadMaxTweakSize* shift) const {
    //max tweak size applies to the embedded par spread curve
    return parSpreads->getName();
}

DoubleArrayConstSP CDSParSpreadsBase::sensMaxTweakSize(ParSpreadMaxTweakSize* shift) const {
    return calculateTweakSizes();
}

/** To support the creation of risky zero curves */
IYieldCurveSP CDSParSpreadsBase::makeRiskyCurve(
	const IYieldCurve& yc, 
	const DateTime*    maturityDate) const
{
    //result depends on whether the CDS curve is switched to work as a CreditSpreadCurve
    IObjectSP riskyCmpt = getRiskyComponent();

    //now do the same checks again, since we could have a CDSParSpreads or a CreditSpreadCurve
    if (CreditSpreadCurve::TYPE->isInstance(riskyCmpt))
    {
        const CreditSpreadCurve* rskCurve = 
            dynamic_cast<const CreditSpreadCurve*>(riskyCmpt.get());

        return yc.makeRiskyCurve(*rskCurve, maturityDate);
    }
    else
    {
        // same as for ICDSParSpreads
        return ICDSParSpreads::makeRiskyCurve(yc, maturityDate);
    }
}

IObjectSP CDSParSpreadsBase::getRiskyComponent() const {
    if (priceViaCreditSpreadCurve) {
        if (!(creditSpreadCurve.get())) {
            // May have been built earlier
            buildCreditSpreadCurve();
        }

        return creditSpreadCurve;
    }
    else {
        CDSParSpreadsBase* aCopy = const_cast<CDSParSpreadsBase*>(this);
        return CDSParSpreadsBaseSP(aCopy);
    }
}


//------------------------------------------
//  IBadDayAdjuster methods
//------------------------------------------
/** Returns "date" bad day adjusted using the bad day convention
 * and holidays in this object */
DateTime CDSParSpreadsBase::badDayAdjust(const DateTime& date) const {
    return parBadDayConv->adjust(date, parHols.get());
}

/** Add a number of business days to a date */
DateTime CDSParSpreadsBase::addBusinessDays(const DateTime& from, int busDays) const {
    return parHols->addBusinessDays(from, busDays);
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
const DefaultRatesSP CDSParSpreadsBase::getCachedDefaultRates(
    const ICDSParSpreads* cds,
    const TypeOfEntry     entryType) const 
{
    static const string method = "CDSParSpreadsBase::getCachedDefaultRates";

    if (!cache.get()) {
        throw ModelException(method, 
                             "Cache is NULL trying to obtain default rates");
    }
    return cache->getEntry(cds, entryType);
}

/** Provides a facility for previously constructed DefaultRates
    *  objects to be added to the cache */
const bool CDSParSpreadsBase::cacheDefaultRates(const ICDSParSpreads* cds,
                                                const TypeOfEntry     entryType,
                                                const DefaultRatesSP  entry)
{
    static const string method = "CDSParSpreadsBase::cacheDefaultRates";

    if (!cache.get()) {
        throw ModelException(method, 
                             "Cache is NULL trying to add default rates");
    }
    return cache->addEntry(cds, entryType, entry);
}



/** Returns a DefaultRates object.
 * This method does NOT use the (potentially availabe) caching mechanism
 * which avoids the slow process of calculating default rates everytime.
 * It is invoked from inside the cache to calculate the default rates
 * for the first time and its direct use from any external classes is
 * strongly DISCOURAGED: the "defaultRates" method should be used instead.
 * Here we are in CDSParSpreadsBase so the way to compute the clean spreads
 * is bootstrapping... so do that */
const DefaultRatesSP CDSParSpreadsBase::computeDefaultRatesForCache() const {
    static const string method = 
        "CDSParSpreadsBase::computeDefaultRatesForCache";

    try {
        if (priceViaCreditSpreadCurve) {
            if (!discount){
                throw ModelException(method, "Discount curve not specified");
            }

            // taken from CorporateBond
            // use credit spread methodology for calculating CREDIT_SPREAD_RHO's
            CDSHelper::CParSpreadDefaultRates psAnnDefRates =
                CDSHelper::CParSpreadDefaultRates(creditSpreadCurve.get(), 
                                                  discount.getSP(),
                                                  valueDate,
                                                  parRecovery);
    
            return psAnnDefRates.createFwdContDefCurve(valueDate);
        }
        else {
            return CDSPricer::bootstrap(ICDSBootstrappableConstSP(this));
        }
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}
    

/** Resets the pointer to the cache - The cache will not be used
 * in this object from this point onwards:
 * This is invoked for the CDSParSpreadsBase clone stored inside the cache, 
 * and is required so that the cache is only pointed at from external 
 * CDSParSpreadsBase. Otherwise the cache object would never get freed */
void CDSParSpreadsBase::resetCache() const {
    cache.reset();
}

/** Hash code function - the cache needs improved performance compared
 * to the default "hashCode" function in CObjet: only the required
 * components are hashed (see comments in equalToOpt below) */
int CDSParSpreadsBase::hashCodeOpt() const {
    int hcode = (size_t) getClass();
    hcode ^= valueDate.hashCode();
    hcode ^= spotOffset; // (integer)
    hcode ^= hash_string(parDCC);
    /* parBadDayConv *may* be parBDC but it may also be set in 
     * setBadDayConvention, so do not trust parBDC (which is optional)*/
    hcode ^= hash_string(parBadDayConv->toString());

    hcode ^= parSwapFreq; // (integer)
    hcode ^= parSpreads->hashCode();
    hcode ^= CDouble::hashCode(parRecovery);
    hcode ^= CBool::hashCode(parAccrueFee);
    hcode ^= discount.get()->zeroCurveHash();

    // prepay curve is an OPTIONAL parameter
    if (prepay.get())
        hcode ^= prepay->hashCode();

    hcode ^= CBool::hashCode(priceViaCreditSpreadCurve);
    if (priceViaCreditSpreadCurve) {
        // If the flag priceViaCreditSpreadCurve is set the creditSpreadCurve
        // will be used for pricing and therefore needs to be cached.
        // It should NOT be cached if the flag is not set.
        if (creditSpreadCurve.get()) {
            hcode ^= creditSpreadCurve->hashCode();
        }
        else {
            ; // The else case should never happen - but even if it did, the
              // hash would be different so there is nothing else to do
        }
    }
    return hcode;
}


/** Comparison function - the cache needs improved performance compared 
 * to the default "equalTo" function in CObjet: only the required
 * components are compared. This means that the objects being compared
 * are not required to be identical: only equal as far as the default
 * rates calculation is concerned (e.g, if the volType is not used in
 * the calculation, the equalToOpt method should not fail if the volTypes
 * are different) */
bool CDSParSpreadsBase::equalToOpt(const ICDSParSpreads* cds) const {
    static const string method = "CDSParSpreadsBase::equalToOpt";
    try {
        if (this == cds) { // Obvious first test
            return true;
        }
        if (!cds || cds->getClass() != getClass()) {
            return false;
        }

        const CDSParSpreadsBase* cdsBase = STATIC_CAST(CDSParSpreadsBase, cds);

        if (!valueDate.equals(cdsBase->valueDate)) {
            return false;
        }
        if (spotOffset != cdsBase->spotOffset) {
            return false;
        }
        if (parDayCountConv->toString() != cdsBase->parDayCountConv->toString()){
            return false;
        }
        if (parBadDayConv->toString() != cdsBase->parBadDayConv->toString()){
            return false;
        }
        if (parSwapFreq != cdsBase->parSwapFreq) {
            return false;
        }
        if (!parHols->equals(cdsBase->parHols.get())) {
            return false;
        }

        DoubleArrayConstSP  spread1 = parSpreads->getParSpreads();
        ExpiryArrayConstSP  expiry1 = parSpreads->getExpiries();
        CDoubleArrayConstSP spread2 = cdsBase->parSpreads->getParSpreads();
        ExpiryArrayConstSP  expiry2 = cdsBase->parSpreads->getExpiries();
        if ((spread1->size() != spread2->size())   ||
            (expiry1->size() != expiry2->size()))
        { 
            return false;
        }

        int len = spread1->size();
        for (int i = 0; i < len; i++) {
            if ((*spread1)[i] != (*spread2)[i]) {
                return false;
            }
            if (!(*expiry1)[i]->equals((*expiry2)[i].get())) {
                return false;
            }
        }

        if(parSpreads->getUpfronts().get())
        {
            DoubleArrayConstSP upfront1 = parSpreads->getUpfronts();
            DoubleArrayConstSP upfront2;
            if(cdsBase->parSpreads->getUpfronts().get())
            {
                upfront2 = cdsBase->parSpreads->getUpfronts();
            } else {
                return false;
            }

            if ((upfront1->size() != upfront2->size()))
            { 
                return false;
            }
            
            int upflen = upfront1->size();
            for (int i = 0; i < upflen; i++) {
                if ((*upfront1)[i] != (*upfront2)[i]) {
                    return false;
                }
            }
        }

        //prepay curve
        if(prepay.get() && !prepay->equalTo(cdsBase->prepay.get())) {
            return false;
        }

            
        if (parRecovery != cdsBase->parRecovery) {
            return false;
        }
        if (parAccrueFee != cdsBase->parAccrueFee) {
            return false;
        }
        if (!discount->zeroCurveEquals(cdsBase->discount.get())) {
            return false;
        }
        
        if (priceViaCreditSpreadCurve != cdsBase->priceViaCreditSpreadCurve) {
            return false;
        }

        if (priceViaCreditSpreadCurve) {
            // If the flag priceViaCreditSpreadCurve is set the 
            // creditSpreadCurve will be used for pricing and therefore needs
            //  to be compared - but not otherwise.
            if (!creditSpreadCurve->equalTo(cdsBase->creditSpreadCurve.get())) {
                return false;
            }
        }

        return true;
    }
    catch(exception &e){
        throw ModelException(e, method);
    }
}

double CDSParSpreadsBase::survivalProb(
                            const DateTime& d1, 
                            const DateTime& d2) const
{
    return this->defaultRates()->calcTotalDefaultPV(d1, d2);
}



/**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
double CDSParSpreadsBase::survivalProb(const DateTime& dt) const
{
    return survivalProb(getValueDate(), dt);
}

/*TODO: this currently return PV with no delay - need to implement fixed date recovery.*/
double CDSParSpreadsBase::protectionPV(const DateTime& paymentDate, 
                                       const DateTime& startDt, 
                                       const DateTime& endDt,
                                       RecoveryType    recTyp,
                                       const DateTime& recoveryDate) const
{
    if (recTyp == IDiscountCurveRisky::RECOVER_0)
    {
        return .0;
    }

    double  protectionValue = (1 - survivalProb(startDt, endDt)) * this->getYieldCurve()->pv(paymentDate, recoveryDate);

    if (recTyp == IDiscountCurveRisky::RECOVER_R)
    {
        protectionValue *= getRecovery();
    }
    else if (recTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
    {
        protectionValue *= (1-getRecovery());
    }

    return protectionValue;
}


/**Returns the value at paymentDate (and conditional on no default before then) 
 * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
 * in case of default between startDate and endDate, and zero otherwise.
 * The protection pay-out is paid at date recoveryTiming.toDate(defaultDate).
 * This allows a recovery delay if the Expiry is a relative date, or an
 * absolute recovery if it's a fixed date (e.g. for use in recovery-at-maturity
 * calculations for unconditional settlement. If the Expiry pointer is null,
 * then no delay.
 * Whether integration is done continuously or discretely and what, 
 * if any approximations,
 * are used is up to the curve implementation. */
double CDSParSpreadsBase::protectionPV(const DateTime& paymentDate, 
                                       const DateTime& startDt, 
                                       const DateTime& endDt,
                                       RecoveryType    recTyp,
                                       double          recoveryDelay) const
{   
    if (recTyp == IDiscountCurveRisky::RECOVER_0)
    {
        return .0;
    }

    DateTimeArray               noDates = DateTimeArray();
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual360());
    double                      protectionValue = .0;
    double                      accrualValue = .0;

    CDSPricer   cdsPricer 
        = CDSPricer(.0,                     // double couponRate,
                    noDates,                // paymentDates,
                    this->defaultRates(),   // defaultRates,
                    1.0,                    // notional,
                    0.0,                    // recovery,
                    false,                  // payAccruedFeed,
                    dcc,                    // swpAccrualDCC
                    paymentDate,            // valueDate,
                    startDt,                // protectionStartDate,
                    endDt,                  // protectionEndDate,
                    startDt,                // accruedStartDate,
                    getYieldCurve(),
                    getPrepayCurve());
    
    if (recoveryDelay == .0)
    {
        cdsPricer.calculateDefaultPayment(  false, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::NO_DELAY_DISCOUNTING);
    }
    else
    {
        cdsPricer.calculateDefaultPayment(  false, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::DELAYED_DISCOUNTING,
                                            recoveryDelay);
    }

    if (recTyp == IDiscountCurveRisky::RECOVER_R)
    {
        protectionValue *= getRecovery();
    }
    else if (recTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
    {
        protectionValue *= (1-getRecovery());
    }

    return protectionValue;
}

/*TODO: this currently return PV with no delay - need to implement fixed date recovery.*/
double CDSParSpreadsBase::annuityPV(
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
        CDSPricer   cdsPricer 
            = CDSPricer(paymentsSP,                // paymentDates,
                        this->defaultRates(),   // defaultRates,
                        1.0,                    // notional,
                        0.0,                    // recovery,
                        true,                   // payAccruedFeed,
                        dcc,                    // swpAccrualDCC
                        paymentDate,            // valueDate,
                        accrueStartDate,            // protectionStartDate,
                        payments[payments.size()-1].date,            // protectionEndDate,
                        accrueStartDate,            // accruedStartDate,
                        getYieldCurve());
        
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

    return couponValue + accrualValue * this->getYieldCurve()->pv(paymentDate, recoveryDate);
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
 * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
 * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
 * this method should return the same value as pv(payments, paymentDate), inherited
 * from IDiscountCurve. */
double CDSParSpreadsBase::annuityPV(const CashFlowArray& payments,
                                    const DateTime&      paymentDate,
                                    RecoveryType         accruedRecTyp,
                                    double               recoveryDelay,
                                    DateTime             accrueStartDate) const
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
        CDSPricer   cdsPricer 
            = CDSPricer(paymentsSP,                // paymentDates,
                        this->defaultRates(),   // defaultRates,
                        1.0,                    // notional,
                        0.0,                    // recovery,
                        true,                   // payAccruedFeed,
                        dcc,                    // swpAccrualDCC
                        paymentDate,            // valueDate,
                        accrueStartDate,            // protectionStartDate,
                        payments[payments.size()-1].date,            // protectionEndDate,
                        accrueStartDate,            // accruedStartDate,
                        getYieldCurve());
        
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
double CDSParSpreadsBase::risklessPV(const DateTime& date1, 
                                     const DateTime& date2) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(date1, date2);
}
    
/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double CDSParSpreadsBase::risklessPV(const DateTime& date) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(date);
}
    
/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double CDSParSpreadsBase::risklessPV(const CashFlowArray& cashFlows,
                                     const DateTime&      baseDate) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(cashFlows, baseDate);
}


double CDSParSpreadsBase::pv(const DateTime& date1, 
                             const DateTime& date2) const 
{
    return this->defaultRates()->calcTotalDefaultPV(date1, date2) * getYieldCurve()->pv(date1, date2);
}
    
/** Compute price for settlement today of a zero-coupon bond
 * maturing at date. Note settlement is to TODAY, not to
 * the curve spot date. This is because some curves
 * may have ambiguous spot-dates - for example should a combined
 * credit and rates curve have spot date T+1 or T+2?
 * @param date To get discount factor/PV for
 * @return Discount factor between today & given date
 */
double CDSParSpreadsBase::pv(const DateTime& date) const
{
    return pv(getValueDate(), date);
}
    
/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double CDSParSpreadsBase::pv(const CashFlowArray& cashFlows,
                             const DateTime&      baseDate) const
{   
    double      value = .0;
        
    for (int i=0; i < cashFlows.size(); ++i)
        {
            const DateTime cfDate = cashFlows[i].date;

            if (cfDate.isGreater(baseDate))
                {
                    value += cashFlows[i].amount * pv(baseDate, cfDate);
                }
        }

    return value;
}


CClassConstSP const CDSParSpreadsBase::TYPE = CClass::registerClassLoadMethod(
    "CDSParSpreadsBase", typeid(CDSParSpreadsBase), CDSParSpreadsBaseHelper::load);

DEFINE_TEMPLATE_TYPE(CDSParSpreadsBaseWrapper);

DRLIB_END_NAMESPACE
