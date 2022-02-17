//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetCVB.cpp
//
//   Description : asset with convertible bond as its underlying
//
//   Date        : 2 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/ConvBond.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/AllStrikes.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE

class CAssetCVB: public CAsset,
                 public virtual DerivativeAsset,
                 public virtual CAsset::IStruckable,
                 public virtual SpotLevelProbability::Shift,
                 public virtual VegaSkewParallel::IShift,
                 public virtual VegaSkewPointwise::IShift,
                 public virtual IAssetFairValue,
                 public virtual INextStrike,
                 public virtual ObjectIteration::IOverride,
                 public virtual DeltaSurface::IShift,
                 public virtual VolRelativeShift::IShift {
private:
    /// fields ////
    string           name;
    double*          SpotPrice;         // spot price of cvb (not stock)
    double           cvbDeltaOverride;  // not supported yet
    CVolBaseWrapper  vol;               // volatility of cvb
    BorrowCurveSP    borrowCurve;       // Borrow curve of cvb
    ConvBondSP       CVBData;           // cvb data
    IModelSP         CVBModel;          // mode for pricing a cvb
    // true = asset fixing (spot) price is clean price, false = spot
    // fixing on dirty price.
    bool             FixOnCleanPrice;
    SettlementSP     settlement; // cvb asset settlement - default T+3
    // transient fields ///
    mutable EquitySP equity; // used for computing forward price of cvb
    mutable bool     resetEquitySpot; // true => must recalc equity spot
    bool             useCvbTheoPrice; // true => use theoretical cvb price
    CControlSP       control; // not registered $unregistered
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        if (name.empty()){
            throw ModelException("CAssetCVB::validatePop2Object()",
                                 "AssetCVB name is empty.");
        }
    }

    /** clone() has to keep the same model for same grid tweaking */
    IObject* clone() const {
        // copy registered fields first
        IObject*  obj = CAsset::clone();
        // shallow copy model to allow same grid tweaking
        CAssetCVB& cvbAsset = dynamic_cast<CAssetCVB&>(*obj);
        cvbAsset.CVBModel = CVBModel;
        cvbAsset.control = control; // shallow copy over
        return obj;
    }

    /** returns the spot price */
    virtual double getSpot() const{
        checkEquity();
        return equity->spot();
    }
    
    /** returns the asset name */
    virtual string getName() const{
        return name;
    }

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const{
        if (!vol) {
            throw ModelException("CAssetCVB::getProcessedVol",
                                 "AssetCVB does not have a vol object.");
        }
        return vol->getProcessedVol(volRequest, this);
    }

    /** Calculates the expected spot price of the asset at the given date.
        Do not use this repeatedly to calculate values over a set of
        dates (poor performance) - instead use other fwdValue method */
    virtual double fwdValue(const DateTime& date) const{
        checkEquity();
        checkFwdDate(date);
        return equity->fwdValue(date);
    }

    /** Default implementation is provided */
    virtual void fwdValue(const DateTimeArray& dateList,
                          CDoubleArray&        result) const{
        checkEquity();
        checkFwdDate(dateList);
        equity->fwdValue(dateList, result);
    }

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const{
        checkEquity();
        checkFwdDate(dateList);
        equity->fwdValue(dateList, algo, result);
    }

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const{
        return CVBData->discountYieldCurveName();
    }

    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const{
        checkEquity();
        checkFwdDate(fwdDate);
        DateTimeArray dateArray(1, fwdDate);
        CDoubleArray  result(1);
        equity->fwdFwd(spotDate, spotDate, spot, dateArray, result);
        return (result[0]);
    }

    /** Array version of fwdFwd */
    virtual void fwdFwd(const DateTime&      spotDate,
                        double               spot, 
                        const DateTimeArray& fwdDates,
                        DoubleArray&         results) const{
        checkEquity();
        checkFwdDate(fwdDates);
        equity->fwdFwd(spotDate, spotDate, spot, fwdDates, results);
    }

   /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const{
        return equity->settles(tradeDate);
    }

    virtual void getSensitiveStrikes(
        const CVolRequest* volRequest,
        OutputNameConstSP outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP sensitiveStrikes) const{
        checkEquity();
        // to do: alter getSensitiveStrikes signature
        if (!!vol) {
            const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
                *(IObject*)(volRequest));
            if (sensStrikeDesc.forwardOnly == false && 
                outputName->equals(vol->getName()) ) {
                request.getSensitiveStrike(equity->spot(), sensitiveStrikes);
            }
            // really should pass onto asset in CVB - to do
        }
    }

    /** return a pdf calculator provided our vol supports it */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const{
        if (!!vol) {
            if (IPDFCalculator::TYPE->isInstance(vol.get())) {
                const IPDFCalculator* pdf = 
                    dynamic_cast< const IPDFCalculator*>(vol.get());
                return pdf->getPDFCalculator(request, this);
            }
            throw ModelException("AssetCVB::pdfCalculator",
                                 "vol of type (" + vol->getClass()->getName() +
                                 ") has no pdf calculator");
        } else {
            throw ModelException("CAssetCVB::getProcessedVol",
                                 "AssetCVB does not have a vol object.");
        }
    }

    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&            sampleDate,
                            const ObservationType*      obsType,
                            const ObservationSource*    source,
                            const FixingType*           fixType,
                            const IObservationOverride* overrides,
                            const SamplingConvention*   sampleRule) const{
        throw ModelException("AssetCVB::pastValue", 
                    "Cannot implement this method for AssetCVB");
    }

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const {
        throw ModelException("AssetCVB::observationDate", 
                    "Cannot implement this method for AssetCVB");
    }

    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const{
        throw ModelException("AssetCVB::isHoliday", 
                    "Cannot implement this method for AssetCVB");
    }

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const {
        throw ModelException("AssetCVB::addPastSampleEvent", 
                    "Cannot implement this method for AssetCVB");
    }

    /** get data from market */
    virtual void getMarket(const IModel* model, const MarketData* market){
        try{
            // call getInstrumentAndModelMarket to fill cvb data and stock data
            // underlying the cvb instrument
            CVBModel->getInstrumentAndModelMarket(market, CVBData.get());
            borrowCurve->getMarket(CVBModel.get(), market);
            settlement->getMarket(CVBModel.get(), market);
            // fill vol 
            if (vol.getName() != "") {
                vol.getData(model, market);
            }
            
            // now build equity - important to build it before the
            // first pricing call to ensure that it gets tweaked (eg
            // the ESW clones the entire instrument)
            buildEquity();
            // flag that we cvb price is not known
            resetEquitySpot = true;
        } catch (exception& e){
            throw ModelException(e, "CAssetCVB::getMarket");
        }
    }
 
    /** stop tweaking for empty equity */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const{
        if (field->getName() == "equity"){
            if(MuParallel::Shift::TYPE->isAssignableFrom(targetClass)   ||
               MuPointwise::IShift::TYPE->isAssignableFrom(targetClass) ||
               MuSpecial::IShift::TYPE->isAssignableFrom(targetClass)){
                // the divs are not real - therefore hide all greeks and
                // scenarios dealing with divs
                return false;
            }
            if (useCvbTheoPrice && 
                (ITweakableWithRespectTo<Spot>::TYPE->isAssignableFrom(targetClass) ||
                 DeltaSurface::Shift::TYPE->isAssignableFrom(targetClass))){
                // for theo, spot price is derived therefore hide delta
                return false;
            }
            // can shift everything else
            return true;
        }
        if (field->getName() == "CVBData"){
            if (!useCvbTheoPrice){
                // never tweak CVBData when doing mtm
                return false;
            }
            // using cvb theo price - must therefore recalculate it if it has
            // been shifted
            // Below is a bit of a hack - we need a common base class for
            // all the shift interfaces to see if out target class is a shift
            if (!DerivativeAsset::TYPE->isAssignableFrom(targetClass) &&
                !FXAsset::TYPE->isAssignableFrom(targetClass)){
                resetEquitySpot = true;
            }
        }
        return true; // tweak everything else as needed
    }
    
    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        asset together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. Note that the struckAsset is indeed the struckAsset cf
        'this' which is the non struck asset */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      struckAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const{
        if (!vol) {
            throw ModelException("CAssetCVB::getProcessedVol",
                                 "AssetCVB does not have a vol object.");
        }
        return vol->getProcessedVol(volRequest, struckAsset,
                                    fxAsset, eqFXCorr);
    }

    /** sets whether the asset should use the theoretical price or the
        mtm price */
    virtual void setUseTheoAssetPrice(bool useTheoAssetPrice){
        if (!equity || useCvbTheoPrice != useTheoAssetPrice){
            useCvbTheoPrice = useTheoAssetPrice;
            resetEquitySpot = true;
            // this is called once we block of pricings so now is a good time
            // to cache the value of cvb spot.
            checkEquity();
        }
    }

    /** sets the control ptr in the asset to the inputted
        control. This method is invoked on each DerivativeAsset just
        before a pricing call is made to the instrument. The supplied
        control is guaranteed to be valid for the duration of the
        pricing call. */
    virtual void setControlAndModel(const CControlSP& ctrl,
                                    const IModelSP&   model){
        CVBModel = model; /* ensure we use the right instance (can vary
                             across greeks eg Delta T+1) */
        if (!control && !resetEquitySpot){
            // control will be null before first pricing. If resetEquitySpot
            // is false then we have done a pricing already (in roll to now)
            // using an empty control. So force a reprice with full control
            resetEquitySpot = true;
        }
        control = ctrl;
    }

    /* Returns the model used to calculate asset's theo price. This is used
       to allow the DerivativeAsset model to manage the model. See 
       setControlAndModel(). */
    virtual IModelSP getModel() const{
        return CVBModel;
    }

    /* Returns the instrument used to calculate asset's theo price. */
    virtual CInstrumentSP getInstrument() {
        return CVBData;
    }

    /* Returns the MTM price of the derivative asset if it exists */
    virtual double getMTM() const {
        if (!SpotPrice) {
            throw ModelException("CAssetCVB::getMTM", "The MTM price has not been set for CVB asset" + name);
        }
        return *SpotPrice;
    }

    /** Returns the name of the stock/asset - used to determine
        whether to shift the object */
    virtual string sensName(SpotLevelProbability* shift) const{
        return name;
    }

    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(SpotLevelProbability* shift){
        checkEquity();
        double spot = shift->spotLevel(borrowCurve->getBaseDate(), this);
        SpotLevel level(spot);
        equity->sensShift(&level);  // set equity spot to our new level
        return false;  // all done;
    }

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewParallel* shift) const{
        if (!!vol) {
            return vol->getName();
        } else {
            return "";
        }
    }

    /** Shifts the object for VegaSkewParallel */
    virtual bool sensShift(VegaSkewParallel* shift){
        checkEquity();
        // store the underlying spot price
        shift->setSpot(equity->spot());
        return true; // continue on and do the vol
    }

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewPointwise* shift) const{
        if (!!vol) {
            return vol->getName();
        } else {
            return "";
        }
    }
    
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const{
        return ExpiryArrayConstSP(); // return null SP - the vol will do this
    }

    virtual bool sensShift(VegaSkewPointwise* shift){
        checkEquity();
        // store the underlying spot price
        shift->setSpot(equity->spot());
        return true; // continue on and do the vol
    }

    /** Returns the name of the asset - used to determine whether to tweak
        the object */
    virtual string sensName(DeltaSurface* shift) const{
        if (vol.get()) {
            // only bother if we have a vol, but report against asset name
            return name;
        }
        return "";
    }
    
    virtual bool sensShift(DeltaSurface* shift){
        checkEquity();
        // store the underlying spot price
        shift->setSpot(equity->spot(), name, vol->getName());
        return true; // continue on and do the vol
    }

    string sensName(VolRelativeShift* shift) const{
        // only bother if we have a vol, but report against asset name
        if (vol.get() && 
            VolRelativeShift::IShift::TYPE->isInstance(vol.get())) {
            return name;
        }
        return "";            
    }
    
    bool sensShift(VolRelativeShift* shift){
        // store the underlying spot price
        shift->setSpot(equity->getValueDate(), this);
        return true; // continue on and do the vol
    }

    /** Returns fair value of stock price */
    virtual double fairValue() const {
        checkEquity();
        return equity->fairValue();
    }

    /** given a current spot level, get the next strike on the vol
        surface where the slope is non-differentiable */
    double getNextStrike(const double& strike,
                         bool          isUp,
                         bool&         offSurface) const {
        double volStrike;
        if (!vol) {
            volStrike  = 0.0;
            offSurface = true;
        } else {
            const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(vol.get());
            if (nextStrike) {
                volStrike = nextStrike->getNextStrike(strike,
                                                      isUp,
                                                      offSurface);
            } else {
                volStrike  = 0.0;
                offSurface = true;
            }
        }
        return volStrike;
    }

    //// public just to shut the compiler up
    CAssetCVB(): CAsset(TYPE), SpotPrice(0), cvbDeltaOverride(0), FixOnCleanPrice(false),
                 resetEquitySpot(true), useCvbTheoPrice(true)  {
        // default our settlement to T + 3
        HolidayConstSP hol(Holiday::weekendsOnly());
        settlement = SettlementSP(new RollingSettlement(3, hol));
    }

    ~CAssetCVB(){
        delete SpotPrice;
    }

private:
    CAssetCVB(const CAssetCVB& rhs);
    CAssetCVB& operator=(const CAssetCVB& rhs);



    /** Validates that no-one is asking for a forward price beyond the length
        of the bond */
    void checkFwdDate(const DateTimeArray& dateList) const {
        if (!dateList.empty()){
            checkFwdDate(dateList.back());
        }
    }

    /** Validates that no-one is asking for a forward price beyond the length
        of the bond */
    void checkFwdDate(const DateTime& date) const {
        if (date.isGreater(CVBData->getBondMaturityDate())){
            throw ModelException("AssetCVB::checkFwdDate",
                                 "requested bond price at date " +
                                 date.toString() +
                                 " which is after bond maturity of " + 
                                 CVBData->getBondMaturityDate().toString());
        }
    }

    /** Builds equity using the mtm value for the spot price */
    void buildEquity(){
        const string method = "AssetCVB::buildEquity";
        try{
            // build equity - which we use to handle fwd
            // price calc create dividend list first
            DateTime valDate = borrowCurve->getBaseDate();
            DateTime firstDate = valDate;
            DateTimeArraySP exDates = CVBData->getExCouponDates();
            if (exDates->size()>0) {
                firstDate = (*exDates)[0].rollDate(-1);
            }
            CashFlowArraySP cf = CVBData->getCoupons(firstDate);
            DividendArray divArr(cf->size());
            if (cf->size() < exDates->size()){
                throw ModelException(method, "ex-coupon dates and "
                                     "coupons have different array sizes");
            }

            double mtmPrice = (!SpotPrice)?100.0:getMTM();
            for (int i = 0; i < cf->size(); i++) {
                DateTime divDate(((*exDates)[i].getDate() <= (*cf)[i].date.getDate())?(*exDates)[i].getDate():(*cf)[i].date.getDate(), 
                                 Dividend::DIVIDEND_EXDIV_TIME);
                divArr[i] = Dividend(divDate,
                                     (*cf)[i].date, 
                                     Dividend::AMOUNT,
                                     (*cf)[i].amount);
            }
            DividendList divList(divArr);
            // MAR: futher work needs to be done regarding fitting the
            // cvb's price into an equity.
            const DateTime& stockDate = valDate; 
            equity = EquitySP(new Equity(name, valDate, stockDate,
                                         mtmPrice, settlement.get(),
                                         CVBData->getDiscount().get(),
                                         &divList, borrowCurve.get()));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
   
    void checkEquity() const{
        if (!equity){
            throw ModelException("AssetCVB::checkEquity", "Null equity");
        }
        if (!SpotPrice && useCvbTheoPrice) {
            throw ModelException("AssetCVB::checkEquity", "The MTM price has not been set for CVB asset" + name);
        }

        if (resetEquitySpot){
            resetEquity();
            resetEquitySpot = false;
        }
    }

    /** calc spot by calling cvb model */
    double priceCVB() const{
        // call cvb model
        CResults   results;
        if (!control){
            // we're yet to start the pricing run (probably in roll to now)
            // so price with empty control
            Control ctrl(SensitivityArrayConstSP(new SensitivityArray()),
                         OutputRequestArrayConstSP(new OutputRequestArray()),
                         false, "");
            ctrl.calculate(CVBModel.get(), CVBData.get(), &results);
        } else {
            CVBModel->Price(CVBData.get(), control.get(), &results);
        }
        if (!FixOnCleanPrice) {
            return results.retrievePrice();
        } else {
            double cleanPrice = results.retrievePrice() - 
                CVBData->getAccruedAtDate(borrowCurve->getBaseDate());
            return cleanPrice;
        }
    }
        
    //// reset cvb spot
    void resetEquity() const {
        const string method = "AssetCVB::resetEquity";
        try{
            if (!useCvbTheoPrice){
                // using mtm - nothing to do. Equity has spot price in it
                // and it will have tweaked it as necessary
            } else {
                // if the equity already exists then just change
                // the spot Note: we are relying on the Equity to
                // do the shifting for the greeks so we must not
                // rebuild it
                double cvbSpot = priceCVB();
                SpotLevel spotLevel(cvbSpot);
                equity->sensShift(&spotLevel);
            }
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    static IObject* defaultAssetCVB(){
        return new CAssetCVB();
    }

    /** Validates the uniqueness of asset names. */
    static void acceptNameCollector(const CAssetCVB*    asset,
                                    AssetNameCollector* collector){
        collector->assetNameValidate(asset->getName());
        // may want to check cvb data etc
    }

    /** validates cross validate start dates/value dates */
    static void acceptValueDateCollector(const CAssetCVB*     asset, 
                                         CValueDateCollector* collector) {
        asset->equity->accept(collector);
        asset->CVBData->accept(collector); // might need a method on cvb
    }

    /** helps determine a good delta shift size */
    static void acceptDeltaShift(const CAssetCVB*    asset, 
                                 ShiftSizeCollector* collector){
//        if (!!vol) {
            const IObject* obj = asset->vol.get();
            const IVolDeltaShiftSize* volDelShift = 
                dynamic_cast<const IVolDeltaShiftSize*>(obj);
            if (volDelShift) {
                volDelShift->adjustDeltaShiftSize(collector,
                                           asset->getName(),
                                           asset->getSpot());
            }
//        }
    }

    /** validate currencies are consistent */
    static void acceptImntCcy(const CAssetCVB*   asset,
                              AssetCcyCollector* collector){
        collector->currencyValidate(asset->equity->getYCIsoCode(),
                                    asset->getName());
    }

    /** retrieves the market holiday from an asset */
    static void acceptHoliday(const CAssetCVB*  asset,
                              HolidayCollector* collector){
        collector->setHoliday(asset->equity->getMarketHolidays());
    }

    /** dollar div and yield div converted to dollar, yield converted
        to amount but the type name is kept as PERCENET for dollar div
        treatment.  better to add a new type AMOUNT_FROM_YIELD in
        dividend */
    static void acceptCriticalDateCollector(const CAssetCVB*       asset, 
                                            CriticalDateCollector* collector){
        asset->equity->accept(collector);
    }
    
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CAssetCVB, clazz);
        SUPERCLASS(CAsset);
        IMPLEMENTS(DerivativeAsset);
        IMPLEMENTS(CAsset::IStruckable);
        IMPLEMENTS(SpotLevelProbability::Shift);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(ObjectIteration::IOverride);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        EMPTY_SHELL_METHOD(defaultAssetCVB);
        FIELD(name, "asset name");
        FIELD(SpotPrice, "spot price of cvb (not stock)");
        FIELD_MAKE_OPTIONAL(SpotPrice);
        FIELD(cvbDeltaOverride, 
                     "delta override for cvb. Not supported yet.");
        FIELD_MAKE_OPTIONAL(cvbDeltaOverride);
        FIELD(vol, "volatility of the convertible bond (not stock)");
        FIELD_MAKE_OPTIONAL(vol);
        FIELD(borrowCurve, "cost of borrowing cvb (not stock)");
        FIELD(CVBData, "cvb instrument data");
        FIELD(CVBModel, "model used to price the cvb");
        FIELD(FixOnCleanPrice, "true = asset fixing (spot) price is "
                     "clean price, "
                     "false(default) = spot fixing on dirty price");
        FIELD_MAKE_OPTIONAL(FixOnCleanPrice);
        
        FIELD(settlement, "cvb asset settlement. Default to T+3 with"
              " weekends as holidays.");
        FIELD_MAKE_OPTIONAL(settlement);
        // hide equity but registered for tweaking and clone
        FIELD_NO_DESC(equity);
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(equity);
        FIELD_NO_DESC(resetEquitySpot);
        FIELD_MAKE_TRANSIENT(resetEquitySpot);
        FIELD_NO_DESC(useCvbTheoPrice);
        FIELD_MAKE_TRANSIENT(useCvbTheoPrice);
        ClassSetAcceptMethod(acceptNameCollector);
        ClassSetAcceptMethod(acceptValueDateCollector);
        ClassSetAcceptMethod(acceptDeltaShift);
        ClassSetAcceptMethod(acceptImntCcy);
        ClassSetAcceptMethod(acceptHoliday);
        ClassSetAcceptMethod(acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our 'dividends' */
    static void acceptDividendCollector(const CAssetCVB*   asset,
                                        DividendCollector* collector){
        asset->checkEquity();
        collector->addDivs(asset->equity->getDivList(), 0 /* not struck */);
    }
};

CClassConstSP const CAssetCVB::TYPE = CClass::registerClassLoadMethod(
    "AssetCVB", typeid(CAssetCVB), load);

// * for class loading (avoid having header file) */
bool AssetCVBLoad() {
    return (CAssetCVB::TYPE != 0);
}

DRLIB_END_NAMESPACE


