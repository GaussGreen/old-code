#include "edginc/config.hpp"
#include "edginc/DerivativeAsset.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/AllStrikes.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FutureAsAsset.hpp"
#include "edginc/BorrowCurve.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/Equity.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Spot.hpp"
#include "edginc/VolTypeSensitiveStrikes.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE
/** Makes an instrument look like an asset. The volatility is provided whilst
    the forward is caclulated using a risk-neutral approach. */
class FutureAsset: public CAsset,
                   public virtual DerivativeAsset {
private:
    /// fields ////
    string                     name;
    CVolBaseWrapper            vol;        // volatility of instrument
    IFutureAsAssetSP           inst;       // the instrument
    IModelSP                   model;      // model for pricing instrument
    SettlementSP               settlement; // instrument settlement
  
    // transient fields ///
    CControlSP                 control; // not registered $unregistered

public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method("FutureAsset::validatePop2Object()");
        if (name.empty()){
            throw ModelException(method, "FutureAsset name is empty.");
        }
        if (!CInstrument::TYPE->isInstance(inst)){
            // this is only because the infrastructure currently works with
            // CInstrument. It needs upgrading to work with IInstrument
            throw ModelException(method, "Instrument doesn't derive from "
                                 "CInstrument. Contact DR");
        }
    }

    /** clone() has to keep the same model for same grid tweaking */
    IObject* clone() const {
        // copy registered fields first
        IObject*  obj = CAsset::clone();
        // shallow copy model to allow same grid tweaking
        FutureAsset& inst = dynamic_cast<FutureAsset&>(*obj);
        inst.model = model;
        inst.control = control; // shallow copy over
        return obj;
    }

    /** returns the spot price */
    virtual double getSpot() const{
        return priceInstrument();
    }
    
    /** returns the asset name */
    virtual string getName() const{
        return name;
    }
    
    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const{
        return inst->discountYieldCurveName();
    }
    
    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const {
        // ISS - don't kow what to do here
        return 0.0;
    }

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const {
        // ISS - don't kow what to do here
        return true;
    }

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const {
        // ISS - don't kow what to do here
        return 0.0;
    }
    
    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const {
        // ISS - don't kow what to do here
        return false;
    }

    /** Calculates the expected spot price of the asset at the given date.
        Do not use this repeatedly to calculate values over a set of
        dates (poor performance) - instead use other fwdValue method */
    virtual double fwdValue(const DateTime& date) const{
        if (date < inst->maturityDate()) {
            return getSpot();
        }
        else {
            return 0.0;
        }
    }
    
    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const{
        result = DoubleArray(dateList.size(),0);
        
        for(int iDate=0; iDate<dateList.size(); iDate++) {
            if (dateList[iDate] < inst->maturityDate()) {
                result[iDate] = getSpot();
            }
        }
    }
        
    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const{
        if (fwdDate < inst->maturityDate()) {
            return spot;
        }
        else {
            return 0.0;
        }
    }
    
    /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const{
        return settlement->settles(tradeDate);
    }
    

    /** for Vega Matrix sensitivity */
    virtual void getSensitiveStrikes(
        const CVolRequest* volRequest,
        OutputNameConstSP outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP sensitiveStrikes) const{
        
        if (!!vol) {
            const CVolRequestLN& request = 
                dynamic_cast<const CVolRequestLN&>(*volRequest);
            
            if (sensStrikeDesc.forwardOnly == false && 
                outputName->equals(vol->getName())) {
                request.getSensitiveStrike(getSpot(), sensitiveStrikes);
            }  
            if (IVolTypeSensitiveStrikes::TYPE->isInstance(vol.get())) {
                // may have to defer to the vol
                const IVolTypeSensitiveStrikes* myVol = 
                        dynamic_cast<const IVolTypeSensitiveStrikes*>(vol.get());
                myVol->getSensitiveStrikes(this, volRequest, outputName, 
                                            sensStrikeDesc, sensitiveStrikes);
            }
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
            throw ModelException("FutureAsset::pdfCalculator",
                                 "vol of type (" + vol->getClass()->getName() +
                                 ") has no pdf calculator");
        }
        else {
            throw ModelException("FutureAsset::getProcessedVol",
                                 "FutureAsset does not have a vol object.");
        }
    }

    /** get data from market */
    virtual void getMarket(const IModel* extModel, const MarketData* market){
        try{
            // call instrument getMarket to fill inst data and stock data
            // underlying the instrument
            const CMarketDataSP src = CMarketDataSP::attachToRef((MarketData*)market);
            
            inst->GetMarket(model.get(), src);
            
            settlement->getMarket(model.get(), market);
            // fill vol 
            if (vol.getName() != "") {
                vol.getData(extModel, market);
            }
            
            // Validate the instrument
            inst->Validate();

        } 
        catch (exception& e){
            throw ModelException(e, "FutureAsset::getMarket");
        }
    }
    
    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const{
        if (!vol) {
            throw ModelException("FutureAsset::getProcessedVol",
                                 "FutureAsset does not have a vol object.");
        }
        return vol->getProcessedVol(volRequest, this);
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
            throw ModelException("FutureAsset::getProcessedVol",
                                 "FutureAsset does not have a vol object.");
        }
        return vol->getProcessedVol(volRequest, struckAsset,
                                    fxAsset, eqFXCorr);
    }

    /*------------------ derivative asset interface --------------------------------*/

    /** sets whether the asset should use the theoretical price or the
        mtm price */
    virtual void setUseTheoAssetPrice(bool useTheoAssetPrice){
        
    }

    /** sets the control ptr in the asset to the inputted
        control. This method is invoked on each DerivativeAsset just
        before a pricing call is made to the instrument. The supplied
        control is guaranteed to be valid for the duration of the
        pricing call. */
    virtual void setControlAndModel(const CControlSP& ctrl,
                                    const IModelSP&   model){
        this->model = model; /* ensure we use the right instance (can vary
                                across greeks eg Delta T+1) */
        control = ctrl;
    }

    /* Returns the model used to calculate asset's theo price. This is used
       to allow the DerivativeAsset model to manage the model. See 
       setControlAndModel(). */
    virtual IModelSP getModel() const{
        return model;
    }

    /* Returns the instrument used to calculate asset's theo price. */
    virtual CInstrumentSP getInstrument() {
        CInstrumentSP theInst(CInstrumentSP::dynamicCast(inst));
        return theInst;
    }

    /* Returns the MTM price of the derivative asset if it exists */
    virtual double getMTM() const {
       
        return getSpot();
    }

    /*---------------------------------------------------------------------------------*/

    // default settlement T+3 with week-end only
    FutureAsset(): CAsset(TYPE)   {
        HolidayConstSP hol(Holiday::weekendsOnly());
        settlement = SettlementSP(new RollingSettlement(3, hol));
    }

    ~FutureAsset(){
    }

private:
    FutureAsset(const FutureAsset& rhs);
    FutureAsset& operator=(const FutureAsset& rhs);
 
    static IObject* defaultFutureAsset(){
        return new FutureAsset();
    }

    /** calc spot by calling instrument */
    double priceInstrument() const{
        CInstrumentSP theInst(CInstrumentSP::dynamicCast(inst));
        CResults   results;
        if (!control){
            // we're yet to start the pricing run (probably in roll to now)
            // so price with empty control
            Control ctrl(SensitivityArrayConstSP(new SensitivityArray()),
                         OutputRequestArrayConstSP(new OutputRequestArray()),
                         false, "");
            ctrl.calculate(model.get(), theInst.get(), &results);
        } else {
            model->Price(theInst.get(), control.get(), &results);
        }
        double price = results.retrievePrice();
      
        return price;
    }

    /** Validates the uniqueness of asset names. */
    static void acceptNameCollector(const FutureAsset*   asset,
                                    AssetNameCollector* collector){
        collector->assetNameValidate(asset->getName());
        // may want to check instrument data etc
    }

    /** validates cross validate start dates/value dates */
    static void acceptValueDateCollector(const FutureAsset* asset, 
                                         CValueDateCollector*   collector) {
      
        asset->inst->accept(collector); // might need a method on instrument
    }

    /** helps determine a good delta shift size */
    static void acceptDeltaShift(const FutureAsset*           asset, 
                                 ShiftSizeCollector*        collector){
        if (asset->vol.get()){
            const IVolDeltaShiftSize* volDelShift = 
                       dynamic_cast<const IVolDeltaShiftSize*>(asset->vol.get());
            if (volDelShift) {
                volDelShift->adjustDeltaShiftSize(collector,
                                            asset->getName(),
                                            asset->getSpot());
            }
        }
    }

    /** validate currencies are consistent */
    static void acceptImntCcy(const FutureAsset*          asset,
                              AssetCcyCollector*         collector){
        collector->currencyValidate(asset->inst->getDiscount()->getCcy(),
                                    asset->getName());
    }

    /** retrieves the market holiday from an asset */
    static void acceptHoliday(const FutureAsset*         asset,
                              HolidayCollector* collector){
        collector->setHoliday(asset->settlement->getMarketHolidays());
    }
    
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(FutureAsset, clazz);
        SUPERCLASS(CAsset);
        IMPLEMENTS(DerivativeAsset);
        EMPTY_SHELL_METHOD(defaultFutureAsset);
        FIELD(name, "asset name");
        FIELD(vol, "volatility of the asset (not stock)");
        FIELD_MAKE_OPTIONAL(vol);
        FIELD(inst, "instrument data within the asset");
        FIELD(model, "model used to price the asset");
        FIELD(settlement, "asset settlement. Default to T+3 with weekends as holidays.");
        FIELD_MAKE_OPTIONAL(settlement);
        ClassSetAcceptMethod(acceptNameCollector);
        ClassSetAcceptMethod(acceptValueDateCollector);
        ClassSetAcceptMethod(acceptDeltaShift);
        ClassSetAcceptMethod(acceptImntCcy);
        ClassSetAcceptMethod(acceptHoliday);
    }
};

CClassConstSP const FutureAsset::TYPE = CClass::registerClassLoadMethod(
    "FutureAsset", typeid(FutureAsset), load);

// * for class loading (avoid having header file) */
bool FutureAssetLoad() {
    return (FutureAsset::TYPE != 0);
}

IFutureAsAsset::~IFutureAsAsset(){}
IFutureAsAsset::IFutureAsAsset(){}

void IFutureAsAsset::load(CClassSP& clazz){
    REGISTER_INTERFACE(IFutureAsAsset, clazz);
    EXTENDS(IInstrument);
}

CClassConstSP const IFutureAsAsset::TYPE =
CClass::registerInterfaceLoadMethod(
    "IFutureAsAsset", typeid(IFutureAsAsset), load);
    

DRLIB_END_NAMESPACE


