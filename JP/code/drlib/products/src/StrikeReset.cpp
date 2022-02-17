//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : StrikeReset.cpp
//
//   Description : Strike reset option.
//                 the strike = max or min of (k0, k1)
//                 where k0 is the initial strike, k1 is moneyness at strike reset date.
//                 We only allowed one reset date for easy vol strike estimate in LN MC.
//                 Note on volStrike: this is the strike for vol interp in LN MC. 
//                 computed as expected value (a call/put option itself).
//
//
//   Date        : 15 May 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/Maths.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/Average.hpp"
#include "edginc/CashSettlePeriod.hpp"

DRLIB_BEGIN_NAMESPACE

/** Striek reset model */
//////// instrument class //////////
class StrikeReset: public Generic1Factor, 
                   virtual public LastSensDate,
                   virtual public IMCIntoProduct{
protected:
    /// fields ////////
    SampleListSP            spotSamples; // average out
    SampleListSP            strikeSamples; // ave in for strike reset
    double                  initStrike; // intial strike
    double                  resetScaling;  // reset strike = resetScaling * avg of strike samples
    bool                    isResetUp; // true max of strike, false = min
    bool                    isCall;

public:
    static CClassConstSP const TYPE;
    friend class StrikeResetProd;

    virtual void Validate(){
        static const string routine("StrikeReset::Validate");
        // NB don't call parent's validate - issues with fwdStart
        if (fwdStarting){
            throw ModelException(routine, "Fwd starting flag on "
                                 "Generic1Factor is not used");
        }

        if (ccyTreatment != CAsset::CCY_TREATMENT_NONE &&
            ccyTreatment != CAsset::CCY_TREATMENT_VANILLA) {
            throw ModelException(routine, "ccy struck or protected is not supported.");
        }

        if (!Maths::isPositive(resetScaling)) {
            throw ModelException(routine, "resetScaling must be positive.");
        }


        AssetUtil::assetCrossValidate(asset.get(),
                                      false, //fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** when to stop tweaking - default implementation assumes product can
        be priced with a MC */
    DateTime endDate(const Sensitivity* sensControl) const{
        DateTime matDate = spotSamples->getDates()[spotSamples->getDates().size()-1];
        DateTime instEnd  = instSettle->settles(matDate, asset.get());
        DateTime assetEnd = asset->settleDate(matDate);
        DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;
    }

private:
    StrikeReset(): Generic1Factor(TYPE){}

    // for reflection
    StrikeReset(const StrikeReset& rhs); // not implemented
    StrikeReset& operator=(const StrikeReset& rhs); // not implemented

    static IObject* defaultStrikeReset(){
        return new StrikeReset();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta){
        // use valueDate before it changes
        spotSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
                         asset.get()); // then roll our past values

        strikeSamples->roll(theta->getUtil(valueDate), 0 /* iAsset */,
                         asset.get()); // then roll our past values

        Generic1Factor::sensShift(theta); // and then call parent's method
        return true; // continue to tweak components which implement Theta
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StrikeReset, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultStrikeReset);
        FIELD(spotSamples, "average out");
        FIELD(strikeSamples, "ave in for strike reset");
        FIELD(initStrike,    "intial strike");
        FIELD(resetScaling,    "scaling applied to reset strike samples");
        FIELD(isResetUp,    "true = strike will be reset up, false = strike will be reset down");
        FIELD(isCall,    "is it a call option");

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

/* MC product class for super rainbow */
//////// product class //////////
class StrikeResetProd : public IMCProduct, virtual public IMCProductLN{
private:
    const StrikeReset*   inst; // reference to original instrument
    IRefLevel::IMCPathSP    avgOutSample;
    double                  callPut;
    int                     numDates;
    mutable double          volStrike;

public:
    
    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    StrikeResetProd(const StrikeReset*         inst,
                const IRefLevelConstSP&   refLevel, // how to 'avg in'IRefLevelSP refLevel,
                const SimSeriesSP&             simSeries,
                 const IPastValuesConstSP& mcPastValues)  : // past values
                IMCProduct(inst->asset.get(),
                          inst->valueDate,
                          inst->discount.get(),
                          refLevel,
                          simSeries,
                          inst->spotSamples,
                          inst->instSettle.get(),
                          simSeries->getLastDate()),
                          inst(inst){

                    // compute expected strike
                    numDates = inst->spotSamples->numDates(0);
                    avgOutSample = 
                        IRefLevel::IMCPathSP(
                            inst->spotSamples->createMCPath(inst->valueDate,
                                                       DoubleArray(1), /* value is
                                                                          irrelevant */
                                                       inst->spotSamples.get()));
                }
        
    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,   // (I)
                IMCPrices&                prices) {  // (O)
        // only do future (or rather avoid doing partial past)
        if (pathGen->end(0) == numDates){
            // first start with 'in' value
            double inValue = pathGen->refLevel(0,0);
            // then do out value
            const double* path = pathGen->Path(0,0); // access path
            path += pathGen->begin(0); // skip over past
            double outValue = avgOutSample->refLevel(0, path);

            // compute strike

            double strike = inst->isResetUp? 
                Maths::max(inValue * inst->resetScaling, inst->initStrike) :
                Maths::min(inValue * inst->resetScaling, inst->initStrike);
            double price = outValue - strike;
            if (!inst->isCall){
                price = -price;
            }
            prices.add(price < 0.0? 0.0: price);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        CVolRequestLNArray   reqarr(1); // one interp level/path per asset here

        DateTime matDate = inst->strikeSamples->getDates()[inst->strikeSamples->getDates().size()-1];

        // we compute expected strike for vol interp
        // vol strike = initStrike + resetScaling * 
        //                   AverageModel(strikeSamples, initStrike/resetScaling)/PV(0, lastStrikeSampleDate)

        CResults results;

        const HolidayConstSP hols(Holiday::noHolidays());
        CashSettlePeriod sett(1);
        AverageSP avg(Average::makeAvgSpot(inst->isResetUp,
                                       matDate,
                                       inst->initStrike/inst->resetScaling,
                                       inst->strikeSamples.get(),
                                       &sett,
                                       inst->premiumSettle.get(),
                                       inst->asset.get(),
                                       inst->ccyTreatment,
                                       inst->discount.get(),
                                       inst->valueDate,
                                       false,
                                       inst->startDate,
                                       true,
                                       inst->notional,
                                       inst->initialSpot));
         CClosedFormLN model;
         model.Price(avg.get(),0, &results);
         double sgn = inst->isResetUp? 1.0 : -1.0;
         volStrike = inst->initStrike + sgn*inst->resetScaling*results.retrievePrice()
                        /inst->discount->pv(matDate);

        reqarr[0] = CVolRequestLNSP(new LinearStrikeVolRequest(
            volStrike,
            inst->valueDate, 
            matDate,
            false));
        
        return reqarr;
    }

    // extra output
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices& ) const
    {
        // store vol strike
        if (control && control->isPricing())
        results->storeScalarGreek(volStrike, Results::DEBUG_PACKET, 
                                      OutputNameSP(new OutputName("VOL_STRIKE")));
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* StrikeReset::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(spotSamples->getAllDates());

    // v simple RefLevel
    IRefLevelConstSP refLevel(IRefLevelConstSP::attachToRef(strikeSamples.get())); 
 
    // similarly for levels
    IPastValuesConstSP pastValues(
        IPastValuesConstSP::attachToRef(spotSamples.get()));
    return new StrikeResetProd(this, refLevel, simSeries, pastValues);
}

CClassConstSP const StrikeReset::TYPE = CClass::registerClassLoadMethod(
    "StrikeReset", typeid(StrikeReset), StrikeReset::load);

// force linker to include this file (avoid having header file) */
bool StrikeResetLoad() {
    return (StrikeReset::TYPE != 0);
}

DRLIB_END_NAMESPACE
