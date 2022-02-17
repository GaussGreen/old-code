//----------------------------------------------------------------------------
//
//   Description : statistics products for monte carlo 
//
//   Author      : 
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/GenericNFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenExpectedEnergyFuture.hpp"
#include "edginc/MonteCarlo.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

// for class loading 
bool MCEnergyFutureLoad();
bool MCEnergyMultiFutureLoad();
bool MCStatisticsProductsLoad() {
    return (
        MCEnergyFutureLoad() && 
        MCEnergyMultiFutureLoad()
        );
}


class MCEnergyFuture: public CInstrument,
    public virtual IMCIntoProduct{

public:
    static CClassConstSP const TYPE; 

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const; // see below

    virtual void GetMarket(const IModel* model, const CMarketDataSP market)
    {
        // if we were to put tie an asset with this product, we would do something like
        futureCurveWrapper.getData(model, market);
        yieldCurveWrapper.getData(model, market);        

        if (!fxWrapper.isEmpty())
            fxWrapper.getData(model, market);
        if (!parentCurveWrapper.isEmpty())
            parentCurveWrapper.getData(model, market);
    }

    virtual void Validate() 
    {
        const string & method = "MCEnergyFuture::Validate";
        EnergyFuturesCurveConstSP futureCurve = futureCurveWrapper.getSP();
        EnergyUnderlyerConstSP underlyer = futureCurve->getEnergyUnderlyer();

        // construct the expiry date array of the future curve:
        //DateTimeArray futureExpiryDates;
        //for (size_t i = 0; i < futureCurve->getSize(); ++i) {
        //    EnergyContractLabel currExpiryLabel(futureCurve->getExpiryLabel(i));
        //    DateTime currExpiryDate = underlyer->expiryDate(currExpiryLabel);
        //    futureExpiryDates.push_back(currExpiryDate);
        //}
        const DateTimeArray& futureExpiryDates = futureCurve->getFutureMaturityDates();

        // check expiry date:
        int expiryIdx = maturityDate.find(futureExpiryDates); 

        // check if fx is needed:
        if (discountYieldCurveName() != underlyer->getPricingCurrency() &&
            fxWrapper.isEmpty())
            throw ModelException(method, "Product domestic ccy is different "
            "from the energy curve pricing ccy --> expect fx but none is "
            "specified by the product");
    }

    void validatePop2Object()
    {
        const string & method = "MCEnergyFuture::validatePop2Object";
        if (maturityDate < measureDate)
            throw ModelException(method, "Observation date is before expiry date!");
    }

    virtual DateTime getValueDate() const { return DateTime(); };

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const {
        return yieldCurveWrapper.getName();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MCEnergyFuture, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultMCEnergyFuture);
        FIELD(futureCurveWrapper, "future curve wrapper");
        FIELD(parentCurveWrapper, "parent future curve wrapper");
        FIELD_MAKE_OPTIONAL(parentCurveWrapper);
        FIELD(yieldCurveWrapper, "yield curve wrapper");
        FIELD(fxWrapper, "fx wrapper");
        FIELD_MAKE_OPTIONAL(fxWrapper);
        FIELD(today, "today");
        FIELD(maturityDate, "maturity date of the future");
        FIELD(measureDate, "measurement date");
        FIELD(doLog, "compute log");
        FIELD(doSquare, "compute price squared");
        FIELD_MAKE_OPTIONAL(doSquare);
        FIELD(dumpPrices, "dump prices to text file");
        FIELD_MAKE_OPTIONAL(dumpPrices);
        FIELD(dumpFileName, "dump file name");
        FIELD_MAKE_OPTIONAL(dumpFileName);
    }

    static IObject* defaultMCEnergyFuture() {
        return new MCEnergyFuture();
    }

private:
    MCEnergyFuture() : 
       CInstrument(TYPE), 
           doSquare(false),
           dumpPrices(false),
           dumpFileName("c:/debugEnergySRM3.txt")
       {}; 

       class MC;
       friend class MC;

       EnergyFuturesCurveWrapper futureCurveWrapper;
       EnergyFuturesCurveWrapper parentCurveWrapper;
       YieldCurveWrapper yieldCurveWrapper;
       FXAssetWrapper fxWrapper;
       DateTime today; 
       DateTime maturityDate; // maturity of the future price that we want to observe
       DateTime measureDate; // measurement date
       bool doLog; // compute log
       bool doSquare; // compute price squared
       bool dumpPrices;
       string dumpFileName;
    };

CClassConstSP const MCEnergyFuture::TYPE = CClass::registerClassLoadMethod(
    "MCEnergyFuture", typeid(MCEnergyFuture), MCEnergyFuture::load);

// for class loading 
bool MCEnergyFutureLoad() {
    return (MCEnergyFuture::TYPE != 0);
}


/* MC view of the MCEnergyFuture product */
class MCEnergyFuture::MC : public MCProductClient /*,
                                                  public IMCStatelessProductClient*/
{
private:
    const MCEnergyFuture* inst; // reference ptr to the original inst

    // state var and state var generators
    SVExpEnergyFutureSP expFp; // future price SV
    SVGenExpectedEnergyFutureSP expFpGen; // future price SV generator
    SVExpEnergyFutureSP expFpParent; // parent SV
    SVGenExpectedEnergyFutureSP expFpParentGen; // parent SV generator

    SVDiscFactorSP domDf; // df state variable
    SVGenDiscFactorSP domDfGen; // df SV generator
    SVDiscFactorSP forDf; 
    SVGenDiscFactorSP forDfGen; 
    MCPath::IStateVarSP fx;
    SVGenSpotSP fxGen;

private:
    DateTimeArray maturityDates;
    long payoffCounter;
    bool needFx;
    bool isT2;
    ofstream dumpFile;

protected:
    /** Override default method on IMCProduct. This method is called every time
    the path generator is changed (which is, at the moment, when the
    past path generator is created, and then when the future path
    generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) 
    {
        expFp = expFpGen->getSVExpEnergyFuture(newPathGen); // create SV
        domDf = domDfGen->getSVDiscFactor(newPathGen);

        if (needFx) {
            forDf = forDfGen->getSVDiscFactor(newPathGen);
            fx = fxGen->getSpotSV(newPathGen);
        }
        if (isT2) {
            expFpParent = expFpParentGen->getSVExpEnergyFuture(newPathGen);
        }
    }

public:
    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const 
    {
        svCollector->append(expFpGen.get()); // collect SV generator
        svCollector->append(domDfGen.get());

        if (needFx) {
            svCollector->append(forDfGen.get());
            svCollector->append(fxGen.get());
        }
        if (isT2) {
            svCollector->append(expFpParentGen.get());
        }
    }

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC( const MCEnergyFuture* inst, 
        const SimSeriesSP& simSeries,
        InstrumentSettlementSP instSettle)
        :
        MCProductClient(
        IMultiMarketFactors::asMulti(inst->futureCurveWrapper.getSP()).get(),
        inst->today, 
        inst->yieldCurveWrapper.get(), 
        IRefLevelSP(IRefLevel::Util::makeZero(inst->today)), // fix!
        simSeries,
        IPastValuesSP(IPastValues::Util::makeTrivial(inst->today, 0.0)), // fix
        instSettle.get(),
        inst->maturityDate),
        payoffCounter(0),
        inst(inst)
        /*MCProductClient(
        IMultiMarketFactors::asMulti(inst->futureCurveWrapper.getSP()).get(),
        inst->today, 
        inst->yieldCurveWrapper.get())
        */
    {
        EnergyFuturesCurveConstSP futureCurve = EnergyFuturesCurveConstSP(inst->futureCurveWrapper.get());
        EnergyUnderlyerConstSP underlyer = futureCurve->getEnergyUnderlyer();
        needFx = !inst->fxWrapper.isEmpty();
        isT2 = !inst->parentCurveWrapper.isEmpty();

        // create maturity date place holder
        maturityDates.clear();
        maturityDates.push_back(inst->maturityDate);

        // create energy future SV generator
        expFpGen = SVGenExpectedEnergyFutureSP(
            new SVGenExpectedEnergyFuture(
            inst->measureDate,
            futureCurve,
            maturityDates,
            inst->doLog)); // compute log

        // create domestic df SV generator
        domDfGen = SVGenDiscFactorSP(
            new SVGenDiscFactor(
            inst->today, 
            YieldCurveConstSP(inst->yieldCurveWrapper.get()),
            instSettle, 
            inst->measureDate));

        if (needFx) {
            // create foreign df SV generator
            forDfGen = SVGenDiscFactorSP(
                new SVGenDiscFactor(
                inst->today, 
                YieldCurveConstSP(underlyer->getYieldCurve()),
                instSettle, 
                inst->measureDate));

            // create fx SV generator
            fxGen = SVGenSpotSP(new SVGenSpot(1, inst->measureDate));
        }

        if (isT2) {
            // create tier1 SV generator
            expFpParentGen = SVGenExpectedEnergyFutureSP(
                new SVGenExpectedEnergyFuture(
                inst->measureDate,
                EnergyFuturesCurveConstSP(inst->parentCurveWrapper.get()),
                maturityDates,
                inst->doLog)); // compute log
        }

        // write row header to dump file
        if (inst->dumpPrices) {
            dumpFile.open(inst->dumpFileName.c_str(), ios_base::app);
            dumpFile.precision(10);
            dumpFile << inst->futureCurveWrapper.get()->getName() << ":" << 
                inst->measureDate.toString() << ":" << 
                inst->maturityDate.toString() << " started: -----" << endl;
        }
    }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator* pathGen,
        IMCPrices& prices)
    {
        // only retreive the price at the measurement date along each path:
        double diffusedFuture = expFp->getExpFuturePrice(0);
        if (!inst->doSquare)
            prices.add(diffusedFuture);
        else
            prices.add(diffusedFuture*diffusedFuture);

        if (inst->dumpPrices) {
            if (!inst->doSquare)
                dumpFile << payoffCounter+1 << "\t" << diffusedFuture << endl;
            else
                dumpFile << payoffCounter+1 << "\t" << diffusedFuture*diffusedFuture << endl;
        }

        ++payoffCounter;
    }
};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* MCEnergyFuture::createProduct(const MonteCarlo* model) const 
{
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); // create empty one
    simSeries->addDates(DateTimeArray(1, maturityDate));

    // this is required
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new MC(this, simSeries, instSettle);    
}



///////////////////////////////////////////////////////////////////////////

class MCEnergyMultiFuture: public GenericNFactor,
    public virtual IMCIntoProduct{

public:
    static CClassConstSP const TYPE; 

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const; // see below

    virtual void Validate() 
    {
        const string & method = "MCEnergyMultiFuture::Validate";
        static CClassConstSP iENGType = CClass::forName("EnergyFuturesCurve");
        static CClassConstSP iFXType = CClass::forName("FXAsset");
        int nbEngs = 0;
        int nbForEngs = 0;
        int nbFxs = 0;

        for (int i = 0; i < assets->NbAssets(); ++i ) 
        {
            IMarketFactorConstSP factor = assets->getFactor(i);
            CClassConstSP type = factor->getClass();

            // make sure all factors are EnergyFuturesCurve:
            if (iENGType->isAssignableFrom(type)) 
            {
                // check maturity date:
                EnergyFuturesCurveConstSP futureCurve = EnergyFuturesCurveConstSP::dynamicCast(factor);
                EnergyUnderlyerConstSP underlyer = futureCurve->getEnergyUnderlyer();

                // construct the expiry date array of the future curve:
                //DateTimeArray futureMaturityDates;
                //for (size_t j = 0; j < futureCurve->getSize(); ++j) {
                //    EnergyContractLabel currExpiryLabel(futureCurve->getExpiryLabel(j));
                //    DateTime currExpiryDate = underlyer->expiryDate(currExpiryLabel);
                //    futureMaturityDates.push_back(currExpiryDate);
                //}
                const DateTimeArray& futureMaturityDates = futureCurve->getFutureMaturityDates();

                // make sure the target maturity is on the curve:
                maturityDates[i].find(futureMaturityDates); 

                ++nbEngs;

                // check if the energy is foreign denominated:
                if (underlyer->getYieldCurve()->getName() != discountYieldCurveName()) {
                    ++nbForEngs;
                }
            }
            else if (iFXType->isAssignableFrom(type)) 
            {
                // check FX asset?
                ++nbFxs;
            }
            else
                throw ModelException(method, "MCEnergyMultiFuture only supports "
                "energy and fx (for quanto) assets");
        }

        // more input checkings:
        if (nbFxs > nbEngs)
            throw ModelException(method, "Cannot have more FX (for quanto) "
            "assets than energys!");

        if (nbFxs != nbForEngs)
            throw ModelException(method, "The number of FX assets and foreign "
            "denominated energy assets must be equal!");
    }

    void validatePop2Object()
    {
        const string & method = "MCEnergyMultiFuture::validatePop2Object";
        if (measureDates.size() != maturityDates.size())
            throw ModelException(method, "Unequal number of measure and maturity dates!");

        for (int i = 0; i < measureDates.size(); ++i) {
            if (maturityDates[i] < measureDates[i])
                throw ModelException(method, "Please make sure each measure date is no later then the corresponding maturity date!");
        }
    }

    virtual DateTime getValueDate() const { return DateTime(); };

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz)
    {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MCEnergyMultiFuture, clazz);
        SUPERCLASS(GenericNFactor);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultMCEnergyMultiFuture);
        FIELD(today, "today");
        FIELD(maturityDates, "maturity date of the future");
        FIELD(measureDates, "measurement date");
        FIELD(doLog, "compute log");
        FIELD(doSquare, "compute price squared");
        FIELD_MAKE_OPTIONAL(doSquare);
        FIELD(doDeltaPricing, "compute delta price (change in price)");
        FIELD_MAKE_OPTIONAL(doDeltaPricing);
        FIELD(nDayLag, "delta price is computed as: Price(measureDate + nDayLag) - Price(measureDate)");
        FIELD_MAKE_OPTIONAL(nDayLag);
        FIELD(dumpPrices, "dump prices to text file");
        FIELD_MAKE_OPTIONAL(dumpPrices);
        FIELD(dumpFileName, "dump file name");
        FIELD_MAKE_OPTIONAL(dumpFileName);
    }

    static IObject* defaultMCEnergyMultiFuture() {
        return new MCEnergyMultiFuture();
    }

private:
    MCEnergyMultiFuture() : 
       GenericNFactor(TYPE), 
           doSquare(false),
           doDeltaPricing(false),
           nDayLag(0),
           dumpPrices(false),
           dumpFileName("c:/debugEnergySRM3.txt")
       {}; 

       class MC;
       friend class MC;

       DateTime today; 
       DateTimeArray maturityDates;    // maturity of the future price that we want to observe
       DateTimeArray measureDates;     // measurement date
       bool doLog;                     // compute log
       bool doSquare;                  // compute price squared
       bool doDeltaPricing;            // compute delta price (change in price) instead of price
       int nDayLag;                    // price delta is computed as:
       // Price(measureDate + nDayLag) - Price(measureDate)
       bool dumpPrices;
       string dumpFileName;

private:
    mutable DoubleArray pastValues; // dummy
    };

CClassConstSP const MCEnergyMultiFuture::TYPE = CClass::registerClassLoadMethod(
    "MCEnergyMultiFuture", typeid(MCEnergyMultiFuture), MCEnergyMultiFuture::load);

// for class loading 
bool MCEnergyMultiFutureLoad() {
    return (MCEnergyMultiFuture::TYPE != 0);
}


/* MC view of the MCEnergyMultiFuture product */
class MCEnergyMultiFuture::MC : public MCProductClient /*,
                                                       public IMCStatelessProductClient*/
{
private:
    const MCEnergyMultiFuture* inst;                // reference ptr to original inst

    // state var and state var generators
    vector<SVExpEnergyFutureSP> expFps;          // future price SV
    vector<SVGenExpectedEnergyFutureSP> expFpGens;  // future price SV generator
    vector<SVExpEnergyFutureSP> expFps2;         // for doing delta pricing, this vec can be 
    vector<SVGenExpectedEnergyFutureSP> expFpGens2; // empty.
    vector<SVDiscFactorSP> forDfs;               // a fx for each foreign denominated energy
    vector<SVGenDiscFactorSP> forDfGens;               //
    SVDiscFactorSP domDf;                        // only 1 domestic df per product
    SVGenDiscFactorSP domDfGen;                        //

    vector<MCPath::IStateVarSP> fxs; 
    //MCPath::IStateVarSP fx; // note: SVGenSpot and MCPath::IStateVar can handle multiple FXs
    SVGenSpotSP fxGen;

private:
    // the product needs to remember these:
    vector<DateTimeArray> maturityVector;
    size_t nbFxs;
    DateTimeArray fxDates;
    long payoffCounter;
    ofstream dumpFile;

protected:
    /** Override default method on IMCProduct. This method is called every time
    the path generator is changed (which is, at the moment, when the
    past path generator is created, and then when the future path
    generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen) 
    {
        expFps.clear();
        expFps2.clear();
        forDfs.clear();
        fxs.clear();

        domDf = domDfGen->getSVDiscFactor(newPathGen);
        /*if (fxGen.get()) {
        fx = fxGen->getSpotSV(newPathGen);
        }
        */
        for (size_t i = 0; i < expFpGens.size(); ++i) {
            expFps.push_back(expFpGens[i]->getSVExpEnergyFuture(newPathGen));
        }        
        for (size_t i = 0; i < expFpGens2.size(); ++i) { // for delta pricing
            expFps2.push_back(expFpGens2[i]->getSVExpEnergyFuture(newPathGen));
        }   
        for (size_t i = 0; i < forDfGens.size(); ++i) {
            forDfs.push_back(forDfGens[i]->getSVDiscFactor(newPathGen));
        }
        for (size_t i = 0; i < nbFxs; ++i) {
            fxs.push_back(fxGen->getSpotSV(newPathGen));
        }
    }

public:
    /** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const 
    {
        svCollector->append(domDfGen.get());
        if (fxGen.get()) {
            svCollector->append(fxGen.get());
        }
        for (size_t i = 0; i < expFpGens.size(); ++i) {
            svCollector->append(expFpGens[i].get());
        }
        for (size_t i = 0; i < expFpGens2.size(); ++i) { 
            svCollector->append(expFpGens2[i].get()); // for delta pricing
        }
        for (size_t i = 0; i < forDfGens.size(); ++i) {
            svCollector->append(forDfGens[i].get());
        }
        /*for (size_t i = 0; i < fxGens.size(); ++i) {
        svCollector->append(fxGens[i].get());
        }*/
    }

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC( const MCEnergyMultiFuture* inst, 
        const SimSeriesSP& simSeries)
        :
        MCProductClient(
        inst->assets.get(), //IMultiMarketFactors::asMulti(inst->futureCurveWrapper.getSP()).get(),
        inst->today, 
        inst->discount.get(), //inst->yieldCurveWrapper.get(), 
        IRefLevelSP(IRefLevel::Util::makeZero(inst->today)), // fix!
        simSeries,
        IPastValuesSP(IPastValues::Util::makeTrivial(inst->today, inst->pastValues)), // fix
        inst->instSettle.get(),
        inst->maturityDates[0]),
        payoffCounter(0),
        inst(inst)
    {
        const string & method = "MCEnergyMultiFuture::MC::MC";
        static CClassConstSP iENGType = CClass::forName("EnergyFuturesCurve");
        static CClassConstSP iFXType = CClass::forName("FXAsset");

        maturityVector.clear();
        expFpGens.clear();
        expFpGens2.clear();
        forDfGens.clear();
        //fxGens.clear();
        fxDates.clear();

        // create domestic df SV generator
        domDfGen = SVGenDiscFactorSP(
            new SVGenDiscFactor(
            inst->today, 
            YieldCurveConstSP(inst->discount.get()),
            inst->instSettle, 
            inst->measureDates[0]));

        nbFxs = 0;
        for (int i = 0; i < mfAsset->NbAssets(); ++i) 
        {
            IMarketFactorConstSP factor = mfAsset->getFactor(i);
            CClassConstSP type = factor->getClass();

            if (iENGType->isAssignableFrom(type)) 
            {
                EnergyFuturesCurveConstSP futureCurve = EnergyFuturesCurveConstSP::dynamicCast(factor);
                EnergyUnderlyerConstSP underlyer = futureCurve->getEnergyUnderlyer();
                maturityVector.push_back(DateTimeArray(1, inst->maturityDates[i]));

                // create energy SV gens:
                expFpGens.push_back(SVGenExpectedEnergyFutureSP(
                    new SVGenExpectedEnergyFuture(
                    inst->measureDates[i],
                    futureCurve,
                    maturityVector.back(),
                    inst->doLog))); // compute log

                // need to create another enrg SV gen if doing delta pricing:
                if (inst->doDeltaPricing) {
                    expFpGens2.push_back(SVGenExpectedEnergyFutureSP(
                        new SVGenExpectedEnergyFuture(
                        inst->measureDates[i].rollDate(inst->nDayLag),
                        futureCurve,
                        maturityVector.back(),
                        inst->doLog))); // compute log
                }

                // create foreign df sv gens, if necessary:
                if (underlyer->getYieldCurve()->getName() != inst->discountYieldCurveName()) {
                    forDfGens.push_back(SVGenDiscFactorSP(
                        new SVGenDiscFactor(
                        inst->today, 
                        YieldCurveConstSP(underlyer->getYieldCurve()),
                        inst->instSettle, 
                        inst->measureDates[i])));
                }
            }
            else if (iFXType->isAssignableFrom(type))
            {
                ++nbFxs;
                fxDates.push_back(inst->measureDates[i]);
            }
            else
                throw ModelException(method, "MCEnergyMultiFuture only supports "
                "energy and fx (for quanto) assets");
        }

        // create fx sv gens:
        if (nbFxs) {
            DateTime::doSortUniq(fxDates);
            //fxGen = SVGenSpotSP(new SVGenSpot(simSeries));
            //fxGen = SVGenSpotSP(new SVGenSpot(nbFxs, fxDates));
            fxGen = SVGenSpotSP(new SVGenSpot(mfAsset->NbAssets(), fxDates));
        }

        // write row header to dump file
        if (inst->dumpPrices) {
            dumpFile.open(inst->dumpFileName.c_str(), ios_base::app);
            dumpFile.precision(10);
            dumpFile << "MultiEnergy:";
            for (size_t i = 0; i < expFpGens.size(); ++i)
                dumpFile << "[" << expFpGens[i]->getName() <<
                ":" << expFpGens[i]->getCalcDate().toString() << 
                ":" << expFpGens[i]->getDates()[0].toString() << "]";
            dumpFile << " started ------------------" << endl;
        }
    }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator* pathGen,
        IMCPrices& prices)
    {
        // only retreive the price at the measurement date along each path:
        double sum = 0.0;
        for (size_t i = 0; i < expFps.size(); ++i) {
            sum += expFps[i]->getExpFuturePrice(0);
        }

        if (!inst->doSquare)
            prices.add(sum);
        else
            prices.add(sum*sum);

        if (inst->dumpPrices) {
            dumpFile << payoffCounter + 1 << "\t";
            if (!inst->doSquare) {
                if (!inst->doDeltaPricing) { // dump prices
                    for (size_t i = 0; i < expFps.size(); ++i)
                        dumpFile << expFps[i]->getExpFuturePrice(0) << "\t";
                }
                else { // dump delta prices
                    for (size_t i = 0; i < expFps.size(); ++i)
                        dumpFile << expFps2[i]->getExpFuturePrice(0) - expFps[i]->getExpFuturePrice(0) << "\t"; 
                }
            }
            else {
                if (!inst->doDeltaPricing) {
                    for (size_t i = 0; i < expFps.size(); ++i)
                        dumpFile << expFps[i]->getExpFuturePrice(0)*expFps[i]->getExpFuturePrice(0) << "\t";
                }
                else {
                    for (size_t i = 0; i < expFps.size(); ++i)
                        dumpFile << (expFps2[i]->getExpFuturePrice(0) - expFps[i]->getExpFuturePrice(0)) * 
                        (expFps2[i]->getExpFuturePrice(0) - expFps[i]->getExpFuturePrice(0)) << "\t";
                }
            }
            dumpFile << endl;
        }

        ++payoffCounter;
    }
};


/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* MCEnergyMultiFuture::createProduct(const MonteCarlo* model) const 
{
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    int nbAssets = assets->NbAssets();
    pastValues = DoubleArray(nbAssets, 0.0);
    SimSeriesSP simSeries(new SimSeries(nbAssets)); // create empty one
    for (int i = 0; i < nbAssets; ++i)
        simSeries->addDates(DateTimeArray(1, maturityDates[i]));

    return new MC(this, simSeries);    
}


DRLIB_END_NAMESPACE
