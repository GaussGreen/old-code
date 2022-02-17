//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskyBondSeries.cpp
//
//   Description : Series of risky bonds. Primarily a test vehicle for SRM3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp" 
#include "edginc/SVGenDiscFactor.hpp" 
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp" 
#include "edginc/SVGenSurvivalDiscFactor.hpp" 
#include "edginc/SVGenDateOfDefault.hpp"
#include "edginc/SVGenPathWeight.hpp"
#include "edginc/CashSettlePeriod.hpp"


DRLIB_BEGIN_NAMESPACE

/** RiskyBondSeries product 
    Series of risky bonds -- pay one if survival and R if default
    valuation at "bondStart" for bond from "bondStart" to "bondMat" 
    valuation conditional on survival until bondStart */
class RiskyBondSeries: public CInstrument,
                       public virtual IMCIntoProduct {
    
public:
    static CClassConstSP const TYPE;

    virtual ~RiskyBondSeries(){}

    void validatePop2Object(){
        
        static const string method("RiskyBondSeries::validatePop2Object");
        
        if (bondStart.size() != bondMat.size()) {
            throw ModelException(method, 
                "Number of bondStart and bondMat dates must be equal");
        }
        for (int i = 0; i < bondStart.size(); i++) {
            if (bondStart[i] > bondMat[i]) {
                throw ModelException(method,
                    "Bond maturity for bond " + Format::toString(i) + " has to be after bond start");
            }
        }
    }

    void GetMarket(const IModel*        model, 
                   const CMarketDataSP  market) { 
        static const string method = "RiskyBondSeries::GetMarket";
        try {
            market->GetReferenceDate(valueDate); // populates valueDate, since optional input

            // risky bond series without fx 
            ICDSParSpreads::getMarketData(model, market.get(), creditYieldCurve.getName(), cdsParSpreads);
            creditYieldCurve.getData(model, market);

            // bad fudge, should be done by engine
            // risky bond series with fx
            if (isFx) {
                discountYieldCurve.getData(model, market);
                CAsset::getAssetMarketData(model, market.get(), "V", discountYieldCurve, fxAsset);
            } 
        } catch (exception &e) {
            throw ModelException(e, method);
        }
    }
    
    virtual void Validate() { // after GetMarket
        static const string method = "RiskyBondSeries::Validate";
        try {
            for (int i = 0; i < bondStart.size(); i++) {
                if (valueDate > bondStart[i]) {
                    throw ModelException(method,
                        "Start of bond " + Format::toString(i) + " is before value date");
                }
            }
            // TODO
        } catch (exception &e) {
            throw ModelException(e, method);
        }
    }
    
    DateTime getValueDate() const {
        return valueDate;
    }


    /** Returns the name of the instrument's discount currency. */
    string discountYieldCurveName() const {
        if (isFx) {
            return discountYieldCurve.getName();
        }
        else {
            return creditYieldCurve.getName();
        }
    }


    /** Implementation of MonteCarlo::IntoProduct interface - the implementation of this is below */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const; 
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(RiskyBondSeries, clazz);
        SUPERCLASS(CInstrument);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultRiskyBondSeries);
        FIELD(valueDate,"valuation date");
        FIELD_MAKE_OPTIONAL(valueDate);
        FIELD(cdsParSpreads, "Credit underlying");
        FIELD(creditYieldCurve,"Credit yield curve");
        FIELD(discountYieldCurve,"Discount yield curve");
        FIELD_MAKE_OPTIONAL(discountYieldCurve);
        FIELD(fxAsset, "fx");
        FIELD_MAKE_OPTIONAL(fxAsset);
        FIELD(useSmoothDefaults, "True means using SurvivalProb, False - simulated date of default");
        FIELD_MAKE_OPTIONAL(useSmoothDefaults);
        FIELD(bondStart, "When each risky bond starts");
        FIELD(bondMat, "When each risky bond ends");
        FIELD(isFx, "Is risky bond series fx or not");
    }

protected:
    RiskyBondSeries(CClassConstSP clazz):
         CInstrument(clazz), useSmoothDefaults(true) {} // extra constructor
        
private:

    friend class RiskyBondSeriesMC;

    /** helpers for MC **/
    
    //// Returns a SVGenDiscFactor for risky bond given by index idx
    SVGenDiscFactor* createDiscountMCDF(int idx) const{
        if (isFx) {
            return new SVGenDiscFactor(valueDate, discountYieldCurve.getSP(), bondStart[idx]); 
        } else {
            return new SVGenDiscFactor(valueDate, creditYieldCurve.getSP(), bondStart[idx]); 
        }
    }
        
    //// Returns a SVGenExpectedDiscFactor for risky bond given by index idx
    SVGenExpectedDiscFactor* createCreditMCDF(int idx) const{
        return new SVGenExpectedDiscFactor(bondStart[idx], bondStart[idx],
                                        creditYieldCurve.getSP(), 
                                        DateTimeArray(1, bondMat[idx]),
                                        false); // computeLog = ?
    }

    // Returns a SVGenSurvivalDiscFactor for risky bond given by index idx
    SVGenSurvivalDiscFactor* createMCSDF(int idx) const{
        return new SVGenSurvivalDiscFactor(valueDate, cdsParSpreads.getSP(), bondStart[idx]); 
    }
    SVGenDateOfDefault* createCreditDoD() const{
        return new SVGenDateOfDefault(valueDate, cdsParSpreads.getSP(), bondStart.back());
    }
    SVGenExpectedSurvivalDiscFactor* createMCExpSDF(int idx) const{
        return new SVGenExpectedSurvivalDiscFactor(bondStart[idx], bondStart[idx],
                                                cdsParSpreads.getSP(), 
                                                DateTimeArray(1, bondMat[idx]),
                                                false); // computeLog = ?
    }

    SVGenPathWeight* createMCPathW(int idx) const{
        return new SVGenPathWeight(DateTimeArray(1, bondStart[idx]));
    }

    // spot FX
    SVGenSpot* createSVGenSpot(int idx) const{
        return new SVGenSpot(1 /* number of assets */, bondStart[idx]);
    }    

    RiskyBondSeries(): CInstrument(TYPE), useSmoothDefaults(true) {}
    RiskyBondSeries(const RiskyBondSeries& rhs);            // not implemented
    RiskyBondSeries& operator=(const RiskyBondSeries& rhs); // not implemented

    static IObject* defaultRiskyBondSeries(){
        return new RiskyBondSeries();
    }
 
    /* fields */
    DateTime                valueDate;
    
    /** in the ideal case, we only have one CDS curve here and one discount curve
        if these two have different currencies, then the engine should ask for an
        additional discount curve and an fxAsset */
    ICDSParSpreadsWrapper   cdsParSpreads;
    YieldCurveWrapper       creditYieldCurve;   // discount curve if non fx
    YieldCurveWrapper       discountYieldCurve; // discount curve if fx         
    CAssetWrapper           fxAsset;
    
    DateTimeArray           bondStart;      // when each risky bond starts
    DateTimeArray           bondMat;        // when each risky bond ends

    bool                    isFx;           // not registered
    bool                    useSmoothDefaults; // false means using DateOfDefault, true - SurvivalProbabilities
};

/* MC product class for super rainbow */
class RiskyBondSeriesMC : public MCProductClient,
                          public IMCStatelessProductClient {
    // a set of state variables for each risky bond
    struct SV{
        SVDiscFactorSP           discFac;          //!< discount factor SV
        SVExpectedDiscFactorSP   expDiscFac;       //!< expected discount factor SV
        SVSurvivalDiscFactorSP   survDiscFac;      //!< survival discount factor SV
        SVExpSurvDiscFactorSP    expSurvDiscFac;   //!< expected survival discount factor SV
        SVGenSpot::IStateVarSP   spotFX;           //!< spot FX SV
        SVPathWeightSP           pathWeight;       //!< path weight SV
        DateTime                 bondStartDate;    // to know what stream this SV refers to
        
    };
    // a set of state variables generators for each set of state variables
    struct GenSV{
        SVGenDiscFactorSP                   discFac;         //!< generator for discFac
        SVGenExpectedDiscFactorSP           expDiscFac;      //!< generator for expDiscFac
        SVGenSurvivalDiscFactorSP           survDiscFac;     //!< generator for survDiscFac
        SVGenExpectedSurvivalDiscFactorSP   expSurvDiscFac;  //!< generator for expSurvDiscFac
        SVGenSpotSP                         spotFX;          //!< generator for spot FX
        SVGenPathWeightSP                   pathWeight;      //!< generator for path weight
    };
    const RiskyBondSeries*  inst;      //!< reference to original inst
    vector< vector<SV> >    sv;        //!< date slot -> vector of state variables (one per cashflow)
    vector<GenSV>           genSV;     //!< state variable generators
    map<DateTime,int> dateToIndex;     // Correspondence betweens key dates and index in sv array
    int maturityIdx;
    SVDateOfDefaultSP        dateOfDefault;    //!< date of default SV (incompatible with survival DF)
    SVGenDateOfDefaultSP     dateOfDefaultGen; //!< generator for dateOfDefault SV

private:
    /** Update our motley collection of state variables */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine("RiskyBondSeriesMC::pathGenUpdated");
        try {
            sv.reserve( dateToIndex.size() );
            for (unsigned int i = 0; i < genSV.size(); i++){
                SV tmpSV;
                tmpSV.discFac = genSV[i].discFac->getSVDiscFactor(newPathGen);
                tmpSV.expDiscFac = genSV[i].expDiscFac->getSVExpectedDiscFactor(newPathGen);                
                tmpSV.expSurvDiscFac = genSV[i].expSurvDiscFac->getSVExpectedSurvivalDiscFactor(newPathGen);
                if (inst->useSmoothDefaults)
                    tmpSV.survDiscFac = genSV[i].survDiscFac->getSVSurvivalDiscFactor(newPathGen);
                if (inst->isFx) {
                    tmpSV.spotFX = genSV[i].spotFX->getSpotSV(newPathGen);
                }
                tmpSV.bondStartDate = inst->bondStart[i];
                int index = dateToIndex[inst->bondStart[i]];

                tmpSV.pathWeight = genSV[i].pathWeight->getSVPathWeight(newPathGen);
                sv[index].push_back( tmpSV );
              
            }
            if ( ! inst->useSmoothDefaults)
                dateOfDefault = dateOfDefaultGen->getSVDateOfDefault(newPathGen);
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    /** Ask for our unusual collection of state variables */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        for (unsigned int i = 0; i < genSV.size(); i++){
            svCollector->append(genSV[i].discFac.get());
            svCollector->append(genSV[i].expDiscFac.get());            
            if (inst->useSmoothDefaults)
                svCollector->append(genSV[i].survDiscFac.get());
            svCollector->append(genSV[i].expSurvDiscFac.get());
            svCollector->append(genSV[i].pathWeight.get());

            if (inst->isFx) {
                svCollector->append(genSV[i].spotFX.get());
            }
        }
        if ( ! inst->useSmoothDefaults)
            svCollector->append(dateOfDefaultGen.get());

    }
    
    /** Need to call parent's constructor - a bit of a mess */
    /** this has to be addressed again ...!!! */
    RiskyBondSeriesMC(const RiskyBondSeries*    inst,
                      SimSeriesSP               simSeries,
                      InstrumentSettlementSP    instSettle,
                      bool                      dummyBoolean):
        MCProductClient(inst->fxAsset.get(), // irrelevant
                        inst->valueDate,
                        inst->discountYieldCurve.get(),
                        IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                        simSeries, // irrelevant - fix
                        IPastValuesSP(IPastValues::Util::
                                      makeTrivial(inst->valueDate, 0.0)), // fix
                        instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst), sv(inst->bondStart.size()), genSV(inst->bondStart.size()){
        for (unsigned int i = 0; i < genSV.size(); i++){
            genSV[i].discFac = SVGenDiscFactorSP(inst->createDiscountMCDF(i));
            genSV[i].expDiscFac = SVGenExpectedDiscFactorSP(inst->createCreditMCDF(i));
            if (inst->useSmoothDefaults)
                genSV[i].survDiscFac = SVGenSurvivalDiscFactorSP(inst->createMCSDF(i));
            genSV[i].expSurvDiscFac = SVGenExpectedSurvivalDiscFactorSP(inst->createMCExpSDF(i));
            if (inst->isFx) {
                genSV[i].spotFX = SVGenSpotSP(inst->createSVGenSpot(i));
            }
            genSV[i].pathWeight.reset(inst->createMCPathW(i));
        }
        if (!inst->useSmoothDefaults)
            dateOfDefaultGen = SVGenDateOfDefaultSP(inst->createCreditDoD());
        // Set up mapping between date and location in sv vector.
        for (int j = 0; j < inst->bondStart.size(); ++j ) {
            const DateTime& date = inst->bondStart[j];
            map<DateTime, int>::iterator I = dateToIndex.find( date );
            if ( I == dateToIndex.end() ) {
        int currentSize = dateToIndex.size();
                dateToIndex[date] = currentSize;
            }
        }
    }
    /** this has to be addressed again ...!!! */
    RiskyBondSeriesMC(const RiskyBondSeries*    inst,
                      SimSeriesSP               simSeries,
                      InstrumentSettlementSP    instSettle):
        MCProductClient(IMultiMarketFactors::asMulti(inst->creditYieldCurve.getSP()).get(),  // irrelevant
                        inst->valueDate,
                        inst->creditYieldCurve.get(),
                        IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)), // fix!
                        simSeries, // irrelevant - fix
                        IPastValuesSP(IPastValues::Util::
                                      makeTrivial(inst->valueDate, 0.0)), // fix
                        instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst), sv(inst->bondStart.size()), genSV(inst->bondStart.size()){
        for (unsigned int i = 0; i < genSV.size(); i++){
            genSV[i].discFac = SVGenDiscFactorSP(inst->createDiscountMCDF(i));
            genSV[i].expDiscFac = SVGenExpectedDiscFactorSP(inst->createCreditMCDF(i));
            if (inst->useSmoothDefaults)
                genSV[i].survDiscFac = SVGenSurvivalDiscFactorSP(inst->createMCSDF(i));
            genSV[i].expSurvDiscFac = SVGenExpectedSurvivalDiscFactorSP(inst->createMCExpSDF(i));
            genSV[i].pathWeight.reset(inst->createMCPathW(i));
        }
        if (!inst->useSmoothDefaults)
            dateOfDefaultGen = SVGenDateOfDefaultSP(inst->createCreditDoD());

        // Set up mapping between date and location in sv vector.
        for (int j = 0; j < inst->bondStart.size(); ++j ) {
            const DateTime& date = inst->bondStart[j];
            map<DateTime, int>::iterator I = dateToIndex.find( date );
            if ( I == dateToIndex.end() ) {
                int currentSize = dateToIndex.size();
                dateToIndex[date] = currentSize;
            }
        }
    }

    /** Called within the simulation loop */
    void payoff(const MCPathGenerator*   pathGen,
                IMCPrices&                 prices) {
        double price = 0.0;
        double rec = inst->cdsParSpreads->getRecovery();
        for (unsigned int t = 0; t < sv.size(); t++) { // one slot for each date 
            price += updatePrice( t, rec );
        }
        prices.add(price);
    }

    // Stateless pricer version

    class HistoricalContext : public IHistoricalContext {
        friend class RiskyBondSeriesMC;
        double priceSoFar;
    public:
        HistoricalContext()
        {
            priceSoFar = 0;
        }

        virtual void deepCopyTo( IHistoricalContext* destination ) const
        {
            HistoricalContext* dest = static_cast<HistoricalContext*>( destination );
            *dest = *this;
        }

    };

    // IMCStatelessProductClient
    virtual IHistoricalContextSP createHistoricalContext()    
    {
        return IHistoricalContextSP(new HistoricalContext());
    }

    virtual IHistoricalContextSP getInitialHistoricalContext()    
    {
        return IHistoricalContextSP(new HistoricalContext());
    }


    virtual DateTimeArray getPastDates()
    {
        return DateTimeArray();
    }

    virtual vector<int> finalize( const DateTimeArray& simDates )
    {
        const DateTimeArray& D = simDates;
        vector<int> keyDate;
        for ( map<DateTime,int>::const_iterator I = dateToIndex.begin();
              I != dateToIndex.end(); ++I ) {
            array<DateTime>::const_iterator J;
            J = std::find( D.begin(), D.end(), I->first );
            if ( J == D.end() ) {
                throw ModelException("Could not find date '" 
                    + I->first.toString() + "' in the timeline"); 
            }
            keyDate.push_back( J - D.begin() );
        }
        maturityIdx = dateToIndex.size() - 1;
        return keyDate;
    }

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices ) 
    {
        HistoricalContext* hContext = 
            static_cast<HistoricalContext*>(history.get()); 

        double rec = inst->cdsParSpreads->getRecovery();
        hContext->priceSoFar += updatePrice( currentDateIdx, rec );

        // register price at maturity
        if ( currentDateIdx == maturityIdx )
            prices.add( hContext->priceSoFar );
    }

    double updatePrice( int t, double rec )
    {
        double price = 0;
        for ( unsigned int i = 0; i < sv[t].size(); ++i ) {// one SV for each risky cashflow
            SV& sv = this->sv[t][i];
            double discFac = sv.discFac->firstDF() * sv.pathWeight->getWeight(0);
            double expDiscFac = sv.expDiscFac->firstDF();
            double survDiscFac;
            if (inst->useSmoothDefaults)
                survDiscFac = sv.survDiscFac->firstSDF();
            else
                survDiscFac = (dateOfDefault->getDateOfDefault() > sv.bondStartDate) ? 1.0 : 0.0 ;

            double expSurvDiscFac = sv.expSurvDiscFac->firstExpSDF();              
            const SVPath& path = sv.expSurvDiscFac->path();
            if (inst->isFx) {
                const SVPath& pathSpotFX = sv.spotFX->path(0); // 1 asset
                for (int iStep = path.begin(); iStep < path.end(); iStep++){
                    price += discFac *                                               // discount using payoff ccy
                        pathSpotFX[iStep] *                                          // convert into payoff ccy
                        (expDiscFac * survDiscFac * expSurvDiscFac +                 // pay one if survival
                        rec * expDiscFac * (1.0 - survDiscFac * expSurvDiscFac));    // pay rec if default
                }
            } else {
                for (int iStep = path.begin(); iStep < path.end(); iStep++){
                    price += discFac *                                                // discount
                        (expDiscFac * survDiscFac * expSurvDiscFac +                  // pay one if survival
                        rec * expDiscFac * (1.0 - survDiscFac * expSurvDiscFac));     // pay rec if default
                }
            }
        }
        return price;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RiskyBondSeries::createProduct(const MonteCarlo* model) const {
     // the sim series here is irrelevant - to do: fix MCProductClient
    SimSeriesSP simSeries(new SimSeries(1));
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    simSeries->addDates(DateTimeArray(1, bondMat.back())); // so getLastDate() works
    if (this->isFx) {
        return new RiskyBondSeriesMC(this, simSeries, instSettle, true);
    } else {
        return new RiskyBondSeriesMC(this, simSeries, instSettle);
    }
}

CClassConstSP const RiskyBondSeries::TYPE = CClass::registerClassLoadMethod(
    "RiskyBondSeries", typeid(RiskyBondSeries), load);

// * for class loading (avoid having header file) */
bool RiskyBondSeriesLoad() {
    return (RiskyBondSeries::TYPE != 0);
}

DRLIB_END_NAMESPACE
