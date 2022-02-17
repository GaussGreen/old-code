//----------------------------------------------------------------------------
//
//   Description : demo product for monte carlo 
//
//   Author      : Jay Z Wang
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/SVGenDemo.hpp"
#include "edginc/MonteCarlo.hpp"

DRLIB_BEGIN_NAMESPACE

class MCDemoProduct: public GenericSimpleIR,
                     public virtual IMCIntoProduct{

public:
    static CClassConstSP const TYPE; 

    /** Implementation of MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const; // see below

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) 
    {
        // delegate to base class 
        GenericSimpleIR::GetMarket(model, market);

        // if we were to put tie an asset with this product, we would do something like
        // asset.getData(model, market);
    };

    virtual void Validate() {};

    virtual DateTime getValueDate() const { return DateTime(); };
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MCDemoProduct, clazz);
        SUPERCLASS(GenericSimpleIR);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultMCDemoProduct);
        FIELD(expDate, "exp date");
    }

    static IObject* defaultMCDemoProduct(){
        return new MCDemoProduct();
    }
    
private:
    MCDemoProduct():GenericSimpleIR(TYPE) {}; 

private:
    class MC;
    friend class MC;

    // CAssetWrapper asset;

    DateTime expDate;
};

CClassConstSP const MCDemoProduct::TYPE = CClass::registerClassLoadMethod(
    "MCDemoProduct", typeid(MCDemoProduct), MCDemoProduct::load);



/* MC product class MCDemoProduct */
class MCDemoProduct::MC : public MCProductClient,
                            public IMCStatelessProductClient
{
private:
    // state var and state var generators
    SVGenDemo::IStateVarSP jsv;
    SVGenDemo::INewStateVarSP jnsv;
    SVGenDemoSP jsvGen;
    DateTime expDate;
    int _numSlices;


protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        // jsv = jsvGen->getDemoStateVar(newPathGen);
        jnsv = jsvGen->getNewDemoStateVar(newPathGen);
    };
public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(jsvGen.get());
    }
    
    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC(const MCDemoProduct*          inst,
       const SimSeriesSP&       simSeries,
       InstrumentSettlementSP   instSettle):
        MCProductClient(IMultiMarketFactors::asMulti(inst->coupon).get(),
                        inst->valueDate,
                        inst->discount.get(),
                        IRefLevelSP(IRefLevel::Util::makeZero(inst->valueDate)),
                        simSeries,
                        IPastValuesSP(IPastValues::Util::makeTrivial(inst->valueDate, 0.0)),
                        instSettle.get(),
                        inst->expDate),
        jsvGen(new SVGenDemo()),
        expDate(inst->expDate)
        {}

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*  pathGen,
                        IMCPrices&                prices) 
    {
        // call method on jsv to get current spot
        double* spot = jsv->getSpot();
        double sum(0);
        for (int i = 0; i < _numSlices; ++i)
            sum += spot[i];
        sum /= _numSlices;
        std::cout << "payoff = " << sum << "\n";
        prices.add(sum);
    }

    class HistoricalContext : public IHistoricalContext {
        friend class MC;
        double sumSoFar;
    public:
        HistoricalContext():sumSoFar(0) {};
        virtual void deepCopyTo( IHistoricalContext* destination ) const
        {
            HistoricalContext* dest = static_cast<HistoricalContext*>( destination );
            *dest = *this;
        }
    };

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
        int numDates = simDates.size();
        _numSlices = numDates;
        vector<int> keyDate;
        keyDate.reserve(numDates);
        for (int i = 0; i < numDates; ++i)
            keyDate.push_back(i);
        return keyDate;
    }

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices )
    {
        HistoricalContext* hContext = 
            static_cast<HistoricalContext*>(history.get()); 
        double spot = jnsv->getSpot();
        hContext->sumSoFar += spot;
        if ( currentDateIdx == _numSlices ) {  // FIXME
            hContext->sumSoFar /= _numSlices;
            // std::cout << "payoff = " << hContext->sumSoFar << "\n";
            prices.add(hContext->sumSoFar);
            hContext->sumSoFar = 0;
        }
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* MCDemoProduct::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(DateTimeArray(1, expDate));
    // this is required
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new MC(this, simSeries, instSettle);
}

// for class loading 
bool MCDemoProductLoad() {
    return (MCDemoProduct::TYPE != 0);
}

DRLIB_END_NAMESPACE
