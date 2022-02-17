//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ForwardOptionSeries.cpp
//
//   Description : IMCPrices a series of options on [asset] forwards.
//                 Primarily a test vehicle for SRM3
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedSpot.hpp"

DRLIB_BEGIN_NAMESPACE

/** ForwardOptionSeries product  - a strip of options on forwards.
    Derived from Generic1Factor just in case
    we ever want to book one - although several fields are a bit meaningless */
class ForwardOptionSeries: public Generic1Factor, 
                           virtual public IMCIntoProduct{
    /// fields ////////
    DateTimeArray           fwdStart;      // when each forward starts
    DateTimeArray           fwdMat;        // when each forward ends
    DoubleArray             strike;        // strike for each forward
    BoolArray               isCall;        // is option a call
public:
    static CClassConstSP const TYPE;
    friend class ForwardOptionSeriesMC;

    // validation
    void validatePop2Object(){
        static const string method("ForwardOptionSeries::validatePop2Object");
        if (fwdStarting){
            throw ModelException(method, "Forward starting not applicable to "
                                 "options on a forward");
        }
        if (!oneContract){
            throw ModelException(method, "Notional based option not supported");
        }
        if (premiumSettle.get()){
            throw ModelException(method, "Premium settlement not supported");
        }
        if (fwdStart.empty()){
            throw ModelException(method, "No fwdStart dates supplied");
        }
        if (fwdStart.size() != fwdMat.size() ||
            fwdStart.size() != strike.size() ||
            fwdStart.size() != isCall.size()){
            throw ModelException(method, "Number of fwdStart and fwdMat dates"
                                 " must both be eqal to the number of strikes "
                                 "and isCall flags");
        }            
    }

    virtual void Validate(){}

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    /// helpers for MC ///
    //// Returns a SVGenDiscFactor for option given by index idx
    SVGenDiscFactor* createMCDF(int idx) const{
        return new SVGenDiscFactor(valueDate, discount.getSP(),
                                instSettle, fwdStart[idx]);
    }
    //// Returns a SVGenExpectedDiscFactor for option given by index idx
    SVGenExpectedDiscFactor* createMCExpDF(int idx) const{
        return new SVGenExpectedDiscFactor(fwdStart[idx], fwdStart[idx],
                                        discount.getSP(),
                                        DateTimeArray(1, fwdMat[idx]),
                                        false);
    }
    //// Returns a SVGenExpectedSpot for option given by index idx
    SVGenExpectedSpot* createMCExpSpot(int idx) const{
        return new SVGenExpectedSpot(0, // asset index
                                  fwdStart[idx],
                                  DateTimeArray(1, fwdMat[idx]));
    }

    //// Returns a SVGenSpot for option given by index idx
    SVGenSpot* createSVGenSpot(int idx) const{
        return new SVGenSpot(1 /* number of assets */, fwdStart[idx]);
    }

    ForwardOptionSeries(): Generic1Factor(TYPE) {} // for reflection
    ForwardOptionSeries(const ForwardOptionSeries& rhs); // not implemented
    ForwardOptionSeries& operator=(
        const ForwardOptionSeries& rhs); // not implemented

    static IObject* defaultForwardOptionSeries(){
        return new ForwardOptionSeries();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ForwardOptionSeries, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultForwardOptionSeries);
        FIELD(fwdStart,         "When each forward starts");
        FIELD(fwdMat,           "When each forward ends");
        FIELD(strike,           "Strike for each forward");
        FIELD(isCall,           "Is option a call");
    }
};

/* MC product class for super rainbow */
class ForwardOptionSeriesMC : public MCProductClient,
                              public IMCStatelessProductClient {
    // a set of state variables for each option
    struct SV{
        SVPathSP         expSpot; //!< expected spot state variable
        SVExpectedDiscFactorSP expDf;   //!< expected df state variable
        SVDiscFactorSP         df;      //!< df state variable
        SVGenSpot::IStateVarSP               spot;    //!< Spot state variable. Can be 0
        int                               i;       //!< Index into inst arrays for ref data
    };
    // a set of state variables generators for each set of state variables
    struct GenSV{
        SVGenExpectedSpotSP          expSpot; //!< Generator for exp spot
        SVGenSpotSP                  spot;    //!< Generator for spot - can be 0
        SVGenExpectedDiscFactorSP    expDf;   //!< Generator for expected df
        SVGenDiscFactorSP            df;      //!< Generator for df
    };
    const ForwardOptionSeries*  inst;      //!< reference to original inst
    vector<vector<SV> >         sv;        //!< state variables
    vector<GenSV>               genSV;     //!< state variable generators
    map<DateTime,int> dateToIndex;     // Correspondence betweens key dates and index in sv array
    int maturityIdx;

private:
    /** Update our motley collection of state variables */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine("ForwardOptionSeriesMC::pathGenUpdated");
        try {
	    sv.clear();
            sv.reserve( dateToIndex.size() );
	    vector<SV> dummy;
            size_t i;
	    for (i = 0; i < dateToIndex.size(); ++i)
	        sv.push_back(dummy);
            for (i = 0; i < genSV.size(); i++){
                SV tmpSV;
                if (genSV[i].spot.get()){
                    tmpSV.spot = genSV[i].spot->getSpotSV(newPathGen);
                } else {
                    tmpSV.expSpot = genSV[i].expSpot->getExpSpotSV(newPathGen);
                    tmpSV.expDf = genSV[i].expDf->getSVExpectedDiscFactor(newPathGen);
                }
                tmpSV.df = genSV[i].df->getSVDiscFactor(newPathGen);
                tmpSV.i = i;
                int index = dateToIndex[inst->fwdStart[i]];
                sv[index].push_back( tmpSV );
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    };

public:
    /** Ask for our unusual collection of state variables */
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        for (unsigned int i = 0; i < genSV.size(); i++){
            if (genSV[i].spot.get()){
                svCollector->append(genSV[i].spot.get());
            } else {
                svCollector->append(genSV[i].expSpot.get());
                svCollector->append(genSV[i].expDf.get());
            }
            svCollector->append(genSV[i].df.get());
        }
    }
    
    /** Need to call parent's constructor - a bit of a mess */
    ForwardOptionSeriesMC(const ForwardOptionSeries* inst,
                          SimSeriesSP                simSeries):
        MCProductClient(inst->asset.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        IRefLevelSP(IRefLevel::Util::
                                    makeZero(inst->valueDate)), // fix!
                        simSeries, // irrelevant - fix
                        IPastValuesSP(IPastValues::Util::
                                      makeTrivial(inst->valueDate, 0.0)), // fix
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst), sv(inst->isCall.size() ), genSV(inst->isCall.size()){
        for (unsigned int i = 0; i < genSV.size(); i++){
            if (inst->fwdStart[i].equals(inst->fwdMat[i])){
                // zero length forward - can do using spot rather than
                // expected spot
                genSV[i].spot = SVGenSpotSP(inst->createSVGenSpot(i));
            } else {
                genSV[i].expSpot = SVGenExpectedSpotSP(inst->createMCExpSpot(i));
                genSV[i].expDf = SVGenExpectedDiscFactorSP(inst->createMCExpDF(i));
            }
            genSV[i].df = SVGenDiscFactorSP(inst->createMCDF(i));
        }
        // Set up mapping between date and location in sv vector.
        for (int j = 0; j < inst->fwdStart.size(); ++j ) {
            const DateTime& date = inst->fwdStart[j];
            map<DateTime, int>::iterator I = dateToIndex.find( date );
            if ( I == dateToIndex.end() ) {
		int currentSize = dateToIndex.size();
                dateToIndex[date] = currentSize;
            }
        }
    }

    /** Called within the simulation loop */
    void payoff(const IPathGenerator*  pathGen,
                IMCPrices&                prices) 
    {
        double price = 0.0;
        for (unsigned int i = 0; i < sv.size(); i++) {
            price += updatePrice( i );
#ifdef DEBUGSRM            
            fprintf(stderr, "i= %d price= %g\n", i, price);
#endif //FIXME            
        }
        prices.add(price);
    }

    class HistoricalContext : public IHistoricalContext {
        friend class ForwardOptionSeriesMC;
        double priceSoFar;
    public:
        HistoricalContext() : priceSoFar(0) {}

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
        maturityIdx = dateToIndex.size() - 1;
        int lastPastDateIdx;
        DateTimeArray date( dateToIndex.size() );
        int k = 0;
        for ( map<DateTime,int>::iterator I = dateToIndex.begin();
            I != dateToIndex.end(); ++I ) 
        {
            date[k] = I->first;
            k += 1;
        }
        return DateTime::createMapping2( 
            simDates, 
            date, 
            simDates[0], 
            lastPastDateIdx);
    }

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices ) 
    {
        HistoricalContext* hContext = 
            static_cast<HistoricalContext*>(history.get()); 
        
        hContext->priceSoFar += updatePrice( currentDateIdx );
 
        // register price at maturity
        if ( currentDateIdx == maturityIdx )
            prices.add( hContext->priceSoFar );
    }

private:
    double updatePrice( int t )
    {
        double price = 0.0;
        for ( unsigned int j = 0; j < sv[t].size(); ++j ) {// one SV for each interval of series
            SV& sv = this->sv[t][j];
            int i = sv.i; // Index into instruments arrays for reference data
            double df = sv.df->firstDF();
            if (!sv.expSpot){
                // price is just a vanilla:
                // ie DF * (MAX(FX - K, 0)    for a call
                const SVPath& path = sv.spot->path(0); // 1 asset
                // at most only one date but possibly zero dates though
                double option = 0.0;
                for (int iStep = path.begin(); iStep < path.end(); ++iStep){
                    option = path[iStep]; // At fwdStart, spot FX at fwdMat
                }
                option = inst->isCall[i]? 
                    option - inst->strike[i]: (inst->strike[i] - option);
                // terse
                price += option > 0.0? option * df: 0.0;
#ifdef DEBUGSRM                
                fprintf(stderr, "updatePrice: spot t= %d j= %d option= %g df= %g price= %g\n", t,j,option, df, price);
#endif //FIXME                

            } else {
                // price is discounted value of a forward contract at fwdStart
                // ie DF * (MAX(Exp FX - K, 0) * Exp DF)    for a call
                double expDF = sv.expDf->firstDF();
                const SVPath& path = *sv.expSpot;
                // at most only one date but possibly zero dates though
                double option = 0.0;
                for (int iStep = path.begin(); iStep < path.end(); ++iStep){
                    option = path[iStep]; /* At fwdStart, expected FX
                                            at fwdMat */
                }
                option = inst->isCall[i]? 
                    option - inst->strike[i]: (inst->strike[i] - option);
    #if 1
                // terse
                price += option > 0.0?
                    option * expDF * df: 0.0;
    #else
                // in full
                option = option > 0.0? option: 0.0;
                double optionTimesExpDf = option * expDF;
                double thisPrice = optionTimesExpDf * df;
                price += thisPrice;
    #endif
#ifdef DEBUGSRM    
                fprintf(stderr, "updatePrice: expSpot t= %d j= %d option= %g expDF= %g df= %g price= %g\n", t,j,option, expDF, df, price);
#endif // FIXME                
            }
        }
        return price;
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* ForwardOptionSeries::createProduct(const MonteCarlo* model) const {
     // the sim series here is irrelevant - to do: fix MCProductClient
    SimSeriesSP simSeries(new SimSeries(1));
    simSeries->addDates(DateTimeArray(1,
                                      fwdMat.back())); // so getLastDate() works
    return new ForwardOptionSeriesMC(this, simSeries);
}

CClassConstSP const ForwardOptionSeries::TYPE = CClass::registerClassLoadMethod(
    "ForwardOptionSeries", typeid(ForwardOptionSeries), load);

// * for class loading (avoid having header file) */
bool ForwardOptionSeriesLoad() {
    return (ForwardOptionSeries::TYPE != 0);
}

DRLIB_END_NAMESPACE
