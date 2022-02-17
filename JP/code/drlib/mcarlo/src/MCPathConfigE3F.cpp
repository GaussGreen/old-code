//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : MCPathConfigE3F.cpp
//
//   Description : MonteCarlo path generator for Energy 3 Factor model.
//                 
//
//   Author      : Spyridon Schismenos
//
//   Date        : June 22, 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/Random.hpp"
#include "edginc/Maths.hpp"
#include "edginc/MarketDataFetcherE3F.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/SVGenEnergyFuturePrice.hpp"
#include "edginc/EnergyInstVolExplicit.hpp"
#include "edginc/EnergyUnderlyer.hpp"
#include "edginc/DependenceGauss.hpp"

#define CAREFUL_RAND_DEBUG 0
#if CAREFUL_RAND_DEBUG
#include <fstream>
#endif

DRLIB_BEGIN_NAMESPACE

class MarketDataFetcherE3F2;
class MCPathConfigE3F : public MCPathConfig 
{

public:

    static CClassConstSP const TYPE;
    
    class Gen;
    friend class Gen;
    class  EnergyFutureSV;    // State Variable
    string     volType;
    EnergyInstVolBaseWrapper         volatility;
// S. Chen. MCProduct became IMCProduct
    virtual MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod );
    virtual MarketDataFetcherSP marketDataFetcher() const;
        
    virtual MCPathGeneratorSP futurePathGenerator(int                      cachingMode,
                                                  int                      numPaths,
                                                  const MCPathGeneratorSP& pastPathGenerator,
                                                  const IMCProduct*         prod,
                                                  Control*                 control, 
                                                  Results*                 results);

    virtual void getMarket(const IModel*        model,
                           const MarketData*    market,
                           CInstrument*   instrument);
    void getComponentMarketData(const IModel*          model,
                                const MarketData*      market,
                                MarketObjectSP         mo,
                                MarketDataFetcherE3F2* mdf) {}

protected:
    MCPathGeneratorSP makePathGenerator(bool                       cachingRequested,
                                        int                        numPaths,
                                        const MCPathGeneratorSP&   pastPathGenerator,
                                        const IMCProduct*           prod,
                                        Control*                   control, 
                                        Results*                   results,
                                        DateTimeArray&             simDates);

    virtual bool carefulRandoms() const { return true; }

private:

    virtual bool vegaMatrixSupported() const { return true; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * Not sure if we will ever want this?  See IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

    MCPathConfigE3F(const string& volType):MCPathConfig(TYPE), volType(volType) {}
    MCPathConfigE3F():MCPathConfig(TYPE) {}

    static void load(CClassSP& clazz)
    {
        REGISTER(MCPathConfigE3F, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(defaultMCPathConfigE3F);
        FIELD(volType, "Volatility Type");
        FIELD(volatility, "Volatility Wrapper");  // Here just for testing reasons. 
        // Definitely not correct approach
        clazz->setPublic();
    }

    static IObject* defaultMCPathConfigE3F() 
    {
        return new MCPathConfigE3F();
    }
};

CClassConstSP const MCPathConfigE3F::TYPE = 
CClass::registerClassLoadMethod("MCPathConfigE3F", typeid(MCPathConfigE3F), load);

class MCPathConfigE3F::EnergyFutureSV: public virtual SVGenEnergyFuturePrice::IStateVar
{

public:

    EnergyFutureSV(DoubleMatrix* spot) :_spot(spot){ };

    double getSpot()
    {
        return  (*_spot)[0][_spot->numRows() -1];
    }

    bool doingPast() const { return false; }

    DoubleMatrix getFullPath()
    {
        return *_spot;
    }

private:

    DoubleMatrix* _spot;

};


class MCPathConfigE3F::Gen: virtual public MCPathGen,
                       virtual public MCPathGenerator,
                       virtual public IStateVariableGen::IStateGen
{

private:

    int             pathIdx;
    IMCRandomSP     randomGen;
    DoubleMatrix    spot;
    DoubleMatrix    lnSpot;
    double          cSpot; //Price at time 0, has to come from market data
    StateVarDBase   svDBase;
    bool            havePast;
    double          alpha;
    double          beta;
    DoubleArray     sigma1;
    DoubleArray     sigma1Bar;
    DoubleArray     sigma2;
    DoubleArray     sigma2Bar;
    DoubleArray     timeSteps;
    DoubleArray     sqrtTimeSteps;
    DoubleArray     drifts;
    DoubleArray     shocks1;
    DoubleArray     shocks2;
    DoubleArray     timeToMaturity;
    DateTimeArray   simTimeline;
    DateTime        maturity;
    const IMultiMarketFactors*   mAsset;   // possibly useful in the future for getting market data

public:

    /** MCpathGenerator methods */
    // Deprecated methods
    virtual int NbSimAssets() const { return 0; }
    virtual const double* Path(int iAsset, int iPath) const { return 0; }
    virtual double refLevel(int iAsset, int iPath) const { return 0; }
    virtual double maxDriftProduct(int iAsset) const { return 0; }
    virtual int begin(int iAsset) const { return 0; }
    virtual int end(int iAsset) const { return 0; }
    
    // live methods
    bool hasPast() const 
    {
        return havePast;
    }

    int getPathIndex() const 
    {
        return pathIdx;
    }

    virtual void generatePath(int pathIdx);

    virtual bool doingPast() const
    {
        return false;
    }
    virtual IStateVariableSP create(const IStateVariableGen* svGen) 
    {
        return svDBase.find(svGen);
    }

    Gen(int numSimPaths, MCPathConfigE3F* mcPathConfig,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient* prodClient,
        EnergyInstVolBaseSP vol);

    void calibrate(); //help method that initializes useful quantities before the
                      // simulation loop starts
};

MCPathConfigE3F::Gen::Gen(int numSimPaths, MCPathConfigE3F* mcPathConfig,
                          const MCPathGeneratorSP& pastPathGenerator,
                          const MCProductClient* prodClient,
                          EnergyInstVolBaseSP vol) :
    pathIdx(-1),
    havePast(pastPathGenerator->hasPast()),
    mAsset(prodClient->getMultiMarketFactors())
{
    StateVariableCollectorSP svCollector(new StateVariableCollector());
    prodClient->collectStateVars(svCollector);
    IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

    std::vector<const SVGenEnergyFuturePrice*> efsv(filterStateVars<SVGenEnergyFuturePrice>(stateVarGenArray));
    for (unsigned iVar = 0; iVar < efsv.size(); ++iVar)
    {
        // create EnergyFuturePrice state vars and put in database
        IStateVariableSP sv(new MCPathConfigE3F::EnergyFutureSV(&spot));
        svDBase.append(efsv[iVar], sv);
    }

    // Only Energy Future State variables are supported. Throw an exception if 
    // another request exists

    if(!stateVarGenArray.empty()) {
        throw ModelException("Unable to recognize all state variable types.");
    }

    // Take the expiry dates, and create the simulation timeline by disregarding the ones that 
    // are later from the maturity date.
    simTimeline=vol->getFutureMaturityDates(); //simTimeline=vol->getOptionExpiryDates();
    maturity=prodClient->getSimSeries()->getLastDate();
    for (int count=simTimeline.size()-1; count>=0; count--)
    {
        if (simTimeline[count].isGreater(maturity)) {
            simTimeline.pop_back();
        }
    }

    // Initialize diffusion quantities.
    spot=DoubleMatrix(1,simTimeline.size());
    lnSpot=DoubleMatrix(1,simTimeline.size());
    spot[0][0]=cSpot=43.6;
    lnSpot[0][0]=log(43.6);

//      Get parameters from volatility
    alpha=(vol->getAlphas())[0];
    beta=(vol->getAlphas())[1];

    int numExpiries = vol->getFutureMaturityDates().size(); //int numExpiries = vol->getOptionExpiryDates().size();
    const DoubleMatrix& locSigmas = vol->getSigmas();
    for ( int i=0; i<numExpiries; i++)
    {
        sigma1[i]=locSigmas[i][0];
        sigma2[i]=locSigmas[i][3];
        sigma1Bar[i]=locSigmas[i][1];
        sigma2Bar[i]=locSigmas[i][2];
    }
    this->calibrate();


    // Create the (identity) correlation matrix
    CDoubleMatrixSP corrs(new CDoubleMatrix(2,2));
    (*corrs)[0][0]=1.0;
    (*corrs)[1][0]=0.0;
    (*corrs)[0][1]=0.0;
    (*corrs)[1][1]=1.0;
    // Declare that we want Normal random variables
    // S.Chen. fromSqr.. not supported anymore,so...   DependenceSP dependence(Gauss::fromSqrtOfCorrelations(*corrs));
    DependenceSP dependence(new Gauss(*corrs));
     
    // Create the random number generator
    randomGen = IMCRandomSP(new MCRandomNoCache(
                                0,
                                dependence,
                                mcPathConfig->getRandomGenerator(),
                                mcPathConfig->carefulRandoms(),
                                timeSteps.size(), // how many steps?  
                                2, // num of factors
                                0)); // number of past dates, does not seem to be used any where

}

MCPathGeneratorSP MCPathConfigE3F::pastPathGenerator(const IMCProduct* prod )
{
    const MCProductClient* prodClient=
        dynamic_cast<const MCProductClient*>(prod);
    return MCPathGeneratorSP(new PastPathGen(prodClient));
}

MCPathGeneratorSP MCPathConfigE3F::makePathGenerator(bool  cachingRequested,
                                                     int numPaths,
                                                     const MCPathGeneratorSP& pastPathGenerator,
                                                     const IMCProduct* prod,
                                                     Control* control, 
                                                     Results* results,
                                                     DateTimeArray&  simDates) 
{
    try 
    {
        const MCProductClient* prodClient = 
            &dynamic_cast<const MCProductClient&>(*prod);
        return MCPathGeneratorSP(new Gen(numPaths, this, pastPathGenerator,
                                         prodClient, volatility.getSP()));
    } 
    catch (exception& e)
    {
        throw ModelException(e, "MCPathConfigE3F::makePathGenerator");
    }
}

MCPathGeneratorSP MCPathConfigE3F::futurePathGenerator(int cachingMode,
                                                       int numPaths,
                                                       const MCPathGeneratorSP& pastPathGenerator,
                                                       const IMCProduct* prod,
                                                       Control* control, 
                                                       Results* results)
{
    static const string method("MCPathConfigE3F::futurePathGenerator");
    if (cachingMode & IMCProduct::CACHE_PRODUCT_BIT)
    {
        throw ModelException(method, "Caching not supported yet");
    }
    // S. Chen. Note, new makePathGenerator takes one more param, simDates
    DateTimeArray dummy;
    return makePathGenerator(false, numPaths, pastPathGenerator, 
                             prod, control, results, dummy);
}

void MCPathConfigE3F::Gen::generatePath(int pathIdx)
{
    this->pathIdx=pathIdx;
    randomGen->generate(pathIdx);
    // get the random numbers
    // randomNums[num of factor][num of step]
    const DoubleMatrix& randomNums = randomGen->getRandomNumbers();
    // actual diffusion logic
    for (int i=0; i< simTimeline.size()-1; i++) {
        lnSpot[0][i+1]=lnSpot[0][i]+drifts[i]+shocks1[i]*randomNums[0][i]+shocks2[i]*randomNums[1][i];
        spot[0][i+1]=exp(lnSpot[0][i+1]);
    }
}

void MCPathConfigE3F::getMarket(const IModel* model,
                                const MarketData*    market,
                                CInstrument*   instrument)
{
    static const string method("MCPathConfigE3F::getMarket");
    try 
    {
    } 
    catch (exception& e)
    {
        throw ModelException(e,method);
    }

}

void MCPathConfigE3F::Gen::calibrate() 
{

    //    Help method that initializes quantities needed in the diffusion
    timeSteps.reserve(simTimeline.size()-1);
    timeToMaturity.reserve(simTimeline.size());
    sqrtTimeSteps.reserve(simTimeline.size()-1);
    shocks1.reserve(simTimeline.size()-1);
    shocks2.reserve(simTimeline.size()-1);
    drifts.reserve(simTimeline.size()-1);
    double drift,shock1,shock2;
    timeToMaturity.push_back((double)maturity.daysDiff(simTimeline[0])/365.0);
    for (int i = 0; i<simTimeline.size()-1 ; i++)
    {   
        drift=shock1=shock2=0.0;
        timeSteps.push_back((double)simTimeline[i+1].daysDiff(simTimeline[i])/365.0);
        sqrtTimeSteps.push_back(sqrt(timeSteps[i]));
        timeToMaturity.push_back((double)maturity.daysDiff(simTimeline[i+1])/365.0);
        shock1=sigma1Bar[i]*sigma1Bar[i]*timeSteps[i];
        shock1+=0.5*sigma1[i]*sigma1[i]*(exp(-2.0*alpha*timeToMaturity[i+1])-exp(-2.0*alpha*timeToMaturity[i]))/alpha;
        shock1+=sigma1[i]*sigma1Bar[i]*(exp(-alpha*timeToMaturity[i+1])-exp(-alpha*timeToMaturity[i]))/alpha;
        shock1+=2.0*sigma1[i]*sigma2[i]*(exp(-(alpha+beta)*timeToMaturity[i+1])-exp(-(alpha+beta)*timeToMaturity[i]))/(alpha+beta);
        shock1+=sigma2[i]*sigma1Bar[i]*(exp(-beta*timeToMaturity[i+1])-exp(-beta*timeToMaturity[i]))/beta;
        shock1+=0.5*sigma2[i]*sigma2[i]*(exp(-2.0*beta*timeToMaturity[i+1])-exp(-2.0*beta*timeToMaturity[i]))/beta;
        shocks1.push_back(sqrt(shock1));
        shock2=sigma2Bar[i]*timeSteps[i]*sigma2Bar[i];
        shocks2.push_back(sigma2Bar[i]*sqrtTimeSteps[i]);
        drift=-0.5*(shock1+shock2);
        drifts.push_back(drift);
    }
}
//// unclear why this class exists - unless still in development
class MarketDataFetcherE3F2: public MarketDataFetcherE3F
{

public:

    MarketDataFetcherE3F2(const string& volType,MCPathConfigE3F* pathconfig):
        MarketDataFetcherE3F(volType)
    {
#if 0
        pathConfig=pathconfig;
#endif
    }
#if 0
    virtual MarketObjectSP fetch(const MarketData *market,
                                 const string& name,
                                 const CClassConstSP& type,
                                 const IModel *model) const {
        static const string   method("MarketDataFetcherE3F::fetch");
        CClassConstSP         typeToUse = type;
        try {

            // Check that the volatility passed is a Energy Vol
            // otherwise throw exception
            if (volClass ==0) {
                volClass = CClass::forName(volType);
                if (!EnergyVolBase::TYPE->isAssignableFrom(volClass)) {
                    throw ModelException(method, "Specified volatility type (" +
                                         volType + ") is not an Energy Volatility");
                } 
            }
            return MarketDataFetcher::fetch(market, name, typeToUse, model);
        }
        catch (exception& e) {
            throw ModelException(e, method, "Failed to get market data " + name +
                                 " of type " + typeToUse->getName()); 
        }
    }

    virtual void getComponentMarketData(const IModel* model, const MarketData* market, MarketObjectSP mo) const
    {
        //just call the parent's method
        MarketDataFetcher::getComponentMarketData(model,market,mo);
    }

private:

    MCPathConfigE3F* pathConfig; // Goal is to pass in the MCPathConfigE3F object in order to populate
                                 // it with the Market Data 
#endif
};

MarketDataFetcherSP MCPathConfigE3F::marketDataFetcher() const
{
    return MarketDataFetcherSP(new MarketDataFetcherE3F2(volType,const_cast<MCPathConfigE3F*>(this)));
}
         
bool MCPathConfigE3FLoad()
{
    return (MCPathConfigE3F::TYPE != 0);
}

DRLIB_END_NAMESPACE
