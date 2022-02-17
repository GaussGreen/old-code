//----------------------------------------------------------------------------
//
//   Description : demo path config for monte carlo 
//
//   Author      : Jay Z Wang
//
//-----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MarketDataFetcherDemo.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/SVGenDemo.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductClient.hpp"
#include "edginc/DependenceGauss.hpp"

#define CAREFUL_RAND_DEBUG 0
#if CAREFUL_RAND_DEBUG
#endif

DRLIB_BEGIN_NAMESPACE

class MCPathConfigDemo : public MCPathConfig {
public:
    static CClassConstSP const TYPE;

    class Gen;
    friend class Gen;

    class NewGen;
    friend class NewGen;

    class JSV;
    class NewJSV;

    // fields, unused, for demo purpose ////
    int x;
    int y;

    virtual MCPathGeneratorSP pastPathGenerator(const IMCProduct* prod)
    {
        const MCProductClient* prodClient =
            dynamic_cast<const MCProductClient*>(prod);
        return MCPathGeneratorSP(new PastPathGen(prodClient));
    }

    virtual MarketDataFetcherSP marketDataFetcher() const 
    {
		return MarketDataFetcherSP(new MarketDataFetcherDemo());
    }

    virtual MCPathGeneratorSP futurePathGenerator(int cachingMode, 
                                                  int numPaths, 
                                                  const MCPathGeneratorSP& pastPathGenerator, 
                                                  const IMCProduct* prod, 
                                                  Control* control, 
                                                  Results* results,
                                                  DateTimeArray& simDates );

private:

    // for reflection
    MCPathConfigDemo(): MCPathConfig(TYPE)
    {
    }

    virtual bool vegaMatrixSupported() const { return true; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this pdf
     *
     * See IModel::wantsRiskMapping().
     */

    IModel::WantsRiskMapping wantsRiskMapping() const {
        return IModel::riskMappingIrrelevant;
    }

protected:

    /** Creates a future path generator */
    MCPathGeneratorSP makePathGenerator(
        bool                     cachingRequested,
        int                      numPaths,
        const MCPathGeneratorSP& pastPathGenerator,
        const IMCProduct*         prod,
        Control*                 control, 
        Results*                 results,
        DateTimeArray&           simDates);
    
    virtual bool carefulRandoms() const { return false; }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(MCPathConfigDemo, clazz);
        SUPERCLASS(MCPathConfig);
        EMPTY_SHELL_METHOD(defaultMCPathConfigDemo);
        FIELD(x, "x");
        FIELD(y, "y");
    }

    static IObject* defaultMCPathConfigDemo() {
        return new MCPathConfigDemo();
    }
};

CClassConstSP const MCPathConfigDemo::TYPE = 
CClass::registerClassLoadMethod("MCPathConfigDemo", typeid(MCPathConfigDemo), load);

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// State Variables
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

class MCPathConfigDemo::JSV: public virtual SVGenDemo::IStateVar{
public:
    JSV(double* spot) : _spot(spot) {};

    double* getSpot()
    {
        return _spot;
    }

    bool doingPast() const { return false; }

private:
    double* _spot;
};

class MCPathConfigDemo::NewJSV: public virtual SVGenDemo::INewStateVar{
public:
    NewJSV(double* data) : _index(0),_data(data) {};

    double getSpot()
    {
        return _data[_index];
    }

    bool doingPast() const { return false; }

    void reset() { _index = 0; }

    void advance() { ++_index; }

private:
    int _index;
    double* _data;
};


class MCPathConfigDemo::Gen: virtual public MCPathGenerator,
                             virtual public MCStatelessPathGen,
                             virtual public IStateVariableGen::IStateGen,
                             public DependenceMakerGauss::Support {

private:
    int pathIdx;
    IMCRandomSP randomGen;
    bool havePast;

    double spot[20]; // max 20 slices diffused state

    StateVarDBase svDBase; // state var

    MCPathConfigDemo::NewJSV* _sv; // XXX HACK


public:

    /** MCpathGenerator methods */
    // Deprecated methods
    virtual int NbSimAssets() const { return 0; }

    virtual const double* Path(int /*iAsset*/, int /*iPath*/) const { return 0; }
    
    virtual double refLevel(int /*iAsset*/, int /*iPath*/) const { return 0; }

    virtual double maxDriftProduct(int /*iAsset*/) const { return 0; }

    virtual int begin(int /*iAsset*/) const { return 0; }

    virtual int end(int /*iAsset*/) const { return 0; }

    // live methods
    bool hasPast() const {
        return havePast;
    }

    int getPathIndex() const {
        return pathIdx;
    }

    virtual void generatePath(int pathIdx){
        this->pathIdx = pathIdx;
        // get random numbers
        randomGen->generate(pathIdx);
        const DoubleMatrix& randomNums = randomGen->getRandomNumbers();
        
        {
            // actual diffusion logic
            const double* randoms = randomNums[0];
            double deltaT = 0.1;
            double r = 0.1;
            double vol = 0.2;

            // iterate through timepoints
            std::cout << "Path " << pathIdx << ": ";
            for (int i = 1; i < 20; ++i) {
                // ds = s * r * dt + vol * s * dw
                spot[i] = spot[i-1] + spot[i-1] * r * deltaT + spot[i-1] * vol * randoms[i];
                std::cout << spot[i] << " ";
            }   
            std::cout << "\n";
        }

        _sv->reset();
    }

    virtual bool doingPast() const {
        return false;
    }
    
    /** constructor*/
    Gen(int                         /*numSimPaths*/,
        MCPathConfigDemo*           mcPathConfig,
        const MCPathGeneratorSP&    pastPathGenerator,
        const MCProductClient*      prodClient):
        pathIdx(-1), havePast(pastPathGenerator->hasPast()) 
    {
        // populate state var database
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

        // we only understand SVGenDemo
        std::vector<const SVGenDemo*> jsv(filterStateVars<SVGenDemo>(stateVarGenArray));
        for (unsigned iVar = 0; iVar < jsv.size(); ++iVar)
        {
            // create individual state vars and put in database
            // for path payoff , create JSV
            // IStateVariableSP sv(new MCPathConfigDemo::JSV(spot));
            // for stateless payoff, create NewJSV
            IStateVariableSP sv(_sv = new MCPathConfigDemo::NewJSV(spot));
            svDBase.append(jsv[iVar], sv);
        }

        // Initialize random number generator for paths
        DependenceSP dependence = mcPathConfig->dependenceMaker->createDependence(this);
        randomGen = IMCRandomSP(new MCRandomNoCache(
            0,
            dependence,
            mcPathConfig->getRandomGenerator(),
            mcPathConfig->carefulRandoms(),
            20, // how many steps?
            1, // num assets
            0)); // number of past dates, does not seem to be used any where

        spot[0] = 30; // S(0)
    }

    void advance() {
        // advance all my state variables
        _sv->advance();
    }

    /** Returns the state variable corresponding to generator.
    Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen){
        return svDBase.find(svGen);
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigDemo::Gen::getGaussData");
        try {
            CDoubleMatrix inputMatrix(1,1);
            inputMatrix[0][0] = 1.0;
            return CDoubleMatrixSP(new DoubleMatrix(inputMatrix));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


class MCPathConfigDemo::NewGen: virtual public MCPathGenerator,
                                virtual public MCStatelessSliceGen,
                                virtual public IStateVariableGen::IStateGen,
                                public DependenceMakerGauss::Support {

private:
    int sliceIdx;
    IMCRandomSP randomGen;
    bool havePast;

    int numPaths;

    double* spot; // diffused state

    StateVarDBase svDBase; // state var

    MCPathConfigDemo::NewJSV* _sv; // XXX HACK

public:

    /** MCpathGenerator methods */
    // Deprecated methods
    virtual int NbSimAssets() const { return 0; }

    virtual const double* Path(int /*iAsset*/, int /*iPath*/) const { return 0; }
    
    virtual double refLevel(int /*iAsset*/, int /*iPath*/) const { return 0; }

    virtual double maxDriftProduct(int /*iAsset*/) const { return 0; }

    virtual int begin(int /*iAsset*/) const { return 0; }

    virtual int end(int /*iAsset*/) const { return 0; }

    // XXX
    virtual void generatePath(int) {}
    virtual bool hasPast() const { return false; }
    virtual int getPathIndex() const { return 0; }

    // live methods
    virtual void generateSlice(int sliceIdx){
        this->sliceIdx = sliceIdx;
        // get random numbers
        randomGen->generate(sliceIdx);
        const DoubleMatrix& randomNums = randomGen->getRandomNumbers();
        
        {
            // actual diffusion logic
            const double* randoms = randomNums[0];
            double deltaT = 0.1;
            double r = 0.1;
            double vol = 0.2;

            // iterate through timepoints
            std::cout << "Slice " << sliceIdx << ": ";
            for (int i = 0; i < numPaths; ++i) {
                // ds = s * r * dt + vol * s * dw
                spot[i]  += spot[i] * r * deltaT + spot[i] * vol * randoms[i];
                std::cout << spot[i] << " ";
            }   
            std::cout << "\n";
        }

        // reset all my state variables
        _sv->reset();

    }

    void advance() {
        // advance all my state variables
        _sv->advance();
    }

    virtual bool doingPast() const {
        return false;
    }

    NewGen(int                      numSimPaths,
        MCPathConfigDemo*        mcPathConfig,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient):
        sliceIdx(-1), havePast(pastPathGenerator->hasPast()) 
    {
        numPaths = numSimPaths;
        spot = new double[numPaths];

        // populate state var database
        StateVariableCollectorSP svCollector(new StateVariableCollector());
        prodClient->collectStateVars(svCollector);
        IElemStateVariableGenArray stateVarGenArray = svCollector->getElemStateVarGens();

        // we only understand SVGenDemo
        std::vector<const SVGenDemo*> jsv(filterStateVars<SVGenDemo>(stateVarGenArray));
        for (unsigned iVar = 0; iVar < jsv.size(); ++iVar)
        {
            // create individual state vars and put in database
            IStateVariableSP sv(_sv = new MCPathConfigDemo::NewJSV(spot));
            svDBase.append(jsv[iVar], sv);
        }

        // Initialize random number generator for paths
        DependenceSP dependence = mcPathConfig->dependenceMaker->createDependence(this);
        randomGen = IMCRandomSP(new MCRandomNoCache(
            0,
            dependence,
            mcPathConfig->getRandomGenerator(),
            mcPathConfig->carefulRandoms(),
            numSimPaths, // how many paths
            1, // num assets
            0)); // number of past dates, does not seem to be used any where

        for (int i = 0; i < numSimPaths; ++i)
            spot[i] = 30; // S(0)

    }

    /** Returns the state variable corresponding to generator.
    Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen){
        return svDBase.find(svGen);
    }

    // implement getGaussData since derivation from DependenceMakerGauss::Support
    virtual CDoubleMatrixConstSP getGaussData() const {
        static const string method("MCPathConfigDemo::NewGen::getGaussData");
        try {
            CDoubleMatrix inputMatrix(1,1);
            inputMatrix[0][0] = 1.0;
            return CDoubleMatrixSP(new DoubleMatrix(inputMatrix));
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    ~NewGen() { delete[] spot; }
};


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


MCPathGeneratorSP MCPathConfigDemo::futurePathGenerator(
    int                      cachingMode,
    int                      numPaths,
    const MCPathGeneratorSP& pastPathGenerator,
    const IMCProduct*         prod,
    Control*                 control, 
    Results*                 results,
    DateTimeArray&           simDates )
{
    static const string method("MCPathConfigDemo::futurePathGenerator");
    if (cachingMode & IMCProduct::CACHE_PRODUCT_BIT){
        throw ModelException(method, "Caching not supported yet");
    }
    return makePathGenerator(false, numPaths, pastPathGenerator, 
                             prod, control, results, simDates);
}

/** Essentially a pass through for futurePathGenerator except that the
    relevant caches are created/updated etc */
MCPathGeneratorSP MCPathConfigDemo::makePathGenerator(
    bool                               /*cachingRequested*/,
    int                                numPaths,
    const MCPathGeneratorSP&           pastPathGenerator,
    const IMCProduct*                   prod,
    Control*                           /*control*/, 
    Results*                           /*results*/,
    DateTimeArray&                     /*simDates*/ ) {
    try {
        const MCProductClient* prodClient = 
            &dynamic_cast<const MCProductClient&>(*prod);

        // for slice simulation, use NewGen
        /*
        return MCPathGeneratorSP(new NewGen(numPaths, this, pastPathGenerator,
                                         prodClient));
        */
        // for path simulation, use Gen
        return MCPathGeneratorSP(new Gen(numPaths, this, pastPathGenerator,
                                         prodClient));

    } catch (exception& e){
        throw ModelException(e, "MCPathConfigSRM::makePathGenerator");
    }
}



// external symbol to allow class to be forced to be linked in
bool MCPathConfigDemoLoad(){
    return (MCPathConfigDemo::TYPE != 0);
}

DRLIB_END_NAMESPACE
