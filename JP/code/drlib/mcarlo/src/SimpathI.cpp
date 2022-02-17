//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimpathI.cpp
//
//   Description : SimpathI model
//

//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/PastPathGenerator.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/SimpathI.hpp"
#include "edginc/Array.hpp"
#include "edginc/MCProduct.hpp"
#include "edginc/MCProductEngineClient.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/IReportResults.hpp"
#include "edginc/MemoryProfiler.hpp"

#include "edginc/Instrument.hpp"
#include "edginc/MarketDataFetcher.hpp" 
#include "edginc/MCPathGenerator.hpp"
#include "edginc/MCPrices.hpp"
#include "edginc/XMLWriter.hpp"
DRLIB_BEGIN_NAMESPACE

SimpathI::~SimpathI()
{ delete instrumentMC; }

SimpathI::SimpathI(CClassConstSP clazz):
    CModel(clazz), instrumentMC(0), nbIter(0) {}

SimpathI::SimpathI():
    CModel(TYPE), instrumentMC(0), nbIter(0){}

void SimpathI::serialize(string fn) {
	XMLWriter xml(fn);
	this->write("SimpathI", &xml);
}
bool SimpathI::Initialize()
{
    static const string method("SimpathI::Initialize");
#if 1
	serialize("simpathi.xml");
#endif
    try {
        MCPathGeneratorSP pastPathGenerator(new DummyPastPathGenerator());

        // load market data

        IInstrumentCollectionSP instrumentCollection(IInstrumentCollection::singleton(instrument));


        instrumentCollection->GetMarket(this, market);

        pathConfig->getMarket(this, market.get(), instrumentCollection);

        const IMCIntoProduct* intoProd =
            dynamic_cast<const IMCIntoProduct*>(instrument.get());
        if (intoProd == 0)
            throw ModelException("Not a monte carlo instrument");
        instrumentMC = intoProd->createProduct(NULL);

        DateTimeArray   dates;

        pathGen = pathConfig->futurePathGenerator(
            0,
            nbIter, 
            pastPathGenerator,
            instrumentMC,
            NULL, 
            NULL, // control & result not used for PathConfigSRM
            dates);

        // Cause all state variables to be collected.
        instrumentMC->pathGenUpdated(pathGen.get());

    } catch (exception& e) {
        throw ModelException(e, method, "SimpathI initialization failed!");
    }


    return true;
}

 IObjectSP SimpathI::GeneratePath(int idx)
{
    IMCPrices * prices = instrumentMC->createOrigPrices(1, 1, 0);
    return GeneratePath(idx, *prices);
}

IObjectSP SimpathI::GeneratePath(int idx, IMCPrices& prices)
{
    static const string method("SimpathI::GeneratePath");

    try
    {
        pathGen->generatePath(idx);
        instrumentMC->payoff(pathGen.get(), prices);
        if (idx == nbIter-1)
            instrumentMC->finalize();
    }
    catch (exception& e)
    {
        throw ModelException(e, method, "SimpathI generate path failed!");
    }

	IGenericPrices * gp = dynamic_cast<IGenericPrices*>(& prices);
	
	return (gp == NULL) ? IObjectSP(CBool::create(true)) : gp->getResult();
	
/*    // extra hacky code to return results suitable for RegTesting framework
    IReportResultsLite * report = dynamic_cast<IReportResultsLite *>(instrumentMC);
    
    // Convert address to a string: find more portable way later.
    std::ostringstream hack;
    hack << report;
    
    return IObjectSP(CString::create(hack.str()));*/
    
}

MarketDataFetcherSP SimpathI::createMDF() const
{
  return pathConfig->marketDataFetcher();
}

/** Invoked when Class is 'loaded' */
void SimpathI::load(CClassSP& clazz){
    REGISTER(SimpathI, clazz);
    SUPERCLASS(CModel);
    IMPLEMENTS(ClientRunnable),
    EMPTY_SHELL_METHOD(defaultSimpathI);
    FIELD(market, "market");
    FIELD(pathConfig, "Path generator");
    FIELD(instrument, "SimpathIInstrument");
    FIELD(nbIter, "NbIter");
    clazz->setPublic();
}

IObject* SimpathI::defaultSimpathI(){
    return new SimpathI();
}

CClassConstSP const SimpathI::TYPE = CClass::registerClassLoadMethod(
    "SimpathI", typeid(SimpathI), SimpathI::load);

int SimpathI::getNbIter() {
    return nbIter;
}

// A sample implementation to make SimpathI testable via models.exe framework
// As we derive from IRegressionTest interface, as runTest() method should be constant, so we copy the object.

IObjectSP SimpathI::runTest() const
{
    SimpathISP obj(copy(this)); // copy object
	obj->serialize("c:/simapthi.xml");
    obj->Initialize();

    // TODO: move the code in a separate function;
    // Add code to print results of diffusions
    IMCPrices * prices = obj->instrumentMC->createOrigPrices(nbIter, 1, 0);
    
    for(int i=0; i+1 < obj->getNbIter(); ++i) {
        IObjectSP res = obj->GeneratePath(i, *prices);
    }
    
    obj->instrumentMC->finalize();
    IGenericPrices * gp = dynamic_cast<IGenericPrices*>(prices);
    return (gp == NULL) ? IObjectSP(CBool::create(true)) : gp->getResult();
    
/*    IReportResults * report = dynamic_cast<IReportResults *>(obj->instrumentMC);
    if (report != NULL)
        return report->getResults();
    
    return IObjectSP(CBool::create(true));*/
}

/////////////////////////////////////////////////////////////////

SimpathIController::SimpathIController(CClassConstSP clazz):
    CObject(clazz), model(0), idx(0) {}

SimpathIController::SimpathIController():
    CObject(TYPE), model(0), idx(0){}

bool SimpathIController::Initialize()
{
    return model->Initialize();
}

IObjectSP SimpathIController::GeneratePath()
{
    return model->GeneratePath(idx);
}

IObjectSP SimpathIController::runTest()
{
	return model->runTest();
}

/** Invoked when Class is 'loaded' */
void SimpathIController::load(CClassSP& clazz){
    REGISTER(SimpathIController, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSimpathIController);
    FIELD(model, "SimpathI model");
    FIELD(idx, "path id");
    FIELD_MAKE_OPTIONAL(idx);
    clazz->setPublic();
    // make visible to ES/spreadsheet
    Addin::registerObjectMethod("GENERATE_PATH",
        Addin::UTILITIES,
        "Generate one path",
        false, //don't bother to read in handle name
        Addin::expandSimple,
        &SimpathIController::GeneratePath);

    Addin::registerBoolMethod(
        "SIMPATHL_INITIALIZE",
        Addin::UTILITIES,
        "Initialize the simulation",
        &SimpathIController::Initialize);

	Addin::registerObjectMethod(
		"GENERATE_PATHS",
		Addin::UTILITIES,
		"Generate a bunch of paths",
		false, //don't bother to read in handle name
		Addin::expandSimple,
		&SimpathIController::runTest);
}

IObject* SimpathIController::defaultSimpathIController(){
    return new SimpathIController();
}

CClassConstSP const SimpathIController::TYPE = CClass::registerClassLoadMethod(
    "SimpathIController", typeid(SimpathIController), SimpathIController::load);

// * for class loading (avoid having header file) */
bool SimpathILoad() {
    return (SimpathI::TYPE != 0 && SimpathIController::TYPE != 0);
}

class SimpathIInitializer : public CObject, public ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~SimpathIInitializer() {};
    IObjectSP run();

protected:
    SimpathIInitializer(CClassConstSP clazz);
    SimpathI* model;

private:
    SimpathIInitializer();
    static IObject* defaultSimpathIInitializer();
    static void load(CClassSP& clazz);
};

SimpathIInitializer::SimpathIInitializer(CClassConstSP clazz):
CObject(clazz), model(0) {}

SimpathIInitializer::SimpathIInitializer():
CObject(TYPE), model(0) {}

IObject* SimpathIInitializer::defaultSimpathIInitializer()
{
    return new SimpathIInitializer();
}

void SimpathIInitializer::load(CClassSP& clazz)
{
    REGISTER(SimpathIInitializer, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultSimpathIInitializer);
    FIELD(model, "SimpathI model");
    clazz->setPublic();
}

CClassConstSP const SimpathIInitializer::TYPE = CClass::registerClassLoadMethod(
    "SimpathIInitializer", typeid(SimpathIInitializer), SimpathIInitializer::load);

// * for class loading (avoid having header file) */
bool SimpathIInitializerLoad() {
    return (SimpathIInitializer::TYPE != 0);
}

IObjectSP SimpathIInitializer::run()
{
    bool result = model->Initialize();
    return IObjectSP(CBool::create(result));
}



class SimpathIPathGenerator : public CObject, public ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~SimpathIPathGenerator() {};
    IObjectSP run();

protected:
    SimpathIPathGenerator(CClassConstSP clazz);
    SimpathI* model;
    int idx;

private:
    SimpathIPathGenerator();
    static IObject* defaultSimpathIPathGenerator();
    static void load(CClassSP& clazz);
};

SimpathIPathGenerator::SimpathIPathGenerator(CClassConstSP clazz):
CObject(clazz), model(0), idx(0) {}

SimpathIPathGenerator::SimpathIPathGenerator():
CObject(TYPE), model(0), idx(0){}

IObject* SimpathIPathGenerator::defaultSimpathIPathGenerator()
{
    return new SimpathIPathGenerator();
}

void SimpathIPathGenerator::load(CClassSP& clazz)
{
    REGISTER(SimpathIPathGenerator, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultSimpathIPathGenerator);
    FIELD(model, "SimpathI model");
    FIELD(idx, "path id");
    FIELD_MAKE_OPTIONAL(idx);
    clazz->setPublic();
}

CClassConstSP const SimpathIPathGenerator::TYPE = CClass::registerClassLoadMethod(
    "SimpathIPathGenerator", typeid(SimpathIPathGenerator), SimpathIPathGenerator::load);

// * for class loading (avoid having header file) */
bool SimpathIPathGeneratorLoad() {
    return (SimpathIPathGenerator::TYPE != 0);
}


IObjectSP SimpathIPathGenerator::run()
{
    IObjectSP result = model->GeneratePath(idx);
    return result;
}

///////////////////////////////////////////////////////////////////////


   

////////////////////////////////////////////////////// SimpathIShell
CClassConstSP const SimpathIShell::TYPE = CClass::registerClassLoadMethod(
    "SimpathIShell", 
    typeid(SimpathIShell), 
    & SimpathIShell::load);

IObjectSP SimpathIShell::run() {
    IObjectSP res = simpathi->run();
	return res;
}

SimpathIShell::SimpathIShell(CClassConstSP clazz) :
        CObject(clazz),
        simpathi(NULL)
{
}

void SimpathIShell::load(CClassSP& clazz)
{
    REGISTER(SimpathIShell, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultSimpathIShell);
    FIELD(simpathi, "handle to SimpathI");
    clazz->setPublic();
    
    Addin::registerObjectMethod(
        "SIMPATHL_EXECUTE",
        Addin::UTILITIES,
        "Execute simpathi->run() and return results",
        false, //don't bother to read in handle name
        Addin::expandSimple,
        &SimpathIShell::run);

}

SimpathIShell::SimpathIShell() :
        CObject(TYPE),
        simpathi(NULL)
{
}

IObject* SimpathIShell::defaultSimpathIShell()
{
    return new SimpathIShell;
}

DRLIB_END_NAMESPACE

