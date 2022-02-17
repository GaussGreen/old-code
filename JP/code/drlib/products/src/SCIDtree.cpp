/**
 * @file SCIDtree.hpp
 */

#include "edginc/config.hpp"
#include "edginc/SCIDtree.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/CInstrumentCollection.hpp"


DRLIB_BEGIN_NAMESPACE

SCIDtree::~SCIDtree() {}

SCIDtree::SCIDtree() : Model(TYPE), MCmethod("DEFAULT") {}

SCIDtree::SCIDtree(CClassConstSP clazz): Model(clazz) {}

IModel::WantsRiskMapping SCIDtree::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** calculate single price and store result in CResult */
void SCIDtree::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "SCIDtree::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support SCIDtree::IntoProduct");
    }
    IProduct*     product = 0;
    try{
        product = intoProd->createProduct(this);
        product->price(this, control, results);
    } catch (exception& e){
        delete product;
        throw ModelException(e, method);
    }
    delete product;
}

void SCIDtree::getMarket(const MarketData*  market,  IInstrumentCollectionSP instruments) {
	SCIDtreeParam.getData(this, market); // get the SCIDtree parameters
    CModel::getMarket(market,instruments);
}

class SCIDtreeHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SCIDtree, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultSCIDtree);

		FIELD(MCmethod,					 "for testing: MC method to use: MARKET, CLOSED FORM, FAST, FULL, FULL1, FULL2, OR HYBRID");
		FIELD_MAKE_OPTIONAL(MCmethod);
		FIELD(freqFastMC,				 "how often do we compute the loss distribution");
		FIELD(freqFullMC,				 "how often do we simulate spread and losses");
		FIELD(seed,						 "seed for the simulation");
		FIELD(timeSteps,                 "time Steps used in the discretization of the CIR process");
 		FIELD(nbPathsFast,               "nbPaths conditional on no jumps and conditional on at least one jump");
  		FIELD(nbPathsFull,               "nbPaths used in the Full MC");
		FIELD(ConvolutionMethod,         "Convolution Method conditional on no jumps and conditional on at least one jump");
 		FIELD(SCIDtreeParam,             "sCID parameters");
		FIELD(test,						 "for testing, to remove");
		FIELD_MAKE_OPTIONAL(test);

    }




    static IObject* defaultSCIDtree() {
        return new SCIDtree();
    }
};

CClassConstSP const SCIDtree::TYPE = CClass::registerClassLoadMethod(
    "SCIDtree", typeid(SCIDtree), SCIDtreeHelper::load);

// for linker
bool   SCIDtreeLoad() {
    return (SCIDtree::TYPE != 0);
   }

CClassConstSP const SCIDtree::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("SCIDtree::IIntoProduct",
                                    typeid(SCIDtree::IIntoProduct), 0);

DRLIB_END_NAMESPACE
