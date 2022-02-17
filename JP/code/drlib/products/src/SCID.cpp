/**
 * @file SCID.hpp
 */

#include "edginc/config.hpp"
#include "edginc/SCID.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/CInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

SCID::~SCID() {}

SCID::SCID() : Model(TYPE), MCmethod("DEFAULT") {}

SCID::SCID(CClassConstSP clazz): Model(clazz) {}

IModel::WantsRiskMapping SCID::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** calculate single price and store result in CResult */
void SCID::Price(CInstrument*  instrument, 
                            CControl*     control, 
                            CResults*     results){
    static const string method = "SCID::Price";
    IIntoProduct* intoProd;
    if (!IIntoProduct::TYPE->isInstance(instrument) ||
        !(intoProd = dynamic_cast<IIntoProduct*>(instrument))){
        throw ModelException(method, "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support SCID::IntoProduct");
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

void SCID::getMarket(const MarketData*  market,  IInstrumentCollectionSP instruments) {
	sCIDparam.getData(this, market); // get the sCID parameters
    CModel::getMarket(market,instruments);
}

class SCIDHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(SCID, clazz);
        SUPERCLASS(Model);
        EMPTY_SHELL_METHOD(defaultSCID);

		FIELD(MCmethod,					 "for testing: MC method to use: MARKET, CLOSED FORM, FAST, FULL, FULL1, FULL2, OR HYBRID");
		FIELD_MAKE_OPTIONAL(MCmethod);
		FIELD(freqFastMC,				 "how often do we compute the loss distribution");
		FIELD(freqFullMC,				 "how often do we simulate spread and losses");
		FIELD(seed,						 "seed for the simulation");
		FIELD(timeSteps,                  "time Steps used in the discretization of the CIR process");
 		FIELD(nbPathsFast,                "nbPaths conditional on no jumps and conditional on at least one jump");
  		FIELD(nbPathsFull,                "nbPaths used in the Full MC");
		FIELD(ConvolutionMethod,          "Convolution Method conditional on no jumps and conditional on at least one jump");
 		FIELD(sCIDparam,                  "sCID parameters");
		FIELD(test,						 "for testing, to remove");
		FIELD_MAKE_OPTIONAL(test);
    }




    static IObject* defaultSCID() {
        return new SCID();
    }
};

CClassConstSP const SCID::TYPE = CClass::registerClassLoadMethod(
    "SCID", typeid(SCID), SCIDHelper::load);

// for linker
bool   SCIDLoad() {
    return (SCID::TYPE != 0);
   }

CClassConstSP const SCID::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("SCID::IIntoProduct",
                                    typeid(SCID::IIntoProduct), 0);

DRLIB_END_NAMESPACE
