
#include "edginc/config.hpp"
#define QLIB_INSTRUMENT_CPP
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE

IInstrument::~IInstrument(){}
IInstrument::IInstrument(){}

void IInstrument::load(CClassSP& clazz){
    REGISTER_INTERFACE(IInstrument, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IInstrument::TYPE = CClass::registerInterfaceLoadMethod(
    "IInstrument", typeid(IInstrument), load);
    

CInstrument::~CInstrument(){}

CInstrument::CInstrument(CClassConstSP clazz): CObject(clazz){}


/** override a control shift (eg for delta on trees) - may return
    null to use original. Default implementation returns null. */
CSensControl* CInstrument::AlterControl(const IModel*       modelParams,
                                        const CSensControl* sensControl) const{
    // use original
    return 0;
}

/** price a dead instrument until settlement - exercised, expired,
    knocked out etc.  returns true if it is dead (and priced), false
    if it is not dead */
bool CInstrument::priceDeadInstrument(CControl* control, 
                                      CResults* results) const{
    // simple default
    return false;
}


class InstrumentHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CInstrument, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IInstrument);
    }
};

CClassConstSP const CInstrument::TYPE = CClass::registerClassLoadMethod(
    "Instrument", typeid(CInstrument), InstrumentHelper::load);

DEFINE_TEMPLATE_TYPE_WITH_NAME("InstrumentArray", CInstrumentArray);


DRLIB_END_NAMESPACE
