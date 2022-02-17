#include "edginc/config.hpp"
#include "edginc/VolRegime.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

void VolRegime::validatePop2Object(){
	const static string method = "VolRegime::validatePop2Object"; 

	try{
		validate();
	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

IObject* VolRegime::defaultCtor(){
    return new VolRegime();
}

/** Invoked when Class is 'loaded' */
void VolRegime::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VolRegime, clazz);
    SUPERCLASS(RegimeFactor);
    EMPTY_SHELL_METHOD(defaultCtor);
}

VolRegime::VolRegime(): 
RegimeFactor(TYPE){}

CClassConstSP const VolRegime::TYPE = CClass::registerClassLoadMethod(
	"VolRegime", typeid(VolRegime), load);

/* external symbol to allow class to be forced to be linked in */
bool VolRegimeLinkIn(){
    return (VolRegime::TYPE != 0);  
}

DRLIB_END_NAMESPACE
