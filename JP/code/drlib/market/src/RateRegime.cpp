#include "edginc/config.hpp"
#include "edginc/RateRegime.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

void RateRegime::validatePop2Object(){
	const static string method = "RateRegime::validatePop2Object";

	try{
		validate();
	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

double	RateRegime::calcRisklessBond(double tradYear) const {
	DoubleMatrix	tmpMatrix(generator);
	DoubleArray		tmpArray(nbStates);
	double			bond = 0.;
	int i = 0;
	for (; i < nbStates; ++i){
		tmpMatrix[i][i] -= states[i];
	}
	tmpMatrix = tmpMatrix.exp(tradYear, -1.0, 7.0);
	tmpMatrix.transpose();
	tmpArray[initStateIdx] = 1.;
	tmpArray = tmpMatrix.mult(tmpArray);
	for (i = 0; i < nbStates; ++i){
		bond += tmpArray[i];
	}
	return bond;
}

double RateRegime::calcYieldCurve(double tradYear) const {
	double bond = calcRisklessBond(tradYear);
	return - log(bond)/tradYear;
}

DoubleArray RateRegime::calcYieldCurve(DoubleArray tradYears) const {
	DoubleArray		yieldCurve(tradYears.size());
	int i = 0;
	for (; i < tradYears.size(); ++i){
		yieldCurve[i] = -log(calcRisklessBond(tradYears[i]))/tradYears[i];
	}
	return yieldCurve;
}

IObject* RateRegime::defaultCtor(){
    return new RateRegime();
}

/** Invoked when Class is 'loaded' */
void RateRegime::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RateRegime, clazz);
    SUPERCLASS(RegimeFactor);
    EMPTY_SHELL_METHOD(defaultCtor);
}

RateRegime::RateRegime(): 
RegimeFactor(TYPE){}

class RateRegimeAddin: public CObject{
    static CClassConstSP const TYPE;

    RateRegimeSP		regime;
	DoubleArray			tradYears;

	static IObjectSP calc(RateRegimeAddin* params) {
        return params->calcYieldCurve();
    }

    IObjectSP calcYieldCurve(){
		DoubleArraySP output(new DoubleArray(tradYears.size()));
        *output = regime->calcYieldCurve(tradYears);
        return IObjectSP(output);
    }

    RateRegimeAddin(): CObject(TYPE) {}

    static void load(CClassSP& clazz){
        REGISTER(RateRegimeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultRateRegimeAddin);

        FIELD(regime, "");
		FIELD(tradYears, "");

        Addin::registerClassObjectMethod("CALC_YIELD_CURVE",
                                          Addin::MARKET,
                                          "returns yield curve",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)calc);
    }

    static IObject* defaultRateRegimeAddin(){
        return new RateRegimeAddin();
    }
};

CClassConstSP const RateRegimeAddin::TYPE = CClass::registerClassLoadMethod(
	"RateRegimeAddin", typeid(RateRegimeAddin), RateRegimeAddin::load);

CClassConstSP const RateRegime::TYPE = CClass::registerClassLoadMethod(
	"RateRegime", typeid(RateRegime), load);

/* external symbol to allow class to be forced to be linked in */
bool RateRegimeLinkIn(){
    return (RateRegime::TYPE != 0);  
}

DRLIB_END_NAMESPACE
