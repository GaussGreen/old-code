#include "edginc/config.hpp"
#include "edginc/RegimeFactor.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

void RegimeFactor::validate(){
	const static string method = "RegimeFactor::validate";

	try{
		nbStates = states.size();
		if (!nbStates){
			throw ModelException(method, "Number of states equals to zero.");
		}

		int nbRows = generator.numRows();
		int nbCols = generator.numCols();
		if(nbRows != nbStates || nbCols != nbStates){
			throw ModelException(method, "Number of rows or columns of generator is not equal to nbStates.");
		}

		if(!jumpsOnLogSpot.empty()){
			nbRows = jumpsOnLogSpot.numRows();
			nbCols = jumpsOnLogSpot.numCols();
			if (nbRows != nbStates || nbCols != nbStates) {
				throw ModelException(method, "Number of rows or columns of jumpsOnLogSpot is not equals to nbStates.");
			}
            for (int iRow = 0; iRow < nbRows-1; iRow++){
                for (int iCol = 0; iCol < nbCols; iCol++){
                    if (iRow == iCol && !Maths::isZero(jumpsOnLogSpot[iCol][iRow])){
                        throw ModelException(method, "jump must be zero when state does not switch");
                    }
                }
            }
		} else {
			jumpsOnLogSpot.resize(nbStates,nbStates);
		}

		if (initStateIdx < 0 || initStateIdx >= nbStates){
			throw ModelException(method, "initStateIdx is not between 0 (included) and nbStates (excluded).");
		}

		transitionProb.resize(nbStates,nbStates);
	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

string RegimeFactor::getName() const {
	return name;
}

void RegimeFactor::getMarket(const IModel* model, const MarketData* market, const string& name) {
    // do nothing
}

void RegimeFactor::calcTransitionProb(double tradYear){
	transitionProb = generator.exp(tradYear, -1.0, 7.0);
}

double RegimeFactor::expect(double tradYear){

	if (Maths::isZero(tradYear)){
		return states[initStateIdx];
	}

	calcTransitionProb(tradYear);

	double sum = 0.;
	for (int i=0; i<nbStates; i++){
		sum += transitionProb[i][initStateIdx] * states[i];
	}
	return sum;
}

int RegimeFactor::getInitStateIdx() const {
    return initStateIdx;
}

int RegimeFactor::getNbStates() const {
    return nbStates;
}

DoubleArray RegimeFactor::getStates() const {
    return states;
}

void RegimeFactor::setInitState(double level){
    states[initStateIdx] = level;
}

void RegimeFactor::setStates(DoubleArray levels){
    states = levels;
}

RegimeFactor::~RegimeFactor(){}

RegimeFactor::RegimeFactor(): 
MarketObject(TYPE){}

RegimeFactor::RegimeFactor(const CClassConstSP& clazz): 
MarketObject(clazz){}

/** Invoked when Class is 'loaded' */
void RegimeFactor::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(RegimeFactor, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IGetMarket);
	EMPTY_SHELL_METHOD(defaultCtor);
	FIELD(name, "name");
	FIELD(states, "states");
	FIELD(generator, "generator");
	FIELD(initStateIdx, "the initial state index");
	FIELD(jumpsOnLogSpot, "jumps on log spot when state switches");
	FIELD_MAKE_OPTIONAL(jumpsOnLogSpot);
	FIELD(nbStates, "number of states");
	FIELD_MAKE_TRANSIENT(nbStates);
	FIELD(transitionProb, "transition probability matrix");
	FIELD_MAKE_TRANSIENT(transitionProb);
}

IObject* RegimeFactor::defaultCtor(){
	return new RegimeFactor();
}

CClassConstSP const RegimeFactor::TYPE = CClass::registerClassLoadMethod(
									"RegimeFactor", typeid(RegimeFactor), load);

/* external symbol to allow class to be forced to be linked in */
bool RegimeFactorLinkIn(){
    return (RegimeFactor::TYPE != 0);  
}

DRLIB_END_NAMESPACE
