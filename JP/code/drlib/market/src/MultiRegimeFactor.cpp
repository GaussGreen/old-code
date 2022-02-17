#include "edginc/config.hpp"
#define QLIB_MULTIREGIMEFACTOR_CPP
#include "edginc/MultiRegimeFactor.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

int MultiRegimeFactor::getRateIdx(int iState)const {
    return rateIdx[iState];
}

int MultiRegimeFactor::getVolIdx(int iState) const {
    return volIdx[iState];
}

int MultiRegimeFactor::getHazardIdx(int iState)const {
    return hazardIdx[iState];
}

HazardRegimeSP MultiRegimeFactor::getHazardRegime()const {
    return hazard;
}

void MultiRegimeFactor::validatePop2Object(){
	const static string method = "MultiRegimeFactor::validatePop2Object";

	try{
        nbFactors = 0;
        nbStates = 1;
        nbRateStates = 0;
        nbVolStates = 0;
        nbHazardStates = 0;
        if (rate.get()){
            nbFactors++;
            rate->validatePop2Object();
            nbRateStates = rate->nbStates;
            initRateStateIdx = rate->initStateIdx;
            nbStates *= nbRateStates;
        }

        if (vol.get()){
            nbFactors++;
            vol->validatePop2Object();
            nbVolStates = vol->nbStates;
            initVolStateIdx = vol->initStateIdx;
            nbStates *= nbVolStates;
        } 
        
        if (hazard.get()){
            nbFactors++;
            hazard->validatePop2Object();
            nbHazardStates = hazard->nbStates;
            initHazardStateIdx = hazard->initStateIdx;
            nbStates *= nbHazardStates;
        }

        if (nbRateStates>0) {
            rateIdx.resize(nbStates);
            rateStates.resize(nbStates);
        }
        if (nbVolStates>0) {
            volIdx.resize(nbStates);
            volStates.resize(nbStates);
        }
        if (nbHazardStates>0) {
            hazardIdx.resize(nbStates);
            hazardStates.resize(nbStates);
        }

        int iState = 0;
        for (int iRate = 0; iRate < nbRateStates || iRate == 0; iRate++){
            for (int iVol = 0; iVol < nbVolStates || iVol == 0; iVol++){
                for (int iHazard = 0; iHazard < nbHazardStates || iHazard == 0; iHazard++){
                    bool isInitState = true;
                    if (nbRateStates>0) {
                        rateIdx[iState] = iRate;
                        if (iRate != initRateStateIdx) isInitState = false;
                    } 
                    if (nbVolStates>0) {
                        volIdx[iState] = iVol;
                        if (iVol != initVolStateIdx) isInitState = false;
                    }
                    if (nbHazardStates>0) {
                        hazardIdx[iState] = iHazard;
                        if (iHazard != initHazardStateIdx) isInitState = false;
                    }
                    if (isInitState) initState = iState;
                    iState++;
                }
            }
        }

		generator.resize(nbStates, nbStates);
		jumpsOnLogSpot.resize(nbStates, nbStates);
		jumpsAtDrift.resize(nbStates, nbStates);
		transitionProb.resize(nbStates, nbStates);
        transitionProb.resize(nbStates,nbStates);

        update();

	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

void MultiRegimeFactor::update(){
	const static string method = "MultiRegimeFactor::update";

	try{
		int		iRow, iCol;
		for (iRow = 0; iRow < nbStates; ++iRow){
            double sum = 0.;
			for (iCol = 0; iCol < nbStates; ++iCol){
				int nbStateJumps = 0;
                generator[iCol][iRow] = 0.;
                jumpsOnLogSpot[iCol][iRow] = 0.;

                if (nbRateStates>0){
                    if (rateIdx[iCol] != rateIdx[iRow]) {
					    nbStateJumps++;
                        generator[iCol][iRow] = rate->generator[rateIdx[iCol]][rateIdx[iRow]];
                    }
                    jumpsOnLogSpot[iCol][iRow] += rate->jumpsOnLogSpot[rateIdx[iCol]][rateIdx[iRow]];
				}

                if (nbVolStates>0) {
                    if (volIdx[iCol] != volIdx[iRow]) {
					    nbStateJumps++;
                        generator[iCol][iRow] = vol->generator[volIdx[iCol]][volIdx[iRow]];
                    }
                    jumpsOnLogSpot[iCol][iRow] += vol->jumpsOnLogSpot[volIdx[iCol]][volIdx[iRow]];
				}

                if (nbHazardStates>0) {
                    if (hazardIdx[iCol] != hazardIdx[iRow]) {
					    nbStateJumps++;
                        generator[iCol][iRow] = hazard->generator[hazardIdx[iCol]][hazardIdx[iRow]];
                    }
                    jumpsOnLogSpot[iCol][iRow] += hazard->jumpsOnLogSpot[hazardIdx[iCol]][hazardIdx[iRow]];
				}
				if (nbStateJumps > 1){
					generator[iCol][iRow] = 0.;
				}
                sum += generator[iCol][iRow];
				jumpsAtDrift[iCol][iRow] = exp(jumpsOnLogSpot[iCol][iRow])-1.;
			}
            generator[iRow][iRow] = -sum;
			if (nbRateStates>0)     rateStates[iRow] = rate->states[rateIdx[iRow]];
			if (nbVolStates>0)      volStates[iRow] = vol->states[volIdx[iRow]];
			if (nbHazardStates>0)   hazardStates[iRow] = hazard->states[hazardIdx[iRow]];
		}
	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

string MultiRegimeFactor::getName() const {
	return name;
}

void MultiRegimeFactor::getMarket(const IModel* model, const MarketData* market, const string& name) {
    // do nothing
}

void MultiRegimeFactor::calcTransitionProb(double tradYear){
	transitionProb = generator.exp(tradYear, -1.0, 7.0);
}

void MultiRegimeFactor::calcTransitionProb(DateTime date){

}

MultiRegimeFactor::~MultiRegimeFactor(){}

MultiRegimeFactor::MultiRegimeFactor(): 
MarketObject(TYPE){}

MultiRegimeFactor::MultiRegimeFactor(const CClassConstSP& clazz): 
MarketObject(clazz){}

/** Invoked when Class is 'loaded' */
void MultiRegimeFactor::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(MultiRegimeFactor, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IGetMarket);
    IMPLEMENTS(CreditSpreadRhoParallel::RestorableShift);
    IMPLEMENTS(IRestorableWithRespectTo<VolParallel>);
	EMPTY_SHELL_METHOD(defaultCtor);

	FIELD(name, "name");

	FIELD(rate, "rate regime");
	FIELD_MAKE_OPTIONAL(rate);

	FIELD(vol, "vol regime");
	FIELD_MAKE_OPTIONAL(vol);

	FIELD(hazard, "hazard regime");
	FIELD_MAKE_OPTIONAL(hazard);

	FIELD(rateStates, "rate states");
	FIELD_MAKE_TRANSIENT(rateStates);

	FIELD(volStates, "vol states");
	FIELD_MAKE_TRANSIENT(volStates);

	FIELD(hazardStates, "hazard states");
	FIELD_MAKE_TRANSIENT(hazardStates);

	FIELD(generator, "generator");
	FIELD_MAKE_TRANSIENT(generator);

	FIELD(jumpsOnLogSpot, "jumps on log spot when state switches");
	FIELD_MAKE_TRANSIENT(jumpsOnLogSpot);

	FIELD(jumpsAtDrift, "drift adjustment due to jumps");
	FIELD_MAKE_TRANSIENT(jumpsAtDrift);

	FIELD(transitionProb, "transition probability matrix");
	FIELD_MAKE_TRANSIENT(transitionProb);

	FIELD(initRateStateIdx, "the initial rate state index");
	FIELD_MAKE_TRANSIENT(initRateStateIdx);

	FIELD(initVolStateIdx, "the initial vol state index");
	FIELD_MAKE_TRANSIENT(initVolStateIdx);

	FIELD(initHazardStateIdx, "the initial hazard state index");
	FIELD_MAKE_TRANSIENT(initHazardStateIdx);

	FIELD(nbRateStates, "number of rate states");
	FIELD_MAKE_TRANSIENT(nbRateStates);

	FIELD(nbVolStates, "number of vol states");
	FIELD_MAKE_TRANSIENT(nbVolStates);

	FIELD(nbHazardStates, "number of hazard states");
	FIELD_MAKE_TRANSIENT(nbHazardStates);

	FIELD(nbStates, "number of states");
	FIELD_MAKE_TRANSIENT(nbStates);

	FIELD(nbFactors, "number of states");
	FIELD_MAKE_TRANSIENT(nbFactors);

    FIELD(rateIdx, "corresponding rate state index for each state");
	FIELD_MAKE_TRANSIENT(rateIdx);

    FIELD(volIdx, "corresponding vol state index for each state");
	FIELD_MAKE_TRANSIENT(volIdx);

    FIELD(hazardIdx, "corresponding hazard state index for each state");
	FIELD_MAKE_TRANSIENT(hazardIdx);
}

IObject* MultiRegimeFactor::defaultCtor(){
	return new MultiRegimeFactor();
}

/** CreditSpreadRhoParallel sensitivity */
/** Returns name identifying yield curve for rho parallel */
string MultiRegimeFactor::sensName(CreditSpreadRhoParallel* shift) const{
    return getName();  
}

/** Shifts the object using given shift */
bool MultiRegimeFactor::sensShift(CreditSpreadRhoParallel* shift){
   static const string method = "MultiRegimeFactor::sensShift";
   try {
       double shiftSize = shift->getShiftSize();
       if (!Maths::isZero(shiftSize)){
          // tweak non-default states
          for (int i = 0; i < hazard->states.size()-1; i++) {
             hazard->states[i] += shiftSize / (1. - hazard->recoveryRate);
             hazard->generator[hazard->states.size()-1][i] = hazard->states[i];
          }
       }
       update();
       return false; // none of our components has a rho type sensitivity
   }
   catch (exception &e) {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void MultiRegimeFactor::sensRestore(CreditSpreadRhoParallel* shift){
   double shiftSize = shift->getShiftSize();
   if (!Maths::isZero(shiftSize)){
      for (int i = 0; i < hazard->states.size()-1; i++) {
         hazard->states[i] -= shiftSize / (1. - hazard->recoveryRate);
         hazard->generator[hazard->states.size()-1][i] = hazard->states[i];
      }
   }
   update();
}

/** Returns name identifying vol for vega parallel */
string MultiRegimeFactor::sensName(const VolParallel*) const{
    return getName();
}

/** Shifts the object using given shift */
TweakOutcome MultiRegimeFactor::sensShift(const PropertyTweak<VolParallel>& tweak){
    
    static const string routine = "MultiRegimeFactor:sensShift<VolParallel>";
    
    if (!Maths::isZero(tweak.coefficient)){
        for (int i = 0; i < vol->states.size(); i++) {
            vol->states[i] += tweak.coefficient;
            if (Maths::isNegative(vol->states[i])) {
                throw ModelException(routine, "tweaked vol level must be positive.");
            }
        }
    }
    update();
    return TweakOutcome(tweak.coefficient,false); 
}

/** Restores the object to its original form */
void MultiRegimeFactor::sensRestore(const PropertyTweak<VolParallel>& tweak){
    if (!Maths::isZero(tweak.coefficient)){
        for (int i = 0; i < vol->states.size(); i++) {
            vol->states[i] -= tweak.coefficient;
        }
    }
    update();
}

CClassConstSP const MultiRegimeFactor::TYPE = CClass::registerClassLoadMethod(
									"MultiRegimeFactor", typeid(MultiRegimeFactor), load);

DEFINE_TEMPLATE_TYPE(MultiRegimeFactorWrapper);

/* external symbol to allow class to be forced to be linked in */
bool MultiRegimeFactorLinkIn(){
    return (MultiRegimeFactor::TYPE != 0);  
}

DRLIB_END_NAMESPACE
