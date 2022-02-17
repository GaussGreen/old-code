//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : AdjCreditSpreadRhoPointwise.cpp
//
//   Description : Adjusted credit spread rho pointwise sensitivity
//
//   Author      : Milan Kovacevic
//
//   Date        : 03 January 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/AdjCreditSpreadRhoPointwise.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/ExpiryResult.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Holiday.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<AdjCreditSpreadRhoPointwise::IRestorableShift> AdjCreditSpreadRhoPointwiseRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
AdjCreditSpreadRhoPointwise::IShift::~IShift(){} // empty
AdjCreditSpreadRhoPointwise::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for nbrho parallel */
const string AdjCreditSpreadRhoPointwise::NAME = "ADJ_CREDIT_SPREAD_RHO_POINTWISE";
const double AdjCreditSpreadRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit IShift size */
AdjCreditSpreadRhoPointwise::AdjCreditSpreadRhoPointwise(double shiftSize):
    VectorShift(TYPE, NAME, shiftSize), origRiskyGrowth(false){
}

/** for reflection */
AdjCreditSpreadRhoPointwise::AdjCreditSpreadRhoPointwise(): 
    VectorShift(TYPE, NAME), origRiskyGrowth(false){
}

/** The supplied object is queried for the expiries array needed
    for doing rho pointwise and this array is returned. The supplied
    object must implement the CreditSpreadRhoPointwise.Shift interface */
IObjectConstSP AdjCreditSpreadRhoPointwise::qualifier(IObjectConstSP obj){
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke qualifier method
    const IShift& rhoObj = 
        dynamic_cast<const IShift&>(*obj);
    return rhoObj.sensExpiries(this);
}


void AdjCreditSpreadRhoPointwise::calculate(TweakGroup*  tweakGroup,
                                            CResults*    results){
    try {
        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } else {
        
            // see if the instrument has a last sens date method
            LastSensDate* lsd = 
                dynamic_cast<LastSensDate*>(tweakGroup->getInstrument());
            DateTime      endDate;
            DateTime      valueDate;

            if (lsd) {
                valueDate = tweakGroup->getInstrument()->getValueDate();
                endDate   = lsd->endDate(this);
            }

            double origPrice;
            for (int idx = 0; idx < names->size(); idx++){
                // store what we want to IShift
                setMarketDataName((*names)[idx]);

                // We need to compute a base price for pricing with risky equity before tweaking
                // We achieve this by shifting and pricing with 0 shiftSize
                double restoreShiftSize = getShiftSize();
                setShiftSize(0.0);           
                
                try {
                    double dummy = 0.0;     // Not used
                    origPrice = shiftAndPrice(tweakGroup, dummy, true);
                }
                catch (exception& e){
                    setShiftSize(restoreShiftSize);
                    throw ModelException(e);
                }
                setShiftSize(restoreShiftSize);
                
                /* skip over blank names or those where result has
                   been calculated already */
                if (!results->exists(this)){
                    try {
                        // now get expiries for which we should tweak
                        cachedExpiries = getExpiries(tweakGroup);
                        // create room for storing the results
                        ExpiryResultArraySP tweaks(new ExpiryResultArray(
                            cachedExpiries->size()));
                        // then loop over the expiries
                        bool expired = false;
                        for (int jdx = 0; jdx < cachedExpiries->size(); jdx++){
                            // to do: ask the instrument when to stop tweaking
                            // store the expiry which we want to tweak
                            expiry = (*cachedExpiries)[jdx];
                            // calculate sens (if not expired)
                            double firstDeriv = 0.0;
                            if (!expired) {
                                double shiftedPrice = shiftAndPrice(tweakGroup, origPrice);
                                double divisor = this->divisor();
                                firstDeriv = (shiftedPrice - origPrice)/divisor;
                            }                           
                            // store result in array
                            (*tweaks)[jdx] = ExpiryResult(expiry, firstDeriv);

                            // do we need to tweak anymore ?
                            if (lsd) {
                                expired = expiry->toDate(valueDate).isGreater(endDate);
                            }
                        }
                        // and store it
                        results->storeGreek(tweaks, this);
                    }
                    catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)), this);
                    }
                }
            }
        }

    } catch (exception& e){
        results->storeGreek(IObjectSP(new Untweakable(e)), getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}


void AdjCreditSpreadRhoPointwise::setRiskyGrowth(bool riskyGrowth) {
    origRiskyGrowth = riskyGrowth;
}

bool AdjCreditSpreadRhoPointwise::getRiskyGrowth() const {
    return origRiskyGrowth;
}

/** Once used to make a IShift, this reports the appropriate divisor
    for this sensitivity */
double AdjCreditSpreadRhoPointwise::divisor() const{
    static const string method = "AdjCreditSpreadRhoPointwise::divisor";
    double shiftSize;
    try{
        // just scale the IShift size for a 1bp move
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(method, "IShift size is zero");
        }
        return (shiftSize/ONE_BASIS_POINT);
    } 
    catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP AdjCreditSpreadRhoPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP AdjCreditSpreadRhoPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    IShift interface */
bool AdjCreditSpreadRhoPointwise::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to AdjCreditSpreadRhoPointwise::IShift and then invoke name method
    const IShift& RhoPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    return name.equals(RhoPointwiseObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    IShift interface */
void AdjCreditSpreadRhoPointwise::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to AdjCreditSpreadRhoPointwise::IShift and then invoke name method
    const IShift& RhoPointwiseObj =
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(RhoPointwiseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool AdjCreditSpreadRhoPointwise::shift(IObjectSP obj) {
    // cast obj to AdjCreditSpreadRhoPointwise::IShift and then invoke IShift method
    IShift& RhoPointwiseObj =
        dynamic_cast<IShift&>(*obj);
    return RhoPointwiseObj.sensShift(this);
}

void AdjCreditSpreadRhoPointwise::restore(IObjectSP obj) {
    // cast obj to AdjCreditSpreadRhoPointwise::IShift and then invoke restore method
    IRestorableShift& RhoPointwiseObj =
        dynamic_cast<IRestorableShift&>(*obj);
    RhoPointwiseObj.sensRestore(this);
}

class AdjCreditSpreadRhoPointwiseHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given IShift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new AdjCreditSpreadRhoPointwise(AdjCreditSpreadRhoPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new AdjCreditSpreadRhoPointwise(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AdjCreditSpreadRhoPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultAdjCreditSpreadRhoPointwise);
        FIELD_NO_DESC(origRiskyGrowth);
        FIELD_MAKE_TRANSIENT(origRiskyGrowth);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(AdjCreditSpreadRhoPointwise::NAME, 
                                    new Factory(), 
                                    new AdjCreditSpreadRhoPointwise(AdjCreditSpreadRhoPointwise::DEFAULT_SHIFT),
                                    AdjCreditSpreadRhoPointwise::IShift::TYPE);
    }

    static IObject* defaultAdjCreditSpreadRhoPointwise(){
        return new AdjCreditSpreadRhoPointwise();
    }
};

CClassConstSP const AdjCreditSpreadRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "AdjCreditSpreadRhoPointwise", typeid(AdjCreditSpreadRhoPointwise), AdjCreditSpreadRhoPointwiseHelper::load);

CClassConstSP const AdjCreditSpreadRhoPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "AdjCreditSpreadRhoPointwise::IShift", typeid(AdjCreditSpreadRhoPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(AdjCreditSpreadRhoPointwise::IRestorableShift, clazz);
    EXTENDS(AdjCreditSpreadRhoPointwise::IShift);
}
    
CClassConstSP const AdjCreditSpreadRhoPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "AdjCreditSpreadRhoPointwise::IRestorableShift", typeid(AdjCreditSpreadRhoPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE

