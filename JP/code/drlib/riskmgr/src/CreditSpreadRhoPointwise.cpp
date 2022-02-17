//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadRhoPointwise.cpp
//
//   Description : Controls calculation of spread Rho pointwise
//
//   Author      : André Segger
//
//   Date        : 16 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditSpreadRhoPointwise.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/ExpiryResult.hpp"

DRLIB_BEGIN_NAMESPACE
typedef smartPtr<CreditSpreadRhoPointwise::IRestorableShift> 
CreditSpreadRhoPointwiseRestorableShiftSP;

// not sure why we have to have these - MSVC gets upset if they're not in the
// header
CreditSpreadRhoPointwise::IShift::~IShift(){} // empty
CreditSpreadRhoPointwise::IRestorableShift::~IRestorableShift(){} // empty

const double ONE_BASIS_POINT = 0.0001;

/** Sens Control for rho pointwise */
const string CreditSpreadRhoPointwise::NAME = "CREDIT_SPREAD_RHO_POINTWISE";
const double CreditSpreadRhoPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;

/** constructor with explicit shift size */
CreditSpreadRhoPointwise::CreditSpreadRhoPointwise(double     shiftSize):
    VectorShift(TYPE, NAME, shiftSize){
}

/** for reflection */
CreditSpreadRhoPointwise::CreditSpreadRhoPointwise(): 
    VectorShift(TYPE, NAME){
}

/** Once used to make a shift, this reports the appropriate divisor
    for this sensitivity */
double CreditSpreadRhoPointwise::divisor() const{
    static const string method = "CreditSpreadRhoPointwise::divisor";
    double shiftSize;
    try{
        // just scale the shift size for a 1bp move
        shiftSize = getShiftSize();
        if (Maths::isZero(shiftSize)){
            throw ModelException(method, "Shift size is zero");
        }
        return (shiftSize/ONE_BASIS_POINT);
    } 
    catch (ModelException& e) {
        throw ModelException(&e, method);
    }
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditSpreadRhoPointwise::shiftInterface() const{
    return IShift::TYPE;
}
 
/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CreditSpreadRhoPointwise::restorableShiftInterface() const{
    return IRestorableShift::TYPE;
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement the
    Shift interface */
bool CreditSpreadRhoPointwise::nameMatches(const OutputName&         name,
                               IObjectConstSP          obj){
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke name method
    const IShift& rhoPointwiseObj = 
        dynamic_cast<const IShift&>(*obj);
    return name.equals(rhoPointwiseObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement the
    Shift interface */
void CreditSpreadRhoPointwise::appendName(OutputNameArray&          namesList,
                              IObjectConstSP          obj){
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke name method
    const IShift& rhoPointwiseObj = 
        dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(rhoPointwiseObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditSpreadRhoPointwise::shift(IObjectSP obj) {
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke shift method
    IShift& rhoObj = 
        dynamic_cast<IShift&>(*obj);
    return rhoObj.sensShift(this);
}

void CreditSpreadRhoPointwise::restore(IObjectSP obj) {
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke restore method
    IRestorableShift& rhoObj = 
        dynamic_cast<IRestorableShift&>(*obj);
    rhoObj.sensRestore(this);
}

/** The supplied object is queried for the expiries array needed
    for doing rho pointwise and this array is returned. The supplied
    object must implement the CreditSpreadRhoPointwise.Shift interface */
IObjectConstSP CreditSpreadRhoPointwise::qualifier(IObjectConstSP obj){
    // cast obj to CreditSpreadRhoPointwise::Shift and then invoke qualifier method
    const IShift& rhoObj = 
        dynamic_cast<const IShift&>(*obj);
    return rhoObj.sensExpiries(this);
}


/** implements a one sided vector derivative for each instance of the
    market data which is sensitive to this SensControl */
void CreditSpreadRhoPointwise::calculate(TweakGroup*  tweakGroup,
                                         CResults*    results)
{
    try{
        // get list of names to calculate result for
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        } 
        else {

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
                // store the name of what we want to shift
                setMarketDataName((*names)[idx]);

                // We need to compute a base price for pricing off credit spreads before tweaking
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

                /* skip over where result has been calculated already */
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

class CreditSpreadRhoPointwiseHelper{
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new CreditSpreadRhoPointwise(CreditSpreadRhoPointwise::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new CreditSpreadRhoPointwise(shiftSize);
        }
    };
        
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditSpreadRhoPointwise, clazz);
        SUPERCLASS(VectorShift);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultCreditSpreadRhoPointwise);
        // no fields
        // register how to build our sensitivity
        SensitivityFactory::addSens(CreditSpreadRhoPointwise::NAME, 
                                    new Factory(), 
                                    new CreditSpreadRhoPointwise(CreditSpreadRhoPointwise::DEFAULT_SHIFT),
                                    CreditSpreadRhoPointwise::IShift::TYPE);
    }

    static IObject* defaultCreditSpreadRhoPointwise(){
        return new CreditSpreadRhoPointwise();
    }
};

CClassConstSP const CreditSpreadRhoPointwise::TYPE = CClass::registerClassLoadMethod(
    "CreditSpreadRhoPointwise", typeid(CreditSpreadRhoPointwise), CreditSpreadRhoPointwiseHelper::load);

CClassConstSP const CreditSpreadRhoPointwise::IShift::TYPE =
CClass::registerInterfaceLoadMethod(
    "CreditSpreadRhoPointwise::IShift", typeid(CreditSpreadRhoPointwise::IShift), 0);


static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CreditSpreadRhoPointwise::IRestorableShift, clazz);
    EXTENDS(CreditSpreadRhoPointwise::IShift);
}
    
CClassConstSP const CreditSpreadRhoPointwise::IRestorableShift::TYPE = 
CClass::registerInterfaceLoadMethod(
    "CreditSpreadRhoPointwise::IRestorableShift",
    typeid(CreditSpreadRhoPointwise::IRestorableShift), 
    restorableShiftLoad);

DRLIB_END_NAMESPACE
