//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ICreditDefaultSens.cpp
//
//   Description : Sensitivities corresponding to 'shifting' credit default, 
//                 with the possibility to override recovery :
//                 - CreditDefaultSens : use current recovery
//                 - CreditDefaultSensWithZeroRecovery : use recovery = 0
//                 - CreditDefaultSensWithSpecifiedRecovery : use recovery = user specified recovery
//
//   Author      : Antoine Gregoire
//
//   Date        : February 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditDefaultSensBase.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/Results.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/TweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

CreditDefaultSensBase::CreditDefaultSensBase(
        const CClassConstSP& clazz,
        const string& outputName) : SensControlPerName(clazz, outputName) {} 

CreditDefaultSensBase::~CreditDefaultSensBase() {}
CreditDefaultSensBase::IShift::~IShift() {}
CreditDefaultSensBase::IRestorableShift::~IRestorableShift() {}

/** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. */
bool CreditDefaultSensBase::discreteShift() const{
    return true;
}

/** Once used to make a shift, this reports the appropriate divisor
for this sensitivity */
double CreditDefaultSensBase::divisor() const {
    return DIVISOR;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditDefaultSensBase::shiftInterface() const {
    return CreditDefaultSensBase::IShift::TYPE;
}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents which is also
    restorable */
CClassConstSP CreditDefaultSensBase::restorableShiftInterface() const {
    return CreditDefaultSensBase::IRestorableShift::TYPE;
}

/** Returns true if the supplied object matches the supplied name
    for this sensitivity.  The object must implement the
    CreditDefaultSensBase::IShift interface */
bool CreditDefaultSensBase::nameMatches(const OutputName& name,
                                        IObjectConstSP    obj) {
    // cast obj to CreditDefaultSensBase::IShift and then invoke name method
    const IShift& creditDefaultableObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(creditDefaultableObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to
    this sensitivity to the supplied list */
void CreditDefaultSensBase::appendName(OutputNameArray&          namesList,
                                       IObjectConstSP            obj) {
    // cast obj to CreditDefaultSensBase::IShift and then invoke name method
    const IShift& creditDefaultableObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(creditDefaultableObj.sensName(this)));
    namesList.push_back(outputName);
}


/** Shifts the object (which supports being tweaked
    by this type of sens control) using given shift. The return value
    indicates whether or not components of this object need to be
    tweaked ie true: infrastructure should continue to recurse through
    components tweaking them; false: the infrastructure shouldn't
    touch any components within this object */
bool CreditDefaultSensBase::shift(IObjectSP obj) {
    // cast obj to CreditDefaultSensBase::IShift and then invoke shift method
    IShift& creditDefaultableObj = dynamic_cast<IShift&>(*obj);
    return creditDefaultableObj.sensShift(this);
}

/** Restores the object (which supports being tweaked
    by this type of sens control) to its original form */
void CreditDefaultSensBase::restore(IObjectSP obj) {
    // cast obj to CreditDefaultSensBase::IShift and then invoke restore method
    IRestorableShift& restorableObj = dynamic_cast<IRestorableShift&>(*obj);
    restorableObj.sensRestore(this);
}

/** Returns the date of the credit event */
DateTime CreditDefaultSensBase::getCreditEventDate(
    const DateTime& valueDate) const 
{
    return (defaultDate.empty() ? valueDate : defaultDate);
}

/** calculates given sensitivity - invoked by calculateSens */
void CreditDefaultSensBase::calculate(TweakGroup*     tweakGroup,
                                      CResults*       results) {
    try {
        // get list of names to calculate result for. Remove blanks
        OutputNameArrayConstSP names(this->names(tweakGroup, results));

        if (names->empty()) {
            results->storeNotApplicable(this);            
        }

        for (int idx = 0; idx < names->size(); idx++){
            // store what we want to shift
            setMarketDataName((*names)[idx]);
            /* skip over where result has been calculated already */
            if (!results->exists(this)){
                try {
                    // calculate sens
                    double firstDeriv = 
                        calcOneSidedFirstDeriv(tweakGroup, results);
                    // and store it
                    results->storeScalarGreek(firstDeriv, this);
                }
                catch (exception& exc) {
                    results->storeGreek(IObjectSP(new Untweakable(exc)),this);
                }
            }
        }
    } catch (exception& e){
        throw ModelException(&e,  "CreditDefaultSensBase::calculate");
    }
}
    
/** Invoked when this class is 'loaded' */
void CreditDefaultSensBase::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditDefaultSensBase, clazz);
    SUPERCLASS(SensControlPerName);
    IMPLEMENTS(Additive);
    FIELD(defaultDate, "Credit event date. Default: today");
    FIELD_MAKE_OPTIONAL(defaultDate);
}

CClassConstSP const CreditDefaultSensBase::TYPE = 
    CClass::registerClassLoadMethod(
        "CreditDefaultSensBase",
        typeid(CreditDefaultSensBase),
        CreditDefaultSensBase::load);
        
CClassConstSP const CreditDefaultSensBase::IShift::TYPE =
    CClass::registerInterfaceLoadMethod(
        "CreditDefaultSensBase::IShift",
        typeid(CreditDefaultSensBase::IShift),
        0);

static void restorableShiftLoad(CClassSP& clazz){
    REGISTER_INTERFACE(CreditDefaultSensBase::IRestorableShift, clazz);
    EXTENDS(CreditDefaultSensBase::IShift);
}
    
CClassConstSP const CreditDefaultSensBase::IRestorableShift::TYPE = 
    CClass::registerInterfaceLoadMethod(
        "CreditDefaultSensBase::IRestorableShift",
        typeid(CreditDefaultSensBase::IRestorableShift), 
        restorableShiftLoad);
        
const double CreditDefaultSensBase::DIVISOR = 1.0;

// -----------------
// CreditDefaultSens
// -----------------

/** Factory class dictates what methods of building CreditDefaultSens
    sensitivity are supported. Here we support merely a default method. */
class CreditDefaultSensFactory: public virtual SensitivityFactory::IDefault {
public:
    virtual Sensitivity* createDefault();
};

class CreditDefaultSens :
    public CreditDefaultSensBase
{
public:
    static CClassConstSP const TYPE;
    const static string NAME;

    CreditDefaultSens() : CreditDefaultSensBase(TYPE, NAME) {}

    /** Return the overriden recovery given the current recovery */
    virtual double getOverriddenRecovery(double currentRecovery) {
        return currentRecovery;
    }

private:
    static IObject* defaultConstructor() {
        return new CreditDefaultSens();
    }
    
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditDefaultSens, clazz);
        SUPERCLASS(CreditDefaultSensBase);
        EMPTY_SHELL_METHOD(defaultConstructor);
    
        // register how to build our sensitivity
        SensitivityFactory::addSens(
            CreditDefaultSens::NAME,
            new CreditDefaultSensFactory(),
            new CreditDefaultSens(),
            CreditDefaultSensBase::IShift::TYPE);
    }
};

Sensitivity* CreditDefaultSensFactory::createDefault(){
    return new CreditDefaultSens();
}

CClassConstSP const CreditDefaultSens::TYPE = 
    CClass::registerClassLoadMethod(
        "CreditDefaultSens",
        typeid(CreditDefaultSens),
        CreditDefaultSens::load);

const string CreditDefaultSens::NAME = "CREDIT_DEFAULT_SENS";


// ---------------------------------
// CreditDefaultSensWithZeroRecovery
// ---------------------------------

/** Factory class dictates what methods of building CreditDefaultSens
    sensitivity are supported. Here we support merely a default method. */
class CreditDefaultSensWithZeroRecoveryFactory: public virtual SensitivityFactory::IDefault {
public:
    virtual Sensitivity* createDefault();
};

class CreditDefaultSensWithZeroRecovery :
    public CreditDefaultSensBase
{
public:
    static CClassConstSP const TYPE;
    const static string NAME;

    CreditDefaultSensWithZeroRecovery() : CreditDefaultSensBase(TYPE, NAME) {}

    /** Return the overriden recovery given the current recovery */
    virtual double getOverriddenRecovery(double currentRecovery) {
        return 0;
    }

private:
    static IObject* defaultConstructor() {
        return new CreditDefaultSensWithZeroRecovery();
    }
    
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditDefaultSensWithZeroRecovery, clazz);
        SUPERCLASS(CreditDefaultSensBase);
        EMPTY_SHELL_METHOD(defaultConstructor);
    
        // register how to build our sensitivity
        SensitivityFactory::addSens(
            CreditDefaultSensWithZeroRecovery::NAME,
            new CreditDefaultSensWithZeroRecoveryFactory(),
            new CreditDefaultSensWithZeroRecovery(),
            CreditDefaultSensBase::IShift::TYPE);
    }
};

Sensitivity* CreditDefaultSensWithZeroRecoveryFactory::createDefault(){
    return new CreditDefaultSensWithZeroRecovery();
}

CClassConstSP const CreditDefaultSensWithZeroRecovery::TYPE = 
    CClass::registerClassLoadMethod(
        "CreditDefaultSensWithZeroRecovery",
        typeid(CreditDefaultSensWithZeroRecovery),
        CreditDefaultSensWithZeroRecovery::load);

const string CreditDefaultSensWithZeroRecovery::NAME =
    "CREDIT_DEFAULT_SENS_WITH_ZERO_RECOVERY";


// --------------------------------------
// CreditDefaultSensWithSpecifiedRecovery
// --------------------------------------

/** Factory class dictates what methods of building CreditDefaultSens
    sensitivity are supported. Here we support merely a default method. */
class CreditDefaultSensWithSpecifiedRecoveryFactory :
    public SensitivityFactory::IScalar
{
public:
    virtual Sensitivity* createScalar(double specifiedRecovery);
};

class CreditDefaultSensWithSpecifiedRecovery :
    public CreditDefaultSensBase
{
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SPECIFIED_RECOVERY;

    CreditDefaultSensWithSpecifiedRecovery(double specifiedRecovery) :
        CreditDefaultSensBase(TYPE, NAME),
        specifiedRecovery(specifiedRecovery) {}

    /** Return the overriden recovery given the current recovery */
    virtual double getOverriddenRecovery(double currentRecovery) {
        return specifiedRecovery;
    }

private:
    static IObject* defaultConstructor() {
        return new CreditDefaultSensWithSpecifiedRecovery(
            DEFAULT_SPECIFIED_RECOVERY);
    }
    
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CreditDefaultSensWithSpecifiedRecovery, clazz);
        SUPERCLASS(CreditDefaultSensBase);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(specifiedRecovery,
            "The recovery used when simulating the default");
    
        // register how to build our sensitivity
        SensitivityFactory::addSens(
            CreditDefaultSensWithSpecifiedRecovery::NAME,
            new CreditDefaultSensWithSpecifiedRecoveryFactory(),
            new CreditDefaultSensWithSpecifiedRecovery(
                DEFAULT_SPECIFIED_RECOVERY),
            CreditDefaultSensBase::IShift::TYPE);
    }
    
    /** The recovery used when simulating the default */
    double specifiedRecovery;
};

Sensitivity* CreditDefaultSensWithSpecifiedRecoveryFactory::createScalar(
    double specifiedRecovery)
{
    return new CreditDefaultSensWithSpecifiedRecovery(specifiedRecovery);
}

CClassConstSP const CreditDefaultSensWithSpecifiedRecovery::TYPE = 
    CClass::registerClassLoadMethod(
        "CreditDefaultSensWithSpecifiedRecovery",
        typeid(CreditDefaultSensWithSpecifiedRecovery),
        CreditDefaultSensWithSpecifiedRecovery::load);

const string CreditDefaultSensWithSpecifiedRecovery::NAME =
    "CREDIT_DEFAULT_SENS_WITH_GIVEN_RECOVERY";

const double CreditDefaultSensWithSpecifiedRecovery::
    DEFAULT_SPECIFIED_RECOVERY = 0.0;

DRLIB_END_NAMESPACE
