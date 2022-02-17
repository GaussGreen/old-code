//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids 
//
//   Description : CDO portfolio name substitution scenario 
//
//   Author      : Linus Thand 
//
//   Date        : 24 May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE

/** A scenario which replaces a name in a CDO portfolio with another, assuming
 * that the name being substituted in exists in the market data. If the name
 * that should be substituted is not found nothing will happen. This is 
 * used by managed trades where clients can request a substitution in the porfolio,
 * and where another variable can be solved for for PV neutrality (e.g a parallel
 * shift in strikes). **/

class RISKMGR_DLL MarketDataNameSubstitution : public CObject,
                                               public virtual IScenarioShift {
  public:
    static CClassConstSP const TYPE;
    virtual bool applyScenario(IObjectSP object);
    virtual bool preapplyScenario(IObjectSP object);
  private:
    static void load (CClassSP& clazz);
    static IObject *defaultConstructor();
    ~MarketDataNameSubstitution();
    MarketDataNameSubstitution();
    MarketDataNameSubstitution(const MarketDataNameSubstitution& rhs); 
    MarketDataNameSubstitution& operator=(const MarketDataNameSubstitution& rhs);

    /* Fields */
    string marketDataType;
    string oldName;
    string newName;
};

FORWARD_DECLARE(MarketDataNameSubstitution);

MarketDataNameSubstitution::~MarketDataNameSubstitution() 
{}

MarketDataNameSubstitution::MarketDataNameSubstitution() : CObject (TYPE) 
{} /// For reflection

/** Apply this scenario shift to the supplied object. The return
    value indicates if anything was actually shifted (true => yes) */
bool MarketDataNameSubstitution::preapplyScenario (IObjectSP object) {
    static const string method = "MarketDataNameSubstitution::preapplyScenario";
    //// Substitute name
    try {
        //// Create action object for iteration
        class Action: public ObjectIteration::IAction {
         public:
            string oldName, newName;
            CClassConstSP baseClass;
            bool shifted;
            Action(const string& marketDataType, 
                   const string& oldName, 
                   const string& newName) : 
                   oldName(oldName), 
                   newName(newName),
                   shifted(false), 
                   baseClass(CClass::forName(marketDataType)) {}

            bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {
                MarketObjectWrapperConstSP candidate = 
                    MarketObjectWrapperConstSP::dynamicCast(obj);
                if (baseClass->isAssignableFrom(candidate->getMOType()) &&
                    (candidate->getName() == oldName) )
                {
                    //// If the object is of the right type, check name and try to substitute
                    MarketObjectWrapperSP currentMOW = 
                        MarketObjectWrapperSP::dynamicCast(state.getObject());
                    shifted |= currentMOW->substituteName(oldName, newName);
                }
                return true;
            }
        };
        Action action(marketDataType, oldName, newName);

        //// Iterate over objects to find name to substitute
        ObjectIteration iteration(MarketObjectWrapper::TYPE);
        iteration.recurse(action, object);

        //// Return true if a name has been substituted       
        return action.shifted;
    } catch (exception & e) {
        throw ModelException (e, method);
    }
}

bool MarketDataNameSubstitution::applyScenario (IObjectSP object) {
    //// Do nothing after market data is fetched
    return false;
}

void MarketDataNameSubstitution::load (CClassSP & clazz) {
    clazz->setPublic ();    // make visible to EAS/spreadsheet
    REGISTER(MarketDataNameSubstitution, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IScenarioShift);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(marketDataType, "market data type to change");
    FIELD(oldName, "name to change from");
    FIELD(newName, "name to change to");
} 

IObject* MarketDataNameSubstitution::defaultConstructor () {
    return new MarketDataNameSubstitution ();
}

CClassConstSP const
MarketDataNameSubstitution::TYPE =
    CClass::registerClassLoadMethod ("MarketDataNameSubstitution",
                                     typeid(MarketDataNameSubstitution),
                                     load);

bool MarketDataNameSubstitutionLinkIn () {
    return MarketDataNameSubstitution::TYPE != NULL;
}

DRLIB_END_NAMESPACE
