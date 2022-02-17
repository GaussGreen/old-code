/**
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/CurrencyBasisSpreadLevelTweak.hpp"
#include "edginc/ScalarPerturbation.hpp"
#include "edginc/RestorableWith.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Sets a Currency Basis spread curve to the level indicated as "shift size".
 *
 * Despite the name and the interface, this is not a greek but a scenario to
 * set the level of the par spread curve of the Currency Basis. The spread 
 * curve will be set flat to the value specified as "shift size".
 *
 */
class CurrencyBasisSpreadLevel : public ScalarPerturbation,
                                 public virtual CurrencyBasisSpreadLevelTweak {
    typedef TweakableWith<CurrencyBasisSpreadLevelTweak> Tweakable;

public:
    /** Reflection type of this class */
    static CClassConstSP const TYPE;

    /** A CurrencyBasisSpreadLevel with given shift size */
    CurrencyBasisSpreadLevel(double level): ScalarPerturbation(TYPE, level) {}

    /** The interface which the objects to be tweaked must implement */
    CClassConstSP shiftInterface() const {
        return Tweakable::TYPE;
    }

    /** Whether @a obj's name is @a name */
    bool nameMatches(const OutputName& name, IObjectConstSP obj) {
        return name.equals(
            dynamic_cast<const Tweakable &>(*obj).sensName(this));
    }

    /** Appends the name(s) of the supplied object with respect to
     *  this sensitivity to the supplied list 
     */
    void appendName(OutputNameArray& namesList, IObjectConstSP obj) {
        namesList.push_back(OutputNameSP(new OutputName(
            dynamic_cast<const Tweakable &>(*obj).sensName(this))));
    }

    /** Shifts the object (which supports being tweaked
     * by this type of sens control) using given shift. The return value
     * indicates whether or not components of this object need to be
     * tweaked ie: true) infrastructure should continue to recurse through
     * components tweaking them; false) the infrastructure shouldn't
     * touch any components within this object 
     */
    bool shift(IObjectSP obj) {
        return dynamic_cast<Tweakable &>(*obj).sensShift(this);
    }

private:
    CurrencyBasisSpreadLevel(): ScalarPerturbation(TYPE) {}

    static IObject* defaultCurrencyBasisSpreadLevel() {
        return new CurrencyBasisSpreadLevel();
    }

    /**
     * Invoked when Class is 'loaded'
     */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CurrencyBasisSpreadLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultCurrencyBasisSpreadLevel);
    }
};

CClassConstSP const CurrencyBasisSpreadLevel::TYPE =
CClass::registerClassLoadMethod("CurrencyBasisSpreadLevel",
                                typeid(CurrencyBasisSpreadLevel),
                                load);


/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the  exe.
 */
bool CurrencyBasisSpreadLevelLinkIn() {
    return CurrencyBasisSpreadLevel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
