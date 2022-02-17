/**
 * @file PerturbationTest.cpp
 */

#include "edginc/config.hpp"
#include "edginc/ScenarioShift.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Very simple test harness for low-level testing of the effect of
 * IPerturbation or IScenarioShift implementations
 *
 * You pass in any old object and a scenario, it returns the object as modified
 * by the scenario.  See for instance
 * testing/flexiblegreeksinp/operator-auto-sig-0_01.xml.
 */

class PerturbationTest: public CObject,
                        public virtual ClientRunnable {
    IScenarioShiftSP scenario;
    IObjectSP object;

    static IObject* defaultOne() {
        return new PerturbationTest();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to XML
        REGISTER(PerturbationTest, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultOne);
        FIELD(scenario, "scenario");
        FIELD(object, "object");
        FIELD_MAKE_OPTIONAL(object); // for NULL testing
    }

public:

    static const CClassConstSP TYPE;

    PerturbationTest(): CObject(TYPE)
    {}

    ~PerturbationTest() {}

    IObjectSP run() {
        try {
            scenario->applyScenario(object);
            return object;
        }
        catch (exception& e) {
            throw ModelException(e, "PerturbationTest::PerturbationTest");
        }
    }
};

CClassConstSP const PerturbationTest::TYPE = CClass::registerClassLoadMethod(
    "PerturbationTest", typeid(PerturbationTest), load);

bool PerturbationTestLinkIn() {
    return PerturbationTest::TYPE != NULL;
}

DRLIB_END_NAMESPACE
