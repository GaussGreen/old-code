//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ShortTermSqueezeParallel.cpp
//
//   Description : Sensitivity to short term correlation spread
//
//   Author      : Eva X Strasser
//
//   Date        : 15 Aug 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ShortTermSqueezeTweak.hpp"
#include "edginc/SensitivityFactory.hpp"


DRLIB_BEGIN_NAMESPACE

/** Sens Control for ShortTermSqueezeParallel */
class ShortTermSqueezeParallel: public Sensitivity,
                                virtual public Additive, 
                                public virtual IScenarioShift,
                                public virtual IPerturbation {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    ShortTermSqueezeParallel(double shiftSize): 
        Sensitivity(TYPE), shiftSize(shiftSize) {}

    /** for reflection */
    ShortTermSqueezeParallel(): 
        Sensitivity(TYPE), shiftSize(DEFAULT_SHIFT) {}

    /** identifies the packet in which the results are stored. ShortTermSqueezeParallel
        results are stored in the instrument packet */
    const string& getPacketName() const{
        return Results::INSTRUMENT_PACKET;
    }    

    /** identifies the name used for storing associated results in the output*/
    const string& getSensOutputName() const{
        return NAME;
    }

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    bool discreteShift() const{
        return false;
    }

    /** calculates given sensitivity - invoked by calculateSens */
    void calculate(TweakGroup*      tweakGroup,
                                         Results*         results){
        static const string method("ShortTermSqueezeParallel::calculate");
        try {
            OutputNameConstSP outputName(new OutputName(getSensOutputName()));
            const string& packetName = getPacketName(); // eg INSTRUMENT
            if (!results->exists(packetName, outputName)) {
                // build a ShortTermSqueezeParallel shift
                ShortTermSqueezeTweakSP phi(new ShortTermSqueezeTweak(shiftSize, true));
                if (!phi.get()) {
                    throw ModelException(method, "Internal error");
                }

                // need to see if there is actually anything to tweak
                OutputNameArrayConstSP names(phi->names(tweakGroup));
                
                if (names->empty()) {
                    results->storeNotApplicable(this);
                } else {
                    // store what we want to shift
                    OutputNameSP nullName; // null => all correlations
                    phi->setMarketDataName(nullName);
                    try {
                        // calculate sens
                        double firstDeriv = phi->calculateOneSidedFirstDeriv(
                            tweakGroup->getModel(),
                            tweakGroup->getInstrument(),
                            control, // from Sensitivity base class
                            results);
                        // and store it
                        results->storeScalarGreek(firstDeriv,
                                                  packetName,
                                                  outputName);
                    }
                    catch (exception& e) {
                        results->storeGreek(IObjectSP(new Untweakable(e)),
                                            packetName, outputName);
                    }
                }
            }
        } catch (exception& e){
            throw ModelException(&e,  "ShortTermSqueezeParallel::calculate");
        }
    }

    virtual bool applyScenario(IObjectSP object) {
        ShortTermSqueezeTweakSP phi(new ShortTermSqueezeTweak(shiftSize, true));
        OutputNameConstSP nullName;
        return phi->findAndShift(object, nullName);
    }
    
    // Nothing to do before market data is retrieved.
    virtual bool preapplyScenario(IObjectSP object) {
        return false;
    }

    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name) {
        return applyScenario(objectToShift);
    }

private:
    friend class ShortTermSqueezeParallelHelper;
    // **** registered fields ****
    double shiftSize;
};

typedef smartPtr<ShortTermSqueezeParallel> ShortTermSqueezeParallelSP;

class ShortTermSqueezeParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new ShortTermSqueezeParallel(ShortTermSqueezeParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new ShortTermSqueezeParallel(shiftSize);
        }
    };

public:

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(ShortTermSqueezeParallel, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        IMPLEMENTS(IScenarioShift);
        IMPLEMENTS(IPerturbation);
        EMPTY_SHELL_METHOD(defaultShortTermSqueezeParallel);
        FIELD(shiftSize, "How big to make the tweak");
        // register how to build our sensitivity
        SensitivityFactory::addSens(ShortTermSqueezeParallel::NAME, 
                                    new Factory(), 
                                    new ShortTermSqueezeParallel(),
                                    ShortTermSqueezeTweak::IShift::TYPE);
    }
    
    static IObject* defaultShortTermSqueezeParallel(){
        return new ShortTermSqueezeParallel();
    }
};

const string ShortTermSqueezeParallel::NAME = "SHORT_TERM_SQUEEZE_PARALLEL";
const double ShortTermSqueezeParallel::DEFAULT_SHIFT = 0.01;


CClassConstSP const ShortTermSqueezeParallel::TYPE = CClass::registerClassLoadMethod(
        "ShortTermSqueezeParallel", typeid(ShortTermSqueezeParallel), ShortTermSqueezeParallelHelper::load);

// to force linker to include file (avoid having header file) */
bool ShortTermSqueezeParallelLinkIn() {
    return (ShortTermSqueezeParallel::TYPE != 0);
}

DRLIB_END_NAMESPACE
