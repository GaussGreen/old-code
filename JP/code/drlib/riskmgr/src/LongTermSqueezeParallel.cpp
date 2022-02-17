//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : LongTermSqueezeParallel.cpp
//
//   Description : Sensitivity to short term correlation spread
//
//   Author      : Eva X Strasser
//
//   Date        : 15 Aug 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/LongTermSqueezeTweak.hpp"
#include "edginc/SensitivityFactory.hpp"


DRLIB_BEGIN_NAMESPACE

/** Sens Control for LongTermSqueezeParallel */
class LongTermSqueezeParallel: public Sensitivity,
                               virtual public Additive,
                               public virtual IScenarioShift,
                               public virtual IPerturbation { 
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    LongTermSqueezeParallel(double shiftSize): 
        Sensitivity(TYPE), shiftSize(shiftSize) {}

    /** for reflection */
    LongTermSqueezeParallel(): 
        Sensitivity(TYPE), shiftSize(DEFAULT_SHIFT) {}

    /** identifies the packet in which the results are stored. LongTermSqueezeParallel
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
        static const string method("LongTermSqueezeParallel::calculate");
        try {
            OutputNameConstSP outputName(new OutputName(getSensOutputName()));
            const string& packetName = getPacketName(); // eg INSTRUMENT
            if (!results->exists(packetName, outputName)) {
                // build a LongTermTweakParallel shift
                LongTermSqueezeTweakSP phi(new LongTermSqueezeTweak(shiftSize, true));
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
            throw ModelException(&e,  "LongTermSqueezeParallel::calculate");
        }
    }

    // Nothing to do before market data is retrieved.
    virtual bool preapplyScenario(IObjectSP object){
        return false;
    }

    virtual bool applyScenario(IObjectSP object) {
        LongTermSqueezeTweakSP phi(new LongTermSqueezeTweak(shiftSize, true));
            OutputNameConstSP nullName;
            return phi->findAndShift(object, nullName);
    }
    
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name) {
        return applyScenario(objectToShift);
    }

private:
    friend class LongTermSqueezeParallelHelper;
    // **** registered fields ****
    double shiftSize;
};

typedef smartPtr<LongTermSqueezeParallel> LongTermSqueezeParallelSP;

class LongTermSqueezeParallelHelper {
    /** Factory class dictates what methods of building this
        sensitivity are supported. Here we support a default + one for
        a given shift size */
    class Factory: public SensitivityFactory::IDefault,
                   public SensitivityFactory::IScalar {
    public:
        virtual Sensitivity* createDefault(){
            return new LongTermSqueezeParallel(LongTermSqueezeParallel::DEFAULT_SHIFT);
        }
        virtual Sensitivity* createScalar(double shiftSize){
            return new LongTermSqueezeParallel(shiftSize);
        }
    };

public:

    /** Invoked when SmallClass is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(LongTermSqueezeParallel, clazz);
        SUPERCLASS(Sensitivity);
        IMPLEMENTS(Additive);
        IMPLEMENTS(IScenarioShift);
        IMPLEMENTS(IPerturbation);
        EMPTY_SHELL_METHOD(defaultLongTermSqueezeParallel);
        FIELD(shiftSize, "How big to make the tweak");
        // register how to build our sensitivity
        SensitivityFactory::addSens(LongTermSqueezeParallel::NAME, 
                                    new Factory(), 
                                    new LongTermSqueezeParallel(),
                                    LongTermSqueezeTweak::IShift::TYPE);
    }
    
    static IObject* defaultLongTermSqueezeParallel(){
        return new LongTermSqueezeParallel();
    }
};

const string LongTermSqueezeParallel::NAME = "LONG_TERM_SQUEEZE_PARALLEL";
const double LongTermSqueezeParallel::DEFAULT_SHIFT = 0.01;


CClassConstSP const LongTermSqueezeParallel::TYPE = CClass::registerClassLoadMethod(
        "LongTermSqueezeParallel", typeid(LongTermSqueezeParallel), LongTermSqueezeParallelHelper::load);

// to force linker to include file (avoid having header file) */
bool LongTermSqueezeParallelLinkIn() {
    return (LongTermSqueezeParallel::TYPE != 0);
}

DRLIB_END_NAMESPACE
