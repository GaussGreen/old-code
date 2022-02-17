//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MonteCarloTurbo.cpp
//
//   Description : Monte Carlo Algorithm with reducde iterations for greeks
//
//   Date        : 7 October 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MCPathConfig.hpp"

DRLIB_BEGIN_NAMESPACE

class MonteCarloTurbo : public MonteCarlo {
public:
    static CClassConstSP const TYPE;
    friend class MonteCarloTurboHelper;

    /** main control - calculates price and sensitivities. */
    virtual CResultsArraySP RunMulti(IInstrumentCollectionSP instruments, 
                                     CControl* control);

protected:
    MonteCarloTurbo(CClassConstSP clazz);
    // fields
    int nbIterForGreeks;
    
private:
    /* for reflection */
    MonteCarloTurbo();
};

CResultsArraySP MonteCarloTurbo::RunMulti(IInstrumentCollectionSP instruments, 
                                          CControl* control) {
    if (nbIterForGreeks > 0 && nbIterForGreeks < nbIter) {
        // clone model and work on that
        smartPtr<MonteCarloTurbo> modelCopy(copy(this));
        modelCopy->nbIter = nbIterForGreeks;
        // then call parent method - may want to split into blocks.
        // If fewIter specified then greeks done using 20 paths
        CResultsArraySP results(modelCopy->RunMulti(instruments, control));
        // now do price at full iterations. 1st take another copy
        modelCopy = smartPtr<MonteCarloTurbo>(copy(this));
        // then switch caching/quick greeks off
        modelCopy->quickGreeks = "NO"; // should be a constant
        modelCopy->cachedQckGreeks = 0;
        modelCopy->cachePaths = "NO"; // should be a constant
        modelCopy->cachedCachePaths =  false;
        // ensure output requests are calculated [again]
        control->reset(); // would prefer to do myControl.calculate()
        // pass original list of sensitivities down so any LR greeks are done
        instruments->Price(modelCopy.get(), control, results);
        return results;
    } else {
        return MonteCarlo::RunMulti(instruments, control);
    }
}

MonteCarloTurbo::MonteCarloTurbo(CClassConstSP clazz): 
    MonteCarlo(clazz), nbIterForGreeks(0) {}

// for reflection
MonteCarloTurbo::MonteCarloTurbo():MonteCarlo(TYPE), nbIterForGreeks(0) {}

class MonteCarloTurboHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MonteCarloTurbo, clazz);
        SUPERCLASS(MonteCarlo);
        EMPTY_SHELL_METHOD(defaultMonteCarloTurbo);
        FIELD(nbIterForGreeks, "Reduced iterations for greeks");
    }

    static IObject* defaultMonteCarloTurbo(){
        return new MonteCarloTurbo();
    }
};

CClassConstSP const MonteCarloTurbo::TYPE = CClass::registerClassLoadMethod(
    "MonteCarloTurbo", typeid(MonteCarloTurbo), MonteCarloTurboHelper::load);

// for class loading 
bool MonteCarloTurboLoad() {
    return (MonteCarloTurbo::TYPE != 0);
}


/** "standard" MonteCarloTurbo that can be captured in Pyramid using 
     current IMS */
class MonteCarloTurboLNDefault: public MonteCarloTurbo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloTurboLNDefault():MonteCarloTurbo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MonteCarloTurboLNDefault, clazz);
        SUPERCLASS(MonteCarloTurbo);
        EMPTY_SHELL_METHOD(defaultMonteCarloTurboLNDefault);
    }

    static IObject* defaultMonteCarloTurboLNDefault(){
        return new MonteCarloTurboLNDefault();
    }
};

CClassConstSP const MonteCarloTurboLNDefault::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloTurboLNDefault", typeid(MonteCarloTurboLNDefault), load);

class MonteCarloTurboImpliedDefault: public MonteCarloTurbo {
public:
    static CClassConstSP const TYPE;
    
private:
    MonteCarloTurboImpliedDefault():MonteCarloTurbo(TYPE) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MonteCarloTurboImpliedDefault, clazz);
        SUPERCLASS(MonteCarloTurbo);
        EMPTY_SHELL_METHOD(defaultMonteCarloTurboImpliedDefault);
    }

    static IObject* defaultMonteCarloTurboImpliedDefault(){
        return new MonteCarloTurboImpliedDefault();
    }
};

CClassConstSP const MonteCarloTurboImpliedDefault::TYPE = 
CClass::registerClassLoadMethod(
    "MonteCarloTurboImpliedDefault", typeid(MonteCarloTurboImpliedDefault),
    load);

DRLIB_END_NAMESPACE
