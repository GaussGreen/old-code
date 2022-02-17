//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : random.cpp
//
//   Description : Provide generation and optional storage of random numbers
//                 Both uniform, normal and correlated normal
//
//   Date        : May 2001
//
//
//---------------------------------------------------------------------------- 
#include "edginc/config.hpp" 
#define QLIB_RANDOM_CPP
#include "edginc/Addin.hpp"  
#include "edginc/Random.hpp" 
#include ext_hash_map

/* Taken from ALIB (and then interface modified) - for which thanks! 
 * Needed so can preserve more state if required to resume a sequence */
#define EDR_MC_NTAB  32
#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
/* This had to be promoted into header to allow structure def #define NTAB  32 */
#define NDIV  (1+IMM1/EDR_MC_NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0-EPS)

DRLIB_BEGIN_NAMESPACE 

// Interface
CClassConstSP const IRandom::TYPE = 
    CClass::registerInterfaceLoadMethod( 
        "IRandom", typeid(IRandom), load);

void IRandom::load(CClassSP& clazz){
    REGISTER_INTERFACE(IRandom, clazz);
    EXTENDS(IObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

// Default implementation of convert and clone methods 
CClassConstSP const RandomImpl::TYPE = CClass::registerClassLoadMethod(
    "RandomImpl", typeid(RandomImpl), load);

RandomImpl::RandomImpl(const CClassConstSP& clazz):
CObject(clazz){}

void RandomImpl::load(CClassSP& clazz){
    REGISTER(RandomImpl, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IRandom); 
}

IRandom* RandomImpl::convert(CClassConstSP requiredType) const{
    CClassConstSP thisType = this->getClass();
    if (requiredType == thisType){
        return copy(this);
    }
    throw ModelException("IRandom::convert", 
                         "Cannot convert a " + thisType->getName()
                         + " object into object of type " + requiredType->getName());
}

IObject* RandomImpl::clone() const{
    // lazy (but generic) way of copying random objects over
    IRandomSP random(dynamic_cast<IRandom*>(CObject::clone()));
    IRandom::StateSP state(getState());
    random->setState(state.get());
    return random.release();
}

// Random
//// holds state info for this Random object
class Random::State: public IRandom::State{
public: // only to this file though
    int      seed;          // for re-init remember orig seed
    int      seed2;
    int      iy;
    int      iv[EDR_MC_NTAB];
    int      iset;
    double   gset;
};

Random::Random():
RandomImpl(TYPE) { 
    seed = 97; 
    IsStored = false;
    init();
}

Random::Random(int seed):
RandomImpl(TYPE),
seed(seed){
    IsStored = false;
    init();
}

void Random::init() {
    int    j;
    int    k;

    if ( seed == 0 )
        MovingSeed = 1; // (long)time( (time_t *)NULL ); hmm...?
    else if (seed < 0 )
        MovingSeed = -(seed);
    else
        MovingSeed = seed;

    seed2 = MovingSeed;
    for( j=EDR_MC_NTAB+7; j>=0; j--) {
        k = (MovingSeed) / IQ1;
        MovingSeed = IA1 * (MovingSeed - k * IQ1) - k * IR1;
        if( MovingSeed < 0 )
            MovingSeed += IM1;
        if( j < EDR_MC_NTAB )
            iv[j] = MovingSeed;
    }
    iy = iv[0];

    // normal rands
    iset = 0;
    gset = 0;
}

void Random::fetch(int     numToFetch,
                   double* rands) {
    for(int i=0;i<numToFetch;i++){
        rands[i] = gasDev();
    }
}

void Random::skip(int  numToSkip) {
    for(int i=0;i<numToSkip;i++){
        gasDev();
    }
}

/** returns a State object capturing the state of the IRandom */
IRandom::State* Random::getState() const{
    Random::State* state = new Random::State();
    state->seed = MovingSeed;
    state->seed2 = seed2;
    state->iy = iy;
    memcpy(state->iv, iv, EDR_MC_NTAB * sizeof(int));
    state->iset = iset;
    state->gset = gset;
    return state;
}

/** restores the state of an IRandom */
void Random::setState(const IRandom::State* state){
    const Random::State& myState = 
        dynamic_cast<const Random::State&>(*state);
    MovingSeed = myState.seed;
    seed2 = myState.seed2;
    iy = myState.iy;
    memcpy(iv, myState.iv, EDR_MC_NTAB * sizeof(int));
    iset = myState.iset;
    gset = myState.gset;
}

double Random::uniform(){
    int    j;
    int    k;
    double temp;
    double rnd;

    k = (MovingSeed) / IQ1;
    MovingSeed = IA1 * (MovingSeed - k * IQ1) - k * IR1;
    if( MovingSeed < 0 )
        MovingSeed += IM1;
    k = seed2 / IQ2;
    seed2 = IA2 * (seed2 - k * IQ2) - k * IR2;
    if( seed2 < 0 )
        seed2 += IM2;
    j = (int)(iy / NDIV);
    iy = iv[j] - seed2;
    iv[j] = MovingSeed;
    if( iy < 1 )
        iy += IMM1;
    if( (temp = AM * iy) > RNMX )
        rnd = RNMX;
    else
        rnd = temp;

    return rnd;
}


double Random::gasDev(){
    double          v1;
    double          v2;
    double          rsq; 
    double          fac;
    double          rnd;

    if (iset == 0)
    {
        do 
        {
            rnd = uniform();
            v1 = 2.0 * rnd - 1.0;
            rnd = uniform();
            v2 = 2.0 * rnd - 1.0; 
            rsq = v1 * v1 + v2 * v2;
        }
        while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        /* NB NT optimisation fails if return statement is without brackets */
        return (v2 * fac);
    }
    else
    {
        iset = 0;
        return gset;
    }
}

typedef hash_map<CClassConstSP, 
                 Random::TIRandomFromRandom*, 
                 CClass::Hash> 
    IRandomFromRandomHash;

class Random::Helper{
public:
    // storage of methods for converting Random objects into other objects
    static IRandomFromRandomHash convertMethods;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Random, clazz);
        SUPERCLASS(RandomImpl);
        EMPTY_SHELL_METHOD(defaultRandom);
        FIELD(seed, "Seed");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultRandom(){
        return new Random();
    }
};

IRandomFromRandomHash Random::Helper::convertMethods;

/** Register a conversion method for a given type */
void Random::registerIRandomFromRandomMethod(CClassConstSP       targetClass,
                                             TIRandomFromRandom* method){
    if (!method){
        throw ModelException("Random::registerIRandomFromRandomMethod",
                             "Null method for class " + targetClass->getName());
    }
    Helper::convertMethods[targetClass] = method;
}


IRandom* Random::convert(CClassConstSP requiredType) const{
    // first check whether required type is same as this type
    // if yes, simply return deep copy
    CClassConstSP thisType = this->getClass();
    if (requiredType == thisType){
        return copy(this);
    }
    // otherwise, locate the required type in the hash table
    IRandomFromRandomHash::const_iterator iter = 
        Helper::convertMethods.find(requiredType);
    if (iter == Helper::convertMethods.end()){
        // fails if not not found
        throw ModelException("Random::convert", "Cannot convert Random object into "
                             "IRandom object of type " + requiredType->getName());
    }
    return iter->second(seed);
}

CClassConstSP const Random::TYPE = CClass::registerClassLoadMethod(
    "Random", typeid(Random), Helper::load);

template<> CClassConstSP const RandUniformDefault::TYPE = CClass::registerClassLoadMethod(
    "RandUniformDefault", typeid(RandUniformDefault), load);

template<> CClassConstSP const RandNormalDefault::TYPE = CClass::registerClassLoadMethod(
    "RandNormalDefault", typeid(RandNormalDefault), load);

template<> CClassConstSP const RandPoissonDefault::TYPE = CClass::registerClassLoadMethod(
    "RandPoissonDefault", typeid(RandPoissonDefault), load);

///////////////////////////////////////////////////////////////////////// 
// Publish the default RNG to the spreadsheet  
/////////////////////////////////////////////////////////////////////////
class GasDevAddin : public CObject {
public:
    static CClassConstSP const TYPE;

    // addin params
    smartPtr<Random>               rand;
    
    static double get(GasDevAddin* params){
        double aRandom;
        params->rand->fetch(1, &aRandom);
        return aRandom;
    }

    GasDevAddin(): CObject(TYPE) {}

    static IObject* defaultGasDevAddin(){
        return new GasDevAddin();
    }

    void validatePop2Object(){
        static const string method("GasDevAddin::validatePop2Object");
        if (!rand) {
            throw ModelException(method, "No Random generator supplied!");
        }
    }

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GasDevAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGasDevAddin);
        FIELD(rand,   "Random number generator");
        Addin::registerInstanceDoubleMethod(
            "GASDEV",
            Addin::UTILITIES,
            "Returns the next normal random number",
            TYPE,
            (Addin::DoubleMethod*)get);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
    
};

CClassConstSP const GasDevAddin::TYPE = CClass::registerClassLoadMethod(
    "GasDevAddin", typeid(GasDevAddin), GasDevAddin::load);


bool RandomLinkIn() {
    return (RandUniformDefault::TYPE != 0
            && RandNormalDefault::TYPE != 0
            && RandPoissonDefault::TYPE != 0);
}

DRLIB_END_NAMESPACE

