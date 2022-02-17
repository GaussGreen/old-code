//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : random.hpp
//
//   Description : Provide generation and optional storage of random numbers
//                 Both uniform, normal and correlated normal
//
//   Date        : May 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/TypeConvert.hpp"

#ifndef EDR_RANDOM_HPP
#define EDR_RANDOM_HPP

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

class IRandom;
typedef smartPtr<IRandom> IRandomSP;
#ifndef QLIB_RANDOM_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<IRandom>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<IRandom>);
#endif

class UTIL_DLL IRandom: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** captures the state of an IRandom ie can use to reproduce the
        same sequence of random numbers */
    class UTIL_DLL State{
    public:
        virtual ~State(){} // for ease and performance
    protected:
        State(){} // for ease and performance
    };
    typedef refCountPtr<State> StateSP;
    typedef refCountPtr<vector<StateSP> > StateArraySP;

    /** Initializes the genarator */
    virtual void init() = 0;

    /** returns a State object capturing the state of the IRandom */
    virtual State* getState() const = 0;

    /** restores the state of an IRandom */
    virtual void setState(const State* state) = 0;

    /** Converts this object to an instance of the requiredType. Throws an
        exception if a conversion to the required Type is not supported */
    virtual IRandom* convert(CClassConstSP requiredType) const = 0;

    // for backward compatibility

    /** Fetches numToFetch random number.
        Beware there is no bound checking */
    virtual void fetch(int     numToFetch,
                       double* rands){
        throw ModelException("IRandom::fetch", "internal error");
    }

    /** Skips numToSkip deviates */
    virtual void skip(int numToSkip){
        throw ModelException("IRandom::skip", "internal error");
    }

private:

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

/** Implementation of the convert and clone methods
    of IRandom */
class UTIL_DLL RandomImpl: public CObject,
                  virtual public IRandom {
public:
    static CClassConstSP const TYPE;

    /** Need to preserve state of random number generator */
    virtual IObject* clone() const;

    /** By default, returns a copy of 'this' object if the required
        type is the same as 'this' type. Fails, otherwise. */
    virtual IRandom* convert(CClassConstSP requiredType) const;

protected:
    RandomImpl(const CClassConstSP& clazz);

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

/** Old default normal random generator */
class UTIL_DLL Random : public RandomImpl {
    // for uniform rands
    int      seed;          // for re-init remember orig seed
    int      MovingSeed; // $unregistered
    int      seed2; // $unregistered
    int      iy; // $unregistered
    int      iv[EDR_MC_NTAB]; // $unregistered

    // for normal rands
    int      iset; // $unregistered
    double   gset; // $unregistered

    // for later
    bool     IsStored; // $unregistered

    //// holds state info for this Random object
    class State;

public:
    static CClassConstSP const TYPE;

    Random(int seed);

    void init();

    void fetch(int     numToFetch,
               double* rands);

    /** returns a State object capturing the state of the IRandom */
    virtual IRandom::State* getState() const;

    /** restores the state of an IRandom */
    virtual void setState(const IRandom::State* state);

    /** Skips numToSkip deviates */
    virtual void skip(int numToSkip);

    /** Support for conversion of objects from Randoms. This defines
        the function that the Random class demands of a type.  */
    typedef IRandom* (TIRandomFromRandom)(int seed);

    /** Register a conversion method for a given type. */
    static void registerIRandomFromRandomMethod(CClassConstSP       targetClass,
                                                TIRandomFromRandom* method);

    /** Converts this object to an instance of the requiredType. Throws an
        exception if a conversion to the required Type is not supported.
        Note that the supplied object is a smart pointer to this. The
        converted object should be stored in the object parameter. */
    virtual IRandom* convert(CClassConstSP requiredType) const;

private:
    double uniform();

    double gasDev();

    class Helper;
    friend class Helper;

    Random();
};

/** Null fetch interface -- should be used when virtuality is not needed */
class UTIL_DLL IRandNullFect{};

/** Default uniform random number generator.
    Gotten from Alib */
template<class IFetch = IRandNullFect>
class RandUniform: public RandomImpl,
                   public virtual IFetch{
public:
    static CClassConstSP const TYPE;

    class IVirtualFetch{
    public:
        /** Fetches numToFetch random number.
            Beware there is no bound checking. */
        // Renamed to avoid confusion with IRandom's virtual fetch method
        virtual void draw(int     numToFetch,
                          double* rands) = 0;
    };

    RandUniform(int seed):
    RandomImpl(TYPE), seed(seed), IsStored(false){}

    RandUniform():
    RandomImpl(TYPE), seed(97), IsStored(false){}

    void validatePop2Object(){
        init();
    }

    /** Initializes the genarator */
    virtual void init() {
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
    }

    /** Fetches numToFetch random number.
        Beware there is no bound checking.
        Can be virtual or not depending on the derived interface */
    void draw(int     numToFetch,
              double* rands) {
        for(int i=0;i<numToFetch;i++){
            rands[i] = uniform();
        }
    }

    /** Skips numToSkip deviates */
    virtual void skip(int numToSkip) {
        for(int i=0;i<numToSkip;i++){
            uniform();
        }
    }

    /** returns a State object capturing the state of the IRandom */
    virtual IRandom::State* getState() const{
        State* state = new State();
        state->seed = MovingSeed;
        state->seed2 = seed2;
        state->iy = iy;
        memcpy(state->iv, iv, EDR_MC_NTAB * sizeof(int));
        return state;
    }

    /** restores the state of an IRandom */
    virtual void setState(const IRandom::State* state){
        const State& myState = dynamic_cast<const State&>(*state);
        MovingSeed = myState.seed;
        seed2 = myState.seed2;
        iy = myState.iy;
        memcpy(iv, myState.iv, EDR_MC_NTAB * sizeof(int));
    }

private:
    class State: public IRandom::State{
    friend class RandUniform<IFetch>;
    private:
        int      seed;          // for re-init remember orig seed
        int      seed2;
        int      iy;
        int      iv[EDR_MC_NTAB];
    };

    double uniform(){
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

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RandUniform, clazz);
        SUPERCLASS(RandomImpl);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(seed, "Seed");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new RandUniform();
    }

    // for uniform rands
    int      seed;          // for re-init remember orig seed
    int      MovingSeed; // $unregistered
    int      seed2; // $unregistered
    int      iy; // $unregistered
    int      iv[EDR_MC_NTAB]; // $unregistered

    // for later
    bool     IsStored; // $unregistered
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class IFetch> CClassConstSP const RandUniform<IFetch>::TYPE =
CClass::templateRegisterClass(typeid(RandUniform<IFetch>));
#endif

typedef RandUniform<> RandUniformDefault;
typedef smartPtr<RandUniformDefault> RandUniformDefaultSP;
#ifndef QLIB_RANDOM_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<RandUniformDefault>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<RandUniformDefault>);
#endif

/** Default normal random number generator.
    Really Numerical Recipes' gasdev */
template<class IFetch = IRandNullFect,
         class RandUni = RandUniformDefault,
         class RandUniSP = smartPtr<RandUni> >
class RandNormal : public RandomImpl,
                   public virtual IFetch {
public:
    static CClassConstSP const TYPE;

    class IVirtualFetch{
    public:
        /** Fetches numToFetch random number.
            Beware there is no bound checking. */
        // Renamed to avoid confusion with IRandom's virtual fetch method
        virtual void draw(int     numToFetch,
                          double* rands) = 0;
    };

    RandNormal():
    RandomImpl(TYPE){}

    RandNormal(const RandUniSP& uniRand):
    RandomImpl(TYPE), uniRand(uniRand){}

    void validatePop2Object(){
        init();
    }

    /** Initializes the genarator */
    virtual void init() {
        uniRand->init();
        iset = 0;
        gset = 0;
    }

    /** Fetches numToFetch random number.
        Beware there is no bound checking.
        Can be virtual or not depending on the derived interface */
    void draw(int     numToFetch,
              double* rands) {
        for(int i=0;i<numToFetch;i++){
            rands[i] = gasDev();
        }
    }

    /** Skips numToSkip deviates */
    virtual void skip(int  numToSkip) {
        for(int i=0;i<numToSkip;i++){
            gasDev();
        }
    }

    /** returns a State object capturing the state of the IRandom */
    virtual IRandom::State* getState() const{
        State* state = new State();
        state->uState = IRandom::StateSP(uniRand->getState());
        state->iset = iset;
        state->gset = gset;
        return state;
    }

    /** restores the state of an IRandom */
    virtual void setState(const IRandom::State* state){
        const State& myState = dynamic_cast<const State&>(*state);
        uniRand->setState(myState.uState.get());
        iset = myState.iset;
        gset = myState.gset;
    }

private:
    //// holds state info for this Random object
    class State: public IRandom::State{
    public: // only to this class though
        IRandom::StateSP     uState;
        int                  iset;
        double               gset;
    };

    double gasDev(){
        double          v1;
        double          v2;
        double          rsq;
        double          fac;
        double          rnd;

        if (iset == 0)
        {
            do
            {
                uniRand->draw(1, &rnd);
                v1 = 2.0 * rnd - 1.0;
                uniRand->draw(1, &rnd);
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

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RandNormal, clazz);
        SUPERCLASS(RandomImpl);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(uniRand,  "Normal random number generator");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new RandNormal();
    }

    RandUniSP   uniRand;
    int         iset; // $unregistered
    double      gset; // $unregistered
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class IFetch, class RandUni, class RandUniSP> CClassConstSP 
const RandNormal<IFetch, RandUni, RandUniSP>::TYPE =
CClass::templateRegisterClass(typeid(RandNormal<IFetch, RandUni, RandUniSP>));
#endif

typedef RandNormal<> RandNormalDefault;
typedef smartPtr<RandNormalDefault> RandNormalDefaultSP;
#ifndef QLIB_RANDOM_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<RandNormalDefault>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<RandNormalDefault>);
#endif

/** Poisson random number generator */
template<class IFetch = IRandNullFect,
         class RandUni = RandUniformDefault,
         class RandUniSP = smartPtr<RandUni> >
class RandPoisson : public RandomImpl,
                    public virtual IFetch {
public:
    static CClassConstSP const TYPE;

    class IVirtualFetch{
    public:
        /** Fetches numToFetch random number.
            Beware there is no bound checking */
        virtual double draw(double  mean) = 0;
    };

    RandPoisson():
    RandomImpl(TYPE){}

    RandPoisson(const RandUniSP& uniRand):
    RandomImpl(TYPE), uniRand(uniRand){}

    void validatePop2Object(){
        init();
    }

    /** Initializes the genarator */
    virtual void init() {
        uniRand->init();
        oldm = (-1.0);
    }

    /** Fetches numToFetch random number.
        Beware there is no bound checking */
    double draw(double  mean) {
        return poidev(mean);
    }

    /** returns a State object capturing the state of the IRandom */
    virtual IRandom::State* getState() const{
        State* state = new State();
        state->uState = IRandom::StateSP(uniRand->getState());
        state->sq = sq;
        state->alxm = alxm;
        state->g = g;
        state->oldm = oldm;
        return state;
    }

    /** restores the state of an IRandom */
    virtual void setState(const IRandom::State* state){
        const State& myState = dynamic_cast<const State&>(*state);
        uniRand->setState(myState.uState.get());
        sq = myState.sq;
        alxm = myState.alxm;
        g = myState.g;
        oldm = myState.oldm;
    }

private:
    //// holds state info for this Random object
    class State: public IRandom::State{
    public: // only to this file though
        IRandom::StateSP uState;
        double    sq;
        double    alxm;
        double    g;
        double    oldm;
    };

    double ran1(long* /*notused*/){
        double rnd;
        uniRand->draw(1, &rnd);
        return rnd;
    }

    /** Returns as a double-point number an integer value that is a random deviate drawn from a
        Poisson distribution of mean xm (Numerical Recipes) */
    double poidev(double xm){
        double em,t,y;
        long bogus;
        long* idum = &bogus;
        if (xm < 12.0) { // Use direct method.
        if (xm != oldm) {
            oldm=xm;
            g=exp(-xm); // If xm is new, compute the exponential.
        }
        em = -1;
        t=1.0;
        do { // Instead of adding exponential deviates it is equivalent
             // to multiply uniform deviates. We never
             // actually have to take the log, merely compare
             // to the pre-computed exponential.
            ++em;
            t *= ran1(idum);
        } while (t > g);
        } else { // Use rejection method.
            if (xm != oldm) { // If xm has changed since the last call, then precompute
                              // some functions that occur below. oldm=xm;
                sq=sqrt(2.0*xm);
                alxm=log(xm);
                g=xm*alxm-gammln(xm+1.0);
                // The function gammln is the natural log of the gamma function, as given in 6.1.
            }
            do {
                do { // y is a deviate from a Lorentzian comparison function.
                    y=tan(Maths::PI*ran1(idum));
                    em=sq*y+xm; // em is y, shifted and scaled.
                } while (em < 0.0); // Reject if in regime of zero probability.
                em=floor(em); // The trick for integer-valued distributions.
                t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                // The ratio of the desired distribution to the comparison function; we accept or
                // reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
                // that t never exceeds 1.
            } while (ran1(idum) > t);
        }
        return em;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RandPoisson, clazz);
        SUPERCLASS(RandomImpl);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(uniRand,  "Poisson random number generator");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultCtor(){
        return new RandPoisson();
    }

    RandUniSP uniRand;
    double    sq; // $unregistered
    double    alxm; // $unregistered
    double    g; // $unregistered
    double    oldm; // $unregistered
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class IFetch, class RandUni, class RandUniSP> CClassConstSP 
const RandPoisson<IFetch, RandUni, RandUniSP>::TYPE =
CClass::templateRegisterClass(typeid(RandPoisson<IFetch, RandUni, RandUniSP>));
#endif

typedef RandPoisson<> RandPoissonDefault;
typedef smartPtr<RandPoissonDefault> RandPoissonDefaultSP;
#ifndef QLIB_RANDOM_CPP
EXTERN_TEMPLATE(class UTIL_DLL_SP smartPtr<RandPoissonDefault>);
#else
INSTANTIATE_TEMPLATE(class UTIL_DLL smartPtr<RandPoissonDefault>);
#endif

// fwd declaration
class MCRandPoisson;
class MCRandNormal;

class UTIL_DLL IMCRandNormal{
public:
    /** Non-virtual reader for access to generated random deviates
        on a 'iAsset/iFactor' basis */
    class UTIL_DLL Path{
    public:
        double operator[](int step) const{
            return deviates[step];
        }

        int size() const{
            return nbDeviates;
        }

    protected:
        Path(int nbDeviates,
             const double* deviates):
        nbDeviates(nbDeviates),
        deviates(deviates){}

    private:
        int nbDeviates;
        const double* deviates;
    };
    // make const ptr + array visible only
    typedef refCountPtr<const Path> PathSP;
    typedef vector<PathSP> PathArray;
};

DRLIB_END_NAMESPACE

#undef EDR_MC_NTAB
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NDIV
#undef EPS
#undef RNMX

#endif




