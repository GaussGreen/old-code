//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCRandom.hpp
//
//   Description : Monte Carlo Random Numbers
//
//   Date        : March 2004
//
//----------------------------------------------------------------------------

#ifndef EDR_MCRANDOM_HPP
#define EDR_MCRANDOM_HPP

#include "edginc/config.hpp"
#include "edginc/MCPathConfig.hpp"
#include "edginc/Random.hpp"
#include "edginc/MCPathBase.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/MCCache.hpp"
#include "edginc/MSVC6Helper.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

/** Interface for random number simulation */
class MCARLO_DLL IMCRandom {
public:
    class MCARLO_DLL Callbacks {
    public:
        /** Configures object for antithetics */
        virtual void configureAntithetics() = 0;

        /** Configures object for nonAntithetics. We need a better name */
        virtual void configureNonAntithetics() = 0;

        /** Virtual destructor */
        virtual ~Callbacks() {}
    };

    /** Destructor */
    virtual ~IMCRandom () {}

    /** Invokes a random number generation of the appropriate dimension */
    virtual void generate(int pathIdx) = 0;

    /** Returns a reference to the simulated random numbers */
    virtual const DoubleMatrix& getRandomNumbers() const = 0;

    /** Returns if the generator is doing antithetics or not */
    virtual bool isAntithetic() const = 0;
};
DECLARE_REF_COUNT(IMCRandom);

/** Object that performs random number generation and correlation
    and prvides support for caching, antithetics and random access */
class MCARLO_DLL MCRandomGenCache{
private:
    class MCARLO_DLL IImpl{
    public:
        virtual void draw(int pathIdx) = 0;
        virtual const DoubleMatrix& getRandomNumbers() const = 0;
        virtual ~IImpl(){}
    };
    typedef refCountPtr<IImpl> IImplSP;
    IImplSP impl;

    template<class RandomCare> class Impl;

public:
    /** Full constructor */
    MCRandomGenCache(IMCRandom::Callbacks* generator,
                     const DependenceSP& dependence,
                     const IRandomSP& rand,
                     const MCRandomCacheSP& randomCache,
                     bool isCarefulRandoms,
                     int numDates,
                     int numFactors,
                     int numPastDates);

    /** Saves/restores state of random number generator as needed. Then
        invokes generateCorrelatedRands or does antithetics as required */
    void generate(int pathIdx);

    /** Returns numbers by reference */
    const DoubleMatrix& getRandomNumbers() const;
};
DECLARE_REF_COUNT(MCRandomGenCache);

/** Object that performs random number generation and correlation
    and prvides support for caching, antithetics and random access */
class MCRandom: virtual public IMCRandom {
private:
    class IImpl{
    public:
        virtual void draw(int pathIdx) = 0;
        virtual const DoubleMatrix& getRandomNumbers() const = 0;
        virtual void blockFillByDate(int blockSize) = 0;
#if CAREFUL_RAND_DEBUG
        virtual void debugCarefulRandoms(int pathIdx) = 0;
#endif
        virtual ~IImpl(){}
    };
    typedef refCountPtr<IImpl> IImplSP;
    IImplSP impl;

    typedef MCPathConfig::RandomCacheSP RandomCacheSP;

    template<class RandomCare, class AntitheticTreatment> class Impl;

public:
    /** Full constructor */
    MCRandom(IMCRandom::Callbacks* generator,
             const DependenceSP& dependence,
             const IRandomSP& rand,
             const RandomCacheSP& randomCache,
             bool isCarefulRandoms,
             int numDates,
             int numFactors,
             int numPastDates,
             bool doAntithetics = true);

    /** Invokes a random number generation of the appropriate dimension */
    virtual void generate(int pathIdx);

    /** Returns a reference to the simulated random numbers */
    virtual const DoubleMatrix& getRandomNumbers() const;

    /** Fetch all the random numbers by looping over dates, then over
        factors, and then over blockSize number of paths. Used for
        matching SRM3 test cases */
    virtual void blockFillByDate(int blockSize);

#if CAREFUL_RAND_DEBUG
    /** Writes random numbers to a file.
        Not virtual function so we can just not compile it if not needed */
    virtual void debugCarefulRandoms(int pathIdx);
#endif

    /** Returns if the generator is doing antithetics or not */
    virtual bool isAntithetic() const;
};
DECLARE_REF_COUNT(MCRandom);

/** Object that performs random number generation and correlation
    without caching but allows careful randoms and antithetics */
class MCRandomNoCache: virtual public IMCRandom {
private:
    class IImpl{
    public:
        virtual void draw(int pathIdx) = 0;
        virtual const DoubleMatrix& getRandomNumbers() const = 0;
        virtual ~IImpl(){}
    };
    typedef refCountPtr<IImpl> IImplSP;
    IImplSP impl;

    template<class RandomCare> class Impl;

public:
    /** Full constructor */
    MCRandomNoCache(IMCRandom::Callbacks* generator,
                    const DependenceSP& dependence,
                    const IRandomSP& rand,
                    bool isCarefulRandoms,
                    int numDates,
                    int numFactors,
                    int numPastDates);

    /** Invokes a random number generation of the appropriate dimension */
    virtual void generate(int pathIdx);

    /** Returns a reference to the simulated random numbers */
    virtual const DoubleMatrix& getRandomNumbers() const;

    /** Returns if the generator is doing antithetics or not */
    virtual bool isAntithetic() const;
};
DECLARE_REF_COUNT(MCRandomNoCache);

/** Object that performs random number generation and correlation
    without caching but allows careful randoms and antithetics */
class MCRandNormalNoCache {
private:
    class IImpl{
    public:
        virtual void draw(int pathIdx) = 0;
        virtual const IMCRandNormal::Path* getPath(int iAsset) = 0;
        // for backward compatibility -- for now
        virtual const DoubleMatrix& getRandomNumbers() const = 0;
        virtual ~IImpl(){}
    };
    typedef refCountPtr<IImpl> IImplSP;
    IImplSP impl;

    template<class RandomCare> class Impl;

public:
    /** Full constructor */
    MCRandNormalNoCache(IMCRandom::Callbacks* generator,
                        const DependenceSP& dependence,
                        const RandNormalDefaultSP& rand,
                        bool isCarefulRandoms,
                        int numDates,
                        int numFactors,
                        int numPastDates);

    /** Invokes a random number generation of the appropriate dimension */
    void generate(int pathIdx);

    const IMCRandNormal::Path* getPath(int iAsset);

    // for backward compatibility -- for now
    /** Returns a reference to the simulated random numbers */
    const DoubleMatrix& getRandomNumbers() const{
        return impl->getRandomNumbers();
    }
};
DECLARE_REF_COUNT(MCRandNormalNoCache);

class MCARLO_DLL MCRandNormal{
public:
    /** Full constructor */
    MCRandNormal(const RandNormalDefaultSP& rand):
    rand(rand){}

    double draw(){
        double rtn;
        rand->draw(1, &rtn);
        return rtn;
    }

private:
    RandNormalDefaultSP rand;
};
typedef refCountPtr<MCRandNormal> MCRandNormalSP;

class MCARLO_DLL MCRandPoisson{
public:
    /** Full constructor */
    MCRandPoisson(const RandPoissonDefaultSP& rand):
    rand(rand){}

    double draw(double mean){
        return rand->draw(mean);
    }

private:
    RandPoissonDefaultSP rand;
};
typedef refCountPtr<MCRandPoisson> MCRandPoissonSP;

DRLIB_END_NAMESPACE

#endif
