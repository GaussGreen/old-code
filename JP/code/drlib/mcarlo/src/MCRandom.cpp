//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCRandom.cpp
//
//   Description : Monte Carlo Random Numbers
//
//   Date        : March 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MCRandom.hpp" 

#include <fstream>

DRLIB_BEGIN_NAMESPACE  

class DiceRecorder {
public:
    DiceRecorder(char* name) : ou(name) {}

    ~DiceRecorder() {
        ou.close();
    }

    void Record(int id, const DoubleMatrix& num) {
        ou << "ID " << id << " " << num.numCols() << " " << num.numRows() << std::endl;
        for (int col = 0; col < num.numCols(); ++col) {
            for (int row = 0 ; row < num.numRows(); ++row)
                ou << " " << num[col][row];
            ou << "\n";
        }
    }

private:
    std::ofstream ou;
};

static DiceRecorder &Recorder() {
    static DiceRecorder it("random.out");
    return it;
}
// Suite of policies which, when composed, make up a Monte Carlo Random class

struct DummyRandom{ 
    template <class Rand> 
    class In{ 
    public:
        typedef smartPtr<Rand> RandSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP& rand,
           int           numDates,
           int           numFactors,
           int           numPastDates):
        rand(rand), 
        numDates(numDates), 
        numFactors(numFactors), 
        numPastDates(numPastDates){}
 
        virtual ~In(){}

        void draw(DoubleMatrix& randoms) {
            // do nthg 
        } 

    protected:
        RandSP  rand;                   //!< Random number generator
        int     numDates;               //!< Number of dates (rows in randoms matrix)
        int     numFactors;             //!< Number of factors (#columns in randoms matrix)
        int     numPastDates;           //!< Number of past dates to skip when retrieving
    };
};

struct CarefulRandom{ 
    template <class Rand> 
    class In{ 
    public:
        typedef smartPtr<Rand> RandSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP& rand,
           int           numDates,
           int           numFactors,
           int           numPastDates):
        rand(rand), 
        numDates(numDates), 
        numFactors(numFactors), 
        numPastDates(numPastDates){}
 
        virtual ~In(){}

        void draw(DoubleMatrix& randoms) {
            // Generate by asset then by date from latest, so passing
            // sample dates drops a contiguous block of randoms "from
            // the back" and can maintain same randoms on same sample
            // dates more easily
            for (int iStep = numDates - 1; iStep >= 0; iStep--) {
                for (int iAsset = 0; iAsset < numFactors; iAsset++) {
                    rand->draw(1, &randoms[iAsset][iStep]);
                }
            }
            /* In order to use the same random numbers from one day to the
               next (eg calculating price each day in prod) we get the
               random number generator to alter its state so that it was as 
               if the numbers had been generated. This ensures consistency
               when the next path is generated */
            rand->skip(numFactors * numPastDates);
        }

    protected:
        RandSP  rand;                   //!< Random number generator
        int     numDates;               //!< Number of dates (rows in randoms matrix)
        int     numFactors;             //!< Number of factors (#columns in randoms matrix)
        int     numPastDates;           //!< Number of past dates to skip when retrieving
    };
};


// For backward compatibility, I need a specialization for IRandom
// since I want to invoke the virtual fetch method (as opposed to draw)
template<> void CarefulRandom::In<IRandom>::draw(DoubleMatrix& randoms) {
    // Generate by asset then by date from latest, so passing
    // sample dates drops a contiguous block of randoms "from
    // the back" and can maintain same randoms on same sample
    // dates more easily
    for (int iStep = numDates - 1; iStep >= 0; iStep--) {
        for (int iAsset = 0; iAsset < numFactors; iAsset++) {
            rand->fetch(1, &randoms[iAsset][iStep]);
        }
    }
    /* In order to use the same random numbers from one day to the
       next (eg calculating price each day in prod) we get the
       random number generator to alter its state so that it was as 
       if the numbers had been generated. This ensures consistency
       when the next path is generated */
    rand->skip(numFactors * numPastDates);
}

struct CarelessRandom{
    template <class Rand>
    class In{
    public:
        typedef smartPtr<Rand> RandSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP& rand,
           int           numDates,
           int           numFactors,
           int           numPastDates):
        rand(rand), 
        numDates(numDates), 
        numFactors(numFactors), 
        numPastDates(numPastDates){}

        virtual ~In(){}

        void draw(DoubleMatrix& randoms) {
            for (int iAsset = 0; iAsset < numFactors; iAsset++) {
                rand->draw(randoms.numRows(), randoms[iAsset]);
            }
        }

    protected:
        RandSP  rand;                   //!< Random number generator
        int     numDates;               //!< Number of dates (rows in randoms matrix)
        int     numFactors;             //!< Number of factors (#columns in randoms matrix)
        int     numPastDates;           //!< Number of past dates to skip when retrieving
    };
};

// For backward compatibility, I need a specialization for IRandom
// since I want to invoke the virtual fetch method (as opposed to draw)
template<> void CarelessRandom::In<IRandom>::draw(DoubleMatrix& randoms) {
    for (int iAsset = 0; iAsset < numFactors; iAsset++) {
        rand->fetch(randoms.numRows(), randoms[iAsset]);
    }
}  

template <class Rand, class RandomCare, class CacheType>
class RawRandCaching: public MSVC6Helper::Apply1<RandomCare, Rand>::Base{
private:
    typedef typename MSVC6Helper::Apply1<RandomCare, Rand>::Base RandCare;

public:
    typedef typename RandCare::RandSP RandSP;
    typedef refCountPtr<CacheType> CacheTypeSP;
 
    RawRandCaching(const RandSP&          rand,
                   int                    numDates,
                   int                    numFactors,
                   int                    numPastDates,
                   const CacheTypeSP&     randomCache):
    RandCare(rand, numDates, numFactors, numPastDates),
    correlated(false),   // uncorrelated rands
    randomCache(randomCache){} 

    void draw(int randPathIdx, DoubleMatrix& randoms) {
        if (isValid()) {
            // read off cache
            read(randPathIdx, randoms);
        } else {
            // draw
            RandCare::draw(randoms);
            // save in cache
            write(randPathIdx, randoms);
        }
    }

    bool hasValidCache() const{
        return isValid();
    }

private:
    const bool  correlated;

    bool isValid() const{
        return randomCache->isValid(correlated);
    }

    void read(int randPathIdx, DoubleMatrix& randoms) const{
        randomCache->read(randPathIdx, correlated, randoms);
    }

    void write(int randPathIdx, DoubleMatrix& randoms){
        if (randomCache->updateAllowed(correlated)) {
            randomCache->write(randPathIdx, correlated, randoms);
        }
    }

protected:
    CacheTypeSP randomCache;         //!< Cache for random numbers
};

struct CorrelatedRandCaching {
    template <class Rand, class RandomCare, class CacheType>
    class In: public RawRandCaching<Rand, RandomCare, CacheType> {
    private:
        typedef RawRandCaching<Rand, RandomCare, CacheType> RawRandCache;

    public:
        typedef typename RawRandCache::RandSP RandSP;
        typedef typename RawRandCache::CacheTypeSP CacheTypeSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP&          rand,
           int                    numDates,
           int                    numFactors,
           int                    numPastDates,
           IMCRandom::Callbacks*  generator,
           const DependenceSP&    dependence,
           const CacheTypeSP&     randomCache):
        RawRandCache(rand, numDates, numFactors, numPastDates, randomCache),
        correlated(true),   // correlated rands
        generator(generator), dependence(dependence){}

        void draw(int randPathIdx, DoubleMatrix& randoms) {
            if (isValid()) {
                // read off cache
                read(randPathIdx, randoms);
            } else {
                // draw
                RawRandCache::draw(randPathIdx, randoms);
                // correlate
                dependence->correlateSeries(randoms, randPathIdx);
                // save in cache
                write(randPathIdx, randoms);
            }
            configurePathGen();
        }

        bool hasValidCache() const{
            return RawRandCache::hasValidCache() || this->isValid();
        }

    private:
        const bool correlated;

        bool isValid() const{
            return this->randomCache->isValid(correlated);
        }

        void read(int randPathIdx, DoubleMatrix& randoms) const{
            this->randomCache->read(randPathIdx, correlated, randoms);
        }

        void write(int randPathIdx, DoubleMatrix& randoms){
            if (this->randomCache->updateAllowed(correlated)) {
                this->randomCache->write(randPathIdx, correlated, randoms);
            }
        }

        void configurePathGen(){
            // Allow pathGen to do extra work
            if (this->generator) {
                this->generator->configureNonAntithetics();
            }
        }

    protected:
        IMCRandom::Callbacks* generator;  //!< Pointer to original pathGen: used to trigger methods
        DependenceSP dependence;          //!< Dependence object for correlating randoms
    };
};


struct NoRandCaching{ 
    template <class Rand, class RandomCare, class CacheType> 
    class In: public MSVC6Helper::Apply1<RandomCare, Rand>::Base{
    private: 
        typedef typename MSVC6Helper::Apply1<RandomCare, Rand>::Base RandCare;
  
    public: 
        typedef typename RandCare::RandSP RandSP; 
        typedef refCountPtr<CacheType> CacheTypeSP;  

        // seems to be needed by gcc !
        In(){}

        In(const RandSP&          rand,
           int                    numDates,
           int                    numFactors,
           int                    numPastDates,
           IMCRandom::Callbacks*  generator,
           const DependenceSP&    dependence,
           const CacheTypeSP&     notused):
        RandCare(rand, numDates, numFactors, numPastDates),
        generator(generator), dependence(dependence){}

        void draw(int notused, DoubleMatrix& randoms) {
            RandCare::draw(randoms);
            dependence->correlateSeries(randoms, notused);
            configurePathGen();
        }

        bool hasValidCache() const{
            return false;
        }

    private:
        void configurePathGen(){
            // Allow pathGen to do extra work
            if(generator) {
                generator->configureNonAntithetics();
            }
        }

    protected:
        IMCRandom::Callbacks* generator;  //!< Pointer to original pathGen: used to trigger methods
        DependenceSP dependence;          //!< Dependence object for correlating randoms
    };
};

struct NullCacheType{};
struct Antithetic{
    template <class Rand, class RandomCare, class Caching, class CacheType = NullCacheType>
    class In: public MSVC6Helper::Apply3<Caching, Rand, RandomCare, CacheType>::Base {
    private:
        typedef typename MSVC6Helper::Apply3<Caching, Rand, RandomCare, CacheType>::Base Cache;

    public:
        typedef typename Cache::RandSP RandSP;
        typedef typename Cache::CacheTypeSP CacheTypeSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP&         rand,
           int                   numDates,
           int                   numFactors,
           int                   numPastDates,
           IMCRandom::Callbacks* generator,
           const DependenceSP&   dependence,
           const CacheTypeSP&     randomCache = CacheTypeSP()):
        Cache(rand, numDates, numFactors,
              numPastDates, generator, 
              dependence, randomCache){ 
            // Initialize transient fields
            randoms = DoubleMatrix(numFactors, numDates);
            lastPathIdx = -1;
        }

        int isAntithetic(int pathIdx) const{
            return (pathIdx % 2);
        }

        void draw(int pathIdx) {
            bool pathsOutOfOrder = pathIdx != lastPathIdx + 1;
            ASSERT(!pathsOutOfOrder || Cache::hasValidCache());
            int randPathIdx = getRandPathIdx(pathIdx);
            if (isAntithetic(pathIdx)) {
                if (pathsOutOfOrder){
                    Cache::draw(randPathIdx, randoms);
                }
                randoms.negate();
                configurePathGen();
            } else {
                Cache::draw(randPathIdx, randoms);
            }
            lastPathIdx = pathIdx;
        }

    private:
        int getRandPathIdx(int pathIdx) const{
            return (pathIdx / 2);
        }

        void configurePathGen(){
            // Allow pathGen to do extra work
            if (this->generator) {
                this->generator->configureAntithetics();
            }
        }

    protected:
        int          lastPathIdx;
        DoubleMatrix randoms;
    };
};

struct NoAntithetic{
    template <class Rand, class RandomCare, class Caching, class CacheType = NullCacheType>
    class In: public MSVC6Helper::Apply3<Caching, Rand, RandomCare, CacheType>::Base {
    private:
        typedef typename MSVC6Helper::Apply3<Caching, Rand, RandomCare, CacheType>::Base Cache;

    public:
        typedef typename Cache::RandSP RandSP;
        typedef typename Cache::CacheTypeSP CacheTypeSP;

        // seems to be needed by gcc !
        In(){}

        In(const RandSP&          rand,
           int                    numDates,
           int                    numFactors,
           int                    numPastDates,
           IMCRandom::Callbacks*  generator,
           const DependenceSP&    dependence,
           const CacheTypeSP&     randomCache = CacheTypeSP(   )):
        Cache(rand, numDates, numFactors,
              numPastDates, generator, 
              dependence, randomCache){ 
            // Initialize transient fields
            randoms = DoubleMatrix(numFactors, numDates);
            lastPathIdx = -1;
        }

        void draw(int pathIdx) {
            bool pathsOutOfOrder = pathIdx != lastPathIdx + 1;
            ASSERT(!pathsOutOfOrder || Cache::hasValidCache());
            Cache::draw(pathIdx, randoms);            
            lastPathIdx = pathIdx;
        }

    private:
        int getRandPathIdx(int pathIdx) const{
            return pathIdx;
        }

        void configurePathGen(){
            // Allow pathGen to do extra work
            if (this->generator) {
                this->generator->configureNonAntithetics();
            }
        }

    protected:
        int          lastPathIdx;
        DoubleMatrix randoms;
    };
};       

// MCRandomGenCache
template<class RandomCare>
class MCRandomGenCache::Impl: public MSVC6Helper::Apply4<Antithetic, 
                                                         IRandom, 
                                                         RandomCare, 
                                                         CorrelatedRandCaching, 
                                                         MCRandomCache>::Base,
                              public IImpl{
private:
    typedef typename MSVC6Helper::Apply4<Antithetic, IRandom, RandomCare, 
                                         CorrelatedRandCaching, MCRandomCache>::Base Anti;

public:
    Impl(IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const IRandomSP& rand,
         const MCRandomCacheSP& randomCache,
         int numDates,
         int numFactors,
         int numPastDates):
    Anti(rand, numDates, numFactors, numPastDates, 
         generator, dependence, randomCache){}
       
private:
    virtual void draw(int pathIdx){
        Anti::draw(pathIdx);
    }

    virtual const DoubleMatrix& getRandomNumbers() const{ 
        return Anti::randoms;
    }
};

MCRandomGenCache::MCRandomGenCache(
         IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const IRandomSP& rand,
         const MCRandomCacheSP& randomCache,
         bool isCarefulRandoms,
         int numDates,
         int numFactors,
         int numPastDates){

    typedef Impl<CarefulRandom> CarefulImpl;
    typedef Impl<CarelessRandom> CarelessImpl;
    if (isCarefulRandoms){
        impl = IImplSP(
            new CarefulImpl(generator,
                            dependence,
                            rand,
                            randomCache,
                            numDates,
                            numFactors,
                            numPastDates));
    }
    else{
        impl = IImplSP(
            new CarelessImpl(generator,
                            dependence,
                            rand,
                            randomCache,
                            numDates,
                            numFactors,
                            numPastDates));
    }
}
         
/** Saves/restores state of random number generator as needed. Then
    invokes generateCorrelatedRands or does antithetics as required */
void MCRandomGenCache::generate(int pathIdx){
    impl->draw(pathIdx);
}
         
/** Returns numbers by reference */
const DoubleMatrix& MCRandomGenCache::getRandomNumbers() const{
    return impl->getRandomNumbers();
}

// MCRandom
template<class RandomCare, class AntitheticTreatment>
class MCRandom::Impl: public MSVC6Helper::Apply4<AntitheticTreatment, 
                                                 IRandom, 
                                                 RandomCare, 
                                                 CorrelatedRandCaching, 
                                                 MCPathConfig::RandomCache>::Base,
                      public IImpl{
private:
    typedef typename MSVC6Helper::Apply4<AntitheticTreatment, IRandom, RandomCare, 
                                         CorrelatedRandCaching, MCPathConfig::RandomCache>::Base Anti;

public:
    Impl(IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const IRandomSP& rand,
         const RandomCacheSP& randomCache,
         int numDates,
         int numFactors,
         int numPastDates):
    Anti(rand, numDates, numFactors, numPastDates, 
         generator, dependence, randomCache){}

private:
    virtual void draw(int pathIdx){
        Anti::draw(pathIdx);
    }

    virtual const DoubleMatrix& getRandomNumbers() const{
        return Anti::randoms;
    }

    virtual void blockFillByDate(int blockSize){
        CarefulRandom* carefulRand = dynamic_cast<CarefulRandom*>(this);
        if (carefulRand){
            throw ModelException("MCRandom::blockFillByDate", "Not supported "
                                 "with careful randoms methodology");
        }
        DummyRandom* dummyRand = dynamic_cast<DummyRandom*>(this);
        if (dummyRand) {
            throw ModelException("MCRandom::blockFillByDate", "Not supported "
                                 "with dummy randoms methdogoldy");
        }
        int totalNumPaths = 0;
        int numPaths = this->randomCache->getNumRandPaths();
        while (totalNumPaths < numPaths){
            for (int iDate = 0; iDate < this->numDates; iDate++){
                for (int iFactor = 0; iFactor < this->numFactors; iFactor++){
                    for (int iPath = 0; iPath < blockSize; iPath++){
                        double randNum;
                        this->rand->fetch(1, &randNum);
                        this->randomCache->write(totalNumPaths+iPath, iFactor, 
                                                 iDate, false, randNum);
                    }
                }
            }
            totalNumPaths += blockSize;
        }
        this->randomCache->configure(false);
    }

#if CAREFUL_RAND_DEBUG
    virtual void debugCarefulRandoms(int pathIdx) {
        ofstream  debugfile("randoms.txt", ios_base::app);
        debugfile << "Iteration : " << pathIdx << "\n";
        for(int j = 0; j < numFactors; j++) {
            debugfile << "Asset : " << j << "\n";
            for (int i = 0; i < numDates; i++) {
                debugfile << randoms[j][i] << "\n";
            }
        }
    } 
#endif
};

/** Full constructor */
MCRandom::MCRandom(IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const IRandomSP& rand,
         const RandomCacheSP& randomCache,
         bool isCarefulRandoms,
         int numDates,
         int numFactors,
         int numPastDates,
         bool doAntithetics){
    // Configure cache here
    randomCache->configure(numFactors, numDates);
    if (!rand) { // LocalCorr
        impl = IImplSP(
            new Impl<DummyRandom, NoAntithetic>(generator,
                          dependence, 
                          rand,
                          randomCache, 
                          numDates, 
                          numFactors, 
                          numPastDates));
    } else if (isCarefulRandoms) { // CAREFUL
        if (doAntithetics) {
            impl = IImplSP(
                new Impl<CarefulRandom, Antithetic>(generator,
                                                    dependence,
                                                    rand,
                                                    randomCache,
                                                    numDates,
                                                    numFactors,
                                                    numPastDates));
        } else {
            impl = IImplSP(
                new Impl<CarefulRandom, NoAntithetic>(generator,
                                                      dependence,
                                                      rand,
                                                      randomCache,
                                                      numDates,
                                                      numFactors,
                                                      numPastDates));
        }
    } else { // CARELESS
        if (doAntithetics) {
            impl = IImplSP(
                new Impl<CarelessRandom, Antithetic>(generator,
                                                     dependence,
                                                     rand,
                                                     randomCache,
                                                     numDates,
                                                     numFactors,
                                                     numPastDates));
        } else {
            impl = IImplSP(
                new Impl<CarelessRandom, NoAntithetic>(generator,
                                                       dependence,
                                                       rand,
                                                       randomCache,
                                                       numDates,
                                                       numFactors,
                                                       numPastDates));
        }
    }
}

/** Invokes a random number generation of the appropriate dimension */
void MCRandom::generate(int pathIdx){
    impl->draw(pathIdx);
    //if (getenv("RECORD_RANDOM"))
    //    Recorder().Record(pathIdx, getRandomNumbers());
}
         
/** Returns a reference to the simulated random numbers */
const DoubleMatrix& MCRandom::getRandomNumbers() const{
    return impl->getRandomNumbers();
}

/** Fetch all the random numbers by looping over dates, then over 
    factors, and then over blockSize number of paths. Used for
    matching SRM3 test cases */
void MCRandom::blockFillByDate(int blockSize){
    impl->blockFillByDate(blockSize);
}

#if CAREFUL_RAND_DEBUG
/** Writes random numbers to a file.
    Not virtual function so we can just not compile it if not needed */
void MCRandom::debugCarefulRandoms(int pathIdx) {
    impl->debugCarefulRandoms(pathIdx);
}
#endif

/** Returns if the generator is doing antithetics or not */
bool MCRandom::isAntithetic() const{
    throw ModelException("MCRandomNoCache::isAntithetic", "not supported");
}

// MCRandomNoCache
template<class RandomCare>
class MCRandomNoCache::Impl: public MSVC6Helper::Apply3<Antithetic, IRandom, RandomCare, NoRandCaching>::Base,
                             public IImpl{
private:
    typedef typename MSVC6Helper::Apply3<Antithetic, IRandom, RandomCare, NoRandCaching>::Base Anti;
public:
    Impl(IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const IRandomSP& rand,
         int numDates,
         int numFactors,
         int numPastDates):
    Anti(rand, numDates, numFactors, numPastDates, 
         generator, dependence){}
       
private:
    virtual void draw(int pathIdx){
        Anti::draw(pathIdx);
    }

    virtual const DoubleMatrix& getRandomNumbers() const{
        return Anti::randoms;
    }
};

MCRandomNoCache::MCRandomNoCache(IMCRandom::Callbacks* generator,
                const DependenceSP& dependence,
                const IRandomSP& rand,
                bool isCarefulRandoms,
                int numDates,
                int numFactors,
                int numPastDates){
    typedef Impl<CarefulRandom> CarefulImpl;
    typedef Impl<CarelessRandom> CarelessImpl;
    if (isCarefulRandoms){
        impl = IImplSP(
            new CarefulImpl(generator, 
                            dependence,
                            rand,
                            numDates,
                            numFactors,
                            numPastDates));
    }
    else{
        impl = IImplSP(
            new CarelessImpl(generator,
                            dependence,
                            rand,
                            numDates,
                            numFactors,
                            numPastDates));
    }
}
         
/** Invokes a random number generation of the appropriate dimension */
void MCRandomNoCache::generate(int pathIdx){
    impl->draw(pathIdx);
}

/** Returns a reference to the simulated random numbers */
const DoubleMatrix& MCRandomNoCache::getRandomNumbers() const{
    return impl->getRandomNumbers();
}

/** Returns if the generator is doing antithetics or not */
bool MCRandomNoCache::isAntithetic() const{
    throw ModelException("MCRandomNoCache::isAntithetic", "not supported");
}

// MCRandNormalNoCache
template<class RandomCare>
class MCRandNormalNoCache::Impl: public MSVC6Helper::Apply3<Antithetic, RandNormalDefault, RandomCare, NoRandCaching>::Base,
                                 public IImpl{
private:
    typedef typename MSVC6Helper::Apply3<Antithetic, RandNormalDefault, RandomCare, NoRandCaching>::Base Anti;

public:
    Impl(IMCRandom::Callbacks* generator,
         const DependenceSP& dependence,
         const RandNormalDefaultSP& rand,
         int numDates,
         int numFactors,
         int numPastDates):
    Anti(rand, numDates, numFactors, numPastDates, 
         generator, dependence){}
       
private:
    class Path: public IMCRandNormal::Path{
    public:
        Path(int nbDeviates,
             const double* deviates):
        IMCRandNormal::Path(nbDeviates, deviates){}
    };

    virtual void draw(int pathIdx){
        Anti::draw(pathIdx);
    }

    virtual const IMCRandNormal::Path* getPath(int iAsset){
        return new Path(Anti::numDates, Anti::randoms[iAsset]);
    }

    virtual const DoubleMatrix& getRandomNumbers() const{
        return Anti::randoms; 
    }
};

MCRandNormalNoCache::MCRandNormalNoCache(IMCRandom::Callbacks* generator,
                                         const DependenceSP& dependence,
                                         const RandNormalDefaultSP& rand,
                                         bool isCarefulRandoms,
                                         int numDates,
                                         int numFactors,
                                         int numPastDates){
    typedef Impl<CarefulRandom> CarefulImpl;
    typedef Impl<CarelessRandom> CarelessImpl;
    if (isCarefulRandoms){
        impl = IImplSP(
            new CarefulImpl(generator, 
                            dependence,
                            rand,
                            numDates,
                            numFactors,
                            numPastDates));
    }
    else{
        impl = IImplSP(
            new CarelessImpl(generator,
                            dependence,
                            rand,
                            numDates,
                            numFactors,
                            numPastDates));
    }
}
         
/** Invokes a random number generation of the appropriate dimension */
void MCRandNormalNoCache::generate(int pathIdx){
    impl->draw(pathIdx);
}

const IMCRandNormal::Path* MCRandNormalNoCache::getPath(int iAsset){
    return impl->getPath(iAsset);
}

DRLIB_END_NAMESPACE
