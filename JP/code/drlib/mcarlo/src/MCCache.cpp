//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCCache.cpp
//
//   Description : 
//
//   Date        : June 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MCCache.hpp"

DRLIB_BEGIN_NAMESPACE

MCCache::KeyUtils::Key MCCache::KeyUtils::getNewKey() {
    // Initialize key once at 0
    static Key key = 0;
    
    // Increment key and return
    key++;
    return key;
}


string MCCache::KeyUtils::keyToString(const Key& key) {
    return Format::toString(key);
}


//////////////////////////////////////////////////////////////////////


MCCacheManager::MCCacheManager() {}


MCCacheManager::~MCCacheManager() {}


int MCCacheManager::storagePerRun() const {
    // Aggregate memory requirements
    int memRequirement = 0;
    MCCacheMap::const_iterator iter = cacheMap.begin();
    for(; iter != cacheMap.end(); ++iter) {
        memRequirement += (*iter).second->storagePerRun();
    }

    return memRequirement;
}


void MCCacheManager::disallowAssets(const IntArray& assets) {
    // Propagate request to components
    MCCacheMap::iterator iter = cacheMap.begin();
    for(; iter != cacheMap.end(); ++iter) {
        MCCacheSP cache = (*iter).second;
        cache->disallowAssets(assets);
    }
}


void MCCacheManager::configureCache(int numRuns) {
    // Propagate request to components
    MCCacheMap::iterator iter = cacheMap.begin();
    for(; iter != cacheMap.end(); ++iter) {
        MCCacheSP cache = (*iter).second;
        cache->configureCache(numRuns);
    }
}


MCCacheSP MCCacheManager::getCache(const Key& key) {
    MCCacheMap::iterator iter = cacheMap.find(key);
    if(iter == cacheMap.end()) {
        // Return Null pointer.
        return MCCacheSP();
    } 

    // Return existing cache
    return (*iter).second;
}


void MCCacheManager::appendCache(const MCCacheSP& cache, const Key& key) {
    static const string routine = "MCCacheManager::appendCache";

    // Ensure that the key is not used
    MCCacheMap::iterator iter = cacheMap.find(key);
    if(iter != cacheMap.end()) {
        throw ModelException(routine, 
            "Cache with key " + MCCache::KeyUtils::keyToString(key) + "already exists." );
    }

    // Append it in the map
    cacheMap[key] = cache;
}


//////////////////////////////////////////////////////////////////////


MCRandomCacheSP MCRandomCache::createCache(MCCacheManager& cacheMgr,
                                           const MCCache::KeyUtils::Key& key,
                                           int numFactors,
                                           int numDates,
                                           bool allowWrite) {
    static const string routine = "MCRandomCache::createCache";

    try {
        MCRandomCacheSP randCache;

        MCCacheSP cache = cacheMgr.getCache(key);
        if(!cache) {
            // Create new random cache and store it in cache manager
            bool storeCorrelatedRands   = allowWrite ? true : false;
            bool storeUncorrelatedRands = allowWrite ? true : false;
            
            randCache = MCRandomCacheSP(new MCRandomCache(
                numFactors, numDates, storeCorrelatedRands, storeUncorrelatedRands, allowWrite));
            cacheMgr.appendCache(randCache, key);
        } else {
            // Recover cache and validate dimensions
            randCache = DYNAMIC_POINTER_CAST<MCRandomCache>(cache);

            if (randCache->numFactors != numFactors ||
                randCache->numDates   != numDates) {
                throw ModelException("Number of assets/dates has changed "
                                     "between greeks");
            }
        }

        return randCache;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


MCRandomCache::MCRandomCache(int numFactors,
                             int numDates,
                             bool storeCorrelatedRands, 
                             bool storeUncorrelatedRands,
                             bool allowWrite):
numFactors(numFactors), numDates(numDates), 
storeCorrelatedRands(storeCorrelatedRands), 
storeUncorrelatedRands(storeUncorrelatedRands), allowWrite(allowWrite),
correlatedRands(0), uncorrelatedRands(0), validCorrelatedRands(false), 
validUncorrelatedRands(false) { 

    static const string routine = "MCRandomCache::MCRandomCache";

    try {
        if(numFactors == 0) {
            throw ModelException("Number of factors cannot be zero.");
        }

        isDummy = !allowWrite;

        if(isDummy) {
            if(storeCorrelatedRands || storeUncorrelatedRands) {
                throw ModelException(
                    "Cannot request cahing of random numbers in a "
                    "read only cache");
            }
            numRandsPerRun = 0;
        } else {
            // A true cache
            numRandsPerRun = numFactors * numDates;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


MCRandomCache::~MCRandomCache() {
    delete[] correlatedRands;
    delete[] uncorrelatedRands;
}


int MCRandomCache::storagePerRun() const {
    int memRequirement = sizeof(double) * numRandsPerRun;
    return memRequirement;
}


void MCRandomCache::disallowAssets(const IntArray& assets) {
    // Nothing for dummy cache
    if(isDummy) {
        return;
    }

    // Allow both
    validUncorrelatedRands = storeUncorrelatedRands;
    validCorrelatedRands   = storeCorrelatedRands;

    // Disallow correlated random numbers if doing Phi
    if(assets.size() == numFactors) {
        validCorrelatedRands = false;
    } else {
        validCorrelatedRands = true;
    }

    // We are in tweaking mode so read only
    allowWrite = false;
}


void MCRandomCache::configureCache(int numRuns) {
    // Nothing for dummy cache
    if(isDummy) {
        return;
    }

    // We only intend to use the cache for non-anithetics so divide by 2.
    // This is not nice as we are assuming here that the cache will only 
    // be used for non-antithetic paths
    int myNumStates = (numRuns + 1) / 2;
    
    // only do anything if we actually caching anything
    if (storeCorrelatedRands || storeUncorrelatedRands) {
        // Allocate memory
        if (numRandsPerRun > 0) {
            int numDbls = numRandsPerRun * myNumStates;
            if (!correlatedRands && storeCorrelatedRands) {
                correlatedRands = new double[numDbls];
                for(int i = 0; i < numDbls; i++) {
                    correlatedRands[i] = 0.0;
                }
            }
            if (!uncorrelatedRands && storeUncorrelatedRands) {
                uncorrelatedRands = new double[numDbls];
                for(int i = 0; i < numDbls; i++) {
                    uncorrelatedRands[i] = 0.0;
                }
            }
        }
    }
}


void MCRandomCache::read(int runIdx, bool correlated, DoubleMatrix& randoms) const {
    // Go to appropriate position
    const double* rand = (correlated? correlatedRands: uncorrelatedRands) + 
        numRandsPerRun * runIdx;
    // Copy over
    for (int i = 0; i < numFactors; i++, rand += numDates) {
        memcpy(randoms[i], rand, sizeof(double) * numDates);
    }
}


void MCRandomCache::write(int runIdx, bool correlated, const DoubleMatrix& randoms) {
    // Go to appropriate position
    double* rand = (correlated? correlatedRands: uncorrelatedRands) + 
        numRandsPerRun * runIdx;
    // Copy over
    for (int i = 0; i < numFactors; i++, rand += numDates) {
        memcpy(rand, randoms[i], sizeof(double) * numDates);
    }
}


bool MCRandomCache::isValid(bool correlated) const {
    return (correlated? validCorrelatedRands: validUncorrelatedRands);
}


bool MCRandomCache::updateAllowed(bool correlated) const {
    return (allowWrite && (correlated? 
                           storeCorrelatedRands: storeUncorrelatedRands));
}

DRLIB_END_NAMESPACE
