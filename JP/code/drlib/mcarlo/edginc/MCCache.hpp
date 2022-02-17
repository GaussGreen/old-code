//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCCache.hpp
//
//   Description : Monte Carlo Cache
//
//   Author      : Manos Venardos
//
//   Date        : 11 Jun 2004
//
//
//----------------------------------------------------------------------------
#ifndef MCCACHE_HPP
#define MCCACHE_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/Format.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

/** Cache Interface */
class MCARLO_DLL MCCache {
public:
    /** Definition of cache Key and utils */
    class MCARLO_DLL KeyUtils {
    public:
        typedef int Key;    //!< Keys for indexing caches

        /** Returns a new Globally Unique Id */
        static Key getNewKey();

        /** Converts key to string */
        static string keyToString(const Key& key);
    };

    /** Returns the number of bytes required for each run */
    virtual int storagePerRun() const = 0;

    /** Specifies which assets will not be read from the
        cache and will therefore be regenerated. */
    virtual void disallowAssets(const IntArray& assets) = 0;

    /** Sets up the number of simulation runs dimension of the cache */
    virtual void configureCache(int numRuns) = 0;

    /** Virtual destructor */
    virtual ~MCCache() {}
};

typedef refCountPtr<MCCache> MCCacheSP;


//////////////////////////////////////////////////////////////////////


/** Stores and manages an array of caches. Behaves like a cache too. */
class MCARLO_DLL MCCacheManager: virtual public MCCache {
private:
    typedef MCCache::KeyUtils::Key  Key;
    typedef map<Key, MCCacheSP>     MCCacheMap;

    MCCacheMap cacheMap;    //!< Collection of indexed MCCaches

public:
    /** Default constructor */
    MCCacheManager();

    /** Destructor */
    ~MCCacheManager();

    /** Returns the number of bytes required for each run */
    virtual int storagePerRun() const;

    /** Specifies which assets will not be read from the
        cache and will therefore be regenerated. */
    virtual void disallowAssets(const IntArray& assets);

    /** Sets up the number of simulation runs dimension of the cache */
    virtual void configureCache(int numRuns);

    /** Returns a cache corresponding to requested key */
    MCCacheSP getCache(const Key& key);

    /** Appends new cache. Fails if key already exists */
    void appendCache(const MCCacheSP& cache, const Key& key);
};

typedef refCountPtr<MCCacheManager> MCCacheManagerSP;


//////////////////////////////////////////////////////////////////////


class MCRandomCache;
typedef refCountPtr<MCRandomCache> MCRandomCacheSP;

/** Random number cache */
class MCARLO_DLL MCRandomCache: virtual public MCCache {
public:
    /** Returns cache from cache manager if it exists. Otherwise, it
        creates it, stores it in the manager and returns it */
    static MCRandomCacheSP createCache(MCCacheManager& cacheMgr,
                                       const MCCache::KeyUtils::Key& key,
                                       int numFactors,
                                       int numDates,
                                       bool allowWrite);

    /** Destructor */
    virtual ~MCRandomCache();

    /** Returns the number of bytes required for each run */
    virtual int storagePerRun() const;

    /** Specifies which assets will not be read from the
        cache and will therefore be regenerated. */
    virtual void disallowAssets(const IntArray& assets);

    /** Sets up the number of simulation runs dimension of the cache */
    virtual void configureCache(int numRuns);

    /** Populates supplied DoubleMatrix with cached randoms for
        specified run.  Either correlated or uncorrelated
        randoms depending on 'correlated'. Results
        undefined if cache not written to or runIdx has value
        outside that specified in constructor or randoms not of
        sufficient size (as specified in configure) */
    void read(int runIdx, bool correlated, DoubleMatrix& randoms) const;

    /** Copies supplied DoubleMatrix into cache for specified run.
        Randoms should either be correlated or uncorrelated randoms
        depending on 'correlated'. configure() must be invoked
        before write can be invoked. Results undefined if
        runIdx has value outside that specified in constructor or
        randoms not of sufficient size (as specified in configure) */
    void write(int runIdx, bool correlated, const DoubleMatrix& randoms);
    
    /** Are the contents of the cache valid */
    bool isValid(bool correlated) const;
    
    /** Are we allowed to write correlatedRands (true) /
        uncorrelated (false) randoms to the cache */
    bool updateAllowed(bool correlated) const;

private:
    /** Full constructor */
    MCRandomCache(int numFactors,
                  int numDates,
                  bool storeCorrelatedRands, 
                  bool storeUncorrelatedRands,
                  bool allowWrite);

    int numFactors;                 //!< Number of factors (usually assets)
    int numDates;                   //!< Number of sampling dates (common dates per factor)
    bool storeCorrelatedRands;      //!< Whether to store correlated randoms
    bool storeUncorrelatedRands;    //!< Whether to store uncorrelated randoms
    bool allowWrite;                //!< False means dummy cache

    bool isDummy;                   //!< Whether the cache will not be any caching at all
    int numRandsPerRun;             //!< Randoms per run = numFactors * numDates
    double* correlatedRands;        //!< Concatenated correlated random numbers
    double* uncorrelatedRands;      //!< Concatenated uncorrelated random numbers
    
    bool validCorrelatedRands;
    bool validUncorrelatedRands;
};


//////////////////////////////////////////////////////////////////////


/** Asset path cache */
template <class T>
class MCPathCache: virtual public MCCache {
public:
    typedef refCountPtr<MCPathCache<T> > MCPathCacheSP;
   
    /** Returns cache from cache manager if it exists. Otherwise, it
        creates it, stores it in the manager and returns it */
    static MCPathCacheSP createCache(MCCacheManager& cacheMgr,
                                     const MCCache::KeyUtils::Key& key,
                                     const IntArray&  numDatesPerAsset, 
                                     bool allowWrite) {
        static const string routine = "MCPathCache::createCache";

        try {
            MCPathCacheSP pathCache;
            
            MCCacheSP cache = cacheMgr.getCache(key);
            if(!cache) {
                // Create new path cache and store it in cache manager
                pathCache = MCPathCacheSP(new MCPathCache<T>(numDatesPerAsset, allowWrite));
                cacheMgr.appendCache(pathCache, key);
            } else {
                // Recover cache and validate dimensions
                pathCache = DYNAMIC_POINTER_CAST<MCPathCache>(cache);
                
                if(pathCache->numDatesPerAsset.size() != numDatesPerAsset.size()) {
                    throw ModelException("Number of assets has changed between greeks");
                }
                for(int iAsset = 0; iAsset < numDatesPerAsset.size(); iAsset++) {
                    if(pathCache->numDatesPerAsset[iAsset] != numDatesPerAsset[iAsset]) {
                        throw ModelException(
                            "Number of dates for asset " +
                            Format::toString(iAsset+1) + 
                            " has changed between greeks");
                    }
                }
            }

            return pathCache;

        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
    
    /** Destructor */
    virtual ~MCPathCache() {
        delete[] paths;
    }
    
    /** Returns the number of bytes required for each run */
    virtual int storagePerRun() const {
        int memRequirement = sizeof(T) * numValuesPerRun;
        return memRequirement;
    }

    /** Specifies which assets will not be read from the
        cache and will therefore be regenerated. */
    virtual void disallowAssets(const IntArray& assets) {
        // Nothing for dummy cache
        if(isDummy) {
            return;
        }
        
        // Mark everything as being valid again
        valid = BoolArray(numAssets, true);
    
        // Disallow the designated assets
        for (int i = 0; i < assets.size(); i++) {
            valid[assets[i]] = false;
        }

        // We are in tweaking mode so read only
        allowWrite = false;
    }

    /** Sets up the number of simulation runs dimension of the cache */
    virtual void configureCache(int numRuns) {
        // Nothing for dummy cache
        if(isDummy) {
            return;
        }

        int numDbls = numRuns * numValuesPerRun;
        if (numDbls > 0) {
            paths = new T[numDbls];
        }
    }

    /** Disallow all assets */
    void invalidateAllAssets() {
        // Disallow all assets
        valid = BoolArray(numAssets, false);
    }
    
    /** Validates that the cache is the right size. Very much a
        defensive piece of programming but it is only called once per
        construction of path generator so very cheap. Derived classes of
        MCPathConfig should call this with the appropriate data: length is
        number of points in cache. */
    void validateSize(int iAsset, int length) const {
        static const string routine = "MCPathCache::validateSize::validateSize";

        try {
            if(isValid(iAsset) || updateAllowed(iAsset)) {
                if (iAsset < 0 || iAsset >= assetOffsets.size() - 1) {
                    // this really shouldn't happen
                    throw ModelException("Asset index out of bounds");
                }
                int assetOffset = assetOffsets[iAsset];
                int cacheLength = assetOffsets[iAsset+1] - assetOffset;
                if (length != cacheLength) {
                    // hopefully this won't happen
                    throw ModelException("Internal error: cache misconfigured");
                }
            }
        } catch(exception& e) {
            throw ModelException(e, routine);
        }
    }
    
    /** Path stored in cache is written to dest. */
    void read(int iAsset, int runIdx, T* dest) const {
        // Go to appropriate position
        int assetOffset = assetOffsets[iAsset];
        int length = assetOffsets[iAsset+1] - assetOffset;
        const T* path = paths + (numValuesPerRun * runIdx) + assetOffset;
        // Copy over
        memcpy(dest, path, sizeof(T) * length);
    }
    
    /** Supplied path is written into cache */
    void write(int iAsset, int runIdx, const T* src) {
        // Go to appropriate position
        int assetOffset = assetOffsets[iAsset];
        int length = assetOffsets[iAsset+1] - assetOffset;
        T* path = paths + (numValuesPerRun * runIdx) + assetOffset;
        // Copy over
        memcpy(path, src, sizeof(T) * length);
    }

    /** Returns true if the cache has a valid path for this asset */
    bool isValid(int iAsset) const {
        return (valid[iAsset]);
    }
    
    /** returns true if the cache should be updated if isCacheValid = false */
    bool updateAllowed(int iAsset) const {
        return allowWrite; // don't support per asset [yet]
    }

private:
    // Mandatory fields
    IntArray  numDatesPerAsset;     //!< Number of dates for each asset
    bool      allowWrite;           //!< Set to false for dummy cache

    // Transient fields
    bool      isDummy;              //!< Whether the cache will not be any caching at all
    int       numAssets;            //!< Number of assets
    T*   paths;                     //!< Concatenated asset paths

    BoolArray valid;                //!< Specifies which assets are allowed    
    IntArray  assetOffsets;         //!< offsets for each asset
    int       numValuesPerRun;      //!< Sum numDatesPerAsset

    /** Full constructor */
    MCPathCache(const IntArray&  numDatesPerAsset, bool allowWrite):
    numDatesPerAsset(numDatesPerAsset), allowWrite(allowWrite), 
    numAssets(numDatesPerAsset.size()), paths(0) {
        static const string routine = "MCPathCache::MCPathCache";

        if(numAssets == 0) {
            throw ModelException(routine, "Number of assets cannot be zero.");
        }

        valid = BoolArray(numAssets, false);
        assetOffsets = IntArray(numAssets + 1);

        numValuesPerRun = 0;
        assetOffsets[0] = 0;
        for (int i = 0; i < numAssets; i++) {
            numValuesPerRun += numDatesPerAsset[i];
            assetOffsets[i+1] = assetOffsets[i] + numDatesPerAsset[i];
        }

        // If dummy override some values
        isDummy = !allowWrite;
        if(isDummy) {
            numValuesPerRun = 0;
        }
    }
};

typedef MCPathCache<double> SVGenSpotPathCache;
typedef SVGenSpotPathCache::MCPathCacheSP SVGenSpotPathCacheSP;

typedef MCPathCache<int> MCHitValuePathCache;
typedef MCHitValuePathCache::MCPathCacheSP MCHitValuePathCacheSP;

DRLIB_END_NAMESPACE

#endif
