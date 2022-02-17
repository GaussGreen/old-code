//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDSParSpreadsCache.hpp
//
//   Description : A cache to store <ICDSParSpreads, DefaultRates>.
//                 Refer to the "Caching of default rates (clean spreads)" 
//                 design document in the QLib DB for details.
//
//   Author      : Jose Hilera
//
//   Date        : 1 Nov 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICDSPARSPREADSCACHE_CPP
#define QLIB_ICDSPARSPREADSCACHE_CPP

#include "edginc/DefaultRates.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include ext_hash_map


DRLIB_BEGIN_NAMESPACE

class LinkedListCache;

/** Class containing the caches used for ICDSParSpreads objects. 
 * Essentially two LinkedListCaches, one to store the clean spreads of "big"
 * curves (quanto curves) and one for small curves (any other curves).
 * Internally, each LinkedListCache contains a hash_map table with:
 * - Key  = A "LinkedListElement *", which contains an ICDSParSpreadsConstSP. 
 *          Effectively the ICDSParSpreadsConstSP is the key used to access the
 *          hash_map - See the "Caching of default rates (clean spreads)"  
 *          design document in the QLib DB for details.
 * - Data = The Clean Spread Curves associated to the ICDSParSpreadsConstSP. */
class MARKET_DLL ICDSParSpreadsCache : virtual public VirtualDestructorBase {
public:
    /* Constructor and destructor */
    ICDSParSpreadsCache();
    ~ICDSParSpreadsCache();

    /** Get the Clean Spread Curves associated to the ICDSParSpreads:
     * - If already in the cache, return the Clean Spread Curves 
     * - Otherwise compute, store and return the Clean Spread Curves
     * The cache to use is determined using the entryType parameter */ 
    const DefaultRatesSP getEntry(const ICDSParSpreads* cds,
                                  const TypeOfEntry     entryType);

    /** Provides a facility for previously constructed DefaultRates
     *  objects to be added to the cache */
    const bool addEntry(const ICDSParSpreads* cds,
                        const TypeOfEntry     entryType,
                        const DefaultRatesSP  entry);
private:
    /* LinkedListCaches to store the curves and their clean spreads, and also 
     * keep track of which entries are most frequently used */
    LinkedListCache* linkedListCacheBig;
    LinkedListCache* linkedListCacheSmall;
};

typedef smartPtr<ICDSParSpreadsCache> ICDSParSpreadsCacheSP;

DRLIB_END_NAMESPACE

#endif
