//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IAssetFilter.cpp
//
//   Description : e.g. drop
//
//   Date        : Apr 04
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/IAssetFilter.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE


void IAssetFilterMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IAssetFilterMaker, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const IAssetFilterMaker::TYPE = CClass::registerInterfaceLoadMethod(
    "IAssetFilterMaker", typeid(IAssetFilterMaker), load);


/*********************************************************************/

class NoFilterMakerHelper{
public:
    static IObject* defaultNoFilterMaker(){
        return new NoFilterMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(NoFilterMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAssetFilterMaker);
        EMPTY_SHELL_METHOD(defaultNoFilterMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const NoFilterMaker::TYPE =
CClass::registerClassLoadMethod("NoFilterMaker", 
                                typeid(NoFilterMaker), NoFilterMakerHelper::load);

/*********************************************************************/
// Drop
class DropFilter : virtual public IAssetFilter {
public:
    DropFilter(const DropFilterMaker*  maker,
               const DateTimeArray&    allSimDates,
               IDoubleArray*           measures):
        maker(maker), 
        activeAssets(measures->size()),
        activeAssetsSoFar(measures->size()),
        activeRemoves(maker->numBestToDrop + maker->numWorstToDrop),
        measures(measures),
        perfs(measures->size()) {

        bool isTrivial;
        dropDateMap = DateTime::createMapping(allSimDates, 
                                              maker->dropDates,
                                              isTrivial);
        // can check we have sensible numbers
        if (maker->dropDates.size() * 
            (maker->numBestToDrop + maker->numWorstToDrop) >= measures->size()) {
            throw ModelException("DropFilter::DropFilter",
                                 "Not enough assets to drop that many!");
        }
        // all active to start with
        for(int i=0; i<activeAssets.size(); i++) {
            activeAssets[i] = i;
            activeAssetsSoFar[i] = i;
        }
    };

    // Will assume moving forward in time. At successive filter dates the
    // set of "active" assets is changed.
    virtual void update(int iStep) {
        if (dropDateMap[iStep]==0) {
            int i;
            int n = activeAssets.size();
            // 1. sort the remaining active assets
            perfs.resize(n);
            for(i=0; i<n; i++) {
                int iAsset = activeAssets[i];
                perfs[i].iAsset = iAsset;
                perfs[i].perf = (*measures)[iAsset];
            }
            sort(perfs.begin(), perfs.end(), OrderPerfs());
            // 2. identify the numBestToDrop and numWorstToDrop
            for (i=0; i<maker->numBestToDrop; i++) {
                activeRemoves[i] = perfs[i].iAsset;
            }
            for(i=0; i<maker->numWorstToDrop; i++) {
                activeRemoves[maker->numBestToDrop+i] = perfs[n-1-i].iAsset;
            }
            // 3. remove them from activeAssets
            // This should leave activeAssets in order.
            if (activeRemoves.size()>0) {
                sort(activeRemoves.begin(), activeRemoves.end());
                // XXX not sure how portable this code is!
                IntArray::iterator elt=activeAssets.begin();
                i = 0;
                while (elt != activeAssets.end() && i<activeRemoves.size()) {
                    if (*elt == activeRemoves[i]) {
                        elt = activeAssets.erase(elt); // even for vector it is more correct to use the result of .erase()
                        // no need to increment elt since the erasure will have
                        // moved the array elements down one
                        i++;
                    } else {
                        elt++;
                    }
                }
#if 0
                for(i=0, elt=activeAssets.begin(); elt<activeAssets.end(); elt++) {
                    if (*elt == activeRemoves[i]) {
                        activeAssets.erase(elt);
                        i++;
                        if (i>=activeRemoves.size()) {
                            // done 
                            break;
                        }
                    }
                }
#endif
            }

            // sanity ...
            if (activeAssets.size()<1) {
                throw ModelException("DropFilter::update",
                                     "Dropped all assets!");
            }

        }
    }

    // Need to know which assets are currently active.
    // The IntArray holds the indexes (as per the "filterMeasures") of measures
    // which are still active
    virtual const IntArray& getActiveAssets(bool doingPast) {
        if (!doingPast) {
            activeAssets = activeAssetsSoFar;
        } 
        return activeAssets;
    }

    // Allow caching of past drops
    virtual void updateForPast() {
        activeAssetsSoFar = activeAssets;
    }

private:
    const DropFilterMaker* maker;
    IntArray               activeAssets;
    IntArray               activeAssetsSoFar; // for past
    IntArray               activeRemoves;
    IntArray               dropDateMap;
    IDoubleArray*          measures;

    // for sorting
    struct IndexedPerfs {
        double perf;
        int    iAsset;
    };
    
    vector<IndexedPerfs>     perfs;
    class OrderPerfs {
    public:
        int operator() (const IndexedPerfs& p1, const IndexedPerfs& p2) {
            return p1.perf > p2.perf;
        }
    };

};

/*********************************************************************/

void DropFilterMaker::validatePop2Object(){
    if (dropDates.size()<1) {
        throw ModelException("DropFilterMaker::validatePop2Object",
                             "Require some drop dates!");
    }
}

const DateTimeArray& DropFilterMaker::getFilterDates() const {
    return dropDates;
}

IAssetFilter* DropFilterMaker::getFilter(const DateTimeArray& allSimDates,
                                         IDoubleArray*        filterMeasures) {
    return new DropFilter(this,
                          allSimDates, 
                          filterMeasures); 

}

class DropFilterMakerHelper{
public:
    static IObject* defaultDropFilterMaker(){
        return new DropFilterMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DropFilterMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAssetFilterMaker);
        EMPTY_SHELL_METHOD(defaultDropFilterMaker);
        FIELD(numBestToDrop,      "Num of best performing assets to drop each drop date");
        FIELD(numWorstToDrop,     "Num of worst performing assets to drop each drop date");
        FIELD(dropDates,          "Dates when drop happens");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const DropFilterMaker::TYPE =
CClass::registerClassLoadMethod("DropFilterMaker", 
                                typeid(DropFilterMaker), DropFilterMakerHelper::load);


/*********************************************************************/
// ConditionalDrop
class ConditionalDropFilter : virtual public IAssetFilter {
public:
    ConditionalDropFilter(const ConditionalDropFilterMaker*  maker,
                          const DateTimeArray&               allSimDates,
                          IDoubleArray*                      measures):
        maker(maker), 
        activeAssets(measures->size()),
        activeAssetsSoFar(measures->size()),
        activeRemoves(maker->numToDrop),
        measures(measures),
        perfs(measures->size()) {

        bool isTrivial;
        dropDateMap = DateTime::createMapping(allSimDates, 
                                              maker->dropDates,
                                              isTrivial);
        // can check we have sensible numbers
        if (maker->dropDates.size() * maker->numToDrop >= measures->size()) {
            throw ModelException("ConditionalDropFilter::ConditionalDropFilter",
                                 "Not enough assets to drop that many!");
        }

        // all active to start with
        for(int i=0; i<activeAssets.size(); i++) {
            activeAssets[i] = i;
            activeAssetsSoFar[i] = i;
        }
    };

    // Will assume moving forward in time. At successive filter dates the
    // set of "active" assets is changed.
    virtual void update(int iStep) {
        if (dropDateMap[iStep]==0) {
            int i;
            int n = activeAssets.size();
            // 1. sort the remaining active assets
            // XXX can we keep perfs around from previous drop date?
            perfs.resize(n);
            for(i=0; i<n; i++) {
                int iAsset = activeAssets[i];
                perfs[i].iAsset = iAsset;
                perfs[i].perf = (*measures)[iAsset];
            }
            sort(perfs.begin(), perfs.end(), OrderPerfs());
            // 2. Decide if any get dropped
            bool conditionMet;
            if (maker->isDropBest) {
                conditionMet = (perfs.front().perf >= maker->conditionalDropLevel);
            } else {
                conditionMet = (perfs.back().perf <= maker->conditionalDropLevel);
            }
            if (conditionMet && maker->numToDrop>0) {
                // 3a. If so, identify the numToDrop 
                for (i=0; i<maker->numToDrop; i++) {
                    int j = maker->isDropBest?i:n-1-i;
                    activeRemoves[i] = perfs[j].iAsset;
                }
                // 3b. and remove them from activeAssets
                // This should leave activeAssets in order.
                sort(activeRemoves.begin(), activeRemoves.end()); // XXX Need this to sort low -> high.
                IntArray::iterator elt=activeAssets.begin();
                i = 0;
                while (elt != activeAssets.end() && i<activeRemoves.size()) {
                    if (*elt == activeRemoves[i]) {
                        elt = activeAssets.erase(elt); // even for vector it is more correct to use the result of .erase()
                        // no need to increment elt since the erasure will have
                        // moved the array elements down one
                        i++;
                    } else {
                        elt++;
                    }
                }
            }

            // sanity ...
            if (activeAssets.size()<1) {
                throw ModelException("ConditionalDropFilter::update",
                                     "Dropped all assets!");
            }

        }
    }

    // Need to know which assets are currently active.
    // The IntArray holds the indexes (as per the "filterMeasures") of measures
    // which are still active
    virtual const IntArray& getActiveAssets(bool doingPast) {
        if (!doingPast) {
            activeAssets = activeAssetsSoFar;
        } 
        return activeAssets;
    }

    // Allow caching of past drops
    virtual void updateForPast() {
        activeAssetsSoFar = activeAssets;
    }

private:
    const ConditionalDropFilterMaker* maker;
    IntArray               activeAssets;
    IntArray               activeAssetsSoFar; // for past
    IntArray               activeRemoves;
    IntArray               dropDateMap;
    IDoubleArray*          measures;

    // for sorting
    struct IndexedPerfs {
        double perf;
        int    iAsset;
    };
    
    vector<IndexedPerfs>     perfs;
    class OrderPerfs {
    public:
        int operator() (const IndexedPerfs& p1, const IndexedPerfs& p2) {
            return p1.perf > p2.perf;
        }
    };

};

/*********************************************************************/

void ConditionalDropFilterMaker::validatePop2Object(){
    if (dropDates.size()<1) {
        throw ModelException("ConditionalDropFilterMaker::validatePop2Object",
                             "Require some drop dates!");
    }
}

const DateTimeArray& ConditionalDropFilterMaker::getFilterDates() const {
    return dropDates;
}

IAssetFilter* ConditionalDropFilterMaker::getFilter(const DateTimeArray& allSimDates,
                                                    IDoubleArray*        filterMeasures) {
    return new ConditionalDropFilter(this,
                                     allSimDates, 
                                     filterMeasures); 

}

class ConditionalDropFilterMakerHelper{
public:
    static IObject* defaultConditionalDropFilterMaker(){
        return new ConditionalDropFilterMaker();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ConditionalDropFilterMaker, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAssetFilterMaker);
        EMPTY_SHELL_METHOD(defaultConditionalDropFilterMaker);
        FIELD(isDropBest,     "Else worst");
        FIELD(numToDrop,      "Num of assets that may be dropped if condition satisfied");
        FIELD(conditionalDropLevel, "Level which must be breached to allow drop");
        FIELD(dropDates,      "Dates when drop happens");
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};

CClassConstSP const ConditionalDropFilterMaker::TYPE =
CClass::registerClassLoadMethod("ConditionalDropFilterMaker", 
                                typeid(ConditionalDropFilterMaker), ConditionalDropFilterMakerHelper::load);

/*********************************************************************/

/* The ultimate wrapping of IAssetFilterMakers mainly for use in Pyramid 
 */
#define ASSET_FILTER_TYPE_NONE        "None"
#define ASSET_FILTER_TYPE_DROP        "Drop"
#define ASSET_FILTER_TYPE_COND_DROP   "ConditionalDrop"

class AssetFilterMakerWrapper : public CObject,
                                virtual public IAssetFilterMaker {
public: // how can I have this protected or private?
    string                              assetFilterMakerType;
    NoFilterMakerSP                     noFilterMaker;
    DropFilterMakerSP                   dropFilterMaker;
    ConditionalDropFilterMakerSP        conDropFilterMaker;

private:
    IAssetFilterMakerSP  realMaker;

public:
    static CClassConstSP const TYPE;

    const DateTimeArray& getFilterDates() const {
        return realMaker->getFilterDates();
    }

    virtual IAssetFilter* getFilter(const DateTimeArray& allSimDates,
                                    IDoubleArray*        filterMeasures){
        return realMaker->getFilter(allSimDates, filterMeasures);
    }

    // validation
    void validatePop2Object(){
        static const string routine = "AssetFilterMakerWrapper::validatePop2Object";

        if (assetFilterMakerType.empty()){
            throw ModelException(routine, "Blank AssetFilter Maker specified!");
        }
        if (assetFilterMakerType==ASSET_FILTER_TYPE_NONE) {
            if (noFilterMaker.get()) {
                realMaker = noFilterMaker;
            } else {
                throw ModelException(routine, "Expected NoFilter Maker but none supplied!");
            }
        } else if (assetFilterMakerType==ASSET_FILTER_TYPE_DROP) {
            if (dropFilterMaker.get()) {
                realMaker = dropFilterMaker;
            } else {
                throw ModelException(routine, "Expected Drop Filter Maker but none supplied!");
            }
        } else if (assetFilterMakerType==ASSET_FILTER_TYPE_COND_DROP) {
            if (conDropFilterMaker.get()) {
                realMaker = conDropFilterMaker;
            } else {
                throw ModelException(routine, "Expected Conditional Drop Filter Maker but none supplied!");
            }
        } else {
            throw ModelException(routine, "Unrecognised Asset Filter Maker " + assetFilterMakerType + 
                                 ". Expected " + ASSET_FILTER_TYPE_NONE + ", " + 
                                 ASSET_FILTER_TYPE_DROP + " or " + ASSET_FILTER_TYPE_COND_DROP);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(AssetFilterMakerWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IAssetFilterMaker);
        EMPTY_SHELL_METHOD(defaultAssetFilterMakerWrapper);
        FIELD(assetFilterMakerType, "Drop or ConditionalDrop");
        FIELD(noFilterMaker,  "No filtering");
        FIELD_MAKE_OPTIONAL(noFilterMaker);
        FIELD(dropFilterMaker,  "Simple drop");
        FIELD_MAKE_OPTIONAL(dropFilterMaker);
        FIELD(conDropFilterMaker,  "Conditional drop");
        FIELD_MAKE_OPTIONAL(conDropFilterMaker);
        FIELD(realMaker, "realMaker");
        FIELD_MAKE_TRANSIENT(realMaker);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
     
    // for reflection
    AssetFilterMakerWrapper(): CObject(TYPE){}

    static IObject* defaultAssetFilterMakerWrapper(){
        return new AssetFilterMakerWrapper();
    }
};

typedef smartPtr<AssetFilterMakerWrapper> AssetFilterMakerWrapperSP;

CClassConstSP const AssetFilterMakerWrapper::TYPE =
CClass::registerClassLoadMethod("AssetFilterMakerWrapper", 
                                typeid(AssetFilterMakerWrapper), load);

DRLIB_END_NAMESPACE

    

