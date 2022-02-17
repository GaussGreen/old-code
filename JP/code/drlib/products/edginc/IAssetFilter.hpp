//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IAssetFilter.hpp
//
//   Description : e.g. drop. No idea if it'll work more generally.
//
//   Date        : Apr 04
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IASSETFILTER_HPP
#define EDR_IASSETFILTER_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/IDoubleArray.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL IAssetFilter {
public:
    // Will assume moving forward in time. At successive filter dates the
    // set of "active" assets is changed.
    virtual void update(int iStep) = 0;

    // Need to know which assets are currently active.
    // The IntArray holds the indexes (as per the "filterMeasures") of measures
    // which are still active
    virtual const IntArray& getActiveAssets(bool doingPast) = 0;

    // Allow caching of past drops
    virtual void updateForPast() = 0;

    virtual ~IAssetFilter() {};
};
typedef refCountPtr<IAssetFilter> IAssetFilterSP;

class PRODUCTS_DLL IAssetFilterMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    // drop/etc at a distinct set of dates
    // Intended use is to add to sim dates
    virtual const DateTimeArray& getFilterDates() const = 0;

    // Somehow need to arrange a context in both asset and time 
    // dimensions - full set of assets we're handling and
    // the set of dates 

    // The allSimDates param allows step index to be used to access the filter
    // rather than needing to use actual dates.
    // "filterMeasures" is a working area of memory common to payoff and filter
    // It's size determines the range of indexes which are considered for filtering.
    virtual IAssetFilter* getFilter(const DateTimeArray& allSimDates,
                                    IDoubleArray*        filterMeasures) = 0;

    virtual ~IAssetFilterMaker() {};
private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IAssetFilterMaker> IAssetFilterMakerSP;

/*********************************************************************/
// No filtering at all - allows dependent structures more feedom 
// (e.g. full rainbow asset basket - see IAggregate::getAggregate
class PRODUCTS_DLL NoFilter : virtual public IAssetFilter {
public:
    NoFilter(IDoubleArray*         measures):
        activeAssets(measures->size()) {

        // all active, always
        for(int i=0; i<activeAssets.size(); i++) {
            activeAssets[i] = i;
        }
    };

    virtual void update(int iStep) {
        // Does nothing.
    }

    virtual const IntArray& getActiveAssets(bool doingPast) {
        return activeAssets;
    }

    virtual void updateForPast() {
        // Does nothing.
    }

private:
    IntArray activeAssets;
};

class PRODUCTS_DLL NoFilterMaker : public CObject,
                      virtual public IAssetFilterMaker {
public:
    static CClassConstSP const TYPE;
    friend class NoFilter;
    friend class NoFilterMakerHelper; // for hiding the registration

    virtual const DateTimeArray& getFilterDates() const{
        static DateTimeArray noDates(0);
        return noDates;
    };

    virtual IAssetFilter* getFilter(const DateTimeArray& allSimDates,
                                    IDoubleArray*        filterMeasures){
        return new NoFilter(filterMeasures); // needed for size
    };

private:
    NoFilterMaker(): CObject(TYPE) {} // for reflection
    NoFilterMaker(const NoFilterMaker& rhs); // not implemented
    NoFilterMaker& operator=(const NoFilterMaker& rhs); // not implemented
};
typedef smartPtr<NoFilterMaker> NoFilterMakerSP;
typedef array<NoFilterMakerSP, NoFilterMaker> NoFilterMakerArray;
typedef smartPtr<NoFilterMakerArray> NoFilterMakerArraySP;

/*********************************************************************/
// Simple Drop
class PRODUCTS_DLL DropFilterMaker : public CObject,
                        virtual public IAssetFilterMaker {
public:
    static CClassConstSP const TYPE;
    friend class DropFilter;
    friend class DropFilterMakerHelper; // for hiding the registration

    void validatePop2Object();

    virtual const DateTimeArray& getFilterDates() const;

    virtual IAssetFilter* getFilter(const DateTimeArray& allSimDates,
                                    IDoubleArray*        filterMeasures);

protected:
    int                  numBestToDrop;
    int                  numWorstToDrop;
    DateTimeArray        dropDates;

private:
    DropFilterMaker(): CObject(TYPE) {} // for reflection
    DropFilterMaker(const DropFilterMaker& rhs); // not implemented
    DropFilterMaker& operator=(const DropFilterMaker& rhs); // not implemented
};
typedef smartPtr<DropFilterMaker> DropFilterMakerSP;
typedef array<DropFilterMakerSP, DropFilterMaker> DropFilterMakerArray;
typedef smartPtr<DropFilterMakerArray> DropFilterMakerArraySP;

/*********************************************************************/
// Conditional drop
// So asset is dropped only if spot asset perf satisfies condition :
//  - up : best spot perf must >= barrier
//  - down : worst spot perf must <= barrier
// XXX Now spot perf is the drop measure, so that's ok, but 
// XXX how do I get hold of the barrier? I think it needs to be given 
// XXX to the maker as a separate param...
class PRODUCTS_DLL ConditionalDropFilterMaker : public CObject,
                                   virtual public IAssetFilterMaker {
public:
    static CClassConstSP const TYPE;
    friend class ConditionalDropFilter;
    friend class ConditionalDropFilterMakerHelper; // for hiding the registration

    void validatePop2Object();

    virtual const DateTimeArray& getFilterDates() const;

    virtual IAssetFilter* getFilter(const DateTimeArray& allSimDates,
                                    IDoubleArray*        filterMeasures);
protected:
    bool                 isDropBest; // only one direction at a time
    int                  numToDrop;  
    double               conditionalDropLevel; //schedule?
    DateTimeArray        dropDates;

private:
    ConditionalDropFilterMaker(): CObject(TYPE) {} // for reflection
    ConditionalDropFilterMaker(const ConditionalDropFilterMaker& rhs); // not implemented
    ConditionalDropFilterMaker& operator=(const ConditionalDropFilterMaker& rhs); // not implemented

};
typedef smartPtr<ConditionalDropFilterMaker> ConditionalDropFilterMakerSP;
typedef array<ConditionalDropFilterMakerSP, ConditionalDropFilterMaker> ConditionalDropFilterMakerArray;
typedef smartPtr<ConditionalDropFilterMakerArray> ConditionalDropFilterMakerArraySP;

DRLIB_END_NAMESPACE

#endif

