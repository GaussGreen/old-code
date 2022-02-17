//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimSeries.hpp
//
//   Description : Provides views of simulation dates for multiple assets
//                 excluding past values
//
//   Author      : Mark A Robson
//
//   Date        : 3 Sep 2001
//
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SIMSERIES_HPP
#define EDR_SIMSERIES_HPP

#include "edginc/DateTime.hpp"
#include "edginc/VirtualDestructorBase.hpp"

DRLIB_BEGIN_NAMESPACE
class IPastValues;

/** SimSeries class - idea is that it holds which dates the MC needs to
    simulate on (although excluding 'reference level' dates). An
    instance of the class is built on the fly through cooperation with
    performance type classes. It needs to support different dates per
    asset - given that upfront you don't know whether the dates will be
    different per asset, there seems little point in having different
    classes. 
    The other important point is that this class does not store past values */
class MARKET_DLL SimSeries: public virtual VirtualDestructorBase{
public:
    virtual ~SimSeries();

    /** Builds an empty SimSeries. Clients should use addDates() to register
        dates with the SimSeries. */
    SimSeries(int numAssets);
    
    /** returns simulation dates (of all assets) which are not
        strictly in the future */
    DateTimeArray getPastDates(const DateTime& valueDate) const;

    /** returns simulation dates (of all assets) which are not
        strictly in the future */
    DateTimeArray getFutureDates(const DateTime& valueDate) const;

    /** returns all the simulation dates (of all assets) */
    const DateTimeArray& getAllDates() const;

    /** Returns all the simulation dates of the given asset */
    const DateTimeArray& getDates(int iAsset) const;

    /** returns the number of future dates for given asset */
    int numFutureDates(const DateTime& valueDate,
                       int             iAsset) const;
    
    /** returns the number of dates (past and future) for given asset */
    int numDates(int iAsset) const;

    /** returns the number of assets */
    int getNumAssets() const;

    /** Returns the last date in the series - will throw an exception if no
        dates have been added */
    const DateTime& getLastDate() const;

    /** Returns the first date in the series - will throw an exception if no
        dates have been added */
    const DateTime& getFirstDate() const;

    /** Returns true if each asset has the same set of simulation dates */
    bool sameDatesPerAsset() const;

    /** Returns an array of integers denoting which assets have a simulation
        date given the index into getAllDates(). The number of the assets
        (0, 1, 2, ..., numAssets-1) is the same as the order in which the
        series is supplied to the constructor */
    const IntArray& assetsOnDate(int index) const;

    /** Creates a map as detailed in DateTime::createMapping ie
        equivalent to
        DateTime::createMapping(getDates(assetIndex), datesForMap, isTrivial)*/
    IntArray  createMap(int                   assetIndex,
                        const DateTimeArray&  datesForMap,
                        bool&                 isTrivial) const;

    /** Returns DoubleArray containing historic values for the given asset.
        Array returned is of length equal to the number of past values */
    DoubleArray getPastValues(const DateTime&    valueDate,
                              int                assetIndex,
                              const IPastValues* pastValues) const;

    /** adds supplied dates to list of dates held for each asset */
    void addDates(const DateTimeArray& dates);

    /** adds supplied cluster of dates to list of dates held on a per
        asset basis */
    void addDates(const DateTimeCluster& dates);

    /** adds supplied dates to list of dates held for specified asset */
    void addDates(int                  iAsset,
                  const DateTimeArray& dates);
    
private:
    // Throws exception if iAsset is out of range [0, numAssets)
    void checkIndex(const string      method,
                    int               iAsset) const;

    void buildCache() const; // modifies mutable fields

    /// fields ///
    mutable DateTimeArray    allDates;       // valid if rebuildCache is false
    int                      numAssets;
    vector<DateTimeArray>    datesPerAsset; 
    mutable vector<IntArray> assetsIndexes;
    mutable bool             areDatesSamePerAsset;
    mutable bool             rebuildCache;  /* => allDates and assetsIndexes
                                             needs rebuilding */

};


typedef smartConstPtr<SimSeries> SimSeriesConstSP;
typedef smartPtr<SimSeries> SimSeriesSP;
#ifndef QLIB_SIMSERIES_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<SimSeries>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<SimSeries>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<SimSeries>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<SimSeries>);
#endif

DRLIB_END_NAMESPACE

#endif

