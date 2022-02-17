//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SimSeries.cpp
//
//   Description : Provides views of simulation dates for multiple assets
//                 excluding past values
//
//   Author      : Mark A Robson
//
//   Date        : 3 Sep 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#define QLIB_SIMSERIES_CPP
#include "edginc/SimSeries.hpp"
#include "edginc/PastValues.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE
/* SimSeries class - idea is that it holds which dates the MC needs to
   simulate on (although excluding 'reference level' dates). An
   instance of the class is built on the fly through cooperation with
   performance type classes. It needs to support different dates per
   asset - given that upfront you don't know whether the dates will be
   different per asset, there seems little point in having different
   classes. 
   The other important point is that this class does not store past values */

SimSeries::~SimSeries(){}

/** Builds an empty SimSeries. Clients should use addDates() to register
    dates with the SimSeries. */
SimSeries::SimSeries(int numAssets): numAssets(numAssets), 
    areDatesSamePerAsset(true), rebuildCache(false) {}

/** returns simulation dates (of all assets) which are not
    strictly in the future */
DateTimeArray SimSeries::getPastDates(const DateTime& valueDate) const{
    return valueDate.getPastDates(getAllDates());
}

/** returns simulation dates (of all assets) which are not
    strictly in the future */
DateTimeArray SimSeries::getFutureDates(const DateTime& valueDate) const{
    return valueDate.getFutureDates(getAllDates());
}

/** returns all the simulation dates (of all assets) */
const DateTimeArray& SimSeries::getAllDates() const{
    if (rebuildCache){
        buildCache();
    }
    return allDates;
}

/** Returns all the simulation dates of the given asset */
const DateTimeArray& SimSeries::getDates(int iAsset) const{
    checkIndex("SimSeries::getDates", iAsset);
    if (rebuildCache){
        buildCache();
    }
    return (sameDatesPerAsset()? allDates: datesPerAsset[iAsset]);
}

/** returns the number of future dates for given asset */
int SimSeries::numFutureDates(const DateTime& valueDate,
                              int             iAsset) const{
    checkIndex("SimSeries::numFutureDates", iAsset);
    return valueDate.numFutureDates(getDates(iAsset));
}
    
/** returns the number of dates (past and future) for given asset */
int SimSeries::numDates(int iAsset) const{
    checkIndex("SimSeries::numDates", iAsset);
    return getDates(iAsset).size();
}
 
/** returns the number of assets */
int SimSeries::getNumAssets() const{
    return numAssets;
}

/** Returns the last date in the series - will throw an exception if no
    dates have been added */
const DateTime& SimSeries::getLastDate() const{
    if (rebuildCache){
        buildCache();
    }
    if (allDates.empty()){
        throw ModelException("SimSeries::getLastDate", "Sim Series contains "
                             "no dates");
    }
    return getAllDates().back();
}

/** Returns the first date in the series - will throw an exception if no
    dates have been added */
const DateTime& SimSeries::getFirstDate() const{
    if (rebuildCache){
        buildCache();
    }
    if (allDates.empty()){
        throw ModelException("SimSeries::getFirstDate", "Sim Series contains "
                             "no dates");
    }
    return getAllDates().front();
}

/** Returns true if each asset has the same set of simulation dates */
bool SimSeries::sameDatesPerAsset() const{
    if (rebuildCache){
        buildCache();
    }
    return areDatesSamePerAsset;
}

/** Returns an array of integers denoting which assets have a simulation
    date given the index into getAllDates(). The number of the assets
    (0, 1, 2, ..., numAssets-1) is the same as the order in which the
    series is supplied to the constructor */
const IntArray& SimSeries::assetsOnDate(int index) const{
    if (rebuildCache){
        buildCache();
    }
    return assetsIndexes[index];
}

/** Creates an array of integers of length numDates(assetIndex). Each
    element in the array is the offset to the index of the array
    which corresponds to the next date in datesForMap. Eg
    Dates for an Asset: A B C D E F G H I J
    dateForMap:           B     E       I
    map:                1 0 2 1 0 3 2 1 0 1
    The output parameter isTrivial is true if datesForMap coincides with
    the dates for assetIndex. An exception is thrown if any of the dates
    in datesForMap are not in the sim series (for that asset). The
    supplied dates must be in order */
IntArray  SimSeries::createMap(int                   assetIndex,
                               const DateTimeArray&  datesForMap,
                               bool&                 isTrivial) const{
    const DateTimeArray& dates = getDates(assetIndex);
    return DateTime::createMapping(dates, datesForMap, isTrivial);
}

/** Returns DoubleArray containing historic values for the given asset.
    Array returned is of length equal to the number of past values */
DoubleArray SimSeries::getPastValues(const DateTime&    valueDate,
                                     int                assetIndex,
                                     const IPastValues* pastValues) const{
    const DateTimeArray& dates = getDates(assetIndex);
    return pastValues->getPastValues(dates, assetIndex, valueDate);
}

void SimSeries::buildCache() const{
    static const string routine("ISimSeries::General:buildCache");
    // tricky piece of code coming up ... currently untested...
    // need to merge date arrays and keep track of which asset's dates
    // are where
        
    int totalNumDates = 0;
    for (int i = 0; i < numAssets; i++){
        totalNumDates += datesPerAsset[i].size();
    }

    // create some storage space
    assetsIndexes = vector<IntArray>();
    assetsIndexes.reserve(totalNumDates);
    allDates.clear();
    allDates.reserve(totalNumDates);
    // set up array holding where we are in each array
    IntArray  posByAsset(numAssets);
    // loop through all dates
    int iDate = 0;
    for (int j = 0; j < totalNumDates; iDate++ /* j incremented below */){
        assetsIndexes.resize(iDate+1); // extend length by 1
        IntArray& assetsIndex = assetsIndexes[iDate];
        assetsIndex.reserve(numAssets); // reserve storage space
        // nextDate will point to the next date or be 0 until date found
        const DateTime* nextDate = 0;
        // loop over assets to see which date is earliest
        for (int iAsset = 0; iAsset < numAssets; iAsset++){
            if (posByAsset[iAsset] < datesPerAsset[iAsset].size()){
                const DateTime& date =
                    datesPerAsset[iAsset][posByAsset[iAsset]];
                if (nextDate == 0 || date.isLess(*nextDate)){
                    // record new date
                    nextDate = &date;
                    // reset assetsIndexes array
                    assetsIndex.resize(1);
                    assetsIndex[0] = iAsset;
                } else if (date.equals(*nextDate)){
                    assetsIndex.push_back(iAsset);
                }
            }
        }
        allDates.push_back(*nextDate); // save the date
        // then increment the indexes which have added a date
        for (int i = 0; i < assetsIndex.size(); i++){
            int iAsset = assetsIndex[i];
            posByAsset[iAsset]++;
            j++; // increment our overall counter
        }
    }
    areDatesSamePerAsset = true;
    for (int iAsset = 0; iAsset < numAssets; iAsset++){
        // sanity check
        if (posByAsset[iAsset] != datesPerAsset[iAsset].size()){
            throw ModelException(routine, "Internal error");
        }
        // have they all merged seamlessly?
        areDatesSamePerAsset = areDatesSamePerAsset &&
            (datesPerAsset[iAsset].size() == allDates.size());
    }
    rebuildCache = false;
}

/** adds supplied dates to list of dates held for each asset */
void SimSeries::addDates(const DateTimeArray& dates){
    if (!dates.empty()){
        DateTime::ensureIncreasing(dates, "SimSeries::addDates", false);
        // ensureIncreasing allows duplicates - which we do not allow here
        for(int j=1; j<dates.size(); j++) {
            if (dates[j-1]==dates[j]) {
                throw ModelException("SimSeries::addDates", 
                                     "The date list supplied contains a duplicate :" + dates[j].toString());
            }
        }

        if (areDatesSamePerAsset){
            // just merge into overall list
            if (allDates.empty()){
                allDates = dates;
            } else {
                allDates = DateTime::merge(allDates, dates);
            }
        } else {
            // merge into each assets dates
            for (unsigned int i = 0; i < datesPerAsset.size(); i++){
                datesPerAsset[i] = DateTime::merge(datesPerAsset[i], dates);
            }
            rebuildCache = true;
        }
    }
}

/** adds supplied dates to list of dates held for specified asset */
void SimSeries::addDates(int                  iAsset,
                         const DateTimeArray& dates){
    checkIndex("SimSeries::addDates", iAsset);
    if (!dates.empty()){
        DateTime::ensureIncreasing(dates, "SimSeries::addDates asset #" + Format::toString(iAsset+1), false);
        // ensureIncreasing allows duplicates - which we do not allow here
        for(int j=1; j<dates.size(); j++) {
            if (dates[j-1]==dates[j]) {
                throw ModelException("SimSeries::addDates", 
                                     "The date list supplied contains a duplicate :" + dates[j].toString());
            }
        }
        if (areDatesSamePerAsset){
            // need to copy allDates to each assets list of dates 
            datesPerAsset = vector<DateTimeArray>(numAssets);
            for (int i = 0; i < numAssets; i++){
                if (i == iAsset){
                    datesPerAsset[i] = allDates.empty()? 
                        dates: DateTime::merge(allDates, dates);
                } else {
                    datesPerAsset[i] = allDates;
                }
            }
            areDatesSamePerAsset = false;
        } else {
            DateTimeArray& datesForAsset = datesPerAsset[iAsset];
            datesForAsset = datesForAsset.empty()? 
                dates: DateTime::merge(datesForAsset, dates);
        }
        rebuildCache = true;
    }
}


/** adds supplied cluster of dates to list of dates held on a per
    asset basis */
void SimSeries::addDates(const DateTimeCluster& dates){
    for (int i = 0; i < dates.size(); i++){
        addDates(i, dates[i]);
    }
}

void SimSeries::checkIndex(const string      method,
                           int               iAsset) const {
    if (iAsset<0 || iAsset>=numAssets)
    {
        throw ModelException(method, 
                             "Asset index " + Format::toString(iAsset) +
                             " out of range [0, " + Format::toString(numAssets) + ")");
    }
}

DRLIB_END_NAMESPACE
