//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCHelperDateTimeCache.hpp
//
//   Description :
//
//   Date        : Wed Apr 12 2006
//
//   Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/smartPtr.hpp"
#include <algorithm>
#include "edginc/QMCHelperDateTimeCache.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////


/** Global cache of all Dates we want to refer to via an iterator */
QMCHelperDateTimeCache::DateTimeSet * QMCHelperDateTimeCache::globalCache()
{
    static QMCHelperDateTimeCache::DateTimeSetSP  g; // make sure the memory will be eventually freed.
        if (g.get() == NULL)
            g = DateTimeSetSP(new DateTimeSet);
        return g.get();
}

// Return megre of all DateTimeArrays submitted directly
DateTimeArray QMCHelperDateTimeCache::extractDirect() const
{
        DateTimeArray r;
        for(IterSet::const_iterator it = directEntries.begin(); it != directEntries.end(); ++it)
                r.push_back(*(*it));
        return r;
}
/** Return megre of all DateTimeArrays submitted via SP */
DateTimeArray QMCHelperDateTimeCache::extractIndirect() const
{
        return DateTime::merge(indirectEntries);
}

/** internal function to invalidate "iterators" array that we use for quick getIdx */
// clears "iterators" vector
void QMCHelperDateTimeCache::invalidate(void)
{
        if (isValid) {
                isValid = false;
                iterators.clear();
                vector<DTSiter>().swap(iterators); // "swap" trick to forcibly release memory
        }
}

/** internal function that makes sure all dates are in cache, and object's dates are in the vector */
// populates "iterators" vector
//TODO: we can do better here: remove from directEntries thouse that are among Indirect ones
void QMCHelperDateTimeCache::initialize(void)
{
        DateTimeSet * cache = globalCache();
        for(vector<DateTimeArrayConstSP>::iterator arr = indirectEntries.begin();
                        arr != indirectEntries.end();
                        ++arr)
            cache->insert((*arr)->begin(), (*arr)->end());

        DateTimeArray dates = getDates();

        iterators.clear();
        iterators.reserve(dates.size());

        for(DateTimeArray::iterator it = dates.begin(); it != dates.end(); ++it)
            iterators.push_back(cache->find(*it));
        isValid = true;
}

/** Add an unshared DateTimeArray; have to copy it right away */
void QMCHelperDateTimeCache::add(const DateTimeArray& dates)
{
        if (!dates.empty()) {
                DateTimeSet * cache = globalCache();
                for(DateTimeArray::const_iterator it = dates.begin(); it != dates.end(); ++it)
                        directEntries.insert(directEntries.begin(), cache->insert(cache->begin(), *it));
                invalidate();
        }
}

/** Add an implicit DateTimeArray; to save memory, just store it untill we need it */
void QMCHelperDateTimeCache::add(DateTimeArrayConstSP dates)
{
        if (dates.get()) {
                indirectEntries.push_back(dates);
                invalidate();
        }
}

DateTimeArray   QMCHelperDateTimeCache::getDates()
{
        return DateTime::merge(extractDirect(), extractIndirect()); // keeps space efficiency
}

// Just like string::npos, indicates the DateTime is Not in cache OR not in getDates() array
//const size_t QMCHelperDateTimeCache::npos = size_t(-1);

int          QMCHelperDateTimeCache::getIdx(const DateTime& date)
{
        // intelligently returns date.find(getDates()) without keeping array
        if (!isValid)
                initialize();

        DateTimeSet * cache = globalCache();
        DTSiter loc = cache->find(date);
//        ASSERT(loc != cache->end());    // FIXME: remove later

        if (loc == cache->end())  // serious error: Date we didn't see yet
            return npos; // should we disable this possibility ?

        // Find loc in the vector using the fact it is a sorted array
        pair<vector<DTSiter>::iterator, vector<DTSiter>::iterator > range =
                equal_range(iterators.begin(),
                            iterators.end(),
                            loc, SetIterComparator());

        if (distance(range.first, range.second) == 1) // check if we found it
            return distance(iterators.begin(), range.first); // return offset
        else
            return npos;
}

DateTime     QMCHelperDateTimeCache::getDate(int idx) //< return getDates()[idx] from cache
{
    ASSERT(idx != npos);

    if (!isValid)
        initialize();
    return *(iterators[idx]);
}

QMCHelperDateTimeCache::QMCHelperDateTimeCache() :
        isValid(false)
{
}

/* Cleanup memory used by iterators. We can always reconstruct them, so we don't care to clear the vector */
/* Should save memory in SimpathI instrument case when all dates are passed via smart pointers and some of the DateTimeCaches are not needed after finalize()
*/

void QMCHelperDateTimeCache::cleanup()
{
    invalidate();
}

#define DUMPSIZE(x) cerr << #x << " p=" << (&(x)) << " " << typeid(*this).name() <<  " size= " << (x.size()) << " * " << sizeof(*(x.begin())) <<  endl

void QMCHelperDateTimeCache::debug() {
    DUMPSIZE( iterators);
    DUMPSIZE( directEntries);
    DUMPSIZE( indirectEntries);
//    DUMPSIZE( (*globalCache()));
}
        
void QMCHelperDateTimeCache::trim(const DateTime& date) 
{
    DateTimeSet * cache = globalCache();
    DTSiter loc = cache->insert(cache->begin(), date);

    if (!isValid) 
        initialize();
    iterators.erase(upper_bound(iterators.begin(), iterators.end(), loc, SetIterComparator()), iterators.end());

}


DRLIB_END_NAMESPACE
