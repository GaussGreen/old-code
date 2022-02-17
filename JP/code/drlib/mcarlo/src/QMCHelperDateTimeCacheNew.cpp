//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCHelperDateTimeCacheNew.hpp
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
#include "edginc/QMCHelperDateTimeCacheNew.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////////////


/** Global cache of all Dates we want to refer to via an iterator */
QMCHelperDateTimeCacheNew::DateTimeSet * QMCHelperDateTimeCacheNew::globalCache()
{
    static QMCHelperDateTimeCacheNew::DateTimeSetSP g; // make sure the memory will be eventually freed.
    if ( g.get() == NULL )
        g = DateTimeSetSP( new DateTimeSet );
    return g.get();
}


/** Add an unshared DateTimeArray; have to copy it right away */
void QMCHelperDateTimeCacheNew::add
    ( const DateTimeArray& dates )
{
    if ( !dates.empty() ) {
        
        if (frozen)
            thaw();

        DateTimeSet * cache = globalCache();
        for ( DateTimeArray::const_iterator it = dates.begin(); it != dates.end(); ++it )
                entries.insert( entries.begin(), cache->insert( cache->begin(), *it ) );
    }
}

/** Add an implicit DateTimeArray; to save memory, just store it untill we need it */
void QMCHelperDateTimeCacheNew::add
    ( DateTimeArrayConstSP dates )
{
    if ( dates.get() ) {
        add( *dates );
    }
}

void QMCHelperDateTimeCacheNew::freeze()
{
    ASSERT( !frozen );
    const size_t esize = entries.size();
    vector<DTSiter> v(entries.begin(), entries.end());
    iterators.swap(v);
    
    IterSet tmp; 
    entries.swap(tmp);
    frozen = true;

    ASSERT(iterators.size() == esize);
    ASSERT(entries.size() == 0);
}

void QMCHelperDateTimeCacheNew::thaw()
{
    ASSERT(frozen);

    const size_t isize = iterators.size();
    // two step process: construct and swap to force memory cleanup
    
    IterSet newEntries(iterators.begin(), iterators.end());
    entries.swap(newEntries);

    vector<DTSiter> empty;
    iterators.swap(empty);
    frozen = false;

    ASSERT(entries.size() == isize);
    ASSERT(iterators.size() == 0);

}

DateTimeArray QMCHelperDateTimeCacheNew::getDates()
{
    if ( ! frozen )
        freeze();
    DateTimeArray res( iterators.size() );
    for ( size_t i = 0; i < iterators.size(); ++i )
        res[ i ] = getDate( i );
    return res;
}

// Just like string::npos, indicates the DateTime is Not in cache OR not in getDates() array
//const size_t QMCHelperDateTimeCacheNew::npos = size_t(-1);

int QMCHelperDateTimeCacheNew::getIdx( const DateTime& date )
{
    if ( !frozen )
        freeze();

    DateTimeSet * cache = globalCache();
    DTSiter loc = cache->find( date );

    if ( loc == cache->end() )   // serious error: Date we didn't see yet
        return npos; // should we disable this possibility ?

    // Find loc in the vector using the fact it is a sorted array
    pair<vector<DTSiter>::iterator, vector<DTSiter>::iterator > range =
        equal_range( iterators.begin(),
                     iterators.end(),
                     loc, SetIterComparator() );

    if ( distance( range.first, range.second ) == 1 )  // check if we found it
        return distance( iterators.begin(), range.first ); // return offset
    else
        return npos;
}

DateTime QMCHelperDateTimeCacheNew::getDate( int idx )  //< return getDates()[idx] from cache
{
    ASSERT( idx != npos );

    if ( !frozen )
        freeze();
    return *( iterators[ idx ] );
}

QMCHelperDateTimeCacheNew::QMCHelperDateTimeCacheNew() :
        frozen( false )
{}

/* Cleanup memory used by iterators. We can always reconstruct them, so we don't care to clear the vector */
/* Should save memory in SimpathI instrument case when all dates are passed via smart pointers and some of the DateTimeCaches are not needed after finalize()
*/

void QMCHelperDateTimeCacheNew::cleanup()
{
    if ( !frozen )
        freeze();
}

#define DUMPSIZE(x) cerr << #x << " p=" << (&(x)) << " " << typeid(*this).name() <<  " size= " << (x.size()) << " * " << sizeof(*(x.begin())) <<  endl

void QMCHelperDateTimeCacheNew::debug()
{
    DUMPSIZE( iterators );
    DUMPSIZE( entries );
//    DUMPSIZE( ( *globalCache() ) );
}

void QMCHelperDateTimeCacheNew::trim(const DateTime& date)
{
    DateTimeSet * cache = globalCache();
    DTSiter loc = cache->insert(cache->begin(), date);

    if (!frozen)
        freeze();
    iterators.erase(upper_bound(iterators.begin(), iterators.end(), loc, SetIterComparator()), iterators.end());
}


DRLIB_END_NAMESPACE
