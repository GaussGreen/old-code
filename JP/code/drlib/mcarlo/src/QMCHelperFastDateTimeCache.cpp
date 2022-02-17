//
// C++ Implementation: QMCHelperFastDateTimeCache
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "edginc/config.hpp"
#include "edginc/QMCHelperFastDateTimeCache.hpp"

DRLIB_BEGIN_NAMESPACE

void QMCHelperFastDateTimeCache::normalize( void )
{
    if ( !isSorted )
        DateTime::doSortUniq( dates );
    isSorted = true;
}

QMCHelperFastDateTimeCache::QMCHelperFastDateTimeCache() : dates(), isSorted( false )
{
}

void QMCHelperFastDateTimeCache::add
    ( const DateTimeArray& _dates )
{
    if ( _dates.size() == 0 )
        return ;
    dates.insert( dates.end(), _dates.begin(), _dates.end() );
    isSorted = false;
    if ( _dates.size() > 1 )
        normalize();
}

DateTimeArray QMCHelperFastDateTimeCache::getDates()
{
    normalize();
    return dates;
}

DateTime QMCHelperFastDateTimeCache::getDate( int idx )
{
    normalize();
    return dates[ idx ];
}

int QMCHelperFastDateTimeCache::getIdx( const DateTime& date )   // logically it is offset of date inside getDates() or npos
{
    normalize();
    return date.find( dates );
}

DRLIB_END_NAMESPACE
