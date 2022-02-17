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

#ifndef QMCHelperFastDateTimeCache_HPP
#define QMCHelperFastDateTimeCache_HPP

//#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp" 
//#include "edginc/IQMCDiffusibleAssetBase.hpp"
#include "edginc/IQMCHelperDateTimeCache.hpp"

DRLIB_BEGIN_NAMESPACE



/** This class allows one to very succintly represent DateTimeArrays.
    Remember that now DateTime is quite a complex object that takes much more memory than needed.
    This class logically represents merge of several DateTimeArrays given explicitly or via smartPtrs. When passed explicitly, dates are accumulated in a global cache and only iterators to this cache are stored; when passed via SP, only SP is copied and no other memory is wasted.
 
    We support two operations: to retrieve all dates (see getDates()) that is to merge all that was passed to add() so far.
    The second one, getIdx, is mostly needed for StateVariable framework where we need to determine indices of Dates inside diffusionDates without spending too much memory (think 10,000 CreditSpreads with 10K time points)
 */
class QMCHelperFastDateTimeCache : public virtual IQMCHelperDateTimeCache
{
    DateTimeArray dates;
    bool isSorted;
    void normalize( void );
public:
    QMCHelperFastDateTimeCache();

    virtual void add( const DateTimeArray& _dates );

    virtual DateTimeArray getDates();
    virtual DateTime getDate( int idx );
    virtual int getIdx( const DateTime& date );
};

DECLARE( QMCHelperFastDateTimeCache );

DRLIB_END_NAMESPACE

#endif
