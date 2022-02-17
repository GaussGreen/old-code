//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCHelperCachingTimeLogic
//
//   Description :
//
//   Date        : Wed Apr 12 2006
//
//   Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>
//
//----------------------------------------------------------------------------

#ifndef QMCHelperCachingTimeLogic_HPP
#define QMCHelperCachingTimeLogic_HPP

#include "edginc/IQMCHelperTimeLogic.hpp" 
//#include "edginc/QMCHelperDateTimeCache.hpp"
//#include "edginc/QMCHelperFastDateTimeCache.hpp"
//#include "edginc/IQMCDiffusibleAssetBase.hpp"
#include "edginc/TemplateIdx.hpp"
#include "edginc/IQMCHelperDateTimeCache.hpp"


DRLIB_BEGIN_NAMESPACE

////////////////////////////////////////////////////////////////////////

/*
Memory requirements analysis for SimpathI
1. we may have up to 10,000 CR
2. Each asset will have essentially 2 copies of timeline (sdfDates and esdfRequestedDates) plus all forward dates (which is like 5*diffusion dates)
3. With Aggregated SVGen, all Gen will get all dates, meaning we keep ~10^4assets * 10^4 dates Given that each DateTime takes ~24bytes (40bytes on 64bit), that quickly runs into more memory than we can afford.
 
Solution:
- Keep dates passed via DateTimeArrayConstSP separately from the "bare" dates
- Use level of inderection when keeping dates:
 - keep _all_ (accross all assets) encountered dates in a set (should be less than ~365*100 as time is not used)
 - point to dates via iterators in this set (which stays valid upod addition)
 - encode dateTimeArray via set of iterators (where each iterator points to the caching set)
 
 TODO: how expensive it is to keep a set of iterators? (potentially overhead is 2 pointers + malloc (i.e. upto 16bytes of overhead), so may be we should use vector and keep it sorted or something even more smarter)
 
 
*/

// Cache all DateTime seen in a set; keep
class QMCHelperCachingTimeLogic : public IQMCHelperTimeLogic
{
private:

    // sets of iterators for dates submitted in arrays
    IQMCHelperDateTimeCacheSP df;
    IQMCHelperDateTimeCacheSP edfReq;
    IQMCHelperDateTimeCacheSP edfFwd;

public:
    QMCHelperCachingTimeLogic( IQMCHelperDateTimeCacheSP _df,
                               IQMCHelperDateTimeCacheSP _edfReq,
                               IQMCHelperDateTimeCacheSP _edfFwd ) :
            df( _df ), edfReq( _edfReq ), edfFwd( _edfFwd )
    {
    }
    
    DateTimeArray getDFDates() const
    {
        return df->getDates();
    }

    DateTimeArray getReqEDFDates() const
    {
        return edfReq->getDates();
    }

    DateTimeArray getFwdEDFDates() const
    {
        return edfFwd->getDates();
    }

    virtual void addAggregatedDates(
        const DateTimeArray& dfDates,
        const DateTimeArray& edfReqDates,
        const DateTimeArray& edfFwdDates )
    {
        df->add
        ( dfDates );
        edfReq->add
        ( edfReqDates );
        edfFwd->add
        ( edfFwdDates );
    }

    virtual void addAggregatedDates(
        DateTimeArrayConstSP dfDates,
        DateTimeArrayConstSP edfReqDates,
        DateTimeArrayConstSP edfFwdDates )
    {
        df->add
        ( dfDates );
        edfReq->add
        ( edfReqDates );
        edfFwd->add
        ( edfFwdDates );
    }

    virtual SpotIdx getDFIdx( const DateTime& date )
    { // returns Idx of a forward date
        return df->getIdx( date );
    }

    virtual SpotIdx getReqEDFIdx( const DateTime& date )
    { // returns Idx of a forward date
        return edfReq->getIdx( date );
    }

    // Given a date, return its FwdIdx in the set of Forward dates
    // Indices are invalidated every time
    virtual FwdIdx getFwdEDFIdx( const DateTime& date )
    { // returns Idx of a forward date
        return edfFwd->getIdx( date );
    }

    virtual SpotIdx getReqEDFIdx( FwdIdx fwdIdx )
    {
        return getReqEDFIdx( edfFwd->getDate( fwdIdx ) );
    }
    virtual DateTime getReqDate( SpotIdx idx )  ///< return date such that getReqEDFIdx(date) == idx
    {
        return edfReq->getDate( idx );
    }

    virtual void cleanup()
    {
        df->cleanup();
        edfReq->cleanup();
        edfFwd->cleanup();
    }
    
    virtual void debug()
    {
        df->debug();
        edfReq->debug();
        edfFwd->debug();
    }
    
    virtual void trim(const DateTime& maxDiffMat, const DateTime& maxCurveMat) {
        df->trim(maxDiffMat);
        edfReq->trim(maxDiffMat);
        edfFwd->trim(maxCurveMat);
    }
};

DECLARE_REF_COUNT(QMCHelperCachingTimeLogic);

DRLIB_END_NAMESPACE

#endif
