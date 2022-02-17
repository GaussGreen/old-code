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
//
//----------------------------------------------------------------------------

#ifndef QMCHelperDateTimeCacheNew_HPP
#define QMCHelperDateTimeCacheNew_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/TemplateIdx.hpp"
#include "edginc/IQMCHelperDateTimeCache.hpp"

#include <set>
#include <vector>

DRLIB_BEGIN_NAMESPACE

// compare iterators (or pointers) via underlying objects
// requirement : 1. T has operator *() 2. (*T) has operator <() (i.e. (*a) < (*b) compiles

template <class T>
struct  PtrComparator : public binary_function<T, T, bool>
{
    bool operator() (const T& l, const T& r) const
    {
        return (*l) < (*r);
    }
};


/** This class allows one to very succintly represent DateTimeArrays.
    Remember that now DateTime is quite a complex object that takes much more memory than needed.
    This class logically represents merge of several DateTimeArrays given explicitly or via smartPtrs. When passed explicitly, dates are accumulated in a global cache and only iterators to this cache are stored; when passed via SP, only SP is copied and no other memory is wasted.

    We support two operations: to retrieve all dates (see getDates()) that is to merge all that was passed to add() so far.
    The second one, getIdx, is mostly needed for StateVariable framework where we need to determine indices of Dates inside diffusionDates without spending too much memory (think 10,000 CreditSpreads with 10K time points)
*/
class MCARLO_DLL QMCHelperDateTimeCacheNew
    : public virtual IQMCHelperDateTimeCache
{
private:
        typedef set<DateTime> DateTimeSet;
        typedef DateTimeSet::iterator DTSiter;
        typedef PtrComparator<DTSiter> SetIterComparator;
        typedef set<DTSiter, SetIterComparator> IterSet;
        typedef refCountPtr<DateTimeSet> DateTimeSetSP;

        bool                frozen;   // while "false" we update "entries"; when "true" entries is moved into iterators.

        IterSet             entries; // things we saw
        vector<DTSiter>     iterators; // vectorized form of directEntries
        
    static DateTimeSet * globalCache(); // singleton of all dates we saw
        void                freeze(); // calculate "iterators" and release "entries"
        void                thaw(); // moves "iterators" back to "entries" and clears "frozen" flag

public:
        QMCHelperDateTimeCacheNew();

        virtual void add(const DateTimeArray& dates);
        virtual void add(DateTimeArrayConstSP dates);

        virtual DateTimeArray       getDates();
        virtual DateTime            getDate(int idx); // logically return getDates()[idx], but from global cache
        virtual int                 getIdx(const DateTime& date); // logically it is offset of date inside getDates() or IQMCHelperDateTimeCache::npos
        virtual void                cleanup();
        virtual void                debug();
        virtual void                trim(const DateTime& date); // remove entries >date
};

DECLARE_REF_COUNT(QMCHelperDateTimeCacheNew);

DRLIB_END_NAMESPACE

#endif
