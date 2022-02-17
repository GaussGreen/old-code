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

#ifndef QMCHelperDateTimeCache_HPP
#define QMCHelperDateTimeCache_HPP

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
class MCARLO_DLL QMCHelperDateTimeCache : public virtual IQMCHelperDateTimeCache
{
private:
        typedef set<DateTime> DateTimeSet;
        typedef DateTimeSet::iterator DTSiter;
        typedef PtrComparator<DTSiter> SetIterComparator;
        typedef set<DTSiter, SetIterComparator> IterSet;
        typedef refCountPtr<DateTimeSet> DateTimeSetSP;

        bool                isValid;   // true iff "iterators" is vecorized directEntries and indirectEntries is empty
        vector<DTSiter>     iterators; // vectorized form of directEntries

        IterSet directEntries; // things submitted directly via DateTimeArray
        vector<DateTimeArrayConstSP> indirectEntries; // things submitted indirectly

        // Global cache will trigger purify/valgrind, but let's worry about it after it works
        static DateTimeSet * globalCache();
        DateTimeArray extractDirect() const;
        DateTimeArray extractIndirect() const;

        void invalidate(void); // called when we add something
        void initialize(void); // populates "iterators" vector
public:
        QMCHelperDateTimeCache();

        virtual void add(const DateTimeArray& dates);
        virtual void add(DateTimeArrayConstSP dates);

        virtual DateTimeArray       getDates();
        virtual DateTime            getDate(int idx); // logically return getDates()[idx], but from global cache
        virtual int                 getIdx(const DateTime& date); // logically it is offset of date inside getDates() or IQMCHelperDateTimeCache::npos
        virtual void                cleanup();
        virtual void                debug();
        virtual void                trim(const DateTime& maxDate); /// more memory tricks
};

DECLARE_REF_COUNT(QMCHelperDateTimeCache);

DRLIB_END_NAMESPACE

#endif
