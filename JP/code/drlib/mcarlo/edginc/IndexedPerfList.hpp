//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IndexedPerfList.hpp
//
//   Description : Class to sort doubles whilst providing access to which
//                 number ended up where
//
//   Author      : Mark A Robson
//
//   Date        : 29 November 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_INDEXED_PERF_LIST_HPP
#define EDR_INDEXED_PERF_LIST_HPP

#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
/** Class to sort doubles whilst providing access to which number
    ended up where */
class MCARLO_DLL IndexedPerfList{
public:
    class MCARLO_DLL Cmpt{
    public:
        int    assetIdx;
        double perf;
    };

    /** Creates IndexedPerfList holding numPerfs Cmpts. The assetIdx are
        initialised */
    static IndexedPerfList* create(int numPerfs);

    /** Deletes components */
    ~IndexedPerfList();

    /** Returns n'th component. Inline for performance */
    Cmpt& operator[](int index){
        return (*(cmpts[(size_t) index]));
    }

    /** Returns n'th component. Inline for performance */
    const Cmpt& operator[](int index) const{
        return (*(cmpts[(size_t) index]));
    }

    /** Returns number of components */
    int size() const;

    /** returns true if size == 0 */
    bool empty() const;

    /** populate existing array of asset performances from a DoubleArray */
    void populate(const DoubleArray&  perfs);

    /** Sorts assets using perf field of components. Highest
        performance first */
    void sortByPerf();

    /** Calculates weighted sum using current order of IndexedPerf */
    double weightedSum(const DoubleArray&  weights) const;

    /** Store smallest absolute difference between each asset and its
        nearest neighbour (in terms of performance) whose performance is
        actually different (assets with equal performance are judged to be
        equally ranked) */
    void storeNearestPerfDiff(DoubleArray&  perfDiff) const;
private:
    // these two not implemented - have to watch out for memory issues
    IndexedPerfList(const IndexedPerfList& rhs);
    IndexedPerfList& operator=(const IndexedPerfList& rhs);

    IndexedPerfList(int numPerfs);

    // fields
    vector<Cmpt*> cmpts;

    static bool sortMethod(const Cmpt* cmpt1, const Cmpt* cmpt2);
    friend class NoFriends; // stop the compiler moaning
};

typedef refCountPtr<IndexedPerfList> IndexedPerfListSP;

DRLIB_END_NAMESPACE
#endif
