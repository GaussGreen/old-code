//
// C++ Interface: TemplateIdx
//
// Description: Aux template to create objects that behave like integers, but that cannot be assigned to each other
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef EDR_TEMPLATEIDX_HPP
#define EDR_TEMPLATEIDX_HPP

#include <edginc/config.hpp>
#include <edginc/smartPtr.hpp>

DRLIB_BEGIN_NAMESPACE


/** this class for strong separation of spot value indexing from
    expected value indexing */
#if !defined (_MSC_VER) || (_MSC_VER >= 1300)

template <int n> class TemplatedIdx
{
    public:
        TemplatedIdx() : idex(0) {}
        TemplatedIdx(int i) : idex(i) {}
        operator int& () { return idex; }

// this part makes certain that templates of two different kind cannot
// be automatically converted to each other
        TemplatedIdx(const TemplatedIdx<n>& b) : idex(b.idex){}
        TemplatedIdx<n>& operator=(const TemplatedIdx<n>& b){idex = b.idex; return *this;}
        static const int npos = -1; // npos means entry was not found when we resolve for index
    private:
        template <int nn> TemplatedIdx(const TemplatedIdx<nn>&);
        template <int nn> TemplatedIdx<n>& operator=(const TemplatedIdx<nn>&);
        int idex;
};

// These are templatized classes for distinguishing indexes of spot
// and forward measurements. they might and should be used in place of int
// for index safety.

// examples:
//  SpotIdx idx_spot = 0;    //<-- ok, since SpotIdx is constructed from int
//  int idx_int = idx_spot.toInt();  //<-- ok, but needs to be explicit
//  FwdIdx  idx_fwd  = idx_spot;//<-- NOT OK, as indexes of the same date
                                // can be different for spot and fwd arrays

typedef TemplatedIdx<0>  SpotIdx;
typedef TemplatedIdx<1>  FwdIdx;

#else // we are in the VC6 world
// unfortunately strong type safety here is impossible for VC6 compiler due to its issues

typedef int  SpotIdx;
typedef int  FwdIdx;
// class SpotIdx {
//     int idex;
// public:
//     SpotIdx(int i=0) : idex(i) {}
//     enum{npos = -1};
//     operator int& () { return idex; }
// };
//
// typedef SpotIdx FwdIdx;

#endif

typedef vector<SpotIdx> SpotIdxArray;
typedef refCountPtr<SpotIdxArray> SpotIdxArraySP;

typedef vector<FwdIdx>  FwdIdxArray;
typedef refCountPtr<FwdIdxArray>   FwdIdxArraySP;

DRLIB_END_NAMESPACE

#endif

