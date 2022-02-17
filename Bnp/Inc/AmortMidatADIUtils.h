// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning(disable : 4786) //"identifier was truncated to '255' characters
                                // in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "crtdbg.h"
#include "math.h"
#include "search.h"
#include "srt_h_all.h"
#include "utmemory.h"

#ifdef _DEBUG
#include "crtdbg.h"
#endif

// macro to swap 2 objects
#define _SWAP_(_object_1_, _object_2_, _object_holder_)                        \
  {                                                                            \
    (_object_holder_) = (_object_1_);                                          \
    (_object_1_) = (_object_2_);                                               \
    (_object_2_) = (_object_holder_);                                          \
  }

#define _INIT_MATRIX_(_matrix_, _row_begin_, _row_end_, _col_begin_,           \
                      _col_end_, _value_, _nr_, _nc_)                          \
  {                                                                            \
    for ((_nr_) = (_row_begin_); (_nr_) < (_row_end_); ++(_nr_)) {             \
      for ((_nc_) = (_col_begin_); (_nc_) < (_col_end_); ++(_nc_)) {           \
        (_matrix_)[(_nr_)][(_nc_)] = (_value_);                                \
      }                                                                        \
    }                                                                          \
  }

#define _SELFMULT_MATRIX_(_matrix_, _row_begin_, _row_end_, _col_begin_,       \
                          _col_end_, _value_, _nr_, _nc_)                      \
  {                                                                            \
    for ((_nr_) = (_row_begin_); (_nr_) < (_row_end_); ++(_nr_))               \
      for ((_nc_) = (_col_begin_); (_nc_) < (_col_end_); ++(_nc_))             \
        (_matrix_)[(_nr_)][(_nc_)] *= (_value_);                               \
  }

#define _SELFADD_MATRIX_(_matrix_, _row_begin_, _row_end_, _col_begin_,        \
                         _col_end_, _value_, _nr_, _nc_)                       \
  {                                                                            \
    for ((_nr_) = (_row_begin_); (_nr_) < (_row_end_); ++(_nr_))               \
      for ((_nc_) = (_col_begin_); (_nc_) < (_col_end_); ++(_nc_))             \
        (_matrix_)[(_nr_)][(_nc_)] += (_value_);                               \
  }

typedef int(__stdcall *PFN_COMPARE_BINARY)(const void * /*pvRight*/,
                                           const void * /*pvLeft */);

// standard interface for unary function
typedef void(__stdcall *PFN_U)(
    const void *, // pIn
    void *,       // pComm      , // common function parameters
    void *        // output
);

// standard interface for binary function
typedef void(__stdcall *PFN_BI)(
    const void *, // pIn1
    const void *, // pIn2
    void *,       // pComm      , // common function parameters
    void *        // output
);
//-------------------------------------------------------------------------------------------------------
//	Description: transform the matrix pOut by calling the binary function
// pfnBI
//
//	Returns : if(n1R)
//						pOut[n_RBegin+i][nC_Begin+j]=pfnBI(pIn1[i]
//    ,pIn2[j]      ,pComm) 					else
//    pOut[n_RBegin+j][nC_Begin+i]=pfnBI(pIn1[i]      ,pIn2[j] ,pComm)
//
//	NB:  pfnBI must take the order of the arguments pIn1 and pIn2
//-------------------------------------------------------------------------------------------------------
void _transform_matrix(
    const void *pIn1,
    unsigned nElem_In1, // # elements of input 1
    unsigned nSize_In1, // size of elements in input1(in bytes)
    const void *pIn2,
    unsigned nElem_In2, // # elements of input 2
    unsigned nSize_In2, // size of elements in input2(in bytes)
    unsigned
        nSize_Out, // size of elements (in bytes) returned by the binary functor
    PFN_BI pfnBI,  // binary functor
    void *pComm,   // common function parameters
    void **pOut,   // output-	pOut[nR_Begin      ,...][nC_Begin      ,...]
    unsigned nR_Begin, unsigned nC_Begin, int n1R);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// function retuns 1 if dX equals machine precision      , 0 otherwise
int is_zero(double dX);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//// exp_integrate_pc first integrates y(x) from x= dLB to x = dUB
//// piece-wise constant      , and then returns the exp of the integration
double exp_integrate_pc(double dLB, // integration lower bound
                        double dUB, // integration upper bound
                        const double *pdX, const double *pdY, int nBegin,
                        int nEnd);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//// exp_minus_integrate_pc first integrates y(x) from x= dLB to x = dUB
//// piece-wise constant      , and then returns the exp of -1. * integration
double exp_minus_integrate_pc(double dLB, // integration lower bound
                              double dUB, // integration upper bound
                              const double *pdX, const double *pdY, int nBegin,
                              int nEnd);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// performs binary search to determine least n s.t. compare(key      ,begin + m)
// <= 0 for all m >= n
const void *lower_bound(const void *pvKey,   // pointer to key
                        const void *pvBegin, // pointer to base element
                        unsigned nElems,     // # elements
                        unsigned nSize,      // size of elements (in bytes)
                        PFN_COMPARE_BINARY pfnCompare // binary compare function
);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// performs binary search to determine n s.t. compare(key      ,begin + m) < 0
// for all m >= n
const void *upper_bound(const void *pvKey,   // pointer to key
                        const void *pvBegin, // pointer to base element
                        unsigned nElems,     // # elements
                        unsigned nSize,      // size of elements (in bytes)
                        PFN_COMPARE_BINARY pfnCompare // binary compare function
);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// populate pOut_Begin with *pKey and return pOut_Begin
const void *_memset(const void *pKey,
                    void *pOut_Begin, // pointer to base element
                    unsigned nElems,  // # elements
                    unsigned nSize    // size of elements (in bytes)
);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
///// returns min(n : pvBegin[n]  > *pvKey      , for 0 <= n <nElems )
///// if no such an index can be found      , return pvBegin+nElems*nSize
const void *
find_min_greater_than(const void *pvKey,   // pointer to key
                      const void *pvBegin, // pointer to base element
                      unsigned nElems,     // # elements
                      unsigned nSize,      // size of elements (in bytes)
                      PFN_COMPARE_BINARY pComp);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// binary compare function
// NB: return value must be:
//	< 0 if right < left;
//	0 if right == left;
//	> 0 if right > left;
int _dcompare(const void *pvL, const void *pvR);

// binary equal function
// NB: return value must be:
//	1 if right != left;
//	0 if right == left;
int _dequal(const void *pvL, const void *pvR);
int _nequal(const void *pvL, const void *pvR);

// most commonly used binaries for double comparison-
int __stdcall dless(const double *pdLeft, const double *pdRight);
int __stdcall lless(const long *pdLeft, const long *pdRight);
int __stdcall dgreater(const double *pdLeft, const double *pdRight);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// look up psubcBegin in pcBegin and return corresponding pvBegin in psubvOut
// NB: nsubElems<=nElems and  pcBegin must contain psubcBegin
void _bmatch(const void *pcBegin,    // linear coordinates
             unsigned ncSize,        // size of coordinates (in bytes)
             const void *pvBegin,    // # elements in values
             unsigned nvSize,        // size of values (in bytes)
             unsigned nElems,        // # elements in the linear coordinates
             const void *psubcBegin, // sub linear coordinates
             unsigned nsubElems,     // # elements in the sub linear coordinates
             void *psubvOut);
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// look up psubcBegin in pcBegin and return corresponding pvBegin in psubvOut
// if elements in psubcBegin are not found in pcBegin      , interpolate
// according to SORT convention
void _bmatch_ex(const void *pcBegin,   // linear coordinates
                unsigned ncSize,       // size of coordinates (in bytes)
                const void *pvBegin,   // # elements in values
                unsigned nvSize,       // size of values (in bytes)
                unsigned nElems,       // # elements in the linear coordinates
                const void *pcInBegin, // sub linear coordinates
                unsigned nInElems, // # elements in the sub linear coordinates
                void *pInvOut);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// combine and then sort pv1Begin and pv2Begin in ascending order
// results returned in pvOut
void _sort2(const void *pv1Begin, // base 1
            unsigned n1Elems,     // # elements in pv1Begin
            const void *pv2Begin, // base 2
            unsigned n2Elems,     // # elements pv2Begin
            unsigned nvSize,      // size of v1/v2 (in bytes)
            void *pvOut);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// combine and then sort pv1Begin      , pv2Begin  and pv3Begin in ascending
// order results returned in pvOut
void _sort3(const void *pv1Begin, // base 1
            unsigned n1Elems,     // # elements in pv1Begin
            const void *pv2Begin, // base 2
            unsigned n2Elems,     // # elements pv2Begin
            const void *pv3Begin, // base 3
            unsigned n3Elems,     // # elements pv2Begin
            unsigned nvSize,      // size of v1/v2/v3 (in bytes)
            void *pvOut);

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+--+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// perform linear interpolation according to the SORT convention : (left open
//     ,righ close] results returned in pv
void sort_interpolate(
    const void *pcBegin,   // linear coordinates
    unsigned ncSize,       // size of coordinates (in bytes)
    const void *pvBegin,   // # elements in values
    unsigned nvSize,       // size of values (in bytes)
    unsigned nElems,       // # elements in the linear coordinates
    const void *pcInBegin, // coordinates to be interpolated
    unsigned ncInElems,    // # elements in the linear coordinates
    void *pv               /// results
);

////////////////////Funtion undocumented
const double *_multiply_vector(const double *pdIn1_Begin,
                               const double *pdIn1_End,
                               const double *pdIn2_Begin, double *pdOut_Begin);

const double *_multiply_vector2(const double *pdIn_Begin,
                                const double *pdIn_End, double dConstant,
                                double *pdOut_Begin);

const double *_add_vector(const double *pdIn1_Begin, const double *pdIn1_End,
                          const double *pdIn2_Begin, double *pdOut_Begin);

const double *_max_dEle(const double *pdBegin, const double *pdEnd);
const int *_max_nEle(const int *pnBegin, const int *pnEnd);

void decomp_sysmetric_tridiagonal(
    /// matrix to be decomposed
    const double *pdDiag, int nDiag_Begin, int nDiag_End,
    const double *pdSubDiag, int nSubDiag_Begin,
    /// result
    double *pdSubDiagL, int nSubDiagL_Begin, double *pdPivot, int nPivot_Begin);

const char *Decode_RefRate(const char *szRefRate, char **pszFeq,
                           char **pszBasis);
double _PayRec(const char *szPayRec);

Err check_null_ptr(void *ptr);
void _free_dvector(double *pd_Begin, int nSize);
void _free_and_zero(const void *);

double _Sum_dVec(const double *pdIn_Begin, const double *pdIn_End);
double _Sum_dMatrix(const double **ppdMatrix, int nR_Begin, int nR_End,
                    int nC_Begin, int nC_End);

Err _fill_df(const char *, long, const long *, const long *, double *);
const double *_fill_tenor(long, const long *, const long *, double *);
const long *_fill_date(long, const double *, const double *, long *);
const char *_fill_swapcashrate_diag(const char *szYC, const long *plStart_Begin,
                                    const long *plStart_End, long lEndDate,
                                    const char *szFreq, const char *szBasis,
                                    const char *szRefRate,
                                    double *pdSwapCashRate_Begin);
int _is_ascending_(const double *, const double *);

const char *_alloc_dvector(const double **ppd_Begin, int nSize);
const char *_alloc_dmatrix(const double ***pppd_Begin, int nR_Begin, int nR_End,
                           int nC_Begin, int nC_End);
void **_alloca_matrix(void **pm, void *pdata, int nR_Begin, int nR_End,
                      int nC_Begin, int nC_End, unsigned nSize);

const char *_free2(const char *szErr, const void *, const void *);

const void *_reverse(void *pvBegin,   // pointer to base element
                     unsigned nElems, // # elements
                     unsigned nSize   // size of elements (in bytes)
);

void memswap(void *pm1, void *pm2, unsigned nSize);

long _ToDate(long lBase, double dT);

double _NormDensity(double dOne_Over_Sqr2Pi, double dX, double dMu,
                    double dVar);

#ifdef _DEBUG
void _print_vector(const double *pdBegin, int nSize);
void _print_nvector(const int *pnBegin, int nSize);
void _print_matrix(const double **ppdMatrix, int nR_Begin, int nR_End,
                   int nC_Begin, int nC_End);
int _non_descending_(const void *pvBegin, // pointer to base element
                     unsigned nElems,     // # elements
                     unsigned nSize       // size of elements (in bytes)
);

#endif
