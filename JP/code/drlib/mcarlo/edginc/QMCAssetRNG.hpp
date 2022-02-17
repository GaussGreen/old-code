
//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QMCAssetRNG.hpp
//
//   Description : simple implementation of per-asset random numbers
//
//----------------------------------------------------------------------------


#ifndef QMCASSETRNG
#define QMCASSETRNG

#include "edginc/IQMCAssetRNG.hpp"
#include "edginc/DoubleMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

class MCARLO_DLL QMCAssetRNG : public virtual IQMCAssetRNG
{
    const DoubleMatrix& m;
    size_t baseIdx;
    IUniformRNGGenSP uni;

public:

    /* Note that we need to know the offset (randomIdx) in the global matrix of correlated RNGs */
    QMCAssetRNG(const DoubleMatrix& _m, size_t randomIdx, IUniformRNGGenSP gen) :
        m(_m),
        baseIdx(randomIdx), // at some point it will be assetIdx, but for now we to know line in the matrix
        uni(gen)
    {
    }

    virtual const double * getCorrelatedRandoms(size_t factor = 0) {
        return & (m[baseIdx+factor][0]);
    }
    virtual IUniformRNGGenSP  getExtra() { return uni;}      // get uniform numbers we can use as we pleased

};

DRLIB_END_NAMESPACE

#endif
