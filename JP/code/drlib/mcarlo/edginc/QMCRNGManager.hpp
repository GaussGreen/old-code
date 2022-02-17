//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : QMCRNGManager.hpp
//
//   Description : simple implementation of per-asset random numbers
//
//----------------------------------------------------------------------------


#ifndef QMC_RNG_MANAGER
#define QMC_RNG_MANAGER

#include "edginc/IQMCRNGManager.hpp"
#include "edginc/IQMCAssetRNG.hpp"
#include "edginc/IRNGGenerator.h" // IUniformRNGGenSP

#include "edginc/DoubleMatrix.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(MCRandom);

class MCARLO_DLL QMCRNGManager : public virtual IQMCRNGManager
{
    const DoubleMatrix * randoms;
    MCRandomSP randomGen;
    ISuperCubeRNGGenSP uniGen; // extra source of randomness per asset; FIXME: implement one stream per asset and another substream per Path

    public:
        QMCRNGManager(MCRandomSP gen, size_t assetsNum);
        virtual void seekToPath(size_t pathIdx);
        virtual IQMCAssetRNGSP  getAssetRNG(size_t randomIdx);
        virtual const double *  getCorrelatedRandoms(size_t assetIndex, size_t factor);
        virtual IUniformRNGGenSP    getSharedGen();
};

DECLARE_REF_COUNT(QMCRNGManager);

DRLIB_END_NAMESPACE

#endif
