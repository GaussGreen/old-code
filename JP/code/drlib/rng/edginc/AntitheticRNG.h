#ifndef Antithetic_RNG_H
#define Antithetic_RNG_H

#include "edginc/coreConfig.hpp"
#include "edginc/IRNG.h"

CORE_BEGIN_NAMESPACE

// AntitheticRNG transformation or RNGGen into RNG
class AntitheticRNG;
DECLARESP(AntitheticRNG);

class RNG_DLL AntitheticRNG : public virtual IUniformRNG
{
    double anti(double x) const
    {
        return 1.0 - x;
    }
    IUniformRNGGenSP rng;
    AntitheticRNG(IUniformRNGGenSP _rng) : rng(_rng)
    {}
public:
    virtual double fetch()
    {
        return anti(rng->fetch());
    }
    virtual IRNGGeneratorSP getGenerator() const
    {
        return rng;
    }
    IRNGSP clone() const
    {
        return AntitheticRNGSP(new AntitheticRNG(DYNAMIC_POINTER_CAST<IUniformRNGGen> (rng->clone())));
    }
    static AntitheticRNGSP create(IUniformRNGGenSP _rng)
    {
        return AntitheticRNGSP(new AntitheticRNG(_rng));
    }

};

CORE_END_NAMESPACE
#endif
