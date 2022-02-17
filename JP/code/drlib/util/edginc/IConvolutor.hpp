//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ICONVOLUTOR_HPP
#define QLIB_ICONVOLUTOR_HPP

#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/IDistribution1D.hpp"
//#include "edginc/ICreditLossConfig.hpp"

DRLIB_BEGIN_NAMESPACE

class DateTime;
FORWARD_DECLARE(ICreditLossConfig);

/**
 * Interface for objects that can perform a convolution algorithm on
 * distributions.
 * We authorise 2 flavours to cope with possible optimisations
 * (eg: Saddle Point, Fourier...):
 * - "convolute"            : "pure" convolution algorithm that returns
 *                            the whole convoluted distribution
 * - "convoluteAndIntegrate": convolution and integration over a "payoff" i.e.
 *                            integral(convolute(D1,D2,...)(x)*payoff(x)dx)
 * */
class UTIL_DLL IConvolutor: public virtual IObject {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** "Pure" convolution algorithm that returns the whole convoluted distribution */
    virtual IDistribution1DConstSP convolute(
        IDistribution1DArrayConstSP distributions) const = 0;

    /**
     * Convolution and integration over a "payoff" i.e.
     * integral(convolute(D1,D2,...)(x)*payoff(x)dx).
     * Need to pass a "timepoint" to allow time-dependent payoffs.
     * */
    virtual double convoluteAndIntegrate(
        IDistribution1DArrayConstSP distributions,
        ICreditLossConfigConstSP lossConfig,
        const DateTime& timepoint) const = 0;
};

DECLARE(IConvolutor);

DRLIB_END_NAMESPACE

#endif /*QLIB_ICONVOLUTOR_HPP*/
