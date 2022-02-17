//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   Date        : 22-Sep-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_DOUBLEARRAYMARKETFACTOR_HPP
#define QLIB_DOUBLEARRAYMARKETFACTOR_HPP

#include "edginc/IMarketFactorValue.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Realisation of the market factor as a "point" in a multi-dimensional space
 * */
class CONVOLUTION_DLL DoubleArrayMarketFactor:
    public IMarketFactorValue // NB: no virtual inheritance to allow static_cast
{
public:

    virtual ~DoubleArrayMarketFactor();

    /** Constructor */
    DoubleArrayMarketFactor(const DoubleArray& value);

    /** Return value */
    const DoubleArray& getValue() const;

private:
    const DoubleArray& value;
};

DECLARE(DoubleArrayMarketFactor);

DRLIB_END_NAMESPACE

#endif /*QLIB_DOUBLEARRAYMARKETFACTOR_HPP*/
