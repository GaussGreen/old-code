#include "edginc/config.hpp"
#include "edginc/DoubleArrayMarketFactor.hpp"

DRLIB_BEGIN_NAMESPACE

DoubleArrayMarketFactor::~DoubleArrayMarketFactor(){}

/** Constructor */
DoubleArrayMarketFactor::DoubleArrayMarketFactor(const DoubleArray& value):
    value(value) {}

/** Return value */
const DoubleArray& DoubleArrayMarketFactor::getValue() const
{
    return value;
}

DRLIB_END_NAMESPACE
