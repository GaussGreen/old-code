//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OneToOneMapping.hpp
//
//   Date        : 20 June 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Range.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

#ifndef EDR_ONE_TO_ONE_MAPPING_HPP
#define EDR_ONE_TO_ONE_MAPPING_HPP

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_REF_COUNT(OneToOneMapping);

class UTIL_DLL OneToOneMapping{
public:
    virtual double operator()(double x) const = 0;    // x in R
    virtual double inverse(double y) const = 0;       // y in (lower, upper)
    // first derivative
    virtual double derivative(double x) const = 0;    // x in R

    virtual ~OneToOneMapping(){}
    static OneToOneMappingConstSP create(const Range& range);
};

DRLIB_END_NAMESPACE

#endif




