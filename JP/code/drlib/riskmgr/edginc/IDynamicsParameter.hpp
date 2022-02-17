/**
 * @file IDynamicsParameter.hpp
 */

#ifndef QLIB_IDynamicsParameter_H
#define QLIB_IDynamicsParameter_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IDynamicsParameter)

/**
 * Tag interface for "vol" objects, like VolSVJ or SRMEQ::Vol, which are really
 * parameters for parametric pdf models
 *
 * This is for RiskMapping: prices computed using models which are
 * parameterised in terms of numbers like SRMEQ::Vol::smileA1 or
 * VolSVJ::meanReversRate look as if they are insensitive to good old fashioned
 * spot/vol/... --- delta and vega come out zero (or, worse, nonzero but
 * misleading).  RiskMapping enables hypotheses about spot/vol/... to be
 * re-expressed in terms of tweaks to opaque dynamics parameters so as to get
 * more meaningful Greeks.
 */

class RISKMGR_DLL IDynamicsParameter: public virtual IObject {

public:

    static CClassConstSP const TYPE;

    IDynamicsParameter();
    ~IDynamicsParameter();
};

DRLIB_END_NAMESPACE

#endif
