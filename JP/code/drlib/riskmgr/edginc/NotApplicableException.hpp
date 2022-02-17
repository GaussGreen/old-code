/**
 * @file NotApplicableException.hpp
 */

#ifndef QLIB_NotApplicableException_H
#define QLIB_NotApplicableException_H

#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Thrown from RiskQuantity::notApplicable(), caught by
 * NamedRiskQuantity::storeResult() to get a "not applicable"
 * into the Results.
 */

class RISKMGR_DLL NotApplicableException: public ModelException {

public:

    NotApplicableException(const string& message);
    NotApplicableException();
    ~NotApplicableException() throw ();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
