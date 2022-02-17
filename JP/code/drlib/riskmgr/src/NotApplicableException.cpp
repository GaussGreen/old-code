/**
 * @file NotApplicableException.cpp
 */

#include "edginc/config.hpp"
#include "edginc/NotApplicableException.hpp"

DRLIB_BEGIN_NAMESPACE

NotApplicableException::NotApplicableException(const string& message):
    ModelException(message)
{}

NotApplicableException::NotApplicableException():
    ModelException("Not applicable")
{}

NotApplicableException::~NotApplicableException() throw () {}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
