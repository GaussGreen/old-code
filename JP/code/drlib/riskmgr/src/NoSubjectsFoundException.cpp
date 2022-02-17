/**
 * @file NoSubjectsFoundException.cpp
 */

#include "edginc/config.hpp"
#include "edginc/NoSubjectsFoundException.hpp"

DRLIB_BEGIN_NAMESPACE

NoSubjectsFoundException::NoSubjectsFoundException(const string& message):
    ModelException(message)
{}

NoSubjectsFoundException::NoSubjectsFoundException() {}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
