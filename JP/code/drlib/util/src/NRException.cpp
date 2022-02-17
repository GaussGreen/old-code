//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : NRException.hpp
//
//   Description : Exception related to Numerical Recipes functions
//
//   Author      : Jose Hilera
//
//   Date        : 27 Sep 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/NRException.hpp"

DRLIB_BEGIN_NAMESPACE

/** Instances of this class are thrown when a Numerical Recipe function
 * fails for any reason.
 * There are no new methods or attributes here - just a different type
 * of exception. */

NRException::NRException(const char* errorMsg) : ModelException(errorMsg)
{}


NRException::NRException(const string& errorMsg) : ModelException(errorMsg)
{}


NRException::NRException(const string& routine, const string& errorMsg) : 
    ModelException(routine, errorMsg)
{}


NRException::NRException(const NRException& e) : ModelException(create(e))
{}


NRException::~NRException() throw()
{}


ModelException* NRException::clone() const {
    NRException* e = new NRException(*this);
    return e;
}


/** indicates whether this exception is derived from ModelException -
    used to drive whether this is stored as a 'cause' when a new
    exception is created */
bool NRException::isDerived() const{
    return true;
}


/** Returns the 'cause' exception if it's a NRException
    otherwise returns null. Equivalent to
    dynamic_cast<NRException*>e.getCause() */
NRException* NRException::getInstance(ModelException& e) {
    return dynamic_cast<NRException*>(e.getCause());
}


DRLIB_END_NAMESPACE
