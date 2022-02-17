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

#ifndef QLIB_NREXCEPTION_HPP
#define QLIB_NREXCEPTION_HPP

#include "edginc/ModelException.hpp"
#include "edginc/DateTime.hpp"


DRLIB_BEGIN_NAMESPACE


/** Instances of this class are thrown when a Numerical Recipe function
 * fails for any reason.
 * There are no new methods or attributes here - just a different type
 * of exception. */

class UTIL_DLL NRException : public ModelException {
public:
    NRException(const char* errorMsg);
    NRException(const string& errorMsg);
    NRException(const string& routine, const string& errorMsg);
    NRException(const NRException& e);

    virtual ~NRException() throw ();
    
    /** creates a [deep] copy of the exception */
    virtual ModelException* clone() const;
    
    /** indicates whether this exception is derived from ModelException -
        used to drive whether this is stored as a 'cause' when a new
        exception is created */
    virtual bool isDerived() const;
    
    /** Returns the 'cause' exception if it's a NRException
        otherwise returns null. Equivalent to
        dynamic_cast<NRException*>e.getCause() */
    static NRException* getInstance(ModelException& e);
};


DRLIB_END_NAMESPACE

#endif
