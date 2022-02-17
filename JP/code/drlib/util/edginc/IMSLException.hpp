//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IMSLException.hpp
//
//   Description : Exception related to IMSL functions (c.f. NRException)
//
//   Author      : Matthias Arnsdorf
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IMSL_EXCEPTION_HPP
#define QLIB_IMSL_EXCEPTION_HPP

#include "edginc/ModelException.hpp"
#include "edginc/imslerror.hpp"


DRLIB_BEGIN_NAMESPACE


/** Instances of this class are thrown when a IMSL function fails. */
class UTIL_DLL IMSLException : public ModelException {
public:
    
	IMSLException(int IMSLErrorCode, string errorMsg);
    IMSLException(const IMSLException& e);

    virtual ~IMSLException() throw ();
    
    /** creates a [deep] copy of the exception */
    virtual ModelException* clone() const;
    
    /** indicates whether this exception is derived from ModelException -
        used to drive whether this is stored as a 'cause' when a new
        exception is created */
    virtual bool isDerived() const;
    
    /** Returns the 'cause' exception if it's a IMSLException
        otherwise returns null. Equivalent to
        dynamic_cast<IMSLException*>e.getCause() */
    static IMSLException* getInstance(ModelException& e);

	/** returns type of IMSL error defined in imslerr.h */
	int getErrorCode() const;


	/** Throws IMSL exception with relevant error code if there is an error.
		Otherwise does nothing. 
		Use this immediately after an IMSL function call.
		Similar to IMSLError::throwExceptionIfError except that throw IMSLException
		instead of ModelException.
	*/
	static void throwIMSLExceptionIfError();

	

private:

	// FIELDS
	int errorCode;
	

};


DRLIB_END_NAMESPACE

#endif

