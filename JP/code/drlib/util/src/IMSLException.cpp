//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IMSLException.cpp
//
//   Description : Exception related to IMSL functions (c.f. NRException)
//
//   Author      : Matthias Arnsdorf
//
//   Date        : March 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IMSLException.hpp"
#include "edginc/Format.hpp"


DRLIB_BEGIN_NAMESPACE

/** Instances of this class are thrown when a IMSL function
 * fails.  */

IMSLException::IMSLException(int IMSLErrorCode,
							 string errorMsg) : ModelException(errorMsg),
												errorCode(IMSLErrorCode)

{
	// specific error msg's
	switch(IMSLErrorCode)
	{
	// imsl_d_quadratic_prog errors
	case IMSL_NO_MORE_PROGRESS : addMsg("IMSL_NO_MORE_PROGRESS"); break ;
	case IMSL_SYSTEM_INCONSISTENT : addMsg("IMSL_SYSTEM_INCONSISTENT"); break ;

	// other errors
	// ADD ERROR MSG'S HERE AS NEEDED
	}
}


IMSLException::IMSLException(const IMSLException& e) : ModelException(create(e)),
errorCode(e.getErrorCode())
{}


IMSLException::~IMSLException() throw()
{}


ModelException* IMSLException::clone() const {
    IMSLException* e = new IMSLException(*this);
    return e;
}


/** indicates whether this exception is derived from ModelException -
    used to drive whether this is stored as a 'cause' when a new
    exception is created */
bool IMSLException::isDerived() const{
    return true;
}


/** Returns the 'cause' exception if it's a IMSLException
    otherwise returns null. Equivalent to
    dynamic_cast<IMSLException*>e.getCause() */
IMSLException* IMSLException::getInstance(ModelException& e) {
    return dynamic_cast<IMSLException*>(e.getCause());
}

/** returns type of IMSL error defined in imslerr.h */
int IMSLException::getErrorCode() const
{
	return errorCode;
}

/** throws IMSL exception with relevant error code and mesg if there is an error.
		Otherwise does nothing */
void IMSLException::throwIMSLExceptionIfError()
{
	int code = imsl_error_code();
	if(code != 0)
	{
		try
		{
			IMSLError::throwExceptionIfError();
		}
		catch (exception& e)
		{
			throw IMSLException(code, e.what());
		}
	}
}


DRLIB_END_NAMESPACE
