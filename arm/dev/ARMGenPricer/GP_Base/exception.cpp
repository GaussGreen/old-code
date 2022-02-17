/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file exception.cpp
 *
 *  \brief file for the exception
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#if defined(__USE_BASE_LIBRARY)

#include "gpbase/exception.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Exception 
///	Routine: ARM_Exception
///	Returns: 
///	Action : Constructor,Copy constructor, assignment operator
////////////////////////////////////////////////////

ARM_Exception::ARM_Exception(long line, const char* file, long code, const string& message )
:	CC_NS(std,exception)(), 
	itsLine(line), 
	itsFile(file), 
	itsCode(code), 
	itsMessage(message)
{}

ARM_Exception::ARM_Exception(const ARM_Exception& rhs)
:	CC_NS(std,exception)(rhs), 
	itsLine(rhs.itsLine), 
	itsFile(rhs.itsFile), 
	itsCode(rhs.itsCode), 
	itsMessage(rhs.itsMessage)
{}

ARM_Exception& ARM_Exception::operator =(const ARM_Exception& rhs)
{
	if( this!= &rhs )
	{
		exception::operator =(rhs), 
		itsLine		= rhs.itsLine;
		itsFile		= rhs.itsFile;
		itsCode		= rhs.itsCode;
		itsMessage	= rhs.itsMessage;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Exception 
///	Routine: GetErrorMessage
///	Returns: 
///	Action : writes the error message based on member data
////////////////////////////////////////////////////

void ARM_Exception::GetErrorMessage(char* msg) const
{
	sprintf(msg, " ERROR: %s, FILE: %s, LINE: %ld, ERROR CODE: %ld", 
		itsMessage.c_str(), itsFile, itsLine, itsCode );
}

////////////////////////////////////////////////////
///	Class  : ARM_Exception 
///	Routine: GetErrorMessage
///	Returns: 
///	Action : return the error message based on member data
////////////////////////////////////////////////////

const char* ARM_Exception::what() const throw()
{
	char* msg = new char[500];
	GetErrorMessage(msg);
	return msg;
}

CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

