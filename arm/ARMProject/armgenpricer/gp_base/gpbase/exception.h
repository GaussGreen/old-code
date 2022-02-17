/*!
 *
 * Copyright (c) IXIS CI January 2005 Paris
 *
 * \file exception.h
 *
 *  \brief file for the exception object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPBASE_EXCEPTION_H
#define _INGPBASE_EXCEPTION_H

#include "gpbase/port.h"

#include <string>
CC_USING_NS(std,string)
#include <exception>

CC_BEGIN_NAMESPACE( ARM )

class ARM_Exception : public CC_NS(std,exception)
{
private:
	long		itsLine;
	const char*	itsFile;
	long		itsCode;
	string		itsMessage;
public:
	ARM_Exception(long line, const char* file, long code, const string& message);
	ARM_Exception(const ARM_Exception& rhs);
	ARM_Exception& operator=(const ARM_Exception& rhs);
	virtual ~ARM_Exception(){};
	virtual const char* what() const throw();

	/// ARM support
	void GetErrorMessage(char* msg) const;

	/// name of the file where it puts full message in case of truncation
	static char* itsLongMessageLog;
};

/// Macro to avoid typing lengthous message!
#define ARM_THROW( errorMessageNb, errorMessage ) \
	throw ARM_Exception(__LINE__, __FILE__, errorMessageNb, errorMessage );

typedef ARM_Exception Exception;

CC_END_NAMESPACE()


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

