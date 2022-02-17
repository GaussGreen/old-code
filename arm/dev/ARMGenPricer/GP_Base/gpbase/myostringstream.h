/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: myostringstream.h,v $
 * Revision 1.1  2003/10/08 16:45:55  ebenhamou
 * Initial revision
 *
 * Revision 1.1  2003/09/30 11:05:24  ebenhamou
 * Initial revision
 *
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file myostringstream.h
 *
 *  \brief brief implementation of ostringstream as
 *	unix seems not to support sstream!
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */

 
/*----------------------------------------------------------------------------*/


/*! \class   oStringStream
 *	\brief  the minimum to support ostringstream for unix!
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPBASE_MYOSTRINGSTREAM_H
#define _INGPBASE_MYOSTRINGSTREAM_H

#include "port.h"
#include <string>
CC_USING_NS( std, string )
#include <stdlib.h>

class oStringStream
{
private:
	string itsString;

public:
	oStringStream(): itsString() {};
	string str() const{ return itsString; }

	oStringStream& operator<<( const string& s ) {
		itsString += s;
		return *this;
	}

	oStringStream& operator<<( char** c ) {
		itsString += string( *c );
		return *this;
	}

	oStringStream& operator<<( char* c ) {
		itsString += string( c );
		return *this;
	}

	oStringStream& operator<<( double d ) {
		char buffer[20];
		sprintf( buffer, "%f", d );
		itsString += string( buffer );
		return *this;
	}
	
	oStringStream& operator<<( int i ) {
		char buffer[20];
		sprintf( buffer, "%d", i );
		itsString += string( buffer );
		return *this;
	}

	oStringStream& operator<<( size_t i ) {
		char buffer[20];
		sprintf( buffer, "%d", i );
		itsString += string( buffer );
		return *this;
	}
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
