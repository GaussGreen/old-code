/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file warning.cpp
 *  \brief object for any warning
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2005
 */

#include "gpbase/warning.h"

#define VC_EXTRALEAN		// Exclude rarely-used stuff from Windows headers
#include <Windows.h>
CC_BEGIN_NAMESPACE( ARM )

bool ARM_Warning::isPopUpOn = false;

const int ARM_Warning::MaxMessageSize = 1000;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_Warning
///	Routine: PopUp
///	Action : does a popup if the flag is on
/////////////////////////////////////////////////////////////////
void ARM_Warning::PopUp() const
{
	if( ARM_Warning::isPopUpOn )
		MessageBox( NULL, itsMessage.c_str(), "Warning", MB_OK );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_Warning
///	Routine: PopUp
///	Action : add a message
/////////////////////////////////////////////////////////////////
void ARM_Warning::AddToMessage( const string& message ) 
{ 
	if (itsMessage.size() < MaxMessageSize)
	{
		itsMessage += message; 
		PopUp(); 
	}
}

CC_END_NAMESPACE()