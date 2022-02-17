/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: eventviewer.cpp,v $
 * Revision 1.1  2003/10/16 07:52:11  ebenhamou
 * Initial revision
 *
 */


/*! \file eventviewer.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/eventviewer.h"
#include "gpbase/singleton.h"

CC_BEGIN_NAMESPACE( ARM )

/// by default verbose mode is on!
bool ARM_EventViewerImp::VerboseIsOn	= true;

/// comment this out to remove debug information!
/// bool ARM_EventViewerImp::DebugIsOn		= true;
bool ARM_EventViewerImp::DebugIsOn		= false;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_EventViewerImp
///	Routine: toString
///	Returns: 
///	Action : stringify
/////////////////////////////////////////////////////////////////
string ARM_EventViewerImp::toString() const
{
	string message;
	if( ARM_EventViewerImp::DebugIsOn  )
	{
		message += "Debug   mode : On\n";
		message += itsDebugMessage;
		message += "\n\n\n";
	}

	message += "Verbose mode : ";
	message += ARM_EventViewerImp::VerboseIsOn ? "On" : "Off";
	message += "\n" + itsMessage;
	return message;
}



/// Creation of the singleton the event viewer!
ARM_SingletonHolder<ARM_EventViewerImp> ARM_TheEventViewer;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

