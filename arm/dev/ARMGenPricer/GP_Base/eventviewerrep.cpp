/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: eventviewerrep.cpp,v $
 * Revision 1.1  2003/10/16 07:52:11  ebenhamou
 * Initial revision
 *
 */


/*! \file eventviewerrep.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpbase/eventviewerrep.h"
#include "gpbase/singleton.h"

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
///	Class  : ARM_EventViewerRep
///	Routine: constructor,copy constructor,assignment operator and destructor
///	Returns: 
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////

ARM_EventViewerRep::ARM_EventViewerRep( )
:	ARM_RootObject() 
{
	CC_ARM_SETNAME(ARM_EVENT_VIEWER); 
}

ARM_EventViewerRep::ARM_EventViewerRep( const ARM_EventViewerRep& rhs ) 
:	ARM_RootObject( rhs )
{}

ARM_EventViewerRep& ARM_EventViewerRep::operator=( const ARM_EventViewerRep& rhs )
{
	if( this != &rhs )
		ARM_RootObject::operator=(rhs);
	return *this;
}

ARM_EventViewerRep::~ARM_EventViewerRep()
{}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_EventViewerRep
///	Routine: Clone
///	Returns: 
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_EventViewerRep::Clone() const
{
	return new ARM_EventViewerRep(*this);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_EventViewerRep
///	Routine: toString
///	Returns: 
///	Action : standard action for derived object
/////////////////////////////////////////////////////////////////
string ARM_EventViewerRep::toString(const string& indent, const string& nextIndent) const
{
	return ARM_TheEventViewer.Instance()->toString();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_EventViewerRep
///	Routine: ResetMessage, AddToMessage
///	Returns: 
///	Action : reset message remove message while add to message appends
///				the new message to the current message
/////////////////////////////////////////////////////////////////
void ARM_EventViewerRep::ResetMessage() const
{
	ARM_TheEventViewer.Instance()->ResetMessage();
}

void ARM_EventViewerRep::AddToMessage( const string& msg ) const
{
	ARM_TheEventViewer.Instance()->AddToMessage(msg);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

