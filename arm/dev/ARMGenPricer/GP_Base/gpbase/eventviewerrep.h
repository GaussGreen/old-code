/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: eventviewerrep.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file eventviewerrep.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPBASE_EVENTVIEWERREP_H
#define _INGPBASE_EVENTVIEWERREP_H

#include "port.h"
#include "eventviewer.h"
#include "rootobject.h"


CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
/// \struct ARM_EventViewerRep (Rep stands for representant!)
/// this is only a representant of the event viewer which is a real singleton!
/////////////////////////////////////////////////////////////////
struct ARM_EventViewerRep : public ARM_RootObject
{
	/// constructor,copy constructor,assignment operator and destructor
	ARM_EventViewerRep( );
	ARM_EventViewerRep( const ARM_EventViewerRep& rhs );
	ARM_EventViewerRep& operator=( const ARM_EventViewerRep& rhs );
	virtual ~ARM_EventViewerRep();

	/// standard ARM Object support
	/// because of the clone implementation in terms of the copy constructor
	/// there is no need to define BitwiseCopy and Copy!
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// eventviewer method
	void ResetMessage() const;
	void AddToMessage( const string& msg ) const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/