/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file eventviewer.h
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef _INGPBASE_WARNING_H
#define _INGPBASE_WARNING_H

#include "port.h"
#include "rootobject.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

class ARM_Warning : public ARM_RootObject
{
public:
	ARM_Warning( const string& message )
	:	ARM_RootObject(), itsMessage( message )	{ PopUp();}
	ARM_Warning( const ARM_Warning& rhs ): ARM_RootObject(rhs), itsMessage(rhs.itsMessage) {}
	ARM_Warning& operator=( const ARM_Warning& rhs )
	{
		if( this != &rhs )
			new (this) ARM_Warning(rhs);
		return *this;
	}

	void AddToMessage( const string& message );
	ARM_Object* Clone() const { return new ARM_Warning(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return itsMessage; }
	virtual string ExportShortName() const { return "LWARN";}
	static void SetPopUp( bool flag ) { ARM_Warning::isPopUpOn = flag; }

	///Accessors
	inline string GetMessage() const {return itsMessage; }
	inline void SetMessage(const string& message) { itsMessage = message; }

private:
	void PopUp() const;
	static bool isPopUpOn;
	string itsMessage;

	static const int MaxMessageSize;
};

CC_END_NAMESPACE()

#endif