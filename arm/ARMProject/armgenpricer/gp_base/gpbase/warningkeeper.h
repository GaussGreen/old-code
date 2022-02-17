/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file warningkeeper.h
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPBASE_WARNINGKERPER_H
#define _INGPBASE_WARNINGKERPER_H

#include "gpbase/port.h"
#include "gpbase/warning.h"
#include "typedef.h"

#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

struct ARM_WarningKeeper
{
	ARM_WarningKeeper(): itsWarning(0), itsKeepCurrentWarning(true) {}
	ARM_WarningKeeper( const ARM_WarningKeeper& rhs );
	ARM_WarningKeeper& operator=( const ARM_WarningKeeper& rhs );
	virtual ~ARM_WarningKeeper();

	inline ARM_Warning* GetWarning() const { if( ARM_WarningKeeper::KeepWarning && itsKeepCurrentWarning ) return itsWarning; return 0; }
	inline void SetKeepCurrentWarning( bool v ) { itsKeepCurrentWarning=v; }
	inline void SetWarning( ARM_Warning* w ) { if( ARM_WarningKeeper::KeepWarning && itsKeepCurrentWarning ) itsWarning=w;} /// no clone!

	void EraseWarningMessage() {itsWarning->SetMessage(string(""));}
	void AddWarningMessage( const string& message )
	{ 
		if( ARM_WarningKeeper::KeepWarning && itsKeepCurrentWarning ) 
			if( itsWarning ) 
				itsWarning->AddToMessage( message ); 
			else 
				itsWarning= new ARM_Warning( message ); 
	}
	static bool KeepWarning;	/// flag to disactivate any warning
private:
	ARM_Warning* itsWarning;
	bool itsKeepCurrentWarning;	/// flag to disactivate the current warning in the current class

};


CC_END_NAMESPACE()

#endif