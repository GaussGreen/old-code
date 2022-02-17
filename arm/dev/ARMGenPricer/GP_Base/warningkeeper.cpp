/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file warningkeeper.cpp
 *  \brief object for any warning
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2005
 */

#include "gpbase/warningkeeper.h"

CC_BEGIN_NAMESPACE( ARM )

bool ARM_WarningKeeper::KeepWarning = true;

ARM_WarningKeeper::ARM_WarningKeeper( const ARM_WarningKeeper& rhs )
:
	itsWarning( rhs.itsWarning? static_cast<ARM_Warning*>(rhs.itsWarning->Clone() ) :NULL ),
	itsKeepCurrentWarning( rhs.itsKeepCurrentWarning )
{}


ARM_WarningKeeper& ARM_WarningKeeper::operator =( const ARM_WarningKeeper& rhs )
{
	if( this != & rhs )
	{
		this->~ARM_WarningKeeper();
		new (this) ARM_WarningKeeper( rhs );
	}
	return *this;
}


ARM_WarningKeeper::~ARM_WarningKeeper()
{ delete itsWarning; itsWarning= 0; }


CC_END_NAMESPACE()
