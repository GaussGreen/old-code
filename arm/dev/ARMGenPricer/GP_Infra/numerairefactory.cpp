/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file numerairefactory.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

#include "gpinfra/numerairefactory.h"
#include "gpinfra/numeraire.h"
#include "gpbase/singleton.h"

CC_BEGIN_NAMESPACE( ARM )




////////////////////////////////////////////////////
///	Class  : ARM_NumeraireFactoryImp
///	Routine: CreateNumeraire
///	Returns: ARM_Numeraire*
///	Action : Create a numeraire
////////////////////////////////////////////////////
ARM_Numeraire* ARM_NumeraireFactoryImp::CreateNumeraire(int numeraireType )
{
	switch( numeraireType )
	{
	case ARM_Numeraire::Cash:
		return new ARM_NumeraireCash(false);
	case ARM_Numeraire::TerminalZc:
		return new ARM_NumeraireTerminalZc;
	case ARM_Numeraire::TerminalEventZc :
		return new ARM_NumeraireTerminalEventZc;
	case ARM_Numeraire::RollingPayment :
		return new ARM_NumeraireRollingPayment;
	case ARM_Numeraire::RollingEvent :
		return new ARM_NumeraireRollingEvent;
	case ARM_Numeraire::RollingCash :
		return new ARM_NumeraireRollingCash;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": unknown numeraire type, allowed TerminalZc, TerminalEventZc and Cash" );
	}
}

ARM_SingletonHolder<ARM_NumeraireFactoryImp> ARM_NumeraireFactory;



CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
