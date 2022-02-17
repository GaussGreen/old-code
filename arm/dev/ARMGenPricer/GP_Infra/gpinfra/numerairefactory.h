/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file numerairefactory.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#ifndef _INGPINFRA_NUMERAIREFACTORY_H
#define _INGPINFRA_NUMERAIREFACTORY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/rootobject.h"
#include "gpbase/port.h"

#include "timeinfo.h"

CC_BEGIN_NAMESPACE( ARM ) 

/// forward declaration
template <typename T> class ARM_SingletonHolder;

class ARM_Numeraire;

struct ARM_NumeraireFactoryImp
{
	ARM_Numeraire* CreateNumeraire(int numeraireType );

private:
	/// to forbid client from using it except for the singleton holder
	ARM_NumeraireFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_NumeraireFactoryImp>;
};

extern ARM_SingletonHolder<ARM_NumeraireFactoryImp> ARM_NumeraireFactory;


CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
