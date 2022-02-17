/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pathschemefactory.h
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPNUMMETHODS_PATHSCHEMEFACTORY_H
#define _INGPNUMMETHODS_PATHSCHEMEFACTORY_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
template <typename T> class ARM_SingletonHolder;
class ARM_PathScheme;

struct ARM_PathSchemeFactoryImp
{
    ARM_PathScheme* CreatePathScheme(
		int pathSchemeType);

private:
	/// to forbid client from using it except for the singleton holder
	ARM_PathSchemeFactoryImp() {};
	friend class ARM_SingletonHolder<ARM_PathSchemeFactoryImp>;
};

extern ARM_SingletonHolder<ARM_PathSchemeFactoryImp> ARM_PathSchemeFactory;


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

