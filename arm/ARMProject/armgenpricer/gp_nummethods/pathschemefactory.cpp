/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pathschemefactory.cpp
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2005
 */


#include "gpnummethods/pathschemefactory.h"

/// gpbase
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"

/// gpnummethods
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/pathscheme.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_PathSchemeFactoryImp
///	Routine: CreatePathScheme
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_PathScheme* ARM_PathSchemeFactoryImp::CreatePathScheme(
        int pathSchemeType)
{
    ARM_PathScheme* pathScheme            = NULL;

    /// Sampler creation
    switch(pathSchemeType)
    {
	case ARM_PathScheme::Incremental:
		pathScheme = new ARM_IncrementalPathScheme();
		break;
    case ARM_PathScheme::BrownianBridge:
		pathScheme = new ARM_BrownianBridgePathScheme();
		break;
	case ARM_PathScheme::IncAdaptative:
		pathScheme = new ARM_IncrementAdaptativePathScheme();
		break;

    default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown path scheme type");
    }

	/// return the result
    return pathScheme;
}


ARM_SingletonHolder<ARM_PathSchemeFactoryImp> ARM_PathSchemeFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/