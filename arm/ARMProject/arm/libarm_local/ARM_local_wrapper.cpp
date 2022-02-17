/*! \file ARM_local_wrapper.cpp
 *
 *  \brief file to factorise the interface
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2003
 */

#include "ARM_local_wrapper.h"
#include "ARM_local_persistent.h"

/*!
 * function to check the global persistance
 * enables to factorise code
 */
bool GlobalPersistanceOk( ARM_result& result)
{
	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return false;
	}
	else
		return true;
}


/// function that return always one (true)
int NOCHECK_RETURNALWAYSTRUE( ARM_Object*, ARM_CLASS_NAME )
{
	return 1;
}
