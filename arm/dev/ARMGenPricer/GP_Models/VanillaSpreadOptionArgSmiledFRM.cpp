/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ARM_VanillaSwaptionArgSFRM.cpp
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"
#include "gpmodels/VanillaSpreadOptionArgSmiledFRM.h"

/// gpbase
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSmiledFRM
///	Routine: constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArgSmiledFRM::ARM_VanillaSpreadOptionArgSmiledFRM()
:	itsSwaptionArg_Long(NULL),
	itsSwaptionArg_Short(NULL),
	itsSwaptionArg_ToTimePayTime(NULL),
	itsSwaptionArg_StartTimePayTime(NULL)
{};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaSpreadOptionArgSmiledFRM::CopyNoCleanUp(const ARM_VanillaSpreadOptionArgSmiledFRM& rhs)
{
	itsSwaptionArg_Long = rhs.itsSwaptionArg_Long ? (ARM_VanillaSwaptionArgSmiledFRM*)rhs.itsSwaptionArg_Long->Clone() :  NULL;
	itsSwaptionArg_Short = rhs.itsSwaptionArg_Short ? (ARM_VanillaSwaptionArgSmiledFRM*)rhs.itsSwaptionArg_Short->Clone() :  NULL;
	itsSwaptionArg_ToTimePayTime = rhs.itsSwaptionArg_ToTimePayTime ? (ARM_VanillaSwaptionArgSmiledFRM*)rhs.itsSwaptionArg_ToTimePayTime->Clone() :  NULL;
	itsSwaptionArg_StartTimePayTime = rhs.itsSwaptionArg_StartTimePayTime ? (ARM_VanillaSwaptionArgSmiledFRM*)rhs.itsSwaptionArg_StartTimePayTime->Clone() :  NULL;
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds semi-empty object
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArgSmiledFRM::ARM_VanillaSpreadOptionArgSmiledFRM( const ARM_VanillaSpreadOptionArgSmiledFRM& rhs )
:	itsSwaptionArg_Long(NULL),
	itsSwaptionArg_Short(NULL),
	itsSwaptionArg_ToTimePayTime(NULL),
	itsSwaptionArg_StartTimePayTime(NULL)
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArgSmiledFRM& ARM_VanillaSpreadOptionArgSmiledFRM::operator=(const ARM_VanillaSpreadOptionArgSmiledFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaSpreadOptionArgSmiledFRM::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: ~ARM_VanillaSwaptionArgSFRM
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArgSmiledFRM::~ARM_VanillaSpreadOptionArgSmiledFRM()
{
	delete itsSwaptionArg_Long;
    itsSwaptionArg_Long = NULL;

	delete itsSwaptionArg_Short;
    itsSwaptionArg_Short = NULL;

	delete itsSwaptionArg_ToTimePayTime;
    itsSwaptionArg_ToTimePayTime = NULL;

	delete itsSwaptionArg_StartTimePayTime;
    itsSwaptionArg_StartTimePayTime = NULL;
}



////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_VanillaSpreadOptionArgSmiledFRM::Clone()
{
	return new ARM_VanillaSpreadOptionArgSmiledFRM(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

