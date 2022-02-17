/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file VanillaDigitalArgSFRM.cpp
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpmodels/VanillaDigitalArgSFRM.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaDigitalArgSFRM::CopyNoCleanUp(const ARM_VanillaDigitalArgSFRM& rhs)
{
	itsLibors			= rhs.itsLibors;
	itsZCPays			= rhs.itsZCPays;
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds empty object
////////////////////////////////////////////////////
ARM_VanillaDigitalArgSFRM::ARM_VanillaDigitalArgSFRM( const ARM_VanillaDigitalArg& arg )
:	
	ARM_VanillaDigitalArg( arg ),
	itsLibors(),
	itsZCPays()
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaDigitalArgSFRM::ARM_VanillaDigitalArgSFRM(const ARM_VanillaDigitalArgSFRM& rhs )
:	
	ARM_VanillaDigitalArg(rhs),
	itsLibors(),
	itsZCPays()
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaDigitalArgSFRM& ARM_VanillaDigitalArgSFRM::operator=(const ARM_VanillaDigitalArgSFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaDigitalArg::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: ~ARM_VanillaDigitalArgSFRM
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_VanillaDigitalArgSFRM::~ARM_VanillaDigitalArgSFRM()
{};



////////////////////////////////////////////////////
///	Struct : ARM_VanillaDigitalArgSFRM
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_VanillaDigitalArgSFRM::Clone()
{
	return new ARM_VanillaDigitalArgSFRM(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

