/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ARM_VanillaCapArgSFRM.cpp
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpmodels/VanillaCapArgSFRM.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaCapArgSFRM::CopyNoCleanUp(const ARM_VanillaCapArgSFRM& rhs)
{
	itsLibors			= rhs.itsLibors;
	itsZCPays			= rhs.itsZCPays;
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds empty object
////////////////////////////////////////////////////
ARM_VanillaCapArgSFRM::ARM_VanillaCapArgSFRM( const ARM_VanillaCapArg& arg )
:	
	ARM_VanillaCapArg( arg ),
	itsLibors(),
	itsZCPays()
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaCapArgSFRM::ARM_VanillaCapArgSFRM(const ARM_VanillaCapArgSFRM& rhs )
:	
	ARM_VanillaCapArg(rhs),
	itsLibors(),
	itsZCPays()
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaCapArgSFRM& ARM_VanillaCapArgSFRM::operator=(const ARM_VanillaCapArgSFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaCapArg::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: ~ARM_VanillaCapArgSFRM
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_VanillaCapArgSFRM::~ARM_VanillaCapArgSFRM()
{};



////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArgSFRM
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_VanillaCapArgSFRM::Clone()
{
	return new ARM_VanillaCapArgSFRM(*this);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

