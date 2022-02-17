/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *	\file VanillaCorridorLegArgSFRM.cpp
 *
 *  \brief this files enables to store various data to fast calibrate
 *      caps for the SFRM model
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date May 2004
 */
#include "gpmodels/VanillaCorridorLegArgSFRM.h"

/// gpbase
#include <glob/expt.h>
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/gpvector.h"

/// gpinfra
#include "gpinfra/typedef.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaCorridorLegArgSFRM::CopyNoCleanUp(const ARM_VanillaCorridorLegArgSFRM& rhs)
{
    itsLibors       = rhs.itsLibors;
    itsZCPays       = rhs.itsZCPays;

    itsRef_Libors   = rhs.itsRef_Libors;;
    itsRef_DfFwds   = rhs.itsRef_DfFwds;;
    itsRef_ZCEnds   = rhs.itsRef_ZCEnds;;
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: default constructor
///	Returns: 
///	Action : builds completely empty object
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM::ARM_VanillaCorridorLegArgSFRM()
:	
	ARM_VanillaCorridorLegArg()
{}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds empty object
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM::ARM_VanillaCorridorLegArgSFRM( const ARM_VanillaCorridorLegArg& arg )
:	
	ARM_VanillaCorridorLegArg( arg ),
	itsLibors(),
	itsZCPays(),

    itsRef_Libors(),
    itsRef_DfFwds(),
    itsRef_ZCEnds()
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM::ARM_VanillaCorridorLegArgSFRM(const ARM_VanillaCorridorLegArgSFRM& rhs )
:	
	ARM_VanillaCorridorLegArg(rhs),
	itsLibors(),
	itsZCPays(),

    itsRef_Libors(),
    itsRef_DfFwds(),
    itsRef_ZCEnds()
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM& ARM_VanillaCorridorLegArgSFRM::operator=(const ARM_VanillaCorridorLegArgSFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaCorridorLegArg::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: ~ARM_VanillaCorridorLegArgSFRM
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM::~ARM_VanillaCorridorLegArgSFRM()
{};



////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArgSFRM
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_VanillaCorridorLegArgSFRM::Clone()
{
	return new ARM_VanillaCorridorLegArgSFRM(*this);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

