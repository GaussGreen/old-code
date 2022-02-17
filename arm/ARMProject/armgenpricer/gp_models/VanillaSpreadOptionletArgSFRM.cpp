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

#include "gpmodels/VanillaSwaptionArgSFRM.h"
#include "gpmodels/VanillaSpreadOptionletArgSFRM.h"

/// gpbase
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaSpreadOptionletArgSFRM::CopyNoCleanUp(const ARM_VanillaSpreadOptionletArgSFRM& rhs)
{
	itsSwaptionArg_Long = rhs.itsSwaptionArg_Long ? (ARM_VanillaSwaptionArgSFRM*)rhs.itsSwaptionArg_Long->Clone() :  NULL;
	itsSwaptionArg_Short = rhs.itsSwaptionArg_Short ? (ARM_VanillaSwaptionArgSFRM*)rhs.itsSwaptionArg_Short->Clone() :  NULL;
	itsSwaptionArg_ToTimePayTime = rhs.itsSwaptionArg_ToTimePayTime ? (ARM_VanillaSwaptionArgSFRM*)rhs.itsSwaptionArg_ToTimePayTime->Clone() :  NULL;
	itsSwaptionArg_StartTimePayTime = rhs.itsSwaptionArg_StartTimePayTime ? (ARM_VanillaSwaptionArgSFRM*)rhs.itsSwaptionArg_StartTimePayTime->Clone() :  NULL;
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds semi-empty object
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM::ARM_VanillaSpreadOptionletArgSFRM( const ARM_VanillaSpreadOptionArg& rhs )
:	ARM_VanillaSpreadOptionArg( rhs ),
	itsSwaptionArg_Long(NULL),
	itsSwaptionArg_Short(NULL),
	itsSwaptionArg_ToTimePayTime(NULL),
	itsSwaptionArg_StartTimePayTime(NULL)
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM::ARM_VanillaSpreadOptionletArgSFRM(const ARM_VanillaSpreadOptionletArgSFRM& arg)
: ARM_VanillaSpreadOptionArg(arg)
{
    CopyNoCleanUp(arg);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM& ARM_VanillaSpreadOptionletArgSFRM::operator=(const ARM_VanillaSpreadOptionletArgSFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaSpreadOptionArg::operator=(rhs);
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
ARM_VanillaSpreadOptionletArgSFRM::~ARM_VanillaSpreadOptionletArgSFRM()
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
ARM_Object* ARM_VanillaSpreadOptionletArgSFRM::Clone()
{
	return new ARM_VanillaSpreadOptionletArgSFRM(*this);
}




////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Constructor
///	Returns: 
///	Action : creates the object!
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM::ARM_VanillaSpreadOptionletArgSFRM(
		const string& curveName,
		double evalTime,
		int CallPut,
		double expiryTime,
		double startTime,
		double endTime,
		std::vector<double>& resetTimes,
		std::vector<double>& payTimes,
		std::vector<double>& payPeriods,
		std::vector<double>& notional,
		std::vector<double>& coeffLong,
		std::vector<double>& coeffShort,
		std::vector<double>& strikes,
		std::vector<double>& swapLongFloatStartTime,
		std::vector<double>& swapLongFloatEndTime,
		ARM_VectorVector swapLongFixPayTimes,
		ARM_VectorVector swapLongFixPayPeriods,
		std::vector<double>& swapShortFloatStartTime,
		std::vector<double>& swapShortFloatEndTime,
		ARM_VectorVector swapShortFixPayTimes,
		ARM_VectorVector swapShortFixPayPeriods)
:
	ARM_VanillaSpreadOptionArg(
		curveName,
		evalTime,
		CallPut,
		expiryTime,
		startTime,
		endTime,
		resetTimes,
		payTimes,
		payPeriods,
		notional,
		coeffLong,
		coeffShort,
		strikes,
		swapLongFloatStartTime,
		swapLongFloatEndTime,
		swapLongFixPayTimes,
		swapLongFixPayPeriods,
		swapShortFloatStartTime,
		swapShortFloatEndTime,
		swapShortFixPayTimes,
		swapShortFixPayPeriods),
	itsSwaptionArg_Long(NULL),
	itsSwaptionArg_Short(NULL),
	itsSwaptionArg_ToTimePayTime(NULL),
	itsSwaptionArg_StartTimePayTime(NULL)
{
	}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

