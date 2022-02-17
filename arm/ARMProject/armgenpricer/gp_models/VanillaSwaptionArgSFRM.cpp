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

/// gpbase
#include "gpbase/gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Copy NoCleanUp
///	Returns: 
///	Action : copies only the new layer of data member
////////////////////////////////////////////////////
void ARM_VanillaSwaptionArgSFRM::CopyNoCleanUp(const ARM_VanillaSwaptionArgSFRM& rhs)
{
	itsAverageShift				= rhs.itsAverageShift;
	itsMu						= rhs.itsMu;
	itsFixAnnuity				= rhs.itsFixAnnuity;
	itsSwapFwd					= rhs.itsSwapFwd;
	itsFixAnnuityWithNominal	= rhs.itsFixAnnuityWithNominal;
	itsSwapFwdWithNominal		= rhs.itsSwapFwdWithNominal;

}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: constructor
///	Returns: 
///	Action : builds semi-empty object
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM::ARM_VanillaSwaptionArgSFRM( const ARM_VanillaSwaptionArg& rhs )
:	ARM_VanillaSwaptionArg( rhs ),
	itsAverageShift(),
	itsMu(),
	itsFixAnnuity(),
	itsSwapFwd(),
	itsFixAnnuityWithNominal(),
	itsSwapFwdWithNominal()
{}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM::ARM_VanillaSwaptionArgSFRM(const ARM_VanillaSwaptionArgSFRM& arg)
: ARM_VanillaSwaptionArg(arg)
{
    CopyNoCleanUp(arg);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: operator=
///	Returns: 
///	Action : assignment operator
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM& ARM_VanillaSwaptionArgSFRM::operator=(const ARM_VanillaSwaptionArgSFRM& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaSwaptionArg::operator=(rhs);
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
ARM_VanillaSwaptionArgSFRM::~ARM_VanillaSwaptionArgSFRM()
{};



////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Clone
///	Returns: 
///	Action : clone the object
////////////////////////////////////////////////////
ARM_Object* ARM_VanillaSwaptionArgSFRM::Clone()
{
	return new ARM_VanillaSwaptionArgSFRM(*this);
}




////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSFRM
///	Routine: Constructor
///	Returns: 
///	Action : creates the object!
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM::ARM_VanillaSwaptionArgSFRM(
	double expiryTime,
	double startTime,
	double endTime,
	double averageShift,
	std::vector<double>& fixPayTimes,
	std::vector<double>& fixPayPeriods,
	std::vector<double>& floatResetTimes,
	std::vector<double>& floatStartTimes,
	std::vector<double>& floatEndTimes,
	std::vector<double>& floatIntTerms,
	const ARM_VectorPtr& mu,
	const ARM_VectorPtr& fixAnnuity,
	const ARM_VectorPtr& swapFwd,
    int fixFrequency,
	const ARM_VectorPtr& fixAnnuityWithNominal,
	const ARM_VectorPtr& swapFwdWithNominal,
	std::vector<double>& fixNominal,
	std::vector<double>& floatNominal)
:
	ARM_VanillaSwaptionArg(
		"",					///		const string& curveName
		0.0,				///		double evalTime
		K_CALL,				///		int CallPut
		new std::vector<double>(1,0.0),
		fixNominal,			///		std::vector<double>& fixNominal	
		floatNominal,		///		std::vector<double>& floatNominal
		expiryTime,			///		double expiryTime
		startTime,			///		double startTime
		endTime,			///		double endTime
		new std::vector<double>(1,0.0),	///		std::vector<double>& strikes	
		fixPayTimes,		///		std::vector<double>& fixPayTimes
		fixPayPeriods,		///		std::vector<double>& fixPayPeriods
		floatResetTimes,	///		std::vector<double>& floatResetTimes
		floatStartTimes,	///		std::vector<double>& floatStartTimes
		floatEndTimes,		///		std::vector<double>& floatEndTimes
		floatIntTerms,		///		std::vector<double>& floatIntTerms
		fixFrequency,		///     int FixFrequency
		K_SEMIANNUAL,		///		Default Parameters for Frequency, DayCount
		K30_360,
		KACTUAL_360
		),
	itsAverageShift(averageShift),
	itsMu(mu),
	itsFixAnnuity(fixAnnuity),
	itsSwapFwd(swapFwd),
	itsFixAnnuityWithNominal(fixAnnuityWithNominal),
	itsSwapFwdWithNominal(swapFwdWithNominal)
{}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

