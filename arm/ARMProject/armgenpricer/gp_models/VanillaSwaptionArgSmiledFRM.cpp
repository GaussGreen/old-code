/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file ARM_VanillaSwaptionArgSmiledFRM.cpp
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSmiledFRM
///	Routine: constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM::ARM_VanillaSwaptionArgSmiledFRM()
:	itsAverageShift(0),
	itsMu(0),
	itsFixAnnuity(0),
	itsSwapFwd(0)
{};


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSmiledFRM
///	Routine: copy constructor
///	Returns: 
///	Action : builds the object
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM::ARM_VanillaSwaptionArgSmiledFRM(const ARM_VanillaSwaptionArgSmiledFRM& rhs)
:	itsAverageShift(rhs.itsAverageShift),
	itsMu(rhs.itsMu),
	itsFixAnnuity(rhs.itsFixAnnuity),
	itsSwapFwd(rhs.itsSwapFwd)
{};


////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSmiledFRM
///	Routine: ~ARM_VanillaSwaptionArgSmiledFRM
///	Returns: 
///	Action : destructor
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM::~ARM_VanillaSwaptionArgSmiledFRM()
{};

////////////////////////////////////////////////////
///	Class   : ARM_VanillaSwaptionArgSmiledFRM
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_VanillaSwaptionArgSmiledFRM::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    return os.str();
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaSwaptionArgSmiledFRM
///	Routine: Constructor
///	Returns: 
///	Action : creates the object!
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM::ARM_VanillaSwaptionArgSmiledFRM(	double averageShift,
																	const ARM_VectorPtr& mu,
																	const ARM_VectorPtr& fixAnnuity,
																	const ARM_VectorPtr& swapFwd)
:	itsAverageShift(averageShift),
	itsMu(mu),
	itsFixAnnuity(fixAnnuity),
	itsSwapFwd(swapFwd)
{};

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

