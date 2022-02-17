/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file VanillaDigitalArgSFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      caps for the SFRM model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPMODEL_VANILLADIGITALARGSFRM_H
#define _INGPMODEL_VANILLADIGITALARGSFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/port.h"
#include "gpcalib/vanilladigital.h"
#include "gpinfra/typedef.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaDigitalArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaDigitalArgSFRM: public ARM_VanillaDigitalArg
{
	ARM_VanillaDigitalArgSFRM( const ARM_VanillaDigitalArg& rhs );

	/// copy constructor
	ARM_VanillaDigitalArgSFRM( const ARM_VanillaDigitalArgSFRM& rhs );
    ARM_VanillaDigitalArgSFRM& operator=(const ARM_VanillaDigitalArgSFRM& rhs);
	virtual ~ARM_VanillaDigitalArgSFRM();
	virtual ARM_Object* Clone();

	/// accessors
	inline void SetLibors( const vector<ARM_VectorPtr>& libors )    { itsLibors = libors; }
	inline void SetZCPays( const vector<ARM_VectorPtr>& ZCPays )    { itsZCPays = ZCPays; }

private :
	friend class ARM_SFRM;
	vector<ARM_VectorPtr>   itsLibors;
	vector<ARM_VectorPtr>   itsZCPays;

// FIXMEFRED: mig.vc8 (25/05/2007 15:44:51): missing return type
	void CopyNoCleanUp(const ARM_VanillaDigitalArgSFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
