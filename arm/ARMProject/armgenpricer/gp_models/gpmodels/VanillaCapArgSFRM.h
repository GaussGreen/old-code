/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file VanillaCapArgSFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      caps for the SFRM model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPMODEL_VANILLACAPARGSFRM_H
#define _INGPMODEL_VANILLACAPARGSFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/port.h"
#include "gpcalib/vanillacap.h"
#include "gpinfra/typedef.h"
#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaCapArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaCapArgSFRM: public ARM_VanillaCapArg
{
	ARM_VanillaCapArgSFRM( const ARM_VanillaCapArg& rhs );

	/// copy constructor
	ARM_VanillaCapArgSFRM( const ARM_VanillaCapArgSFRM& rhs );
    ARM_VanillaCapArgSFRM& operator=(const ARM_VanillaCapArgSFRM& rhs);
	virtual ~ARM_VanillaCapArgSFRM();
	virtual ARM_Object* Clone();

	/// accessors
	inline void SetLibors( const vector<ARM_VectorPtr>& libors )    { itsLibors = libors; }
	inline void SetZCPays( const vector<ARM_VectorPtr>& ZCPays )    { itsZCPays = ZCPays; }

private :
	friend class ARM_SFRM;
	vector<ARM_VectorPtr>   itsLibors;
	vector<ARM_VectorPtr>   itsZCPays;

    void CopyNoCleanUp(const ARM_VanillaCapArgSFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
