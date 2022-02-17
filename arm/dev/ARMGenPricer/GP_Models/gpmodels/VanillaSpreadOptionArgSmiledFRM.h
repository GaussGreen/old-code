/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: VanillaCapArgSFRM.h,v $
 * Revision 1.1  2004/03/31 07:51:19  ebenhamou
 * Initial revision
 *
 */

/*! \file VanillaSwaptionArgSFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      swaption for the SFRM model
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date March 2004
 */


#ifndef _INGPMODEL_VANILLASPREADOPTIONLETARGSFRM_H
#define _INGPMODEL_VANILLASPREADOPTIONLETARGSFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/port.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpinfra/typedef.h"

#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSwaptionArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
class ARM_VanillaSpreadOptionArgSmiledFRM: public ARM_VanillaSwaptionArgSmiledFRM
{
public:
	ARM_VanillaSpreadOptionArgSmiledFRM();
	ARM_VanillaSpreadOptionArgSmiledFRM( const ARM_VanillaSpreadOptionArgSmiledFRM& rhs );
    ARM_VanillaSpreadOptionArgSmiledFRM& operator=(const ARM_VanillaSpreadOptionArgSmiledFRM& rhs);
    virtual ~ARM_VanillaSpreadOptionArgSmiledFRM();
	virtual ARM_Object* Clone();

	/// accessor
	inline ARM_VanillaSwaptionArgSmiledFRM* GetSwaptionArg_Long() const { return itsSwaptionArg_Long; }
	inline ARM_VanillaSwaptionArgSmiledFRM* GetSwaptionArg_Short() const { return itsSwaptionArg_Short; }
	inline ARM_VanillaSwaptionArgSmiledFRM* GetSwaptionArg_ToTimePayTime() const { return itsSwaptionArg_ToTimePayTime; }
	inline ARM_VanillaSwaptionArgSmiledFRM* GetSwaptionArg_StartTimePayTime() const { return itsSwaptionArg_StartTimePayTime; }

	inline void SetSwaptionArg_Long(const ARM_VanillaSwaptionArgSmiledFRM& swaptionArg) const {((ARM_VanillaSwaptionArgSmiledFRM*)itsSwaptionArg_Long)= new ARM_VanillaSwaptionArgSmiledFRM(swaptionArg); } 
	inline void SetSwaptionArg_Short(const ARM_VanillaSwaptionArgSmiledFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSmiledFRM*)itsSwaptionArg_Short) = new ARM_VanillaSwaptionArgSmiledFRM(swaptionArg); }
	inline void SetSwaptionArg_ToTimePayTime(const ARM_VanillaSwaptionArgSmiledFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSmiledFRM*)itsSwaptionArg_ToTimePayTime) = new ARM_VanillaSwaptionArgSmiledFRM(swaptionArg); }
	inline void SetSwaptionArg_StartTimePayTime(const ARM_VanillaSwaptionArgSmiledFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSmiledFRM*)itsSwaptionArg_StartTimePayTime) = new ARM_VanillaSwaptionArgSmiledFRM(swaptionArg); }


private :

	ARM_VanillaSwaptionArgSmiledFRM* itsSwaptionArg_Long;
	ARM_VanillaSwaptionArgSmiledFRM* itsSwaptionArg_Short;
	ARM_VanillaSwaptionArgSmiledFRM* itsSwaptionArg_ToTimePayTime;
	ARM_VanillaSwaptionArgSmiledFRM* itsSwaptionArg_StartTimePayTime;

// FIXMEFRED: mig.vc8 (30/05/2007 16:31:59):no return type
	void CopyNoCleanUp(const ARM_VanillaSpreadOptionArgSmiledFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
