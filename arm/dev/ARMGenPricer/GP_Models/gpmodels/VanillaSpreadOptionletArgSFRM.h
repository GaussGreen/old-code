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
#include "gpmodels/VanillaSwaptionArgSFRM.h"

#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;


///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaSwaptionArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaSpreadOptionletArgSFRM: public ARM_VanillaSpreadOptionArg
{
	ARM_VanillaSpreadOptionletArgSFRM(
		const string& curveName,
		double evalTime,
		int CallPut,
		double expiryTime,
		double startTime,
		double endTime,
		ARM_GP_Vector* resetTimes,
		ARM_GP_Vector* payTimes,
		ARM_GP_Vector* payPeriods,
		ARM_GP_Vector* notional,
		ARM_GP_Vector* coeffLong,
		ARM_GP_Vector* coeffShort,
		ARM_GP_Vector* strikes,
		ARM_GP_Vector* swapLongFloatStartTime,
		ARM_GP_Vector* swapLongFloatEndTime,
		ARM_VectorVector swapLongFixPayTimes,
		ARM_VectorVector swapLongFixPayPeriods,
		ARM_GP_Vector* swapShortFloatStartTime,
		ARM_GP_Vector* swapShortFloatEndTime,
		ARM_VectorVector swapShortFixPayTimes,
		ARM_VectorVector swapShortFixPayPeriods);

	ARM_VanillaSpreadOptionletArgSFRM( const ARM_VanillaSpreadOptionArg& rhs );
	ARM_VanillaSpreadOptionletArgSFRM( const ARM_VanillaSpreadOptionletArgSFRM& rhs );
    ARM_VanillaSpreadOptionletArgSFRM& operator=(const ARM_VanillaSpreadOptionletArgSFRM& rhs);
    virtual ~ARM_VanillaSpreadOptionletArgSFRM();
	virtual ARM_Object* Clone();

	/// accessor
	inline ARM_VanillaSwaptionArgSFRM* GetSwaptionArg_Long() const { return itsSwaptionArg_Long; }
	inline ARM_VanillaSwaptionArgSFRM* GetSwaptionArg_Short() const { return itsSwaptionArg_Short; }
	inline ARM_VanillaSwaptionArgSFRM* GetSwaptionArg_ToTimePayTime() const { return itsSwaptionArg_ToTimePayTime; }
	inline ARM_VanillaSwaptionArgSFRM* GetSwaptionArg_StartTimePayTime() const { return itsSwaptionArg_StartTimePayTime; }

	inline void SetSwaptionArg_Long(const ARM_VanillaSwaptionArgSFRM& swaptionArg) const {((ARM_VanillaSwaptionArgSFRM*)itsSwaptionArg_Long)= new ARM_VanillaSwaptionArgSFRM(swaptionArg); } 
	inline void SetSwaptionArg_Short(const ARM_VanillaSwaptionArgSFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSFRM*)itsSwaptionArg_Short) = new ARM_VanillaSwaptionArgSFRM(swaptionArg); }
	inline void SetSwaptionArg_ToTimePayTime(const ARM_VanillaSwaptionArgSFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSFRM*)itsSwaptionArg_ToTimePayTime) = new ARM_VanillaSwaptionArgSFRM(swaptionArg); }
	inline void SetSwaptionArg_StartTimePayTime(const ARM_VanillaSwaptionArgSFRM& swaptionArg) const { ((ARM_VanillaSwaptionArgSFRM*)itsSwaptionArg_StartTimePayTime) = new ARM_VanillaSwaptionArgSFRM(swaptionArg); }


private :

	ARM_VanillaSwaptionArgSFRM* itsSwaptionArg_Long;
	ARM_VanillaSwaptionArgSFRM* itsSwaptionArg_Short;
	ARM_VanillaSwaptionArgSFRM* itsSwaptionArg_ToTimePayTime;
	ARM_VanillaSwaptionArgSFRM* itsSwaptionArg_StartTimePayTime;

	void CopyNoCleanUp(const ARM_VanillaSpreadOptionletArgSFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
