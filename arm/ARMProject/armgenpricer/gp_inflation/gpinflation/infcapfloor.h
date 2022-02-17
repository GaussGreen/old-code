/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: infcapfloor.h,v $
 * Revision 1.10  2003/12/15 18:21:00  ebenhamou
 * change for new unix compiler
 *
 * Revision 1.9  2003/10/24 14:53:50  ebenhamou
 * add GetConstantWithNoSpread
 *
 * Revision 1.8  2003/09/29 16:29:51  ebenhamou
 * more in the view and not overwritting of the forwards!
 *
 * Revision 1.7  2003/09/11 10:56:06  ebenhamou
 * storeVolInfo with argument of renormalisation.
 *
 * Revision 1.6  2003/09/10 15:15:59  ebenhamou
 * remove pricing strike to the model
 *
 * Revision 1.5  2003/09/09 11:54:54  ebenhamou
 * added volInfo for view method
 *
 * Revision 1.4  2003/09/08 10:34:47  ebenhamou
 * dos2unix
 *
 * Revision 1.3  2003/09/05 07:28:17  ebenhamou
 * ersion working with bsmodel
 *
 * Revision 1.2  2003/09/02 17:25:06  ebenhamou
 * version with some code
 *
 * Revision 1.1  2003/07/18 09:02:04  ebenhamou
 * Initial revision
 *
 *
 */

    
/*----------------------------------------------------------------------------*/

/*! \file inflcapfloor.h
 *
 *  \brief Vanilla inflation cap and floor
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */

 
/*----------------------------------------------------------------------------*/


/*! \class   ARM_InfCapFloor
 *	\brief  defines an inflation cap floor object. Inherits form ARM_CapFloor
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */


#ifndef _INGPINFLATION_INFCAPFLOOR_H
#define _INGPINFLATION_INFCAPFLOOR_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

/// inflation
#include "infleg.h"

/// kernel
#include <inst/capfloor.h>
#include <glob/dates.h>

CC_USING_NS(std,map)



CC_BEGIN_NAMESPACE( ARM )

class ARM_InfIdx;

/// very simple class to get the context of each caplet!
class ARM_InfCapFloorContext : public ARM_Object
{
private:
	ARM_Date itsNumDate;
 	ARM_Date itsDenomDate;
	int itsOptionType;
	double itsRenormalisationFactor;

public:
	ARM_InfCapFloorContext(
		const ARM_Date& numDate,
 		const ARM_Date& denomDate,
		int optionType,
		double renormalisationFactor )
	:
		itsNumDate( numDate ),
 		itsDenomDate( denomDate ),
		itsOptionType( optionType ),
		itsRenormalisationFactor( renormalisationFactor )
	{}
	
	inline const ARM_Date& GetNumDate() const {	return itsNumDate; }
	inline const ARM_Date& GetDenomDate() const { return itsDenomDate; }
	inline int GetOptionType() const { return itsOptionType; }
	inline double GetRenormalisationFactor() const { return itsRenormalisationFactor; }
};


class ARM_InfCapFloor : public ARM_CapFloor
{
private:
	ARM_Matrix* itsVolInfo;
	ARM_GP_Vector* itsCFValues;

	/// copy mehtod that only copies the member data above
	void CopyNoCleanUp( const ARM_InfCapFloor& infLeg );
	
	void CleanUp();
	void Init();
	double GetConstantWithNoSpread(int swapType) const;

public:

	ARM_InfCapFloor( const ARM_InfCapFloor& rhs );
	ARM_InfCapFloor& operator=( const ARM_InfCapFloor& rhs );

	ARM_InfCapFloor(
	const ARM_Date& startDate,
	const ARM_Date& endDate,
	const string& indexName,
	int capOrFloor,
	double strike,
	double leverage,
	double spread,
	int swapType,
	int rcvOrPay,
	int interpType					= K_CPILINEAR,
	int resetFreq					= K_DEF_FREQ,
	int dayCount					= KACTUAL_ACTUAL,
	const char* resetCalendar		= "INF",				/// calendar used for reset dates
	int fwdRule						= K_MOD_FOLLOWING,		/// whether fwds are with adjusted dates
	int intRule						= K_UNADJUSTED,			/// whether period of interest are adjusted or not
	int stubRule					= K_SHORTSTART,			/// ability to have K_SHORTSTART etc
	int resetGap					= GETDEFAULTVALUE,		/// reset gap default is -spotDay of calendar
	int payFreq						= GETDEFAULTVALUE,		/// payment frequency
	int payGap						= GETDEFAULTVALUE,		/// pay gap default is 0
	const char* payCalendar			= GETDEFAULTVALUESTR,	/// calendar used for payment
	int adjFirstDate				= 1,					/// we adjust for first date
	double firstReset				= GETDEFAULTVALUE,		/// ability to overwritte the first reset
	const ARM_Currency*	discountCcy = NULL );

	//// constructor with a swapleg
	ARM_InfCapFloor( ARM_SwapLeg* swapLeg, int capFloor, double strike = 0.0);
	//// constructor with a refValue
	ARM_InfCapFloor( ARM_SwapLeg* swapLeg, int capFloor, ARM_ReferenceValue *strike);

	/// desctructor
	virtual ~ARM_InfCapFloor();

	/// ARM_Object stuff
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// ARM pricing support
	virtual double ComputePrice(int mode = 0);		
	virtual void PropagateModel(ARM_Model *model);
	virtual void CptCashFlowValues(void);

	void StoreVol( double vol, double volLookupStrike, double strike, double timeToStart, 
		double tenor, int k );
};

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
