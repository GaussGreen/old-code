/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file LN_Eq.cpp
 *
 *  \brief Lognormal Equity 1 factor
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

/// gpinfra
#include "gpinfra/pricingmodeltype.h"

/// gpmodel
#include "gpmodels/LN_Eq.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Bs_ModelParams.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_LN_Eq
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_LN_Eq::ARM_LN_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsBS_Eq* modelParam)
:	ARM_EqFxBase(zc,modelParam)
{
	if(	!dynamic_cast<const ARM_ModelParamsBS_Eq*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Equity model param" );
}


////////////////////////////////////////////////////
///	Class   : ARM_LN_Eq
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_LN_Eq::GetSettlementCalendar(const string& modelName) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( GetZeroCurve() == ARM_ZeroCurvePtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no discount curve set!" );
	if( !GetZeroCurve()->GetCurrencyUnit() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no ccy unit set on the discount curve set!" );
#endif

	return GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
}


////////////////////////////////////////////////////
///	Class   : ARM_LN_Eq
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_LN_Eq::GetSettlementGap(const string& modelName) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( GetZeroCurve() == ARM_ZeroCurvePtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no discount curve set!" );
	if( !GetZeroCurve()->GetCurrencyUnit() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no ccy unit set on the discount curve set!" );
#endif

	return GetZeroCurve()->GetCurrencyUnit()->GetSpotDays();
}


////////////////////////////////////////////////////
///	Class  : ARM_LN_Eq
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_LN_Eq::GetType() const
{
	return MT_EQUITY_MODEL;
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

