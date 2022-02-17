/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file LN_Eq.cpp
 *
 *  \brief Lognormal Equity 1 factor
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date December 2004
 */

/// gpinfra
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingstates.h"

/// gpmodel
#include "gpmodels/SABR_Eq.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SABR_Eq
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SABR_Eq::ARM_SABR_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsSABR_Eq* modelParams)
:	ARM_EqFxBase(zc,modelParams)
{
	if(	!dynamic_cast<const ARM_ModelParamsSABR_Eq*>( modelParams ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Equity model param" );
}


////////////////////////////////////////////////////
///	Class   : ARM_SABR_Eq
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_SABR_Eq::GetSettlementCalendar(const string& modelName) const
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
///	Class   : ARM_SABR_Eq
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_SABR_Eq::GetSettlementGap(const string& modelName) const
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
///	Class  : ARM_SABR_Eq
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_SABR_Eq::GetType() const
{
	return MT_EQUITY_MODEL;
}


////////////////////////////////////////////////////
///	Class  : ARM_SABR_Eq
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SABR_Eq::Forward(
	const string& curveName, 
	double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t stateSize = states->size();
		double forwardValue = ComputeFwdAtTime(settlementTime);
		return ARM_VectorPtr( new ARM_GP_Vector(stateSize,forwardValue) );
	}
	else
	{
		size_t modelNb = GetModelNb();
		size_t stateSize = states->size();
		ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize, 0.0 ));
		double forwardEnd = ComputeFwdAtTime( settlementTime );
		double forwardStart = ComputeFwdAtTime( evalTime );
		double spot;

		for( size_t i=0; i<stateSize; ++i )
		{
			spot = states->GetModelState(i,modelNb);
			(*fwdVector )[i] = spot* forwardEnd/forwardStart;
		}
		return fwdVector;
	}
};


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

