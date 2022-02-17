/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1F_Eq.cpp
 *
 *  \brief Q model Equity 1 factor
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

/// gpinfra
#include "gpinfra/pricingmodeltype.h"

/// gpmodel
#include "gpmodels/Q1F_Eq.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsQ1F.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_QModel1F_Eq
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_QModel1F_Eq::ARM_QModel1F_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsQ1F_Eq* modelParam)
:	ARM_EqFxBase(zc,modelParam)
{
	if (!zc.IsNull())
	{
		itsSettlementCalendar = ComputeSettlementCalendar( GetModelName() );
		itsSettlementGap = ComputeSettlementGap( GetModelName() );
	}

	if(	!dynamic_cast<const ARM_ModelParamsQ1F_Eq*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Equity model param" );
}



////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Eq
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_QModel1F_Eq::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F Q Equity Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}



////////////////////////////////////////////////////
///	Class   : ARM_QModel1F_Eq
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_QModel1F_Eq::ComputeSettlementCalendar(const string& modelName) const
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
///	Class   : ARM_QModel1F_Eq
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_QModel1F_Eq::ComputeSettlementGap(const string& modelName) const
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
///	Class  : ARM_QModel1F_Eq
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_QModel1F_Eq::GetType() const
{
	return MT_EQUITY_MODEL;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

