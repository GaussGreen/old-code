/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillacap.cpp
 *
 *  \brief vanilla cap is a vanilla interest rate cap
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

#include "gpcalib/vanillacap.h"

/// gpbase
#include "gpbase/gpvector.h"


/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingfunctionir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapDigitalArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaCapDigitalArg::CopyNoCleanUp(const ARM_VanillaCapDigitalArg& rhs)
{
#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "Nominals",		rhs.itsNominals   );
	ThrowErrorOnNullObject( "ResetTimes",	rhs.itsResetTimes );
	ThrowErrorOnNullObject( "StartTimes",	rhs.itsStartTimes );
	ThrowErrorOnNullObject( "EndTimes",		rhs.itsEndTimes   );
	ThrowErrorOnNullObject( "Strikes",		rhs.itsStrikes    );
	ThrowErrorOnNullObject( "PayTimes",	 	rhs.itsPayTimes   ); 
	ThrowErrorOnNullObject( "PayPeriods",	rhs.itsPayPeriods ); 
#endif
	itsNominals     = static_cast<ARM_GP_Vector*>( rhs.itsNominals->Clone() );
	itsResetTimes   = static_cast<ARM_GP_Vector*>( rhs.itsResetTimes->Clone());
	itsStartTimes   = static_cast<ARM_GP_Vector*>( rhs.itsStartTimes->Clone());
	itsEndTimes     = static_cast<ARM_GP_Vector*>( rhs.itsEndTimes->Clone());
	itsStrikes      = static_cast<ARM_GP_Vector*>( rhs.itsStrikes->Clone());
	itsPayTimes     = static_cast<ARM_GP_Vector*>( rhs.itsPayTimes->Clone());
	itsPayPeriods   = static_cast<ARM_GP_Vector*>( rhs.itsPayPeriods->Clone());
	itsAtTheMoneyFlag = rhs.itsAtTheMoneyFlag;

}


ARM_VanillaCapDigitalArg::ARM_VanillaCapDigitalArg(const ARM_VanillaCapDigitalArg& rhs)
:	
	ARM_VanillaArg(rhs),
	itsNominals(NULL),
	itsResetTimes(NULL), 
	itsStartTimes(NULL), 
	itsEndTimes(NULL),
	itsStrikes(NULL), 
	itsPayTimes(NULL), 
	itsPayPeriods(NULL)
{
	CopyNoCleanUp(rhs);
}


ARM_VanillaCapDigitalArg& ARM_VanillaCapDigitalArg::operator=(const ARM_VanillaCapDigitalArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaCapDigitalArg::~ARM_VanillaCapDigitalArg()
{
	CleanUp();
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapDigitalArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaCapDigitalArg::CleanUp()
{
	delete itsNominals;
	delete itsResetTimes;
	delete itsStartTimes;
    delete itsEndTimes;
    delete itsStrikes;
    delete itsPayTimes; 
	delete itsPayPeriods;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapDigitalArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaCapDigitalArg::ImpliedVol(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
	double impliedVol;
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
	
    if(IRModel)
		impliedVol = IRModel->ImpliedVol(*this);
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price caplet, please advise");
	
    return impliedVol;

}
////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapDigitalArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaCapDigitalArg::Price(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );
	double OpletPrice;
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
	
    double totalPrice = 0.0;
    if(IRModel)
    {
        for(size_t i=0;i<itsResetTimes->size();++i)
        {
			/// in order to communicate the current index of a closed form
			/// we use the dumStates to dump the current index!
			dumStates->SetModelState(0,0,i);

			/// Set the strike At the money at the value of the Libor calculated by the IRModel
			if (itsAtTheMoneyFlag)
			{
				double period = ((*(itsEndTimes))[i] - (*(itsStartTimes))[i])/360	;			
				(*(itsStrikes))[i]=(*(IRModel->Libor(GetCurveName(),0,(*(itsStartTimes))[i],(*(itsEndTimes))[i],	(*(itsPayPeriods))[i],(*(itsStartTimes))[i],	(*(itsEndTimes))[i], ARM_PricingStatesPtr(NULL))))[0];
			}
			
			OpletPrice = PriceOplet(
				IRModel, 
				GetCurveName(),
				GetEvalTime(),
				(*(itsPayTimes))[i],
				(*(itsPayPeriods))[i],
				(*(itsNominals))[i],
				(*(itsResetTimes))[i],
				(*(itsStartTimes))[i],
				(*(itsEndTimes))[i],
				(*(itsPayPeriods))[i], // index period = payment period
				(*(itsStrikes))[i],
				GetCallPut(),
				dumStates);
			
			totalPrice += OpletPrice;
        }
    }
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price caplet, please advise");
	
    return totalPrice;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArg
///	Routine: PriceOplet
///	Returns: double
///	Action : price a caplet
////////////////////////////////////////////////////

double ARM_VanillaCapArg::PriceOplet(
	ARM_PricingFunctionIR* IRModel,
	const string& curveName, 
	double evalTime,
	double payTime,
	double period,
	double payNotional,
	double fwdResetTime,	/// used for volatility computation
	double fwdStartTime,
	double fwdEndTime,
	double fwdPeriod,
	double strike,
	int capFloor,
	const ARM_PricingStatesPtr& states ) const
{
	static ARM_VectorPtr result;

	result = IRModel->VanillaCapletScalar(
		curveName,
		evalTime,
		payTime,
		period,
		payNotional,
		fwdResetTime,
		fwdStartTime,
		fwdEndTime,
		fwdPeriod,
		strike,
		capFloor,
		states);
	return (*result)[0];
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaCapArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaCapArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

