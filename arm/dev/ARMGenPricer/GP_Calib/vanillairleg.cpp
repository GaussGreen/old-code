/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairleg.cpp
 *
 *  \brief vanilla interest rates leg
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#include "gpcalib/vanillairleg.h"

/// gpbase
#include "gpbase/gpvector.h"


/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaIRSwaplegArg::CopyNoCleanUp(const ARM_VanillaIRSwaplegArg& rhs)
{
#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "Nominals",		rhs.itsNominals   );
	ThrowErrorOnNullObject( "ResetTimes",	rhs.itsResetTimes );
	ThrowErrorOnNullObject( "StartTimes",	rhs.itsStartTimes );
	ThrowErrorOnNullObject( "EndTimes",		rhs.itsEndTimes   );
	ThrowErrorOnNullObject( "PayTimes",	 	rhs.itsPayTimes   ); 
	ThrowErrorOnNullObject( "PayPeriods",	rhs.itsPayPeriods ); 
#endif
	itsNominals     = static_cast<ARM_GP_Vector*>( rhs.itsNominals->Clone() );
	itsResetTimes   = static_cast<ARM_GP_Vector*>( rhs.itsResetTimes->Clone());
	itsStartTimes   = static_cast<ARM_GP_Vector*>( rhs.itsStartTimes->Clone());
	itsEndTimes     = static_cast<ARM_GP_Vector*>( rhs.itsEndTimes->Clone());
	itsPayTimes     = static_cast<ARM_GP_Vector*>( rhs.itsPayTimes->Clone());
	itsPayPeriods   = static_cast<ARM_GP_Vector*>( rhs.itsPayPeriods->Clone());
}


ARM_VanillaIRSwaplegArg::ARM_VanillaIRSwaplegArg(const ARM_VanillaIRSwaplegArg& rhs)
:	
	ARM_VanillaArg(rhs),
	itsNominals(NULL),
	itsResetTimes(NULL), 
	itsStartTimes(NULL), 
	itsEndTimes(NULL),
	itsPayTimes(NULL), 
	itsPayPeriods(NULL)
{
	CopyNoCleanUp(rhs);
}


ARM_VanillaIRSwaplegArg& ARM_VanillaIRSwaplegArg::operator=(const ARM_VanillaIRSwaplegArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaIRSwaplegArg::~ARM_VanillaIRSwaplegArg()
{
	CleanUp();
}

ARM_Object* ARM_VanillaIRSwaplegArg::Clone() const
{
	return new ARM_VanillaIRSwaplegArg(*this);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaIRSwaplegArg::CleanUp()
{
	delete itsNominals;
	delete itsResetTimes;
	delete itsStartTimes;
    delete itsEndTimes;
    delete itsPayTimes; 
	delete itsPayPeriods;
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaIRSwaplegArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaIRSwaplegArg::Price(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );
	double swaplegPrice;
	ARM_PricingModelIR* IRModel = dynamic_cast<ARM_PricingModelIR*>(model);
	
    double totalPrice = 0.0;
    if(IRModel)
    {
        for(size_t i=0;i<itsResetTimes->size();++i)
        {
			/// in order to communicate the current index of a closed form
			/// we use the dumStates to dump the current index!
			dumStates->SetModelState(0,0,i);
			
			swaplegPrice = PriceOplet(
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
				0.0,
				GetCallPut(),
				dumStates);
			
			totalPrice += swaplegPrice;
        }
    }
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingModelIR. So cannot price caplet, please advise");
	
    return totalPrice;
}




////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: PriceOplet
///	Returns: double
///	Action : price a swapleg
////////////////////////////////////////////////////

double ARM_VanillaIRSwaplegArg::PriceOplet(
	ARM_PricingModelIR* IRModel,
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
	static ARM_GP_VectorPtr result;

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
string ARM_VanillaIRSwaplegArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaCapArg"; 
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

/*
#include "gpcalib/vanillairleg.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"


/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_VanillaIRSwaplegArg::CopyNoCleanUp(const ARM_VanillaIRSwaplegArg& rhs)
{
    itsNominal          =   rhs.itsNominal;
	itsStartTime        =   rhs.itsStartTime;
    itsEndTime          =   rhs.itsEndTime;
    itsFixFrequency     =   rhs.itsFixFrequency;


#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "FixPayTimes",		rhs.itsFixPayTimes );
	ThrowErrorOnNullObject( "FixPayPeriods",	rhs.itsFixPayPeriods );
	ThrowErrorOnNullObject( "FloatResetTimes",	rhs.itsFloatResetTimes );
	ThrowErrorOnNullObject( "FloatStartTimes",	rhs.itsFloatStartTimes );
	ThrowErrorOnNullObject( "FloatEndTimes",	rhs.itsFloatEndTimes );
	ThrowErrorOnNullObject( "FloatIntTerms",	rhs.itsFloatIntTerms );
    ThrowErrorOnNullObject( "Strikes",	        rhs.itsStrikes );
#endif


	/// for fast access no test of NULL pointor in release!	
	itsFixPayTimes      =  static_cast<ARM_GP_Vector*>( rhs.itsFixPayTimes->Clone());
	itsFixPayPeriods    =  static_cast<ARM_GP_Vector*>( rhs.itsFixPayPeriods->Clone());
	itsFloatResetTimes  =  static_cast<ARM_GP_Vector*>( rhs.itsFloatResetTimes->Clone());
	itsFloatStartTimes  =  static_cast<ARM_GP_Vector*>( rhs.itsFloatStartTimes->Clone());
	itsFloatEndTimes    =  static_cast<ARM_GP_Vector*>( rhs.itsFloatEndTimes->Clone());
	itsFloatIntTerms    =  static_cast<ARM_GP_Vector*>( rhs.itsFloatIntTerms->Clone());
    itsStrikes          =  static_cast<ARM_GP_Vector*>( rhs.itsStrikes->Clone());


#if defined(__GP_STRICT_VALIDATION)
    CC_NS(ARM_Check,CheckSameArgSize)(*itsStrikes, *itsFixPayTimes, "itsStrikes", "FixPayTimes",__LINE__,__FILE__);
#endif
}


ARM_VanillaIRSwaplegArg::ARM_VanillaIRSwaplegArg(const ARM_VanillaIRSwaplegArg& arg)
:	ARM_VanillaArg(arg),
	itsFixPayTimes( NULL ),
	itsFixPayPeriods( NULL ),  
	itsFloatResetTimes( NULL ),
	itsFloatStartTimes( NULL ),
	itsFloatEndTimes( NULL ),
	itsFloatIntTerms( NULL ),
    itsStrikes(NULL),
    itsFixFrequency(K_SEMIANNUAL)
{
    CopyNoCleanUp(arg);
}

ARM_VanillaIRSwaplegArg& ARM_VanillaIRSwaplegArg::operator=(const ARM_VanillaIRSwaplegArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaIRSwaplegArg::~ARM_VanillaIRSwaplegArg()
{
	delete itsFixPayTimes;
	delete itsFixPayPeriods;
	delete itsFloatResetTimes;
	delete itsFloatStartTimes;
	delete itsFloatEndTimes;
	delete itsFloatIntTerms; 
    delete itsStrikes;

};

ARM_Object* ARM_VanillaIRSwaplegArg::Clone() const
{
	return new ARM_VanillaIRSwaplegArg(*this);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: Price
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaIRSwaplegArg::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
    ARM_VectorPtr Price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
    if(IRModel)
        Price = IRModel->VanillaSwaptionScalar(
		GetCurveName(),
		GetEvalTime(),
		GetExpiry(),
		itsNominal,
		itsStartTime,
		itsEndTime,
		*itsFixPayTimes,
		*itsFixPayPeriods,
		*itsStrikes,
		GetCallPut(),
		dumStates);
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price swaption, please advise");
  
    return (*Price)[0];
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwaplegArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaIRSwaplegArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaIRSwaplegArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

