/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairswap.cpp
 *
 *  \brief vanilla interest rates swap
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#include "gpcalib/vanillairswap.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"


/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwapArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaIRSwapArg::CopyNoCleanUp(const ARM_VanillaIRSwapArg& rhs)
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


ARM_VanillaIRSwapArg::ARM_VanillaIRSwapArg(const ARM_VanillaIRSwapArg& arg)
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

ARM_VanillaIRSwapArg& ARM_VanillaIRSwapArg::operator=(const ARM_VanillaIRSwapArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaIRSwapArg::~ARM_VanillaIRSwapArg()
{
	CleanUp();
}

ARM_Object* ARM_VanillaIRSwapArg::Clone() const
{
	return new ARM_VanillaIRSwapArg(*this);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwapArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaIRSwapArg::CleanUp()
{
	delete itsFixPayTimes;
	delete itsFixPayPeriods;
	delete itsFloatResetTimes;
	delete itsFloatStartTimes;
	delete itsFloatEndTimes;
	delete itsFloatIntTerms; 
    delete itsStrikes;

};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwapArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaIRSwapArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRSwapArg
///	Routine: Price
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaIRSwapArg::Price(ARM_PricingModel* model) const
{
	/// force to use a real dumStates that is not null!
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);
    ARM_VectorPtr Price;
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    
	ARM_GP_Vector FixNominalVec( itsFixPayTimes->size(),itsNominal);
    if(IRModel)
        Price = IRModel->VanillaSwaptionScalar(
		GetCurveName(),
		GetEvalTime(),
		GetExpiry(),
		FixNominalVec,
		FixNominalVec,
		itsStartTime,
		itsEndTime,
		*itsFloatResetTimes,
		*itsFloatStartTimes,
		*itsFloatEndTimes,
		*itsFloatIntTerms,
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
///	Struct : ARM_VanillaIRSwapArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaIRSwapArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaIRSwapArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

