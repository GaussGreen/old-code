/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillairfixedleg.cpp
 *
 *  \brief vanilla interest rates fixed leg
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#include "gpcalib/vanillairfixedleg.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/checkarg.h"


/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaIRFixedLegArg::CopyNoCleanUp(const ARM_VanillaIRFixedLegArg& rhs)
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


ARM_VanillaIRFixedLegArg::ARM_VanillaIRFixedLegArg(const ARM_VanillaIRFixedLegArg& arg)
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

ARM_VanillaIRFixedLegArg& ARM_VanillaIRFixedLegArg::operator=(const ARM_VanillaIRFixedLegArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaIRFixedLegArg::~ARM_VanillaIRFixedLegArg()
{
	CleanUp();
}
	
////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaIRFixedLegArg::CleanUp()
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
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: Clone
///	Returns: 
///	Action : Clone the object
////////////////////////////////////////////////////

ARM_Object* ARM_VanillaIRFixedLegArg::Clone() const
{
	return new ARM_VanillaIRFixedLegArg(*this);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaIRFixedLegArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: Price
///	Returns: 
///	Action : price a swaption with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaIRFixedLegArg::Price(ARM_PricingModel* model) const
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
///	Struct : ARM_VanillaIRFixedLegArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaIRFixedLegArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaIRFixedLegArg"; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

