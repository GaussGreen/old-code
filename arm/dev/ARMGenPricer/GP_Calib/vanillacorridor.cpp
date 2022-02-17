/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillacorridor.cpp
 *
 *  \brief vanilla corridor 
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpcalib/vanillacorridor.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodelir.h"

/// ARM Kernel
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaCorridorLegArg::CopyNoCleanUp(const ARM_VanillaCorridorLegArg& rhs)
{
#if defined(__GP_STRICT_VALIDATION)
	ThrowErrorOnNullObject( "ResetTimes",			rhs.itsResetTimes           );
	ThrowErrorOnNullObject(  "StartTimes",			rhs.itsStartTimes           );
	ThrowErrorOnNullObject(  "EndTimes",		    rhs.itsEndTimes             );
	ThrowErrorOnNullObject(  "PayTimes",	 	    rhs.itsPayTimes             );
    ThrowErrorOnNullObject(  "FwdPaymentPeriod",	rhs.itsFwdPaymentPeriod     ); 
    ThrowErrorOnEmptyVector( "DownBarrier",			rhs.itsDownBarrier          );
    ThrowErrorOnEmptyVector( "UpBarrier",			rhs.itsUpBarrier            );
    ThrowErrorOnNullObject(  "Nominals",		    rhs.itsNominals             );
    ThrowErrorOnNullObject(  "CouponMargin",	    rhs.itsCouponMargin         );
    ThrowErrorOnEmptyVector( "RefIdxResettimes",	rhs.itsRefIdxResettimes     );
    ThrowErrorOnEmptyVector( "RefIdxStarttimes",	rhs.itsRefIdxStarttimes     );
    ThrowErrorOnEmptyVector( "RefIdxEndtimes",		rhs.itsRefIdxEndtimes       );
    ThrowErrorOnEmptyVector( "RefFwdPeriods",		rhs.itsRefFwdPeriods        );
    ThrowErrorOnEmptyVector( "RefIndexWeight",		rhs.itsRefIndexWeight       );
#endif

	itsResetTimes       =   (ARM_GP_Vector*)   rhs.itsResetTimes->Clone();
	itsStartTimes       =   (ARM_GP_Vector*)   rhs.itsStartTimes->Clone();
	itsEndTimes         =   (ARM_GP_Vector*)   rhs.itsEndTimes->Clone();
    itsPayTimes         =   (ARM_GP_Vector*)   rhs.itsPayTimes->Clone();
    itsIndexPaymentType =   rhs.itsIndexPaymentType;
    itsFwdPaymentPeriod =   (ARM_GP_Vector*)   rhs.itsFwdPaymentPeriod->Clone();
    
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsDownBarrier,itsDownBarrier);
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsUpBarrier,  itsUpBarrier);

	itsCouponMargin     =   (ARM_GP_Vector*)   rhs.itsCouponMargin->Clone();
  	itsNominals         =   (ARM_GP_Vector*)   rhs.itsNominals->Clone();

    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsRefIdxResettimes,itsRefIdxResettimes);
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsRefIdxStarttimes,itsRefIdxStarttimes);
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsRefIdxEndtimes,  itsRefIdxEndtimes);
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsRefFwdPeriods,   itsRefFwdPeriods);
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(rhs.itsRefIndexWeight,  itsRefIndexWeight);
}

ARM_VanillaCorridorLegArg::ARM_VanillaCorridorLegArg(const ARM_VanillaCorridorLegArg& arg)
: 
	ARM_VanillaArg(arg),
	itsResetTimes(NULL), 
	itsStartTimes(NULL), 
	itsEndTimes(NULL),
	itsPayTimes(NULL), 
    itsIndexPaymentType(IDXFIXED),
    itsFwdPaymentPeriod(NULL),

    itsDownBarrier(0), 
	itsUpBarrier(0), 
    itsNominals(NULL), 
	itsCouponMargin(NULL), 

    itsRefIdxStarttimes(0),
    itsRefIdxResettimes(0),
    itsRefIdxEndtimes(0),
    itsRefFwdPeriods(0),
    itsRefIndexWeight(0)
{
    CopyNoCleanUp(arg);
}

ARM_VanillaCorridorLegArg& ARM_VanillaCorridorLegArg::operator=(const ARM_VanillaCorridorLegArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaCorridorLegArg::~ARM_VanillaCorridorLegArg()
{
	CleanUp();
}
	
////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaCorridorLegArg::CleanUp()
{
	delete itsResetTimes;
	delete itsStartTimes;
    delete itsEndTimes;
    delete itsPayTimes; 
	delete itsFwdPaymentPeriod;
    
    DeletePointorVector<ARM_GP_Vector>(itsDownBarrier);
    DeletePointorVector<ARM_GP_Vector>(itsUpBarrier);
    delete itsCouponMargin;
    delete itsNominals;

    DeletePointorVector<ARM_GP_Vector>(itsRefIdxStarttimes);
    DeletePointorVector<ARM_GP_Vector>(itsRefIdxResettimes);
    DeletePointorVector<ARM_GP_Vector>(itsRefIdxEndtimes);
    DeletePointorVector<ARM_GP_Vector>(itsRefFwdPeriods);
    DeletePointorVector<ARM_GP_Vector>(itsRefIndexWeight);
};

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArg
///	Routine: ImpliedVol
///	Returns: Exeption
///	Action :  No Implied Volatility 
////////////////////////////////////////////////////
double ARM_VanillaCorridorLegArg::ImpliedVol(ARM_PricingModel* model) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : No formula is valid to calculate Implied Volatilty";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaCorridorLegArg::Price(ARM_PricingModel* model) const
{
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );
	ARM_VectorPtr Price;
	ARM_PricingFunctionIR* IRModel = dynamic_cast<ARM_PricingFunctionIR*>(model);

    double capPrice = 0.0;
    if(IRModel)
    {
        for(size_t i=0; i<itsResetTimes->size(); ++i)
        {
			/// in order to communicate the current index of a closed form
			/// we use the dumStates to dump the current index!
			/// currently inactivated!

			dumStates->SetModelState(0,0,i);
            Price =IRModel->VanillaCorridorletScalar(
				GetCurveName(),
		        GetEvalTime(),
		        (*itsPayTimes)[i],
                (*itsResetTimes)[i],
                (*itsStartTimes)[i],
		        (*itsEndTimes)[i],
                itsIndexPaymentType,
                (*itsFwdPaymentPeriod)[i],
                *itsRefIdxResettimes[i],
                *itsRefIdxStarttimes[i],                
                *itsRefIdxEndtimes[i],
                *itsRefFwdPeriods[i],
                *itsRefIndexWeight[i],
                (*itsCouponMargin)[i],
                *itsDownBarrier[i],
                *itsUpBarrier[i],
		        (*itsNominals)[i], 		        
		        GetCallPut(),
		        dumStates);

            capPrice += (*Price)[0];
        }
    }
    else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Model is not an interest rate model: not derived from ARM_PricingFunctionIR. So cannot price CorridorLeg, please advise");
	
    return capPrice;
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaCorridorLegArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaCorridorLegArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaCorridorLegArg"; 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

