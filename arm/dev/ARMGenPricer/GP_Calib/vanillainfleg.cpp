/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file vanillainfleg.cpp
 *
 *  \brief vanilla interest rates leg
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#include "gpcalib/vanillainfleg.h"

/// gpbase
#include "gpbase/gpvector.h"


/// gpinfra
#include "gpinfra/pricingfunctioninflation.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Struct : ARM_VanillaInfSwaplegArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_VanillaInfSwaplegArg::CopyNoCleanUp(const ARM_VanillaInfSwaplegArg& rhs)
{
	itsSpread = rhs.itsSpread;
	itsNumDateStrip = rhs.itsNumDateStrip;
	itsDenomDateStrip = rhs.itsDenomDateStrip;
	itsSwapType = rhs.itsSwapType;
}


ARM_VanillaInfSwaplegArg::ARM_VanillaInfSwaplegArg(	const string& curveName,const string& infCurveName, double evalTime,int swapType,int CallPut,const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double expiry,double tenor,	double Spread,int DecompFreq )
: ARM_VanillaArg(curveName, evalTime, CallPut, expiry, 0, Spread, DecompFreq), itsSpread(Spread), itsNumDateStrip(numDateStrip), itsDenomDateStrip(denomDateStrip), itsSwapType(swapType), itsIRCurveName(curveName), itsInfCurveName(infCurveName)
{
}


ARM_VanillaInfSwaplegArg::ARM_VanillaInfSwaplegArg(const ARM_VanillaInfSwaplegArg& rhs)
:	
	ARM_VanillaArg(rhs), itsSpread(0), itsNumDateStrip(NULL), itsDenomDateStrip(NULL), itsIRCurveName(rhs.itsIRCurveName), itsInfCurveName(rhs.itsInfCurveName)
{
	CopyNoCleanUp(rhs);
}


ARM_VanillaInfSwaplegArg& ARM_VanillaInfSwaplegArg::operator=(const ARM_VanillaInfSwaplegArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
		CleanUp();
        CopyNoCleanUp(rhs);
	}
	return *this;
}

ARM_VanillaInfSwaplegArg::~ARM_VanillaInfSwaplegArg()
{
	CleanUp();
}

ARM_Object* ARM_VanillaInfSwaplegArg::Clone() const
{
	return new ARM_VanillaInfSwaplegArg(*this);
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaInfSwaplegArg
///	Routine: CleanUp
///	Returns: 
///	Action : Delete the various data members
////////////////////////////////////////////////////
void ARM_VanillaInfSwaplegArg::CleanUp()
{
}

////////////////////////////////////////////////////
///	Struct : ARM_VanillaInfSwaplegArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaInfSwaplegArg::Price(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
	size_t modelFactors = model->FactorCount();
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,modelFactors,0) );

	ARM_PricingFuncInflation* InfModel = dynamic_cast<ARM_PricingFuncInflation*>(model);
	ARM_GP_VectorPtr result(NULL);

	if( !InfModel )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Inflation Swap Legs can be priced with InflationModels only!");

	switch( itsSwapType )
	{
	case K_YEARTOYEAR_LEG:
		{
			result = InfModel->YoYSwap( itsIRCurveName, itsInfCurveName, 0, 0, itsSpread, itsNumDateStrip, itsDenomDateStrip, itsDenomDateStrip, itsSpread, dumStates );
			break;
		}
	case K_ZEROCOUPON_LEG:
		{
			result = InfModel->OATSwap( itsIRCurveName, itsInfCurveName, 0, 0, itsSpread, itsNumDateStrip, itsDenomDateStrip, itsDenomDateStrip, itsSpread, dumStates );
			break;
		}
	case K_OATTYPE_LEG:
		{
			result = InfModel->OATSwap( itsIRCurveName, itsInfCurveName, 0, 0, itsSpread, itsNumDateStrip, itsDenomDateStrip, itsDenomDateStrip, itsSpread, dumStates );
			break;
		}
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Unsupported Inflation Leg Type!");
	}

	return (*result)[0];
}



////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaInfSwaplegArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaInfSwaplegArg"; 
}

///////////////////////////////////////////////////////////////////////////////////
//// InflationCap
///////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Struct : ARM_VanillaInfCapArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor. mean stuff
///	Returns: 
///	Action :
////////////////////////////////////////////////////

ARM_VanillaInfCapArg::ARM_VanillaInfCapArg(	const string& curveName,const string& infCurveName, double evalTime,int swapType,int CallPut,double Strike,const ARM_DateStripPtr& numDateStrip,const ARM_DateStripPtr& denomDateStrip,double expiry,double tenor,	double Spread,int DecompFreq )
:   ARM_VanillaInfSwaplegArg(curveName,infCurveName, evalTime,swapType,CallPut,numDateStrip,denomDateStrip,expiry,tenor,Spread,DecompFreq), itsStrike(Strike)
{
}


ARM_VanillaInfCapArg::ARM_VanillaInfCapArg(const ARM_VanillaInfCapArg& rhs)
:	ARM_VanillaInfSwaplegArg(rhs), itsStrike(rhs.itsStrike)
{
}


ARM_VanillaInfCapArg& ARM_VanillaInfCapArg::operator=(const ARM_VanillaInfCapArg& rhs)
{
	if( this != & rhs )
	{
		ARM_VanillaArg::operator=(rhs);
	}
	return *this;
}

ARM_VanillaInfCapArg::~ARM_VanillaInfCapArg()
{
}

ARM_Object* ARM_VanillaInfCapArg::Clone() const
{
	return new ARM_VanillaInfCapArg(*this);
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaInfCapArg
///	Routine: Price
///	Returns: 
///	Action : price a cap with a model checking that it
///				is derived from an interest rate model
////////////////////////////////////////////////////
double ARM_VanillaInfCapArg::Price(ARM_PricingModel* model) const
{
    /// create a dumState to avoid dummy price..
	size_t modelFactors = model->FactorCount();
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,modelFactors,0) );

	ARM_PricingFuncInflation* InfModel = dynamic_cast<ARM_PricingFuncInflation*>(model);
	ARM_GP_VectorPtr result(NULL);

	if( !InfModel )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Inflation Swap Legs can be priced with InflationModels only!");

	switch( getSwapType() )
	{
	case K_YEARTOYEAR_LEG:
		{
			result = InfModel->YoYCapFloor( getIRCurveName(), getInfCurveName(), 0, itsStrike, getSpread(), GetCallPut(), getNumDateStrip(), getDenomDateStrip(), getSpread(), dumStates );
			break;
		}
	case K_OATTYPE_LEG:
		{
			result = InfModel->OATCapFloor( getIRCurveName(), getInfCurveName(), 0, itsStrike, getSpread(), GetCallPut(), getNumDateStrip(), getDenomDateStrip(), getSpread(), dumStates );
			break;
		}
	case K_ZEROCOUPON_LEG:
		{
			result = InfModel->OATCapFloor( getIRCurveName(), getInfCurveName(), 0, itsStrike, getSpread(), GetCallPut(), getNumDateStrip(), getDenomDateStrip(), getSpread(), dumStates );
			break;
		}
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Unsupported Inflation Leg Type!");
	}

	return (*result)[0];
}



////////////////////////////////////////////////////
///	Struct : ARM_VanillaCapArg
///	Routine: toString
///	Returns: string
///	Action : stringify the object to give details about it
////////////////////////////////////////////////////
string ARM_VanillaInfCapArg::toString(const string& indent, const string& nextIndent) const
{ 
	return "ARM_VanillaInfCapArg"; 
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/