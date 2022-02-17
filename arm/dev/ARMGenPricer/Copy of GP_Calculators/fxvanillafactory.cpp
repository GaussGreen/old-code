/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file fxVanillaFactory.cpp
 *
 *  \brief fx vanilla factory
 *
 *	\author  K. Belkheir
 *	\version 1.0
 *	\date february 2007
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/singleton.h"

#include "gpcalculators/fxvanillafactory.h"
#include "gpcalculators/forexvanilla.h"
#include "gpmodels/fxname.h"
#include "gpcalculators/typedef.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class   : ARM_FXVanillaFactoryImp
///	Routine : CreateFXVanillaModel
///	Returns : void
///	Action  : creates the corresponding FXVanilla and 
///		GaussReplic
////////////////////////////////////////////////////

ARM_GaussReplic2D* ARM_FXVanillaFactoryImp::CreateFXVanillaAndGaussReplic( const string& payName,
		const string& fx1Name,
		int callPut,
		double strike,
		ARM_VanillaType vanillaType,
		const string& fx2Name ,
		double rho ,
		double alpha ,
		double beta )
{
	//Get the currencies for the first FX
	string leftCcy1	= fx1Name.substr(0,3);
	string rightCcy1 = fx1Name.substr(3,3);

	//Get the currencies for the second FX
	string leftCcy2		= fx2Name.substr(0,3);
	string rightCcy2 = fx2Name.substr(3,3);

	// Set the comCcy, the DoubleIntegrationType, and the index of the FX in the integration ("FORDOM" or "DOMFOR")
	string comCcy, Gindex1, Gindex2;
	ARM_FXVanilla2D::FXVanilla2D  fxvanilla2DType;
	if( leftCcy1==leftCcy2 )
	{
		comCcy = leftCcy1;
		fxvanilla2DType = ARM_FXVanilla2D::InvFx1_InvFx2;
		Gindex1 = rightCcy1 + leftCcy1;
		Gindex2 = rightCcy2 + leftCcy2;
	}
	else if( leftCcy1==rightCcy2 )
	{
		comCcy=leftCcy1;
		fxvanilla2DType = ARM_FXVanilla2D::InvFx1_Fx2;
		Gindex1=rightCcy1+leftCcy1;
		Gindex2=leftCcy2+rightCcy2;
	}
	else if( rightCcy1==leftCcy2 )
	{
		comCcy=rightCcy1;
		fxvanilla2DType = ARM_FXVanilla2D::Fx1_InvFx2;
		Gindex1=leftCcy1+rightCcy1;
		Gindex2=rightCcy2+leftCcy2;
	}
	else if( rightCcy1==rightCcy2 )
	{
		comCcy=rightCcy1;
		fxvanilla2DType = ARM_FXVanilla2D::Fx1_Fx2;
		Gindex1=leftCcy1+rightCcy1;
		Gindex2=leftCcy2+rightCcy2;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+": The two Fx must have a common currency");
	
	if( (leftCcy1==leftCcy2)&(rightCcy1==rightCcy2) )
		comCcy=rightCcy1;	//Particular case where there are only two currencies, ie the 2 FX are the same

	ARM_FXName fxName1(fx1Name),fxName2(fx2Name);
	string mktFx1Name = fxName1.GetMktName();
	string mktFx2Name = fxName2.GetMktName();

	//Determination of the type of quanto(Quanto1, Quanto2 or None) 	
	ARM_GaussReplic2D::QuantoType quantoType;
	if( payName != comCcy)
	{
		string quanto = payName + comCcy;
		ARM_FXName quantoName(quanto);
		string mktFxQuantoName = quantoName.GetMktName();
		
		if( mktFxQuantoName==mktFx1Name )
			quantoType = ARM_GaussReplic2D::Quanto1;
		else if( mktFxQuantoName==mktFx2Name )
			quantoType = ARM_GaussReplic2D::Quanto2;
	}
	else
		quantoType = ARM_GaussReplic2D::None;

	//FXVanillaType
	ARM_FXVanilla* opt = NULL ;
	ARM_FXName fxName(Gindex1);
	bool isInvG1 = fxName.GetIsInvMkt();
	fxName = ARM_FXName(Gindex2);
	bool isInvG2 = fxName.GetIsInvMkt();

	if(isInvG1 != isInvG2) rho = -rho;
	if ( vanillaType == ARM_FXVanillaType::spread )
		opt =  new ARM_FXSpread( strike, callPut,fxvanilla2DType,alpha, beta,isInvG1, isInvG2);
	else if ( vanillaType == ARM_FXVanillaType::vanilla )//vanilla in the final version
		opt = new ARM_FXCall( strike, callPut, fxvanilla2DType, isInvG1, isInvG2);
	else if ( vanillaType == ARM_FXVanillaType::basket )
	{
		opt =  new ARM_FXSpread( strike, callPut,fxvanilla2DType,alpha, beta,isInvG1, isInvG2);
		/*if (fxvanilla2DType == ARM_FXVanilla2D::InvFx1_InvFx2)
			opt = new ARM_FXBasketInvFx1InvFx2(mktFx1Name, mktFx2Name, comCcy, strike, callput, quantoType, isInvG1, isInvG2);
		else if (fxvanilla2DType == ARM_FXVanilla2D::InvFx1_Fx2)
			opt =  new ARM_FXBasketInvFx1Fx2(mktFx1Name, mktFx2Name, comCcy, strike, callput, quantoType, isInvG1, isInvG2);
		else if (fxvanilla2DType == ARM_FXVanilla2D::Fx1_InvFx2)
			opt =  new ARM_FXBasketFx1InvFx2(mktFx1Name, mktFx2Name, comCcy, strike, callput, quantoType, isInvG1, isInvG2);
		else if (fxvanilla2DType == ARM_FXVanilla2D::Fx1_Fx2)
			opt =  new ARM_FXBasket(mktFx1Name, mktFx2Name, comCcy, strike, callput, quantoType, isInvG1, isInvG2);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Error in fxvanilla2DType");*/
	}
	else if ( vanillaType == ARM_FXVanillaType::digit )
		opt = new ARM_FXDigital( strike, callPut, fxvanilla2DType, isInvG1, isInvG2);
	else if ( vanillaType == ARM_FXVanillaType::quotient )
		opt =  new ARM_FXQuotient( strike, callPut,fxvanilla2DType,isInvG1, isInvG2);
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Vanilla non computed yet");

	//Gauss Legendre
	ARM_GaussReplic2D* gReplic = new ARM_GaussReplic2D(opt, rho, quantoType);
	delete opt;
	return gReplic;
}

ARM_SingletonHolder<ARM_FXVanillaFactoryImp> ARM_FXVanillaFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

