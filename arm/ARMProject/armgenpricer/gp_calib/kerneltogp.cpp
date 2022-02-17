/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file kerneltogp.cpp
 *
 *  \brief file to convert kernel to gp objects
 *
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/kerneltogp.h"

/// gpbase
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/autocleaner.h"
#include "gpbase/datestrip.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/countedptr.h"

/// gpinfra
#include "gpinfra/pricingmodelir.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/vanillacorridor.h"
#include "gpcalib/vanilladigital.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/vanillafxoption.h"
#include "gpcalib/vanillaequityoption.h"
#include "gpcalib/vanillairleg.h"
#include "gpcalib/vanillaportfolio.h"
#include "gpcalib/vanillainfleg.h"
#include "gpcalib/vanillasumopt.h"
#include "gpcalib/vanillaswaptionsmile.h"
#include "gpcalib/typedef.h"
#include "gpcalib/vanillairfxswaption.h"
#include "gpcalib/vanillastrip.h"

/// gpinflation
//#include "gpinflation/infleg.h"
//#include "gpinflation/infcapfloor.h"
//
/////kernel
//#include <inst/swaption.h>
//#include <inst/armdigital.h>
//#include <inst/corridorleg.h>
//#include <inst/option.h>
//#include <inst/forex.h>
//#include <inst/portfolio.h>
//#include <inst/sumopt.h>
//#include <inst/swaption_smile.h>
//#include <inst/optionportfolio.h>
//#include <inst/stripoption.h>
//#include <inst/stripdigitaloption.h>
//
////include for spreadOption
//#include <inst/spreadoption.h>
#include "gpcalib/vanillaspreadoption.h"
#include "gpbase/countedptr.h"
#include "gpinfra/typedef.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingadviser.h"

#include <ccy/currency.h>


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertSecuritytoArgObject
///	Returns: 
///	Action : convert a security into an object for the generic pricer
////////////////////////////////////////////////////
ARM_VanillaArg* ARM_ConverterFromKernel::ConvertSecuritytoArgObject( ARM_Security* Security, double asOfDate, const string& modelName)
{
	/*/// test the type of the security!
	/// support only 
	if(ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(Security))
    {
		if(ARM_SpreadOption* spreadoption = dynamic_cast<ARM_SpreadOption*>(Security))
		{
			return ARM_ConverterFromKernel::ConvertVanillaSpreadOption( spreadoption, asOfDate, modelName );
		}
		else if(ARM_InfCapFloor* infCapFloor = dynamic_cast<ARM_InfCapFloor*>(Security))
			return ARM_ConverterFromKernel::ConvertVanillaInfCapFloor( infCapFloor, asOfDate, modelName );

		return ARM_ConverterFromKernel::ConvertVanillaCapFloor( capFloor, asOfDate, modelName );
	}
	/// Has to come first because smiledSwaption derives from ARM_Swaption 
	else if (ARM_SmiledSwaption* smiledSwaption = dynamic_cast<ARM_SmiledSwaption*>(Security))
	{
		return ARM_ConverterFromKernel::ConvertVanillaSmiledSwaption( smiledSwaption, asOfDate, modelName );
	}

    else if(ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(Security))
    {
		return ARM_ConverterFromKernel::ConvertVanillaSwaption( swaption, asOfDate, modelName );
    }
    else if(ARM_Digital* digital = dynamic_cast<ARM_Digital*>(Security))
    {
        return ARM_ConverterFromKernel::ConvertVanillaDigital( digital, asOfDate, modelName );    
    }
    else if(ARM_CorridorLeg* corridorleg = dynamic_cast<ARM_CorridorLeg*>(Security))
    {
        return ARM_ConverterFromKernel::ConvertVanillaCorridorLeg( corridorleg, asOfDate, modelName );    
    }
	else if( ARM_Option* option = dynamic_cast<ARM_Option*>(Security))
	{
		ARM_Security* underlying = option->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			return ARM_ConverterFromKernel::ConvertToFXOption( option, asOfDate );
		}
		else
		{
			return ARM_ConverterFromKernel::ConvertToEqOption( option, asOfDate, modelName );
		}
	}
	else if( ARM_InfLeg* infleg = dynamic_cast<ARM_InfLeg*>(Security) )
	{
		return ARM_ConverterFromKernel::ConvertVanillaInfSwapLeg( infleg, asOfDate, modelName );
	}
	else if(ARM_SwapLeg* swapleg = dynamic_cast<ARM_SwapLeg*>(Security))
    {
		return ARM_ConverterFromKernel::ConvertVanillaIRSwapleg( swapleg, asOfDate, modelName );
    }
	else if(ARM_Portfolio* portfolio = dynamic_cast<ARM_Portfolio*>(Security))
    {
		return ARM_ConverterFromKernel::ConvertPortoflio( portfolio, asOfDate, modelName );
    }
	else if(ARM_OptionPortfolio* optPortfolio = dynamic_cast<ARM_OptionPortfolio*>(Security))
    {
		ARM_Portfolio* portfolio = optPortfolio->GetPtf();
		if( portfolio->size()>=2 &&
			dynamic_cast<ARM_Option*>(portfolio->GetAsset(0)) &&
			dynamic_cast<ARM_Swaption*>(portfolio->GetAsset(portfolio->size()-1)) )
			/// Here is a IR/FX hybrid portfolio describing a IR/FX swaption
			return ARM_ConverterFromKernel::ConvertIrFxSwaption( optPortfolio, asOfDate, modelName );
		else
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"No option portfolio currently supported!");
		}

    }
	else if (ARM_SumOpt* sumOpt = dynamic_cast<ARM_SumOpt*>(Security))
	{
		return ARM_ConverterFromKernel::ConvertVanillaSumOpt( sumOpt, asOfDate, modelName );
	}
	else if (ARM_StripOption* stripOption = dynamic_cast<ARM_StripOption*>(Security))
	{
		if (ARM_StripDigitalOption* stripDigitalOption = dynamic_cast<ARM_StripDigitalOption*>(Security))
		{
			return ARM_ConverterFromKernel::ConvertToStripDigitalOption( stripDigitalOption, asOfDate, modelName );
		}
		else
		{
			return ARM_ConverterFromKernel::ConvertToStripOption( stripOption, asOfDate, modelName );
		}
	}
	else 
    {
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Only capfloor,swaption,digital, corridorleg and swapleg currently supported!");
	}*/
	return NULL;
}

void ARM_ConverterFromKernel::TransferToSecurity(ARM_VanillaArg* vanillaArg, ARM_Security* sec)
{
	/// test the type of the security!
	/// support only 
	//if(ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(sec))
 //   {
	//	if(ARM_SpreadOption* spreadoption = dynamic_cast<ARM_SpreadOption*>(sec))
	//	{
	//	}
	//	else if(ARM_InfCapFloor* infCapFloor = dynamic_cast<ARM_InfCapFloor*>(sec))
	//	{
	//	}
	//}
 //   else if(ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(sec))
 //   {
 //   }
 //   else if(ARM_Digital* digital = dynamic_cast<ARM_Digital*>(sec))
 //   {
 //   }
 //   else if(ARM_CorridorLeg* corridorleg = dynamic_cast<ARM_CorridorLeg*>(sec))
 //   {
 //   }
	//else if( ARM_Option* option = dynamic_cast<ARM_Option*>(sec))
	//{
	//	ARM_Security* underlying = option->GetUnderlying();
	//	if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
	//	{
	//	}
	//	else
	//	{
	//	}
	//}
	//else if( ARM_InfLeg* infleg = dynamic_cast<ARM_InfLeg*>(sec) )
	//{
	//}
	//else if(ARM_SwapLeg* swapleg = dynamic_cast<ARM_SwapLeg*>(sec))
 //   {
 //   }
	//else if(ARM_Portfolio* portfolio = dynamic_cast<ARM_Portfolio*>(sec))
 //   {
 //   }
	//else if (ARM_SumOpt* sumOpt = dynamic_cast<ARM_SumOpt*>(sec))// && ())
	//{
	//	if (ARM_VanillaSumOptArg* sumOptArg = dynamic_cast<ARM_VanillaSumOptArg*>(vanillaArg))
	//		ARM_ConverterFromKernel::TranferToSumOptSec( sumOptArg, sumOpt );
	//}
	//else if(ARM_StripOption* portfolio = dynamic_cast<ARM_StripOption*>(sec))
 //   {
 //   }
	//else if(ARM_StripDigitalOption* portfolio = dynamic_cast<ARM_StripDigitalOption*>(sec))
 //   {
 //   }
	//else 
 //   {
	//	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
	//	"Only capfloor,swaption,digital, corridorleg and swapleg currently supported!");
	//}
}



////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaCapFloor
///	Returns: 
///	Action : transform an CapFloor into an ARM_VanillaCapArg
////////////////////////////////////////////////////
ARM_VanillaCapArg* ARM_ConverterFromKernel::ConvertVanillaCapFloor( ARM_CapFloor* capFloor, double asOfDate, const string& modelName )
{
	return NULL;
	/// fills the specific ARM_VanillaCapArg arguments
    /*size_t size		= capFloor->GetFlowEndDates()->GetSize();
	int CallPut		= capFloor->GetCapFloorType();
	double Spread	= capFloor->GetSwapLeg()->GetSpread();
	int DecompFreq	= capFloor->GetSwapLeg()->GetDecompFreq();
	std::vector<double>& ResetTimes    = To_pARM_GP_Vector( capFloor->GetResetDates() );
	std::vector<double>& StartTimes    = To_pARM_GP_Vector( capFloor->GetSwapLeg()->GetFwdRateStartDates() );
	std::vector<double>& EndTimes      = To_pARM_GP_Vector( capFloor->GetSwapLeg()->GetFwdRateEndDates() );
	std::vector<double>& PayTimes      = To_pARM_GP_Vector( capFloor->GetPaymentDates() );
	std::vector<double>& PayPeriods    = To_pARM_GP_Vector( capFloor->GetSwapLeg()->GetInterestTerms() );
	std::vector<double>& Nominals      = new std::vector<double>(size);
	std::vector<double>& Strikes       = new std::vector<double>(size);

	bool isAtTheMoneyFlag = false;
	for(size_t i=0;i<size;++i)
	{
		(*Strikes)[i] = (( capFloor->GetStrikes()!=NULL
			? capFloor->GetStrikes()->CptReferenceValue( (*ResetTimes)[i] )
			: capFloor->GetStrike() ) - capFloor->GetSwapLeg()->ComputeSpread(i))
			/ CC_NS(ARM_Constants,rateBase);
		(*ResetTimes)[i] -= asOfDate;
	
		(*StartTimes)[i] -= asOfDate;
		(*EndTimes)[i] -= asOfDate;
		
		(*Nominals)[i] = capFloor->GetAmount()->CptReferenceValue( (*PayTimes)[i] );
		(*PayTimes)[i] -= asOfDate;
	}
	if ((*Strikes)[0]==(-0.01))
		isAtTheMoneyFlag = true;
	/// part on the curve
    string CurveName = modelName;
    double evalTime = 0.0;
	double mktPrice = -1.0;

	return new ARM_VanillaCapArg(
		CurveName,
		evalTime,
		CallPut,
		Nominals,
		ResetTimes,
		StartTimes,
		EndTimes,
		Strikes,
		PayTimes,
		PayPeriods,
		Spread,
		DecompFreq,
		isAtTheMoneyFlag);*/
	return NULL;
}
 
////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaDigital
///	Returns: 
///	Action : transform an ARM_Digital into an ARM_VanillaDigitalArg
////////////////////////////////////////////////////
ARM_VanillaDigitalArg* ARM_ConverterFromKernel::ConvertVanillaDigital( ARM_Digital* digital, double asOfDate, const string& modelName )
{
	return NULL;
	/// fills the specific ARM_VanillaCapArg arguments
    /*size_t size						    = digital->GetSwapLeg()->GetFlowEndDates()->GetSize();
	int CallPut		= digital->GetCapFloor()->GetCapFloorType();
	std::vector<double>& ResetTimes    = To_pARM_GP_Vector( digital->GetResetDates() );
	std::vector<double>& StartTimes    = To_pARM_GP_Vector( digital->GetSwapLeg()->GetFwdRateStartDates() );
	std::vector<double>& EndTimes      = To_pARM_GP_Vector( digital->GetSwapLeg()->GetFwdRateEndDates() );
	std::vector<double>& PayTimes      = To_pARM_GP_Vector( digital->GetSwapLeg()->GetPaymentDates() );
	std::vector<double>& PayPeriods    = To_pARM_GP_Vector( digital->GetSwapLeg()->GetInterestTerms() );
	std::vector<double>& Nominals      = new std::vector<double>(size);
	std::vector<double>& Strikes       = new std::vector<double>(size);

	for(size_t i=0;i<size;++i)
	{
		(*(Strikes))[i] = digital->GetStrikes()->CptReferenceValue((*(ResetTimes))[i])
            / CC_NS(ARM_Constants,rateBase);
		(*(ResetTimes))[i] -= asOfDate;
	
		(*(StartTimes))[i] -= asOfDate;
		(*(EndTimes))[i] -= asOfDate;

        double paytime = (*(PayTimes))[i];
		(*(Nominals))[i] = digital->GetAmount()->CptReferenceValue(paytime) * digital->GetPayoff()->CptReferenceValue(paytime)/100.0;
		(*(PayTimes))[i] -= asOfDate;
	}

	/// part on the curve
    string CurveName = modelName;
    double EvalTime = 0.0;

	return new ARM_VanillaDigitalArg( 
		CurveName,
		EvalTime,
		CallPut,
		Nominals,
		ResetTimes,
		StartTimes,
		EndTimes,
		Strikes,
		PayTimes,
		PayPeriods	);*/
	return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaCorridorLeg
///	Returns: 
///	Action : transform an ARM_CorridorLeg into an ARM_VanillaCorridorLegArg
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArg* ARM_ConverterFromKernel::ConvertVanillaCorridorLeg(ARM_CorridorLeg* corridorleg, double asOfDate, const string& modelName )
{
	
	/// fills the specific ARM_VanillaCapArg arguments
 //   size_t size						= corridorleg->GetPaymentEndDates()->GetSize();
	//int Callput		                = corridorleg->GetRcvOrPay();

	//std::vector<double>& ResetTimes       = To_pARM_GP_Vector( corridorleg->GetResetDates() );
	//std::vector<double>& StartTimes       = To_pARM_GP_Vector( corridorleg->GetPaymentStartDates() );
	//std::vector<double>& EndTimes         = To_pARM_GP_Vector( corridorleg->GetPaymentEndDates() );
	//std::vector<double>& PayTimes         = To_pARM_GP_Vector( corridorleg->GetPaymentDates() );
 //   int indexPaymentType            = corridorleg->GetPaymentIndex()->GetIndexType();
 //   std::vector<double>& fwdPaymentPeriod = To_pARM_GP_Vector( corridorleg->GetPaymentInterestterms() );

 //   ARM_VectorVector DownBarriers;
	//ARM_VectorVector UpBarriers  ;

	//std::vector<double>& Nominals      = new std::vector<double>(size);
	//std::vector<double>& CouponMargin  = new std::vector<double>(size);

 //   size_t nbRefFlows = corridorleg->GetRefIndexFixingDates()->size();

 //   ARM_VectorVector VectRefIdxResettimes;
 //   ARM_VectorVector VectRefIdxStarttimes;
 //   ARM_VectorVector VectRefIdxEndtimes;
 //   ARM_VectorVector VectRefFwdPeriods;
 //   ARM_VectorVector VectRefIndexWeight;

 //   size_t i, k=0;
 //   for(i=0; i< size; ++i )
 //   {
 //       std::vector<double> DownBarrier(0);
	//    std::vector<double> UpBarrier(0)  ;

 //       std::vector<double> RefIdxResettimes(0);
 //       std::vector<double> RefIdxStarttimes(0);
 //       std::vector<double> RefIdxEndtimes(0);
 //       std::vector<double> RefFwdPeriods(0);
 //       std::vector<double> RefIndexWeight(0);
 //       
 //       double epsilon = 0.001;         
 //       if ( corridorleg->GetRefIndex()->GetResetTiming() == K_ADVANCE )
 //          epsilon = 0.0;  // le taux de ref est resette en avance, 
 //                          //on ne prend pas le taux qui fixe a la date de paiement
 //       while ( (*(corridorleg->GetRefIndexFixingDates()))[k] < corridorleg->GetFlowEndDates()->Elt(i) + epsilon )
 //       {  
 //           double resetIndexdate = (*(corridorleg->GetRefIndexFixingDates()))[k];
 //           double barrier = corridorleg->GetDownBarriers()->CptReferenceValue(resetIndexdate)
 //                                                                                   / CC_NS(ARM_Constants,rateBase);
 //           DownBarrier.push_back(barrier);

 //           barrier = corridorleg->GetUpBarriers()->CptReferenceValue(resetIndexdate)
 //                                                                                   / CC_NS(ARM_Constants,rateBase);
 //           UpBarrier.push_back(barrier);

 //           resetIndexdate -= asOfDate;
 //           RefIdxResettimes.push_back(resetIndexdate);

 //           double startIndexdate = (*(corridorleg->GetRefIndexStartDates()))[k] - asOfDate;
 //           RefIdxStarttimes.push_back(startIndexdate);

 //           double endIndexdate = (*(corridorleg->GetRefIndexEndDates()))[k] - asOfDate;
 //           RefIdxEndtimes.push_back(endIndexdate);

 //           double Refperiod = (*(corridorleg->GetRefIndexInterestTerms()))[k] ;
 //           RefFwdPeriods.push_back(Refperiod);

 //           double Weight = corridorleg->GetRefIndexWeight()[k];
 //           RefIndexWeight.push_back(Weight);  
 //           ++k;
 //           if(k > nbRefFlows -1 )
 //               break;
 //       }
 //       DownBarriers.push_back((std::vector<double>&)DownBarrier.Clone());
 //       UpBarriers.push_back((std::vector<double>&)UpBarrier.Clone());

 //       VectRefIdxResettimes.push_back((std::vector<double>&)RefIdxResettimes.Clone());
 //       VectRefIdxStarttimes.push_back((std::vector<double>&)RefIdxStarttimes.Clone());
 //       VectRefIdxEndtimes.push_back((std::vector<double>&)RefIdxEndtimes.Clone());
 //       VectRefFwdPeriods.push_back((std::vector<double>&)RefFwdPeriods.Clone());
 //       VectRefIndexWeight.push_back((std::vector<double>&)RefIndexWeight.Clone());
 //   }

	//for(i=0;i<size;++i)
	//{
 //       double startdate = (*StartTimes)[i];
 //       double paydate = (*PayTimes)[i];

	//	if(corridorleg->GetSpreads())
 //           (*CouponMargin)[i] = corridorleg->GetSpreads()->CptReferenceValue(startdate)
 //                                                                                   / CC_NS(ARM_Constants,rateBase);
 //       else
 //           (*CouponMargin)[i] = corridorleg->GetSpread()/ CC_NS(ARM_Constants,rateBase);

 //       (*Nominals)[i] = corridorleg->GetAmount()->CptReferenceValue(paydate);

	//	(*ResetTimes)[i]  -= asOfDate;	
	//	(*StartTimes)[i]  -= asOfDate;
	//	(*EndTimes)[i]    -= asOfDate;        
	//	(*PayTimes)[i]    -= asOfDate;
	//}

	///// part on the curve
 //   double EvalTime = 0.0;

 //   /// creates an empty ARM_VanillaCapArg
	//ARM_VanillaCorridorLegArg* argVanillaCorridorLeg= new ARM_VanillaCorridorLegArg(
 //           modelName               ,
 //           EvalTime                ,
 //           Callput                 ,
 //           ResetTimes              ,
 //           StartTimes              ,
 //           EndTimes                ,
 //           PayTimes                ,
 //           indexPaymentType        ,
 //           fwdPaymentPeriod        ,
 //           VectRefIdxResettimes    ,
 //           VectRefIdxStarttimes    ,
 //           VectRefIdxEndtimes      ,
 //           VectRefFwdPeriods       ,
 //           VectRefIndexWeight      ,
 //           DownBarriers            , 
 //           UpBarriers              ,             
 //           CouponMargin            ,
 //           Nominals                ); 

	return NULL;//argVanillaCorridorLeg;
}


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaSwaption
///	Returns: 
///	Action : transform an ARM_Swaption into an ARM_VanillaSwaptionArg
////////////////////////////////////////////////////
ARM_VanillaSwaptionArg* ARM_ConverterFromKernel::ConvertVanillaSwaption( ARM_Swaption* swaption, double asOfDate, const string& modelName )
{
	return NULL;
	/// convention is to say that a swaption is a call or put on the swap rate
	/// hence a payer swap is a call on the swap rate while a receiver is a put
	/// beware that if you price the swaption as a bond option, the convention is the opposite
	//int CallPut		= swaption->IsPayer() - swaption->IsReceiver();
	//
	//double ExpiryTime   = swaption->GetExpiryDate().GetJulian() - asOfDate;
	//
	///// --> old : unadjusted start & end dates
	///// double StartTime    = swaption->GetFloatLeg()->GetStartDateNA().GetJulian() - asOfDate;
 //   /// if( StartTime < ExpiryTime)
 //   ///    StartTime    = swaption->GetFloatLeg()->GetStartDate().GetJulian() - asOfDate;
	/////
	///// double EndTime      = swaption->GetFloatLeg()->GetEndDateNA().GetJulian() - asOfDate;
	//
	///// --> new : we now take adjusted start & end dates
	//double StartTime = swaption->GetFloatLeg()->GetFlowStartDates()->Elt(0) - asOfDate;
	//	
	//double EndTime 
	//	= swaption->GetFloatLeg()->GetFlowEndDates()->Elt(swaption->GetFloatLeg()->GetFlowEndDates()->size()-1)
	//	 - asOfDate;

 //   int FixFrequency = swaption->GetFixedLeg()->GetPaymentFreq();
	//int FlaotFrequency = swaption->GetFloatLeg()->GetPaymentFreq();
	//int FixDayCount = swaption->GetFixedLeg()->GetDayCount();
	//int FloatDayCount = swaption->GetFloatLeg()->GetDayCount();


	//
	//std::vector<double>& fixPayTimes = To_pARM_GP_Vector( swaption->GetFixedLeg()->GetFlowEndDates() );
	//if( fixPayTimes )
	//	(*fixPayTimes) -=asOfDate;

	//std::vector<double>& fixPayPeriods = To_pARM_GP_Vector( swaption->GetFixedLeg()->GetInterestTerms() );
	//
	//std::vector<double>& floatResetTimes = To_pARM_GP_Vector( swaption->GetFloatLeg()->GetResetDates() );
	//if( floatResetTimes)
	//	(*floatResetTimes) -=asOfDate;
	//
	//std::vector<double>& floatStartTimes = To_pARM_GP_Vector( swaption->GetFloatLeg()->GetFwdRateStartDates() );
	//if( floatStartTimes )
	//	(*floatStartTimes) -=asOfDate;
	//
	//std::vector<double>& floatEndTimes   = To_pARM_GP_Vector( swaption->GetFloatLeg()->GetFwdRateEndDates() );
	//if( floatEndTimes )
	//	(*floatEndTimes) -=asOfDate;
	//
	//std::vector<double>& floatInterestTerms = To_pARM_GP_Vector( swaption->GetFloatLeg()->GetInterestTerms() );

 //   ARM_SwapLeg* fixLeg = swaption->GetFixedLeg();

	//ARM_Vector* fixResetDates = fixLeg->GetResetDates();
 //   size_t sizefix = fixResetDates? fixResetDates->size(): 0;

	//std::vector<double>& Strikes;
	//if(!fixLeg->GetVarCoupons())
	//	Strikes = new std::vector<double>(sizefix, swaption->GetStrike()/CC_NS(ARM_Constants,rateBase));
 //   else
 //   {
 //       Strikes = new std::vector<double>(sizefix);
 //       for(int i=0; i<sizefix; ++i)
 //           (*Strikes)[i]=(fixLeg->GetVarCoupons()->CptReferenceValue((*fixResetDates)[i])/CC_NS(ARM_Constants,rateBase));
 //   }

	//bool isAtTheMoneyFlag = false;
	//if ((*Strikes)[0]==(-10000))
	//	isAtTheMoneyFlag = true;

	//string CurveName = modelName;
 //   double EvalTime	= 0.0;

	////// Handle Variable Notional
	////// Calculate a Notional Vector for the floatLeg and the FixLeg
	//int nbFxFlows = fixPayTimes->size();
	//int nbFloatFlows = floatResetTimes->size();

	//std::vector<double>&  FixNotional =new std::vector<double>(nbFxFlows);
	//std::vector<double>&  FloatNotional = new std::vector<double>(nbFloatFlows);
	//for (int i = 0; i < nbFxFlows; i++) 
 //   {
 //       (*FixNotional)[i] = swaption->GetFixedLeg()->GetAmount()->CptReferenceValue(swaption->GetFixedLeg()->GetPaymentDates()->Elt(i));        
 //   }
	//for (i = 0; i < nbFloatFlows; i++) 
 //   {
 //       (*FloatNotional)[i] = swaption->GetFloatLeg()->GetAmount()->CptReferenceValue(swaption->GetFloatLeg()->GetPaymentDates()->Elt(i));        
 //   }

	//bool isConstantNominal = true;
	//bool isConstantStrike = true;
	//bool isConstantSpread = true;

	//if (!(swaption->GetFixedLeg()->GetAmount()->IsConstant()) || !(swaption->GetFloatLeg()->GetAmount()->IsConstant()))
	//	isConstantNominal = false;
	//if ((swaption->GetStrikes()) && (!(swaption->GetStrikes()->IsConstant())))
	//	isConstantStrike = false;
	//

	//std::vector<double>& MarginVector ;
	//if ((swaption->GetFloatLeg()->GetVariableSpread()) && !(swaption->GetFloatLeg()->GetVariableSpread()->IsConstant()))
	//{
	//	isConstantSpread = false;
	//	MarginVector = new std::vector<double>(floatResetTimes->size());
	//	for (i = 0; i < nbFloatFlows; i++) 
	//	{
	//		(*MarginVector)[i] = swaption->GetFloatLeg()->GetVariableSpread()->CptReferenceValue((*swaption->GetFloatLeg()->GetFwdRateStartDates())[i])/CC_NS(ARM_Constants,rateBase);
	//	}
	//}
	//else
	//	MarginVector = new std::vector<double>(floatResetTimes->size(),swaption->GetFloatLeg()->GetSpread()/CC_NS(ARM_Constants,rateBase));

	//return new ARM_VanillaSwaptionArg(
	//	CurveName,
	//	EvalTime,
	//	CallPut,
	//	MarginVector,
	//	FixNotional,
	//	FloatNotional,
	//	ExpiryTime,
	//	StartTime,
	//	EndTime,
	//	Strikes,
	//	fixPayTimes,
	//	fixPayPeriods,
	//	floatResetTimes,
	//	floatStartTimes,
	//	floatEndTimes,
	//	floatInterestTerms,
 //       FixFrequency,
	//	FlaotFrequency,
	//	FixDayCount,
	//	FloatDayCount,
	//	asOfDate,
	//	isConstantNominal,
	//	isConstantSpread,
	//	isConstantStrike,
	//	swaption->GetGenSecString(),
	//	isAtTheMoneyFlag);
return NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToFXOption
///	Returns: 
///	Action : transform an Option into an ARM_VanillaFXOption
////////////////////////////////////////////////////
ARM_VanillaFxOption* ARM_ConverterFromKernel::ConvertToFXOption( ARM_Option* option,double asOfDate)
{
	/*double maturityTime = option->GetMaturity().GetJulian()-asOfDate;
	double fwdTime = option->GetForwardDate().GetJulian()-asOfDate;
    double payTime = fwdTime;
    if(!option->IsStdPayDate())
        payTime = option->GetPayDate().GetJulian()-asOfDate;
	double strike = option->GetStrike();
    int callPut = option->IsCall()? K_CALL : K_PUT;

	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(option->GetUnderlying());
	string forCcy = forex->GetMainCurrency()->GetCcyName();
	string domCcy = forex->GetMoneyCurrency()->GetCcyName();
	string name(forCcy + domCcy);

	return new ARM_VanillaFxOption(name, maturityTime, fwdTime, payTime, 
                                   strike, callPut, option->GetFixedFx());*/
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToStripOption
///	Returns: 
///	Action : transform an Strip Option into an ARM_VanillaStripOption
////////////////////////////////////////////////////
ARM_VanillaStripArg* ARM_ConverterFromKernel::ConvertToStripOption( ARM_StripOption* stripOption,double asOfDate, const string& modelName )
{
	/*vector<ARM_Option*> strip = stripOption->GetStrip();
	ARM_GP_T_Vector<ARM_VanillaArg*> vanillas;
	std::vector<double> coeffs(strip.size());

	std::vector<double>& startDates = stripOption->GetSchedule()->GetFlowStartDates();
	std::vector<double>& payDates   = stripOption->GetSchedule()->GetPaymentDates();
	std::vector<double>& terms      = stripOption->GetSchedule()->GetInterestTerms();

	size_t i;

	for(i = 0; i < strip.size(); ++i)
	{
		ARM_Security* underlying = strip[i]->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			ARM_VanillaFxOption* fxOption = ARM_ConverterFromKernel::ConvertToFXOption( strip[i], asOfDate );

            fxOption->SetDiscPricingMode(stripOption->GetDiscPricingMode());


			ARM_Currency* mainCcy = (((ARM_Forex *) strip[i]->GetUnderlying())->GetMainCurrency());
			
			if (*mainCcy == stripOption->GetPaymentCcy())
				fxOption->SetPayFgnCcy(true);

			vanillas.push_back(fxOption);
		}
		else
		{
			vanillas.push_back(ARM_ConverterFromKernel::ConvertToEqOption( strip[i], asOfDate, modelName ));
		}

		int PorSale = stripOption->GetPorS();

		coeffs[i] += (*terms)[i] * stripOption->GetNotional().Interpolate((*payDates)[i]-asOfDate);
		coeffs[i] *= PorSale * stripOption->GetLeverage().Interpolate((*startDates)[i]-asOfDate);
	}

    ARM_VanillaStripArg* createdVanillaStrip = new ARM_VanillaStripArg(vanillas, coeffs);
    
    createdVanillaStrip->SetDiscPricingMode(stripOption->GetDiscPricingMode());

	return(createdVanillaStrip);;*/
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToStripDigitalOption
///	Returns: 
///	Action : transform an Strip Digital Option into an ARM_VanillaStrip
////////////////////////////////////////////////////
ARM_VanillaStripArg* ARM_ConverterFromKernel::ConvertToStripDigitalOption( ARM_StripDigitalOption* stripDigitalOption,double asOfDate, const string& modelName )
{
	/*vector<ARM_Option*> stripUp = stripDigitalOption->GetStripUp();
	vector<ARM_Option*> stripDown = stripDigitalOption->GetStripDown();
	ARM_GP_T_Vector<ARM_VanillaArg*> vanillas;
	std::vector<double> coeffs(stripUp.size());
	std::vector<double>& resetDates   = stripDigitalOption->GetSchedule()->GetResetDates();
	std::vector<double>& startDates   = stripDigitalOption->GetSchedule()->GetFlowStartDates();
	std::vector<double>& payDates   = stripDigitalOption->GetSchedule()->GetPaymentDates();
	std::vector<double>& terms      = stripDigitalOption->GetSchedule()->GetInterestTerms();

	size_t i;
	for (i = 0; i < stripUp.size(); i++)
	{
		coeffs[i] += stripDigitalOption->GetStrikesUp().Interpolate((*startDates)[i]-asOfDate);
		coeffs[i] -= stripDigitalOption->GetStrikesDown().Interpolate((*startDates)[i]-asOfDate);
		coeffs[i] = 1/coeffs[i];
		coeffs[i] *= stripDigitalOption->GetPayOff().Interpolate((*startDates)[i]-asOfDate)*0.01;
		coeffs[i] *= stripDigitalOption->GetNotional().Interpolate((*payDates)[i]-asOfDate);

		int PorSale = stripDigitalOption->GetPorS();

		coeffs[i] *= PorSale*(*terms)[i];
		coeffs[i] *= stripDigitalOption->GetLeverage().Interpolate((*startDates)[i]-asOfDate);
	}

	std::vector<double> optCoeffs;
	for(i = 0; i < stripUp.size(); ++i)
	{
		ARM_Security* underlying = stripUp[i]->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			ARM_VanillaFxOption* fxOption = ARM_ConverterFromKernel::ConvertToFXOption( stripUp[i], asOfDate );

            fxOption->SetDiscPricingMode(stripDigitalOption->GetDiscPricingMode());

			ARM_Currency* mainCcy = (((ARM_Forex *) stripUp[i]->GetUnderlying())->GetMainCurrency());
			
			if (*mainCcy == stripDigitalOption->GetPaymentCcy())
				fxOption->SetPayFgnCcy(true);

			vanillas.push_back(fxOption);
		}
		else
		{
			vanillas.push_back(ARM_ConverterFromKernel::ConvertToEqOption( stripUp[i], asOfDate, modelName ));
		}
		int optType=(stripUp[i]->IsCall()?1:-1);
		optCoeffs.push_back(-coeffs[i]*optType);
	}

	for(i = 0; i < stripDown.size(); ++i)
	{
		ARM_Security* underlying = stripDown[i]->GetUnderlying();
		if( ARM_Forex* forex = dynamic_cast<ARM_Forex*>(underlying) )
		{
			ARM_VanillaFxOption* fxOption = ARM_ConverterFromKernel::ConvertToFXOption( stripDown[i], asOfDate );

            fxOption->SetDiscPricingMode(stripDigitalOption->GetDiscPricingMode());

			ARM_Currency* mainCcy = (((ARM_Forex *) stripDown[i]->GetUnderlying())->GetMainCurrency());
			
			if (*mainCcy == stripDigitalOption->GetPaymentCcy())
				fxOption->SetPayFgnCcy(true);

			vanillas.push_back(fxOption);
		}
		else
		{
			vanillas.push_back(ARM_ConverterFromKernel::ConvertToEqOption( stripDown[i], asOfDate, modelName ));
		}
		int optType=(stripDown[i]->IsCall()?1:-1);
		optCoeffs.push_back(+coeffs[i]*optType);
	}

    ARM_VanillaStripArg* createdVanillaDigitalStrip = new ARM_VanillaStripArg( vanillas, optCoeffs);

    createdVanillaDigitalStrip->SetDiscPricingMode(stripDigitalOption->GetDiscPricingMode());

	return(createdVanillaDigitalStrip);*/
	return NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToFXOption
///	Returns: 
///	Action : transform an Option into an ARM_VanillaFXOption
////////////////////////////////////////////////////
ARM_VanillaEqOption* ARM_ConverterFromKernel::ConvertToEqOption( ARM_Option* option, double asOfDate, const string& modelName ){
	return NULL;
}/*
{
	double maturityTime = option->GetMaturity().GetJulian()-asOfDate;
	double fwdTime = option->GetForwardDate().GetJulian()-asOfDate;
    double payTime = fwdTime;
    if(!option->IsStdPayDate())
        payTime = option->GetPayDate().GetJulian()-asOfDate;
	double strike = option->GetStrike();
    int callPut = option->IsCall()? K_CALL : K_PUT;

	return new ARM_VanillaEqOption( modelName, maturityTime, fwdTime, payTime, strike, callPut );
}*/


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaSpreadOption
///	Returns: 
///	Action : transform an ARM_SpreadOption into an ARM_VanillaSpreadOptionArg
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArg* ARM_ConverterFromKernel::ConvertVanillaSpreadOption( ARM_SpreadOption* spreadoption, double asOfDate, const string& modelName )
{
	return NULL;
}
//{
//	int dayCount = spreadoption->GetSwapLeg()->GetDayCount();
//	int nbCapLet = spreadoption->GetNumFlows();
//	int CallPut = spreadoption->IsCap()-spreadoption->IsFloor();
//	string curveName = modelName;
//	double evalTime	= 0.0;
//	double expiryTime = evalTime;
//	double startTime = spreadoption->GetStartDate().GetJulian()-asOfDate;
//	double endTime = spreadoption->GetSwapLeg()->GetEndDate().GetJulian()-asOfDate;
//	double coeffShortValue = spreadoption->GetWeight1();
//	double coeffLongValue = spreadoption->GetWeight2();
//	std::vector<double>& coeffShort = new std::vector<double>(nbCapLet,coeffShortValue);
//	std::vector<double>& coeffLong = new std::vector<double>(nbCapLet,coeffLongValue);
//	
//	std::vector<double>& resetTimes		= To_pARM_GP_Vector(spreadoption->GetSwapLeg()->GetResetDates());
//	std::vector<double>& payTimes			= To_pARM_GP_Vector(spreadoption->GetSwapLeg()->GetPaymentDates()); 
//	std::vector<double>& swapShortFloatStartTime	= To_pARM_GP_Vector(spreadoption->GetSpreadLeg()->GetFirstLeg()->GetFwdRateStartDates());
//	std::vector<double>& swapLongFloatStartTime	= To_pARM_GP_Vector(spreadoption->GetSpreadLeg()->GetSecondLeg()->GetFwdRateStartDates());
//	
//					
//	ARM_Vector* theopaymentDates = spreadoption->GetSwapLeg()->GetTheoPayDates();
//	ARM_Vector* flowStartDates =  spreadoption->GetSwapLeg()->GetFlowStartDates();
//	ARM_Vector* flowEndDates = spreadoption->GetSwapLeg()->GetFlowEndDates();
//	ARM_Vector* flowPayDates = spreadoption->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates();
//	
//	std::vector<double>& payPeriods = new std::vector<double>(nbCapLet);
//	std::vector<double>& strikes    = new std::vector<double>(nbCapLet);
//	std::vector<double>& notional   = new std::vector<double>(nbCapLet);
//	std::vector<double>& swapShortFloatEndTime  = new std::vector<double>(nbCapLet);
//	std::vector<double>& swapLongFloatEndTime   = new std::vector<double>(nbCapLet);
//    
//	int tenor_1(0), tenor_2(0), tenor_p(0);
//	bool shortIndexIsLibor (false);
//	bool longIndexIsLibor  (false);
//			
//	if (spreadoption->GetSpreadLeg()->GetFirstLeg()->GetName() == ARM_SWAPLEG)
//		shortIndexIsLibor = true;
//
//	if (spreadoption->GetSpreadLeg()->GetSecondLeg()->GetName() == ARM_SWAPLEG)
//		longIndexIsLibor = true;
//		
//	if (!shortIndexIsLibor)
//// FIXMEFRED: mig.vc8 (22/05/2007 18:07:17): floor(int) doesnt mean anything
//		tenor_1 = ((ARM_CMSLeg*) (((ARM_SpreadLeg*)spreadoption->GetSwapLeg())->GetFirstLeg()))->GetSwapYearTerm();
//
//	if (!longIndexIsLibor)
//// FIXMEFRED: mig.vc8 (22/05/2007 18:07:17): floor(int) doesnt mean anything
//		tenor_2 = ((ARM_CMSLeg*) (((ARM_SpreadLeg*)spreadoption->GetSwapLeg())->GetSecondLeg()))->GetSwapYearTerm();
//	
//	
//	ARM_VectorVector swapLongFixPayTimes;
//	ARM_VectorVector swapLongFixPayPeriods;
//	ARM_VectorVector swapShortFixPayTimes;
//	ARM_VectorVector swapShortFixPayPeriods;	
//
//	if (!longIndexIsLibor)
//	{
//		swapLongFixPayTimes		= ARM_VectorVector(nbCapLet,NULL);
//		swapLongFixPayPeriods	= ARM_VectorVector(nbCapLet,NULL);
//	}
//	if (!shortIndexIsLibor)
//	{
//		swapShortFixPayTimes	= ARM_VectorVector(nbCapLet,NULL);
//		swapShortFixPayPeriods	= ARM_VectorVector(nbCapLet,NULL);
//	}
//
//	int daycount   = spreadoption->GetSpreadLeg()->GetFirstLeg()->GetIRIndex()->GetDayCount();
//	
//	int k;
//    for (k = 0; k < nbCapLet; k++) 
//	{
//		double startDatek = flowStartDates->Elt(k);
//		double endDatek = flowEndDates->Elt(k);
//				
//		(*payPeriods)[k] = CountYears(dayCount,startDatek,endDatek); 
//		
//		/// INDEX SHORT
//		/// CMS index
//		if (!shortIndexIsLibor)
//		{
//			ARM_Date startDate_1 = swapShortFloatStartTime->Elt(k);
//			ARM_Date endDate_1 = startDate_1;
//			endDate_1.AddYears(tenor_1);
//			(*swapShortFloatEndTime)[k]= endDate_1.GetJulian()-asOfDate;
//			(*swapShortFloatStartTime)[k]= startDate_1.GetJulian()-asOfDate;
//		}
//		/// libor index
//		else
//		{
//			ARM_Date start = spreadoption->GetSpreadLeg()->GetFirstLeg()->GetFwdRateStartDates()->Elt(k);
//			ARM_Date end   = spreadoption->GetSpreadLeg()->GetFirstLeg()->GetFwdRateEndDates()->Elt(k);
//			int daycount   = spreadoption->GetSpreadLeg()->GetFirstLeg()->GetIRIndex()->GetDayCount();
//			double period  = CountYears(daycount, start, end);
//
//			(*swapShortFloatStartTime)[k]   = start.GetJulian() - asOfDate;
//			(*swapShortFloatEndTime)[k]     = end.GetJulian()   - asOfDate;
//			
//			swapShortFixPayTimes.push_back(new std::vector<double>(1,(*swapShortFloatEndTime)[k]));
//			swapShortFixPayPeriods.push_back(new std::vector<double>(1,period));
//		}
//		
//		/// INDEX LONG
//		/// CMS index
//		if (!longIndexIsLibor)
//		{
//			ARM_Date startDate_2 = swapLongFloatStartTime->Elt(k);
//			ARM_Date endDate_2 = startDate_2;
//			endDate_2.AddYears(tenor_2);
//			(*swapLongFloatEndTime)[k]= endDate_2.GetJulian()-asOfDate;
//			(*swapLongFloatStartTime)[k]= startDate_2.GetJulian()-asOfDate;
//		}
//		/// libor  index
//		else
//		{
//			ARM_Date start = spreadoption->GetSpreadLeg()->GetSecondLeg()->GetFwdRateStartDates()->Elt(k);
//			ARM_Date end   = spreadoption->GetSpreadLeg()->GetSecondLeg()->GetFwdRateEndDates()->Elt(k);
//			int daycount   = spreadoption->GetSpreadLeg()->GetSecondLeg()->GetIRIndex()->GetDayCount();
//			double period  = CountYears(daycount, start, end);
//
//			(*swapLongFloatStartTime)[k]   = start.GetJulian() - asOfDate;
//			(*swapLongFloatEndTime)[k]     = end.GetJulian()   - asOfDate;
//			
//			swapLongFixPayTimes.push_back(new std::vector<double>(1,(*swapLongFloatEndTime)[k]));
//			swapLongFixPayPeriods.push_back(new std::vector<double>(1,period));
//		}
//
//
//		ARM_Date theodate = theopaymentDates->Elt(k);
//
//		ARM_ReferenceValue* not = spreadoption->GetAmount();
//		if (not)
//			(*notional)[k] = not->CptReferenceValue(theodate.GetJulian());
//		else
//			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ConvertVanillaSpreadOption : unresolved error");
//		
//		(*strikes)[k]  = spreadoption->GetStrikes()->CptReferenceValue(theodate.GetJulian())/100;
//		
//		ARM_Date payDate_k = payTimes->Elt(k);
//		(*payTimes)[k]= payDate_k.GetJulian()-asOfDate;
//
//		ARM_Date resetDate_k = resetTimes->Elt(k);
//		(*resetTimes)[k]= resetDate_k.GetJulian()-asOfDate;
//	}
//
//
//	/// -------------------------------------------------------
//	/// this part of code is specific to 
//	///	the CMS SPREAD CORRIDOR (fixed or variable payments)
//	/// -------------------------------------------------------
//	///	therefore we have to extract:
//	///		o type of payments
//	///		o payment rate characteristics if variable payments
//	///		o fixValues if fixed payments
//	///		o spread1 & spread2
//	///	we also have to overwrite
//	///		o payPeriods (to have sum weights = 1)
//	///		o notional 
//	/// -------------------------------------------------------
//	bool isDigital = false;
//	int payIndexType = K_FIXED;
//	std::vector<double>& swapPayFloatStartTime	= NULL;
//	std::vector<double>& swapPayFloatEndTime		= NULL;
//	std::vector<double>& fixValues				= NULL;
//	std::vector<double>& payIndexLeverages		= NULL;
//	double spread1 (0.0);
//	double spread2 (0.0);
//	ARM_VectorVector swapPayFixPayTimes;
//	ARM_VectorVector swapPayFixPayPeriods;	
//
//	ARM_IntVector* periodIndexes			= NULL;
//	std::vector<double>& payIndexResetTimes		= NULL;
//
//	if (spreadoption->IsDigital())
//	{
//		isDigital = true;
//		spread1 = spreadoption->GetSpread1()/100.0;
//		spread2 = spreadoption->GetSpread2()/100.0;
//
//		for (int l = 0; l < nbCapLet; l++)
//		{
//			ARM_Date theodate = theopaymentDates->Elt(l);
//			notional[l] = spreadoption->GetPayOffs()->CptReferenceValue(theodate.GetJulian());
//		}
//	}
//	
//	if( spreadoption->IsCorridorSpread() || spreadoption->IsDigitalFLT() )
//	{
//		/// whe now consider the number of payment periods
//		size_t nbPeriods = spreadoption->GetPayIndexLeg()->GetPaymentDates()->size();
//
//		swapPayFloatEndTime    = new std::vector<double>(nbPeriods);
//		swapPayFloatStartTime  = new std::vector<double>(nbPeriods);
//
//		ARM_Vector* periodResetDates = spreadoption->GetPayIndexLeg()->GetResetDates();
//		ARM_Vector* periodStartDates = spreadoption->GetPayIndexLeg()->GetFlowStartDates();
//		ARM_Vector* periodEndDates   = spreadoption->GetPayIndexLeg()->GetFlowEndDates();
//		ARM_Vector* periodPayDates   = spreadoption->GetPayIndexLeg()->GetPaymentDates();
//
//		spread1 = spreadoption->GetSpread1()/100.0;
//		spread2 = spreadoption->GetSpread2()/100.0;
//
//		payIndexType			= spreadoption->GetPayIndexLeg()->GetIRIndex()->GetIndexType();
//		fixValues				= new std::vector<double>(nbCapLet);
//		
//		double tenor_p;
//
//		payIndexLeverages	= new std::vector<double>(nbPeriods);
//		payIndexResetTimes	= new std::vector<double>(nbPeriods);
//		periodIndexes		= new ARM_IntVector(nbCapLet);
//			
//		if (spreadoption->GetPayIndexWeight() == 0.0 || payIndexType == IDXFIXED)
//			payIndexType = K_FIXED;
//		else
//		{
//			/// payIndexType == K_LIBOR
//			if ( spreadoption->GetPayIndexLeg()->GetName() == ARM_SWAPLEG )
//			{
//				payIndexType = K_LIBOR;
//				tenor_p = spreadoption->GetPayIndexLeg()->GetIRIndex()->GetYearTerm();
//			}
//			else
//			{	/// payIndexType == K_CMS
//				payIndexType = K_CMS;
//				if ( spreadoption->GetPayIndexLeg()->GetName() == ARM_CMSLEG )
//					tenor_p = ((ARM_CMSLeg *) spreadoption->GetPayIndexLeg())->GetSwapYearTerm();
//				else if (spreadoption->GetPayIndexLeg()->GetName() == ARM_CMTLEG ) 
//					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
//											"ConvertVanillaSpreadOption : CMT payment not supported");
//			}					
//		}
//
//				
//
//		for (k = 0; k < nbCapLet; k++) 
//		{
//			double startDatek		= flowStartDates->Elt(k);
//			double endDatek			= flowEndDates->Elt(k);
//			double payDatek			= flowPayDates->Elt(k);
//			int periodIndex			= periodPayDates->find(payDatek,ARM_VanillaSpreadOptionArg::PayDateToleranceInDays);
//			periodIndexes->Elt(k)   = periodIndex;
//			
//			if (periodIndex == -1)
//				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ConvertVanillaSpreadOption : unresolved error");
//				
//			double periodResetDate	= periodResetDates->Elt(periodIndex);
//			double periodStartDate	= periodStartDates->Elt(periodIndex);
//			double periodEndDate	= periodEndDates->Elt(periodIndex);
//			
//			(*payPeriods)[k]		= CountYears(KACTUAL_360, startDatek, endDatek)/CountYears(KACTUAL_360, periodStartDate, periodEndDate)*CountYears(daycount, periodStartDate, periodEndDate);
//			(*fixValues)[k]			= spreadoption->GetPayIndexMargins()->CptReferenceValue(periodResetDate);
//			
//			// Take care the kernel is vicious
//			ARM_Date theodate = theopaymentDates->Elt(k);
//			ARM_ReferenceValue* not = spreadoption->GetPayIndexLeg()->GetAmount();
//			if (not)
//				(*notional)[k] = not->CptReferenceValue(theodate.GetJulian());
//			else
//				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "ConvertVanillaSpreadOption : unresolved error");
//
//			(*strikes)[k]  = spreadoption->GetStrikes()->CptReferenceValue(periodResetDate)/100;
//		}
//
//		//
//		// payment rate caracteristics
//		// 
//		for (k = 0; k<nbPeriods; k++)
//		{			
//			(*payIndexLeverages)[k]	 = spreadoption->GetPayIndexWeight();
//			(*payIndexResetTimes)[k] = periodResetDates->Elt(k) - asOfDate;
//
//						   			
//			/// CMS index
//			if (payIndexType == K_CMS)
//			{
//				ARM_Date startDate_p = spreadoption->GetPayIndexLeg()->GetFwdRateStartDates()->Elt(k);
//				ARM_Date endDate_p   = startDate_p;
//				endDate_p.AddYears(tenor_p);
//				(*swapPayFloatEndTime)[k]= endDate_p.GetJulian()-asOfDate;
//				(*swapPayFloatStartTime)[k]= startDate_p.GetJulian()-asOfDate;
//			}
//			/// libor  index
//			else if (payIndexType == K_LIBOR)
//			{
//				ARM_Date start = spreadoption->GetPayIndexLeg()->GetFwdRateStartDates()->Elt(k);
//				ARM_Date end   = spreadoption->GetPayIndexLeg()->GetFwdRateEndDates()->Elt(k);;
//				int daycount   = spreadoption->GetPayIndexLeg()->GetIRIndex()->GetDayCount();
//				double period  = CountYears(daycount, start, end);
//
//				(*swapPayFloatStartTime)[k]   = start.GetJulian() - asOfDate;
//				(*swapPayFloatEndTime)[k]     = end.GetJulian()   - asOfDate;
//				
//				swapPayFixPayTimes.push_back(new std::vector<double>(1,(*swapPayFloatEndTime)[k]));
//				swapPayFixPayPeriods.push_back(new std::vector<double>(1,period));
//			}
//			else
//			{
//				(*swapPayFloatStartTime)[k]   = 0.0;
//				(*swapPayFloatEndTime)[k]     = 0.0;
//				
//				swapPayFixPayTimes.push_back(new std::vector<double>(1, 0.0));
//				swapPayFixPayPeriods.push_back(new std::vector<double>(1, 0.0));
//			}
//
//		}
//
//	}	
//	
//	
//	return new ARM_VanillaSpreadOptionArg(
//						curveName,
//						evalTime,
//						CallPut,
//						expiryTime,
//						startTime,
//						endTime,
//						resetTimes,
//						payTimes,
//						payPeriods,
//						notional,
//						coeffLong,
//						coeffShort,
//						strikes,
//						swapLongFloatStartTime,
//						swapLongFloatEndTime,
//						swapLongFixPayTimes,
//						swapLongFixPayPeriods,
//						swapShortFloatStartTime,
//						swapShortFloatEndTime,
//						swapShortFixPayTimes,
//						swapShortFixPayPeriods,
//						longIndexIsLibor,
//						shortIndexIsLibor,
//						isDigital,
//						fixValues,
//						spread1,
//						spread2,
//						payIndexType,
//						payIndexLeverages,
//						swapPayFloatStartTime,
//						swapPayFloatEndTime,
//						swapPayFixPayTimes,
//						swapPayFixPayPeriods,
//						payIndexResetTimes,
//						periodIndexes);
//	
//}

////////////////////////////////////////////////////
///	Routine: ConvertVanillaIRSwapleg
///	Returns: 
///	Action : transform a SwapLeg into an ARM_VanillaIRSwaplegArg
////////////////////////////////////////////////////
ARM_VanillaIRSwaplegArg* ARM_ConverterFromKernel::ConvertVanillaIRSwapleg( ARM_SwapLeg* swapleg, double asOfDate, const string& modelName ){
	return NULL;
}/*
{
	/// fills the specific ARM_VanillaIRSwaplegArg arguments
    size_t size		= swapleg->GetFlowEndDates()->GetSize();
	double Spread	= swapleg->GetSpread();
	int DecompFreq	= swapleg->GetDecompFreq();
	std::vector<double>& ResetTimes    = To_pARM_GP_Vector( swapleg->GetResetDates() );
	std::vector<double>& StartTimes    = To_pARM_GP_Vector( swapleg->GetFwdRateStartDates() );
	std::vector<double>& EndTimes      = To_pARM_GP_Vector( swapleg->GetFwdRateEndDates() );
	std::vector<double>& PayTimes      = To_pARM_GP_Vector( swapleg->GetPaymentDates() );
	std::vector<double>& PayPeriods    = To_pARM_GP_Vector( swapleg->GetInterestTerms() );
	std::vector<double>& Nominals      = new std::vector<double>(size);
	std::vector<double>& Strikes       = new std::vector<double>(size);

	for(size_t i=0;i<size;++i)
	{
		(*Strikes)[i] = 0.0;
		(*ResetTimes)[i] -= asOfDate;
	
		(*StartTimes)[i] -= asOfDate;
		(*EndTimes)[i] -= asOfDate;
		
		(*Nominals)[i] = swapleg->GetAmount()->CptReferenceValue( (*PayTimes)[i] );
		(*PayTimes)[i] -= asOfDate;
	}

	/// part on the curve
    string CurveName = modelName;
    double evalTime = 0.0;
	double mktPrice = -1.0;

	return new ARM_VanillaIRSwaplegArg(
		CurveName,
		evalTime,
		K_CALL,
		Nominals,
		ResetTimes,
		StartTimes,
		EndTimes,
		Strikes,
		PayTimes,
		PayPeriods,
		Spread,
		DecompFreq	);
}*/


////////////////////////////////////////////////////
///	Routine: ConvertVanillaInfSwapleg
///	Returns: 
///	Action : transform an infLeg into an ARM_VanillaInfSwaplegArg
////////////////////////////////////////////////////
ARM_VanillaInfSwaplegArg* ARM_ConverterFromKernel::ConvertVanillaInfSwapLeg	( ARM_InfLeg* infleg, double asOfDate, const string& modelName ){
	return NULL;
}/*
{
	std::vector<double>& numResetDates = infleg->GetNumResetDates();
	std::vector<double>& numFlowStartDates = infleg->GetNumPublicationDates();
	std::vector<double>& numFlowPaymentDates = To_pARM_GP_Vector(infleg->GetPaymentDates() );
	std::vector<double>& numInterestTerms = To_pARM_GP_Vector( infleg->GetInterestTerms() );
	std::vector<double>& numInterestDays = To_pARM_GP_Vector( infleg->GetInterestDays() );
	std::vector<double>& denomResetDates = infleg->GetDenomResetDates();

	ARM_DateStripPtr numDateStrip( new ARM_DateStrip( numFlowStartDates, numFlowPaymentDates, numFlowStartDates, numFlowPaymentDates, numResetDates, numFlowPaymentDates, numInterestDays, numInterestTerms ) );
	ARM_DateStripPtr denomDateStrip( new ARM_DateStrip( numFlowStartDates, numFlowPaymentDates, numFlowStartDates, numFlowPaymentDates, denomResetDates, numFlowPaymentDates, numInterestDays, numInterestTerms ) );
	string irModelName = infleg->getPayCurrencyName();
	string infModelName = infleg->getIndexName();
	int infLegType = infleg->GetSwapType();
	double spread = ( infLegType == K_YEARTOYEAR_LEG )? (1+infleg->GetConstant()/infleg->GetMultiple()) : ( infleg->GetConstant()/infleg->GetMultiple() );

	if ( infLegType == K_OATTYPE_LEG && infleg->GetNxFlag() != K_NX_NONE )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NX_NONE is the only final notional type supported");

	return new ARM_VanillaInfSwaplegArg( irModelName, infModelName, asOfDate, infleg->GetSwapType(), K_CALL, numDateStrip, denomDateStrip,0,0, spread );
}*/

////////////////////////////////////////////////////
///	Routine: ConvertVanillaInfCapFloor
///	Returns: 
///	Action : transform an ARM_InfCapFloor into an ConvertVanillaInfCapFloor
////////////////////////////////////////////////////
ARM_VanillaInfCapArg* ARM_ConverterFromKernel::ConvertVanillaInfCapFloor( ARM_InfCapFloor* infCapFloor, double asOfDate, const string& modelName )
{
	return NULL;
}
/*
{
	ARM_SwapLeg* underlyingSwapLeg = infCapFloor->GetSwapLeg();
	ARM_InfLeg* infleg = dynamic_cast<ARM_InfLeg*> (underlyingSwapLeg);

	if (!infleg )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Trying to price something that is not an InfCap as an InfCap. Please Advise.");

	std::vector<double>& numResetDates = infleg->GetNumResetDates();
	std::vector<double>& numFlowStartDates = infleg->GetNumPublicationDates();
	std::vector<double>& numFlowPaymentDates = To_pARM_GP_Vector(infleg->GetPaymentDates() );
	std::vector<double>& numInterestTerms = To_pARM_GP_Vector( infleg->GetInterestTerms() );
	std::vector<double>& numInterestDays = To_pARM_GP_Vector( infleg->GetInterestDays() );
	std::vector<double>& denomResetDates = infleg->GetDenomResetDates();

	ARM_DateStripPtr numDateStrip( new ARM_DateStrip( numFlowStartDates, numFlowPaymentDates, numFlowStartDates, numFlowPaymentDates, numResetDates, numFlowPaymentDates, numInterestDays, numInterestTerms ) );
	ARM_DateStripPtr denomDateStrip( new ARM_DateStrip( numFlowStartDates, numFlowPaymentDates, numFlowStartDates, numFlowPaymentDates, denomResetDates, numFlowPaymentDates, numInterestDays, numInterestTerms ) );
	string irModelName = infleg->getPayCurrencyName();
	string infModelName = infleg->getIndexName();
	int infLegType = infleg->GetSwapType();
	double spread = ( infLegType == K_YEARTOYEAR_LEG )? (1+infleg->GetConstant()/infleg->GetMultiple()) : ( infleg->GetConstant()/infleg->GetMultiple() );
	double strike = infCapFloor->GetStrike()*0.01/infleg->GetMultiple();

	if ( infLegType == K_OATTYPE_LEG && infleg->GetNxFlag() != K_NX_NONE )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NX_NONE is the only final notional type supported");


	return new ARM_VanillaInfCapArg( irModelName, infModelName, asOfDate, infleg->GetSwapType(), infCapFloor->GetCapFloorType(),strike, numDateStrip, denomDateStrip, 0, 0, spread );
}*/


////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToFXOption
///	Returns: 
///	Action : transform an Option into an ARM_VanillaFXOption
////////////////////////////////////////////////////
ARM_VanillaPortfolio* ARM_ConverterFromKernel::ConvertPortoflio( ARM_Portfolio* port, double asOfDate, const string& modelName )
{
	return NULL;
}/*
{
	size_t portSize=port->size();
	if( portSize )
	{
		std::vector<double>& weights		= To_pARM_GP_Vector( port->GetWeights()	);
		ARM_AutoCleaner<std::vector<double>> HoldWeights(weights);
		std::vector<double>& mktprices	= To_pARM_GP_Vector( port->GetMktPrices() );
		ARM_AutoCleaner<std::vector<double>> HoldMktPrices(mktprices);
		std::vector<double>& vegas		= To_pARM_GP_Vector( port->GetPrecision() );
		ARM_AutoCleaner<std::vector<double>> HoldVegas(vegas);
		vector<ARM_VanillaArgPtr> assets(portSize);

		for( size_t i=0; i<portSize; ++i )
			assets[i] = ARM_VanillaArgPtr( ARM_ConverterFromKernel::ConvertSecuritytoArgObject( port->GetAsset(i), asOfDate, modelName) );
		return new ARM_VanillaPortfolio( assets, *weights, *mktprices, *vegas );
	}
	else
	{
		std::vector<double> weights(0);
		std::vector<double> mktprices(0);
		std::vector<double> vegas(0);
		vector<ARM_VanillaArgPtr> assets(0);
		return new ARM_VanillaPortfolio( assets, weights, mktprices, vegas );
	}
}*/

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertVanillaSumOpt
///	Returns: 
///	Action : transform an Option into an ARM_VanillaFXOption
////////////////////////////////////////////////////
ARM_VanillaSumOptArg* ARM_ConverterFromKernel::ConvertVanillaSumOpt( ARM_SumOpt* sumOpt, double asOfDate, const string& modelName)
{
	return NULL;
}/*
{

	string curveName = modelName;
    double evalTime = 0.0;

	ARM_GP_VectorPtr FwdResetDates(static_cast<ARM_GP_Vector*>(sumOpt->GetFwdResetDates()->Clone()));
	ARM_GP_VectorPtr FwdStartDates(static_cast<ARM_GP_Vector*>(sumOpt->GetFwdStartDates()->Clone()));
	ARM_GP_VectorPtr FwdEndDates(static_cast<ARM_GP_Vector*>(sumOpt->GetFwdEndDates()->Clone()));

	size_t i;

	for (i = 0; i < FwdResetDates->size(); ++i)
	{
		(*FwdResetDates)[i] -= asOfDate;
		(*FwdStartDates)[i] -= asOfDate;
		(*FwdEndDates)[i] -= asOfDate;
	}

	double payDate = sumOpt->GetPayDate().GetJulian();
	payDate -= asOfDate;

	return new ARM_VanillaSumOptArg( 
		curveName,
		evalTime,
		sumOpt->GetCapFloor(),
		sumOpt->GetCoeffs(),
		FwdResetDates,
		FwdStartDates,
		FwdEndDates,
		payDate,
		sumOpt->GetFwdPeriods(),
		sumOpt->GetStrike(),
		sumOpt->GetVolatilityRatio());
}*/

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertToSwaptionSmile
///	Returns: 
///	Action : transform an Option into an ARM_VanillaFXOption
////////////////////////////////////////////////////
ARM_VanillaSmiledSwaption* ARM_ConverterFromKernel::ConvertVanillaSmiledSwaption	( ARM_SmiledSwaption* opt,double asOfDate, const string& modelName)
{
	return NULL;
}
/*
{
	ARM_VanillaSwaptionArg* baseSwaption = ConvertVanillaSwaption(opt,asOfDate,modelName);
	std::vector<double> data = opt->GetDatas();
	ARM_VanillaSmiledSwaption* result = new ARM_VanillaSmiledSwaption(*baseSwaption,data);
	delete baseSwaption;
	return result;
}*/

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: TranferToSumOptSec
///	Returns: 
///	Action : Transfer information from the vanilla arg 
/// to the security
////////////////////////////////////////////////////
void ARM_ConverterFromKernel::TranferToSumOptSec(ARM_VanillaSumOptArg* sumOptArg, ARM_SumOpt* sumOptSec)
{
}	
/*
{
	if ((sumOptArg != NULL) && (sumOptSec != NULL))
	{
		sumOptSec->SetSumFwd(sumOptArg->GetSumFwd());
		sumOptSec->SetSumVol(sumOptArg->GetSumVol());
	}
}*/

////////////////////////////////////////////////////
///	Class  : ARM_ConverterFromKernel
///	Routine: ConvertIrFxSwaption
///	Returns: 
///	Action : transform a hybrid portfolio into an ARM_VanillaIrFxSwaption
////////////////////////////////////////////////////
ARM_VanillaIrFxSwaption* ARM_ConverterFromKernel::ConvertIrFxSwaption(ARM_OptionPortfolio* irFxSwaption,double asOfDate, const string& modelName)
{
	return NULL;
}
//{
//	ARM_Portfolio* irFxBasket = irFxSwaption->GetPtf();
//	ARM_VanillaFxOption* fxOption;
//
//	int nbFwdFxs = irFxBasket->size()-1;
//	std::vector<double> fxResetTimes(nbFwdFxs);
//	std::vector<double> fxSettlementTimes(nbFwdFxs);
//	std::vector<double> fxPayTimes(nbFwdFxs);
//	std::vector<double> fxNominals(nbFwdFxs);
//
//	/// Convert Fx options
//	for(size_t i=0;i<nbFwdFxs;++i)
//	{
//		fxOption = ARM_ConverterFromKernel::ConvertToFXOption(static_cast<ARM_Option*>(irFxBasket->GetAsset(i)),asOfDate);
//		fxResetTimes[i]		= fxOption->GetExpiry();
//		fxSettlementTimes[i]= fxOption->GetFwdTime();
//		fxPayTimes[i]		= fxOption->GetPayTime();
//		fxNominals[i]		= irFxBasket->GetWeight(i);
//	}
//
//
//	/// Convert variable notional swaption
//	ARM_VanillaSwaptionArg* irSwaption = ARM_ConverterFromKernel::ConvertVanillaSwaption(static_cast<ARM_Swaption*>(irFxBasket->GetAsset(nbFwdFxs)),asOfDate,modelName);
//
//	double irFxSwaptionExpiryTime = irFxSwaption->GetExpiryDate().GetJulian()-asOfDate;
//	int callPut = irFxSwaption->IsCall() ? K_CALL : K_PUT;
//	double strike = 0.01 * irFxSwaption->GetStrike()->CptReferenceValue(asOfDate); // Kernel uses % !
//	return new ARM_VanillaIrFxSwaption(modelName,irFxSwaptionExpiryTime,callPut,strike,fxResetTimes,fxSettlementTimes,fxPayTimes,fxNominals,ARM_CountedPtr<ARM_VanillaSwaptionArg>(irSwaption));
//}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

