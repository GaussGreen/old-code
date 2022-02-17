/*!
 *
 * Copyright (c) IXIS-CIB March 2006
 *
 *	\file mktdata.cpp
 *  \brief file for the mkt datas container
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date March 2006
 */

#include "xxxproject/vanillapricer.h"
#include "xxxproject/mktdatas.h"
#include "xxxproject/observator.h"
#include "xxxproject/argconvdefault.h"

#include <pricer\ipricer.h>

#include <inst/capfloor.h>
#include <inst/armdigital.h>
#include <inst/swaption.h>
#include <inst/corridorleg.h>
#include <inst/option.h>
#include <inst/forex.h>
#include <inst/spreadoption.h>


#include <gpcalculators/pricerfactory.h>
#include <gpcalculators/gencalculator.h>
#include <gpinfra/gensecurity.h>
#include <gpinfra/pricingmodel.h>
#include <gpinfra/pricingadviser.h>
#include <gpinfra/dealdescription.h>

#include <gpinfra/mktdatamanagerrep.h>
#include <gpinfra\mktdatamanager.h>
#include <gpcalib\calibmethod.h>
#include <mod/bssmiled.h>

#include <gpinflation/infcorridorLeg.h>
#include <gpinflation/infbsmodel.h>

#include <crv\zerointspreaded.h>

CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
///	Class  : ARM_VanillaPricer
///	Routine: PriceWithModel
///	Returns: 
///	Action : Price a vanilla security with an ARM model
/////////////////////////////////////////////////////////////////

double ARM_VanillaPricer::PriceWithModel(ARM_Model* model )
{
	double price = 0.0;

    if (dynamic_cast<ARM_PricingModel*>(model))
	{
		CC_NS(std,pair)<bool,double> pricingResult = ARM_PricerFactory::Price(itsSecurity, model);
		if (pricingResult.first)
        {
		   price = pricingResult.second;
        }
	}
	else
	{
		ARM_IFPricer pricer(itsSecurity, model);
		price = pricer.Price();
	}

	return price;
}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_CAP_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

ARM_CAP_Pricer::ARM_CAP_Pricer	( ARM_CapFloor* capFloor):	ARM_VanillaPricer	( capFloor ){ }

void ARM_CAP_Pricer::SetMkt( ARM_MktData * mkt ){

	string  key		= GetSecurity()->GetCurrencyUnit()->GetCcyName();
	V_Sub	sub		=  mkt->itsSub;	

	for( int i = 0; i< sub.size();i++){
		if( dynamic_cast<ARM_CAP_Subject<>* > ( sub[i] ) ){
			ARM_Object* tmp = &* sub[i]->GetModel(key)->GetObject( );
			itsModel = dynamic_cast<ARM_Model*>(tmp );
			break;
		}
	}
}

double ARM_CAP_Pricer::Price(  ){	return PriceWithModel(&*itsModel);	}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_OSW_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


ARM_OSW_Pricer::ARM_OSW_Pricer	( ARM_Swaption* OswFloor):	ARM_VanillaPricer	( OswFloor ){ }

void ARM_OSW_Pricer::SetMkt( ARM_MktData * mkt ){

	string  key		= GetSecurity()->GetCurrencyUnit()->GetCcyName();
	V_Sub	sub		=  mkt->itsSub;

	for( int i = 0; i< sub.size();i++){
		if( dynamic_cast<ARM_OSW_Subject<>* > ( sub[i] ) ){
			ARM_Object* tmp = &* sub[i]->GetModel(key)->GetObject( );
			itsModel = dynamic_cast<ARM_Model*>(tmp );
			break;
		}
	}
}

double ARM_OSW_Pricer::Price		(  ){	return PriceWithModel(&*itsModel);	}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_INF_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


ARM_INF_Pricer::ARM_INF_Pricer	( ARM_InfCorridorLeg* infOpt):	ARM_VanillaPricer	( infOpt ){ }


void ARM_INF_Pricer::SetMkt( ARM_MktData * mkt ){
	string			key		=  dynamic_cast<ARM_InfCorridorLeg*> (GetSecurity())->GetInfLeg()->getIndexName();
	V_Sub			sub		=  mkt->itsSub;

	for( int i = 0; i< sub.size();i++){
		if( dynamic_cast<ARM_IBS_Subject<>* > ( sub[i] ) ){
			ARM_Object* tmp = &* sub[i]->GetModel(key)->GetObject( );
			itsModel = dynamic_cast<ARM_Model*>(tmp );
			break;
		}
	}
}

double ARM_INF_Pricer::Price		(  ){	return PriceWithModel(&*itsModel);	}

					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_FX_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
/*
ARM_FX_Pricer::ARM_FX_Pricer ( ARM_Option* fxOption ):	ARM_VanillaPricer	( fxOption ){ }

void ARM_FX_Pricer::	SetMkt ( ARM_MktData * mkt ){
	
	int test =0;
	ARM_MktData * tmp = NULL;
	
	ARM_Forex* forex= dynamic_cast<ARM_Forex*>( ((ARM_Option*) itsSecurity)->GetUnderlying());
	string	forCcy	= forex->GetMoneyCurrency()->GetCcyName();
	string	domCcy	= forex->GetMainCurrency()->GetCcyName();

	tmp = new ARM_MktSubject< ARM_VFX_Subject<> >( *mkt );
	if ( &*tmp->GetObs( forCcy+"/"+domCcy ) ){
		itsObservator	= tmp->GetObs( forCcy+"/"+domCcy );
		test += 1;
	}

	tmp = new ARM_MktSubject< ARM_MIX_Subject<> >( *mkt );
	if ( &*tmp->GetObs( forCcy+"/"+domCcy ) ){
		itsObservator	= tmp->GetObs( forCcy+"/"+domCcy );
		test += 1;
	}
	if ( test>0) 
		itsModel		= ARM_ModelPtr (dynamic_cast<ARM_Model*>(itsObservator->GetObject()->Clone() ) );
	else
		throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "The corresponding Mkt is missing." );

}

double ARM_FX_Pricer::Price	( )	{ return PriceWithModel(&*itsModel);	}
*/

					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*					Class  : ARM_CalculatorPricer				*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
/*
ARM_CalculatorPricer::ARM_CalculatorPricer ( ARM_GenCalculator* genCalculator ): itsGenCalculator( genCalculator ){ }


void ARM_CalculatorPricer::	SetMkt ( ARM_MktData * mkt ){

	M_Obs tmpObs = mkt->itsObs;
	It_Mbs it;
	ARM_MarketData_ManagerRepPtr MktManager =	itsGenCalculator->GetMktDataManager();
	ARM_MarketData_ManagerImp*   MktMap		=	MktManager->GetMktDataManager();
	ARM_MarketData_ManagerImp::iterator itm;

	for ( itm = MktMap->begin(); itm != MktMap->end(); itm++){
		it = tmpObs.find(itm->first);
		if ( it  != tmpObs.end() ){
			itsObservator.insert(pair<string,ARM_ObservPtr> (it->first,ARM_ObservPtr ( it->second) ) );
			MktManager->RegisterData(it->first,it->second->ReduceObject() );
		}
	}
	itsGenCalculator->SetMktDataManager( ARM_MarketData_ManagerRepPtr( MktManager ) );
	itsGenCalculator->Update();
	
}

double ARM_CalculatorPricer::Price	( )	{	

	ARM_MarketData_ManagerRepPtr MktManager =	itsGenCalculator->GetMktDataManager();
	ARM_MarketData_ManagerImp*   MktMap		=	MktManager->GetMktDataManager();

	for ( It_Mbs it = itsObservator.begin(); it !=itsObservator.end(); it++){
		MktManager->RegisterData(it->first,it->second->ReduceObject() );
	}

	itsGenCalculator->SetMktDataManager( ARM_MarketData_ManagerRepPtr( MktManager ) );
	itsGenCalculator->Update();	
	return itsGenCalculator->Price();	
}

*/
CC_END_NAMESPACE( )

