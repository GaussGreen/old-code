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


#include <mod/bssmiled.h>

CC_BEGIN_NAMESPACE( ARM )

/////////////////////////////////////////////////////////////////
///	Class  : ARM_VanillaPricer
///	Routine: PriceWithModel
///	Returns: 
///	Action : Price a vanilla security with an ARM model
/////////////////////////////////////////////////////////////////

double ARM_VanillaPricer::PriceWithModel(ARM_Model* model ){
	ARM_IFPricer pricer(itsSecurity, model);
	return pricer.Price();
}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_CAP_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

ARM_CAP_Pricer::ARM_CAP_Pricer	( ARM_CapFloor* capFloor):	ARM_VanillaPricer	( capFloor ){ }

void ARM_CAP_Pricer::SetMkt( ARM_MktData * mkt ){

	string forCcy	= GetSecurity()->GetCurrencyUnit()->GetCcyName();
	
	ARM_MktSubject< ARM_CAP_Subject<> > tmp( *mkt );

	itsObservator	= tmp.GetObs( forCcy );
	itsModel		= dynamic_cast<ARM_Model*> ( itsObservator->GetObject() );
}

double ARM_CAP_Pricer::Price( ){	return PriceWithModel(itsModel);	}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_OSW_Pricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/


ARM_OSW_Pricer::ARM_OSW_Pricer	( ARM_Swaption* OswFloor):	ARM_VanillaPricer	( OswFloor ){ }

void ARM_OSW_Pricer::SetMkt( ARM_MktData * mkt ){

	string forCcy	= GetSecurity()->GetCurrencyUnit()->GetCcyName();

	ARM_MktSubject< ARM_OSW_Subject<> > tmp( *mkt );

	itsObservator	= tmp.GetObs( forCcy );
	itsModel		= dynamic_cast<ARM_Model*> (itsObservator->GetObject() );
}

double ARM_OSW_Pricer::Price		(	 ){	return PriceWithModel(itsModel);	}

					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_FxPricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

ARM_FX_Pricer::ARM_FX_Pricer ( ARM_Option* fxOption ):	ARM_VanillaPricer	( fxOption ){ }

void ARM_FX_Pricer::	SetMkt ( ARM_MktData * mkt ){

	ARM_MktData * tmp =NULL;

	ARM_Forex* forex= dynamic_cast<ARM_Forex*>( ((ARM_Option*) itsSecurity)->GetUnderlying());

	string	forCcy	= forex->GetMoneyCurrency()->GetCcyName();
	string	domCcy	= forex->GetMainCurrency()->GetCcyName();

	int n = ARM_ArgConv_MktCurType.GetNumber( domCcy);
	switch( n ){
	
	case 0:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<0.0> >( *mkt );
		break;
	case 1:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<1.0> >( *mkt );
		break;
	case 2:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<2.0> > ( *mkt );
		break;
	case 3:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<3.0> > ( *mkt );
		break;
	case 4:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<4.0> > ( *mkt );
		break;
	}
	itsObservator	= tmp->GetObs( forCcy );
	itsModel		= dynamic_cast<ARM_Model*> ( itsObservator->GetObject() );
}

double ARM_FX_Pricer::Price	(	)	{	return PriceWithModel(itsModel);	}


					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/
					/*																*/
					/*						Class  : ARM_TarnFxPricer					*/
					/*																*/		
					/*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*	*/

/*ARM_FX_Pricer::ARM_FX_Pricer ( ARM_PricingModel* prModel ):	ARM_VanillaPricer	( prModel ){ }

void ARM_FX_Pricer::	SetMkt ( ARM_MktData * mkt ){

	ARM_MktData * tmp =NULL;

	ARM_Forex* forex= dynamic_cast<ARM_Forex*>( ((ARM_Option*) itsSecurity)->GetUnderlying());

	string	forCcy	= forex->GetMoneyCurrency()->GetCcyName();
	string	domCcy	= forex->GetMainCurrency()->GetCcyName();

	int n = ARM_ArgConv_MktCurType.GetNumber( domCcy);
	switch( n ){
	
	case 0:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<0.0> >( *mkt );
		break;
	case 1:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<1.0> >( *mkt );
		break;
	case 2:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<2.0> > ( *mkt );
		break;
	case 3:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<3.0> > ( *mkt );
		break;
	case 4:
		tmp = new ARM_MktSubject< ARM_VFX_Subject<4.0> > ( *mkt );
		break;
	}
	itsObservator	= tmp->GetObs( forCcy );
	itsModel		= dynamic_cast<ARM_Model*> ( itsObservator->GetObject() );
}

double ARM_FX_Pricer::Price	(	)	{	return PriceWithModel(itsModel);	}
*/

CC_END_NAMESPACE( )

