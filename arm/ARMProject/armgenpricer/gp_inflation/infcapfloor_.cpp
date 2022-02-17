/*----------------------------------------------------------------------------
 *
 *
/*! \file inflcapfloor.h
 *
 *  \brief Vanilla inflation cap and floor
 *	\Copyright (c) NATIXIS July 2007 Paris
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date July 2007
 */

/*----------------------------------------------------------------------------*/

#include "gpinflation/infcapfloor_.h"

/// gpinflation
#include "gpinflation/infdata.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infbssmiledmodel.h"

/// gpbase
#include "gpbase/gpvector.h"
/// To conversion
#include "gpbase/gplinalgconvert.h"		/// for conversion

CC_BEGIN_NAMESPACE( ARM )



InfCapFloor::InfCapFloor(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow, 
							int capFloor, double strike,
							//for backward compatibility : to remove
							int interpType,	double firstReset )
: InfLeg(vCahFlow,interpType,	firstReset), itsCapFloor(capFloor), itsStrike(strike)
{
	int nbCapLet = vCahFlow.size() ;
	itsVolInfo   = new ARM_Matrix( nbCapLet, 5, 0. );
	itsCFValues	 = new ARM_GP_Vector( nbCapLet, 0. );
}


InfCapFloor::InfCapFloor(const InfCapFloor& rhs)
:	InfLeg(rhs.itsLeg, rhs.itsInterpType, rhs.itsFirstReset), 
	itsCapFloor(rhs.itsCapFloor), itsStrike(rhs.itsStrike),
	itsVolInfo(rhs.itsVolInfo), itsCFValues(rhs.itsCFValues)	
{
}
InfCapFloor::~InfCapFloor()
{
	delete itsVolInfo ;
	delete itsCFValues ;
}


InfCapFloor& InfCapFloor::operator = (const InfCapFloor &rhs )
{
	if( this !=	 &rhs )
	{
		InfLeg::operator = ( rhs );
		itsVolInfo = rhs.itsVolInfo ;
		itsCFValues = rhs.itsCFValues ;	
		itsCapFloor = rhs.itsCapFloor ;
		itsStrike = rhs.itsStrike ;	
	}

	return *this;
}
/*
void InfCapFloor::PropagateModel(ARM_Model *model) 
{
	Instrument::SetModel(ARM_CountedPtr<ARM_Model>(model)) ;
}

void  InfCapFloor::SetModel(ARM_CountedPtr<ARM_Model> model) 
{
	itsInstrumentModel = model ;

}
*/
void InfCapFloor::PropagateModel(ARM_Model* model )
{
	//SetModel( model );
}

void InfCapFloor::performCalculation() 
{
	CptExpectedFwdRates() ;
    CptCashFlowValues();
}

double InfCapFloor::ComputePrice(int mode)
{
	//calculate() ;

	CptExpectedFwdRates() ;
    CptCashFlowValues();
	return GetPrice() ;
}
void InfCapFloor::CptCashFlowValues()
{
	//ARM_CountedPtr<ARM_Model> model(GetModel()) ;

	ARM_InfBSSmiledModel* pricingMod = dynamic_cast<ARM_InfBSSmiledModel*> (GetModel()/*&*itsInstrumentModel */);
	if( ! pricingMod )
		ARMTHROW(ERR_INVALID_ARGUMENT," wrong model type in InfCapFloor::ComputePrice");	
	ARM_InfIdx* infidx = &* GetIndex()->GetIndex() ;
	//ARM_CountedPtr<	InfCapFloor> capfloor (this) ;

	double price = 0., optLet, yoy=0., timeToStart = 0., tenor = 1., discountfactor = 1., payTime = 0.;
	double timeToExpiry, remainingTenor, strike=itsStrike, interestTerm ;
	double constant	= 1.; //TODO : check this out bro !

	int nbCapLet	= itsLeg.size(), k;
	int callPut		= itsCapFloor;
	int infRoP		= GetRcvOrPay();
		

	ARM_Date numDate, denomDate;
	ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> > vCoupon = itsLeg ;
	ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >::const_iterator coupon = vCoupon.begin() ;

	for (k = 0; k<nbCapLet; k++) 
	{
		numDate			= ARM_Date( GetNumJulianDates()->Elt(k) );
		denomDate		= ARM_Date( GetDemJulianDates()->Elt(k) );

		//discountfactor	= GetDiscountFactors()[k];
		discountfactor	= (*itsDiscountFactors)[k];
		yoy				= (*itsExpectedFwdRates)[k] +1;
		payTime			= (*itsPaymentTimes)[k];
		interestTerm	= (*itsAccrualTimes)[k];

		timeToStart		= pricingMod->GetModelTimeWPublishLag( denomDate, infidx);
		
		timeToExpiry	= pricingMod->GetModelTimeWPublishLag( numDate, infidx );
		tenor			= timeToExpiry	- timeToStart;
		remainingTenor	= tenor;

	
		//ARM_CountedPtr<StoreVolInfoWithNewCF> StoringFunc( new StoreVolInfoWithNewCF(strike, timeToStart, remainingTenor, k, capfloor ));
		StoreVolInfoWithNewCF StoringFunc( strike, timeToStart, remainingTenor, k, this);

		/// to handle weird case of yoy negative ! //PFFFFFFFFFFFF
		if( yoy < 0. )
		{
			strike		= -strike;
			callPut		= -itsCapFloor;
			yoy			= -yoy;
		}

		//ARM_CountedPtr<InfCapFloorContext> context( new InfCapFloorContext(numDate, denomDate ));
		InfCapFloorContext context( numDate, denomDate );
		optLet				= pricingMod->SingleAssetOptionPrice( yoy, infRoP*strike, (infRoP*callPut), discountfactor, infidx, &context, &StoringFunc );
		itsCFValues->Elt(k)	= optLet/discountfactor*interestTerm;
		price				+= optLet*interestTerm;
	}

	SetPrice(price);
}

void InfCapFloor::StoreVol( double vol, double volLookupStrike, double strike, double timeToStart, double tenor, int k )
{
	itsVolInfo->Elt(k,0 ) = vol;
	itsVolInfo->Elt(k,1 ) = volLookupStrike;
	itsVolInfo->Elt(k,2 ) = strike;
	itsVolInfo->Elt(k,3 ) = timeToStart;
	itsVolInfo->Elt(k,4 ) = tenor;
}

/*
void InfCapFloor::View(char* id, FILE* ficOut)
{
	ARM_SwapLeg* tmpInfLeg = GetSwapLeg();
	ARM_InfLeg* infLeg = dynamic_cast<ARM_InfLeg*>(tmpInfLeg);
	if( !infLeg )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"inflation cap/floor use a leg that is not an inflation leg... Please advise");

	
	FILE* fOut;
	char fOutName[200];
	
	/// first determine that the file is not already opened
	if ( ficOut == NULL )
	{
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
	}
	else
		fOut = ficOut;
	
	/// since the view method is common
	/// for all inflationswap
	/// get the type of the swap
	fprintf(fOut, "\n\n =======> INFLATION CAP AND FLOOR  <====== \n" );
	
	fprintf(fOut, "\n\n Corresponding leg\n" );
	
	GetSwapLeg()->View(id, fOut );

	/// Pricing information
	fprintf(fOut, "\n\n =======> PRICING INFORMATION  <====== \n\n" );

	fprintf(fOut, "Denom Dates\t Num Dates\t Tenor\t Raw Strike\t Lookup Strike\t Forward\t Vol Level\t ForwardOption\n" );

	int nbCapLet = GetNumFlows();
	int i;
	char d1[20];
	char d2[20];
	double numDate, denomDate;

	for (i = 0; i < nbCapLet; ++i )
	{
		numDate			= infLeg->GetNumResetDates()->Elt(i);
		denomDate		= infLeg->GetDenomResetDates()->Elt(i);
		((ARM_Date) denomDate).JulianToStrDate(d1);
		((ARM_Date) numDate).JulianToStrDate(d2);
		fprintf(fOut, "%s\t %s\t %.3lf\t %.3lf\t\t %.6lf\t %.6lf\t %.6lf\t %.6lf\n", d1, d2, 
			itsVolInfo->Elt(i,4),
			itsVolInfo->Elt(i,2), 
			itsVolInfo->Elt(i,1), 
			infLeg->GetCashFlowValues()->Elt(i),
			itsVolInfo->Elt(i,0),
			itsCFValues->Elt(i) 
		);
	}
	
	/// to allow to have nested view
	if ( ficOut == NULL )
		fclose(fOut);
}*/


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

