/*!
 *
 * Copyright (c) NATIXIS May 2007 Paris
 *
 *	\file infComposedLeg.cpp
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date May 2007
 */


/// gpinflation
#include "gpinflation/infComposedLeg.h"
#include "gpinflation/infPayOff.h"
#include "gpinflation/infFloatingCoupon.h"
/// gpbase
#include "gpbase/gplinalgconvert.h"
#include "gpbase/checkarg.h"




CC_BEGIN_NAMESPACE( ARM )

InfComposedLeg::InfComposedLeg(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow,
				//for backward compatibility : to remove
				int interpType,	double firstReset): 
Instrument<InfComposedCoupon<YOYPayOff> >(vCahFlow),	
ARM_SwapLeg(), 
itsInterpType(interpType), itsFirstReset(firstReset)
{

	Init();


	int nFlows = itsLeg.size(),i, j=0, dayCount;

	double time = 0 ;
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;

	ARM_CountedPtr<InfPerformance > yoy ;
	ARM_CountedPtr<ARM_InfIdx > cpi ;
	ARM_CountedPtr<ARM_InfCurv> cpiCurve ;

	ARM_Date resDate, payDate, today, tmp, demResetDate ;
	Period tenor ;
	
	coupon	=	*itsLeg.begin() ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
	int nResets = coupon->GetResDates().size() * nFlows;

	ARM_Vector startDates(nResets);
	ARM_Vector intDays(nResets);
	ARM_Vector intTerms(nResets);
	ARM_Vector endDates(nResets);
	ARM_Vector paymentDates(nResets); 
	ARM_Vector resetDates(nResets); 

	itsResetTimes		= new ARM_GP_Vector(nResets) ;
	itsPaymentTimes		= new ARM_GP_Vector(nResets) ;
	itsAccrualTimes		= new ARM_GP_Vector(nResets) ;
	itsNumJulianDates	= new ARM_GP_Vector(nResets) ;
	itsDemJulianDates	= new ARM_GP_Vector(nResets) ;

	ARM_GP_T_Vector<ARM_Date> vResetDates ;
	ARM_GP_T_Vector<ARM_Date>::const_iterator resetDate ;
	
	for( i=0; i<nFlows; ++i )
	{
		coupon	=	itsLeg[i] ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfComposedLeg::InfComposedLeg() NULL coupon	");
		yoy			= coupon->GetIndex() ;
		if (yoy.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfComposedLeg::InfComposedLeg() NULL YOY performance	");
		cpi			= yoy->GetIndex() ;
		if (cpi.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfComposedLeg::InfComposedLeg() NULL CPI index	");
		cpiCurve  = cpi->GetInfFwdCurve(); 
		if (!cpiCurve.IsNull()) 
			today					= cpi->GetInfFwdCurve()->GetAsOf() ;
		/* 
			s'il n'y a pas de courbe infl today = today()  (cf. ARM_Date)
			TODO : exiger que ARM_InfIdx contienne une courbe infl
		*/
		tenor					= yoy->GetTenor() ;
		
		vResetDates		= coupon->GetResDates() ;
		resetDate		= vResetDates.begin() ;
		
		for (; resetDate!=vResetDates.end(); ++resetDate,++j) 
		{
			/**/
			startDates[j]			= coupon->GetStartDate().GetJulian();
			endDates[j]				= (coupon->GetEndDate().GetJulian()) ;
			payDate					= coupon->GetPayDate() ;
			paymentDates[j]			= (payDate.GetJulian()) ;
			dayCount				= cpi->GetDayCount() ; 
			
			time					= CountYears(dayCount, today, payDate);
			(*itsPaymentTimes)[j]	= time;
			(*itsAccrualTimes)[j]	= coupon->GetAccrualTime() ;
			intTerms[j]				= coupon->GetAccrualTimeBetweenResets() ;
			intDays[j]				= coupon->GetAccrualDaysBetweenResets() ;	
			/**/
			resDate					= *resetDate ;
			tmp						= *resetDate ;
			demResetDate			= tmp.AddPeriodNoAdj(tenor.GetUnit(), -tenor.GetLength());
			time					= CountYears(dayCount, today, resDate ) ;
			resetDates[j]			= (*resetDate).GetJulian() ;
			(*itsNumJulianDates)[j]	= (*resetDate).GetJulian() ;
			(*itsDemJulianDates)[j]	= demResetDate.GetJulian() ;
			(*itsResetTimes)[j]		= time ;

		}
	}

	//for backward compatibility : to remove
	SetFlowStartDates(&startDates);
	SetFlowEndDates(&endDates);
	SetPaymentDates( new ARM_Vector(paymentDates) );
	SetResetDates(&resetDates);
	SetInterestTerms(&intTerms) ; 
	SetInterestDays(&intDays) ;
}	

InfComposedLeg::InfComposedLeg(	const InfComposedLeg& InfComposedLeg): 
Instrument<InfComposedCoupon<YOYPayOff> >(InfComposedLeg.itsLeg),	ARM_SwapLeg(), itsNumJulianDates(InfComposedLeg.itsNumJulianDates), 
itsDemJulianDates(InfComposedLeg.itsDemJulianDates), itsPaymentTimes(InfComposedLeg.itsPaymentTimes),
itsAccrualTimes(InfComposedLeg.itsAccrualTimes),itsResetTimes(InfComposedLeg.itsResetTimes)
{
	//for backward compatibility : to remove
	SetFlowStartDates( InfComposedLeg.GetFlowStartDates() ) ;

}

InfComposedLeg&	
InfComposedLeg::operator=(	const InfComposedLeg& InfComposedLeg)
{
	if( this !=	 &InfComposedLeg )
	{
		ARM_SwapLeg::operator = ( InfComposedLeg );
		itsLeg = InfComposedLeg.itsLeg ;
		itsDemJulianDates	= InfComposedLeg.itsDemJulianDates ;
		itsNumJulianDates	= InfComposedLeg.itsNumJulianDates ;
		itsPaymentTimes		= InfComposedLeg.itsPaymentTimes ;
		itsAccrualTimes		= InfComposedLeg.itsAccrualTimes ;
		itsResetTimes		= InfComposedLeg.itsResetTimes ;

		//for backward compatibility : to remove
		SetFlowStartDates( new ARM_Vector(InfComposedLeg.GetFlowStartDates() ) ) ;

	}
	return *this;
}
	
/* 
 TODO : obliger à utiliser setModel() avant les calculate() ... 
 comme dans Instrument::Computeprice(Model)
*/
void InfComposedLeg::performCalculation() 
{
    CptExpectedFwdRates();
    CptCashFlowValues();
}
/*
void  InfComposedLeg::SetModel(ARM_CountedPtr<ARM_Model> model) 
{
	itsInstrumentModel = model ;
}
*/
double InfComposedLeg::ComputePrice(int mode)
{
	calculate() ;
	return GetPrice() ;
}

void InfComposedLeg::CptExpectedFwdRates()
{
	/// FwdCPIRatio returns in the following order
	/// CPIRatio[0] is the ratio of CPI
	/// CPIRatio[1] is the numerator CPI
	/// CPIRatio[2] is the denominator CPI

	InfFwdModel* pricingMod = dynamic_cast<InfFwdModel*> (GetModel());

	if( ! pricingMod ) ARM_THROW( ERR_INVALID_ARGUMENT, " the model you use cannot price inflation derivatives");

	int nFlows = itsLeg.size();
	double firstReset = 0., tmp = 0., spread = 0.;

	ARM_GP_Vector CPIRatio;

	int i;
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;
	ARM_CountedPtr<ARM_InfIdx > cpi ;
	ARM_CountedPtr<InfPerformance > yoy ;

	ARM_Date payDate;
	Period tenor ;
	ARM_Model* model	= GetModel();

	coupon	=	*itsLeg.begin() ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
	int nResets = coupon->GetResDates().size() * nFlows;

	if (itsExpectedFwdRates ) 
		delete itsExpectedFwdRates ;
	itsExpectedFwdRates = new ARM_GP_Vector( nResets, 0.) ;
	if (itsDemCPIs ) 
		delete itsDemCPIs ;
	itsDemCPIs = new ARM_GP_Vector( nResets, 0.);
	if (itsNumCPIs ) 
		delete itsNumCPIs ;
	itsNumCPIs = new ARM_GP_Vector( nResets, 0.);
	if (itsDiscountFactors) 
		delete itsDiscountFactors ;
	itsDiscountFactors	= new ARM_GP_Vector(nFlows, 0.) ;	
	ARM_Date demResetDate, resDate ;
	ARM_GP_T_Vector<ARM_Date> resetDates ;
	ARM_GP_T_Vector<ARM_Date>::const_iterator resetDate ;
	
	for( i=0; i<nFlows; ++i )
	{
		coupon	=	itsLeg[i] ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		yoy			= coupon->GetIndex() ;
		if (yoy.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		tenor		= yoy->GetTenor() ;
		cpi			= yoy->GetIndex() ;
		payDate		= coupon->GetPayDate() ;
		
		resetDates	= coupon->GetResDates() ;
		resetDate	= resetDates.begin() ;
		
		if( (*itsPaymentTimes)[i] < 0. )	
			(*itsDiscountFactors)[i]	= 0.;
		else
			(*itsDiscountFactors)[i]= model->ZeroPrice( 0.0, (*itsPaymentTimes)[i], GetDiscountingYC() );
		
		for (; resetDate!=resetDates.end(); ++resetDate) 
		{
			resDate					= *resetDate ;
			demResetDate			= resDate.AddPeriodMult(tenor.GetUnit(), -tenor.GetLength());
			CPIRatio				= pricingMod->FwdCPIRatio(	*resetDate,	demResetDate, payDate,
												1,0.,itsInterpType,	itsFirstReset,	cpi.operator->() );
			
			tmp						= GetRcvOrPay() ;
			tmp						*= GetInterestTerms()->Elt(i) ; 
			tmp						*= CPIRatio[0] ;
			(*itsExpectedFwdRates)[i]	=  tmp ;
			(*itsNumCPIs)[i]			= CPIRatio[1];
			(*itsDemCPIs)[i]			= CPIRatio[2];

		}
		
	}

	ARM_Vector tmpFwds2 = To_ARM_Vector(*itsExpectedFwdRates);
	SetFwdRates(&tmpFwds2);

}

void InfComposedLeg::CptCashFlowValues()
{
	/*
		RMQ : s'il y a plusieurs resets par paiement, on paie la moyenne arithmethique des fixings
	*/
	double fwd = 0., nominal = 0., price = 0. ;
	int nFlows = itsLeg.size(), i;
	ARM_GP_Vector tmpFlowValues( nFlows, 0.); 
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;
	
	ARM_GP_Vector::const_iterator fwdIte = (*itsExpectedFwdRates).begin() ;
	
	ARM_GP_T_Vector<ARM_Date> resetDates ;
	ARM_GP_T_Vector<ARM_Date>::const_iterator resetDate ;
	
	for( i=0; i<nFlows; ++i, ++fwdIte )
	{
		coupon		=	itsLeg[i] ;
		nominal		=	coupon->GetNominal() ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		if( (*itsPaymentTimes)[i] < 0. )	
		{
			tmpFlowValues[i]	= 0.;
		}
		else
		{
			resetDates	= coupon->GetResDates() ;
			resetDate	= resetDates.begin() ;
			for (; resetDate!=resetDates.end(); ++resetDate, ++fwdIte) 
			{
				fwd += *fwdIte ;
			}
			
			tmpFlowValues[i]		= fwd*nominal*(*itsAccrualTimes)[i]*(*itsDiscountFactors)[i];
			tmpFlowValues[i]		/= (resetDates.end()-resetDates.begin());
		}
		price += tmpFlowValues[i] ;


	}

	SetDecompRates(GetFwdRates());		
	ARM_Vector* tmpFlowValues2 = To_pARM_Vector(&tmpFlowValues);
    SetCashFlowValues(tmpFlowValues2);	
	SetPrice(price) ;

}

///////////////////////////////////////////
// needed interface for ARM_InfIrIndex
///////////////////////////////////////////

ARM_GP_Vector* InfComposedLeg::GetNumJulianDates() const 
{
	return itsNumJulianDates  ;
}

ARM_GP_Vector* InfComposedLeg::GetDemJulianDates()  const 
{
	return itsDemJulianDates ;
}

string InfComposedLeg::GetIndexName() const 
{
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	
	if (coupon.IsNull()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "InfComposedLeg::GetIndexName() null coupon");
	
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;
	return cpiperf->GetIndex()->GetIndexName() ;
}


ARM_CountedPtr<InfPerformance> InfComposedLeg::GetIndex() const 
{
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "InfComposedLeg::GetIndexName() null coupon");
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;	
	return cpiperf ;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/


