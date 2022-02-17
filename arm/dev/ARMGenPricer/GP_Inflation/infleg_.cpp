/*!
 *
 * Copyright (c) NATIXIS May 2007 Paris
 *
 *	\file infleg.cpp
 *
 *	\author  Francois Poitou
 *	\version 1.0
 *	\date May 2007
 */


/// gpinflation
#include "gpinflation/infleg_.h"
#include "gpinflation/infPayOff.h"
#include "gpinflation/infFloatingCoupon.h"
/// gpbase
#include "gpbase/gplinalgconvert.h"




CC_BEGIN_NAMESPACE( ARM )

InfLeg::InfLeg(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow,
				//for backward compatibility : to remove
				int interpType,	double firstReset): 
Instrument(vCahFlow),	ARM_SwapLeg(), itsInterpType(interpType), itsFirstReset(firstReset)
{
	int nFlows = itsLeg.size(),i, dayCount;
	double time = 0 ;
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;
	ARM_CountedPtr<InfPerformance > yoy ;
	ARM_CountedPtr<ARM_InfIdx > cpi ;

	ARM_Date demResetDate, resDate, payDate, today, tmp ;
	Period tenor ;
	ARM_GP_T_Vector<ARM_Date> vResetDates ;
	ARM_GP_T_Vector<ARM_Date>::const_iterator resetDate ;

	ARM_Vector startDates, intDays, intTerms, endDates, paymentDates, resetDates; //to remove

	for( i=0; i<nFlows; ++i )
	{
		coupon	=	itsLeg[i] ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		yoy			= coupon->GetIndex() ;
		if (yoy.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		cpi			= yoy->GetIndex() ;

		tenor		= yoy->GetTenor() ;
		vResetDates	= coupon->GetResDates() ;
		resetDate	= vResetDates.begin() ;
		startDates	.push_back(coupon->GetStartDate()	.GetJulian()) ;
		endDates	.push_back(coupon->GetEndDate()		.GetJulian()) ;

		payDate		= coupon->GetPayDate() ;
		paymentDates.push_back(payDate		.GetJulian()) ;
		
		dayCount	= cpi->GetDayCount() ; 
		today		= cpi->GetSettlement() ;
		time		= CountYears(dayCount, today, payDate) ;
		itsPaymentTimes.push_back(time) ;
		time		= CountYears(dayCount, coupon->GetStartDate(), coupon->GetEndDate()) ;
		itsAccrualTimes.push_back(time) ;

		intTerms	.push_back(coupon->GetAccrualTime()) ;
		intDays		.push_back(coupon->GetAccrualDays()) ;
		
		for (; resetDate!=vResetDates.end(); ++resetDate) 
		{
			resDate				= *resetDate ;
			tmp					= *resetDate ;
			demResetDate		= tmp.AddPeriodNoAdj(tenor.GetUnit(), -tenor.GetLength());
			resetDates			.push_back((*resetDate).GetJulian()) ;
			itsNumJulianDates	.push_back( (*resetDate).GetJulian() ) ;
			itsDemJulianDates	.push_back( demResetDate.GetJulian()  ) ;
			time				= CountYears(dayCount, today, resDate ) ;
			itsResetTimes		.push_back(time) ;
		
		}
	}

	//for backward compatibility : to remove
	SetFlowStartDates( new ARM_Vector(startDates) ) ;
	SetFlowEndDates( new ARM_Vector(endDates) ) ;
	SetPaymentDates( new ARM_Vector(paymentDates) );
	SetResetDates( new ARM_Vector(resetDates) );
	SetInterestTerms(new ARM_Vector(intTerms)) ; 
	SetInterestDays(new ARM_Vector(intDays)) ;
}	

InfLeg::InfLeg(	const InfLeg& infLeg): 
Instrument(infLeg.itsLeg),	ARM_SwapLeg(), itsNumJulianDates(infLeg.itsNumJulianDates), 
itsDemJulianDates(infLeg.itsDemJulianDates), itsPaymentTimes(infLeg.itsPaymentTimes),
itsAccrualTimes(infLeg.itsAccrualTimes),itsResetTimes(infLeg.itsResetTimes)

{
	//for backward compatibility : to remove
	SetFlowStartDates( new ARM_Vector(infLeg.GetFlowStartDates() ) ) ;

}

InfLeg&	
InfLeg::operator=(	const InfLeg& infLeg)
{
	if( this !=	 &infLeg )
	{
		ARM_SwapLeg::operator = ( infLeg );
		itsLeg = infLeg.itsLeg ;
		itsDemJulianDates	= infLeg.itsDemJulianDates ;
		itsNumJulianDates	= infLeg.itsNumJulianDates ;
		itsPaymentTimes		= infLeg.itsPaymentTimes ;
		itsAccrualTimes		= infLeg.itsAccrualTimes ;
		itsResetTimes		= infLeg.itsResetTimes ;

		//for backward compatibility : to remove
		SetFlowStartDates( new ARM_Vector(infLeg.GetFlowStartDates() ) ) ;

	}
	return *this;
}
	
/* 
 TODO : obliger à utiliser setModel() avant les calculate() ... 
 comme dans Instrument::Computeprice(Model)
*/
void InfLeg::performCalculation() 
{
    CptExpectedFwdRates();
    CptCashFlowValues();
}
void  InfLeg::SetModel(ARM_CountedPtr<ARM_Model> model) 
{
	itsInstrumentModel = model ;
}
void InfLeg::ComputePrice()
{
	calculate() ;
}

void InfLeg::CptExpectedFwdRates()
{
	ARM_CountedPtr<InfFwdModel>  pricingMod  = itsInstrumentModel ;

	if( pricingMod.IsNull() ) ARM_THROW( ERR_INVALID_ARGUMENT, " the model you use cannot price inflation derivatives");

	int nFlows = itsLeg.size();
	double firstReset = 0. ;

	ARM_GP_Vector itsExpectedFwdRates( nFlows, 0.);
	ARM_GP_Vector CPIRatio;

	int i;
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;
	ARM_CountedPtr<ARM_InfIdx > cpi ;
	ARM_CountedPtr<InfPerformance > yoy ;

	ARM_Date payDate, demResetDate, resDate ;
	Period tenor ;
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
		resetDates	= coupon->GetResDates() ;
		resetDate	= resetDates.begin() ;
		payDate		= coupon->GetPayDate() ;

		for (; resetDate!=resetDates.end(); ++resetDate) 
		{
			/// FwdCPIRatio returns in the following order
			/// CPIRatio[0] is the ratio of CPI
			/// CPIRatio[1] is the numerator CPI
			/// CPIRatio[2] is the denominator CPI

			resDate					= *resetDate ;
			demResetDate			= resDate.AddPeriodMult(tenor.GetUnit(), -tenor.GetLength());
			CPIRatio				= pricingMod->FwdCPIRatio(	*resetDate,	demResetDate, payDate,
												1,0.,itsInterpType,	itsFirstReset,	cpi.operator->() );
			itsExpectedFwdRates[i]	= GetRcvOrPay() * (CPIRatio[0] * CC_NS( ARM_Constants, rateBase ) + GetSpread())* GetInterestTerms()->Elt(i);
			itsNumCPIs[i]			= CPIRatio[1];
			itsDemCPIs[i]			= CPIRatio[2];

		}
	}

	ARM_Vector tmpFwds2 = To_ARM_Vector(itsExpectedFwdRates);
	SetFwdRates(&tmpFwds2);

}

void InfLeg::CptCashFlowValues()
{
	/*
		RMQ : s'il y a plusieurs resets par paiement, on paie la moyenne arithmethique des fixings
	*/
	double fwd= 0., nominal = 0., price = 0. ;
	int nFlows = itsLeg.size(), i;
	ARM_GP_Vector tmpFlowValues( nFlows, 0.0); 
	ARM_GP_Vector::const_iterator fwdIte = (*itsFwdRates).begin() ;
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon ;

	ARM_GP_T_Vector<ARM_Date> resetDates ;
	ARM_GP_T_Vector<ARM_Date>::const_iterator resetDate ;
	itsDiscountFactors = ARM_GP_Vector(nFlows, 0.) ;
	ARM_Model* model			 = GetModel();
	if(! model->GetZeroCurve()->GetDiscountFactors() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");

	for( i=0; i<nFlows; ++i )
	{
		coupon		=	itsLeg[i] ;
		nominal		=	coupon->GetNominal() ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		resetDates	= coupon->GetResDates() ;
		resetDate	= resetDates.begin() ;
		if( itsPaymentTimes[i] < 0. )	{
			itsDiscountFactors[i]	= 0.0;
			tmpFlowValues[i]	= 0.0;
		}
		else
			itsDiscountFactors[i]	= model->ZeroPrice( 0.0, itsPaymentTimes[i], GetDiscountingYC() );
		
		for (; resetDate!=resetDates.end(); ++resetDate, ++fwdIte) 
		{
			fwd += *fwdIte ;
		}
		tmpFlowValues		= fwd*nominal*itsAccrualTimes[i]*itsDiscountFactors[i];
		tmpFlowValues		/= (resetDates.end()-resetDates.begin());
		price				+= tmpFlowValues[i] ;

	}

	////PFFFFF...
	SetDecompRates(GetFwdRates());		
	ARM_Vector* tmpFlowValues2 = To_pARM_Vector(&tmpFlowValues);
	//delete tmpFlowValues;
    SetCashFlowValues(tmpFlowValues2);	
	SetPrice(price) ;
	////fin PFFFFF...

}

///////////////////////////////////////////
// needed interface for ARM_InfIrIndex
///////////////////////////////////////////

ARM_GP_Vector* InfLeg::GetNumJulianDates() const 
{
	return new ARM_GP_Vector(itsNumJulianDates ) ;
}

ARM_GP_Vector* InfLeg::GetDemJulianDates()  const 
{
	return new ARM_GP_Vector(itsDemJulianDates) ;
}

string InfLeg::GetIndexName() const 
{
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;
	return cpiperf->GetIndex()->GetIndexName() ;
}
ARM_CountedPtr<InfPerformance> InfLeg::GetIndex() const 
{
	ARM_CountedPtr<InfComposedCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;
	return cpiperf ;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/


