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
#include "gpbase/checkarg.h"




CC_BEGIN_NAMESPACE( ARM )

InfLeg::InfLeg(	const ARM_GP_T_Vector<ARM_CountedPtr<CashFlow> >& vCahFlow,
				//for backward compatibility : to remove
				int interpType,	double firstReset): 
Instrument<InfFloatingCoupon<YOYPayOff> >(vCahFlow),	
ARM_SwapLeg(), 
itsInterpType(interpType), itsFirstReset(firstReset)
{
	Init();

	int nFlows = itsLeg.size(),i, dayCount;
	double time = 0 ;
	ARM_CountedPtr<InfFloatingCoupon<YOYPayOff> > coupon ;

	ARM_CountedPtr<InfPerformance> yoy ;
	ARM_CountedPtr<ARM_InfIdx> cpi ;
	ARM_CountedPtr<ARM_InfCurv> cpiCurve ;


	ARM_Date resDate, payDate, today, tmp ;
	Period tenor ;
	ARM_Vector startDates(nFlows);
	ARM_Vector intDays(nFlows);
	ARM_Vector intTerms(nFlows);
	ARM_Vector endDates(nFlows);
	ARM_Vector paymentDates(nFlows); 
    ARM_Vector resetDates(nFlows); 

	itsResetTimes = new ARM_GP_Vector(nFlows) ;
	itsPaymentTimes = new ARM_GP_Vector(nFlows) ;
	itsAccrualTimes = new ARM_GP_Vector(nFlows) ;
	itsNumJulianDates = new ARM_GP_Vector(nFlows) ;
	itsDemJulianDates = new ARM_GP_Vector(nFlows) ;
	

	for( i=0; i<nFlows; ++i )
	{
		coupon	=	itsLeg[i] ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfLeg::InfLeg() NULL coupon	");
		yoy			= coupon->GetIndex() ;
		if (yoy.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfLeg::InfLeg() NULL YOY performance	");
		cpi			= yoy->GetIndex() ;
		if (cpi.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR : InfLeg::InfLeg() NULL CPI index	");
		cpiCurve  = cpi->GetInfFwdCurve(); 
		if (!cpiCurve.IsNull()) 
			today					= cpi->GetInfFwdCurve()->GetAsOf() ;
		/* 
			s'il n'y a pas de courbe infl today = today()  (cf. ARM_Date)
			TODO : exiger que ARM_InfIdx contienne une courbe infl
		*/

		tenor					= yoy->GetTenor() ;
		startDates[i]			= coupon->GetStartDate().GetJulian();
		endDates[i]				= (coupon->GetEndDate().GetJulian()) ;
		payDate					= coupon->GetPayDate() ;
		paymentDates[i]			= (payDate.GetJulian()) ;
		dayCount				= cpi->GetDayCount() ; 
		time					= CountYears(dayCount, today, payDate);
		(*itsPaymentTimes)[i]	= time;
		time					= CountYears(dayCount, coupon->GetStartDate(), coupon->GetEndDate()) ;		
		(*itsAccrualTimes)[i]	= time ;
		intTerms[i]				= (coupon->GetAccrualTime()) ;
		intDays[i]				= (coupon->GetAccrualDays()) ;

		(*itsNumJulianDates)[i] =  coupon->GetNumDate().GetJulian()  ;
		(*itsDemJulianDates)[i] = coupon->GetDemDate().GetJulian()  ;
		
		resDate		= coupon->GetNumDate() ;
		time		= CountYears(dayCount, today, resDate ) ;
		(*itsResetTimes)[i] = time ;

		resetDates[i] = (resDate.GetJulian()) ;
	}

	//for backward compatibility : to remove
	SetFlowStartDates(&startDates);
	SetFlowEndDates(&endDates);
	SetPaymentDates( new ARM_Vector(paymentDates) );
	SetResetDates(&resetDates);
	SetInterestTerms(&intTerms) ; 
	SetInterestDays(&intDays) ;
}	

InfLeg::InfLeg(const InfLeg& infLeg): 
Instrument<InfFloatingCoupon<YOYPayOff> >(infLeg.itsLeg),	ARM_SwapLeg(), 
itsNumJulianDates(infLeg.itsNumJulianDates ? (ARM_GP_Vector *) infLeg.itsNumJulianDates->Clone() : NULL), 
itsDemJulianDates(infLeg.itsDemJulianDates ? (ARM_GP_Vector *) infLeg.itsDemJulianDates->Clone() : NULL), 
itsPaymentTimes(infLeg.itsPaymentTimes ? (ARM_GP_Vector *) infLeg.itsPaymentTimes->Clone() : NULL),
itsAccrualTimes(infLeg.itsAccrualTimes ? (ARM_GP_Vector *) infLeg.itsAccrualTimes->Clone() : NULL),
itsResetTimes(infLeg.itsResetTimes ? (ARM_GP_Vector *) infLeg.itsResetTimes->Clone() : NULL)

{
	//for backward compatibility : to remove
	SetFlowStartDates(infLeg.GetFlowStartDates()) ;

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
		SetFlowStartDates(infLeg.GetFlowStartDates()) ;

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
/*
void  InfLeg::SetModel(ARM_CountedPtr<ARM_Model> model) 
{
	//itsInstrumentModel = model ;
	itsInstrumentModel = ARM_CountedPtr<ARM_Model>(new ARM_Model(model));
}
*/
double InfLeg::ComputePrice(int mode)
{
	calculate() ;
	return GetPrice() ;
}

void InfLeg::CptExpectedFwdRates()
{
	/// FwdCPIRatio returns in the following order
	/// CPIRatio[0] is the ratio of CPI
	/// CPIRatio[1] is the numerator CPI
	/// CPIRatio[2] is the denominator CPI

	InfFwdModel* pricingMod = dynamic_cast<InfFwdModel*> (GetModel());

	if( ! pricingMod) ARM_THROW( ERR_INVALID_ARGUMENT, " the model you use cannot price inflation derivatives");

	int nFlows = itsLeg.size();
	double firstReset = 0., tmp = 0., spread = 0.;

	if (itsExpectedFwdRates ) 
		delete itsExpectedFwdRates ;
	itsExpectedFwdRates = new ARM_GP_Vector( nFlows, 0.) ;
	if (itsDemCPIs ) 
		delete itsDemCPIs ;
	itsDemCPIs = new ARM_GP_Vector( nFlows, 0.);
	if (itsNumCPIs ) 
		delete itsNumCPIs ;
	itsNumCPIs = new ARM_GP_Vector( nFlows, 0.);
	if (itsDiscountFactors) 
		delete itsDiscountFactors ;
	itsDiscountFactors	= new ARM_GP_Vector(nFlows, 0.) ;

	ARM_GP_Vector CPIRatio;

	int i;
	ARM_CountedPtr<InfFloatingCoupon<YOYPayOff> > coupon ;
	ARM_CountedPtr<ARM_InfIdx> cpi ;
	ARM_CountedPtr<InfPerformance> yoy ;

	ARM_Date payDate;
	Period tenor ;
	ARM_Model* model	= GetModel();

	for( i=0; i<nFlows; ++i )
	{
		coupon	=	itsLeg[i] ;
		if (coupon.IsNull()) 
		   ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");		
		yoy			= coupon->GetIndex() ;
		if (yoy.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		tenor		= yoy->GetTenor() ;
		cpi			= yoy->GetIndex() ;
		if (cpi.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR");
		payDate		= coupon->GetPayDate() ;
		
		if( (*itsPaymentTimes)[i] < 0. )	
			(*itsDiscountFactors)[i]	= 0.;
		else
			(*itsDiscountFactors)[i]= model->ZeroPrice( 0.0, (*itsPaymentTimes)[i], GetDiscountingYC() );

		CPIRatio				= pricingMod->FwdCPIRatio(	coupon->GetNumDate(),	
															coupon->GetDemDate(), 
															payDate,1,0.,itsInterpType,	
															itsFirstReset,	&*cpi );
			
		tmp						= GetRcvOrPay() ;
		tmp						*= GetInterestTerms()->Elt(i) ; 
		tmp						*= CPIRatio[0] ;
		(*itsExpectedFwdRates)[i]	=  tmp ;
		(*itsNumCPIs)[i]			= CPIRatio[1];
		(*itsDemCPIs)[i]			= CPIRatio[2];
	}

	ARM_Vector tmpFwds2 = To_ARM_Vector(*itsExpectedFwdRates);
	SetFwdRates(&tmpFwds2);

}

void InfLeg::CptCashFlowValues()
{
	double fwd = 0., nominal = 0., price = 0. ;
	int nFlows = itsLeg.size(), i;
	ARM_GP_Vector tmpFlowValues( nFlows, 0.); 
	ARM_CountedPtr<InfFloatingCoupon<YOYPayOff> > coupon ;
	
	ARM_GP_Vector::const_iterator fwdIte = (*itsExpectedFwdRates).begin() ;

	for( i=0; i<nFlows; ++i, ++fwdIte )
	{
		coupon		=	itsLeg[i] ;
		nominal		=	coupon->GetNominal() ;
		if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "ERR nul coupon");
		if( (*itsPaymentTimes)[i] < 0. )	
		{
			tmpFlowValues[i]	= 0.;
		}
		else 
		{
			fwd					= *fwdIte ;
			tmpFlowValues[i]		= fwd*nominal*(*itsAccrualTimes)[i]*(*itsDiscountFactors)[i];
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

ARM_GP_Vector* InfLeg::GetNumJulianDates() const 
{
	return itsNumJulianDates  ;
}

ARM_GP_Vector* InfLeg::GetDemJulianDates()  const 
{
	return itsDemJulianDates ;
}

string InfLeg::GetIndexName() const 
{
	ARM_CountedPtr<InfFloatingCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "InfLeg::GetIndexName() null coupon");
	//InfPerformance* cpiperf = (InfPerformance*) coupon->GetIndex() ;
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;
	return cpiperf->GetIndex()->GetIndexName() ;
}
ARM_CountedPtr<InfPerformance> InfLeg::GetIndex() const 
{
	ARM_CountedPtr<InfFloatingCoupon<YOYPayOff> > coupon = itsLeg[0] ;
	if (coupon.IsNull()) ARM_THROW( ERR_INVALID_ARGUMENT, "InfLeg::GetIndexName() null coupon");
	ARM_CountedPtr<InfPerformance> cpiperf = coupon->GetIndex() ;

	//	InfPerformance* cpiperf = (InfPerformance*)coupon->GetIndex() ;
	return cpiperf ;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/


