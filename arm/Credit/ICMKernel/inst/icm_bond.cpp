
#include "ARMKernel\glob\firsttoinc.h"
#include "ICMKernel\inst\icm_bond.h"
#include "ARMKernel\mod\bondmath.h"
#include "ARMKernel\util\fromto.h"

void ICM_Bond::BitwiseCopy(const ARM_Object* srcleg)
{
    ICM_Bond* leg = (ICM_Bond *) srcleg;
}

void ICM_Bond::Copy(const ARM_Object* srcleg)
{
     ICM_Leg::Copy(srcleg);
 
     BitwiseCopy(srcleg);
}

ARM_Object* ICM_Bond::Clone(void)
{
     ICM_Bond* theClone = new ICM_Bond();

     theClone->Copy(this);
 
     return(theClone);
}

void ICM_Bond::Init(void)
{
	SetName(ICM_BOND);
}

void ICM_Bond::Set(double CouponRate,
			 const ARM_Date &Int_Accrual_Date,
			 const ARM_Date &Maturity, 
			 const ARM_Date* First_Period_Reference_Date,
			 int	  frequency,
			 double NotionalAmount , 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
			 const std::string& discountCcy,
			 const std::string& payCalendar,
 			 int	DayCount , 
			 int	AccruedDayCount , 
			 int    SettlementGap   , 
			 int	stubrule ,
			 double RedemptionValue)
{

	ICM_Leg::Set(Int_Accrual_Date, 
			Maturity, 
			First_Period_Reference_Date,
			0,
			CouponRate,
			AccruedOnDefault,  //AccruedOnDefault
			AccruedDayCount,
			0., // LastIndexFixing
			K_RCV, // rcvOrPay
			frequency, 
			DayCount, // dayCount
			K_COMP_PROP, // decompFreq
			K_ARREARS, // payTiming
			K_ADJUSTED, //intRule
			stubrule, //stubRule
			discountCcy, // discountCcy
			payCalendar, //payCalName
			K_NX_END,
			EXCLUDE_MATURITY,// 			 const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
			 K_ADJUSTED, // const int& adjStartDate /* = K_ADJUSTED*/ ,
			 qRunning_Leg, // const qCredit_Leg_Type& /* LegType = qRunning_Leg*/ ,
			 CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
			 ISSUER_UNDEFINE // const string& name /* = ISSUER_UNDEFINE*/ );
			);  

	SetName(ICM_BOND);

	// le SettlementGap est setté au niveau d'arm_security
	SetSettlementGap(SettlementGap);

	ARM_ReferenceValue* refval = new ARM_ReferenceValue(NotionalAmount , 1 /* price */, 0 /* K_CONSTANT */);

	SetAmount(refval,RedemptionValue);

	ARM_Vector* dummyFwd=new ARM_Vector(GetNumFlows(),0.);
	SetFwdRates(dummyFwd);

	if (refval)
		delete refval;
	refval = NULL;

}

// *********************************************************************************
// Methode Set pour ICM_Bond avec ARM_ReferenceValue* NotionalSchedule (amortization)
// *********************************************************************************

void ICM_Bond::Set(double CouponRate,
			 const ARM_Date &Int_Accrual_Date,
			 const ARM_Date &Maturity, 
			 const ARM_Date* First_Period_Reference_Date ,
			 int	  frequency,
			 double NotionalAmount , 
			 qPAYMENT_PREMIUM_LEG AccruedOnDefault,
			 const std::string& discountCcy,
			 const std::string& payCalendar,
			 int	DayCount , /*< Si non affecté il sera attribué par défaut*/
			 int	AccruedDayCount , /*< Si non affecté il sera attribué par défaut de la devise */
			 int    SettlementGap   ,
			 int	stubrule , /*<  Si non affecté il sera attribué par défaut de la devise */ 
			 ARM_ReferenceValue* NotionalSchedule)
{

	ICM_Leg::Set(Int_Accrual_Date, 
				 Maturity, 
				 First_Period_Reference_Date ,
				 0,
				 CouponRate,														
				 AccruedOnDefault,  //AccruedOnDefault
				 AccruedDayCount,
				 0., // LastIndexFixing
				 K_RCV, // rcvOrPay
				 frequency, 
				 DayCount, // dayCount
				 K_COMP_PROP, // decompFreq
				 K_ARREARS, // payTiming
				 K_UNADJUSTED, //intRule
				 stubrule, //stubRule
				 discountCcy, // discountCcy
				 payCalendar, //payCalName
				 K_NX_END,  //nxChange
			EXCLUDE_MATURITY,// 			 const bool& includematurity /* = EXCLUDE_MATURITY*/ ,
			 K_ADJUSTED, // const int& adjStartDate /* = K_ADJUSTED*/ ,
			 qRunning_Leg, // const qCredit_Leg_Type& /* LegType = qRunning_Leg*/ ,
			 CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/ ,
			 ISSUER_UNDEFINE // const string& name /* = ISSUER_UNDEFINE*/ 
			 );

				 
	SetName(ICM_BOND);

	// le SettlementGap est setté au niveau d'arm_security
	SetSettlementGap(SettlementGap);

	ARM_ReferenceValue* refval = new ARM_ReferenceValue(NotionalAmount , 1 /* price */, 0 /* K_CONSTANT */);

	SetAmount(refval,100.);

	SetAmortization(NotionalSchedule);

	ARM_Vector* dummyFwd=new ARM_Vector(GetNumFlows(),0.);
	SetFwdRates(dummyFwd);

	if (refval)
		delete refval;
	refval = NULL;
}




/*----------------------------------------------------------------------------*
    Compute the coupon for bond on a given period
*----------------------------------------------------------------------------*/ 
double ICM_Bond::CouponPeriod(int index)
{
	double dt = 1/(double)GetPaymentFreq();
	double Tol = 2.;

	ARM_Vector* CouponDates = GetFlowEndDates();
	int indxfin = CouponDates->GetSize()-1;
	ARM_Date LastCouponDate = CouponDates->Elt(indxfin);

	if (index<0) return 0.; //After CashFlowEndDates

	double Notional = GetAmount()->CptReferenceValue(LastCouponDate.GetJulian());
	double coupon = Notional*GetFixedRate()/100.;
	double yt_period = CountYears(GetDayCount(), (ARM_Date)(GetFlowStartDates()->Elt(index)), (ARM_Date)(GetFlowEndDates()->Elt(index)));
	double coupon_period = coupon* yt_period;

	if ((index>0) && (index<indxfin)) 
	{
		coupon_period = coupon*	dt;	
	}
	else if (fabs(GetInterestDays()->Elt(index) - 365.*dt)<= Tol) 
	{
		coupon_period = coupon*	dt;	
	}

	return (coupon_period);
}


/*----------------------------------------------------------------------------*
    Compute flows (from forward rates) to be paid and fill the cash flow values
*----------------------------------------------------------------------------*/ 
void ICM_Bond::CptCashFlowValues(void)
{
	ARM_SwapLeg::CptCashFlowValues();

    int i, nFlows;
    double rcvPay;

    double payJulDate;

    rcvPay = (double) GetRcvOrPay();

    nFlows = GetNumFlows();

    ARM_Vector* flowValues = new ARM_Vector(GetCashFlowValues());

    double* paymentDatesElt = GetPaymentDates()->GetElt();
    double* flowValuesElt   = flowValues->GetElt();
	ARM_Vector* decompRates = GetDecompRates();

    for (i = 0; i < nFlows; i++) 
    {
        payJulDate = paymentDatesElt[i];
        
		flowValuesElt[i] = CouponPeriod(i);

    }
	
    SetCashFlowValues(flowValues);
}


/*----------------------------------------------------------------------------*
    PriceToYield computes the yield of a bond for a given price and AsOfDate
	Uses a simple dichotomy
*----------------------------------------------------------------------------*/ 
double ICM_Bond::PriceToYield(ARM_Date& settlement, double price)
{
	if(GetAmount()->GetCalcMethod()!=K_CONSTANT_REF) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield on amortizing bond");
	
	double NotionalAmount=GetAmount()->CptReferenceValue(settlement);

	if (settlement < ARM_SwapLeg::GetStartDate())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield for settlement inferior to bond's start date");

	if (settlement > GetMaturity())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield for settlement superior to bond's maturity date");

	ARM_Vector* vector = new ARM_Vector(GetNumFlows(),0.);
	SetDecompRates(vector);
	CptCashFlowValues();
	delete vector;

	int i,index;
    double	r,
			dt(0.),
			presentValue(0.);

	// search for active period
	index = PeriodIndex(settlement);

	if (index<0) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot find period, check bond's construction!");
	int indxfin = GetPaymentDates()->GetSize() - 1;

    double refPV=(price+Accrued(settlement))/NotionalAmount;

	//price computed with yield=0
	presentValue=(GetRedemption()/NotionalAmount);
	for(i=indxfin;i>index;i--)
	{
		presentValue+=(GetCashFlowValues()->Elt(i)/NotionalAmount)/GetPaymentFreq();
	};
	presentValue+=(GetCashFlowValues()->Elt(index)/NotionalAmount)/GetPaymentFreq();

	if(presentValue==refPV) return 0.0;

	//check wether yield will be positive or negative
	double refcomp(1.0);
	if(presentValue<refPV) refcomp=-1.0;
	
	double yield_min(0.0),
		yield_max=refcomp;
	
	r=0.01*(yield_max/GetPaymentFreq());
	// present value	
	presentValue=(GetRedemption()/NotionalAmount);
	for(i=indxfin;i>index;i--)
	{
		dt=CountYears(GetDayCount(),GetFlowEndDates()->Elt(i-1),GetFlowEndDates()->Elt(i))*GetPaymentFreq();
		presentValue+=(GetCashFlowValues()->Elt(i)/NotionalAmount)/GetPaymentFreq();
		presentValue*=pow(1./(1.+r),dt);
	};
	dt=CountYears(GetDayCount(),settlement.GetJulian(),GetFlowEndDates()->Elt(index))*GetPaymentFreq();
	presentValue+=(GetCashFlowValues()->Elt(index)/NotionalAmount)/GetPaymentFreq();
	presentValue*=pow(1./(1.+r),dt);
	//look for integer upper_bound for yield (or lower_bound if negative)
	while((refcomp*presentValue)>(refcomp*refPV))
	{
		yield_min+=refcomp;
		yield_max+=refcomp;
		if(yield_max>100.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield, price is too low, yield>100");
		if(yield_max<-99.0) throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield, price is too high, yield<(-99)");
		
		r=0.01*(yield_max/GetPaymentFreq());
		// present value	
		presentValue=(GetRedemption()/NotionalAmount);
		for(i=indxfin;i>index;i--)
		{
			dt=CountYears(GetDayCount(),GetFlowEndDates()->Elt(i-1),GetFlowEndDates()->Elt(i))*GetPaymentFreq();
			presentValue+=(GetCashFlowValues()->Elt(i)/NotionalAmount)/GetPaymentFreq();
			presentValue*=pow(1./(1.+r),dt);
		};
		dt=CountYears(GetDayCount(),settlement.GetJulian(),GetFlowEndDates()->Elt(index))*GetPaymentFreq();
		presentValue+=(GetCashFlowValues()->Elt(index)/NotionalAmount)/GetPaymentFreq();
		presentValue*=pow(1./(1.+r),dt);
	}
	
	//dichotomy within integer boundaries
	double yield;
	while(fabs(yield_max-yield_min)>1.0e-5)
	{
		yield=(yield_min+yield_max)/2.0;

		r=0.01*(yield/GetPaymentFreq());
		// present value	
		presentValue=(GetRedemption()/NotionalAmount);
		for(i=indxfin;i>index;i--)
		{
			dt=CountYears(GetDayCount(),GetFlowEndDates()->Elt(i-1),GetFlowEndDates()->Elt(i))*GetPaymentFreq();
			presentValue+=(GetCashFlowValues()->Elt(i)/NotionalAmount)/GetPaymentFreq();
			presentValue*=pow(1./(1.+r),dt);
		};
		dt=CountYears(GetDayCount(),settlement.GetJulian(),GetFlowEndDates()->Elt(index))*GetPaymentFreq();
		presentValue+=(GetCashFlowValues()->Elt(index)/NotionalAmount)/GetPaymentFreq();
		presentValue*=pow(1./(1.+r),dt);

		if((refcomp*presentValue)>(refcomp*refPV))
		{
			yield_min=yield;
		}
		else
		{
			yield_max=yield;
		}
	}

    return((yield_min+yield_max)/2.0);
}


/*----------------------------------------------------------------------------*
    YieldToPrice computes the price of a bond for a given yield and AsOfDate
*----------------------------------------------------------------------------*/ 
double ICM_Bond::YieldToPrice(ARM_Date& settlement, double yield)
{
	if(yield<-99.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Yield must be superior to -99");
	
	if(yield>100.0)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Yield must be inferior to 100");

	if(GetAmount()->GetCalcMethod()!=K_CONSTANT_REF) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield on amortizing bond");
	
	CptCashFlowValues();
	
	double NotionalAmount=GetAmount()->CptReferenceValue(settlement);

	if (settlement < ARM_SwapLeg::GetStartDate())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield for settlement inferior to bond's start date");

	if (settlement > GetMaturity())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot compute yield for settlement superior to bond's maturity date");

	int index = 0;
    double	dt(0.),
			presentValue(0.);

	// search for active period
	index = PeriodIndex(settlement);
	if (index<0) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
			"Cannot find period, check bond's construction!");
	int indxfin = GetPaymentDates()->GetSize() - 1;

	// present value	
	presentValue=(GetRedemption()/NotionalAmount);
	double r=0.01*(yield/GetPaymentFreq());

	for(int i=indxfin;i>index;i--)
	{
		dt=CountYears(GetDayCount(),GetFlowEndDates()->Elt(i-1),GetFlowEndDates()->Elt(i))*GetPaymentFreq();
		presentValue+=(GetCashFlowValues()->Elt(i)/NotionalAmount)/GetPaymentFreq();
		presentValue*=pow(1./(1.+r),dt);
	};

	dt=CountYears(GetDayCount(),settlement.GetJulian(),GetFlowEndDates()->Elt(index))*GetPaymentFreq();
	presentValue+=(GetCashFlowValues()->Elt(index)/NotionalAmount)/GetPaymentFreq();
	presentValue*=pow(1./(1.+r),dt);

    return ((presentValue*NotionalAmount)-Accrued(settlement));
}


/*----------------------------------------------------------------------------*
    Accrued computes the accrued of a bond for a given AsOfDate
*----------------------------------------------------------------------------*/ 
double ICM_Bond::Accrued(ARM_Date & AsOfDate)
{
    double previousSettl =0., accrued=0.;

	ARM_Date IssueDate, previousCouponDate;

	ARM_Date Settlement = AsOfDate;
	//Settlement.GapBusinessDay(GetSettlementGap());

    int dayCount = KACTUAL_360;//GetAccruedDayCount();
	int payfreq = GetPaymentFreq();

	ARM_Vector* CouponDate = GetFlowEndDates();
	int NumFlows = GetNumFlows();

	ARM_Date FirstCouponDate = CouponDate->Elt(0);
	ARM_Date LastCouponDate = CouponDate->Elt(NumFlows-1);


	IssueDate = ARM_SwapLeg::GetStartDate();
	int index = PeriodIndex(Settlement);
	double Coupon = CouponPeriod(index);
	    
	if ( Coupon < K_DOUBLE_TOL )
			return(0.0);
		
	if (( Settlement < FirstCouponDate ) && (Settlement > IssueDate))
	    previousSettl = CountYears(dayCount, IssueDate, Settlement);
	else if (( Settlement > FirstCouponDate ) && ( Settlement < LastCouponDate ))
	{
		previousCouponDate = CouponDate->Elt(index-1);
		previousSettl = CountYears(dayCount, previousCouponDate, Settlement);
	}
	else
		previousSettl = 0.;
	
    accrued = Coupon * previousSettl;

	return(accrued);
}


