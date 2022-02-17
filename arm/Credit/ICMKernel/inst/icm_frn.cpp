#pragma warning (disable:4786)
#include "ICMKernel\inst\ICM_frn.h"

void ICM_Frn::BitwiseCopy(const ARM_Object* srcleg)
{
    ICM_Frn* leg = (ICM_Frn *) srcleg;
}

void ICM_Frn::Copy(const ARM_Object* srcleg)
{
     ICM_Leg::Copy(srcleg);
 
     BitwiseCopy(srcleg);
}

ARM_Object* ICM_Frn::Clone(void)
{
     ICM_Frn* theClone = new ICM_Frn();

     theClone->Copy(this);
 
     return(theClone);
}

void ICM_Frn::Init(void)
{
	SetName(ICM_FRN);
}

void ICM_Frn:: Set( double Spread,
					const ARM_Date& Int_Accrual_Date,  
					const ARM_Date& Maturity, 
					const ARM_Date* First_Period_Reference_Date,
					// ARM_Date First_Period_Reference_Date, //Reference Date
					ARM_IRIndex* irindex,
					double InitialRate , // InitialRate	
					double LastIndexFixing,
					qPAYMENT_PREMIUM_LEG AccruedOnDefault,
					double Notional ,
					int AccruedDayCount ,
					int SettlementGap ,
					int DayCount , // dayCount
					int stubrule , //stubRule
					const std::string& Devise , // discountCcy
					const std::string&  resetCalName ,
					const std::string& payCalName ,
					int nxChange)
{

	ICM_Leg::Set(Int_Accrual_Date, 
				 Maturity, 
				 First_Period_Reference_Date,
				 0,
				 irindex,
				 // First_Period_Reference_Date, //Reference Date
				 Spread,
				 AccruedOnDefault,  //AccruedOnDefault
				 AccruedDayCount,
				 InitialRate, // InitialRate	
				 LastIndexFixing, // LastIndexFixing
				 K_RCV, // rcvOrPay
				 DayCount, // dayCount
				 K_COMP_PROP, // decompFreq
				 stubrule, //stubRule
				 10000, //resetgap
				 resetCalName, //resetCalName
			     Devise, // discountCcy
				 payCalName, //payCalName
				 nxChange,
				 EXCLUDE_MATURITY, // const bool& includematurity /* = EXCLUDE_MATURITY*/,
				 K_ADJUSTED, // const int& adjStartDate /* = K_ADJUSTED*/,
				 qRunning_Leg, // const qCredit_Leg_Type& LegType /* = qRunning_Leg*/,
				 CREDIT_DEFAULT_VALUE, // const double& Binary /* = CREDIT_DEFAULT_VALUE*/,
				 ISSUER_UNDEFINE // const string& name /* = ISSUER_UNDEFINE*/ );
				 ); 


	SetName(ICM_FRN);

	ARM_ReferenceValue* refval = new ARM_ReferenceValue(Notional , 1 /* price */, 0 /* K_CONSTANT */);

	SetAmount(refval,100. /*Redemption Value*/);

	if (refval)
		delete refval;
	refval = NULL;

}


/*----------------------------------------------------------------------------*
    Compute the coupon for frn a given period
*----------------------------------------------------------------------------*/ 
double ICM_Frn::CouponPeriod(int index)
{
	double coupon_period = GetCashFlowValues()->Elt(index);

	return (coupon_period);
}


/*----------------------------------------------------------------------------*
    Compute flows (from forward rates) to be paid and fill the cash flow values
*----------------------------------------------------------------------------*/ 
void ICM_Frn::CptCashFlowValues(void)
{
	ARM_SwapLeg::CptCashFlowValues();
}

/*----------------------------------------------------------------------------*
    Accrued computes the accrued of a bond for a given AsOfDate
*----------------------------------------------------------------------------*/ 
double ICM_Frn::Accrued(ARM_Date & AsOfDate)
{
    double previousSettl =0., accrued=0.;

	ARM_Date IssueDate, FirstCouponDate,previousCouponDate;

	ARM_Date Settlement = AsOfDate;

    int dayCount = KACTUAL_360; //GetAccruedDayCount();
	int payfreq = GetPaymentFreq();

	ARM_Vector* CouponDate = GetFlowEndDates();
	FirstCouponDate = CouponDate->Elt(0);
	int NumFlows = GetNumFlows();

	ARM_Date LastCouponDate = CouponDate->Elt(NumFlows-1);

	if (GetModel())
	{
		//accrued = AccruedCoupon(AsOfDate,0);
	}
	else
	{
		IssueDate = ARM_SwapLeg::GetStartDate();
		int index = PeriodIndex(Settlement);
		double Coupon = CouponPeriod(index);
	    
		if ( Coupon < K_DOUBLE_TOL )
				return(0.0);
		
		if ( Settlement < FirstCouponDate )
		    previousSettl = CountYears(dayCount, IssueDate, Settlement);
		else if (( Settlement > FirstCouponDate ) && ( Settlement < LastCouponDate ))
		{
			previousCouponDate = CouponDate->Elt(index-1);
			previousSettl = CountYears(dayCount, previousCouponDate, Settlement);
		}
		else
			previousSettl = 0.;
	
	    accrued = Coupon * previousSettl;
	}

	return(accrued);
}


