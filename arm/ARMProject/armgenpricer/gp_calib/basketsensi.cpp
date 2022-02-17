/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

/// gpcalib
#include "gpcalib/basketsensi.h"
#include "gpcalib/stripper.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/vanillaspreadoption.h"


/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"

// kernel
//#include <inst/spreadoption.h>



CC_BEGIN_NAMESPACE( ARM )

const double SENSI_CMS		= 0.00001;

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensi
///	Routine: FillCpnIndex
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_BasketSensi::FillCpnIndex(const std::vector<double>& exerciseDates,const std::vector<double>& startDates)
{
	double days	= 7.;

	int i = 0;
	size_t k;

	for (k=0; k<exerciseDates.size(); k++)
	{
		if (itsStartDates[itsStartDates.size()-1]+days>=startDates[k])
		{

			while (startDates[k]>itsStartDates[i]+days) 
				i++;

			if(exerciseDates[k]>itsResetDates[i]+days)
			{
				char msg[500];
				sprintf(msg, "ARM_BasketSensi : call date[%d] (%.0f) between Reset[%d] (%.0f) and Start[%d] (%.0f). Not allowed !", 
						k, exerciseDates[k], i, itsResetDates[i], i, itsStartDates[i]);
				ARM_THROW( ERR_INVALID_ARGUMENT, msg );
			}
			itsCpnIndex.push_back(i);
		}
		else
			itsCpnIndex.push_back(itsStartDates.size());
	}
	itsCpnIndex.push_back(itsStartDates.size());
		
	
	int prevIndex;
	itsRompu.resize(exerciseDates.size());
	itsFees.resize(exerciseDates.size());

	if (itsCpnIndex[exerciseDates.size()-1]==itsStartDates.size())
	{
		SetRompu(exerciseDates.size()-1,true);
		prevIndex = itsCpnIndex[exerciseDates.size()-1]-1;
		if (prevIndex<0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BasketSensi::FillCpnIndex : c'est le drame" );
		SetFee(exerciseDates.size()-1,ComputeFee(prevIndex,startDates[exerciseDates.size()-1]));
	}

	for (k=exerciseDates.size()-1; k>0; k--)
	{
		if (itsCpnIndex[k-1]>=itsCpnIndex[k])
		{
			SetRompu(k-1,true);
			prevIndex = itsCpnIndex[k-1]-1;
			if (prevIndex<0)
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BasketSensi::FillCpnIndex : c'est le drame" );
			SetFee(k-1,ComputeFee(prevIndex,startDates[k-1]));
		}
	}

	for (k=0;k<exerciseDates.size();k++)
		if (IsRompu(k))
			itsCpnIndex[k]--;


}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensi
///	Routine: UpdateRompu
///	Returns: void
///	Action : computes stub
////////////////////////////////////////////////////
double ARM_BasketSensi::ComputeFee(int i,double startDate)
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BasketSensi::ComputeFee : Freq (call) > Freq (Coupon) is forbidden" );
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiCSO
///	Routine: constructor
///	Returns: void
///	Action : computes sensis of SO coupon
////////////////////////////////////////////////////
ARM_BasketSensiCSO::ARM_BasketSensiCSO(ARM_ZeroCurve* curve,ARM_SpreadOption* so,double w, double p)
	:	ARM_BasketSensi(w,p)
{
	//double asOfDate	  = 0;//curve->GetAsOfDateJul();
	//ARM_Currency* ccy = NULL;//curve->GetCurrencyUnit();

	//ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)ARM_ConverterFromKernel::ConvertSecuritytoArgObject(so,asOfDate);
	//soVanillaArg->ComputeIndexSchedulesAndAdjustDates(ccy, asOfDate);
	//itsArg = soVanillaArg;

	//ARM_GP_Vector * cpnCoverages = soVanillaArg->GetPayPeriods();
	//ARM_Vector* flowPayDates = so->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates();
	//size_t cpnSize = flowPayDates->size();

	//std::vector<double>& fixValues = soVanillaArg->GetFixValues();
	//std::vector<double>& payMults   = soVanillaArg->GetPayIndexLeverages();

	//ARM_Vector* periodResetDates =  so->GetPayIndexLeg()->GetResetDates();
	//ARM_Vector* periodPayDates   = so->GetPayIndexLeg()->GetPaymentDates();
	//ARM_Vector* periodStartDates =  so->GetPayIndexLeg()->GetFlowStartDates();

	//
	//itsEndDate = periodPayDates->Elt(cpnSize-1);
	//itsRealEndDate = soVanillaArg->GetSwapLongFloatEndTime()->Elt(cpnSize-1) + asOfDate;
	//
	//itsResetDates.resize(cpnSize);
	//itsStartDates.resize(cpnSize);
	//itsPayDates.resize(cpnSize);
	//itsNotional.resize(cpnSize);
	//
	//itsFixValues.resize(cpnSize);
	//itsFixedCpnValue.resize(cpnSize);
	//itsFixedCpnLongSensi.resize(cpnSize);
	//itsFixedCpnShortSensi.resize(cpnSize);
	//
	//itsPayIndexMult.resize(cpnSize);
	//itsPayIndexFwd.resize(cpnSize);
	//itsVarCpnValue.resize(cpnSize);
	//itsVarCpnLongSensi.resize(cpnSize);
	//itsVarCpnShortSensi.resize(cpnSize);

	//itsCpnCoverages.resize(cpnSize);

	//itsLong0.resize(cpnSize);
	//itsShort0.resize(cpnSize);
	//itsPay0.resize(cpnSize);
	//itsDF.resize(cpnSize);
	//
	//const double NSTDEV_NO_IMPLIED_VOL = 6.0;
	//		
	//for (size_t i (0); i<cpnSize; i++)
	//{
	//	itsResetDates[i]	= (*periodResetDates)[i];
	//	itsStartDates[i]	= (*periodStartDates)[i];
	//	itsPayDates[i]		= (*periodPayDates)[i];
	//	itsNotional[i]		= so->GetPayIndexLeg()->GetAmount()->CptReferenceValue(itsPayDates[i]);

	//	itsCpnCoverages[i]  = (*cpnCoverages)[i];
	//	itsPayIndexMult[i]  = (*payMults)[i];
	//	itsFixValues[i]		= (*fixValues)[i];
	//	
	//	double payDate=(*flowPayDates)[i];
	//	int periodIndex			= periodPayDates->find(payDate);

	//	double shortFwd	  = 0.01 * so->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(i);
	//	double longFwd    = 0.01 * so->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(i);
	//	double correl	  = 0.01 * so->GetCorrelVector()->Elt(i);
	//	double coeffShort = so->GetWeight1();
	//	double coeffLong  = so->GetWeight2();
	//	double optMat     = (periodResetDates->Elt(i) - asOfDate)/K_YEAR_LEN;

	//	if (optMat < 0)
	//		continue;
	//	
	//	double strike = 0.01 * so->GetStrikes()->CptReferenceValue((*periodResetDates)[periodIndex]);
	//	double fixedRate  = so->GetPayIndexMargins()->CptReferenceValue((*periodResetDates)[periodIndex]);

	//	double spread1 = 0.01*so->GetSpread1();
	//	double spread2 = 0.01*so->GetSpread2();

	//	double longTenor  = (soVanillaArg->GetSwapLongFloatEndTime()->Elt(i) - soVanillaArg->GetSwapLongFloatStartTime()->Elt(i)) / K_YEAR_LEN;
	//	double shortTenor = (soVanillaArg->GetSwapShortFloatEndTime()->Elt(i) - soVanillaArg->GetSwapShortFloatStartTime()->Elt(i)) / K_YEAR_LEN;
	//				
	//	double atmShortVol  =  0.01 * so->GetVol1ATM(optMat, shortTenor);
	//	double atmLongVol   =  0.01 * so->GetVol2ATM(optMat, longTenor);
	//	atmShortVol	*= coeffShort * shortFwd ;
	//	atmLongVol	*= coeffLong  * longFwd;
	//	double atmSpreadVol =  sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
	//	double atmSpreadStdev = atmSpreadVol * sqrt(optMat);
	//	
	//	double forwardSpread = coeffLong * longFwd - coeffShort * shortFwd;

	//	double payDf		= curve->DiscountPrice((so->GetSwapLeg()->GetPaymentDates()->Elt(i)-asOfDate)/K_YEAR_LEN);
	//	double payNotional	= so->GetPayIndexLeg()->GetAmount()->CptReferenceValue(so->GetSwapLeg()->GetPaymentDates()->Elt(i));
	//	double payPeriod	= itsCpnCoverages[i];

	//	/// take vols ATM
	//	double VolLeft  = atmSpreadVol;
	//	double VolRight = atmSpreadVol;
	//	
	//	double VolLeftVar  = atmSpreadVol;
	//	double VolRightVar = atmSpreadVol;

	//	double t = periodStartDates->Elt(i);
	//	bool isVariable = (so->GetPayIndexLeg()->GetIndexType()!=K_FIXED);
	//	itsIsVariable = isVariable;

	//	double forwardSpreadVar (0.0);
	//	
	//	if (isVariable)
	//	{
	//		double shortFwdVar	= 0.01 * so->GetFirstFwdRateFLT()->Elt(i);
	//		double longFwdVar	= 0.01 * so->GetSecondFwdRateFLT()->Elt(i);
	//		forwardSpreadVar	= coeffLong * longFwdVar- coeffShort * shortFwdVar;
	//	}

	//	/// overwrite only if not too far from ATM
	//	if (fabs(strike - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
	//	{
	//		int callPut = so->IsCap()-so->IsFloor();
	//		double targetPriceLeft_N = so->GetCap1Prices()->Elt(i);
	//		VolLeft = VanillaImpliedVol_N (forwardSpread, targetPriceLeft_N, strike+spread1, optMat, callPut, &atmSpreadVol);
	//		double targetPriceRight_N = so->GetCap2Prices()->Elt(i);
	//		VolRight = VanillaImpliedVol_N (forwardSpread, targetPriceRight_N, strike+spread2, optMat, callPut, &atmSpreadVol);

	//		if (isVariable)
	//		{
	//			targetPriceLeft_N = so->GetCap1PricesFLT()->Elt(i);
	//			VolLeftVar = VanillaImpliedVol_N (forwardSpreadVar, targetPriceLeft_N,  strike+spread1, optMat, callPut, &atmSpreadVol);
	//			targetPriceRight_N = so->GetCap2PricesFLT()->Elt(i);
	//			VolRightVar = VanillaImpliedVol_N (forwardSpreadVar, targetPriceRight_N, strike+spread2, optMat, callPut, &atmSpreadVol);
	//		}
	//	}


	//	///---------------------------------
	//	/// initial prices
	//	///---------------------------------
	//	/// it is safer to recompute them instead of taking the price
	//	/// already computed in the portfolio
	//	double sensi;
	//	double price;
	//	
	//	price  = (VanillaOption_N (forwardSpread, VolLeft, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread, VolRight, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//	
	//	
	//	///---------------------------------
	//	/// sensi w.r.t short CMS rate
	//	///---------------------------------
	//	sensi  = (VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, VolLeft, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, VolRight, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//	
	//	sensi-=price;

	//	sensi /= SENSI_CMS;

	//	itsFixedCpnShortSensi[i] = sensi;


	//	///---------------------------------
	//	/// sensi w.r.t long CMS rate
	//	///---------------------------------
	//	sensi  = (VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, VolLeft, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, VolRight, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//	
	//	sensi-=price;

	//	sensi /= SENSI_CMS;

	//	itsFixedCpnLongSensi[i] = sensi;


	//	///---------------------------------
	//	/// Cpn value
	//	///---------------------------------
	//	itsFixedCpnValue[i] = price;
	//			
	//	if (isVariable)
	//	{
	//		price  = (VanillaOption_N (forwardSpreadVar, VolLeftVar, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVar, VolRightVar, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//		
	//		
	//		///---------------------------------
	//		/// sensi w.r.t short CMS rate
	//		///---------------------------------
	//		sensi  = (VanillaOption_N (forwardSpreadVar- coeffShort * SENSI_CMS, VolLeftVar, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVar - coeffShort * SENSI_CMS, VolRightVar, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//		
	//		sensi-=price;

	//		sensi /= SENSI_CMS;

	//		itsVarCpnShortSensi[i] = sensi;


	//		///---------------------------------
	//		/// sensi w.r.t long CMS rate
	//		///---------------------------------
	//		sensi  = (VanillaOption_N (forwardSpreadVar+ coeffLong * SENSI_CMS, VolLeftVar, strike+spread1, optMat, K_CALL)-VanillaOption_N (forwardSpreadVar + coeffLong * SENSI_CMS, VolRightVar, strike+spread2, optMat, K_CALL))/(spread2-spread1);
	//		
	//		sensi-=price;

	//		sensi /= SENSI_CMS;

	//		itsVarCpnLongSensi[i] = sensi;


	//		///---------------------------------
	//		/// Cpn value
	//		///---------------------------------
	//		itsVarCpnValue[i] = price;


	//		///---------------------------------
	//		/// Fwd of payment index
	//		///---------------------------------
	//		itsPayIndexFwd[i] = 0.01 * so->GetPayFwdRateFLT()->Elt(i);

	//	}
	//	else
	//	{
	//		itsVarCpnValue[i]      = 0.0;
	//		itsVarCpnLongSensi[i]  = 0.0;
	//		itsVarCpnShortSensi[i] = 0.0;
	//		itsPayIndexFwd[i]      = 0.0;
	//	}
	//}
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiCSO
///	Routine: GetCpn
///	Returns: double
///	Action : returns structured coupon
////////////////////////////////////////////////////
double ARM_BasketSensiCSO::GetCpn(int i)
{
	double cpn;
	cpn  = itsCpnCoverages[i] * itsFixedCpnValue[i] * itsFixValues[i]/ARM_Constants::rateBase;
	cpn += itsCpnCoverages[i] * itsVarCpnValue[i]   * itsPayIndexMult[i] * itsPayIndexFwd[i];
	return cpn;
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiCSO
///	Routine: GetDCpn
///	Returns: double
///	Action : computes sensi of structured coupon from sensi of building rates
////////////////////////////////////////////////////
double ARM_BasketSensiCSO::GetDCpn(int i, const ARM_STRIPPER& stripper)
{
	double dcpn,Blong,Bshort,Bpay;
	Blong	= itsCpnCoverages[i] * itsFixedCpnLongSensi[i]  * itsFixValues[i]/ARM_Constants::rateBase;
	Bshort	= itsCpnCoverages[i] * itsFixedCpnShortSensi[i] * itsFixValues[i]/ARM_Constants::rateBase;
	Blong	+= itsCpnCoverages[i] * itsVarCpnLongSensi[i]  * itsPayIndexMult[i] * itsPayIndexFwd[i];
	Bshort	+= itsCpnCoverages[i] * itsVarCpnShortSensi[i] * itsPayIndexMult[i] * itsPayIndexFwd[i];
	Bpay	= itsCpnCoverages[i] * itsVarCpnValue[i] * itsPayIndexMult[i];

	double longStartTime = itsArg->GetSwapLongFloatStartTime()->Elt(i);
	double longEndTime   = itsArg->GetSwapLongFloatEndTime()->Elt(i);
	ARM_GP_Vector* longPayTimes    = itsArg->GetSwapLongFixPayTimes()[i];
	ARM_GP_Vector* longPayPeriods  = itsArg->GetSwapLongFixPayPeriods()[i];
	double Long = stripper.swapRate(longStartTime, longEndTime, longPayTimes->GetValues(), longPayPeriods->GetValues());

	double shortStartTime = itsArg->GetSwapShortFloatStartTime()->Elt(i);
	double shortEndTime   = itsArg->GetSwapShortFloatEndTime()->Elt(i);
	ARM_GP_Vector* shortPayTimes    = itsArg->GetSwapShortFixPayTimes()[i];
	ARM_GP_Vector* shortPayPeriods  = itsArg->GetSwapShortFixPayPeriods()[i];
	double Short = stripper.swapRate(shortStartTime, shortEndTime, shortPayTimes->GetValues(), shortPayPeriods->GetValues());

	dcpn	=	Blong * (Long  - itsLong0[i]);
	dcpn	+=	Bshort * (Short - itsShort0[i]);
	return dcpn;
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiCSO
///	Routine: Init
///	Returns: void
///	Action : computes value of building rates
////////////////////////////////////////////////////
void ARM_BasketSensiCSO::Init(int i, const ARM_STRIPPER& stripper)
{
	double asof = stripper.asOf();
	itsDF[i] = stripper.df(itsPayDates[i]-asof);
	double longStartTime = itsArg->GetSwapLongFloatStartTime()->Elt(i);
	double longEndTime   = itsArg->GetSwapLongFloatEndTime()->Elt(i);
	ARM_GP_Vector* longPayTimes    = itsArg->GetSwapLongFixPayTimes()[i];
	ARM_GP_Vector* longPayPeriods  = itsArg->GetSwapLongFixPayPeriods()[i];
	itsLong0[i] = stripper.swapRate(longStartTime, longEndTime, longPayTimes->GetValues(), longPayPeriods->GetValues());

	double shortStartTime = itsArg->GetSwapShortFloatStartTime()->Elt(i);
	double shortEndTime   = itsArg->GetSwapShortFloatEndTime()->Elt(i);
	ARM_GP_Vector* shortPayTimes    = itsArg->GetSwapShortFixPayTimes()[i];
	ARM_GP_Vector* shortPayPeriods  = itsArg->GetSwapShortFixPayPeriods()[i];
	itsShort0[i] = stripper.swapRate(shortStartTime, shortEndTime, shortPayTimes->GetValues(), shortPayPeriods->GetValues());

	itsPay0[i] = 0.0;

	if (IsVariableCpn(i))
	{
		ARM_IntVector* periodIndexes = itsArg->GetPeriodIndex();
		if (i==0 || ( (i>0) && ((*periodIndexes)[i]!=(*periodIndexes)[i-1]) ) )
		{
			int periodIndex = (*periodIndexes)[i];
			double payStartTime = itsArg->GetSwapPayFloatStartTime()->Elt(periodIndex);
			double payEndTime   = itsArg->GetSwapPayFloatEndTime()->Elt(periodIndex);
			ARM_GP_Vector* payTimes    = itsArg->GetSwapPayFixPayTimes()[periodIndex];
			ARM_GP_Vector* payPeriods  = itsArg->GetSwapPayFixPayPeriods()[periodIndex];
			itsPay0[i] = stripper.swapRate(payStartTime, payEndTime, payTimes->GetValues(), payPeriods->GetValues());
		}
		else
			itsPay0[i] = 0.0;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFunding
///	Routine: constructor
///	Returns: void
///	Action : computes sensis of SO coupon
////////////////////////////////////////////////////
ARM_BasketSensiFunding::ARM_BasketSensiFunding(ARM_ZeroCurve* curve,ARM_SwapLeg* leg,double w, double p)
	:	ARM_BasketSensi(w,p)
{
	/*double asOfDate	  = curve->GetAsOfDateJul();
	ARM_Currency* ccy = curve->GetCurrencyUnit();


	ARM_Vector* periodResetDates =  leg->GetResetDates();
	ARM_Vector* periodStartDates =  leg->GetFlowStartDates();
	ARM_Vector* periodEndDates =  leg->GetFlowEndDates();
	ARM_ReferenceValue* fundingSpread =  leg->GetSpreads();
	ARM_Vector* cpncov =  leg->GetInterestTermValues();
	
	size_t cpnSize = periodResetDates->size();

	itsResetDates.resize(cpnSize);
	itsStartDates.resize(cpnSize);
	itsPayDates.resize(cpnSize);
	itsNotional.resize(cpnSize);

	itsFundingSpreads.resize(cpnSize);
	itsCpnCoverages.resize(cpnSize);
	itsFwdRates.resize(cpnSize);
	itsDF.resize(cpnSize);

	for (size_t i (0); i<cpnSize; i++)
	{
		itsResetDates[i] = (*periodResetDates)[i];
		itsStartDates[i] = (*periodStartDates)[i];
		itsPayDates[i] = (*periodEndDates)[i];
		itsNotional[i] = leg->GetAmount()->CptReferenceValue(itsResetDates.Elt(i));
		itsFundingSpreads[i] = (fundingSpread?fundingSpread->CptReferenceValue(itsResetDates.Elt(i)):leg->GetSpread());
		itsCpnCoverages[i] = (*cpncov)[i]*100/itsNotional[i];
	}

	itsEndDate = itsPayDates.Elt(cpnSize-1);
	itsRealEndDate = itsEndDate;*/

}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFunding
///	Routine: GetCpn
///	Returns: double
///	Action : returns structured coupon
////////////////////////////////////////////////////
double ARM_BasketSensiFunding::GetCpn(int i)
{
	double cpn=0.;
	cpn  = itsCpnCoverages[i] * (itsFwdRates[i] + itsFundingSpreads[i]/100.);
	return cpn;
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFunding
///	Routine: GetDCpn
///	Returns: double
///	Action : computes sensi of structured coupon from sensi of building rates
////////////////////////////////////////////////////
double ARM_BasketSensiFunding::GetDCpn(int i, const ARM_STRIPPER& stripper)
{
	double asof		= stripper.asOf();
	double newfwd	= (stripper.df(itsStartDates[i]-asof)/stripper.df(itsPayDates[i]-asof)-1.)/itsCpnCoverages[i];
	double dcpn		= itsCpnCoverages[i] * (newfwd - itsFwdRates[i]);
	return dcpn;
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFunding
///	Routine: Init
///	Returns: void
///	Action : computes value of building rates
////////////////////////////////////////////////////
void ARM_BasketSensiFunding::Init(int i, const ARM_STRIPPER& stripper)
{
	double asof		= stripper.asOf();
	itsDF[i]		= stripper.df(itsPayDates[i]-asof);
	itsFwdRates[i]	= (stripper.df(itsStartDates[i]-asof)/itsDF[i]-1.)/itsCpnCoverages[i];
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFunding
///	Routine: UpdateRompu
///	Returns: void
///	Action : computes stub
////////////////////////////////////////////////////
double ARM_BasketSensiFunding::ComputeFee(int i,double startDate)
{
	double fee;
	fee = itsNotional[i] * itsCpnCoverages[i] * (startDate-itsStartDates[i])/(itsPayDates[i]-itsStartDates[i]) * itsFwdRates[i]/100.;
	return fee;
}



////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFixLeg
///	Routine: constructor
///	Returns: void
///	Action : computes sensis of SO coupon
////////////////////////////////////////////////////
ARM_BasketSensiFixLeg::ARM_BasketSensiFixLeg(ARM_ZeroCurve* curve,ARM_SwapLeg* leg,double w, double p)
	:	ARM_BasketSensi(w,p)
{
	/*double asOfDate	  = curve->GetAsOfDateJul();
	ARM_Currency* ccy = curve->GetCurrencyUnit();


	ARM_Vector* periodResetDates =  leg->GetResetDates();
	ARM_Vector* periodStartDates =  leg->GetFlowStartDates();
	ARM_Vector* payDates =  leg->GetPaymentDates();
	ARM_Vector* cpncov =  leg->GetInterestTermValues();
	ARM_Vector* fixRates =  leg->GetFwdRates();
	
	size_t cpnSize = payDates->size();

	itsResetDates.resize(cpnSize);
	itsStartDates.resize(cpnSize);
	itsPayDates.resize(cpnSize);
	itsNotional.resize(cpnSize);

	itsCpnCoverages.resize(cpnSize);
	itsFixRates.resize(cpnSize);
	itsDF.resize(cpnSize);
	
	for (size_t i (0); i<cpnSize; i++)
	{
		itsResetDates[i] = (*periodResetDates)[i];
		itsStartDates[i] = (*periodStartDates)[i];
		itsPayDates[i] = (*payDates)[i];
		itsNotional[i] = leg->GetAmount()->CptReferenceValue(itsPayDates[i]);
		itsCpnCoverages[i] = (*cpncov)[i]*100/itsNotional[i];
		itsFixRates[i] = (*fixRates)[i];
	}

	itsEndDate = itsPayDates.Elt(cpnSize-1);
	itsRealEndDate = itsEndDate;*/

}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFixLeg
///	Routine: GetCpn
///	Returns: double
///	Action : returns structured coupon
////////////////////////////////////////////////////
double ARM_BasketSensiFixLeg::GetCpn(int i)
{
	double cpn=0.;
	cpn  = itsCpnCoverages[i] * itsFixRates[i]/100.;
	return cpn;
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFixLeg
///	Routine: Init
///	Returns: void
///	Action : computes value of building rates
////////////////////////////////////////////////////
void ARM_BasketSensiFixLeg::Init(int i, const ARM_STRIPPER& stripper)
{
	double asof		= stripper.asOf();
	itsDF[i]		= stripper.df(itsPayDates[i]-asof);
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiFixLeg
///	Routine: UpdateRompu
///	Returns: void
///	Action : computes stub
////////////////////////////////////////////////////
double ARM_BasketSensiFixLeg::ComputeFee(int i,double startDate)
{
	double fee;
	fee = itsNotional[i] * itsCpnCoverages[i] * (startDate-itsStartDates[i])/(itsPayDates[i]-itsStartDates[i]) * itsFixRates[i]/100.;
	return fee;
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiSO
///	Routine: constructor
///	Returns: void
///	Action : computes sensis of SO coupon
////////////////////////////////////////////////////
ARM_BasketSensiSO::ARM_BasketSensiSO(ARM_ZeroCurve* curve,ARM_SpreadOption* so,double w, double p)
	:	ARM_BasketSensi(w,p)
{
//	double asOfDate	  = curve->GetAsOfDateJul();
//	ARM_Currency* ccy = curve->GetCurrencyUnit();
//
//	ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)ARM_ConverterFromKernel::ConvertSecuritytoArgObject(so,asOfDate);
//	soVanillaArg->ComputeIndexSchedulesAndAdjustDates(ccy, asOfDate);
//	itsArg = soVanillaArg;
//
//	ARM_Vector* periodPayDates = so->GetSpreadLeg()->GetFirstLeg()->GetPaymentDates();
//	ARM_Vector* periodResetDates = so->GetSpreadLeg()->GetFirstLeg()->GetResetDates();
//	ARM_Vector* periodStartDates = so->GetSpreadLeg()->GetFirstLeg()->GetFlowStartDates();
//	ARM_GP_Vector* cpnCoverages = soVanillaArg->GetPayPeriods();
//
//	size_t cpnSize = periodPayDates->size();
//
//	itsEndDate = periodPayDates->Elt(cpnSize-1);
//	itsRealEndDate = soVanillaArg->GetSwapLongFloatEndTime()->Elt(cpnSize-1) + asOfDate;
//	
//	itsResetDates.resize(cpnSize);
//	itsStartDates.resize(cpnSize);
//	itsPayDates.resize(cpnSize);
//	itsNotional.resize(cpnSize);
//
//	itsCpnValue.resize(cpnSize);
//	itsCpnLongSensi.resize(cpnSize);
//	itsCpnShortSensi.resize(cpnSize);
//
//	itsCpnCoverages.resize(cpnSize);
//
//	itsLong0.resize(cpnSize);
//	itsShort0.resize(cpnSize);
//	itsDF.resize(cpnSize);
//	
////	size_t soIdx (0);
//	double fixRate;
//
//	const double NSTDEV_NO_IMPLIED_VOL = 6.0;
//		
//	for (size_t i (0); i<cpnSize; i++)
//	{
//		itsResetDates[i]	= (*periodResetDates)[i];
//		itsStartDates[i]	= (*periodStartDates)[i];
//		itsPayDates[i]		= (*periodPayDates)[i];
//		itsNotional[i]		= soVanillaArg->GetNotional()->Elt(i);
//
//		itsCpnCoverages[i]  = (*cpnCoverages)[i];
//		
//		double optMat     = (periodResetDates->Elt(i) - asOfDate)/K_YEAR_LEN;
//
//		if(optMat <0)
//			continue;
//	
//		if (1==0)
//		{
//			//cas fixe à traiter
//			itsCpnShortSensi[i] = 0.0;
//			itsCpnLongSensi[i]  = 0.0;
//			itsCpnValue[i]	    = fixRate;
//		}
//		else
//		{
////			if (soIdx>=sizeSO || ( (i==cpnSize-1) && (soIdx !=sizeSO-1) ) )
////				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOSensitivities : unresolved problem" );
//
//			double shortFwd	 = 0.01 * so->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(i);
//			double longFwd   = 0.01 * so->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(i);
//			double correl	 = 0.01 * so->GetCorrelVector()->Elt(i);
//			double coeffShort= so->GetWeight1();
//			double coeffLong = so->GetWeight2();			
//			
//			double strike = 0.01 * so->GetStrikes()->CptReferenceValue(itsResetDates[i]);
//			
//			double longTenor  = (soVanillaArg->GetSwapLongFloatEndTime()->Elt(i) - soVanillaArg->GetSwapLongFloatStartTime()->Elt(i)) / K_YEAR_LEN;
//			double shortTenor = (soVanillaArg->GetSwapShortFloatEndTime()->Elt(i) - soVanillaArg->GetSwapShortFloatStartTime()->Elt(i)) / K_YEAR_LEN;
//
//			double atmShortVol  =  0.01 * so->GetVol1ATM(optMat, shortTenor);
//			double atmLongVol   =  0.01 * so->GetVol2ATM(optMat, longTenor);
//			atmShortVol	*= coeffShort * shortFwd ;
//			atmLongVol	*= coeffLong  * longFwd;
//			double atmSpreadVol =  sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
//			double atmSpreadStdev = atmSpreadVol * sqrt(optMat);
//			
//			double forwardSpread = coeffLong * longFwd - coeffShort * shortFwd;
//
//			double payDf		= curve->DiscountPrice(soVanillaArg->GetPayTimes()->Elt(i)/K_YEAR_LEN);
//			double payNotional	= itsNotional[i];
//			double payPeriod	= itsCpnCoverages[i];
//
//			/// take vols ATM
//			double Vol = atmSpreadVol;
//						
//			/// overwrite only if not too far from ATM
//			if (fabs(strike - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
//			{
//				int callPut = so->IsCap()-so->IsFloor();
//				double targetPrice_N = so->GetCashFlowValues()->Elt(i);
//				targetPrice_N /= (payDf * payNotional * payPeriod);
//				Vol = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strike, optMat, callPut, &atmSpreadVol);
//			}
//
//			///---------------------------------
//			/// initial prices
//			///---------------------------------
//			/// it is safer to recompute them instead of taking the price
//			/// already computed in the portfolio
//			double sensi;
//			double price0;
//			
//			price0 = VanillaOption_N (forwardSpread, Vol, strike, optMat, K_CALL);
//			
//			
//			///---------------------------------
//			/// sensi w.r.t short CMS rate
//			///---------------------------------
//			sensi  = VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, Vol, strike, optMat, K_CALL);
//			sensi -= price0;
//
//			sensi /= SENSI_CMS;
//
//			itsCpnShortSensi[i] = sensi;
//
//
//			///---------------------------------
//			/// sensi w.r.t long CMS rate
//			///---------------------------------
//			sensi  = VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, Vol, strike, optMat, K_CALL);
//			sensi -= price0;
//
//			sensi /= SENSI_CMS;
//
//			itsCpnLongSensi[i] = sensi;
//
//
//			///---------------------------------
//			/// Cpn value
//			///---------------------------------
//			itsCpnValue[i] = price0;
//
//			///---------------------------------
//			/// don't forget to soIdx ++
//			///---------------------------------
////			soIdx ++;
//		}
//	}
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiSO
///	Routine: GetCpn
///	Returns: double
///	Action : returns structured coupon
////////////////////////////////////////////////////
double ARM_BasketSensiSO::GetCpn(int i)
{
	double cpn;
	cpn  = itsCpnCoverages[i] * itsCpnValue[i];
	return cpn;
}


////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiSO
///	Routine: GetDCpn
///	Returns: double
///	Action : computes sensi of structured coupon from sensi of building rates
////////////////////////////////////////////////////
double ARM_BasketSensiSO::GetDCpn(int i, const ARM_STRIPPER& stripper)
{
	double dcpn,Blong,Bshort;
	Blong	= itsCpnCoverages[i] * itsCpnLongSensi[i];
	Bshort	= itsCpnCoverages[i] * itsCpnShortSensi[i];
	
	double longStartTime = itsArg->GetSwapLongFloatStartTime()->Elt(i);
	double longEndTime   = itsArg->GetSwapLongFloatEndTime()->Elt(i);
	ARM_GP_Vector* longPayTimes    = itsArg->GetSwapLongFixPayTimes()[i];
	ARM_GP_Vector* longPayPeriods  = itsArg->GetSwapLongFixPayPeriods()[i];
	double Long = stripper.swapRate(longStartTime, longEndTime, longPayTimes->GetValues(), longPayPeriods->GetValues());

	double shortStartTime = itsArg->GetSwapShortFloatStartTime()->Elt(i);
	double shortEndTime   = itsArg->GetSwapShortFloatEndTime()->Elt(i);
	ARM_GP_Vector* shortPayTimes    = itsArg->GetSwapShortFixPayTimes()[i];
	ARM_GP_Vector* shortPayPeriods  = itsArg->GetSwapShortFixPayPeriods()[i];
	double Short = stripper.swapRate(shortStartTime, shortEndTime, shortPayTimes->GetValues(), shortPayPeriods->GetValues());

	dcpn	=	Blong * (Long  - itsLong0[i]);
	dcpn	+=	Bshort * (Short - itsShort0[i]);
	return dcpn;
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketSensiSO
///	Routine: Init
///	Returns: void
///	Action : computes value of building rates
////////////////////////////////////////////////////
void ARM_BasketSensiSO::Init(int i, const ARM_STRIPPER& stripper)
{
	double asof = stripper.asOf();
	itsDF[i] = stripper.df(itsPayDates[i]-asof);
	double longStartTime = itsArg->GetSwapLongFloatStartTime()->Elt(i);
	double longEndTime   = itsArg->GetSwapLongFloatEndTime()->Elt(i);
	ARM_GP_Vector* longPayTimes    = itsArg->GetSwapLongFixPayTimes()[i];
	ARM_GP_Vector* longPayPeriods  = itsArg->GetSwapLongFixPayPeriods()[i];
	itsLong0[i] = stripper.swapRate(longStartTime, longEndTime, longPayTimes->GetValues(), longPayPeriods->GetValues());

	double shortStartTime = itsArg->GetSwapShortFloatStartTime()->Elt(i);
	double shortEndTime   = itsArg->GetSwapShortFloatEndTime()->Elt(i);
	ARM_GP_Vector* shortPayTimes    = itsArg->GetSwapShortFixPayTimes()[i];
	ARM_GP_Vector* shortPayPeriods  = itsArg->GetSwapShortFixPayPeriods()[i];
	itsShort0[i] = stripper.swapRate(shortStartTime, shortEndTime, shortPayTimes->GetValues(), shortPayPeriods->GetValues());
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

