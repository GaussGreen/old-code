// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 8/11/99 Afshin Bayrooti
//

#include "kplatform.h"
#include "Magnet/Magnet.h"

#include "smileaddin.h"
#include "genOption.h"
#include "kratecurve.h"
#include "defswap.h"
#include "fstream.h"
#include "cashflow.h"
#include "kdate.h"
#include "ddmap.h"

#include "capital.h"


#if defined(_WIN32)
    #pragma warning( disable: 4190 )
#endif

using namespace CM;

//
//   F i r s t A d d i n
//

MAGNET_FUNCTION0( String, Version )
{
	return HY_version();
   // return "HY V2.4";
}

//
//   A d d
//

MAGNET_FUNCTION2( double, Add, double a, double b )
{
    return a + b;
}


//
//   F r a c t i o n
//

// How to read an ALIB data type.

MAGNET_FUNCTION8( 
	double, optionsmile,
		String optionTypeC,
		double fwd,
		double strike,
		double yearsToExpiry,
		double yearsToPayment,
		double volatility, 
		double rate, 
		BaseFunction* inputProcess)
{
BaseOption* option;

	if(CompareNoCase( optionTypeC, "C" ) == 0)
	{
		option = new CallOption(strike,yearsToExpiry, yearsToPayment);
	}
	else
	{
		option = new PutOption(strike,yearsToExpiry, yearsToPayment);
	}
    
    return GenOptionPricer(fwd,
						   option,
						   volatility,
						   rate,
						   inputProcess);
}

MAGNET_FUNCTION2(double, Inverse, double y, CM::SharedPointer<BaseFunction> funcPtr)
{
	return funcPtr->inverse(y);
}

MAGNET_FUNCTION2(double, Derivative, double x, CM::SharedPointer<BaseFunction> funcPtr)
{
	return funcPtr->deriv(x);
}

MAGNET_FUNCTION0(CM::SharedPointer<LogNormalFunction>, Log_Normal)
{
	return CM::Raw2SharedPointer( new LogNormalFunction() );
}

MAGNET_FUNCTION0(CM::SharedPointer<NormalLogNormalFunction>, Normal_Log_Normal)
{
	return CM::Raw2SharedPointer( new NormalLogNormalFunction() );
}


MAGNET_FUNCTION0(CM::SharedPointer<NormalLogNormalFunction2>, Normal_Log_Normal2)
{
	return CM::Raw2SharedPointer( new NormalLogNormalFunction2() );
}



MAGNET_FUNCTION3(CM::SharedPointer<DualLogNormalFunction>, Dual_Log_Normal, double turnPoint, double factor, double a)
{
	return CM::Raw2SharedPointer( new DualLogNormalFunction(turnPoint, factor,a) );
}


MAGNET_FUNCTION5(CM::SharedPointer<AnnuityIndex>, Annuity_Index, double maturity, double coupon, double turnPoint, double factor, double a)
{
	return CM::Raw2SharedPointer(new AnnuityIndex(maturity, coupon,turnPoint, factor,a) );
}



MAGNET_FUNCTION3(CM::SharedPointer<CallOption>, Call_Option, double strike,double yearsToExpiry, double yearsToPayment)
{
	return CM::Raw2SharedPointer( new CallOption(strike,yearsToExpiry,yearsToPayment) );
}

MAGNET_FUNCTION3(CM::SharedPointer<PutOption>, Put_Option, double strike,double yearsToExpiry, double yearsToPayment)
{
	return CM::Raw2SharedPointer( new PutOption(strike,yearsToExpiry,yearsToPayment) );
}


MAGNET_FUNCTION3(CM::SharedPointer<AssetToStockMapping1>, Asset_To_Stock, double x1,double level, double debtPerShare)
{
	AssetToStockMapping1* atsm;
	try
	{
		atsm = new AssetToStockMapping1(x1,level,debtPerShare);
	}

	catch (KException& e){ throw RuntimeError(e.what());}
	
	return CM::Raw2SharedPointer( atsm );
}

MAGNET_FUNCTION3(double, AssetVol_To_StockVol, double stock, double assetVol, CM::SharedPointer<AssetToStockMapping1> assetToStock)
{
	double  stockVol = assetToStock->assetVolToStockVol(stock, assetVol);
	return stockVol;
}

MAGNET_FUNCTION3(double, StockVol_To_AssetVol, double stock,double stockVol, CM::SharedPointer<AssetToStockMapping1> assetToStock)
{
	double  assetVol = assetToStock->stockVolToAssetVol(stock, stockVol);
	return assetVol;
}



MAGNET_FUNCTION0(CM::SharedPointer<DummyMapping>, Dummy)
{
	return CM::Raw2SharedPointer(new DummyMapping() );
}



MAGNET_FUNCTION1(CM::SharedPointer<CashFlowList>, CashFlow_List, String bondSpec)
{
	return CM::Raw2SharedPointer(new CashFlowList(bondSpec) );
}

MAGNET_DELETE_RESULT
MAGNET_FUNCTION3(DateMap<double>*, Get_CashFlow, Date currDate, CM::SharedPointer<CashFlowList> cf, 
				 CM::SharedPointer<KRateCurve> zc)
{
	DateMap<double>* cashflow = new DateMap<double>;
	CashFlowList::iterator Iter;
	TDate tdate = currDate;
	KDate kdate = tdate;
	double couponPay;
	for (Iter = cf->lower_bound(kdate); Iter != cf->end(); Iter++)
	{
		couponPay =Iter->second.get_cashflow(zc.get());
		TDate temp = Iter->first;
		Date  temp2 = TDate2Date(temp);
		cashflow->insert(temp2,couponPay);
	}
	return cashflow;
}

MAGNET_FUNCTION4(double, Get_Payment, Date date, CM::SharedPointer<KRateCurve> zc, String interpTypeC,
				 CM::SharedPointer<CashFlowList> cfl)
{
	TDate tdate = date;
	KDate kdate = tdate;
	long interpType;
	char* ptr = &interpTypeC[0];

	GtoHYInterpTypeCToI(ptr, &interpType);

	double payment = cfl->get_Payment(kdate,zc.get(),interpType);
	return payment;

}

MAGNET_X_FUNCTION2(CM::SharedPointer<DDMap>, DD_Map, String spec,MAGNET_MANDATORY, Date endDate, = Date() )
{
	TDate tdate;
	KDate kdate;
	if(endDate == Date())
	{
		kdate =0;
	}
	else
	{
		tdate = endDate;
		kdate = tdate;
	}
	DDMap* ddmap;
	try{
	 ddmap = new DDMap(spec,kdate);
	}
	catch(KException &e)
	{
		std::ofstream	err("C:\\temp\\error.log");
		err<<e.what();
		throw RuntimeError(e.what());
	}

	return CM::Raw2SharedPointer( ddmap );
}

MAGNET_DELETE_RESULT
MAGNET_FUNCTION1(DateMap<double>*, Get_DDMap, CM::SharedPointer<DDMap> ddmap)
{
	DateMap<double>* dMap = new DateMap<double>;

	DDMap::iterator Iter;
	for (Iter = ddmap->begin(); Iter != ddmap->end(); Iter++)
	{
		TDate temp = Iter->first;
		Date  temp2 = TDate2Date(temp);
		dMap->insert(temp2,Iter->second);
	}
	return dMap;
}


 MAGNET_FUNCTION3(double, Get_Ave_Amount, Date date1,Date date2, CM::SharedPointer<DDMap> ddmap)
{
	TDate tdate1 = date1, tdate2 = date2;
	KDate kdate1 = tdate1, kdate2 = tdate2;

	double amount = ddmap->get_average_amount(kdate1, kdate2);
	return amount;
}


MAGNET_X_FUNCTION2(CM::SharedPointer<CallOptionT>, Call_OptionT, 
				   CM::SharedPointer<KRateCurve> strike,MAGNET_MANDATORY, String isamerican, = "E")
{
	return CM::Raw2SharedPointer(new CallOptionT(strike.get(),isamerican) );
}

MAGNET_X_FUNCTION2(CM::SharedPointer<PutOptionT>, PUT_OptionT, CM::SharedPointer<KRateCurve> strike, MAGNET_MANDATORY,String isamerican, = "E")
{
	return CM::Raw2SharedPointer( new PutOptionT(strike.get(),isamerican)) ;
}

MAGNET_X_FUNCTION5(CM::SharedPointer<DefaultSwap>, Default_Protection, 
				   CM::SharedPointer<KRateCurve> strike1,MAGNET_MANDATORY, 
				   CM::SharedPointer<KRateCurve> strike2,MAGNET_MANDATORY,
				   double recovery,MAGNET_MANDATORY,
				   String feeSpec, ="",
				   CM::SharedPointer<OptionContext> option, = CM::SharedPointer<OptionContext>() )
{
	return CM::Raw2SharedPointer ( new DefaultSwap( strike1.get(), strike2.get(), recovery,feeSpec,option.get()) );
}

MAGNET_FUNCTION1(CM::SharedPointer<StrikesClass>, Strikes, String strikeSpec)
{
	return CM::Raw2SharedPointer (new StrikesClass(strikeSpec) );
}

MAGNET_X_FUNCTION3(CM::SharedPointer<OptionContext>, Option_Context, 
				   String optionSpec,MAGNET_MANDATORY, 
				   CM::SharedPointer<const StrikesClass> strikes,MAGNET_MANDATORY,
				   CM::SharedPointer<const StrikesClass> soft_strikes, = CM::SharedPointer<const StrikesClass>() )
{
	return CM::Raw2SharedPointer ( new OptionContext(optionSpec, strikes.get(),soft_strikes.get()) );
}

MAGNET_X_FUNCTION7(CM::SharedPointer<Bond>, Bond_Note, 
				   CM::SharedPointer<KRateCurve> strike1,MAGNET_MANDATORY, 
				   CM::SharedPointer<KRateCurve> strike2,MAGNET_MANDATORY,
				   double recovery,MAGNET_MANDATORY,
				   String bondSpec,MAGNET_MANDATORY,
				   double resetRate, = 0.0,
				   CM::SharedPointer<OptionContext> option, = CM::SharedPointer<OptionContext>(),
				   CM::SharedPointer<OptionContext> cvOption, = CM::SharedPointer<OptionContext>())
{
	return CM::Raw2SharedPointer (new Bond( strike1.get(), strike2.get(), recovery,bondSpec, resetRate, 
			option.get(),cvOption.get()) );
}


MAGNET_FUNCTION5(CM::SharedPointer<KRateCurve>,  Rate_Curve, 
				 Date baseDate, 
				 Array<Date>& dates, 
				 Array<double>& rates, 
				 double basis, 
				 long dayCount)
{
	TDate  valueDate = baseDate;
	if(dates.size() != rates.size())
	{
		KException("dates size and rates size are not equal.");
		return CM::SharedPointer<KRateCurve>();  //NULL
	}
	else
	{
		TDate* t = new TDate[dates.size()];
		double* r= new double[ rates.size()];

		for( size_t i = 0; i < dates.size(); ++i )
		{
			t[i] = dates[i];
			r[i] = rates[i];
		}
		return CM::Raw2SharedPointer( new KRateCurve(valueDate,t,r,rates.size(),basis,dayCount) );
	}
}


MAGNET_DELETE_RESULT

MAGNET_FUNCTION14(Array<double>*,  Gen_Pricer, double spotStockPrice,CM::SharedPointer<DDMap> dividentYield, 
				 CM::SharedPointer<KRateCurve> repoCurve,CM::SharedPointer<KRateCurve> irCurve,
				 CM::SharedPointer<KRateCurve> volCurve, CM::SharedPointer<KRateCurve> volCurveShift,
				 CM::SharedPointer<Instrument> instrument, 
				 CM::SharedPointer<BaseFunction> assetToStockMapping, 
				 CM::SharedPointer<BaseFunction> assetProcess, double ppy, double beta, 
				 double vollim, Date settleDate, bool isCVOption)
{
	Array<double>* result = new Array<double>;
	TDate tSettleDate = settleDate;
	KDate kSettleDate = tSettleDate;
	
	try
	{
		KValarray<double> p = GeneralPricer(spotStockPrice,
											   dividentYield.get(),
											   repoCurve.get(),
											   irCurve.get(),
											   volCurve.get(),
											   volCurveShift.get(),			//HY3.4v
											   instrument.get(),
											   assetToStockMapping.get(),
											   assetProcess.get(),
											   ppy,
											   beta,
											   vollim,						//HY4.1v
											   kSettleDate,
											   isCVOption);					//HY3.4v
		for(int i =0; i<p.size(); i++)
		{
			result->push_back( p[i]);
		}
	}
	catch (KException& e){
		std::ofstream	err("C:\\temp\\error.log");
		err<<e.what();
	};

	return result;
}

MAGNET_DELETE_RESULT

MAGNET_FUNCTION14(Array<double>*,  Gen_Pricer2, double spotStockPrice,CM::SharedPointer<DDMap> dividentYield, 
				 CM::SharedPointer<KRateCurve> repoCurve,CM::SharedPointer<KRateCurve> irCurve,
				 CM::SharedPointer<KRateCurve> volCurve, CM::SharedPointer<KRateCurve> volCurveShift,
				 CM::SharedPointer<Instrument> instrument, 
				 CM::SharedPointer<BaseFunction> assetToStockMapping, 
				 CM::SharedPointer<BaseFunction> assetProcess, double ppy, double beta,
				 double vll, Date settleDate, bool isCVOption)
{
	Array<double>* result = new Array<double>;
	TDate tSettleDate = settleDate;
	KDate kSettleDate = tSettleDate;
	int i;
	
	double spotPrice[2];
	spotPrice[0] = 1;
	spotPrice[1] = spotStockPrice;

	double divRefSpot[2];
	divRefSpot[0] = 1;
	divRefSpot[1] = spotStockPrice;

	double volRefSpot[2];
	volRefSpot[0] = 1;
	volRefSpot[1] = spotStockPrice;

	TDate   *divDates = NULL;
	double *divYields = NULL;
	int    divSize = dividentYield->size();
	if(divSize == 0)
	{
		divSize = 1;
		divDates = new TDate[1+1];
		divYields = new double[1+1];
		divDates[1] = tSettleDate;
		divYields[1] = 0;
	}
	else
	{
		divDates = new TDate[divSize+1];
		divYields = new double[divSize+1];
		
		DDMap::const_iterator iter = dividentYield->begin();
		for(i=1; i<= divSize; i++)
		{
			divDates[i] = iter->first;
			divYields[i] = iter->second;
			if (divRefSpot[1] != 0.0)
			{
				divYields[i] = divYields[i]*divRefSpot[1];
			}
			iter++;
		}
	}
	divDates[0] = divSize;
	divYields[0] = divSize;

	double dps[2];
	dps[0]=1;
	AssetToStockMapping1  *atsPt = dynamic_cast<AssetToStockMapping1*>(assetToStockMapping.get());
	dps[1]=atsPt->get_debtPerShare();

	double x[2];
	x[0] = 1;
	x[1] = atsPt->get_x1();

	double lim[2];
	lim[0] = 1;
	lim[1] = atsPt->get_lim();

	int repoSize = repoCurve->size();
	TDate  *repoDates = new TDate[repoSize+1];
	double  *repoRates = new double[repoSize+1];
	repoDates[0] = repoSize;
	repoRates[0] = repoSize;
	for(i=1; i<= repoSize; i++)
	{
		repoDates[i] = repoCurve->iDate(i-1);
		repoRates[i] = repoCurve->iRate(i-1);
	}

	int swapSize = irCurve->size();
	TDate  *swapDates = new TDate[swapSize+1];
	double  *swapRates = new double[swapSize+1];
	swapDates[0] = swapSize;
	swapRates[0] = swapSize;
	for(i=1; i<= swapSize; i++)
	{
		swapDates[i] = irCurve->iDate(i-1);
		swapRates[i] = irCurve->iRate(i-1);
	}

	int volSize = volCurve->size();
	TDate  *volDates = new TDate[volSize+1];
	double  *volRates = new double[volSize+1];
	volDates[0] = volSize;
	volRates[0] = volSize;
	for(i=1; i<= volSize; i++)
	{
		volDates[i] = volCurve->iDate(i-1);
		volRates[i] = volCurve->iRate(i-1);
		
		if(volRefSpot[1] != 0.0)
			volRates[i] = atsPt->assetVolToStockVol(volRefSpot[1],volRates[i]);
	}
	
	double volShift[2];
	volShift[0] = 1;
	volShift[1] = 0;

	long instType[2];
	instType[0] =1;
	double notional[2];
	notional[0] = 1;
	notional[1] =1.0;
	double recoveryRate[2];
	recoveryRate[0] = 1;
	TDate  *instCFAccStart = NULL;
	TDate  *instCFAccEnd = NULL;
	TDate  *instCFDates = NULL;
	double *instCFCoupons = NULL;
	double *instCFAmorts = NULL;
	TDate  *claimCFDates = NULL;
	double *claimCFAmounts = NULL;
	TDate *exerStartDates=NULL;
	TDate *exerEndDates=NULL;
	double *exerStartStrikes=NULL;
	double *exerEndStrikes=NULL;
	long *optionDirections=NULL;
	long *exerTypes = NULL;
	double  *option_barrier = NULL;
	long  *optionType = NULL;

	
	if(DefaultSwap *pt1 = dynamic_cast<DefaultSwap*>(instrument.get()))
	{
		
		recoveryRate[1] = pt1->get_recoveryRate();
		CashFlowList cfList = pt1->get_cashflowlist();
		CashFlowList::const_iterator cfIter = cfList.begin();
		int cfSize = cfList.size();
		instCFAccStart = new TDate[cfSize+1];
		instCFAccStart[0] = cfSize;
		instCFAccEnd = new TDate[cfSize+1];
		instCFAccEnd[0] = cfSize;
		instCFDates = new TDate[cfSize+1];
		instCFDates[0] = cfSize;
		instCFCoupons = new double[cfSize+1];
		instCFCoupons[0] = cfSize;
		instCFAmorts = new double[cfSize+1];
		instCFAmorts[0] = cfSize;
		claimCFDates = new TDate[3];
		claimCFDates[0] = 2;
		claimCFAmounts = new double[3];
		claimCFAmounts[0] = 2;
		
		for(i=1; i<= cfSize; i++)
		{
			instCFAccStart[i] = cfIter->second.get_accruStartDate();
			instCFAccEnd[i] = cfIter->second.get_accruEndDate();
			instCFDates[i] = cfIter->first;
			instCFCoupons[i] = -cfIter->second.get_couponPay(irCurve.get());
			instCFAmorts[i] = cfIter->second.get_amortPay();
			cfIter++;
		}
		claimCFDates[1] = instCFAccStart[1];
		claimCFDates[2] = instCFAccEnd[cfSize];
		claimCFAmounts[1] = 1.0;
		claimCFAmounts[2] = 1.0;
		OptionContext   *optionPt = pt1->get_option_pt();
		if(optionPt == NULL)
		{
			instType[1] = CDS;
			exerStartDates = new TDate[1+1];
			exerStartDates[0] = 1;
			exerStartDates[1] = tSettleDate;
			exerEndDates = new TDate[1+1];
			exerEndDates[0] = 1;
			exerEndDates[1] = tSettleDate;
			exerStartStrikes = new double[1+1];
			exerStartStrikes[0] = 1;
			exerStartStrikes[1] = 0;
			exerEndStrikes = new double[1+1];
			exerEndStrikes[0] = 1;
			exerEndStrikes[1] = 0;
			option_barrier = new double[1+1];
			option_barrier[0] = 1;
			option_barrier[1] = 0;
			optionType = new long[1+1];
			optionType[0] = 1;
			optionType[1] = HY_CAPITAL_CALL;
			optionDirections = new long[1+1];
			optionDirections[0] = 1;
			optionDirections[1] = 0;
			exerTypes = new long[1+1];
			exerTypes[0] = 1;
			exerTypes[1] = 0;

			
		}
		else
		{
			int strikeSize = optionPt->strikesSize();
			exerStartDates = new TDate[strikeSize+1];
			exerStartDates[0] = strikeSize;
			exerEndDates = new TDate[strikeSize+1];
			exerEndDates[0] = strikeSize;
			exerStartStrikes = new double[strikeSize+1];
			exerStartStrikes[0] = strikeSize;
			exerEndStrikes = new double[strikeSize+1];
			exerEndStrikes[0] = strikeSize;
			optionDirections = new long[strikeSize+1];
			optionDirections[0] = strikeSize;
			exerTypes = new long[strikeSize+1];
			exerTypes[0] = strikeSize;
			option_barrier = new double[strikeSize+1];
			option_barrier[0] = strikeSize;
			optionType = new long[strikeSize+1];
			optionType[0] = strikeSize;
			
			long optDirection = optionPt->get_option_direction();
			for(i=1; i<=strikeSize; i++)
			{
				exerStartDates[i] = optionPt->iDate(i-1);
				exerEndDates[i] = optionPt->iDate(i-1);
				exerStartStrikes[i] = optionPt->iRate(i-1);
				exerEndStrikes[i] = optionPt->iRate(i-1);
				optionDirections[i] = optDirection;
				if(optDirection == -1)
				{
					optionType[i] = HY_CAPITAL_CALL;
				}
				else
				{
					optionType[i] = HY_CAPITAL_PUT;
				}
				exerTypes[i] = optionPt->get_option_isamerican();
				option_barrier[i] = 0;
			}

			instType[1] = CDSOPT;
		
			
			
		}


	}
	else if(Bond *pt2 = dynamic_cast<Bond*>(instrument.get()))
	{
		recoveryRate[1] = pt2->get_recoveryRate();
		CashFlowList cfList = pt2->get_cashflowlist();
		CashFlowList::const_iterator cfIter = cfList.begin();
		int cfSize = cfList.size();
		instCFAccStart = new TDate[cfSize+1];
		instCFAccStart[0] = cfSize;
		instCFAccEnd = new TDate[cfSize+1];
		instCFAccEnd[0] = cfSize;
		instCFDates = new TDate[cfSize+1];
		instCFDates[0] = cfSize;
		instCFCoupons = new double[cfSize+1];
		instCFCoupons[0] = cfSize;
		instCFAmorts = new double[cfSize+1];
		instCFAmorts[0] = cfSize;
		claimCFDates = new TDate[3];
		claimCFDates[0] = 2;
		claimCFAmounts = new double[3];
		claimCFAmounts[0] = 2;
		
		for(i=1; i<= cfSize; i++)
		{
			instCFAccStart[i] = cfIter->second.get_accruStartDate();
			instCFAccEnd[i] = cfIter->second.get_accruEndDate();
			instCFDates[i] = cfIter->first;
			instCFCoupons[i] = cfIter->second.get_couponPay(irCurve.get());
			instCFAmorts[i] = cfIter->second.get_amortPay();
			cfIter++;
		}
		claimCFDates[1] = instCFAccStart[1];
		claimCFDates[2] = cfList.get_claimParDate();
		claimCFAmounts[1] = cfList.get_issuePrice();
		claimCFAmounts[2] = 1.0;
		OptionContext   *optionPt = pt2->get_option_pt();
		OptionContext   *optionCVPt = pt2->get_option_cv_pt();

		if(optionPt == NULL && optionCVPt == NULL)
		{
			instType[1] = BOND;
			
			exerStartDates = new TDate[1+1];
			exerStartDates[0] = 1;
			exerStartDates[1] = tSettleDate;
			exerEndDates = new TDate[1+1];
			exerEndDates[0] = 1;
			exerEndDates[1] = tSettleDate;
			exerStartStrikes = new double[1+1];
			exerStartStrikes[0] = 1;
			exerStartStrikes[1] = 0;
			exerEndStrikes = new double[1+1];
			exerEndStrikes[0] = 1;
			exerEndStrikes[1] = 0;
			option_barrier = new double[1+1];
			option_barrier[0] = 1;
			option_barrier[1] = 0;
			optionType = new long[1+1];
			optionType[0] = 1;
			optionType[1] = HY_CAPITAL_CALL;
			optionDirections = new long[1+1];
			optionDirections[0] = 1;
			optionDirections[1] = 0;
			exerTypes = new long[1+1];
			exerTypes[0] = 1;
			exerTypes[1] = 0;
			
		}
		else
		{
			int optionStrikeSize = (optionPt == NULL)? 0:optionPt->strikesSize();
			int softStrikeSize = (optionPt == NULL)? 0:optionPt->softStrikesSize();
			int convertStrikeSize = (optionCVPt == NULL)? 0:optionCVPt->strikesSize();

			int strikeSize = optionStrikeSize+softStrikeSize+convertStrikeSize;
			exerStartDates = new TDate[strikeSize+1];
			exerStartDates[0] = strikeSize;
			exerEndDates = new TDate[strikeSize+1];
			exerEndDates[0] = strikeSize;
			exerStartStrikes = new double[strikeSize+1];
			exerStartStrikes[0] = strikeSize;
			exerEndStrikes = new double[strikeSize+1];
			exerEndStrikes[0] = strikeSize;
			optionDirections = new long[strikeSize+1];
			optionDirections[0] = strikeSize;
			exerTypes = new long[strikeSize+1];
			exerTypes[0] = strikeSize;
			option_barrier = new double[strikeSize+1];
			option_barrier[0] = strikeSize;
			optionType = new long[strikeSize+1];
			optionType[0] = strikeSize;
			
			long optDirection = (optionPt==NULL)? 0:optionPt->get_option_direction();
			long optCVDirection = (optionCVPt==NULL)? 0:optionCVPt->get_option_direction();
			for(i=1; i<=strikeSize; i++)
			{
				if(i<=optionStrikeSize)
				{
					exerStartDates[i] = optionPt->iDate(i-1);
					exerEndDates[i] = optionPt->iDate(i-1);
					exerStartStrikes[i] = optionPt->iRate(i-1);
					exerEndStrikes[i] = optionPt->iRate(i-1);
					optionDirections[i] = optDirection;
					
					optionType[i] = HY_CAPITAL_CLEAN;
					exerTypes[i] = optionPt->get_option_isamerican();				
				}
				else if(i<=optionStrikeSize+softStrikeSize)
				{
					exerStartDates[i] = optionPt->iSDate(i-1-optionStrikeSize);
					exerEndDates[i] = exerStartDates[i];
					exerStartStrikes[i] = optionPt->iSRate(i-1-optionStrikeSize);
					exerEndStrikes[i] = exerStartStrikes[i];
					optionDirections[i] = optDirection;
					
					optionType[i] = HY_CAPITAL_SOFTCALL;
					exerTypes[i] = optionPt->get_option_isamerican();
				}
				else
				{
					if(optionCVPt!=NULL)
					{
						exerStartDates[i] = optionCVPt->iDate(i-1-optionStrikeSize-softStrikeSize);
						exerEndDates[i] = exerStartDates[i];
						exerStartStrikes[i] = optionCVPt->iRate(i-1-optionStrikeSize-softStrikeSize);
						exerEndStrikes[i] = exerStartStrikes[i];
						optionDirections[i] = optCVDirection;
						
						optionType[i] = HY_CAPITAL_SHARES;
						exerTypes[i] = optionCVPt->get_option_isamerican();
					}
				}

				
				option_barrier[i] = 0;
			}

			instType[1] = BONDOPT;
			
		}


	}
	else
	{
		instType[1] = EQUITYOPT;
		recoveryRate[1] = 0;
		
		instCFAccStart = new TDate[2];
		instCFAccStart[0] = 1;
		instCFAccStart[1] = settleDate;
		instCFAccEnd = new TDate[2];
		instCFAccEnd[0] = 1;
		instCFAccEnd[1] = settleDate;
		instCFDates = new TDate[2];
		instCFDates[0] = 1;
		instCFDates[1] = settleDate;
		instCFCoupons = new double[2];
		instCFCoupons[0] = 1;
		instCFCoupons[1] = 0;
		instCFAmorts = new double[2];
		instCFAmorts[0] = 1;
		instCFAmorts[1] = 0;
		claimCFDates = new TDate[3];
		claimCFDates[0] = 2;
		claimCFAmounts = new double[3];
		claimCFAmounts[0] = 2;
		
		claimCFDates[1] = settleDate;
		claimCFDates[2] = settleDate;
		claimCFAmounts[1] = 0;
		claimCFAmounts[2] = 0;

		

		KRateCurve  *strikePt =  instrument->get_strikes();
		int strikeSize = strikePt->size();
		exerStartDates = new TDate[strikeSize+1];
		exerStartDates[0] = strikeSize;
		exerEndDates = new TDate[strikeSize+1];
		exerEndDates[0] = strikeSize;
		exerStartStrikes = new double[strikeSize+1];
		exerStartStrikes[0] = strikeSize;
		exerEndStrikes = new double[strikeSize+1];
		exerEndStrikes[0] = strikeSize;

		optionDirections = new long[strikeSize+1];
		optionDirections[0] = strikeSize;
		exerTypes = new long[strikeSize+1];
		exerTypes[0] = strikeSize;
		option_barrier = new double[strikeSize+1];
		option_barrier[0] = strikeSize;
		optionType = new long[strikeSize+1];
		optionType[0] = strikeSize;

		long tempOptionType, tempExerType;
		if(CallOptionT *pt3 = dynamic_cast<CallOptionT*>(instrument.get()) )
		{		
			tempOptionType = HY_CAPITAL_CALL;			
			tempExerType = pt3->get_isamerican();
		}
		else
		{
			tempOptionType = HY_CAPITAL_PUT;
			PutOptionT *pt4 = dynamic_cast<PutOptionT*>(instrument.get()); 
			tempExerType = pt4->get_isamerican();
		}


		for(i=1; i<=strikeSize; i++)
		{
			exerStartDates[i] = strikePt->iDate(i-1);
			exerEndDates[i] = strikePt->iDate(i-1);
			exerStartStrikes[i] = strikePt->iRate(i-1);
			exerEndStrikes[i] = strikePt->iRate(i-1);
			optionDirections[i] = 1;
			exerTypes[i] = tempExerType;
			option_barrier[i] = 0;
			optionType[i] = tempOptionType;
		}
	}


	double lim1[2];
	lim1[0] = 1;
	lim1[1] = instrument->get_strikes()->iRate(0);	\

	double lim2[2];
	lim2[0] = 1;
	lim2[1] = 1;
	
	double vollim[2];
	vollim[0] = 1;
	vollim[1] = vll;

	double betaa[2];
	betaa[0] = 1;
	betaa[1] = beta;

	long assetProcessType[2];
	assetProcessType[0] = 1;
	assetProcessType[1] = 1;

	long ppya[2];
	ppya[0] =1;
	ppya[1] = ppy;

	TDate valueDates[3];
	valueDates[0] = 2;
	valueDates[1] = irCurve->valueDate();
	valueDates[2] = settleDate;

	char  outputStrings[250];
	double  outputNumbers[40];
	
	if(HYMCapitalWrapper(spotPrice,
						 divRefSpot,					 							
						 divDates,
						 divYields,
						 dps,
						 repoDates,
						 repoRates,
						 swapDates,
						 swapRates,
						 volRefSpot,
						 volDates,
						 volRates,
						 volShift,
						 instType,
						 notional,
						 recoveryRate,
						 instCFAccStart,
						 instCFAccEnd,
						 instCFDates,
						 instCFCoupons,
						 instCFAmorts,
						 claimCFDates,
						 claimCFAmounts,
						 exerStartDates,
						 exerEndDates,
						 exerStartStrikes,
						 exerEndStrikes,
						 optionDirections,
						 optionType,
						 option_barrier,
						 exerTypes,
						 lim1,
						 lim2,
						 vollim,
						 x,
						 lim,
						 betaa,
						 assetProcessType,
						 ppya,
						 valueDates,
						 outputStrings,
						 outputNumbers) == FAILURE)
	{
		GtoErrMsg("function failed.\n");
	}



	result->push_back( outputNumbers[1]);   //value
	result->push_back( outputNumbers[3]);   //delta
	result->push_back( outputNumbers[4]);   //gamma
	result->push_back( outputNumbers[5]);   //vega
	result->push_back( outputNumbers[13]);   //annuity
	result->push_back(0);   //nothing
	result->push_back( outputNumbers[7]);   //accruel
	result->push_back( outputNumbers[18]);   //asset
	
	free(divDates);
	free(divYields);
	free(repoDates);
	free(repoRates);
	free(swapDates);
	free(swapRates);
	free(volDates);
	free(volRates);
	free(instCFAccStart);
	free(instCFAccEnd);
	free(instCFDates);
	free(instCFCoupons);
	free(instCFAmorts);
	free(claimCFDates);
	free(claimCFAmounts );
	free(exerStartDates);
	free(exerEndDates);
	free(exerStartStrikes);
	free(exerEndStrikes);
	return result;
}

