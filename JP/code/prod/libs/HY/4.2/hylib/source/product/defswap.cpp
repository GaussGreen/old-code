//defswap.cpp

#include "defswap.h"

// extern std::ofstream out("c:\\yy.dat");




KValarray<double> GeneralPricer(double  spotStockPrice,
							   DDMap*  dividentYield,
							   KRateCurve*   repoCurve,
							   KRateCurve*   irCurve,
							   KRateCurve*   volCurve,
							   KRateCurve*   volCurveShift,			//HY3.4v
							   Instrument*   instrument,
							   BaseFunction* assetToStockMapping,
							   BaseFunction* assetProcess,
							   double        ppy,
							   double        beta,
							   double		 vollim,				//HY4.1v
							   KDate         settleDate,
							   bool			 isCVOption)			//HY3.4v
{
	KDate	curDate, treeStartDate;
	KValarray<double>  result(10);									//HY3.4v
//	KValarray<double>  result(9);

	SettleInfo settlement(settleDate);								//
	instrument->set_settlement(settlement);

	KDate	maturity = instrument->get_expiryDate();
	KDate	valueDate = irCurve->valueDate();
	
	AssetTree 	tree(spotStockPrice,
					 *irCurve,
					 *volCurve,
					 *volCurveShift,				//HY3.4v
					 assetToStockMapping,
					 assetProcess,
					 *repoCurve,
					 beta,
					 vollim,				//HY4.1v
					 dividentYield,
					 isCVOption);

//	out<<(*irCurve);

	//set critical point, i.e., all payment date
	tree.insert_tp(valueDate);
	KVector(KDate) dates = instrument->get_datelist();
	for (int i= 0; i<dates.size();i++)
	{
		if(GetDaysDiff(valueDate,dates[i])>0)		//HY4.1v
		{
			tree.insert_tp(dates[i]);
		}
	}

	double assetVol0 = volCurve->get_rate(GTO_LINEAR_INTERP,valueDate);				//HY4.1v
	double timeShift = vollim * vollim/assetVol0/assetVol0 * 365.0;					//HY4.1v
	treeStartDate = valueDate.dFwd(-timeShift);										//HY4.1v
	tree.insert_tp(treeStartDate);													//HY4.1v

	tree.set_last_date(maturity);
	// build tree
	tree.init(ppy);
	
//	tree.initDebtDrift();				//HY4.2, moved to init()

//	tree.printTimeLine(out);

	//set up instrument
	instrument->set_tree(&tree);

//	CashFlowList fcCopy = instrument->get_cashflowlist();		//old annuity calculation
	
//	out<<fcCopy;

	// set rate to 1 to calculate pvbp
//	fcCopy.reset_rate(1.0);				//old annuity calculation
//	fcCopy.reset_rateType(HY_FIX);		//old annuity calculation
//	fcCopy.reset_amortPay(0.0);			//old annuity calculation

	KRateCurve* strike1 = instrument->get_strikes();
	KRateCurve* strike2;

	int AnnuCalFlag = 0;
	if(DefaultSwap *pt = dynamic_cast<DefaultSwap*>(instrument))
	{
		strike2 = pt->get_strike2();
		AnnuCalFlag = -1;
	}
	else if(Bond *pt1 = dynamic_cast<Bond*>(instrument))
	{
		strike2 = pt1->get_strike2();
		AnnuCalFlag = 1;
	}
	else
	{
		strike2 = strike1;
	}
	
//	DefaultSwap   dsAverageLife(strike1,strike2,1.0,fcCopy);		//old annuity calculation
//	dsAverageLife.set_tree(&tree);									//old annuity calculation

	Instrument* instrumentCopy;														//
	if(AnnuCalFlag !=0)																	//
	{
		instrumentCopy = instrument->clone();							//
		CashFlowList* cashflowlistCopy = instrumentCopy->get_cashflowlistRef();		//
		cashflowlistCopy->shift_rate(0.0001);										// 1bp
		instrumentCopy->set_tree(&tree);											//
	}
	
	while(tree.roll_back() == true)
	{     
		instrument->update();
//		dsAverageLife.update();			//old annuity calculation
		if(AnnuCalFlag !=0)
		{
			instrumentCopy->update();											//
		}
	}

	double settleFactor = irCurve->get_pv(settleDate);

//	double temp = instrument->get_spot_ts_value();

	result[0] = (instrument->get_spot_ts_value() + instrument->get_option_ts_value())/settleFactor;		//HY3.4v

//	result[1] = instrument->get_delta_value()/settleFactor;		//delta
//	result[2] = instrument->get_gamma_value()/settleFactor;		// Gamma
//	result[3] = instrument->get_vega_value()/settleFactor;		// Vega

	result[1] = (instrument->get_delta_value() + instrument->get_option_delta_value())/settleFactor;	//delta   HY3.4v

	result[2] = (instrument->get_gamma_value() + instrument->get_option_gamma_value())/settleFactor;	// Gamma  HY3.4v

	result[3] = (instrument->get_vega_value() + instrument->get_option_vega_value())/settleFactor;		// Vega   HY3.4v
	
//	result[4] = -dsAverageLife.get_spot_ts_value();				//pvbp, old annuity calculation

	double asset = tree.get_spotAsset();
	double stock = tree.get_spotStock();
	double assetVol = volCurve->get_rate(maturity);

	result[5] =((AssetToStockMapping1*)assetToStockMapping)->assetVolToStockVol(stock,assetVol);  //stock vol

    result[6] = instrument->get_accrual(settleDate,irCurve);
//  result[6] = instrument->get_accrual(valueDate,irCurve);				//test settlement days

	result[7] =asset;
	result[8] =assetToStockMapping->deriv(asset);

	result[9] = instrument->get_option_vega_value()/settleFactor;		// optVega   HY3.4v

	if(AnnuCalFlag)
	{
//		temp = instrumentCopy->get_spot_ts_value();
		double temp = instrumentCopy->get_spot_ts_value() + instrumentCopy->get_option_ts_value();	//HY3.4v
		double temp1 = instrumentCopy->get_accrual(settleDate,irCurve);
		result[4] =ABS(temp/settleFactor-AnnuCalFlag*temp1 - (result[0]-AnnuCalFlag*result[6]))*10000.0;		//
	}
	else
	{
		result[4] =0;		//
	}

//	out<<"do vega\n"<<std::endl;
		// do vega
//	double dvol = 0.00001;
//	KRateCurve  shiftVolCurve(*volCurve);
//	shiftVolCurve.shiftCurve(dvol);
	// set new vol curve in the tree
//	tree.set_vol_curve(shiftVolCurve);
	
	// build tree without rebuild timeline
//	tree.init(ppy,false);
//	tree.printTimeLine(out);
//	instrument->set_tree(&tree);

//	while(tree.roll_back() == true)
//	{     
//		instrument->update();
//	}

//	result[3] = (instrument->get_spot_ts_value()-temp)/dvol;      //Vega

	return result;
}
