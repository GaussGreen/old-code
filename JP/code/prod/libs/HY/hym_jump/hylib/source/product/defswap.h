//defswap.cpp

#ifndef	_DEFSWAP_H_
#define	_DEFSWAP_H_

#include "tree_asset.h"
#include "instruments.h"
#include "kratecurve.h"
#include "kvalarray.h"
#include "settleInfo.h"			//
#include "macros.h"


KValarray<double> NYDefaultProtection(double  spotStockPrice,
							   DDMap*  dividentYield,
							   KRateCurve*   repoCurve,
							   KRateCurve*   irCurve,
							   KRateCurve*   volCurve,
							   KRateCurve*   lowerStrike,
							   KRateCurve*   upperStrike,
							   double        recoveryRate,
							   BaseFunction* assetToStockMapping,
							   BaseFunction* assetProcess,
							   double        ppy,
							   double        beta);

KValarray<double> GeneralPricer(double  spotStockPrice,
							   DDMap*  dividentYield,
							   KRateCurve*   repoCurve,
							   KRateCurve*   irCurve,
							   KRateCurve*   volCurve,
							   KRateCurve*   volCurveShift,
							   Instrument*   instrument,
							   BaseFunction* assetToStockMapping,
							   BaseFunction* assetProcess,
							   double        ppy,
							   double        beta,
							   double		 vollim,					//HY4.1v	
							   double		 jumpProb,					//jump version		
							   KDate         settleDate,
							   bool			 isCVOption = false);		//HY3.4v



#endif



