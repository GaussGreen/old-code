/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file basisconvert.cpp
 *
 *  \brief
 *
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date March 2006
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/gpvector.h"
#include "gpbase/numericconstant.h"
#include "gpbase/globalconstant.h"
#include "gpbase/curve.h"
#include "gpbase/curveconvert.h"  

/// gpinfra
#include "gpinfra/argconvdefault.h"

/// kernel
//#include <inst/swapleg.h>
//#include < mod/y2cmodel.h>

#include <iomanip>	

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_BasisConverter
///	Routine:  copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BasisConverter::ARM_BasisConverter(const ARM_BasisConverter& rhs)
:
	ARM_RootObject(rhs),
	itsDomesticDatestrip(rhs.itsDomesticDatestrip),
	itsFundingDatestrip(rhs.itsFundingDatestrip),
	itsForeignDatestrip(rhs.itsForeignDatestrip),
	itsDomesticNotional(rhs.itsDomesticNotional),
	itsForeignNotional(rhs.itsForeignNotional),
	itsForeignSpread(rhs.itsForeignSpread),
	itsDomesticZeroCurve(rhs.itsDomesticZeroCurve),
	itsForeignZeroCurve(rhs.itsForeignZeroCurve),
	itsDomesticDiscountZeroCurve(rhs.itsDomesticDiscountZeroCurve),
	itsForeignDiscountZeroCurve(rhs.itsForeignDiscountZeroCurve),
	itsForex(rhs.itsForex),
	itsDomesticCcy(rhs.itsDomesticCcy),
	itsForeignCcy(rhs.itsForeignCcy)
{
}
	///////////////////////////////////////////////////
///	Class  : ARM_BasisConverter
///	Routine:  consturctor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_BasisConverter::ARM_BasisConverter( const ARM_Currency& domesticCcy,
		const ARM_Currency& foreignCcy,
		const ARM_DateStripPtr&	domdatestrip,
		const ARM_DateStripPtr&	fordatestrip,
		const ARM_DateStripPtr&	funddatestrip,
		ARM_ZeroCurve* domesticZeroCurve,
		ARM_ZeroCurve* foreignZeroCurve,
		ARM_ZeroCurve* domesticDiscountZeroCurve,
		ARM_ZeroCurve* foreignDiscountZeroCurve,
		const ARM_Forex&        forex, 
		const std::vector<double>&	domesticNotional,
		const std::vector<double>&	foreignNotional,
		const std::vector<double>&	foreignSpread)
:ARM_RootObject()
	//itsDomesticDatestrip(domdatestrip),
	//itsFundingDatestrip(funddatestrip),
	//itsForeignDatestrip(fordatestrip),
	////itsAsOfDate(domesticZeroCurve->GetAsOfDate()),
	//itsDomesticNotional(domesticNotional),
	//itsForeignNotional(foreignNotional),
	//itsForeignSpread(foreignSpread),
	//itsDomesticZeroCurve(domesticZeroCurve),
	//itsForeignZeroCurve(foreignZeroCurve),
	//itsDomesticDiscountZeroCurve(domesticDiscountZeroCurve),
	//itsForeignDiscountZeroCurve(foreignDiscountZeroCurve),
	//itsDomesticCcy(domesticCcy),
	////itsForeignCcy(foreignCcy),
	////itsForex(forex)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_BasisConverter
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
std::vector<double> ARM_BasisConverter::ComputeDomMargin(ARM_SwapLeg* fundingLeg)
{
	/*std::vector<double>&  domResetDates		= itsDomesticDatestrip->GetResetDates(); 
	std::vector<double>&  domStartDates		= itsDomesticDatestrip->GetFlowStartDates();
	std::vector<double>&  domEndDates			= itsDomesticDatestrip->GetFlowEndDates();
	std::vector<double>&  dompayDates			= itsDomesticDatestrip->GetPaymentDates();

	std::vector<double>&  domfwdStartDates	= itsDomesticDatestrip->GetFwdRateStartDates();
	std::vector<double>&  domfwdEndDates		= itsDomesticDatestrip->GetFwdRateEndDates();

	std::vector<double>&  fundStartDates		= itsFundingDatestrip->GetFlowStartDates();
	std::vector<double>&  fundEndDates		= itsFundingDatestrip->GetFlowEndDates();
	std::vector<double>&  fundpayDates		= itsFundingDatestrip->GetPaymentDates();

	std::vector<double>&  fundfwdStartDates	= itsFundingDatestrip->GetFwdRateStartDates();
	std::vector<double>&  fundfwdEndDates		= itsFundingDatestrip->GetFwdRateEndDates();

	std::vector<double>&  forStartDates		= itsForeignDatestrip->GetFlowStartDates();
	std::vector<double>&  forEndDates			= itsForeignDatestrip->GetFlowEndDates();
	std::vector<double>&  forpayDates			= itsForeignDatestrip->GetPaymentDates();

	std::vector<double>&  forfwdStartDates	= itsForeignDatestrip->GetFwdRateStartDates();
	std::vector<double>&  forfwdEndDates		= itsForeignDatestrip->GetFwdRateEndDates();

	size_t index = 0;
	bool isFundingLegInputed = false;
	if(fundingLeg)
	{
		if (fundingLeg->GetPaymentDates()->size() < forpayDates->size())
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : funding leg: inconsistency between sizes, please advise" );
		index = fundingLeg->GetPaymentDates()->size() - forpayDates->size();
		isFundingLegInputed = true;
		ARM_Y2CModel Model(itsForeignZeroCurve, itsForeignDiscountZeroCurve);
		fundingLeg->SetModelVariable(NULL);
		fundingLeg->SetModel(&Model);
		double x = fundingLeg->ComputePrice();
	}

	size_t domsize	= itsDomesticDatestrip->size();
	int ccyDayCount = itsForeignCcy.GetLiborIndexDayCount();	
	double  spotValue = itsForex.GetMarketPrice();
	std::vector<double> domesticMargin(domsize);
	double startDomExchangeNotional,
		endDomExchangeNotional,
		startForExchangeNotional,
		endForExchangeNotional;

	double IT,fwdperiod,periodRatio,zcPay,
		zcStart,zcEnd,margin,nominal;

	/// ATTENTION ATTENTION:We suppose that both foreign first startDate
	/// and domestic first startDate are synchronized.
	for(int i(0),j(0),k(0); i < domsize; ++i)
	{
		ARM_Date startDate((*domStartDates)[i]); 
		ARM_Date endDate((*domEndDates)[i]); 

		/// Exchange Notional computing
		double startlag = startDate-itsAsOfDate;
		if(startlag < ARM_NumericConstants::ARM_TOLERENCE) continue;
		double endlag = endDate-itsAsOfDate;
		/// We suppose that the notionals are indexed on resetDates and 
		/// the interpolation method is STEPUP_LEFT
		startDomExchangeNotional = itsDomesticNotional[i] * itsDomesticDiscountZeroCurve->DiscountPrice(startDate);
		endDomExchangeNotional = itsDomesticNotional[i] *itsDomesticDiscountZeroCurve->DiscountPrice(endDate);

		startForExchangeNotional = itsForeignNotional[j]*itsForeignDiscountZeroCurve->DiscountPrice(startDate);	

		double forprice = 0.0, domprice = 0.0;	
		size_t forsize = itsForeignDatestrip->size();

		double dateTraget =  i == domsize -1 ? (*domEndDates)[i] : (*domResetDates)[i+1];
		while( j < forsize && (*forStartDates)[j] <= dateTraget )
		{
			if((*forStartDates)[j] > (*domResetDates)[i]- 3.0*ARM_GlobalConstant::ARM_SEVENDAYS_LAG){

				/// Foreign Leg
				ARM_Date fwdStartDate((*forfwdStartDates)[j]);
				ARM_Date fwdEndDate((*forfwdEndDates)[j]);
				IT = (*itsForeignDatestrip->GetInterestTerms())[j];
				fwdperiod = CountYearsWithoutException(ccyDayCount,fwdStartDate, fwdEndDate);
				periodRatio = IT/fwdperiod;

				
				zcPay = itsForeignDiscountZeroCurve->DiscountPrice(ARM_Date((*forpayDates)[j]));
				zcStart = itsForeignZeroCurve->DiscountPrice(fwdStartDate);
				zcEnd = itsForeignZeroCurve->DiscountPrice(fwdEndDate);
				margin = itsForeignSpread[j]*fwdperiod;
				nominal = itsForeignNotional[j];
				double flowPrice = nominal*periodRatio*zcPay*((zcStart/zcEnd -1.0)+margin);	
				if(isFundingLegInputed)
					flowPrice = (*fundingLeg->GetFlowsPV())[index+j];
				forprice += flowPrice;
			}
			++j;
		}
		
		endForExchangeNotional = itsForeignNotional[j-1]*itsForeignDiscountZeroCurve->DiscountPrice(endDate);
		
		/// notional foreign exchange
		forprice-= (startForExchangeNotional - endForExchangeNotional);	
		
		/// convert to domestic currency
		forprice*= spotValue;
		
		/// notional domestic exchange
		forprice+= (startDomExchangeNotional - endDomExchangeNotional);


		double O1 = 0.0;
		size_t fundsize  = itsFundingDatestrip->size();
		while( k < fundsize && (*fundfwdStartDates)[k] <= (*domEndDates)[i] - ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
		{
			if((*fundfwdStartDates)[k] > (*domStartDates)[i] - ARM_GlobalConstant::ARM_SEVENDAYS_LAG){

				/// domestic Leg
				ARM_Date fundfwdStartDate((*fundfwdStartDates)[k]);
				ARM_Date fundfwdEndDate((*fundfwdEndDates)[k]);
				IT = (*itsFundingDatestrip->GetInterestTerms())[k];
				fwdperiod = CountYearsWithoutException(ccyDayCount,fundfwdStartDate, fundfwdEndDate);
				periodRatio = IT/fwdperiod;

				zcPay = itsDomesticDiscountZeroCurve->DiscountPrice(ARM_Date((*fundpayDates)[k]));
				zcStart = itsDomesticZeroCurve->DiscountPrice(fundfwdStartDate);
				zcEnd = itsDomesticZeroCurve->DiscountPrice(fundfwdEndDate);
				nominal = itsDomesticNotional[i];
				double flowPrice = nominal*periodRatio*zcPay*(zcStart/zcEnd -1.0);				
				domprice += flowPrice;

				/// O1
				O1 += IT*zcPay*nominal;
			}
			++k;
		}
		
		/// value in real value like 0.5% 
		double dommargin = (forprice-domprice)/O1;

		domesticMargin[i] = dommargin;

	}*/

	return std::vector<double>(0/*domesticMargin*/);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

