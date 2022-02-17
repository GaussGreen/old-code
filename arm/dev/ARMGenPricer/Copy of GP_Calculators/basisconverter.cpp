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
#include <inst/swapleg.h>
#include < mod/y2cmodel.h>

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
		const ARM_GP_Vector&	domesticNotional,
		const ARM_GP_Vector&	foreignNotional,
		const ARM_GP_Vector&	foreignSpread)
:ARM_RootObject(),
	itsDomesticDatestrip(domdatestrip),
	itsFundingDatestrip(funddatestrip),
	itsForeignDatestrip(fordatestrip),
	itsAsOfDate(domesticZeroCurve->GetAsOfDate()),
	itsDomesticNotional(domesticNotional),
	itsForeignNotional(foreignNotional),
	itsForeignSpread(foreignSpread),
	itsDomesticZeroCurve(domesticZeroCurve),
	itsForeignZeroCurve(foreignZeroCurve),
	itsDomesticDiscountZeroCurve(domesticDiscountZeroCurve),
	itsForeignDiscountZeroCurve(foreignDiscountZeroCurve),
	itsDomesticCcy(domesticCcy),
	itsForeignCcy(foreignCcy),
	itsForex(forex)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_BasisConverter
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GP_Vector ARM_BasisConverter::ComputeDomMargin()
{
	ARM_GP_Vector*  domResetDates		= itsDomesticDatestrip->GetResetDates(); 
	ARM_GP_Vector*  domStartDates		= itsDomesticDatestrip->GetFlowStartDates();
	ARM_GP_Vector*  domEndDates			= itsDomesticDatestrip->GetFlowEndDates();
	ARM_GP_Vector*  dompayDates			= itsDomesticDatestrip->GetPaymentDates();

	ARM_GP_Vector*  domfwdStartDates	= itsDomesticDatestrip->GetFwdRateStartDates();
	ARM_GP_Vector*  domfwdEndDates		= itsDomesticDatestrip->GetFwdRateEndDates();

	ARM_GP_Vector*  fundStartDates		= itsFundingDatestrip->GetFlowStartDates();
	ARM_GP_Vector*  fundEndDates		= itsFundingDatestrip->GetFlowEndDates();
	ARM_GP_Vector*  fundpayDates		= itsFundingDatestrip->GetPaymentDates();

	ARM_GP_Vector*  fundfwdStartDates	= itsFundingDatestrip->GetFwdRateStartDates();
	ARM_GP_Vector*  fundfwdEndDates		= itsFundingDatestrip->GetFwdRateEndDates();

	ARM_GP_Vector*  forStartDates		= itsForeignDatestrip->GetFlowStartDates();
	ARM_GP_Vector*  forEndDates			= itsForeignDatestrip->GetFlowEndDates();
	ARM_GP_Vector*  forpayDates			= itsForeignDatestrip->GetPaymentDates();

	ARM_GP_Vector*  forfwdStartDates	= itsForeignDatestrip->GetFwdRateStartDates();
	ARM_GP_Vector*  forfwdEndDates		= itsForeignDatestrip->GetFwdRateEndDates();


	size_t domsize	= itsDomesticDatestrip->size();

	int ccyDayCount = itsForeignCcy.GetLiborIndexDayCount();
	
	double  spotValue = itsForex.GetMarketPrice();

	ARM_GP_Vector domesticMargin(domsize);
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

	}

	return ARM_GP_Vector(domesticMargin);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

