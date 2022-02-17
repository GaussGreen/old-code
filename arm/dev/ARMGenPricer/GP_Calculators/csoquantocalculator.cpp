#include "gpbase/removeidentifiedwarning.h"

#include "gpcalculators/csoquantocalculator.h"

// gpinfra
#include "gpinfra/gramnode.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/gensecurity.h"

#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpcalib/calibmethod.h"

CC_BEGIN_NAMESPACE(ARM)

//------------------------------------------------------------------------------

const string CallableQuantoSpreadOptionCreator::CQSOColNamesTable [] = {
	"RD",
	"SD",
	"ED",
	"IndexFund",
	"DF",
	"ItFund",
	"ItCpn",
	"Fees",
	"N",
	"Margin",
	"K",
	"L",
	"F",
	"C",
	"Index",
	"Coupon",
	"Funding",
	"Callable",
	"Exotic"
}; 


//! Constructor
CallableQuantoSpreadOptionCreator::CallableQuantoSpreadOptionCreator(
								  const ARM_Date& startDate, 
								  const ARM_Date& endDate			)
								  :itsStartDate(startDate),
								  itsEndDate(endDate) 
{
}

//! Fill schedule information
void CallableQuantoSpreadOptionCreator::setScheduleInformation(
						const ScheduleInformation& scheduleInfo)
{
	itsScheduleInformation = scheduleInfo;

}
	
//! Fill underlying information
void CallableQuantoSpreadOptionCreator::setUnderlyings(
						const ARM_Currency& cpnCcy, 
						ARM_INDEX_TYPE shortIndex,
						ARM_INDEX_TYPE longIndex,
						const ARM_Currency& fundCcy,
						ARM_INDEX_TYPE fundIndex)
{
	itsCpnCcy		= cpnCcy;
	itsShortIndex	= shortIndex;
	itsLongIndex	= longIndex;
	itsFundCcy		= fundCcy;
	itsFundIndex	= fundIndex;
}

//! Fill payoff information
void CallableQuantoSpreadOptionCreator::setPayoffInformation(
							const ARM_Curve&	nominal,
							const ARM_Curve&	cpnStrikes,
							const ARM_Curve&	cpnLeverageShort,
							const ARM_Curve&	cpnLeverageLong,
							const ARM_Curve&	cpnMin,
							const ARM_Curve&	cpnMax,
							const ARM_Curve&	fundMargin)
{
	itsNominal			= nominal; 
	itsCpnStrikes		= cpnStrikes; 
	itsCpnLeverageShort	= cpnLeverageShort;
	itsCpnLeverageLong	= cpnLeverageLong;
	itsCpnMin			= cpnMin;
	itsCpnMax			= cpnMax;
	itsFundMargin		= fundMargin;
}

//! Construct Callability information
void CallableQuantoSpreadOptionCreator::setCallabilityInformation(
								int	exerciseFreq,
								int	noticeGap,
								const ARM_Curve&	fees) 
{
	itsFees = fees; 
	// TO DO use other arguments to modify fees
}


//! True work: create a GenericSecurity
ARM_GenSecurityPtr CallableQuantoSpreadOptionCreator::create() const
{

	// CreateDateStrip
	ARM_Date maturity = itsEndDate;
	maturity.AddPeriodMult(itsScheduleInformation.frequency,1);
	
	ARM_DateStrip schedule(itsStartDate,
							maturity,
							itsScheduleInformation.frequency,
							itsScheduleInformation.dayCounter,
							itsScheduleInformation.resetCalendar.c_str(),
							K_MOD_FOLLOWING,							// fwdRule
							itsScheduleInformation.isAdjusted, 
							K_SHORTSTART,								// stubRule
							GETDEFAULTVALUE,							// resetGap
							itsScheduleInformation.frequency,			// payFreq	
							GETDEFAULTVALUE,							// payGap	
							itsScheduleInformation.paymentCalendar.c_str(),
							K_ADVANCE,									// resetTiming
							K_ARREARS,									// payTiming 
							1,											// adjFirstdate 	
							GETDEFAULTVALUESTR,							// refDateChar
							GETDEFAULTVALUE,							// stdSpotDays 
							GETDEFAULTVALUE,							// indexTerm 
							K_ACCRUED,									// accruedOrFull 
							9999										// firstDateFwdRule  
						   ); 


	// DealDescription arguments
	size_t rowsNb = schedule.size()+1;  
	size_t colsNb = 19;		//sizeof(CQSOColNamesTable)/sizeof(CQSOColNamesTable[0]);
	ARM_StringVector txt(rowsNb*colsNb);
	vector< ARM_GP_VALUE_TYPE > format(rowsNb*colsNb); 
	ARM_StringVector pricedColumns= ARM_StringVector(4);
	pricedColumns[0] = CQSOColNamesTable[15];
	pricedColumns[1] = CQSOColNamesTable[16];
	pricedColumns[2] = CQSOColNamesTable[17];
	pricedColumns[3] = CQSOColNamesTable[18];

	// Fill DealDesc
	for(size_t j=0; j<colsNb; ++j) 
	{
		txt[j]		= CQSOColNamesTable[j]; 
		format[j]	= ARM_STRING; 
	}

	int fundingDayCounter = itsFundCcy.GetLiborIndexDayCount();
	bool isSetInAdvance = (itsScheduleInformation.resetType==K_ADVANCE); 

	// Get tenors
	string tenor1, tenor2; 
	switch(itsLongIndex) 
	{
	case CMS1 : tenor1 = "1y";  break;
	case CMS2 : tenor1 = "2y";  break;
	case CMS3 : tenor1 = "3y";  break;
	case CMS4 : tenor1 = "4y";  break;
	case CMS5 : tenor1 = "5y";  break;
	case CMS6 : tenor1 = "6y";  break;
	case CMS7 : tenor1 = "7y";  break;
	case CMS8 : tenor1 = "8y";  break;
	case CMS9 : tenor1 = "9y";  break;
	case CMS10: tenor1 = "10y"; break;
	case CMS11: tenor1 = "11y"; break;
	case CMS12: tenor1 = "12y"; break;
	case CMS13: tenor1 = "13y"; break;
	case CMS14: tenor1 = "14y"; break;
	case CMS15: tenor1 = "15y"; break;
	case CMS16: tenor1 = "16y"; break;
	case CMS17: tenor1 = "17y"; break;
	case CMS18: tenor1 = "18y"; break;
	case CMS19: tenor1 = "19y"; break;
	case CMS20: tenor1 = "20y"; break;
	case CMS21: tenor1 = "21y"; break;
	case CMS22: tenor1 = "22y"; break;
	case CMS23: tenor1 = "23y"; break;
	case CMS24: tenor1 = "24y"; break;
	case CMS25: tenor1 = "25y"; break;
	case CMS26: tenor1 = "26y"; break;
	case CMS27: tenor1 = "27y"; break;
	case CMS28: tenor1 = "28y"; break;
	case CMS29: tenor1 = "29y"; break;
	case CMS30: tenor1 = "30y"; break;
	default:
		ARM_THROW(ERR_INVALID_INPUT, "Unrecognized 1st index in CallableQuantoSpreadOptionCreator!"); 
	}
	switch(itsShortIndex) 
	{
	case CMS1 : tenor2 = "1y";  break;
	case CMS2 : tenor2 = "2y";  break;
	case CMS3 : tenor2 = "3y";  break;
	case CMS4 : tenor2 = "4y";  break;
	case CMS5 : tenor2 = "5y";  break;
	case CMS6 : tenor2 = "6y";  break;
	case CMS7 : tenor2 = "7y";  break;
	case CMS8 : tenor2 = "8y";  break;
	case CMS9 : tenor2 = "9y";  break;
	case CMS10: tenor2 = "10y"; break;
	case CMS11: tenor2 = "11y"; break;
	case CMS12: tenor2 = "12y"; break;
	case CMS13: tenor2 = "13y"; break;
	case CMS14: tenor2 = "14y"; break;
	case CMS15: tenor2 = "15y"; break;
	case CMS16: tenor2 = "16y"; break;
	case CMS17: tenor2 = "17y"; break;
	case CMS18: tenor2 = "18y"; break;
	case CMS19: tenor2 = "19y"; break;
	case CMS20: tenor2 = "20y"; break;
	case CMS21: tenor2 = "21y"; break;
	case CMS22: tenor2 = "22y"; break;
	case CMS23: tenor2 = "23y"; break;
	case CMS24: tenor2 = "24y"; break;
	case CMS25: tenor2 = "25y"; break;
	case CMS26: tenor2 = "26y"; break;
	case CMS27: tenor2 = "27y"; break;
	case CMS28: tenor2 = "28y"; break;
	case CMS29: tenor2 = "29y"; break;
	case CMS30: tenor2 = "30y"; break;
	default:
		ARM_THROW(ERR_INVALID_INPUT, "Unrecognized 2nd index in CallableQuantoSpreadOptionCreator!"); 
	}

	for(size_t i=0; i<rowsNb-1; i++) 
	{

		double rd = (*(schedule.GetResetDates()))[i]; 
		CC_Ostringstream rdDesc;
		rdDesc << CC_NS(std,fixed) << rd;
		txt[(i+1)*colsNb   + 0]	= rdDesc.str(); 
		format[(1+i)*colsNb+ 0]	= ARM_DATE_TYPE; 

		double sd = (*(schedule.GetFlowStartDates()))[i]; 
		CC_Ostringstream sdDesc;
		sdDesc << CC_NS(std,fixed) << sd;
		txt[(i+1)*colsNb   + 1]	= sdDesc.str(); 
		format[(1+i)*colsNb+ 1]	= ARM_DATE_TYPE; 

		double ed = (*(schedule.GetFlowEndDates()))[i]; 
		CC_Ostringstream edDesc;
		edDesc << CC_NS(std,fixed) << ed;
		txt[(i+1)*colsNb   + 2]	= edDesc.str(); 
		format[(1+i)*colsNb+ 2]	= ARM_DATE_TYPE; 

		txt[(i+1)*colsNb   + 3]	= "LIBOR(YC_DOM,SD[i],ED[i])"; 
		format[(1+i)*colsNb+ 3]	= ARM_STRING; 

		txt[(i+1)*colsNb   + 4]	= "DF(YC_DOM,SD[i])"; 
		format[(1+i)*colsNb+ 4]	= ARM_STRING; 

		CC_Ostringstream itFundDesc;
		itFundDesc << CC_NS(std,fixed) 
			<< CountYears(fundingDayCounter,  (*(schedule.GetFlowStartDates()))[i], (*(schedule.GetFlowEndDates()))[i]);
		txt[(i+1)*colsNb   + 5]	= itFundDesc.str(); 
		format[(1+i)*colsNb+ 5]	= ARM_DOUBLE; 

		CC_Ostringstream itCpnDesc;
		itCpnDesc << CC_NS(std,fixed) 
			<< CountYears(itsScheduleInformation.dayCounter,  (*(schedule.GetFlowStartDates()))[i], (*(schedule.GetFlowEndDates()))[i]);
		txt[(i+1)*colsNb   + 6]	= itCpnDesc.str(); 
		format[(1+i)*colsNb+ 6]	= ARM_DOUBLE; 

		CC_Ostringstream feeDesc;
		feeDesc << CC_NS(std,fixed) << itsFees.Interpolate(rd);
		txt[(i+1)*colsNb   + 7]	= feeDesc.str(); 
		format[(1+i)*colsNb+ 7]	= ARM_DOUBLE; 

		CC_Ostringstream nomDesc;
		nomDesc << CC_NS(std,fixed) << itsNominal.Interpolate(rd);
		txt[(i+1)*colsNb   + 8]	= nomDesc.str(); 
		format[(1+i)*colsNb+ 8]	= ARM_DOUBLE; 

		CC_Ostringstream margDesc;
		margDesc << CC_NS(std,fixed) << itsFundMargin.Interpolate(rd);
		txt[(i+1)*colsNb   + 9]	= margDesc.str(); 
		format[(1+i)*colsNb+ 9]	= ARM_DOUBLE; 

		CC_Ostringstream kDesc;
		kDesc << CC_NS(std,fixed) << itsCpnStrikes.Interpolate(rd);
		txt[(i+1)*colsNb   +10]	= kDesc.str(); 
		format[(1+i)*colsNb+10]	= ARM_DOUBLE; 

		double w1 = itsCpnLeverageLong.Interpolate(rd); 
		double w2 = itsCpnLeverageShort.Interpolate(rd); 
		double leverage = 0.5*(w1+w2);
		
		CC_Ostringstream levDesc;
		levDesc << CC_NS(std,fixed) << leverage;
		txt[(i+1)*colsNb   +11]	= levDesc.str(); 
		format[(1+i)*colsNb+11]	= ARM_DOUBLE; 
		
		CC_Ostringstream cpnMinDesc;
		cpnMinDesc << CC_NS(std,fixed) << itsCpnMin.Interpolate(rd);
		txt[(i+1)*colsNb   +12]	= cpnMinDesc.str(); 
		format[(1+i)*colsNb+12]	= ARM_DOUBLE; 

		CC_Ostringstream cpnMaxDesc;
		cpnMaxDesc << CC_NS(std,fixed) << itsCpnMax.Interpolate(rd);
		txt[(i+1)*colsNb   +13]	= cpnMaxDesc.str(); 
		format[(1+i)*colsNb+13]	= ARM_DOUBLE; 

		if (leverage!=0.) 
		{
			w1/=leverage;
			w2/=leverage; 
		}

		CC_Ostringstream indexDesc;
		indexDesc << "SpreadCMS(LOCAL,SD[i]," << tenor1 << "," << tenor2 <<",";  
		indexDesc << CC_NS(std,fixed) << w1 << ","; 
		indexDesc << CC_NS(std,fixed) << w2 << ")"; 
		txt[(i+1)*colsNb   +14]	= indexDesc.str(); 
		format[(1+i)*colsNb+14]	= ARM_STRING; 

		txt[(i+1)*colsNb   +15]	= "N[i-1]*ItCpn[i-1]*Min(Max(F[i-1],L[i-1]*(Index[i";
		if(isSetInAdvance)
			txt[(i+1)*colsNb+15]+= "-1";
		txt[(i+1)*colsNb   +15]+= "]-K[i-1]) ),C[i-1])*DF[i]";
		format[(1+i)*colsNb+15]	= ARM_STRING; 

		txt[(i+1)*colsNb   +16]	= "N[i-1]*ItFund[i-1]*(IndexFund[i-1]+Margin[i-1])*DF[i]"; 
		format[(1+i)*colsNb+16]	= ARM_STRING; 

		txt[(i+1)*colsNb   +17]	= "Exercise(Funding[i]-Coupon[i],-Fees[i],Callable[i+1],Funding[i]-Coupon[i])"; 
		format[(1+i)*colsNb+17]	= ARM_STRING; 

		txt[(i+1)*colsNb   +18]	= "0"; 
		format[(1+i)*colsNb+18]	= ARM_DOUBLE; 
	}
	
	txt[i*colsNb   +17]	= "Funding[i]-Coupon[i]"; 
	format[i*colsNb+17]	= ARM_STRING; 
	
	txt[1*colsNb   +15]	= "0"; 
	format[1*colsNb+15]	= ARM_DOUBLE; 

	txt[1*colsNb   +16]	= "0"; 
	format[1*colsNb+16]	= ARM_DOUBLE; 

	txt[1*colsNb   +17]	= "Exercise(Funding[i]-Coupon[i],-Fees[i],Callable[i+1],IndexFund[i])"; 
	format[1*colsNb+17]	= ARM_STRING; 

	txt[1*colsNb   +18]	= "Callable[i]"; 
	format[1*colsNb+18]	= ARM_STRING; 

	ARM_DealDescriptionPtr dealDescription(
		new ARM_DealDescription(txt, format, rowsNb, colsNb, pricedColumns)
											); 

	// GenSecurity arguments
	ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL); 
	string paymodelName ="YC_DOM"; 
	bool exercBoundaryResetFlag = true; 
	bool otherPayoffsFlag = true; 
	bool ivFlag =true; 
	ARM_GenSecurityPtr security(
		new ARM_GenSecurity(dealDescription, 
							paymodelName,
							cstManager, 
							exercBoundaryResetFlag, 
							otherPayoffsFlag, 
							ivFlag )
								); 
	
	if(security.IsNull()) {
		ARM_THROW( ERR_PRICING_PB, ARM_USERNAME+" : error in CallableQuantoSpreadOptionCreator: security not created" );
	}
	
	return security; 
}



CC_END_NAMESPACE()
							