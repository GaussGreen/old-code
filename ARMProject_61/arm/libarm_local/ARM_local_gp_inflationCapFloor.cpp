/*! \file ARM_local_gp_inflationCapFloor.h
 *
 *
 *	\author  F. Poitou
 *
 */


/*
mising classes :
	Variant
	PeriodParser
	StringFormatter
*/
/*
GenericDate::GenericDate(const Variant& vData){
	if (vData.isDate()) {
		itsIds = ARM_CountedPtr<DoubleDate>(new DoubleDate(vData.asDate()));
	}
	else if (vData.isString()){
		itsIds = ARM_CountedPtr<StringDate>(new StringDate(StringFormatter::toUpperCase(vData.asString())));
	}
	else if (vData.isDouble()) {
		itsIds = ARM_CountedPtr<DoubleDate>(new DoubleDate(vData.asDouble()));
	}
	else ARMTHROW(ERR_INVALID_ARGUMENT,"Date must an XL date or a string describing a period");


}
//! -----------------------------------------------------------------------
ARM_Date DateFactory::create(const ARM_CountedPtr<InputDateStructure>& dt)
{
	ARM_CountedPtr<DoubleDate> ddt = idt ;
	if (!ddt.IsNull())
		return ARM_Date(ddt->dDate());
	else
	{
		ARM_CountedPtr<StringDate> sdt = idt;
		if (!sdt.IsNull())
		{
			if (sdt->sDate().find("/") != std::string::npos)
			{
				return DateParser::parse(sdt->sDate(),"dd/mm/yyyy");
			} else {
				Period lpPer = PeriodParser::parse(sdt->sDate());
				return DateFactory::today().AddPeriodMult(lpPer.unit(), lpPer.length());
			}
		}
		else ARMTHROW(ERR_INVALID_ARGUMENT,"Date must be an XL date or a string describing a period");
	}
}

//! -----------------------------------------------------------------------
ARM_Date DateFactory::create(const ARM_CountedPtr<InputDateStructure>& dt,
						ARM_Date refDate,
						char * calendar,
						int settlementDays,
						int rollingConvention)
{
	ARM_CountedPtr<DoubleDate> ddt = idt ;
	if (!ddt.IsNull())
		return ARM_Date(ddt->dDate());
	else
	{
		ARM_CountedPtr<StringDate> sdt = idt;
		if (!sdt.IsNull())
		{
			const string & strDate = sdt->sDate();
			Period lpPer = PeriodParser::parse(strDate);

			ARM_Date ldDate = refDate.AddPeriodMult(lpPer.unit(), lpPer.length(), calendar) ;
			if (lpPer.units() != Days )
				ldDate = ldDate.AddPeriodMult(K_DAYLY, settlementDays, calendar) ;
			return ldDate;
		}
		else ARMTHROW(ERR_INVALID_ARGUMENT,"Date must be an XL date or a string describing a period");
	}
}

//! -----------------------------------------------------------------------
ARM_Date DateFactory::create(const ARM_CountedPtr<InputDateStructure>& dt, ARM_Date adRefDate, int settlementDays)
{
	ARM_CountedPtr<DoubleDate> ddt = idt ;
	if (!ddt.IsNull())
		return ARM_Date(ddt->dDate());
	else
	{
		ARM_CountedPtr<StringDate> sdt = idt;
		if (!sdt.IsNull())
		{
			const string & strDate = sdt->sDate();
			Period lpPer = PeriodParser::parse(strDate);

			ARM_Date ldDate = adRefDate.AddPeriodMult(lpPer.unit(), lpPer.length());
			if (lpPer.units() != Days )
				ldDate = ldDate.AddPeriodMult(K_DAYLY, settlementDays);

			return ldDate;
		}
		else ARMTHROW(ERR_INVALID_ARGUMENT,"Date must be an XL date or a string describing a period");
	}
}

//! -----------------------------------------------------------------------
ARM_Date DateFactory::create(const ARM_CountedPtr<InputDateStructure>& dt, ARM_Date refDate)
{
	ARM_CountedPtr<DoubleDate> ddt = idt ;
	if (!ddt.IsNull())
		return ARM_Date(ddt->dDate());
	else
	{
		ARM_CountedPtr<StringDate> sdt = idt;
		if (!sdt.IsNull())
		{
			const string & strDate = sdt->sDate();
			Period lpPer = PeriodParser::parse(strDate);
			return ARM_Date.AddPeriodMult(lpPer.unit(), lpPer.length());
		}
		else ARMTHROW(ERR_INVALID_ARGUMENT,"Date must be an XL date or a string describing a period");
	}
}
ARM_Date DateFactory::today()
{
	ARM_Date d(0.);
	d.Today()
	return d ;
}
//! -----------------------------------------------------------------------
*/

/* !
 * function to create a inflation capfloor
 *
 */
#include "ARM_local_gp_inflationCapFloor.h"

#include "util\fromto.h"

#include <gpinfra\argconvdefault.h>
#include "ARM_local_wrapper.h"
#include "ARM_local_persistent.h"
#include "ARM_local_glob.h"
#include "GP_Inflation\gpinflation\infdata.h"
#include "gpinflation/infPayOff.h"


using namespace ARM ;

extern long ARMLOCAL_InfCapFloor_Load(
	double startDate,
	double endDate,
	const CCString& indexName,
	int capOrFloor,
	double strike,
	int interpType,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	double firstReset,
	ARM_result&	result,
	long objId 	)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CountedPtr<InfCapFloor> infCapFloor ;

	try
	{

		/// dateConversion
		char  startDateChar[11];
		char  endDateChar[11];

		Local_XLDATE2ARMDATE( startDate, startDateChar );
		Local_XLDATE2ARMDATE( endDate, endDateChar );

		char* resetCal	= resetCalendar.c_str();
		char* payCal	= payCalendar.c_str();



		ARM_CountedPtr<ARM_Currency>	ccy(new ARM_Currency	( InfData::GetCurrency( indexName.c_str() ) ));
		int				tmpDayCount		= ARM_ArgConv_LgNameDayCount.GetNumber("30/360");
		int				tmpAdjFirstRule	= ARM_ArgConv_IntRule.GetNumber("UNADJ");


		ARM_GP_Vector nominals(1,1.);
		bool monthBegin(true) ;
		Period numGap(resetGap,K_DAILY) ;
		Period demGap(resetGap,K_DAILY) ;
		ARM_CountedPtr <YOYPayOff> payoff( new  YOYPayOff());
		int resetRollConv = fwdRule, payRollConv = fwdRule ;//check this out
		Period tenor(1, K_ANNUAL) ;
		ARM_CountedPtr<ARM_InfIdx> cpi(new ARM_InfIdx(indexName.c_str()));

		ARM_CountedPtr<InfPerformance> cpiPerf ( new InfPerformance(cpi,tenor));



		ARM_GP_T_Vector<ARM_CountedPtr < CashFlow > >
		vCahFlows = InfFloatingCouponVector(nominals,payoff, cpiPerf, (ARM_Date) startDateChar,
							    (ARM_Date) endDateChar, resetFreq,
							    payFreq, dayCount,
							    fwdRule, stubRule,
								intRule, resetRollConv,
								resetCal, payRollConv,
								payCal, monthBegin,demGap,numGap);

		infCapFloor = ARM_CountedPtr<InfCapFloor>( new InfCapFloor( vCahFlows, capOrFloor, strike, interpType, firstReset ));



		delete resetCal;
		delete payCal;

		/// assign object
		if( !assignObject( infCapFloor.operator->(), result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}

}

