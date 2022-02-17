/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file marketconvention.cpp
 *
 *  \brief file for the market convention object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#if defined(__USE_BASE_LIBRARY)

#include "gpbase/marketconvention.h"
#include "gpbase/definecst.h"

#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"

/// ARM Kernel
#include <glob/paramview.h>
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )


ARM_MarketConventionData currencyData[] =
{
///	Name	FwdRule				LiborTerm		LiborIndexDayCount	SpotDaysMMDayCount	FixedPayFreq	FixedDayCount	CashPayFreq		CashDayCount
	"ATS",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"AUD",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_365,		1,	KACTUAL_365,	K_SEMIANNUAL,	KACTUAL_365	,	K_QUARTERLY	,	KACTUAL_365,
	"BEF",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		2,	KACTUAL_360,	K_ANNUAL	,	KACTUAL_365	,	K_ANNUAL	,	K30_360,
	"CAD",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		0,	KACTUAL_360,	K_ANNUAL	,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
	"CHF",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"CZK",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	KACTUAL_360	,	K_SEMIANNUAL,	K30_360,
	"DEM",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"DKK",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		0,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"ESP",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"EUR",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	KACTUAL_FEB29,
	"FIM",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_SEMIANNUAL,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"FRF",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_360,		1,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	KACTUAL_FEB29,
	"GBP",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		0,	KACTUAL_365,	K_SEMIANNUAL,	KACTUAL_365	,	K_SEMIANNUAL,	KACTUAL_365,
	"GRD",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		2,	KACTUAL_360,	K_SEMIANNUAL,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
	"HKD",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_365,		0,	KACTUAL_365,	K_QUARTERLY	,	KACTUAL_365	,	K_QUARTERLY	,	KACTUAL_365,
	"HUF",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_360,		2,	KACTUAL_360,	K_QUARTERLY	,	KACTUAL_360	,	K_SEMIANNUAL,	K30_360,
	"IEP",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_SEMIANNUAL,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"ITL",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"JPY",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_SEMIANNUAL,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
	"NLG",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"NOK",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"NZD",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_SEMIANNUAL,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"PHP",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		1,	KACTUAL_360,	K_SEMIANNUAL,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
	"PLN",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		2,	KACTUAL_360,	K_ANNUAL	,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
	"PTE",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"SEK",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"SGD",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_365,		2,	KACTUAL_365,	K_SEMIANNUAL,	KACTUAL_360	,	K_SEMIANNUAL,	KACTUAL_360,
	"USD",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_360,		2,	KACTUAL_360,	K_SEMIANNUAL,	K30_360		,	K_SEMIANNUAL,	K30_360,
	"XEU",	K_MOD_FOLLOWING,	K_SEMIANNUAL,	KACTUAL_360,		2,	KACTUAL_360,	K_ANNUAL	,	K30_360		,	K_ANNUAL	,	K30_360,
	"ZAR",	K_MOD_FOLLOWING,	K_QUARTERLY	,	KACTUAL_365,		0,	KACTUAL_360,	K_QUARTERLY	,	KACTUAL_365	,	K_SEMIANNUAL,	K30_360,
};

const size_t CcyTableSize = sizeof(currencyData)/sizeof(currencyData[0]);


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketConvention
///	Routine: CreateCcyDataTable
///	Returns: MapStringToCurrencyData*
///	Action : create the static table of all currency data
/////////////////////////////////////////////////////////////////

ARM_MarketConvention::MapStringToCurrencyData* ARM_MarketConvention::CreateCcyDataTable()
{
	MapStringToCurrencyData* result = new MapStringToCurrencyData;
	string name;

	for( size_t i=0; i<CcyTableSize; ++i )
	{
		name = stringGetUpper( currencyData[i].itsName );
		ARM_MarketConventionData data = { 
			currencyData[i].itsName,
			currencyData[i].itsFwdRule,
			currencyData[i].itsLiborTerm,
			currencyData[i].itsLiborIndexDayCount,
			currencyData[i].itsSpotDays,
			currencyData[i].itsMMDayCount,
			currencyData[i].itsFixedPayFreq,
			currencyData[i].itsFixedDayCount,
			currencyData[i].itsCashPayFreq,
			currencyData[i].itsCashDayCount };
		if( !result->insert ( CC_NS(std,pair)< const string, ARM_MarketConventionData >( name, data ) ).second )
		{
			CC_Ostringstream os;
			os << "could not insert data " << i;
			ARM_THROW(ERR_INVALID_DATA, os.str() );
		}
	}
	return result;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketConvention
///	Routine: ReleaseTheCcyDataTable
///	Returns: void 
///	Action : release the static table of all currency data
/////////////////////////////////////////////////////////////////

void ARM_MarketConvention::ReleaseTheCcyDataTable() 
{
	delete ARM_MarketConvention::TheCcyDataTable;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketConvention
///	Initialise TheCcyDataTable
/////////////////////////////////////////////////////////////////

ARM_MarketConvention::MapStringToCurrencyData* ARM_MarketConvention::TheCcyDataTable = ARM_MarketConvention::CreateCcyDataTable();


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketConvention
///	Routine: ARM_MarketConvention
///	Returns: Constructor
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_MarketConvention::ARM_MarketConvention( const string& name )
{
	MapStringToCurrencyData* data = ARM_MarketConvention::TheCcyDataTable;
	MapStringToCurrencyData::iterator iter;

	if( (iter = data->find( stringGetUpper( name ) ) ) != data->end() )
		itsData = (*iter).second;
	else
		ARM_THROW( ERR_INVALID_DATA, " could not find the ccy " + name );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_MarketConvention
///	Routine: toString
///	Returns: 
///	Action : stringify the object
/////////////////////////////////////////////////////////////////
string ARM_MarketConvention::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << "\t ======> Currency Infos <======\n\n";
	os << "\t Currency Name      : " << itsData.itsName << "\n";
	os << "\t Forward Rule       : " << ARM_ParamView::GetMappingName( S_FORWARD_RULES, itsData.itsFwdRule ) << "\n";
	os << "\t Libor Term         : " << ARM_ParamView::GetMappingName( S_FREQUENCY,  itsData.itsLiborTerm ) << "\n";
	os << "\t Libor DayCount     : " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsData.itsLiborIndexDayCount ) << "\n";
	os << "\t Spot Days          : " << itsData.itsSpotDays << "\n";
	os << "\t MM DayCount        : " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsData.itsMMDayCount ) << "\n";
	os << "\t Fixed PayFrequency : " << ARM_ParamView::GetMappingName( S_FREQUENCY, itsData.itsFixedPayFreq ) << "\n";
	os << "\t Fixed DayCount     : " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsData.itsFixedDayCount ) << "\n";
	os << "\t Cash PayFrequency  : " << ARM_ParamView::GetMappingName( S_FREQUENCY, itsData.itsCashPayFreq ) << "\n";
	os << "\t Cash DayCount      : " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsData.itsCashDayCount ) << "\n";
	os << "\n\t ======> End of Currency Infos <======\n";
	return os.str();
}


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
