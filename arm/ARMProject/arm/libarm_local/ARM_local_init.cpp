#include "ARM_local_init.h"

#include <ARM/libarm_frometk/arm_local_xgiga.h>
#include "tech_macro.h"	// listener : util>tech_macro


ARMLOCAL_Init::ARMLOCAL_Init (const char* pFileName) : CInitReader (pFileName)
{
}



void ARMLOCAL_Init::Init ()
{
	Read ();

	SetSection ("VIEWER");

	default_currency    = GetString("DEFAULT_CURRENCY");
	data_folder         = GetString("DATA_FOLDER");

	SetSection ("FX");

	string currencies	= GetString("CURRENCIES");
	string indirects	= GetString("INDIRECT");

	long fxPos=0;
	long indPos=0;
	long fxPosOld=0;
	long indPosOld=0;

	while ((fxPos != -1) && (indPos != -1))
	{
		fxPos = currencies.find_first_of(",",fxPosOld);
		indPos = indirects.find_first_of(",",indPosOld);

		string currency = string(currencies,fxPosOld,fxPos-fxPosOld);
		string indirect = string(indirects,indPosOld,indPos-indPosOld);
		currencyList.push_back(currency);
		indirectList.push_back(atol(indirect.c_str()));

		fxPosOld = fxPos+1;
		indPosOld = indPos+1;
	}

	try
	{
		SetSection ("WS");

		wsetkprod		= GetString("PROD_WSETK_WSDL_URL");
		wsetkrec		= GetString("REC_WSETK_WSDL_URL");
		wspricingprod	= GetString("PROD_WSPRICING_WSDL_URL");
		wspricingrec	= GetString("REC_WSPRICING_WSDL_URL");
	}
	catch(...)
	{
	}

	try
	{
		SetSection ("GIGASPACES");
		string useGigaSpaces = GetString("USEGIGASPACES");
		string spaceURL	= GetString("SPACEURL");

		if (useGigaSpaces == "YES")
			ARM_XGigaToolKit::init(spaceURL);
	}
	catch(...)
	{
	}

	try
	{
		SetSection ("AUDIT");
		bool isActivated = GetString("ACTIVATED")=="YES" ? 1 : 0;
		string listenerIp = GetString("LISTENERIP");
		bool useBuffer = GetString("BUFFER")=="YES" ? 1 : 0;
		int bufferSize = GetInt("BUFFERSIZE");

		if(isActivated)
			Logger::instance().enable();
		else
			Logger::instance().disable();
		Logger::instance().setConfig(listenerIp, useBuffer, bufferSize);
	}
	catch(...)
	{
		Logger::instance().disable();
	}

}




string ARMLOCAL_Init::ToString ()
{
	return "";
}
