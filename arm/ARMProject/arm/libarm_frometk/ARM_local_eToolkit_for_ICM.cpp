#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include <ARM\libarm_local\firstToBeIncluded.h>

#include <comdef.h>
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
#include <ARM\libarm_local\undef_va_vars.h>

#include "ARM\libarm_frometk\ARM_local_paesexml_Calypso.h"

// A cause de windows.h
#ifdef GetObject
#undef GetObject
#endif

#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libarm_frometk\ARM_local_etoolkit_for_ICM.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_local\ARM_local_persistent.h>

#include <ARM\libarm_frometk\eToolkitX.h>

#include <libCCatl\CCatl.h>

#include <glob\dates.h>

#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;




long etoolkit_GetDefProbCurve (const CCString& issuer,
							   ARM_Date asof,
							   const CCString& CurveId,
							   CCString& xmlResponse,
							   CCString& messageList)
{

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString commande("c_curves:GetDefaultProbCurve");
	CCString requete = (CCString)"<Request><Issuer>"+issuer+(CCString)"</Issuer><AsOfDate>"+sDate+(CCString)"</AsOfDate><CurveId>"+CurveId+"</CurveId><Ext>ON</Ext></Request>";

	try
	{
		etoolkit_setCurveId(CurveId);

		etoolkit_execute(commande,requete,xmlResponse,messageList);

		if (xmlResponse == "") return ARM_KO;

		return ARM_OK;
	}
	catch (...)
	{}

	return ARM_KO;
}


long etoolkit_GetVolFromSummit (ARM_Date asof,
							   const CCString& IssuerId,
							   const CCString& curveName,
							   const CCString& currency,
							   CCString& xmlResponse,
							   CCString& messageList)
{

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString commande("c_curves:GetDefaultProbCurve");
	CCString requete = (CCString)"<Request><Issuer>"+IssuerId+(CCString)"</Issuer><AsOfDate>"+sDate+(CCString)"</AsOfDate><CurveId>"+curveName+"</CurveId></Request>";

	try
	{
		etoolkit_setCurveId(curveName);

		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{}

	return ARM_KO;
}


long etoolkit_GetCorrFromSummit (ARM_Date asof,
							   const CCString& Name1,
							   const CCString& Name2,
							   const CCString& CurveId,
							   CCString& xmlResponse,
							   CCString& messageList)
{

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());


	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString) "<Type>ISSRCORR</Type>");
	    xmlReq.Replace("<Name/>",(CCString) "<Name>ISSRCORR/" + Name1 + (CCString)"/NCOPUL/1D</Name>");
	    xmlReq.Replace("<Name2/>",(CCString) "<Name2>"+Name2+ (CCString)"</Name2>");
	    xmlReq.Replace("<Id/>", (CCString)"<Id>" + CurveId + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		if (xmlResponse == "")
		{
		    xmlReq.Replace((CCString) "<Name>ISSRCORR/" + Name1 + (CCString)"/NCOPUL/1D</Name>",(CCString) "<Name>ISSRCORR/" + Name2 + (CCString)"/NCOPUL/1D</Name>");
		    xmlReq.Replace((CCString) "<Name2>"+Name2+ (CCString)"</Name2>",(CCString) "<Name2>"+Name1+ (CCString)"</Name2>");
		}

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{}

	return ARM_KO;
}

long etoolkit_GetRecoveryCurve (const CCString& SecId,
							   ARM_Date asof,
							   const CCString& CurveId,
							   CCString& xmlResponse,
							   CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString commande("c_curves:GetRecoveryRates");
	CCString requete = (CCString)"<Request>" +
					   (CCString)"<TradeType>BOND</TradeType>"+
					   (CCString)"<SecId>"+SecId+(CCString)"</SecId>"+
					   (CCString)"<MMCcy/>"+
					   (CCString)"<CurveId>"+CurveId+(CCString)"</CurveId>"+
					   (CCString)"<MMType/>"+
					   (CCString)"<Issuer/>"+
					   (CCString)"<AsOfDate>"+sDate+(CCString)"</AsOfDate>"+
					   (CCString)"</CommsetColumn>"+
					   (CCString)"</Request>";
				
	try
	{
		etoolkit_setCurveId(CurveId);

		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{}

	return ARM_KO;

}

long etoolkit_getXML_ZC_DP_FromSummit(const CCString& Issuer,
										  const CCString& currency,
										  const CCString& cvName,
										  ARM_Date aSdate,
										  CCString& xmlResponse,
										  CCString& messageList)
{

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", aSdate.GetYear(), aSdate.GetMonth(), aSdate.GetDay());
	
	CCString myMarket("c_curves:BuildDefProbCurve");
	CCString myRequest;

	myRequest = (CCString)"<Request>";
	myRequest = myRequest + (CCString)"<Issuer>" + Issuer + (CCString)"</Issuer>";
	myRequest = myRequest + (CCString)"<AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate>";
	myRequest = myRequest + (CCString)"<CurveId>" + cvName + (CCString)"</CurveId>";
	myRequest = myRequest + (CCString)"</Request>";

	try
	{
		etoolkit_setCurveId(cvName);

		if (etoolkit_execute(myMarket,myRequest,xmlResponse,messageList)==ARM_KO)
		{
			throw Exception(__LINE__, __FILE__,ERR_OBJECT_NULL,"Erreur lors de l'acces etoolkit");
		}

		return ARM_OK;
	}
	catch (...)
	{}

	return ARM_KO;

}