#pragma warning(disable : 4786)
#pragma warning(disable : 4541)

#include <ARM\libarm_local\firstToBeIncluded.h>

#include <comdef.h>
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
#else
#	include <ARM\libarm_local\undef_va_vars.h>
#endif

// A cause de windows.h
#ifdef GetObject
#undef GetObject
#endif

#include <ARM\libarm_frometk\arm_local_xgiga.h>

#include <ARM\libarm_local\ARM_local_glob.h>

#include <ARM\libarm_frometk\ARM_local_wsetk.h>
#include <ARM\libarm_frometk\ARM_local_etoolkitX.h>

#include <ARM\libarm_frometk\ARM_local_etoolkit.h>
#include <ARM\libarm_frometk\ARM_local_parsexml.h>
#include <ARM\libarm_local\ARM_local_persistent.h>
#include <ARM\libarm_frometk\arm_local_paesexml_calypso.h> 
#include <ARM\libarm_frometk\ARM_local_etoolkit.h>

#include <glob\dates.h>
#include "VariantTools.h"
#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;
// Global definitions

CCString RESCOMMSET = "<Request><Entity><COMMSET><dmCommodityId/><Id/><Type/><Name/><Name2/><AsOfDate/><DRuleId/><PYNature/><ShiftFormat/><Value/><AskValue/><ts>19800101 00:00:00</ts><SmileForm/><UnitsTime/><DistribType/><InterpMethod/><InterpOrder/><CommData/><FxVolChar><FXVOLCHAR><SmileAdj/><FXVolType/><Ccy1/><Ccy2/><CcyAlt1/><CcyAlt2/></FXVOLCHAR></FxVolChar><RecRateExpr/><VolNature/><BidValue/><DefaultLevel/></COMMSET></Entity></Request>";

int DATARETRIEVER = FFRETRIEVER;

ARM_EtkSoapClient* eToolkitSoapPtr = NULL;
ARM_Etoolkit* eToolkitPtr = NULL;


void checkETK()
{
	if (IsETKVersion())
	{
		if ( eToolkitPtr == NULL )
		{
			eToolkitPtr = new ARM_Etoolkit();
		}

		if ( eToolkitPtr->isConnecte() == 0 )
		{
			eToolkitPtr->Connect();
		}
	}
	else
	{
		if ( eToolkitSoapPtr == NULL )
		{
			eToolkitSoapPtr = new ARM_EtkSoapClient();
		}
	}
}

void set_etoolkit(const CCString& username,
				  const CCString& password,
				  const CCString& databasecontext,
				  const CCString& itConfigDomainsDir,
				  const CCString& itDomainName)


{
	if (IsETKVersion())
	{
		eToolkitPtr = new ARM_Etoolkit();

		eToolkitPtr->Init(username,password,databasecontext,itConfigDomainsDir,itDomainName);
	}
}


void connection_etoolkit()
{
	if (IsETKVersion())
	{
		if (eToolkitPtr)
			eToolkitPtr->Connect();
	}
}


void connection_etoolkit(const CCString& username,
						 const CCString& password,
						 const CCString& databasecontext,
						 const CCString& itConfigDomainsDir,
						 const CCString& itDomainName)
{
	if (IsETKVersion())
	{
		set_etoolkit(username,password,databasecontext,itConfigDomainsDir,itDomainName);

		connection_etoolkit();
	}

}


void connection_wsetoolkit(const CCString& pProd,
						   const CCString& pBase,
						   const CCString& pUserName)
{
	if (GetDataRetrieverVersion () == WSETKRETRIEVER)
	{
		eToolkitSoapPtr = new ARM_EtkSoapClient(pProd,pBase,pUserName);
	}
}



long etoolkit_getasof() 
{	
	CCString XMLResponse,messageList;

	CCString commande("s_base:GetAsOfDate");
	CCString requete("<Request></Request>");

	try
	{
		return etoolkit_execute(commande,requete,XMLResponse,messageList);
	}
    catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_execute(CCString command,
					  CCString xmlRequest,
					  CCString & xmlResponse,
					  CCString & messageList )
{
	try
	{
		checkETK();

		if (IsETKVersion())
			eToolkitPtr->Execute(command,xmlRequest,xmlResponse,messageList);
		else
			eToolkitSoapPtr->Execute(command,xmlRequest,xmlResponse,messageList);

		return ARM_OK;
	}

    catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_connecte()
{
	if (IsETKVersion())
	{
		return(eToolkitPtr->isConnecte());
	}

	return ARM_OK;
}



long deconnection_etoolkit()
{
	if (IsETKVersion())
	{
		if (eToolkitPtr!=NULL)
		{
			eToolkitPtr->Deconnect();
		}
	}
	return ARM_OK;
}


long shutdown_etoolkit()
{
	if (IsETKVersion())
	{
		if (eToolkitPtr!=NULL)
		{
			eToolkitPtr->Shutdown();

			delete eToolkitPtr;

			eToolkitPtr = NULL;
		}
	}
	return ARM_OK;
}


long kill_etoolkit()
{
	if (IsETKVersion())
	{
		if ( eToolkitPtr != NULL )
		{
			delete eToolkitPtr;

			eToolkitPtr = NULL;
		}
	}

	return(ARM_OK);
}



long etoolkit_getlastdatewarm(CCString& xmlResponse)
{
	CCString messageList;

	CCString commande("c_base:GetLastDateWarm");
	CCString requete=(CCString)"<Request></Request>";

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}

long etoolkit_getcommsetname(const CCString& index,
							 const CCString& ccy,
							 const CCString& cvname,
							 const CCString& type,
							 const CCString& impOrHist,
							 ARM_Date asof,
							 CCString& xmlResponse,
							 CCString& messageList,
							 const char* matuIndex,
							 const char* CallorPut)
{
	CCString ccName;
	CCString ccDate;

	// String de 6 caracteres pour passer le type
	CCString sType;

	sType = type;
	while ( sType.GetLen() < 6 )
		sType += " ";

	if (type == "ISSUER")
		ccName = (CCString)"ISSRVOL/" + index + "/" + ccy + "%";
	else if ( matuIndex == NULL)
	   ccName = (CCString)"IRFWDVOL/" + ccy + "/" + index + "/" + type + "%";
	else
	   ccName = (CCString)"SMILE/" + ccy + "/" + index + "/" + sType + "/" + matuIndex + "%/" + CallorPut + "/";
	
	char sDate[11];

    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());
	

	try
	{
		checkETK();

		if (IsETKVersion())
		{
			eToolkitPtr->GetCommsetName(ccName,"",
										sDate,
										impOrHist,cvname,xmlResponse,messageList);

			return ARM_OK;
		}
		else
		{
			eToolkitSoapPtr->GetCommsetName(ccName,"",
										sDate,
										impOrHist,cvname,xmlResponse,messageList);
			
			return ARM_OK;
		}
	}

	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_newgetcommsetname(const CCString& index,
								const CCString& ccy,
								const CCString& cvname,
								const CCString& type,
								const CCString& impOrHist,
								ARM_Date asof,
								CCString& xmlResponse,
								CCString& messageList,
								const char* matuIndex,
								const char* CallorPut)
{
	CCString commande("c_curves:GetCommset");

	CCString ccName;
	CCString ccDate;

	// String de 6 caracteres pour passer le type
	CCString sType;

	if (type == "ISSUER")
		ccName = (CCString)"ISSRVOL/" + index + "/" + ccy + "%";
	else if ( matuIndex == NULL)
	   ccName = (CCString)"IRFWDVOL/" + ccy + "/" + index + "/" + type;
	else
	   ccName = (CCString)"SMILE/" + ccy + "/" + index + "/" + type + "%/" + matuIndex + "%/" + CallorPut;
	
	char sDate[11];

    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString requete=(CCString)"<Request><Id>" + cvname + (CCString)"</Id><Type>" + impOrHist + (CCString)"</Type><Vol>";
	requete += ccName + (CCString)"</Vol><AsOfDate>" + sDate + (CCString)"</AsOfDate></Request>";

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_getrefrate(const CCString& source,
						 const CCString& ccy,
						 const CCString& index,
						 const CCString& tenor,
						 CCString& xmlResponse,
						 CCString& messageList)
{	
	try
	{
		checkETK();


		if (IsETKVersion())
		{
			eToolkitPtr->GetRefRate(source,
									ccy,
									index,
									tenor,
									xmlResponse,
									messageList);

			return ARM_OK;
		}
		else
		{			
			eToolkitSoapPtr->GetRefRate(source,
										ccy,
										index,
										tenor,
										xmlResponse,
										messageList);
			return ARM_OK;
		}
	}

	catch (...)
	{
		return ARM_KO;
	}
}



long etoolkit_setCurveId(const CCString& cvname)
{
	CCString xmlresponse;
	CCString messageList;

	CCString commande("s_market:FreeStaticCurve");
	CCString requete=(CCString)"<Request/>";

	CCString commande2("s_market:SetCurveId");
	CCString requete2=(CCString)"<Request><CurveId>" + cvname + (CCString)"</CurveId></Request>";

	try
	{
		etoolkit_execute(commande,requete,xmlresponse,messageList);
		etoolkit_execute(commande2,requete2,xmlresponse,messageList);

		return ARM_OK;
	}
	catch (...)
	{
		return ARM_KO;
	}
}



long etoolkit_getVolCurveByTenor(const CCString& sortRequete,
								 const CCString& cvname,
								 ARM_Date asof,
								 const CCString& typeForGVC,
								 const CCString& impOrHist,
								 CCString& xmlResponse,
								 CCString& messageList)
{
	char sDate[11];

    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString xmlReq = RESCOMMSET;
	CCString commande;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");

		if (impOrHist == "ISSRVOL")
			xmlReq.Replace("<Type/>", (CCString)"<Type>ISSRVOL</Type>");
		else if (impOrHist == "HISTVOL")
			xmlReq.Replace("<Type/>", (CCString)"<Type>HISTVOL</Type>");
		else
			xmlReq.Replace("<Type/>", (CCString)"<Type>IRFWDVOL</Type>");

	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + sortRequete + (CCString)"</Name>");

		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvname + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}



long etoolkit_getStringVolCurveByTenor(const CCString& sortRequete,
									   const CCString& cvname,
									   ARM_Date asof,
									   const CCString& typeForGVC,
									   CCString& xmlResponse,
									   CCString& messageList)
{
	CCString ccName;
	CCString ccName2;
	CCString ccDate;

	// String de 6 caracteres pour passer le type
	CCString sType;

	sType = "IRFWDVOL";

	ccName = sortRequete;
	ccName2 = "";

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString commande("s_market:GetVolCurve");
	CCString requete = (CCString)"<Request><Commodity>" + sortRequete + (CCString)"</Commodity><Type>" + typeForGVC + "</Type><AsOfDate>" + sDate + "</AsOfDate></Request>";

	try
	{
		checkETK();

		if (IsETKVersion())
		{
			eToolkitPtr->GetCommsetName(ccName,ccName2,sDate,sType,cvname,xmlResponse,messageList);

			return ARM_OK;
		}
		else
		{
			eToolkitSoapPtr->GetCommsetName(ccName,ccName2,sDate,sType,cvname,xmlResponse,messageList);
			
			return ARM_OK;
		}
	}
	catch (...)
	{
		return ARM_KO;
	}
}



long etoolkit_getSmileCurveByStrike(const CCString& strike,
									const CCString& request,
									const CCString& cvname,
									ARM_Date asof,
									const CCString& type,
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
	    xmlReq.Replace("<Type/>", (CCString)"<Type>" + type + (CCString)"</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + request + (CCString)"</Name>");
		
		if (strcmp(type,"SMILE") == 0)
		   xmlReq.Replace("<Name2/>", (CCString)"<Name2>" + strike + (CCString)"</Name2>");

		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvname + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}



long etoolkit_getCorrelCurveByMatu(const CCString& sortRequete,
								   const CCString& ccy2,
								   const CCString& index2,
								   const CCString& cvname,
								   ARM_Date asof,
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
	    xmlReq.Replace("<Type/>", "<Type>INDEXCOR</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + sortRequete + (CCString)"</Name>");
	    xmlReq.Replace("<Name2/>", (CCString)"<Name2>ZERO/" + ccy2 + (CCString)"/" + index2 + (CCString)"/1D</Name2>");
	    xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvname + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_getfxvol(const CCString& ccy1,
					   const CCString& ccy2,
					   const CCString& cvName,
					   ARM_Date date,
					   const CCString& impOrHist,
					   CCString& xmlResponse,
					   CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>" + impOrHist + (CCString)"</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>FXVOL/" + ccy1 + (CCString)"/" + ccy2 + (CCString)"</Name>");
	    xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}
	catch(...)
	{
		return ARM_KO;
	}
}



//	----------------------------------------------------------------------------------
CCString etoolkit_getXMLZCFromSummit(const CCString& index,
									 const CCString& currency,
									 const CCString& cvName,
									 ARM_Date aSdate)
{
	CCString xmlResponse;

	CCString msgList;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetZC((const char*)index,(const char*)currency,"",(const char*)cvName,aSdate,xmlOutput);

	if (xmlOutput != "")
	{
		return (CCString)(xmlOutput.c_str());
	}
	else
	{
		CCString commande("s_market:FreeStaticCurve");
		CCString requete=(CCString)"<Request/>";

		if (etoolkit_execute(commande,requete,xmlResponse,msgList)==ARM_KO)
		{
			// throw Exception(__LINE__, __FILE__,ERR_OBJECT_NULL,"Erreur lors de l'acces etoolkit");
			ARMTHROW(ERR_OBJECT_NULL,"Erreur lors de l'acces etoolkit");
		}

		char sDate[11];
		sprintf(sDate, "%04d%02d%02d", aSdate.GetYear(), aSdate.GetMonth(), aSdate.GetDay());
		
		CCString myMarket("s_market:GetZeroCurve");
		CCString myRequest;
		myRequest = (CCString)"<Request><Ccy>" + currency + (CCString)"</Ccy><Index>" + index + (CCString)"</Index><AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate><CurveId>" + cvName + (CCString)"</CurveId></Request>";

		if (etoolkit_execute(myMarket,myRequest,xmlResponse,msgList)==ARM_KO)
		{
			throw Exception(__LINE__, __FILE__,ERR_OBJECT_NULL,"Erreur lors de l'acces etoolkit");
		}

		return xmlResponse;
	}
}
//	----------------------------------------------------------------------------------
void etoolkit_getXMLZCFromCalypso(const std::string& index,
								  const std::string ccy,
								  const std::string& term,
								  const std::string& forceCurveName,
								  const std::string& pricingEnv,
								  const ARM_Date& date,
								  const std::string& xmlFileName,
								  std::string& xmlOutput)
{
	ARM_CalypsoToolkit::GetCurveZero(index,ccy,term,forceCurveName,pricingEnv,date,xmlFileName,xmlOutput); 
}

CCString etoolkit_getXMLMYAndZCFromSummit(const CCString& index,
										  const CCString& currency,
										  const CCString& cvName,
										  ARM_Date aSdate)
{
	CCString xmlResponse;

	CCString msgList;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetYLD((const char*)index,(const char*)currency,"",(const char*)cvName,aSdate,xmlOutput);

	if (xmlOutput != "")
	{
		return (CCString)(xmlOutput.c_str());
	}
	else
	{
		char sDate[11];

		sprintf(sDate, "%04d%02d%02d", aSdate.GetYear(), aSdate.GetMonth(), aSdate.GetDay());

		CCString myMarket("s_base:EntityRead");
		CCString myRequest;
		myRequest = (CCString)"<Request><Entity><MKTDATA><Id>" + cvName + (CCString)"</Id><Ccy>" + currency + (CCString)"</Ccy><dmIndex>" + index + (CCString)"</dmIndex><AsOfDate>" + (CCString) sDate + (CCString)"</AsOfDate></MKTDATA></Entity></Request>";

		if (etoolkit_execute(myMarket,myRequest,xmlResponse,msgList)==ARM_KO)
		{
			throw Exception(__LINE__, __FILE__,ERR_OBJECT_NULL,"Erreur lors de l'acces etoolkit");
	//		return NULL;
		}

		return xmlResponse;
	}
}



CCString etoolkit_getXMLObjectFromSummit(const CCString& idSummit,
										 const CCString& typeId)
{
	CCString xmlResponse("");
	CCString xmlResponse2("");
	CCString xmlMerge("");

	CCString tradeType;

	if (   (typeId == "SWOPT") 
		|| (typeId == "PRCS")
		|| (typeId == "PRCS2")
		|| (typeId == "ACCRUALOPTION")
		|| (typeId == "CRF")
		|| (typeId == "CRA")
		|| (typeId == "CDSOPT")
		)
		tradeType = "SWAPTION";
	else if (  (typeId == "IRG")
			|| (typeId == "SPDOPT")
			|| (typeId == "MATURITYCAP")
			|| (typeId == "RFTARN")
			|| (typeId == "MEMORYCAP")
			|| (typeId == "MEMORYSO")
			|| (typeId == "FXOPTSTRIP")
			|| (typeId == "FXSTRIP")
			|| (typeId == "CFXSTRIP")
			)
		tradeType = "CAPTR";
	else if (typeId == "FXOPT")
		tradeType = "FXOPT_TR";
	else if (  (typeId == "NTD")
			|| (typeId == "CDS")
			|| (typeId == "CDO2")
			|| (typeId == "CDO")
			|| (typeId == "EXOTIC"))
		tradeType = "EXOTIC";
	else if (typeId == "SWAP")
		tradeType = "SWAP";
	else if (typeId == "RNG_DOUBLE")
		tradeType = "SWAP";
	else
		return xmlResponse;

	CCString xmlReq;
	CCString messageList;

	CCString commande="s_base:EntityCreate";
	CCString requete=((CCString)"<Request><EntityName>" + tradeType + "</EntityName></Request>");

	try
	{
		etoolkit_execute(commande,requete,xmlResponse,messageList);

		xmlReq = xmlResponse;
		xmlReq.Replace("<TradeId/>", (CCString)"<TradeId>" + idSummit + (CCString)"</TradeId>");
		xmlReq.Replace("<Response>", (CCString)"<Request>");
		xmlReq.Replace("</Response>", (CCString)"</Request>");

		commande = "s_base:EntityRead";
		requete=xmlReq;

		etoolkit_execute(commande,requete,xmlResponse,messageList);
/*
		if (typeId == "IRG")
		{
			xmlResponse2 = etoolkit_getXMLCapCashFlowFromSummit(idSummit,typeId);
			xmlMerge = (CCString)"<ROOT>" + xmlResponse + xmlResponse2 + (CCString)"</ROOT>";
			xmlResponse = xmlMerge;
		}
*/	}
	catch (...)
	{
		return ("");
	}

	xmlResponse.Replace("<Entity>", "");
	xmlResponse.Replace("</Entity>", "");
	return xmlResponse;
}


long etoolkit_getfxcorrel(const CCString& ccy1,
						  const CCString& index,
						  const CCString& term,
						  const CCString& ccy2,
						  const CCString& cvName,
						  ARM_Date date,
						  CCString& xmlResponse,
						  CCString& messageList)
{

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetFXCORREL((const char*)index,(const char*)ccy1,(const char*)ccy1,(const char*)ccy2,(const char*)term,(const char*)cvName,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
		return ARM_OK;
	}
	else
	{
		ARM_XGigaToolKit::doGetFXCORREL((const char*)index,(const char*)ccy1,(const char*)ccy2,(const char*)ccy1,(const char*)term,(const char*)cvName,date,xmlOutput);
		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
			return ARM_OK;
		}
		else
		{
			char sDate[11];
			sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

			CCString commande;

			CCString tmpTerm;
			tmpTerm = term;
			while (tmpTerm.GetLen() < 4)
				tmpTerm += " ";

			CCString name2 = (CCString)"IRFWDVOL/" + ccy1 + (CCString)"/" + index + (CCString)"/IRG   /" + tmpTerm + (CCString)"/ /";

			CCString xmlReq = RESCOMMSET;

			try
			{
				xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
				xmlReq.Replace("<Type/>", "<Type>CORRELAT</Type>");
				xmlReq.Replace("<Name/>", (CCString)"<Name>FXVOL/" + ccy1 + (CCString)"/" + ccy2 + (CCString)"</Name>");
				xmlReq.Replace("<Name2/>", (CCString)"<Name2>" + name2 + (CCString)"</Name2>");
				xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

				commande = "s_base:EntityRead";

				try
				{
					etoolkit_execute(commande,xmlReq,xmlResponse,messageList);
				}
				catch(...)
				{
					xmlReq.Replace(ccy1 + (CCString)"/" + ccy2,ccy2 + (CCString)"/" + ccy1);

					try
					{
						etoolkit_execute(commande,xmlReq,xmlResponse,messageList);
					}
					catch(...)
					{
						return ARM_KO;
					}

					if (xmlResponse == "")
						return ARM_KO;
				}

				if (xmlResponse == "")
				{
					xmlReq.Replace(ccy1 + (CCString)"/" + ccy2,ccy2 + (CCString)"/" + ccy1);

					try
					{
						etoolkit_execute(commande,xmlReq,xmlResponse,messageList);
					}
					catch(...)
					{
						return ARM_KO;
					}

					if (xmlResponse == "")
						return ARM_KO;
				}

				return ARM_OK;
			}

			catch (...)
			{
				return ARM_KO;
			}
		}
	}
}


long etoolkit_getcorrel(const CCString& ccy1,
						const CCString& index1,
						const CCString& ccy2,
						const CCString& index2,
						const CCString& cvName,
						ARM_Date date,
						CCString& xmlResponse,
						CCString& messageList)
{
	CCString ccName;
	CCString ccName2;
	CCString ccDate;

	// String de 6 caracteres pour passer le type
	CCString sType;

	sType = "INDEXCOR";

	ccName = (CCString)"ZERO/" + ccy1 + (CCString)"/" + index1 + (CCString)"/1D/%/RATES";
	ccName2 = (CCString)"ZERO/" + ccy2 + (CCString)"/" + index2 + (CCString)"/1D";

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	try
	{
		checkETK();

		if (IsETKVersion())
		{
			eToolkitPtr->GetCommsetName(ccName,ccName2,sDate,sType,cvName,xmlResponse,messageList);

			return ARM_OK;
		}
		else
		{
			eToolkitSoapPtr->GetCommsetName(ccName,ccName2,sDate,sType,cvName,xmlResponse,messageList);
			
			return ARM_OK;
		}
	}
	catch (...)
	{
		return ARM_KO;
	}
}

long etoolkit_getMeanRev(const CCString& ccy,
						 const CCString& index,
						 const CCString& cvName,
						 const CCString& C_NumFactor,
						 ARM_Date date,
						 CCString& xmlResponse,
						 CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	if (strcmp((const char*)C_NumFactor,"2F") == 0)
	// 2 Facteur -> "MC2"
	{
		Name=(CCString)"MODELFAC/MC2/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/REVERSION";
	}
	else		// 3 Facteurs  -> FXHD
	{
		Name=(CCString)"MODELFAC/FXHD/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/REVERSION";
	}


	// ** GIGASPACE
	std::string xmlOutput;
	if (strcmp((const char*)C_NumFactor,"2F") == 0)
		ARM_XGigaToolKit::doGetMODELFACTOR((const char*)ccy,(const char*)index,"MC2","SWOPT","REVERSION",(const char*)cvName,date,xmlOutput);
	else
		ARM_XGigaToolKit::doGetMODELFACTOR((const char*)ccy,(const char*)index,"FXHD","SWOPT","REVERSION",(const char*)cvName,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
		return ARM_OK;
	}
	else
	{

		CCString commande;
		CCString xmlReq = RESCOMMSET;

		try
		{
			xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
			xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
			xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
			xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

			commande = "s_base:EntityRead";

			long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

			return retcode;
		}

		catch (...)
		{
			return ARM_KO;
		}
	}
}



long etoolkit_getCutOff(const CCString& ccy,
						const CCString& index,
						const CCString& cvName,
						const CCString& NumFactor,
						ARM_Date date,
						CCString& xmlResponse,
						CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	if (strcmp((const char*)NumFactor,"2F") == 0)
	// 2 Facteur -> "MC2"
	{
		Name = (CCString)"MODELFAC/MC2/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/FACTOR1";
	}
	else		// 3 Facteurs  -> FXHD
	{
		Name = (CCString)"MODELFAC/FXHD/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/FACTOR1";
	}

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "c_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_getFilterByStructId(const CCString& structId,
								  const CCString& desk,
								  CCString& xmlResponse,
								  CCString& messageList)
{
	CCString commande="s_trade:GetTradeListIter";
	
	CCString requete=(CCString)"<Request><FilterEntity><FILTER><TradeTypes><FLTTTYPES><TradeType>SWOPT</TradeType></FLTTTYPES><FLTTTYPES><TradeType>EXOTIC</TradeType></FLTTTYPES></TradeTypes>";
	requete = requete + (CCString)"<TradeFilterData><TFLTDATA><TradeStatus>+DONE+VER</TradeStatus>";
	requete = requete + (CCString)"<Desk>" + desk + (CCString)"</Desk>";
	requete = requete + (CCString)"<StructureId>" + structId + (CCString)"</StructureId>";
	requete = requete + (CCString)"</TFLTDATA></TradeFilterData></FILTER></FilterEntity></Request>";

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_getFilter(const CCString& filter,
						CCString& xmlResponse,
						CCString& messageList)
{
	CCString commande="s_trade:GetTradeListIter";
	
	CCString requete=(CCString) "<Request><FilterEntity><Filter>" + (CCString)filter + (CCString)"</Filter></FilterEntity></Request>";
	
	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_iterateTradeList(const CCString& handle,
							   const CCString& action,
							   CCString& xmlResponse,
							   CCString& messageList)
{
	CCString commande="s_trade:IterateTradeList";
	CCString requete=((CCString)"<Request><Iterator>" + handle + (CCString)"</Iterator><Action>" + action + (CCString)"</Action></Request>");

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}

long etoolkit_releaseIterateTradeList(const CCString& handle,
									  CCString& xmlResponse,
									  CCString& messageList)
{
	CCString commande="s_trade:ReleaseTradeListIter";

	CCString requete=((CCString)"<Request><Iterator>" + handle + (CCString)"</Iterator></Request>");

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{
		return ARM_KO;
	}
}


long etoolkit_getQmodParam(const CCString& ccy,
						   const CCString& index,
						   const CCString& cvName,
						   ARM_Date date,
						   CCString& xmlResponse,
						   CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	Name = (CCString)"MODELFAC/QMOD/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/FACTOR1";

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}

}


long etoolkit_getQFXParam(const CCString& ccy,
						  const CCString& index,
						  const CCString& cvName,
						  ARM_Date date,
						  CCString& xmlResponse,
						  CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	Name = (CCString)"MODELFAC/QMOD/SWOPT/" + ccy + (CCString)"/" + index + (CCString)"/FACTOR2";

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}

long etoolkit_getFxSmileByStrike(const CCString& strike,	// ex "1000.000"
								const CCString& request,	//"SMILE/" & Dev1 & "/" & Dev2 & "/FXOPT/C/"
								const CCString& cvname,		//"MO"
								ARM_Date asof,
								const CCString& type,		//"SMILE"
								CCString& xmlResponse,
								CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	
	CCString commande;
	CCString xmlReq = RESCOMMSET;
	CCString xmlReq2 = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>" + type + (CCString)"</Type>");

		CCString newRequest = request;
		newRequest.Replace("%","C");

		xmlReq.Replace("<Name/>", (CCString)"<Name>" + newRequest + (CCString)"</Name>");
		
		if (strcmp(type,"SMILE") == 0)
		   xmlReq.Replace("<Name2/>", (CCString)"<Name2>" + strike + (CCString)"</Name2>");

		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvname + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		if (strcmp(xmlResponse,"") == 0)
		{
			xmlReq2.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
			xmlReq2.Replace("<Type/>", (CCString)"<Type>" + type + (CCString)"</Type>");

			newRequest = request;
			newRequest.Replace("%","P");

			xmlReq2.Replace("<Name/>", (CCString)"<Name>" + newRequest + (CCString)"</Name>");
			
			if (strcmp(type,"SMILE") == 0)
			   xmlReq2.Replace("<Name2/>", (CCString)"<Name2>" + strike + (CCString)"</Name2>");

			xmlReq2.Replace("<Id/>", (CCString)"<Id>" + cvname + (CCString)"</Id>");

			commande = "s_base:EntityRead";

			etoolkit_execute(commande,xmlReq2,xmlResponse,messageList);	
		}	
		return ARM_OK;
	}

	catch (...)
	{
		return ARM_KO;
	}
}

long etoolkit_getFxSmile(const CCString& request,	//"SMILE/" & Dev1 & "/" & Dev2 & "/FXOPT/C/"
						 const CCString& cvname,
						 ARM_Date asof,
						 CCString& xmlResponse,
						 CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	try
	{
		checkETK();
	
		if (IsETKVersion())
		{
			eToolkitPtr->GetCommsetName(request,"",
										sDate,"SMILE",cvname,xmlResponse,messageList);
		}
		else
		{
			eToolkitSoapPtr->GetCommsetName(request,"",
										sDate,"SMILE",cvname,xmlResponse,messageList);
		}
		return ARM_OK;
	}

	catch (...)
	{
		return ARM_KO;
	}
}


bool IsETKVersion()
{
	return ( ( GetDataRetrieverVersion() == ETKRETRIEVER ) 
			|| (GetFallBackDataRetrieverVersion() == ETKRETRIEVER ) );
}


int GetFallBackDataRetrieverVersion()
{
	if (DATARETRIEVER < 10)
		return FFRETRIEVER;
	else
		return DATARETRIEVER % 10;
}

int GetDataRetrieverVersion()
{
	if (DATARETRIEVER < 10)
		return DATARETRIEVER;
	else
		return (DATARETRIEVER / 10) - 1;
}

void switchToETK(int withFallBack)
{
	if (withFallBack == 1)
		DATARETRIEVER = FFETKRETRIEVER;
	else
		DATARETRIEVER = ETKRETRIEVER;
}


void switchToWSETK()
{
	DATARETRIEVER = WSETKRETRIEVER;
}


void switchToFLATFILE()
{
	DATARETRIEVER = FFRETRIEVER;
}


CCString etoolkit_getXMLCapCashFlowFromSummit(const CCString& idSummit,
											const CCString& typeId)
{
	CCString xmlResponse("");

	if (!(typeId == "IRG"))
		return xmlResponse;

	CCString xmlReq;
	CCString messageList;

	CCString commande="c_trade:GetFlows";
	CCString requete=((CCString)"<Request><Entity><TradeId>" + idSummit + "</TradeId><TradeType>" + typeId + "</TradeType></Entity></Request>");

	try
	{
		etoolkit_execute(commande,requete,xmlResponse,messageList);
	}
	catch (...)
	{
		return ("");
	}

	return xmlResponse;
}


long etoolkit_getFixing(const CCString& source,
						const CCString& index,
						const CCString& term,
						const CCString& ccy,
						ARM_Date asof,
						CCString& xmlResponse,
						CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asof.GetYear(), asof.GetMonth(), asof.GetDay());

	CCString commande="s_base:EntityRead";
	CCString requete=((CCString)"<Request><Entity><REFRATE><Source>" + source + (CCString)"</Source><Ccy>" + ccy + (CCString)"</Ccy><dmIndex>" + index + (CCString)"</dmIndex><Term>" + term + (CCString)"</Term><Date>" + sDate + (CCString)"</Date></REFRATE></Entity></Request>");

	try
	{
		long retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		return retcode;
	}
	catch (...)
	{}

	return ARM_KO;
}



long etoolkit_getFactor(const CCString& ccy,
						const CCString& index,
						const CCString& cvName,
						const CCString& Type,
						const CCString& NumFactor,
						ARM_Date date,
						CCString& xmlResponse,
						CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	Name = (CCString)"MODELFAC/FXHD/" + Type + (CCString)"/" + ccy + (CCString)"/" + index + (CCString)"/" + NumFactor ;

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}



CCString etoolkit_getXMLSeasonMgrFromSummit(const CCString& index,
											const CCString& ccy,
											const CCString& cvname,
											ARM_Date date)
{
	CCString xmlResponse;
	CCString messageList;

	CCString requete = (CCString)"IRFWDVOL/" + ccy + (CCString)"/" + index + (CCString)"/IRG   /3M  / /" ;

	long retCode = etoolkit_getVolCurveByTenor(requete,
											   cvname,
											   date,
											   "HISTVOL",
											   "HISTVOL",
											   xmlResponse,
											   messageList);

	return xmlResponse;

}

VECTOR<string> etoolkit_getTradeListByStructId(const CCString& structureId, const CCString& desk)
{
	CCString requete=(CCString)"<Request><FilterEntity><FILTER><TradeTypes><FLTTTYPES><TradeType>SWOPT</TradeType></FLTTTYPES><FLTTTYPES><TradeType>EXOTIC</TradeType></FLTTTYPES></TradeTypes>";
	requete = requete + (CCString)"<TradeFilterData><TFLTDATA><TradeStatus>+DONE+VER</TradeStatus>";
	requete = requete + (CCString)"<Desk>" + desk + (CCString)"</Desk>";
	requete = requete + (CCString)"<StructureId>" + structureId + (CCString)"</StructureId>";
	requete = requete + (CCString)"</TFLTDATA></TradeFilterData></FILTER></FilterEntity></Request>";

	checkETK();

	if (GetDataRetrieverVersion () == WSETKRETRIEVER)
	{
		return eToolkitSoapPtr->GetTradeList(requete);
	}
	else
	{
		VECTOR<string> listXMLResponse;
		long retcode;

		CCString commande="s_trade:GetTradeListIter";
		CCString messageList;
		CCString xmlResponse;

		retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		HRESULT hr;
		VARIANT_BOOL bOK;
	
		MSXML2::IXMLDOMDocument *XMLDoc = NULL;
		MSXML2::IXMLDOMNodeList * resultList = NULL, * resultListRefValue = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * listItemAsset2 = NULL;

		try
		{
			hr = CoInitialize(NULL); 
//			SUCCEEDED(hr) ? 0 : throw hr;
			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = xmlResponse;
			_bstr_t tmpChaine ; 
			VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 

			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating XML document for getting PRCS \n" + xmlResponse);

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		CCString tradeListHandle;
		int nbTrades;

		try
		{
			XMLDoc->selectSingleNode(_bstr_t((const char *)"Response/Length"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				nbTrades = atoi((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			XMLDoc->selectSingleNode(_bstr_t((const char *)"Response/Iterator"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				tradeListHandle = ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}
		}
		catch(...)
		{		
			if (theNode) theNode->Release();
			if (resultList) resultList->Release();
			if (XMLDoc) XMLDoc->Release();

			hr = S_FALSE;

			CCString msg((CCString)"Error in XML parsing for getting first Iterator in PRCS \n" + xmlResponse);

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		long retCode = etoolkit_iterateTradeList(tradeListHandle,
											"BEGINNING",
											xmlResponse,
											messageList);

		for (int i = 0; i < nbTrades; i++)
		{
			retCode = etoolkit_iterateTradeList(tradeListHandle,
												"NEXT",
												xmlResponse,
												messageList);

			listXMLResponse.push_back((const char*)xmlResponse);
		}

		retCode = etoolkit_releaseIterateTradeList(tradeListHandle,
												   xmlResponse,
												   messageList);

		return listXMLResponse;
	}

}


VECTOR<string> etoolkit_getTradeListByBook(const CCString& book)
{
	CCString requete=(CCString)"<Request><FilterEntity><FILTER><TradeTypes><FLTTTYPES><TradeType>EXOTIC</TradeType></FLTTTYPES><FLTTTYPES><TradeType>EXOTIC</TradeType></FLTTTYPES></TradeTypes>";
	requete = requete + (CCString)"<TradeFilterData><TFLTDATA><TradeStatus>+DONE+VER</TradeStatus>";
	requete = requete + (CCString)"<Book>" + book + (CCString)"</Book>";
	requete = requete + (CCString)"</TFLTDATA></TradeFilterData></FILTER></FilterEntity></Request>";

	checkETK();

	if (GetDataRetrieverVersion () == WSETKRETRIEVER)
	{
		return eToolkitSoapPtr->GetTradeList(requete);
	}
	else
	{
		VECTOR<string> listXMLResponse;
		long retcode;

		CCString commande="s_trade:GetTradeListIter";
		CCString messageList;
		CCString xmlResponse;

		retcode = etoolkit_execute(commande,requete,xmlResponse,messageList);

		HRESULT hr;
		VARIANT_BOOL bOK;
	
		MSXML2::IXMLDOMDocument *XMLDoc = NULL;
		MSXML2::IXMLDOMNodeList * resultList = NULL, * resultListRefValue = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * listItemAsset2 = NULL;

		try
		{
			hr = CoInitialize(NULL); 
//			SUCCEEDED(hr) ? 0 : throw hr;
			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = xmlResponse;
			_bstr_t tmpChaine ; 
			VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 

			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating XML document for getting PRCS \n" + xmlResponse);

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		CCString tradeListHandle;
		int nbTrades;

		try
		{
			XMLDoc->selectSingleNode(_bstr_t((const char *)"Response/Length"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				nbTrades = atoi((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			XMLDoc->selectSingleNode(_bstr_t((const char *)"Response/Iterator"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				tradeListHandle = ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}
		}
		catch(...)
		{		
			if (theNode) theNode->Release();
			if (resultList) resultList->Release();
			if (XMLDoc) XMLDoc->Release();

			hr = S_FALSE;

			CCString msg((CCString)"Error in XML parsing for getting first Iterator in PRCS \n" + xmlResponse);

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		long retCode = etoolkit_iterateTradeList(tradeListHandle,
											"BEGINNING",
											xmlResponse,
											messageList);

		for (int i = 0; i < nbTrades; i++)
		{
			retCode = etoolkit_iterateTradeList(tradeListHandle,
												"NEXT",
												xmlResponse,
												messageList);

			listXMLResponse.push_back((const char*)xmlResponse);
		}

		retCode = etoolkit_releaseIterateTradeList(tradeListHandle,
												   xmlResponse,
												   messageList);

		return listXMLResponse;
	}
}


CCString etoolkit_getXMLResetMgrFromSummit(const CCString& index,
										   const CCString& source,
										   const CCString& ccy,
										   const CCString& term)
{
	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetRESET((const char*)ccy,(const char*)index,(const char*)source,(const char*)term,xmlOutput);
	if (xmlOutput != "")
	{
		return (CCString)(xmlOutput.c_str());
	}
	else
	{
		CCString xmlResponse;
		CCString messageList;

		long retCode = etoolkit_getrefrate(source,
										   ccy,
										   index,
										   term,
										   xmlResponse,
										   messageList);
		return xmlResponse;
	}
}


CCString etoolkit_getListTenors(const CCString& index,
								const CCString& ccy,
								const CCString& cvname,
								ARM_Date date,
								const CCString& cvtype)
{
	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetSMILETenorList((const char*)index,(const char*)ccy,(const char*)cvname,(const char*)cvtype,date,xmlOutput);

	if (xmlOutput != "")
	{
		return (CCString)(xmlOutput.c_str());
	}
	else
	{
		char sDate[11];
		sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

		CCString commande;
		CCString xmlReq;

		CCString messageList;
		CCString xmlResponse;

		try
		{
			xmlReq = (CCString)"<Request><Id>" + cvname + (CCString)"</Id><Type>SMILE</Type><Vol>SMILE/";
			xmlReq += ccy + (CCString)"/" + index + (CCString)"/" + cvtype + (CCString) "</Vol><AsOfDate>" + sDate + (CCString)"</AsOfDate><PorC>C</PorC></Request>";

			commande = "c_curves:GetTenor";

			long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

			return xmlResponse;
		}
		catch (...)
		{
			return "";
		}
	}
}

long etoolkit_getModelParam(ARM_Date date,
							const CCString& model,
							const CCString& type,
							const CCString& factorName,
							const CCString& ccy,
							const CCString& index,
							const CCString& cvName,
							CCString& xmlResponse,
							CCString& messageList)
{
	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", date.GetYear(), date.GetMonth(), date.GetDay());

	CCString Name;

	Name = (CCString)"MODELFAC/" + 
			model + (CCString)"/" + 
			type + (CCString)"/" + 
			ccy + (CCString)"/" + 
			index + (CCString)"/" +
			factorName;

	CCString commande;
	CCString xmlReq = RESCOMMSET;

	try
	{
		xmlReq.Replace("<AsOfDate/>", (CCString)"<AsOfDate>" + sDate + (CCString)"</AsOfDate>");
	    xmlReq.Replace("<Type/>", (CCString)"<Type>MODELFAC</Type>");
	    xmlReq.Replace("<Name/>", (CCString)"<Name>" + Name + (CCString)"</Name>");
		xmlReq.Replace("<Id/>", (CCString)"<Id>" + cvName + (CCString)"</Id>");

		commande = "s_base:EntityRead";

		long retcode = etoolkit_execute(commande,xmlReq,xmlResponse,messageList);

		return retcode;
	}

	catch (...)
	{
		return ARM_KO;
	}
}

/*-----------------------------------------------------------------------------------*/
/*---- End Of File ----*/
