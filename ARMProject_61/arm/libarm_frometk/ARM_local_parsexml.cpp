/*
First Version  For Parse CRF For Calypso 8004
  */


#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_frometk\arm_local_xgiga.h>

#include <ARM\libarm_local\arm_local_glob.h>
#include <ARM\libarm_frometk\arm_local_parsexml.h>
#include <ARM\libarm_frometk\arm_local_parsexml_nt_security.h>
#include <ARM\libarm_frometk\PaserManagerUtilities.h>
#include <ARM\libarm_frometk\arm_local_parsexml_credit.h>

#include <ARM\libarm_local\ARM_local_swap.h>
#include <ARM\libarm_local\ARM_local_irindex.h>
#include <ARM\libarm_local\ARM_local_ccy.h>
#include <ARM\libarm_local\ARM_local_zccurve.h>

#include <ARM\local_xlarm\ARM_local_interglob.h>
/// headers to remove definition of va_start and va_end as this is redefined later on!
/// handle this with care after sorting out why this is so!!
#include <ARM\libarm_local\undef_va_vars.h>

#include <glob\linalg.h>
#include <glob\dates.h>
#include <crv\zeroint.h>
#include <crv\volint.h>
#include <crv\volspline.h>
#include <crv\volcube.h>
#include <crv\zerointspreaded.h>
#include <ccy\currency.h>
#include <inst\swapleg.h>
#include <inst\swaption.h>
#include <inst\option.h>
#include <inst\forex.h>
#include <inst\powrev.h>
#include <inst\irindex.h>
#include <inst\fixleg.h>
#include <inst\flexaccretswaption.h>
#include <inst\optionportfolio.h>
#include <inst\corridorleg.h>
#include <inst\spreadoption.h>
#include <ICMKernel\crv\icm_defaultcurve.h>
#include <util\fromto.h>
#include <inst\stripoption.h>
#include <inst\stripdigitaloption.h>
#include <inst\corridorDblCondition.h>

#include <gpbase\gpvector.h>
#include <gpbase\curve.h>
#include <gpbase\curveconvert.h>
#include <gpbase\globalportfolio.h>

/// gpcalculators
#include <GP_Calculators\gpcalculators\typedef.h>
#include <GP_Calculators\gpcalculators\argconvdefault.h>
#include <gpcalculators\crfcalculator.h>
#include <gpcalculators\maturitycapcalculator.h>
#include <gpcalculators\tarncalculator.h>
#include <gpcalculators\bermudaswaptioncalculator.h>
#include <gpcalculators\callablesnowballcalculator.h>
#include <gpcalculators\captioncalculator.h>
#include <gpcalculators\fxvanillacalculator.h>

using ARM::ARM_CRFCalculator;
using ARM::ARM_MaturityCapCalculator;
using ARM::ARM_TARNCalculator;
using ARM::ARM_BermudaSwaptionCalculator;
using ARM::ARM_CallableSnowBallCalculator;
using ARM::ARM_CaptionCalculator;
using ARM::ARM_FXVanillaCalculator;

/// gpbase
using ARM::ARM_GlobalPortfolio;

/// gpinfra
#include <gpinfra\pricingmodel.h>
#include <gpinfra\gensecurity.h>
#include <gpcalib\calibmethod.h>
#include <gpinfra\pricingadviser.h>
#include <gpinfra\modelparams.h>
#include <gpinfra\modelparamtype.h>
#include <gpinfra\curvemodelparam.h>
using ARM::ARM_ModelParam;
using ARM::ARM_CurveModelParam;
using ARM::ARM_ModelParamType;
using ARM::ARM_GP_Vector;

#include <gpinflation\seasonmanager.h>
#include <gpinflation\resetmanager.h>
#include <gpinflation\infcurv.h>
using ARM::ARM_SeasonalityManager;
using ARM::ARM_ResetManager;
using ARM::ARM_InfCurv;

#include <gpmodels\enummodel.h>
using ARM::ARM_PricingModelType ;
using ARM::ARM_SigmaCalibType ;


#include <wtypes.h>
#include <time.h>

#import <msxml3.dll> raw_interfaces_only
using namespace MSXML2;

#include <ARM\libarm_frometk\arm_local_etoolkit.h>
#include <ARM\libarm_frometk\arm_local_parsexml_common.h> 

#include <ARM\libarm_frometk\arm_local_parsexml_util.h> 
#include <ARM\libarm_frometk\ParserManager.h>
#include <ARM\libarm_frometk\CalypsoSwapLeg.h>
#include <ARM\libarm_frometk\CalypsoOptionSchedule.h>

#include "XMLTools.h"
#include "VariantTools.h"

using namespace etoolkit;
//reading from file
#include <iostream>
#include <fstream>
using namespace std;
// end reading from file



// A cause de windows.h
#ifdef GetMessage
#undef GetMessage
#endif



// For Retrieving A Summit Flow

void RetrieveSummitFlow(const CCString& tradeId, const CCString& tradeType,
						const CCString& curveId,
						const ARM_Date& inDate, 
						ARM_SummitFlow* outFlow)
{
	CCString request;
	CCString xmlReq;

	CCString xmlResFlows;
	CCString msgList;
	CCString command;


	long retCode = ARM_OK;

	if ( curveId != CCString((const char *) "" ) )
	{
	   retCode = etoolkit_setCurveId(curveId);

	   if ( retCode == ARM_KO )
	   {
		  throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			       (char *) "ETK Error: RetrieveSummitFlow:etoolkit_setCurveId(curveId) : Failed!!");
	   }
    }
 
	command = "c_trade:GetFlows";

	xmlReq = CCString((const char *)"<Request><Entity><TradeId>")+CCString((const char *) tradeId)
		     +CCString((const char *)"</TradeId><TradeType>")+CCString((const char *) tradeType)
			 +CCString((const char *)"</TradeType></Entity></Request>");

 
	retCode = etoolkit_execute(command, xmlReq, xmlResFlows, msgList);

    if ( retCode == ARM_KO )
	{
	   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		               (char *) "ETK Error: RetrieveSummitFlow:etoolkit_execute: Failed!!");   
    }
	
	/*  //->test
	char* chaineXMLOut = xmlResFlows;
	FILE* pFile1;
	pFile1 = fopen ("c:/test/armOut.xml","w");
	fprintf (pFile1, "%s\n", chaineXMLOut);
	fclose (pFile1);
	*/ //<-test

    // Parse xmlResFlows and select the the relevant flow
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
		
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
			                  CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument),
                              (void **) &XMLDoc);

		SUCCEEDED(hr) ? 0 : throw hr;


        wchar_t* wcharStr = constchar2wchar((const char *) xmlResFlows);

        _bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
		   XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString) "Pb in creating XML document for getting flows in RetrieveSummitFlow");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode* theNode = NULL;
	
	long nbNodes;

	CCString tmpChaine;

	tmpChaine = (CCString)"*/*/*/Flows/Flow";

	if ( XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes == 0 )
		{
			hr = S_FALSE;

			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}
			
	        if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}

			CCString msg((CCString)"Invalid XML string for getting flows");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		for (int indexNode = 0; indexNode < nbNodes; indexNode++)
		{
			theNode = NULL;

			hr = resultList->get_item(indexNode, &theNode);

			if ( hr == S_OK && theNode != NULL )
			{
				ARM_Date fixing = GetDateFromXMLNode (theNode, SUMMIT_NAMES::FIXINGDATE);
				
				if ( (inDate.GetJulian() <= fixing.GetJulian()) && (fixing.GetYear() <= 3000) )
				{
					ARM_Currency* ccy = NULL;

					outFlow->fStartDate		= GetDateFromXMLNode (theNode, SUMMIT_NAMES::STARTDATE);
					outFlow->fEndDate		= GetDateFromXMLNode (theNode, SUMMIT_NAMES::ENDDATE);
					outFlow->fFixingDate	= GetDateFromXMLNode (theNode, SUMMIT_NAMES::FIXINGDATE);
					outFlow->fIntDays		= GetIntFromXMLNode  (theNode, SUMMIT_NAMES::INTDAYS);
					outFlow->fRate			= 0.01*GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::RATE);
					outFlow->fFxRate		= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::FXRATE);
					outFlow->fForward		= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::FORWARD);
					outFlow->fSpread		= 0.0001*GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::SPREAD);
					outFlow->fDecompRate	= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::DECOMPRATE);
					outFlow->fNotional		= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::NOTIONAL);
					outFlow->fPayDate		= GetDateFromXMLNode  (theNode, SUMMIT_NAMES::PAYDATE);
					outFlow->fInterimInterest=GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::INTERIMINTEREST);
					outFlow->fFlows			= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::FLOWS);
					outFlow->fStrike		= 0.01*GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::STRIKE);
					outFlow->fType			= GetStringFromXMLNode(theNode, SUMMIT_NAMES::TYPE);					
					outFlow->fDays			= GetIntFromXMLNode(theNode, SUMMIT_NAMES::DAYS);
					ccy						= GetCcy(theNode);
					outFlow->fCcy			= ccy->GetCcyName();
					outFlow->fZeroRate		= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::ZERORATE);
					outFlow->fDiscFactor	= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::DISCFACTOR);
					outFlow->fPV			= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::PV);
					outFlow->fAICRate		= GetDoubleFromXMLNode(theNode, SUMMIT_NAMES::AICRATE);

					delete ccy;

					break;
				}

				if (theNode)
					theNode->Release();
				theNode = NULL;
			}
			else
			{
				if (theNode)
				{
					theNode->Release();
					theNode = NULL;
				}

				if (resultList)
				{
				   resultList->Release();
				
				   resultList = NULL;
				}

		        if (XMLDoc)
				{
					XMLDoc->Release();
					XMLDoc = NULL;
				}

				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							    "RetrieveSummitFlow: Probleme in get_item");
			}
		}

		if ( indexNode > nbNodes )
		{
			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}

		    if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}
			
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "RetrieveSummitFlow: First Reset Date after AsOf not Found");
		}
		else
		{
			if (resultList)
			{
			   resultList->Release();
			
			   resultList = NULL;
			}

		    if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}		
		}
	}
	else
	{
		if (resultList)
		{
		   resultList->Release();
		
		   resultList = NULL;
		}

        if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}
		
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "RetrieveSummitFlow: Parsing problem");
	}
}



int ARMLOCAL_ParseXMLForInterpolator(const char* chaineXML)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

        hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, 
                                       __uuidof(MSXML2::IXMLDOMDocument), (void **) &XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr; 
	}

	catch(...)
	{
		if (XMLDoc)
           XMLDoc->Release();
		
        hr = S_FALSE;

		CCString msg("Pb in creating XML document for Getting Interpolator");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNode * theNode = NULL;

	int interpMethodId;

	try
	{
		CCString tmpChaine = "/Response/Entity/MKTDATA/InterpMethod";
		if (XMLDoc->selectSingleNode(_bstr_t((const char *)(tmpChaine)), &theNode) == S_OK)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff2(resultat,false);
			char * ff3=(char *)ff2;

			if (strcmp(ff3,"CONT") == 0)
				interpMethodId = K_CONTINUOUS;
			else
				interpMethodId = K_LINEAR;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}
		if (XMLDoc) XMLDoc->Release();
		return interpMethodId;
	}
	catch(...)
	{
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for Getting interpolator");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}

ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForZC(const char* chaineXML, 
                                          ARM_Date aSdate, const char* sCcy, int interp)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

        hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, 
                                       __uuidof(MSXML2::IXMLDOMDocument), (void **) &XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr; 
	}

	catch(...)
	{
		if (XMLDoc)
           XMLDoc->Release();
		
        hr = S_FALSE;

		CCString msg("Pb in creating XML document for ZC");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode* listItem = NULL, *theNode = NULL;
	long nbNodes;

	try
	{
		CCString tmpChaine = "/Response/CurveData/Row";

		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if ( nbNodes == 0 )
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for ZC");
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			pdRate = new double[nbNodes];
			pdMatu = new double[nbNodes];

			for (long indexNode = 0 ; indexNode < nbNodes ; indexNode++)
			{
				hr = resultList->get_item(indexNode, &listItem);

				if ( hr == S_OK && listItem != NULL )
				{
					ARM_Date endDate = GetDateFromXMLNode(listItem,"Date");
					pdMatu[indexNode] = (endDate - aSdate) / 365.;
					
					pdRate[indexNode] = GetDoubleFromXMLNode(listItem,"Rate")*100.0;

					listItem->Release();
					listItem=NULL;
				}
			}
		}
		
        if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		
        if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{
		if (pdRate)
			delete[] pdRate;
		
        if (pdMatu)
			delete[] pdMatu;
		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for ZC");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Vector* mktMatu = NULL;

	ARM_Currency* aCcy = NULL;

	ARM_ZeroLInterpol* newZcLin = NULL;

	try
	{
		mktData = new ARM_Vector(nbNodes, pdRate);
		mktMatu = new ARM_Vector(nbNodes, pdMatu);

		newZcLin = new ARM_ZeroLInterpol(aSdate, mktMatu, mktData, 0, 0, interp);

		aCcy = new ARM_Currency((char *) sCcy);

		newZcLin->SetCurrencyUnit(aCcy);

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		delete[] pdRate;
		delete[] pdMatu;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (mktMatu)
			delete mktMatu;
		mktMatu = NULL;

		return newZcLin;
	}

	catch(...)
	{
		delete[] pdRate;
		delete[] pdMatu;

		if (mktData)
			delete mktData;
		mktData = NULL;

		if (mktMatu)
			delete mktMatu;
		mktMatu = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		hr = S_FALSE;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
               "Error in getting ZC Curve");
	}
}



ARM_ZeroLInterpol* ARMLOCAL_ParseXMLForMYAndZC(const char* chaineXML,
                                               ARM_Date aSdate, const char* sCcy,
											   long adjOrNotId,
											   long rawId,
											   long swapFrqId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	VECTOR<double> pdRate;
	VECTOR<CCString> pdMatu;

	try
	{
		hr = CoInitialize(NULL); 
		
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void **) &XMLDoc);
		
        SUCCEEDED(hr) ? 0 : throw hr;


        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();
		
        hr = S_FALSE;

		CCString msg((CCString)"Error in Loading XML string for MY and ZC");

		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						(const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNodeList* resultList2 = NULL;
	MSXML2::IXMLDOMNodeList* resultListSegment = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * listItemSegment = NULL;
	MSXML2::IXMLDOMNode * listItem2 = NULL;
	long nbNodes;
	long nbNodesSegment;

	int priorityMMFut;
	int priorityFutSwap;
	int interpMethodId;
	int methodId;

	long isBS = 0;

	try
	{
		CCString tmpChaine = "/Response/Entity/MKTDATA";
		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for MY and ZC");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			if (listItem->selectSingleNode(_bstr_t((const char *)"MktSpecs/MKTSPEC/MktType"), &theNode) == S_OK)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"SPR") == 0)
					isBS = 1;
				else
					isBS = 0;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (hr==S_OK && listItem!=NULL)
			{
				BSTR resultat = NULL;

				if (isBS == 0)
				{
					hr=listItem->selectSingleNode(_bstr_t((const char *)"GenOptions"), &theNode);

					// Récupération des paramètres de priorité
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					int myOptions = atoi(ff1);

					int maskMM = 020;
					int maskFut = 002;

					if ( (myOptions & maskMM) != 0)
						priorityMMFut = K_MM;
					else if ( (myOptions & maskFut) != 0)
						priorityMMFut = K_FUT;
					else
					{
						CCString msg ((CCString)"Invalid priority bentween MM and Fut");

						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
										 (const char *) msg);
					}

					int maskSSwap = 004;

					if ( (myOptions & maskSSwap) != 0)
						priorityFutSwap = K_SWAP;
					else
						priorityFutSwap = K_FUT;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				// Récupération de la méthode d'interpolation
				hr = listItem->selectSingleNode(_bstr_t((const char *)"InterpMethod"), &theNode);

				resultat = NULL;
				theNode->get_text(&resultat);

			    _bstr_t ff2(resultat,false);
				char * ff3=(char *)ff2;

				if (strcmp(ff3,"CONT") == 0)
					interpMethodId = K_CONTINUOUS;
				else
					interpMethodId = K_LINEAR;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);

				// Récupération de la méthode de génération
				hr = listItem->selectSingleNode(_bstr_t((const char *)"Method"), &theNode);

				resultat = NULL;
				theNode->get_text(&resultat);

			    _bstr_t ff4(resultat,false);
				char * ff5=(char *)ff4;

				if (strcmp(ff5,"RAW") == 0)
					methodId = K_RAW;
				else if ( (strcmp(ff5,"FWD") == 0) || (strcmp(ff5,"cPATY") == 0) )
					methodId = K_FORWARD;
				else
					methodId = K_PAR;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			if (listItem->selectNodes(_bstr_t((const char *)"MktSpecs/MKTSPEC"), &resultList2) == S_OK)
			{
				resultList2->get_length(&nbNodes);

				if (nbNodes == 0)
				{
					hr = S_FALSE;

					CCString msg ((CCString)"Invalid XML string for MY and ZC");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}
					
				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList2->get_item(indexNode, &listItem2);

					if (listItem2->selectNodes(_bstr_t((const char *)"MktPoints/MKTPOINT"), &resultListSegment) == S_OK)
					{
						resultListSegment->get_length(&nbNodesSegment);

						for (long indexNodeSegment=0; indexNodeSegment<nbNodesSegment;indexNodeSegment++)
						{
							hr=resultListSegment->get_item(indexNodeSegment, &listItemSegment);
							if (hr==S_OK && listItemSegment!=NULL)
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);
									
									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									pdMatu.push_back(ff1);

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
								if (isBS == 0)
									listItemSegment->selectSingleNode(_bstr_t((const char *)"MidRate"),&theNode);
								else
									listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);

								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									double rate = atof(ff1);
									if (isBS == 0)
										rate *= 100.;	//printf("rate : %.13lf\n", rate);
									else
									{
										// Les spreads sont en bp
										// Contrairement au GetInitial(), ici on multiplie par -1
										rate *= -10000.;
									}

									pdRate.push_back(rate);
									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}

								if ( (isBS == 0) && (adjOrNotId == K_YES) )
								{
									listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);

									if (theNode!=NULL)
									{
										BSTR resultat = NULL;
										theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

										_bstr_t ff(resultat,false);
										char * ff1=(char *)ff;

										double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);

										pdRate[pdRate.size()-1] += rate * 100.;

										theNode->Release();
										theNode=NULL;
										if (resultat) SysFreeString(resultat);
									}
								}

								listItemSegment->Release();
								listItemSegment=NULL;
							}
						}
					}
					resultListSegment->Release();
					resultListSegment = NULL;

					listItem2->Release();
					listItem2=NULL;
				}
				if (resultList2)
				{
					resultList2->Release();
					resultList2 = NULL;
				}
			}

			if (listItemSegment)
			{
				listItemSegment->Release();
				listItemSegment = NULL;
			}

			if (listItem)
			{
				listItem->Release();
				listItem = NULL;
			}

			if (theNode)
			{
				theNode->Release();
				theNode = NULL;
			}

			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}

			if (resultListSegment)
			{
				resultListSegment->Release();
				resultListSegment = NULL;
			}

			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}

			hr = S_OK;
		}
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Invalid XML string for MY and ZC");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	ARM_Vector* mktData = NULL;
	ARM_Currency* aCcy = NULL;
	double* dRate = NULL;

	ARM_ZeroLInterpol* createdZcLin = NULL;
	
	try
	{
		char tmp[20];
		
        strcpy(tmp, sCcy);
		
	    aCcy = new ARM_Currency(tmp);


		dRate = new double[pdRate.size()];
		
        for (int i = 0; i < pdRate.size(); i++)
		{
			dRate[i] = pdRate[i];
		}

		char psMatu[ARM_NB_TERMS][ARM_NB_MAX_CHAR_TERMS];
		
        for (int j = 0; j < ARM_NB_TERMS; j++)
        {
			strcpy(psMatu[j], "X");
        }

		for (j = 0; j < pdRate.size(); j++)
		{
		    strcpy(psMatu[j], (const char *) pdMatu[j]);
		}

		mktData = new ARM_Vector(pdRate.size(),dRate);

		if ( isBS == 0 )
		{
		   if ( rawId != K_DEFAULT_CURVEMOD )
			  methodId = rawId;

			createdZcLin = new ARM_ZeroLInterpol(aSdate, psMatu,
												 mktData, priorityMMFut, priorityFutSwap, 
												 methodId, interpMethodId, 0, aCcy, swapFrqId);
		}
		else
		{
			createdZcLin = new ARM_ZeroLInterpol(aSdate, mktData, psMatu,
												 -1, 0, interpMethodId, aCcy);
		}

		if (mktData)
		   delete mktData;
		mktData = NULL;

		if (aCcy)
		   delete aCcy;
		aCcy = NULL;

		if (dRate)
		   delete [] dRate;
		dRate = NULL;

	}

    catch(...)
	{
		if (mktData)
			delete mktData;
		mktData = NULL;

		if (aCcy)
			delete aCcy;
		aCcy = NULL;

		if (dRate)
			delete [] dRate;
		dRate = NULL;

		hr = S_FALSE;

		CCString msg ((CCString)"Pb in Creating ZC");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return createdZcLin;
}


ARM_InfCurv* ARMLOCAL_ParseXMLForInfZC(const char* chaineXML,
									   double aSdate,
									   const CCString& index)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 
		
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void **) &XMLDoc);
		
        SUCCEEDED(hr) ? 0 : throw hr;


        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();
		
        hr = S_FALSE;

		CCString msg((CCString)"Error in Loading XML string for Inflation ZC");

		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						(const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNodeList* resultListSegment = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * listItemSegment = NULL;
	long nbNodes;
	long nbNodesSegment;

	char buffer[50];

	vector<string> MktTerms;
	vector<double> MktValues;

	ARM_Date CPIIndexDate;
	double CPIIndex;

	ARM_InfCurv* infCurv = NULL;

	try
	{
		CCString tmpChaine = "/Response/Entity/MKTDATA";
		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for Inflation ZC");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			if (listItem->selectNodes(_bstr_t((const char *)"MktSpecs/MKTSPEC"), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
					hr = S_FALSE;

					CCString msg ((CCString)"Invalid XML string for MY and ZC");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}
					
				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{
					hr=resultList->get_item(indexNode, &listItem);

					if (listItem->selectNodes(_bstr_t((const char *)"MktPoints/MKTPOINT"), &resultListSegment) == S_OK)
					{
						resultListSegment->get_length(&nbNodesSegment);

						for (long indexNodeSegment=0; indexNodeSegment<nbNodesSegment;indexNodeSegment++)
						{
							hr=resultListSegment->get_item(indexNodeSegment, &listItemSegment);
							if (hr==S_OK && listItemSegment!=NULL)
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);
									
									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									if (strlen((const char*)ff1) == 8)
									{
										ARM_Date tmpDate(ff1,"YYYYMMDD");
										sprintf( buffer, "%f", tmpDate.GetJulian() );
										MktTerms.push_back(CCSTringToSTLString(buffer));

										if ( (indexNodeSegment == 0) && (indexNode == 0) )
											CPIIndexDate = tmpDate;
									}
									else
									{
										MktTerms.push_back(CCSTringToSTLString((const char*)ff1));
									}

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
								listItemSegment->selectSingleNode(_bstr_t((const char *)"MidRate"),&theNode);

								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									double rate = atof(ff1)*100.;
									MktValues.push_back(rate);
									if ( (indexNodeSegment == 0) && (indexNode == 0) )
										CPIIndex = rate;

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}

								listItemSegment->Release();
								listItemSegment=NULL;
							}
						}
					}
					resultListSegment->Release();
					resultListSegment = NULL;

					listItem->Release();
					listItem=NULL;
				}
			}

			if (listItemSegment)
			{
				listItemSegment->Release();
				listItemSegment = NULL;
			}

			if (listItem)
			{
				listItem->Release();
				listItem = NULL;
			}

			if (theNode)
			{
				theNode->Release();
				theNode = NULL;
			}

			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}

			if (resultListSegment)
			{
				resultListSegment->Release();
				resultListSegment = NULL;
			}

			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}

			hr = S_OK;
		}
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Invalid XML string for Inflation ZC");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	char sAsOfDate[11];
	Local_XLDATE2ARMDATE(aSdate,sAsOfDate);
	ARM_Date asOfARMDate( sAsOfDate );

	infCurv = new ARM_InfCurv(asOfARMDate,
							  CCSTringToSTLString(index),
							  CPIIndex,
							  CPIIndexDate,
							  MktTerms,
							  MktValues);

	return infCurv;
}





ARM_ZeroLInterpol* ARMLOCAL_CreateZC(const CCString& index,
									 const CCString& currency,
									 const CCString& cvName,
									 ARM_Date aSdate,
									 const CCString& raw,
									 long adjOrNotId,
									 long swapFrqId)
{
	CCString xmlResponse = etoolkit_getXMLMYAndZCFromSummit(index,currency,cvName,aSdate);

	long rawId = ARM_ConvCvMethod (raw);

    ARM_ZeroLInterpol* createdZcLin = ARMLOCAL_ParseXMLForMYAndZC(xmlResponse, aSdate, 
                                                                  currency, adjOrNotId,
                                                                  rawId, swapFrqId);
    string	vType("ZERO");

    if( index == "BS" )
    {
		vType = string("BASIS USD");
    }

	string	vIndex((const char *) index);
	string	vCurrency((const char *) currency);
	string	vCrvId((const char *) cvName);
    createdZcLin->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);												

	return(createdZcLin);
}




ARM_BasisCurve* ARMLOCAL_CreateZCSpreaded(const CCString& index,
										  const CCString& currency,
										  const CCString& cvName,
										  ARM_Date aSdate,
										  long adjOrNotId,
										  const CCString& raw,
										  long swapFrqId,
										  long mmFrqId,
										  long interpId,
										  ARM_ZeroCurve* ZeroCurve)
{
	ARM_ZeroLInterpol*	zc = NULL;
	ARM_BasisCurve*		newZc = NULL;
	ARM_ZeroCurve*		ZCSpread = NULL;
	ARM_ZeroCurve*		ZCInit = NULL;
	ARM_Currency*		ccy = NULL;

    try
    {
		long	vSummitSwapFreq = -1;

		ZCSpread = ARMLOCAL_CreateZC("BS", currency, cvName, aSdate, raw, adjOrNotId, vSummitSwapFreq);
		if( ZCSpread == NULL )
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Can not find BS curve");
		}

        // Update Mkt data characteristics
		string	vType("BASIS USD");
		string	vIndex((const char*) index);
		string	vCurrency((const char*) currency);
		string	vCrvId((const char*) cvName);

        ZCSpread->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

		if (ZeroCurve)
			ZCInit = ZeroCurve;
		else
		{
			ZCInit = ARMLOCAL_CreateZC(index, currency, cvName, aSdate, raw, adjOrNotId, vSummitSwapFreq);
		
			if( ZCInit == NULL )
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Can not find ZC curve");
			}
		}
		
		if( currency == "DEFAULT" )
			ccy = ARM_DEFAULT_CURRENCY;
		else
			ccy = new ARM_Currency((const char*) currency);

		// if no frequencies input, depend on currency
		if (mmFrqId == K_DEF_FREQ)
			mmFrqId = K_MONTHLY;

		if (swapFrqId == K_DEF_FREQ)
		{
			if ((vCurrency == "EUR") || (vCurrency == "SEK"))
				swapFrqId = K_ANNUAL;
			else
				swapFrqId = K_SEMIANNUAL;
		}

        newZc = new ARM_BasisCurve(	aSdate, ZCSpread, ZCInit, mmFrqId, swapFrqId, ccy);

		return newZc;
    }
	catch (Exception& e)
	{
		throw e;
	}
    catch(...)
    {
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Error in creating ZC spreaded");
    }
}



long ARMLOCAL_ParseXMLForMY(const char* chaineXML,
							long adjOrNotId,
							VECTOR<CCString>* matu,
							VECTOR<double>* yield)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine ; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	
    catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for MY");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultListSegment = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * theNode1 = NULL, * theNode2 = NULL, * listItemSegment = NULL;
	long nbNodes;
	long nbNodesSegment;
	long IsSpreadCurve = 0;

	char segment[20];

	try
	{
		string sCcy;
		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/Entity/MKTDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;
				CCString msg((CCString)"Invalid XML string for MY");
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,(const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);
			listItem->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				sCcy = ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->Release();
			listItem=NULL;
		}

		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/Entity/MKTDATA/MktSpecs/MKTSPEC"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for MY");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"MktType"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					strcpy(segment,ff1);

					if (strcmp(ff1,"SPR") == 0)
						IsSpreadCurve = 1;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				if (listItem->selectNodes(_bstr_t((const char *)"MktPoints/MKTPOINT"), &resultListSegment) == S_OK)
				{
					resultListSegment->get_length(&nbNodesSegment);

					double resVal (0.);

					for (long indexNodeSegment=0; indexNodeSegment<nbNodesSegment;indexNodeSegment++)
					{
						hr=resultListSegment->get_item(indexNodeSegment, &listItemSegment);
						if (hr==S_OK && listItemSegment!=NULL)
						{
							if (strcmp(segment,"FWD") == 0)
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode1);
								listItemSegment->selectSingleNode(_bstr_t((const char *)"ToDate"), &theNode2);

								if ((theNode1!=NULL) && (theNode2!=NULL))
								{
									BSTR resultat = NULL;
									theNode1->get_text(&resultat);
									_bstr_t f1(resultat,false);
									char * ff1=(char *)f1;

									string day1 = string(ff1).substr(6, 2);
									string month1 = string(ff1).substr(4, 2);
									string year1 = string(ff1).substr(2, 2);

									theNode2->get_text(&resultat);
									_bstr_t f2(resultat,false);
									char * ff2=(char *)f2;

									string day2 = string(ff2).substr(6, 2);
									string month2 = string(ff2).substr(4, 2);
									string year2 = string(ff2).substr(2, 2);

									char str[22];
									if (strcmp(sCcy.c_str(), ARM_DEFAULT_CURRENCY->GetCcyName()) == 0)
										sprintf(str, "%s/%s/%s-%s/%s/%s", day1.c_str(), month1.c_str(), year1.c_str(), day2.c_str(), month2.c_str(), year2.c_str());
									else
										sprintf(str, "%s/%s/%s-%s/%s/%s", month1.c_str(), day1.c_str(), year1.c_str(), month2.c_str(), day2.c_str(), year2.c_str());
									
									matu->push_back(str);

									theNode1->Release();
									theNode1=NULL;
									theNode2->Release();
									theNode2=NULL;
									if (resultat) SysFreeString(resultat);
								}
							}
							else if (strcmp(segment,"PRI") == 0)
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);

								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);
									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									string day = string(ff1).substr(6, 2);
									string month = string(ff1).substr(4, 2);
									string year = string(ff1).substr(2, 2);

									char str[22];
									if (strcmp(sCcy.c_str(), ARM_DEFAULT_CURRENCY->GetCcyName()) == 0)
										sprintf(str, "%s/%s/%s", day.c_str(), month.c_str(), year.c_str());
									else
										sprintf(str, "%s/%s/%s", month.c_str(), day.c_str(), year.c_str());
									
									matu->push_back(str);

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
							}
							else
							{
								if (strcmp(segment,"BOND") == 0)
								{
									listItemSegment->selectSingleNode(_bstr_t((const char *)"Sec"), &theNode);
								}
								else
								{
									listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
								}
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);
									
									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									matu->push_back(ff1);

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
							}
							if (IsSpreadCurve == 0)
								listItemSegment->selectSingleNode(_bstr_t((const char *)"MidRate"),&theNode);
							else
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);
							
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);
								if (IsSpreadCurve == 1)
									resVal = rate * 10000.;
								else if ( (strcmp(segment,"FXS") == 0) || (strcmp(segment,"FXF") == 0) )
									resVal = rate;
								else
									resVal = rate * 100.;

								theNode->Release();
								theNode=NULL;
								if (resultat) SysFreeString(resultat);
							}

							if ( (IsSpreadCurve == 0) && (adjOrNotId == K_YES) )
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);

								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);

									resVal += rate * 100.;

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
							}

							yield->push_back(resVal);

							listItemSegment->Release();
							listItemSegment=NULL;
						}
					}
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItemSegment)
		{
			listItemSegment->Release();
			listItemSegment = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (resultListSegment)
		{
			resultListSegment->Release();
			resultListSegment = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for MY");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}


void ARMLOCAL_ParseXMLForMYCalypso(const std::string& chaineXML,
							bool doAdjust,
							std::vector<std::string>&matus ,
							std::vector<double>& yields)
{
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 

	xmlDoc=XMLTools::CreateDOMDocument30(); 
									   
	_bstr_t b;
	VariantTools::convert(chaineXML,b); 
	xmlDoc->loadXML(b,&bOK); 

	#ifdef _DEBUG
	_variant_t v_chtmp; 
	v_chtmp.SetString("c:\\temp\\Calypso_zerocurve.xml"); 
	// VariantTools::convert("c:\\temp\\Calypso_zerocurve.xml",v_chtmp); 
	xmlDoc->save(v_chtmp);
	#endif


	MSXML2::IXMLDOMNodeListPtr resultList; 
	resultList = XMLTools::selectNodes(xmlDoc,"/ZeroCurveList/CurveZero/Underlyings/Underlying") ;
	long nbNodes; 
	resultList->get_length(&nbNodes);
	if (nbNodes==0) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"ARMLOCAL_ParseXMLForMYCalypso: no underlyings found"); 
	
	matus.resize(nbNodes); 
	yields.resize(nbNodes); 

	for(long i=0;i<nbNodes;i++) 
	{
		MSXML2::IXMLDOMNodePtr item = XMLTools::get_item(resultList,i); 
		MSXML2::IXMLDOMNodePtr item2 = XMLTools::selectSingleNode(item,"Maturity"); 
		XMLTools::convert(item2,matus[i]); 
		item2 = XMLTools::selectSingleNode(item,"Value"); 
		XMLTools::convert(item2,yields[i]);  
		std::string type ; 
		item2 = XMLTools::selectSingleNode(item,"Type"); 
		XMLTools::convert(item2,type); 
		if (type=="MM" || type=="Swap") yields[i]*=100; 
		if (doAdjust)
		{
			item2 = XMLTools::selectSingleNode(item,"FutureConvexity"); 
			double value; XMLTools::convert(item2,value); 
			yields[i] += value*100; 
		}
	}

/**
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		_bstr_t tmpChaine = chaineXML;

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	
    catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for MY");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNodeList * resultListSegment = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * listItemSegment = NULL;
	long nbNodes;
	long nbNodesSegment;
	long IsSpreadCurve = 0;

	char segment[20];

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/SPECS/Spec"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for MY");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"SpecType"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					strcpy(segment,ff1);

					if (strcmp(ff1,"SPR") == 0)
						IsSpreadCurve = 1;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				if (listItem->selectNodes(_bstr_t((const char *)"MKTPOINTS/MktPoint"), &resultListSegment) == S_OK)
				{
					resultListSegment->get_length(&nbNodesSegment);

					double resVal (0.);

					for (long indexNodeSegment=0; indexNodeSegment<nbNodesSegment;indexNodeSegment++)
					{
						hr=resultListSegment->get_item(indexNodeSegment, &listItemSegment);
						if (hr==S_OK && listItemSegment!=NULL)
						{
							if (strcmp(segment,"BOND") == 0)
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"SecId"), &theNode);
							}
							else
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
							}
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);
								
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								matu->push_back(ff1);

								theNode->Release();
								theNode=NULL;
								if (resultat) SysFreeString(resultat);
							}

							if (IsSpreadCurve == 0)
								listItemSegment->selectSingleNode(_bstr_t((const char *)"MidRate"),&theNode);
							else
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);
							
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);
								if (IsSpreadCurve == 1)
									resVal = rate * 10000.;
								else
									resVal = rate * 100.;
								
								theNode->Release();
								theNode=NULL;
								if (resultat) SysFreeString(resultat);
							}

							if ( (IsSpreadCurve == 0) && (adjOrNotId == K_YES) )
							{
								listItemSegment->selectSingleNode(_bstr_t((const char *)"Spread"),&theNode);

								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);	//printf("resultat : %S\n", resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									double rate = atof(ff1);	//printf("rate : %.13lf\n", rate);

									resVal += rate * 100.;

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
							}

							yield->push_back(resVal);

							listItemSegment->Release();
							listItemSegment=NULL;
						}
					}
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItemSegment)
		{
			listItemSegment->Release();
			listItemSegment = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (resultListSegment)
		{
			resultListSegment->Release();
			resultListSegment = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for MY");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK; **/ 
}



VECTOR<CCString> ARMLOCAL_GetListTenorsFromXML(const char* chaineXML, VECTOR<CCString>* listRequetes, CCString NomParam)
{
	VECTOR<CCString> res;

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Tenors for Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/Commdata/Name"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Tenors for Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)NomParam), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					res.push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->selectSingleNode(_bstr_t((const char *)"Val"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);
					
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					listRequetes->push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Tenors for Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return res;
}


VECTOR<CCString> ARMLOCAL_GetListTenorsForVol(const char* chaineXML,bool isSmile)
{
	VECTOR<CCString> res;

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Tenors for Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/COMMSET_LIST/COMMSET"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Tenors for Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				if (isSmile)
					listItem->selectSingleNode(_bstr_t("Name2"), &theNode);
				else
					listItem->selectSingleNode(_bstr_t("Date"), &theNode);

				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					res.push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Tenors for Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return res;
}



VECTOR<CCString> ARMLOCAL_GetListStrikesFromXML(const char* chaineXML, VECTOR<CCString>* listRequetes)
{
	VECTOR<CCString> res;

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Strikes for Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)"Response/Commdata/Name"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Strikes for Smile");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*)msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Val"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					listRequetes->push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/Commdata/Name2"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Strikes for Smile");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Val"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					res.push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Strikes for Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return res;
}


// Fonction de conversion des chaines (2D, 1M, Contrats Futur) en fraction d'année
double convPlotInYearTerm(const CCString& inPlot, ARM_Date AsOf, ARM_Date Settle, 
                          char* calenDar)
{
	long isDate = 0;
	long Nb;
	char cMatu;
	long freqId;

    char plot[50];

    char* currency = calenDar;

    strcpy(plot, (const char *) inPlot);

	// Contrat
	if ( strlen((const char *) plot) == 5 )
	{
		int month, year;
		ARM_Date matDate;

		GetMonthYearFromExpiryDate(plot, &month, &year);
		matDate.ChangeDate(1, month, year);

		matDate.PiborDelivery();
		return ( (matDate.GetJulian() - AsOf.GetJulian()) /365.);
	}
	else
	{
		if ( strlen((const char *) plot) >= 8 )
		{
			isDate = 1;
		}
		else
		{
			sscanf(plot, "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
				isDate = 1;
		}
		
		if ( isDate == 1 )
		{
			if ( (IsETKVersion()) || (GetDataRetrieverVersion () == WSETKRETRIEVER) )
			{
				if(plot[2] == '/')
				{
					ARM_Date tmpDate(plot, "DD/MM/YYYY");

					return ((tmpDate-AsOf)/365.);
				}
				else
				{
					ARM_Date tmpDate(plot, "YYYYMMDD");

					return ((tmpDate-AsOf)/365.);
				}
			}
			else
			{
				if ( strcmp(ARM_DEFAULT_COUNTRY, "USD") == 0 )
				{
					ARM_Date tmpDate(plot, "MM/DD/YYYY");
										
                    return (tmpDate - AsOf) /365.;
				}
				else
				{
					ARM_Date tmpDate(plot);
					
                    return (tmpDate - AsOf) /365.;
				}
			}
		}
		else
		{
			if ( freqId == K_DAILY )
			{
				Settle.AddDays(Nb);
				Settle.AdjustToBusDate(currency, K_FOLLOWING);
			}
			else
			{
				Settle.AddPeriodMult(freqId, (long) Nb, currency);
				Settle.AdjustToBusDate(currency, K_FOLLOWING);
			}

			return ( (Settle - AsOf) /365.);
		}
	}	
}




long ARMLOCAL_GetDatesAndDatasForVol(const char* chaineXML,
									 ARM_Date asof,
									 long index,
									 ARM_Matrix*& Volatility,
									 ARM_Vector*& YearTerms,
									 const char* currency,
									 long nbStrikes,
									 VECTOR<CCString>* listYearTerms)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	ARM_Matrix* tmpVol;
	ARM_Currency sCCY(currency);

	long spotDays = sCCY.GetSpotDays();

    char* payCalTmp = sCCY.GetPayCalName(sCCY.GetVanillaIndexType());

    char payCal[30];

    strcpy(payCal, payCalTmp);
    delete payCalTmp;

	ARM_Date settleDate = asof;
	settleDate.NextBusinessDay(spotDays, payCal);

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();

        hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode *listItem = NULL, *theNode = NULL;
	
    long nbNodes;

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/Entity/COMMSET/CommData/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			if ( Volatility == NULL )
			{
			   YearTerms = new ARM_Vector(nbNodes);
				
               tmpVol = new ARM_Matrix(nbNodes,nbStrikes);
			}

			for (long indexNode = 0; indexNode<nbNodes; indexNode++)
			{
				hr = resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				
                if ( theNode != NULL )
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char* ff1 = (char *) ff;

					double vol = atof(ff1);
					if ( Volatility == NULL )
					   tmpVol->Elt(indexNode,index) = vol * 100.;
					else
					   Volatility->Elt(indexNode,index) = vol * 100.;

					theNode->Release();
					theNode = NULL;
					
                    if (resultat) 
                       SysFreeString(resultat);
				}

				if ( Volatility == NULL )
				{
					listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
					
                    if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						YearTerms->Elt(indexNode) = convPlotInYearTerm(ff1, asof, settleDate, 
																			payCal);

                        if (listYearTerms)
						   listYearTerms->push_back(ff1);

						theNode->Release();
						theNode = NULL;
						
                        if (resultat) 
                           SysFreeString(resultat);
					}
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		if (Volatility == NULL)
		{
			Volatility = tmpVol;
		}

		return ARM_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}



long ARMLOCAL_NewGetDatesAndDatasForVol(const char* chaineXML,
										ARM_Date asof,
										ARM_Matrix*& Volatility,
										ARM_Vector*& YearTerms,
										const char* currency,
										VECTOR<CCString>* listYearTerms)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	ARM_Currency sCCY(currency);

	long spotDays = sCCY.GetSpotDays();

    char* payCalTmp = sCCY.GetPayCalName(sCCY.GetVanillaIndexType());

    char payCal[30];

    strcpy(payCal, payCalTmp);
    delete payCalTmp;

	ARM_Date settleDate = asof;
	settleDate.NextBusinessDay(spotDays, payCal);

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();

        hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode *listItem = NULL, *listItem2 = NULL, *theNode = NULL;
	
    long nbTenors;
    long nbExpiries;

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/COMMSET_LIST/COMMSET"), &resultList) == S_OK)
		{
			resultList->get_length(&nbTenors);

			if (nbTenors == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode = 0; indexNode<nbTenors; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				if (listItem->selectNodes((_bstr_t)((const char *)"COMMDATA"), &resultList2) == S_OK)
				{
					resultList2->get_length(&nbExpiries);

					if (nbExpiries == 0)
					{
						hr = S_FALSE;

						CCString msg((CCString)"Invalid XML string for getting Volatility");

						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
										 (const char *) msg);
					}

					if ( Volatility == NULL )
					{
					   YearTerms = new ARM_Vector(nbExpiries);
						
					   Volatility = new ARM_Matrix(nbExpiries,nbTenors);
					}

					for (long indexNode2 = 0; indexNode2 < nbExpiries; indexNode2++)
					{
						hr=resultList2->get_item(indexNode2, &listItem2);
						
						if (indexNode == 0)
						{
							listItem2->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
							
							if ( theNode != NULL )
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);

								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								YearTerms->Elt(indexNode2) = convPlotInYearTerm(ff1, asof, settleDate, 
																					payCal);

								if (listYearTerms)
								   listYearTerms->push_back(ff1);

								theNode->Release();
								theNode = NULL;
								
								if (resultat) 
								   SysFreeString(resultat);
							}
						}

						listItem2->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);

						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char* ff1 = (char *) ff;

							double vol = atof(ff1);

							Volatility->Elt(indexNode2,indexNode) = vol * 100.;

							theNode->Release();
							theNode = NULL;
							
							if (resultat) 
							   SysFreeString(resultat);
						}
						listItem2->Release();
						listItem2=NULL;
					}
					resultList2->Release();
					resultList2=NULL;
				}
				listItem->Release();
				listItem=NULL;
			}
			resultList->Release();
			resultList=NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return ARM_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}


long ARMLOCAL_GetMonthsAndDatasForSeasonMgr(const char* chaineXML,
											ARM_GP_Vector& seasonSpreadList,
											vector<string>& monthList)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();

        hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode *listItem = NULL, *theNode = NULL;
	
    long nbNodes;

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/Entity/COMMSET/CommData/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			bool alimMonthList(false);

			if ( monthList.size() == NULL )
			{
				alimMonthList = true;
			}

			for (long indexNode = 0; indexNode<nbNodes; indexNode++)
			{
				hr = resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				
                if ( theNode != NULL )
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char* ff1 = (char *) ff;

					double vol = atof(ff1);
					seasonSpreadList.push_back(vol * 100.);

					theNode->Release();
					theNode = NULL;
					
                    if (resultat) 
                       SysFreeString(resultat);
				}

				if ( alimMonthList)
				{
					listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
					
                    if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						int NumMonth = atoi(ff1);

						monthList.push_back(ShortMonthsCapitalLetter[NumMonth-1]);

						theNode->Release();
						theNode = NULL;
						
                        if (resultat) 
                           SysFreeString(resultat);
					}
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return ARM_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}



long ARMLOCAL_GetStringAndDatasForVol(const char* chaineXML,
									  ARM_Date asof,
									  long index,
									  ARM_Matrix*& Volatility,
									  VECTOR<CCString> &YearTerms,
									  long nbStrikes)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Matrix* tmpVol;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/Entity/COMMSET/CommData/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			if (Volatility == NULL)
			{
				tmpVol = new ARM_Matrix(nbNodes,nbStrikes);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					double vol = atof(ff1);
					if (Volatility == NULL)
						tmpVol->Elt(indexNode,index) = vol * 100.;
					else
						Volatility->Elt(indexNode,index) = vol * 100.;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				if (Volatility == NULL)
				{
					listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						YearTerms.push_back((const char*) ff1);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					listItem->Release();
					listItem=NULL;
				}
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		if (Volatility == NULL)
		{
			Volatility = tmpVol;
		}

		return ARM_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}



long ARMLOCAL_GetDatesAndDatasForCurve(const char* chaineXML,
									   const CCString& currency,
									   ARM_Date asof,
									   long index,
									   ARM_Matrix* &Volatility,
									   ARM_Vector* &vYearTerms,
									   long nbStrikes,
									   long isSmile,
									   VECTOR<CCString>* listYearTerms,
									   const CCString& impOrHist)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	CCString cal;
	CCString xmlResponse;
	CCString msgList;

    ARM_Currency sCCY((const char*) currency);

    char* payCalTmp = sCCY.GetPayCalName(sCCY.GetVanillaIndexType());

    char payCal[30];

    strcpy(payCal, payCalTmp);
    delete [] payCalTmp;

	ARM_Matrix* tmpVol;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	ARM_Date settleDate;

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/Entity/COMMSET/CommData/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Smile");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			if (Volatility == NULL)
			{
				vYearTerms = new ARM_Vector(nbNodes);
				tmpVol = new ARM_Matrix(nbNodes,nbStrikes);
				
				long spotDays = sCCY.GetSpotDays();
			
                settleDate = asof;
			
				settleDate.NextBusinessDay(spotDays, payCal);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					double vol = atof(ff1);
					
					if (isSmile == 1)
						vol *= 100.;

					if (Volatility == NULL)
						tmpVol->Elt(indexNode,index) = vol;
					else
						Volatility->Elt(indexNode,index) = vol;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				if (Volatility == NULL)
				{
					listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (listYearTerms)
							listYearTerms->push_back(ff1);

						if (strlen((const char*) ff1) < 8)
						{

							if (impOrHist=="AIMM")
							{
								// ARM_Date date = AddPeriod(asof,ff1,&sCCY,false,qCredit_Adjust20).GetJulian();
								ARM_Date date = AddPeriod(asof,ff1,sCCY.GetCcyName(),false,qCredit_Adjust20).GetJulian();
								vYearTerms->Elt(indexNode) = (date - asof)/365.;
							}
							else
							vYearTerms->Elt(indexNode) = convPlotInYearTerm(ff1, asof, 
                                                                            settleDate, 
                                                                            payCal);
						}
						else
						{
							ARM_Date tmpDate(ff1,"YYYYMMDD");
							vYearTerms->Elt(indexNode) = (tmpDate - asof) /365.;
						}

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		if (Volatility == NULL)
		{
			Volatility = tmpVol;
		}

		return ARM_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_OK;
}




CCString ARMLOCAL_GetCcyCalFromXML(const char* chaineXML)
{
	char res[20];

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Summit Calendar");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNode* theNode = NULL;

	try
	{
		if (XMLDoc->selectSingleNode((_bstr_t)((const char *)"Response/Cal"), &theNode) == S_OK)
		{
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(res,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return CCString(res);
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Summit Calendar");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


CCString ARMLOCAL_ParseXMLForGetDateLastWarm()
{
	CCString response;
	
	long retCode = etoolkit_getlastdatewarm(response);
	
	char res[20];

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar((const char*)response);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr; 
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Last Warm Date");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNode* theNode = NULL;

	try
	{
		if (XMLDoc->selectSingleNode((_bstr_t)((const char *)"Response"), &theNode) == S_OK)
		{
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(res,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return CCString(res);
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting LastWarmDate");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


ARM_Date ARMLOCAL_GetFwdDate(const char* chaineXML)
{
	char res[20];

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Date resDate;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Fwd Date");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNode* theNode = NULL;

	try
	{
		if (XMLDoc->selectSingleNode((_bstr_t)("Response/FwdDate"), &theNode) == S_OK)
		{
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(res,ff1);

				ARM_Date tmpDate(res,"YYYYMMDD");
				resDate = tmpDate;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return resDate;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Summit Calendar");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



double ARMLOCAL_ParseFxCurve(const char* chaineXML)
{
	double val = 0.0;

	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * theNode2 = NULL;
	long nbNodes;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting FX Rate");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)("Response/COMMSET_LIST/COMMSET/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FX Rate");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					const char * ff1=(const char *)ff;

					if (strcmp(ff1,"SPOT") == 0)
					{
						listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode2);
						if (theNode2!=NULL)
						{
							BSTR resultat2 = NULL;
							theNode2->get_text(&resultat2);

							_bstr_t ff(resultat2,false);
							char * ff1=(char *)ff;

							val = atof(ff1)*100.;

							theNode2->Release();
							theNode2=NULL;
							if (resultat2) SysFreeString(resultat2);

							indexNode = nbNodes;
						}
					}
					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (theNode2)
		{
			theNode2->Release();
			theNode2 = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return val;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (theNode2) theNode2->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Fx Rate");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


long ARMLOCAL_TestGetFxCurve(const CCString& ccy1,
							 const CCString& ccy2,
							 ARM_Date asOf,
							 const CCString& cvName,
							 bool withUSD,
							 double &val)
{
	CCString xmlResponse;
	CCString msgList;

	char sDate[11];
    sprintf(sDate, "%04d%02d%02d", asOf.GetYear(), asOf.GetMonth(), asOf.GetDay());

	int isInverse = 0;

	CCString commande;
	CCString xmlReq;

	try
	{
		// ** GIGASPACE
		std::string xmlOutput ;
		ARM_XGigaToolKit::doGetFX((const char*)ccy1,(const char*)ccy2,(const char*)cvName,asOf,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			isInverse=1;
			ARM_XGigaToolKit::doGetFX((const char*)ccy2,(const char*)ccy1,(const char*)cvName,asOf,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				// SUMMIT
				try
				{
					isInverse = 0;
					
					xmlReq = (CCString)"<Request><Id>";
					xmlReq = xmlReq + cvName + (CCString)"</Id><Type>FXRATE</Type><Vol>FXRATE/";
					xmlReq = xmlReq + ccy1 + (CCString)"/" + ccy2 + (CCString)"</Vol><AsOfDate>";
					xmlReq = xmlReq + sDate + (CCString)"</AsOfDate></Request>";
					
					commande = "c_curves:GetCommset";
					
					etoolkit_execute(commande,xmlReq,xmlResponse,msgList);
					
					if (xmlResponse == "<Response><COMMSET_LIST/></Response>")
					{
						isInverse = 1;
						
						xmlReq.Replace((CCString)"FXRATE/" +  ccy1 + (CCString)"/" + ccy2, (CCString)"FXRATE/" +  ccy2 + (CCString)"/" + ccy1);
						
						etoolkit_execute(commande,xmlReq,xmlResponse,msgList);
					}
				}
				catch (...)
				{
					try
					{
						if (isInverse == 0)
						{
							isInverse = 1;
							
							xmlReq.Replace((CCString)"FXRATE/" +  ccy1 + (CCString)"/" + ccy2, (CCString)"FXRATE/" +  ccy2 + (CCString)"/" + ccy1);

							etoolkit_execute(commande,xmlReq,xmlResponse,msgList);
						}
					}
					catch(...)
					{
						if (withUSD)
						{
							CCString msg("Pb in getting FX Spot");

							throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								(const char *) msg);
						}
						else
							return ARM_KO;
					}
				}
			}
		}

		val = ARMLOCAL_ParseFxCurve(xmlResponse);

		if (isInverse == 1)
			val = 1.0 / val;

		return ARM_OK;
	}
	catch (...)
	{
		return ARM_KO;
	}
}

double ARMLOCAL_GetFxCurve(const CCString& ccy1,
						   const CCString& ccy2,
						   ARM_Date asOf,
						   double amount,
						   const CCString& cvName)
{
	long retCode;
	double val;

	try
	{
		if ( !(ccy1 == "USD") && !(ccy2 == "USD") )
		{
			retCode = ARMLOCAL_TestGetFxCurve(ccy1,
											  ccy2,
											  asOf,
											  cvName,
											  false,
											  val);
		}
		else
		{
			retCode = ARMLOCAL_TestGetFxCurve(ccy1,
											  ccy2,
											  asOf,
											  cvName,
											  true,
											  val);
		}

		if ( (retCode == ARM_KO) && !(ccy1 == "USD") && !(ccy2 == "USD") )
		{
			double val1;
			retCode = ARMLOCAL_TestGetFxCurve(ccy1,
											  "USD",
											  asOf,
											  cvName,
											  true,
											  val1);
			if (retCode == ARM_KO)
			{
				CCString msg("Pb in getting FX Spot");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}


			double val2;
			retCode = ARMLOCAL_TestGetFxCurve("USD",
											  ccy2,
											  asOf,
											  cvName,
											  true,
											  val2);

			if (retCode == ARM_KO)
			{
				CCString msg("Pb in getting FX Spot");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			val = val1*val2;
		}
		else if (retCode == ARM_KO)
		{
			CCString msg("Pb in getting FX Spot");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		return val * amount;
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		CCString msg("Pb in getting FX Spot");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



long ARMLOCAL_GetDatesAndDatasForFxVol(const char* chaineXML,
									   ARM_Date asOfDate,
									   const char* crossCal,
									   ARM_Matrix* &Volatility,
									   ARM_Vector* &vYearTerms,
									   VECTOR<CCString>* &vMaturity)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting FX Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		if (XMLDoc->selectNodes((_bstr_t)("Response/Entity/COMMSET/CommData/COMMDATA"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FX Volatility");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			vYearTerms = new ARM_Vector(nbNodes);
			Volatility = new ARM_Matrix(nbNodes,1);

			long Nb;
			char cMatu;
			long freqId;

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					vMaturity->push_back(ff1);

					sscanf(ff1, "%d%c", &Nb, &cMatu);

					cMatu = toupper(cMatu);

					if ( cMatu == 'D' ) // Ex : "1D"
						freqId = K_DAILY;
					else if ( cMatu == 'W' )
						freqId = K_WEEKLY;
					else if ( cMatu == 'M' ) 
						freqId = K_MONTHLY;
					else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
						freqId = K_ANNUAL;
					else
					{
						hr = S_FALSE;

						CCString msg((CCString)"Invalid XML string for getting FxVol");

						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
										 (const char *) msg);
					}

					ARM_Date tmpDate (asOfDate);
					tmpDate.AddPeriodMult(freqId,Nb,(char*)crossCal);

					vYearTerms->Elt(indexNode) = (tmpDate - asOfDate) / 365.;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Volatility->Elt(indexNode,0) = atof(ff1) * 100.;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				listItem->Release();
				listItem=NULL;
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return ARM_OK;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Fx Volatility");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	return ARM_KO;
}




ARM_Object* ARMLOCAL_ParseObject(const char* chaineXML, const char* typeDeal, const ARM_Date& date, const char* filter, CCString& bookName, CCString& structureId, CCString& custId, CCString& dealId)
{
	try
	{
		double prime;
		ARM_Date primeDate;
		CCString primeCcy = "";
		ARM_ReferenceValue rPremium;

		structureId = "";
		custId = "";
		dealId = "";
		bookName = "";

		if (strcmp(typeDeal,"SWOPT") == 0)
			return ARMLOCAL_ParseSwaption(chaineXML, bookName, structureId, custId, dealId, rPremium);
		else if (strcmp(typeDeal,"SWAP") == 0)
		{
			VECTOR<string> listAssetId;
			ARM_Security* sec = ARMLOCAL_ParseSwap(chaineXML, date, bookName, custId, dealId,listAssetId);

			return sec;
		}
		else if (strcmp(typeDeal,"PRCS") == 0)
		{
			ARM_Security* sec = ARMLOCAL_ParseGenPRCS(chaineXML, date, bookName, structureId, custId, dealId);

			return sec;
		}
		else if (strcmp(typeDeal,"FXOPT") == 0)
		{
			ARM_Security* sec = ARMLOCAL_ParseFxOption(chaineXML, bookName, structureId, primeDate, prime, primeCcy, custId, dealId);

			return sec;
		}
		else if (strcmp(typeDeal,"ACCRUALOPTION") == 0)
		{
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseFlexAcc(chaineXML,date, bookName, custId, dealId,1);

			return sec;
		}
		else if (strcmp(typeDeal,"CRF") == 0)
		{
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseCRF(chaineXML,date, bookName, custId, dealId,1);

			return sec;
		}
		else if (strcmp(typeDeal,"CRA") == 0)
		{
			ARM_Vector* CraPricing = new ARM_Vector();
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseCRA(chaineXML, date, bookName, custId, dealId, CraPricing, 1);
			delete CraPricing;
			CraPricing=NULL;
			return sec;
		}
		else if (strcmp(typeDeal,"SPDOPT") == 0)
		{
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseSpreadOption(chaineXML,date, bookName, custId, dealId, rPremium,1);

			return sec;
		}
		else if (strcmp(typeDeal,"IRG") == 0)
		{
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseCap(chaineXML,date, bookName, custId, dealId, rPremium,1);

			return sec;
		}
		else if (strcmp(typeDeal,"EXOTIC") == 0)
		{
			VECTOR<string> listAssetId;
			ARM_Object* sec = ARMLOCAL_ParseExotic(chaineXML, date, filter, structureId, bookName, custId, dealId,listAssetId);

			return sec;
		}
		else if(strcmp(typeDeal,"EXOTIC.RA") == 0)
		{
			VECTOR<string> listAssetId;
			ARM_Object* sec = ARMLOCAL_ParseCorridorSpreadOption(chaineXML, date, bookName, filter, custId, dealId,listAssetId, 1, 0);
			return sec;
		}
		else if(strcmp(typeDeal,"SWAP.RA") == 0)
		{
			VECTOR<string> listAssetId;
			ARM_Object* sec = ARMLOCAL_ParseCorridorSpreadOption(chaineXML, date, bookName, filter, custId, dealId, listAssetId, 1, 1);
			return sec;
		}

		else if (strcmp(typeDeal,"MATURITYCAP") == 0)
		{
			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseMaturityCap(chaineXML,date, bookName, custId, dealId, rPremium,1);
			
			return sec;
		}
		else if (strcmp(typeDeal,"RFTARN") == 0)
 		{
 			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseTarn(chaineXML,date, bookName, custId, dealId, rPremium,1);
 			
 			return sec;
 		}
		// Bermuda swaption Calculator
		else if (strcmp(typeDeal,"BERM") == 0)
 		{
 			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseBermudaSwaption(chaineXML, date, bookName, custId, dealId, rPremium,1);
 			
			return sec;
 		}
		else if (strcmp(typeDeal,"MEMORYSO") == 0)
 		{
 			ARM_Object* sec = (ARM_Object*) ARMLOCAL_ParseMemorySO(chaineXML, date, bookName, custId, dealId, 1);
 			
			return sec;
 		}
		// Caption Calculator
		else if (strcmp(typeDeal,"ALMCAPTION") == 0)
 		{
 			ARM_Security* sec = (ARM_Security*) ARMLOCAL_ParseCaption(chaineXML, date, bookName, custId, dealId, rPremium,1);
 			
			return sec;
 		}
		else if ( strcmp(typeDeal, "FXOPTSTRIP") == 0 )
		{
			ARM_Object*	sec = (ARM_Object*) ARMLOCAL_ParseFxOptionStrip(chaineXML, date, bookName, structureId, custId, dealId);
			return	sec;
		}
		// FXVanilla Calculator
		else if ( strcmp(typeDeal, "FXSTRIP") == 0 )
		{
			ARM_Object*	sec = (ARM_Object*) ARMLOCAL_ParseFxStrip(chaineXML, date, bookName, structureId, custId, dealId);
			return	sec;
		}
		else if (strcmp(typeDeal, "CFXSTRIP") == 0 )
		{
			ARM_Object* sec = (ARM_Object*) ARMLOCAL_ParseFxStripCalculator(chaineXML, date, bookName, structureId, custId, dealId);
			return sec;
		}
		// RA double Calculator
		else if ( strcmp(typeDeal, "RNG_DOUBLE") == 0 )
		{
			VECTOR<string> listAssetId;
			ARM_Object*	sec = (ARM_Object*) ARMLOCAL_ParseCorridorDblCondition(chaineXML, date, bookName,  custId, dealId,1);
			return	sec;
		}
		
		else if (strcmp(typeDeal,"CDS") == 0) 
			return ARMLOCAL_XML_CDS(chaineXML,bookName, custId, dealId,0);
		else if (strcmp(typeDeal,"NTD") == 0) 
			return ARMLOCAL_XML_NTD(chaineXML,bookName, custId, dealId,0);
//		else if (strcmp(typeDeal,"CDO") == 0) 
//			return ARMLOCAL_XML_CDO(chaineXML ,bookName, custId, dealId,0);
// 		else if (strcmp(typeDeal,"CDO2") == 0) 
// 			return ARMLOCAL_XML_CDO2(chaineXML ,bookName, custId, dealId,0);
		else if (strcmp(typeDeal,"FTD") == 0)
			return ARMLOCAL_XML_NTD(chaineXML,bookName, custId, dealId,0);
		else if (strcmp(typeDeal,"CDSOPT") == 0)
			return ARMLOCAL_XML_CDS_OPTION(chaineXML,bookName, custId, dealId,NULL);
		else
			return NULL;
	}
	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Error Unknown in XML parsing for getting Object");

	}
}

ARM_Object* ARMLOCAL_ParseExotic(const char* chaineXML, const ARM_Date& date, const char* filter, 
								 CCString& structureId, CCString& bookName, CCString& custId, CCString& dealId, 
								 VECTOR<string>& listAssetId, int aResetFreqForCorridorOptim)
{
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr); 

		XMLDoc->loadXML(tmpChaine, &bOK);
        
        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in ARMLOCAL_ParseExotic");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	CCString tmpFilter(filter);
	tmpFilter.toUpper();
	ParserManager::Init();
	string sBookName, sStructureId,sCustId,sDealId;

	ARM_Object* sec = NULL;

	try
	{
		MSXML2::IXMLDOMNodeList * resultList = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

		if (XMLDoc->selectNodes(_bstr_t((const char *)((CCString)"Response/EXOTIC/Assets")), &resultList) == S_OK)
		{
			long nbNodes = 0;
			resultList->get_length(&nbNodes);

			for (int i=0; i<nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);
				listItem->selectSingleNode(_bstr_t((const char *)"ASSET/ProductName"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strncmp((const char*)ff1,"FX_OPT_STRIP",12) == 0)
					{
						ARM_Object* prcs = ARMLOCAL_ParseUnderlyingPRCS(chaineXML, date, bookName, structureId, custId, dealId);
						sec = new ARM_StdPortfolio();
						((ARM_StdPortfolio*)sec)->AddInstrument(dynamic_cast<ARM_Security*>(prcs), 1.0, 1.0, 0.0);
						break;
					}
				}
			}
			if (sec == NULL)
			{
				sec = ParserManager::BuildInstrument("EXOTIC", XMLDoc, date, (const char*)tmpFilter,
																 sBookName, sStructureId,
																 sCustId, sDealId, listAssetId,
																 aResetFreqForCorridorOptim);

				bookName = sBookName.c_str();
				structureId = sStructureId.c_str();
				custId = sCustId.c_str();
				dealId = sDealId.c_str();
			}
		}
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating Exotic Instrument from Summit");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	if (XMLDoc)
		XMLDoc->Release();

	return sec;
}

ARM_Object* ARMLOCAL_ParseGlobalExotic(	const char* chaineXML, const ARM_Date& date, const char* filter, 
										CCString& structureId, CCString& bookName, CCString& custId, CCString& dealId, 
										VECTOR<string>& listAssetId, int aResetFreqForCorridorOptim)
{
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr); 

		XMLDoc->loadXML(tmpChaine, &bOK);
        
        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in ARMLOCAL_ParseExotic");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	CCString tmpFilter(filter);
	tmpFilter.toUpper();
	ParserManager::Init();
	string sBookName, sStructureId,sCustId,sDealId;

	ARM_Object* portfolio = NULL;

	try
	{
		MSXML2::IXMLDOMNodeList * resultList = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

		if (XMLDoc->selectNodes(_bstr_t((const char *)((CCString)"Response/EXOTIC/Assets")), &resultList) == S_OK)
		{
			long nbNodes = 0;
			resultList->get_length(&nbNodes);

			for (int i=0; i<nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);
				listItem->selectSingleNode(_bstr_t((const char *)"ASSET/ProductName"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strncmp((const char*)ff1,"FX_OPT_STRIP",12) == 0)
					{
						vector<ARM_Object*> assets;
						ARM_GP_Vector weights;
						ARM_GP_Vector prices;

						ARM_Object* prcs = ARMLOCAL_ParseUnderlyingPRCS(chaineXML, date, bookName, structureId, custId, dealId);

						if (prcs)
						{
							assets.push_back(dynamic_cast<ARM_Object*>(prcs));
							weights.push_back(1.0);
							prices.push_back(1.0);

							portfolio = new ARM_GlobalPortfolio(assets, weights, prices);
							
							break;
						}				
					}
				}
			}
			if (portfolio == NULL)
			{
				portfolio = ParserManager::BuildInstrument("GLOBLAEXOTIC", XMLDoc, date, (const char*)tmpFilter,
															sBookName, sStructureId,
															sCustId, sDealId, listAssetId,
															aResetFreqForCorridorOptim);

				bookName = sBookName.c_str();
				structureId = sStructureId.c_str();
				custId = sCustId.c_str();
				dealId = sDealId.c_str();
			}
		}
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating Global Exotic Instrument from Summit");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	if (XMLDoc)
		XMLDoc->Release();

	return portfolio;
}


ARM_Object* ARMLOCAL_ParseMemorySO(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, long isETK)
{
		HRESULT hr;
		VARIANT_BOOL bOK;
		MSXML2::IXMLDOMDocument *XMLDoc = NULL;
		CCString structureId("");

		try
		{
			hr = CoInitialize(NULL); 

			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

            wchar_t* wcharStr = constchar2wchar(chaineXML);

			_bstr_t tmpChaine(wcharStr); 

			XMLDoc->loadXML(tmpChaine, &bOK);
            
            if (wcharStr)
               delete wcharStr;
		}

		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in ARMLOCAL_ParseMemorySO");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		ParserManager::Init();
		string sBookName, sStructureId,sCustId,sDealId;
		VECTOR<string> listAssetId;

		ARM_Object* sec = NULL;

		try
		{
			sec = ParserManager::BuildInstrument("MEMORYSO", XMLDoc, date, "",
															 sBookName, sStructureId,
															 sCustId, sDealId, listAssetId);
		}
		catch (Exception& e)
		{
			throw e;
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating Exotic Instrument from Summit");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		bookName = sBookName.c_str();
		structureId = sStructureId.c_str();
		custId = sCustId.c_str();
		dealId = sDealId.c_str();

		if (XMLDoc)
			XMLDoc->Release();

		return sec;
}



ARM_Object* ARMLOCAL_ParseMemoryCAP(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, long isETK)
{
		HRESULT hr;
		VARIANT_BOOL bOK;
		MSXML2::IXMLDOMDocument *XMLDoc = NULL;
		CCString structureId("");

		try
		{
			hr = CoInitialize(NULL); 

			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

            wchar_t* wcharStr = constchar2wchar(chaineXML);

			_bstr_t tmpChaine(wcharStr); 

			XMLDoc->loadXML(tmpChaine, &bOK);
            
            if (wcharStr)
               delete wcharStr;
		}

		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in ARMLOCAL_ParseMemoryCap");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		ParserManager::Init();
		string sBookName, sStructureId,sCustId,sDealId;
		VECTOR<string> listAssetId;

		ARM_Object* sec = NULL;

		try
		{
			sec = ParserManager::BuildInstrument("MEMORYCAP", XMLDoc, date, "",
															 sBookName, sStructureId,
															 sCustId, sDealId, listAssetId);
		}
		catch (Exception& e)
		{
			throw e;
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating Exotic Instrument from Summit");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		bookName = sBookName.c_str();
		structureId = sStructureId.c_str();
		custId = sCustId.c_str();
		dealId = sDealId.c_str();

		if (XMLDoc)
			XMLDoc->Release();

		return sec;
}


ARM_Swaption* ARMLOCAL_ParseSwaption(const char* chaineXML, CCString& bookName, CCString& structureId, CCString& custId, CCString& dealId, ARM_ReferenceValue& rPremium, long isEtk)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, *item = NULL;
	MSXML2::IXMLDOMNode * floatingNode = NULL;
	MSXML2::IXMLDOMNode * fixedNode = NULL;
	long nbNodes;

	ARM_Date startDate;
	ARM_Date endDate;
	int rcvPay;
	double strike;
	char index[20];
	char term[20];
    ARM_Date maturity;
	int liborType;
    double swapYearTerm = -1000000.0;
	int resetFreq = K_DEF_FREQ;
	int payFreq = K_DEF_FREQ;
	string optionType;
	char resetCal[4];
	char payCal[4];

	double dPorS = 1.0;

	ARM_ReferenceValue* notional = NULL;

	ARM_Swaption* newSwaption = NULL;
	ARM_Swaption* Swaption = NULL;
	ARM_SwapLeg* fixLeg = NULL;
	ARM_SwapLeg* varLeg = NULL;
	ARM_Swap* swap = NULL;
	ARM_Currency* currency = NULL;

	VECTOR<double> PrimeDate;
	VECTOR<double> PrimeVal;
	VECTOR<double> PrimeCcy;
	VECTOR<double> dExerDates;
	VECTOR<double> dFees;
	ARM_Vector* exerDates = NULL;
	ARM_ExerciseStyle* exerStyle = NULL;
	ARM_ExerciseStyle* CStyle    = NULL;
	ARM_ReferenceValue* spreads = NULL;
	ARM_ReferenceValue* strikes = NULL;
	string statusExercised;
	string status;

	structureId = "";

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Swaption");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		CCString tmpChaine;

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Back/BACK";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Back/BACK";

		if (XMLDoc->selectSingleNode(_bstr_t((const char *)(tmpChaine)), &item) == S_OK)
		{
			status = GetStringFromXMLNode(item, "TermAssignStatus");
			item->Release();
		}

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Env/ENV";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Env/ENV";

		// Récupération du customer id
		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Swaption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);
			custId = GetCustumerId(listItem);
			dealId = GetDealId(listItem);
			bookName = GetBook(listItem);
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		// Recuperation de la partie optionnelle
		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Option/OPTION";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Option/OPTION";

		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;
				CCString msg((CCString)"Invalid XML string for getting Swaption");
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &listItem);

			strike = XML_doubleNodeTreating(listItem,"Strike") * 100.;
			statusExercised = GetStringFromXMLNode(listItem, "Exercised");
			dPorS = GetPorS(listItem);

			listItem->selectSingleNode(_bstr_t((const char *)"PorC"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp(ff1,"C") == 0)
					rcvPay = K_PAY;
				else
					rcvPay = K_RCV;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			//Check the option style (BERM/EURO)
			optionType = GetStringFromXMLNode(listItem, "Style");
			maturity = GetDateFromXMLNode(listItem,"ExpDate");
			CStyle = GetExerciceStyle(listItem);

			if (listItem->selectNodes(_bstr_t((const char *)"OpEvents/OPEVENT"), &resultList2) == S_OK)
			{
				long nbOpEvents;
				resultList2->get_length(&nbOpEvents);

				for (int i = 0; i < nbOpEvents; i++)
				{
					hr=resultList2->get_item(i, &item);

					item->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (resultat) SysFreeString(resultat);
						theNode->Release();
						theNode=NULL;

						//Premium Event:
						if ( (strcmp((const char*) ff1, "PRM") == 0) || (strcmp((const char*) ff1, "F/M") == 0) )
						{
							PrimeDate.push_back(GetDateFromXMLNode(item,"Date").GetJulian()); //or ADate ? 

							item->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode);
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								double amount = atof((const char*) ff1);
								PrimeVal.push_back(amount);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							item->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								PrimeCcy.push_back((double)ARM_GetCountrySymbol(ff1));

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							if (item) item->Release();
						}
						//Bermuda exercise type
						else if (strcmp((const char*) ff1, "CL") == 0) 
						{
							dExerDates.push_back(GetDateFromXMLNode(item,"ADate").GetJulian()); //or ADate ? 
							//Doit-on ranger en order croissant les dates d'exercise?
							item->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode);
							if (theNode!=NULL)
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								double amount = atof((const char*) ff1);
								dFees.push_back(amount);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}
						}
					}
				}
			
				if (resultList2)
				{
					resultList2->Release();
					resultList2 = NULL;
				}
			}

			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		// Recuperation de la partie swap
		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Assets/ASSET";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Assets/ASSET";

		if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 2)
			{
				hr = S_FALSE;
				CCString msg((CCString)"Invalid XML string for getting Swaption");
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strcmp(ff1,"FLO") == 0)
						floatingNode = listItem;
					else
						fixedNode = listItem;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
			}

			startDate = GetStartDate(floatingNode);
			endDate = GetEndDate(floatingNode);

			floatingNode->selectSingleNode(_bstr_t((const char *)"INTEREST_dmIndex"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(index,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			floatingNode->selectSingleNode(_bstr_t((const char *)"INTEREST_Term"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(term,ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			liborType = FromStringToIndexType(term,index);

			//Bermuda swaption
			if  (strcmp(optionType.c_str(), "BERM") == 0)
			{
				//fixed leg
				double fixRate = GetCapStrike(fixedNode);
				GetPayCalendar(fixedNode, payCal);
				GetResetCalendar(fixedNode, resetCal);

				fixLeg = new ARM_SwapLeg(startDate, 
										 endDate, 
										 fixRate, 
										 rcvPay,
										 GetPayFreq(fixedNode), 
										 GetAssetDayCount(fixedNode), 
										 GetDecompFreq(fixedNode),
										 GetPayTiming(fixedNode),
										 GetPayIntRule(fixedNode),
										 GetStubRule(fixedNode),
										 GetCcy(fixedNode),
										 payCal,
										 GetNotionalExchangeFlag(fixedNode),
										 NULL, //ref date
										 1);   //adjust start

				//float leg
				GetPayCalendar(floatingNode, payCal);
				GetResetCalendar(floatingNode, resetCal);

				varLeg = new ARM_SwapLeg(startDate,
										 endDate,
										 (ARM_INDEX_TYPE) liborType,
										 rcvPay,
										 GetSpread(floatingNode), 
										 GetResetFreq(floatingNode),
										 GetPayFreq(floatingNode),
										 GetResetTiming(floatingNode),
										 GetPayTiming(floatingNode),
										 GetCcy(floatingNode),
										 GetIntRule(floatingNode),
										 GetResetGap(floatingNode),
										 resetCal,
										 payCal,
										 GetDecompFlag(floatingNode),
										 GetNotionalExchangeFlag(floatingNode),
										 GetStubRule(floatingNode),
										 NULL, //ref date
										 1,    //adjust start
										 GetAssetDayCount(floatingNode),
										 GetFwdRule(floatingNode),
										 GetPayGap(floatingNode));   

				ARM_ReferenceValue* stepUpSpreadOnStarts = GetSwaplegStepUpSpread(floatingNode,true);
				if (stepUpSpreadOnStarts)
				{
					ARM_ReferenceValue* stepUpSpread = ARM_FromStartRefToResetRef(stepUpSpreadOnStarts, varLeg);
					delete stepUpSpreadOnStarts;

					varLeg->SetVariableSpread(stepUpSpread);
					delete stepUpSpread;
				}

				swap = new ARM_Swap(fixLeg, varLeg);
			
				//Exercise
				exerDates = new ARM_Vector(dExerDates);
				exerStyle = new ARM_ExerciseStyle(exerDates);
				strikes = GetStepUpFixCoupon(fixedNode);

				newSwaption = new ARM_Swaption(swap, 
											   rcvPay,
											   CStyle,//exerStyle,
											   strikes,
											   0);

				ARM_Vector* exerFees = new ARM_Vector(dFees);
				ARM_ReferenceValue* fees = new ARM_ReferenceValue(exerDates, exerFees);
				newSwaption->SetFee(fees);

				delete exerDates;
				delete exerFees;
				delete exerStyle;
				delete CStyle;
				delete strikes;
				delete fixLeg;
				delete varLeg;
				delete swap;
			}
			//European swaption
			else
			{
				newSwaption = new ARM_Swaption(startDate,
											   endDate,
											   rcvPay,
											   K_EUROPEAN,
											   strike,
											   maturity,
											   (ARM_INDEX_TYPE)liborType,
											   0.0,
											   -1000000.0,
											   K_DEF_FREQ,
											   K_DEF_FREQ,
											   GetCcy(floatingNode));
			}
			
			ARM_ReferenceValue* notional = GetNotional(floatingNode);
			newSwaption->SetAmount(notional);

			newSwaption->SetPorS(dPorS);
		
			if ( (status == "TERM") || (statusExercised == "EXER") )
				newSwaption->SetFullterminated(true);

			// Remontée des primes
			if (PrimeDate.size() > 0)
			{
				ARM_Vector* vPrimeDate = new ARM_Vector(PrimeDate.size());
				ARM_Vector* vPrimeVal = new ARM_Vector(PrimeDate.size());
				ARM_Vector* vPrimeCcy = new ARM_Vector(PrimeDate.size());

				for (int i = 0; i < PrimeDate.size(); i++)
				{
					vPrimeDate->Elt(i) = PrimeDate[i];
					vPrimeVal->Elt(i) = PrimeVal[i];
					vPrimeCcy->Elt(i) = PrimeCcy[i];
				}

				ARM_ReferenceValue rPrime((ARM_Vector*)vPrimeDate->Clone(),(ARM_Vector*)vPrimeVal->Clone(),(ARM_Vector*)vPrimeCcy->Clone());

				if (vPrimeDate)
					delete vPrimeDate;
				vPrimeDate = NULL;

				if (vPrimeVal)
					delete vPrimeVal;
				vPrimeVal = NULL;

				if (vPrimeCcy)
					delete vPrimeCcy;
				vPrimeCcy = NULL;

				rPremium = rPrime;
			}

			if (notional)
				delete notional;
			notional = NULL;
		}

		floatingNode->Release();
		floatingNode = NULL;
		fixedNode->Release();
		fixedNode = NULL;
       
		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		return newSwaption;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"ParseSwaption : Error in XML parsing for getting ExoSwaption");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}




ARM_Option* ARMLOCAL_ParseFxOption(const char* chaineXML, CCString& bookName, CCString& structureId, ARM_Date& primeDate, double& prime, CCString& primeCcy, CCString& custId, CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	char ccy1[4];
	char ccy2[4];
	int isCallOrPut;
	double mainAmount;
	double dPorS = 1.0;

    ARM_Date maturity;
	double strike;
	ARM_Option* newFxOption = NULL;

	structureId = "";

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting FxOption");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		// Récupération du customer id
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/FXOPT_TR/Env/ENV")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FxOption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			custId = GetCustumerId(listItem);

			dealId = GetDealId(listItem);
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		// Recuperation de la partie optionnelle
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/FXOPT_TR/Option/OPTION")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FxOption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			maturity = GetOptionExpiry(listItem);

			primeDate = GetDateFromXMLNode(listItem,"PREMDATA_Date");

			strike = XML_doubleNodeTreating(listItem,"Strike");

			prime = XML_doubleNodeTreating(listItem,"PREMDATA_Premium");

			listItem->selectSingleNode(_bstr_t((const char *)"PREMDATA_Ccy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				primeCcy = ff1;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			isCallOrPut = GetPorC(listItem);

			dPorS = GetPorS(listItem);

			listItem->Release();
			listItem=NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/FXOPT_TR/FxOptList/FXOPTION")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				//if (XMLDoc) XMLDoc->Release();
				//if (resultList) resultList->Release(); //Non release fait dans le catch
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FxOption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t((const char *)"MainCcy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				strcpy(ccy1,(char *)ff);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			listItem->selectSingleNode(_bstr_t((const char *)"MoneyCcy"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				strcpy(ccy2,(char *)ff);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			mainAmount = XML_doubleNodeTreating(listItem,"MainAmount");
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		// Recuperation du book
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/FXOPT_TR/Env/ENV")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FxOption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			bookName = GetBook(listItem);

			listItem->Release();
			listItem=NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		ARM_Currency Ccy1(ccy1);
		ARM_Currency Ccy2(ccy2);
		ARM_Forex newForex(&Ccy1,&Ccy2);

		newFxOption = new ARM_Option(&newForex,maturity,strike,isCallOrPut,K_EUROPEAN,K_YIELD,(ARM_Date)ARM_DEFAULT_DATE);
		ARM_ReferenceValue amount(mainAmount);
		newFxOption->SetAmount(&amount);
		newFxOption->SetPorS(dPorS);

		return newFxOption;
	}
	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting FX option");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



ARM_FlexAccretSwaption* ARMLOCAL_ParseFlexAcc(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, long isEtk)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *xmlDoc = NULL;
	
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&xmlDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		xmlDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr; 
	}
	catch(...)
	{
		if (xmlDoc) xmlDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Accrual Option");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
	
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode* item = NULL, * theNode = NULL, * fixItem = NULL, * floatItem = NULL;
	long nbAssets = 0;
	BSTR resultat = NULL;

	double floatNotional, floatSpread, fixRate;
	int payTime;
	long frequence = K_ANNUAL, dayCount;
	int stubRule = K_SHORTSTART;
	int nbCP4A = 0;
	ARM_INDEX_TYPE armIndex = LIBOR3M;
	ARM_Date startDate, endDate;
	ARM_Currency * currency = NULL;
	ARM_ExerciseStyle * exerciseDates = NULL;
	ARM_Vector * exDates = NULL;
	ARM_FlexAccretSwaption * security = NULL;
	ARM_Vector * echeancier = NULL, * paymentDates = NULL;
	ARM_SwapLeg * fixedLeg = NULL;
	long nbNodes;
	
	try
	{
		CCString tmpChaine;

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Env/ENV";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Env/ENV";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Swaption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &item);

			custId = GetCustumerId(item);

			dealId = GetDealId(item);

			bookName = GetBook(item);
		}

		if (resultList) resultList->Release();
		if (item) item->Release();

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Assets/ASSET";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Assets/ASSET";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList)== S_OK)
		{
			resultList->get_length(&nbAssets);

			if (nbAssets != 2)
			{
				CCString msg((CCString)"Option XML string is not valid");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			// Récupération des 2 chaines XML correspondant à chaque asset, puis détermination de quel est le fix, quel est le float
			for (long indexAsset=0 ; indexAsset<nbAssets ; indexAsset++)
			{
				hr=resultList->get_item(indexAsset, &item);
				if (hr==S_OK && item!=NULL)
				{
					item->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
					if (theNode!=NULL)
					{
						theNode->get_text(&resultat);
						theNode->Release();
						theNode = NULL;
						if (strcmp((const char*)(_bstr_t)resultat,"FLO") == 0)	// Jambe variable
						{
							item->selectSingleNode(_bstr_t((const char *)"Notional"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								floatNotional = atof((const char*) ff1);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							item->selectSingleNode(_bstr_t((const char *)"INTEREST_Spread"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								floatSpread = atof((const char*) ff1) * 100.;

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							char index[5];
							item->selectSingleNode(_bstr_t((const char *)"INTEREST_dmIndex"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								strcpy(index,(char*)(_bstr_t)resultat);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							char term[4];

							item->selectSingleNode(_bstr_t((const char *)"INTEREST_Term"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								strcpy(term,(char*)(_bstr_t)resultat);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}
							try
							{
								armIndex = (ARM_INDEX_TYPE)FromStringToIndexType(term,index);
							}
							catch (...)
							{
								armIndex = LIBOR3M;
							}

							startDate = GetStartDate(item);

							endDate = GetEndDate(item);

							item->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);	// Rem : les 2 jambes ont la même devise
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								currency = new ARM_Currency((char*)(_bstr_t)resultat);

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}
						}
						else	// Jambe fixe
						{
							item->selectSingleNode(_bstr_t((const char *)"INTEREST_Rate"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								fixRate = atof((const char*) ff1) * 100.;

								if (resultat) SysFreeString(resultat);
								theNode->Release();
								theNode=NULL;
							}

							frequence = GetPayFreq(item);

							stubRule = GetStubRule(item);

							dayCount = GetAssetDayCount(item);

							payTime = GetPayTiming(item);							
						}

						if (resultat) SysFreeString(resultat);
					}
					item->Release();
				}
			}
		}

		if (resultList) resultList->Release();

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Option/OPTION/OpEvents/OPEVENT";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Option/OPTION/OpEvents/OPEVENT";
		
		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			long nbOpEvents;
			resultList->get_length(&nbOpEvents);
			long tmpNbOpeEvents = nbOpEvents;

			if (nbOpEvents == 0)
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
                                 "Invalid XML string (no OPEVENT)");
			}

			double amount = 0., theDate = 0., theADate, lastExpDate = MAXDATE, lastExpDate4Opt = MAXDATE, lastExercisedDate = MAXDATE;
			// 1) Récupérer la dernière date d'exercice passée ( <= asofdate )
			for (int i=0, j=0; i<nbOpEvents; i++)
			{
				hr=resultList->get_item(i, &item);

				item->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
				if (theNode!=NULL)
				{
					theNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode=NULL;

					if (strcmp((const char*) ff1, "CL") == 0)
					{
						theDate = GetDateFromXMLNode(item,"Date").GetJulian();

						theADate = GetDateFromXMLNode(item,"ADate").GetJulian();

						item->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode);
						if (theNode!=NULL)
						{
							theNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							amount = atof((const char*) ff1);

							if (resultat) SysFreeString(resultat);
							theNode->Release();
							theNode=NULL;
						}

						if (item) item->Release();

						// Dates passées uniquement
						if (theADate <= date.GetJulian())
						{
							// Si exercice (si amount = 1)
							if (amount > 0.9)
							{
								j = i + 1;
								lastExercisedDate = theDate;
							}
						}

						lastExpDate = theADate;
					}
					else
						tmpNbOpeEvents--;
				}
			}

			lastExpDate4Opt = lastExpDate;

			// 2) On commence la liste des exercices futurs après la dernière trouvée (exclue)
			// Le tableau ARM doit contenir exactement le bon nombre de dates.
			exDates = new ARM_Vector(tmpNbOpeEvents - j);

			for (i=0, j=0 ; i<nbOpEvents ; i++)
			{
				hr=resultList->get_item(i, &item);

				item->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
				if (theNode!=NULL)
				{
					theNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode=NULL;

					if (strcmp((const char*) ff1, "CL") == 0)
					{
						theDate = GetDateFromXMLNode(item,"Date").GetJulian();

						theADate = GetDateFromXMLNode(item,"ADate").GetJulian();

						item->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode);
						if (theNode!=NULL)
						{
							theNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							amount = atof((const char*) ff1);

							if (resultat) SysFreeString(resultat);
							theNode->Release();
							theNode=NULL;
						}

						// Cas où aucun exercice à la 1ère expiry date ou bien la date d'exercice event(i)->Date est une date future (strictement),
						// ie: > lastExercisedDate.
						if ((lastExercisedDate == MAXDATE)
							||
							(theDate > lastExercisedDate))
						{
							if (theADate == MAXDATE)	// full term : prix nul donc on sort
							{
								throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Full term");
							}

							exDates->Elt(j++) = theADate;
						}
					}
				}
				if (item) item->Release();
			}

			// 3) Création de l'objet arm (expiry date, exercice (0/1) )
			exerciseDates = new ARM_ExerciseStyle(exDates);

			// ---------------------------------------------------------------
			// Gestion des exercices précédents
			// ---------------------------------------------------------------

			// ---------------------------------------------------------------
			// 1) Calcul de nbCP4A : Nb Current Period For Accrued, c'est le
			// nombre de périodes entière passée sans exercice, par ex. si on a
			// exercé à la dernière expiry date, nbCP4A := 0.

			lastExpDate = startDate.GetJulian();

			for (i = 0; i < nbOpEvents ; i++)
			{
				hr=resultList->get_item(i, &item);

				item->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
				if (theNode!=NULL)
				{
					theNode->get_text(&resultat);
					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (resultat) SysFreeString(resultat);
					theNode->Release();
					theNode=NULL;

					if (strcmp((const char*) ff1, "CL") == 0)
					{
						theADate = GetDateFromXMLNode(item,"ADate").GetJulian();

						item->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode);
						if (theNode!=NULL)
						{
							theNode->get_text(&resultat);
							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							amount = atof((const char*) ff1);

							if (resultat) SysFreeString(resultat);
							theNode->Release();
							theNode=NULL;
						}

						double expDate = theADate;
						bool isEx = false;
						if (amount > 0.01)
							isEx = true;

						// expiry date dans le passé (<0) ou le présent (==0)
						if (expDate <= date.GetJulian())
						{
							if (isEx)
							{
								lastExpDate = expDate;
								nbCP4A = 0;
							}
							else
							{
								nbCP4A++;
							}
						}
					}
				}
				if (item) item->Release();
			}

			// -------------------------------------------------------------------
			// 2) On récupère la date de paiement (taux fixe) immédiatement
			// supérieure à la dernière date d'exercice, ça nous done la start
			// date pour ARM (itsLastExpDate). Si aucun exercice passé, on a la
			// start date du trade.

			// Création de l'échéancier ([startDate, paymentDates, endDate])

			fixedLeg = new ARM_SwapLeg(startDate, endDate, fixRate, K_PAY, frequence, dayCount, K_COMP_PROP, payTime, K_UNADJUSTED, stubRule, ARM_DEFAULT_CURRENCY, NULL, K_NX_NONE, NULL);

			// on prend les end dates car elles ne sont pas ajustees
			paymentDates = fixedLeg->GetFlowEndDates();

			// On met la startDate au début de l'échéancier, les dates de paiement ensuite
			if (paymentDates->Elt(0) != startDate.GetJulian())
			{
				echeancier = new ARM_Vector(paymentDates, paymentDates->GetSize()+1, 0, paymentDates->GetSize()-1, 1);
				echeancier->Elt(0) = startDate.GetJulian();
			}

			// On rajoute la endDate à la fin de l'échéancier
			if (paymentDates->Elt(paymentDates->GetSize()-1) < endDate.GetJulian())
			{
				echeancier->insert(endDate.GetJulian());
			}
				
			// Parcours de l'échéancier
			for (i = 0; i < tmpNbOpeEvents ; i++)
			{
				if (ARM_Date(echeancier->Elt(i)).GetJulian() >= lastExpDate)
				{
					lastExpDate = ARM_Date(echeancier->Elt(i)).GetJulian();
					break;
				}
			}
			startDate = ARM_Date(lastExpDate);
		}
		if (fixedLeg) delete fixedLeg;
		if (echeancier) delete echeancier;
		if (resultList) resultList->Release(); resultList = NULL;
		
		security = new ARM_FlexAccretSwaption(startDate, endDate, fixRate, nbCP4A, K_PAY, frequence, armIndex, floatSpread, exerciseDates, currency);

		if (currency) delete currency;
		if (exerciseDates) delete exerciseDates;
		if (exDates) delete exDates;

		if (xmlDoc) xmlDoc->Release(); xmlDoc = NULL;

	}
	catch (...)
	{
		if (fixedLeg) delete fixedLeg;
		if (echeancier) delete echeancier;
		if (resultat) SysFreeString(resultat);
		if (theNode) theNode->Release();
		if (resultList) resultList->Release();
		if (currency) delete currency;
		if (exerciseDates) delete exerciseDates;
		if (exDates) delete exDates;
		if (security) delete security;
		if (xmlDoc) xmlDoc->Release(); xmlDoc = NULL;

		try
		{
			throw;
		}
		catch (Exception&)
		{
			throw;
		}
		catch (...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Global error creating security");
		}
	}

	return security;
}




ARM_CRFCalculator* ARMLOCAL_ParseCRF(const char* chaineXML,
									 const ARM_Date& date,
									 CCString& bookName,
									 CCString& custId,
									 CCString& dealId,
                                     long isEtk,
                                     double FxSpot)
{
	ARM_CRFCalculator* security = NULL;

	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument* xmlDoc = NULL;
	
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&xmlDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		xmlDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr; 
	}

	catch(...)
	{
		if (xmlDoc) xmlDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Accrual Option");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}
	
	MSXML2::IXMLDOMNodeList* resultList = NULL, *resultList2 = NULL;
	MSXML2::IXMLDOMNode* item = NULL, * theNode = NULL, * theNode2 = NULL, * listItem2 = NULL;
	BSTR resultat = NULL;
	long nbNodes;
	long nbAssets;

	ARM_Date startDate, endDate, fixEndDate;
	VECTOR<double> spreadDate;
	VECTOR<double> spreadVal;
	int fundFreq;
	int fundDaycount;
	VECTOR<double> notionalDate;
	VECTOR<double> notionalVal;

	ARM_Currency* fixCurrency   = NULL;
    ARM_Currency* floatCurrency = NULL;

	double dPorS = 1.0;
	int fixedDaycount;
	int cpnFreq;
	char cpnPayCal[4];
	char cpnResetCal[4];
	int cpnDaycount;
	int cpnIndexDaycount;
	int cpnTiming;
	char IndexTerm[15];
	int cpnResetGap;
	int deal_PorS;
	int stubRule;

	ARM_ReferenceValue* rLeverage    = NULL;
	ARM_ReferenceValue* rCpnMin      = NULL;
	ARM_ReferenceValue* rCpnMax      = NULL;
	ARM_ReferenceValue* rSpread      = NULL;
	ARM_ReferenceValue* rNotional    = NULL;
	ARM_ReferenceValue* strike       = NULL;
	ARM_ReferenceValue* exerFee      = NULL;
	ARM_ReferenceValue* fundNotional = NULL;

	double meanRev;
	double meanRev_default = 10000.0;

	CCString tmpChaine;

	string statusExercised;

	try
	{
		ARM_PricingModelType::ModelType	vModelType = ARM_PricingModelType::HWM1F;

        if (isEtk != 1)
        {
            xmlDoc->selectSingleNode(_bstr_t("/ARMNT_CALL/OPSPEC"), &item);

            string pModel = GetStringFromXMLNode(item, "PModel");

            item->Release();

            if( pModel == "MC" || pModel == "cQGM" )
            {
				vModelType = ARM_PricingModelType::QGM1F;
            }
        }

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Back/BACK";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Back/BACK";

		xmlDoc->selectSingleNode((_bstr_t)((const char*)tmpChaine), &item);
		string status = GetStringFromXMLNode(item, "TermAssignStatus");
		item->Release();

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Env/ENV";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Env/ENV";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Swaption");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr=resultList->get_item(0, &item);

			custId = GetCustumerId(item);
			dealId = GetDealId(item);
			bookName = GetBook(item);

			if (item) item->Release();
		}

		if (resultList) resultList->Release();


		if ( isEtk == 1 )
		   tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Assets/ASSET";
		else
		   tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Assets/ASSET";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbAssets);

			if (nbAssets != 2)
			{
				CCString msg((CCString)"Option XML string is not valid");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			// Récupération des 2 chaines XML correspondant à chaque asset, puis détermination de quel est le fix, quel est le float
			for (long indexAsset=0 ; indexAsset<nbAssets ; indexAsset++)
			{
				hr = resultList->get_item(indexAsset, &item);

				if (hr == S_OK && item != NULL)
				{
					item->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
					if ( theNode!= NULL )
					{
						theNode->get_text(&resultat);
						theNode->Release();
						theNode = NULL;

						if ( strcmp((const char*)(_bstr_t)resultat,"FLO") == 0 )	// Jambe variable
						{
							startDate=GetStartDate(item);
							endDate = GetEndDate(item);
							rSpread = GetSpreadVariable(item);
							if (!rSpread)
								rSpread = new ARM_ReferenceValue(GetSpread(item));
							else
								*rSpread *= 0.01;
							fundFreq = GetPayFreq(item);
							fundDaycount = GetAssetDayCount(item);
							stubRule = GetStubRule(item);

							fundNotional = GetNotional(item);

							floatCurrency = GetCcy(item);
						}
						else	// Jambe fixe
						{
							dPorS = GetPorS(item);
							fixedDaycount = GetAssetDayCount(item);
							cpnFreq = GetPayFreq(item);
							GetPayCalendar(item,cpnPayCal);
							meanRev = GetMeanRev(item);

							rNotional = GetNotional(item);
							strike = GetCRFStrike(item);

                            fixCurrency = GetCcy(item);
						}

						if (resultat) SysFreeString(resultat);
					}

					if (item) item->Release();
				}
			}
		}

		if (resultList) 
			resultList->Release();
		resultList = NULL;

		if ( isEtk == 1 )
		   tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Option/OPTION";
		else
		   tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Option/OPTION";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbAssets);

			if ( nbAssets != 1 )
			{
				CCString msg((CCString)"Option XML string is not valid");

				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			hr = resultList->get_item(0, &item);
			
            cpnIndexDaycount = GetCRFIndexDayCount(item);
			
			cpnDaycount = GetCRFDayCount(item);
			if (cpnDaycount == KNOBASE)
				cpnDaycount = fixedDaycount;

			statusExercised = GetStringFromXMLNode(item, "Exercised");

			cpnTiming = GetCRFTiming(item);

			GetCRFIndexTerm(item, fixCurrency, IndexTerm);

			cpnResetGap = GetCRFResetGap(item);

            fixEndDate = GetCRFFixEndDate(item);

			GetCRFResetCal(item,cpnResetCal);

			// s'il y a une liste Schedule/cCRFLST_I on la parcourt
			item->selectNodes(_bstr_t((const char *)"ProdData/cCRFARM_I/Schedule/cCRFLST_I"), &resultList2);
			
            if ( resultList2 != NULL )
			{
				rLeverage = GetCRFLeverage(item);

				rCpnMin = GetCRFCpnMin(item);
				rCpnMax = GetCRFCpnMax(item);
			}
			else
			// sinon on prend les 3 elements au niveau du Schedule et fixEndDate = startDate
			{
				rLeverage = GetCRFOneLeverage(item);
				rCpnMin = GetCRFOneCpnMin(item);
				rCpnMax = GetCRFOneCpnMax(item);
			}

			exerFee = GetExerFees(item);

			// Purchase or Sale (on inverse le nominal si Sale)
			deal_PorS = GetPorS(item);

			if (item) item->Release();
			if (resultList2) resultList2->Release();

		}

		if (resultList) resultList->Release();
		resultList = NULL;
/*
		ARM_DateStrip FixLegSched(startDate,fixEndDate,cpnFreq,fixedDaycount,
				cpnResetCal,
				K_MOD_FOLLOWING,K_MOD_FOLLOWING,stubRule,
				- (fixCurrency->GetSpotDays()),
				cpnFreq,GETDEFAULTVALUE,
				cpnResetCal,
				K_ADVANCE,K_ARREARS);

		ARM_DateStrip RFLegSched(fixEndDate,endDate,cpnFreq,cpnDaycount,
				cpnResetCal,
				K_MOD_FOLLOWING,K_MOD_FOLLOWING,stubRule,
				cpnResetGap,
				cpnFreq,GETDEFAULTVALUE,
				cpnResetCal,
				cpnTiming,K_ARREARS);

		ARM_DateStrip ExerSched(startDate,endDate,cpnFreq,cpnDaycount,
				cpnResetCal, // if exerCal has not been setted, it's the payCal
				K_MOD_FOLLOWING,K_MOD_FOLLOWING,stubRule,
				- (fixCurrency->GetSpotDays()),
				cpnFreq,GETDEFAULTVALUE,
				cpnResetCal,
				K_ADVANCE,K_ARREARS);
*/
		security = new ARM_CRFCalculator(/*FixLegSched,
										 RFLegSched,
										 ExerSched,*/
										 date,
										 startDate,
										 endDate,
										 *strike,
										 (int) dPorS,
										 fixEndDate,
										 fixedDaycount,
										 cpnDaycount,
										 cpnFreq,
										 cpnTiming,
										 IndexTerm,
										 //"3M",
										 cpnIndexDaycount,
										 cpnResetCal,
										 cpnPayCal,
										 stubRule,
										 cpnResetGap,
										 *rLeverage,
										 *rCpnMin,
										 *rCpnMax,
										 *rSpread,
										 fundFreq,
										 fundDaycount,
										 *rNotional,
										 *exerFee,
										 *fixCurrency,
										 *floatCurrency,
										 *fixCurrency,
                                         *fixCurrency,
                                         *floatCurrency,
										 *fundNotional,
										 vModelType);

		security->SetPorS(deal_PorS);

		if ( (status == "TERM") || (statusExercised == "EXER") )
			security->SetFullterminated(true);

		if ( meanRev != meanRev_default )
		{
			ARM_CurveModelParam MRParam(ARM_ModelParamType::MeanReversion,meanRev,"MRS");
			security->SetMRS(&MRParam);
		}

		// get flags from CONFIG_VECTOR if SUMMIT
		if (isEtk != 1)
		{
			tmpChaine = (CCString)"ARMNT_CALL/CONFIG_VECTOR";
			if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{		
					CCString msg((CCString)"Invalid XML string for getting Config Vector");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
				}

				hr=resultList->get_item(0, &item);

				bool sigmaflag = GetConfigVectorFlag(item, CRF_OSW_CFG_VECT_IDX);
				string sigmaTypeStr = true ?  string("EQUIVALENT"):string("UNKNOWN");
				ARM_SigmaCalibType sigmaType = (ARM_SigmaCalibType) ARM::ARM_ArgConv_SigmaCalibType.GetNumber(sigmaTypeStr);
				security->SetOSWCalibFlag(sigmaType);
				security->SetCapCalibFlag(GetConfigVectorFlag(item, CRF_CAP_CFG_VECT_IDX));
				security->SetFloorCalibFlag(GetConfigVectorFlag(item, CRF_FLOOR_CFG_VECT_IDX));

				if (item) item->Release();
			}

			if (resultList)
				resultList->Release();
			resultList = NULL;

		}

		if (xmlDoc) 
			xmlDoc->Release();
		xmlDoc = NULL;	

		if (rLeverage)
			delete rLeverage;
		rLeverage = NULL;

		if (rNotional)
			delete rNotional;
		rNotional = NULL;

		if (rCpnMin)
			delete rCpnMin;
		rCpnMin = NULL;

		if (rCpnMax)
			delete rCpnMax;
		rCpnMax = NULL;

		if (rSpread)
			delete rSpread;
		rSpread = NULL;

		if (strike)
			delete strike;
		strike = NULL;

		if (fixCurrency)
			delete fixCurrency;
		fixCurrency = NULL;

        if (floatCurrency)
			delete floatCurrency;
		fixCurrency = NULL;

		if (exerFee)
			delete exerFee;
		exerFee = NULL;

		if (fundNotional)
			delete fundNotional;
		fundNotional = NULL;
	}

	catch (Exception e)
	{
		if (resultat) SysFreeString(resultat);
		if (theNode) theNode->Release();
		if (resultList) resultList->Release();		
		if (security) delete security;
		if (xmlDoc) xmlDoc->Release(); xmlDoc = NULL;
		throw e;
	}

	catch (...)
	{
        try
        {
		   if (resultat) SysFreeString(resultat);
		   if (theNode) theNode->Release();
		   if (resultList) resultList->Release();
		   if (security) delete security;
		   if (xmlDoc) xmlDoc->Release(); xmlDoc = NULL;
        }

		catch (...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Global error in creating CRF security");
		}
	}
	
	return(security);
}



ARM_OptionPortfolio* ARMLOCAL_ParseCRA(const char* chaineXML,
									   const ARM_Date& date,
									   CCString& bookName,
									   CCString& custId,
									   CCString& dealId,
									   ARM_Vector* CraPricing,
									   long isEtk)
{
	ARM_OptionPortfolio* security = NULL;
	ARM_SwapLeg* swapleg=NULL;
	ARM_CorridorLeg* corridor=NULL;
	ARM_Vector* vCraPricing= new ARM_Vector(8);
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument* xmlDoc = NULL;


	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), 
                              (void**) &xmlDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		bstr_t tmpChaine(wcharStr); 

		xmlDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (xmlDoc) xmlDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting CRA");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
	
	MSXML2::IXMLDOMNodeList* resultList = NULL, *resultList2 = NULL;
	MSXML2::IXMLDOMNode* item = NULL, * theNode = NULL, * theNode2 = NULL, * listItem2 = NULL;
	BSTR resultat = NULL;
	long nbNodes;
	long nbAssets;

	ARM_Date CstartDate, CendDate;
	VECTOR<double> CspreadDate;
	VECTOR<double> CspreadVal;
	int CpayFreq;
	int CresetFreq;
	int CDaycount;
	ARM_Currency* Ccurrency = NULL;
	double CPorS = 1.0;
	double PorC = 1.0;
	double PorS = 1.0;
	double Cspread = 0.0;
	int CresetPayTiming;
	int CresetRefTiming;
	int Cstub;
	int CdecompPricingFlag = 0; // M.A: Forced to zero because this is the general case
	int CRefIndexType;
	int CPayIndexType;
	char Ccust[20];

	VECTOR<double> CstrikeDate;
	VECTOR<double> CstrikeVal;

	ARM_ReferenceValue* Cstrike = NULL;
	ARM_ReferenceValue* CRefValSpread = NULL;
	ARM_ReferenceValue* CRefBarUp = NULL;
	ARM_ReferenceValue* CRefBarDown = NULL;
	ARM_ReferenceValue* CRefValNotional = NULL;
	ARM_Vector* CvExerciseDates = NULL;
	ARM_ExerciseStyle* CStyle=NULL;

	ARM_Date SstartDate, SendDate;
	double SPorS = 1.0;
	int SDaycount;
	int SpayFreq;
	int SpayTiming;
	int Sstub;
	ARM_Currency* Scurrency = NULL;
	double Sspread=0.0;
	ARM_ReferenceValue* SspreadVariable = NULL;
	int SintRule;
	int SIndexType;
	int SresetFreq;
	int SresetTiming;
	int SresetGap;
	int SpayGap;
	char CpayCal[7];
	char CresetCal[7];
	char SpayCal[7];
	char SresetCal[7];
	int SnotionalEx;
	int SdecompFreq;
	ARM_ReferenceValue* SRefValNotional = NULL;
	char Scust[20];

	string statusExercised;

	int mode_fictif=1;

	CCString tmpChaine;

	try
	{
		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Back/BACK";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Back/BACK";

		xmlDoc->selectSingleNode(_bstr_t((const char *)tmpChaine), &item);
		string status = GetStringFromXMLNode(item, "TermAssignStatus");
		item->Release();

		if ( isEtk == 1 )
		   tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Env/ENV";
		else
		   tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Env/ENV";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting CRA");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &item);
			custId= GetCustumerId(item);
			dealId = GetDealId(item);
			bookName = GetBook(item);

			if (item) item->Release();
		}

		if (resultList) resultList->Release();
		resultList=NULL;

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Option/OPTION";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Option/OPTION";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;
				CCString msg((CCString)"Invalid XML string for getting CRA");
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &item);
			PorS = GetPorS(item);
			ARM_Currency* Cccy = GetOptionCcy(item);
			CPayIndexType=GetBLOBIdx(item,Cccy,2);

			// In the corridorLeg, the payment is always made at the end of the period.
			// Reset timing of the payment index: used for a	Libor pay index.
			CresetPayTiming = GetBLOBTiming(item,2);
			// Reset timing of the reference index.
			CresetRefTiming = GetBLOBTiming(item,1);

			CRefValSpread=GetCRASpread(item);
			CRefBarDown=GetCRABarrierDown(item);
			CRefBarUp=GetCRABarrierUp(item);
			CRefIndexType=GetBLOBIdx(item,Cccy,1);

			if (Cccy)
				delete Cccy;
			Cccy = NULL;

			CvExerciseDates=GetExerciseDates(item);
			CStyle= GetExerciceStyle(item);
			Cstrike=GetStrikes(item);

			statusExercised = GetStringFromXMLNode(item, "Exercised");

			//First Period Rank 
			vCraPricing->Elt(0) = GetIntFromXMLNode(item, "ProdData/cBLOB_I/Num1");
			//Nb period Max 
			vCraPricing->Elt(1) = GetIntFromXMLNode(item, "ProdData/cBLOB_I/Num2");

			//A min  
			vCraPricing->Elt(2) = GetDoubleFromXMLNode(item, "ProdData/cBLOB_I/Rate1");
			//A max  
			vCraPricing->Elt(3) = GetDoubleFromXMLNode(item, "ProdData/cBLOB_I/Rate2");
			//Beta min  
			vCraPricing->Elt(4) = GetDoubleFromXMLNode(item, "ProdData/cBLOB_I/Rate3");
			//Beta max  
			vCraPricing->Elt(5) = GetDoubleFromXMLNode(item, "ProdData/cBLOB_I/Rate4");
			//Nb Step 
			vCraPricing->Elt(6) = GetIntFromXMLNode(item, "ProdData/cBLOB_I/Num3");
	
			if (item) item->Release();
		}

		if (resultList) resultList->Release();
		resultList=NULL;

		// Recherche de la prochaine NoticeDates
		int		compteur(0);
		int		vExerciseStartDatesSize = CStyle->GetExerciseStartDates()->GetSize();
		bool	vFound = false;

		if (!date.isMinDate())
		{
			while( (compteur < vExerciseStartDatesSize) && !vFound )
			{
				if( date.GetJulian() <  CStyle->GetExerciseStartDates()->Elt(compteur) )
					vFound = true;
				else
					compteur++;
			}
		}

		if(vFound)
		{
			// Changement des dates d'exercice...
			ARM_Vector* newNoticeDates = new ARM_Vector(CvExerciseDates->GetSize()-compteur);
			ARM_Vector* newExpiryDates = new ARM_Vector(CvExerciseDates->GetSize()-compteur);
			// ... et des fees !
			ARM_Vector* newFees = new ARM_Vector(Cstrike->GetDiscreteValues()->GetSize()-compteur);

			for (int i = 0; i < newNoticeDates->GetSize(); i++)
			{
				newNoticeDates->Elt(i) = CStyle->GetExerciseStartDates()->Elt(i+compteur);
				newExpiryDates->Elt(i) = CStyle->GetExerciseEndDates()->Elt(i+compteur);
				newFees->Elt(i) = Cstrike->GetDiscreteValues()->Elt(i+compteur);
			}
			CStyle->SetExerciseStartDates(newNoticeDates);
			CStyle->SetExerciseEndDates(newExpiryDates);
			CStyle->SetExerciseDatesSize(newNoticeDates->GetSize());
			CStyle->SetExerciseType(K_BERMUDAN);

			// les dates sont clonees, pas les valeurs
			Cstrike->SetDiscreteDates(newNoticeDates);
			Cstrike->SetDiscreteValues(newFees);

			delete newNoticeDates;
			delete newExpiryDates;
//			delete newFees;
		}
		else
			compteur--;

		if (isEtk == 1)
			tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Assets/ASSET";
		else
			tmpChaine = (CCString)SUMMITSTRING + (CCString)"SWAPTION/Assets/ASSET";

		if (xmlDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
		{
			resultList->get_length(&nbAssets);

			if (nbAssets != 2)
			{
				CCString msg((CCString)"Option XML string is not valid");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
			}

			// Récupération des 2 chaines XML correspondant à chaque asset, puis détermination de quel est le fix, quel est le float
			for (long indexAsset=0 ; indexAsset<nbAssets ; indexAsset++)
			{
				hr=resultList->get_item(indexAsset, &item);
				if (hr==S_OK && item!=NULL)
				{
					item->selectSingleNode(_bstr_t((const char *)"INTEREST_dmIndex"), &theNode);
					if (theNode!=NULL)
					{
						theNode->get_text(&resultat);
						theNode->Release();
						theNode = NULL;
						// on selectionne la patte CorridorLeg
						if (strcmp((const char *)_bstr_t(resultat),"FORM") == 0)	// Jambe Corridor
						{
							CstartDate = GetStartDate(item);
							CendDate = GetEndDate(item);
							CPorS = GetPorS(item);
							CpayFreq = GetPayFreq(item);
							CresetFreq = GetResetFreq(item);
							Ccurrency = GetCcy(item);
							CRefValNotional = GetNotional(item);
							CDaycount = GetAssetDayCount(item);
							Cstub = GetStubRule(item);
							GetPayCalendar(item,CpayCal);
							GetResetCalendar(item,CresetCal);
							GetCustom(item, Ccust);

							// take first forward start date
							char tmpDay[15];
							GetPayRollDate(item,tmpDay);
							int payDay = atoi(tmpDay);
							ARM_Date nextStartDate(CvExerciseDates->Elt(compteur));
							nextStartDate.UnadjustDate(payDay);

							if (CstartDate < nextStartDate)
							{
								CstartDate = nextStartDate;
								CstartDate.GoodBusinessDay(K_MOD_FOLLOWING, CpayCal);
							}
						}
						else
						// il s'agit de la patte swapleg
						{
							item->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
							if (theNode!=NULL)
							{
								theNode->get_text(&resultat);
								theNode->Release();
								theNode = NULL;
								
								SstartDate = GetStartDate(item);
								SendDate=GetEndDate(item);
								SPorS = GetPorS(item);
								SDaycount = GetAssetDayCount(item);
								SpayFreq = GetPayFreq(item);
								SpayTiming = GetPayTiming(item);
								Sstub = GetStubRule(item);
								Scurrency = GetCcy(item);
								Sspread=GetSpread(item)*100.0;
								GetPayCalendar(item,SpayCal);
								SintRule=GetIntRule(item);
								SnotionalEx=GetNotionalExchangeFlag(item);
								SRefValNotional= GetNotional(item);
								GetCustom(item, Scust);

								// take first forward start date
								char tmpDay[15];
								GetPayRollDate(item,tmpDay);
								int payDay = atoi(tmpDay);
								ARM_Date nextStartDate(CvExerciseDates->Elt(compteur));
								nextStartDate.UnadjustDate(payDay);

								if (SstartDate < nextStartDate)
								{
									SstartDate = nextStartDate;
									SstartDate.GoodBusinessDay(K_MOD_FOLLOWING, SpayCal);
								}

								// cas d'une jambe fixe
								if (strcmp((const char*)_bstr_t(resultat),"FIX") == 0)	
								{
									SIndexType=K_FIXED;
									SdecompFreq= GetDecompFreq(item);
									
									swapleg= new ARM_SwapLeg(SstartDate, 
															SendDate, Sspread, 
															SPorS, SpayFreq, 
															SDaycount, 
															SdecompFreq,
															SpayTiming,
															SintRule,
															Sstub,
															Scurrency,
															SpayCal,
															SnotionalEx,
															NULL,
															1);
								}
								else
								// cas d'une jambe flottante
								if (strcmp((const char*)_bstr_t(resultat),"FLO") == 0)	
								{
									SIndexType=GetIndexType(item);
									SresetFreq = GetResetFreq(item);
									SresetTiming = GetResetTiming(item);
									SresetGap=GetResetGap(item);
									SpayGap=GetPayGap(item);
									GetResetCalendar(item,SresetCal);
									SspreadVariable = GetSpreadVariable(item);

									swapleg= new ARM_SwapLeg(SstartDate, SendDate, 
															(ARM_INDEX_TYPE) SIndexType, 
															SPorS, 
															Sspread, 
															SresetFreq, SpayFreq, 
															SresetTiming, SpayTiming,
															Scurrency,
															SintRule,
															SresetGap,
															SresetCal,
															SpayCal,
															1,
															SnotionalEx,
															Sstub,
															NULL,
															1);

									if (SspreadVariable)
									{
										// Summit donne le spread variable indexé sur la startdate 
										// or on a besoin de la resetdate correspondante
										VECTOR<double> ResetDates;
										VECTOR<double> Spreads;
										ARM_SwapLeg* fictifswapleg= NULL;
										int size= SspreadVariable->size();
										if (size > 0)
										{
											if (mode_fictif)
											{	
												fictifswapleg=new ARM_SwapLeg(GetStartDate(item), SendDate, 
															(ARM_INDEX_TYPE) SIndexType, 
															SPorS, 
															Sspread, 
															SresetFreq, SpayFreq, 
															SresetTiming, SpayTiming,
															Scurrency,
															SintRule,
															SresetGap,
															SresetCal,
															SpayCal,
															1,
															SnotionalEx,
															Sstub,
															NULL,
															1);
											}
											else
												fictifswapleg = swapleg;

											ARM_Vector* tmpResetDates = fictifswapleg->GetResetDates();
											ARM_Vector* tmpStartDates = fictifswapleg->GetFlowStartDates();
											int istart=0;
											while ( SspreadVariable->GetDiscreteDates()->Elt(istart) <= tmpStartDates->Elt(0) )
											{
												istart++;
												if (istart==size)
													break;	
											}
											if (istart==0) istart++;

											for (int i = istart-1; i < size; i++)
											{
												int j = 0;
												while (tmpStartDates->Elt(j) < SspreadVariable->GetDiscreteDates()->Elt(i))
												j++;
												Spreads.push_back(SspreadVariable->GetDiscreteValues()->Elt(i));
												ResetDates.push_back(tmpResetDates->Elt(j));
											}
										}
										ARM_Vector* VResetDates= CreateARMVectorFromVECTOR(ResetDates);
										ARM_Vector* VSpreads= CreateARMVectorFromVECTOR(Spreads);
										ARM_ReferenceValue* newRefVal = new ARM_ReferenceValue(VResetDates, VSpreads,K_YIELD,1);
										newRefVal->SetCalcMethod(K_STEPUP_LEFT);

										swapleg->SetVariableSpread(newRefVal);
										delete newRefVal;
										delete SspreadVariable;
										if (mode_fictif)
											delete fictifswapleg;

									}
								}				
							}
						}
					}
				}
				if (item) item->Release();
			}

			ARM_IRIndex* CPaymentIndex=NULL; 
			// cas d'un indice fixe
			if (CPayIndexType == K_FIXED)
			{
				CPaymentIndex = new ARM_IRIndex(Ccurrency->GetCcyName(), CDaycount);
				vCraPricing->Elt(7) = 0;
			}
			else // cas d'un indice type Libor
			{
				CPaymentIndex= new ARM_IRIndex((ARM_INDEX_TYPE) CPayIndexType,K_DEF_FREQ,K_DEF_FREQ,Ccurrency,CDaycount);
				vCraPricing->Elt(7) = 1;
			}

			ARM_IRIndex* CRefIndex =NULL;
			// on ne prend pas en compte de CDaycount dans l'indice de reference!!!!!!!!
			// cas d'un indice fixe
			if (CRefIndexType == K_FIXED)
				CRefIndex= new ARM_IRIndex(Ccurrency->GetCcyName(), CDaycount);
			else // cas d'un indice type Libor
				CRefIndex= new ARM_IRIndex((ARM_INDEX_TYPE) CRefIndexType,K_DEF_FREQ,K_DEF_FREQ,Ccurrency, CDaycount);

			corridor= new ARM_CorridorLeg(CstartDate, CendDate,
										CPorS, 
										CPaymentIndex,
										CpayFreq,
										CRefValSpread,
										CRefIndex, CresetFreq,
										CresetPayTiming, CresetRefTiming,
										Cstub,
										CRefBarDown, K_STD,
										CRefBarUp, K_STD,
										Ccurrency,
										K_DEF_FREQ, K_LINEAR,
										K_DIGITALE,
										CdecompPricingFlag,
										CresetCal,
										CpayCal);

			if (CPaymentIndex)
				delete CPaymentIndex;
			CPaymentIndex=NULL;

			if (CRefIndex)
				delete CRefIndex;
			CRefIndex=NULL;

			if (CRefValSpread)
				delete CRefValSpread;
			CRefValSpread=NULL;

			if (CRefBarDown)
				delete CRefBarDown;
			CRefBarDown=NULL;

			if (CRefBarUp)
				delete CRefBarUp;
			CRefBarUp=NULL;

			if (Ccurrency)
				delete Ccurrency;
			Ccurrency = NULL;
		}

		if (resultList) 
			resultList->Release();
		resultList = NULL;

		if (xmlDoc) 
			xmlDoc->Release();
		xmlDoc = NULL;

        if (CvExerciseDates)
		   delete CvExerciseDates;

		corridor->SetAmount(CRefValNotional);

		if (CRefValNotional)
			delete CRefValNotional;
		CRefValNotional = NULL;

		swapleg->SetAmount(SRefValNotional);

		if (SRefValNotional)
			delete SRefValNotional;
		SRefValNotional = NULL;

		int nbsec=2;
		double* weights = new double[nbsec];
		double* mktprices = new double[nbsec];

		weights[0] = 1.0;
		weights[1] = 1.0;
		mktprices[0] = 1.0;
		mktprices[1] = 1.;
				
		ARM_StdPortfolio* Portfolio= new ARM_StdPortfolio(nbsec, weights, mktprices);

		// les ARM_Security ne sont pas clonés:
		Portfolio->SetAsset((ARM_Security *) corridor, 0);
		Portfolio->SetAsset((ARM_Security *) swapleg, 1);

		// verif qu'il existe au moins une option sur le portefeuille
		if (CStyle->GetExerciseStartDates()->GetSize()>0)
		{
			security = new ARM_OptionPortfolio(Portfolio,
											   CStyle,
											   Cstrike,
											   K_CALL);
			security->SetPorS(PorS);

			if ( (status == "TERM") || (statusExercised == "EXER") )
				security->SetFullterminated(true);

			if (CStyle->GetExerciseStartDates()->Elt(CStyle->GetExerciseStartDates()->GetSize() - 1) <= date.GetJulian() )
				security->SetFullterminated(true);

			// TMP : Pour hedge d'un ARM_Object créé sous Excel
			// Il faudra voir si ces paramètres correspondent à ceux de ARM_OptionPortfolio::SetPFCalibParams()
			security->SetCraPricing((ARM_Vector*)vCraPricing);

			security->SetAmount(corridor->GetAmount());
		}

		security->SetCurrencyUnit(swapleg->GetCurrencyUnit());
	
		if (vCraPricing != NULL)
			*CraPricing = *vCraPricing;
		delete vCraPricing;
		vCraPricing=NULL;
		
		if (Cstrike)
			delete Cstrike;
		Cstrike = NULL;

		if (CStyle)
			delete CStyle;
		CStyle = NULL;

		if (Portfolio)
			delete Portfolio;
		Portfolio = NULL;

		delete [] weights;
		weights = NULL;

		delete [] mktprices;
		mktprices = NULL;
	}
	
	catch (Exception &e)
	{
		if (resultat) 
           SysFreeString(resultat);
		
        if (theNode) 
           theNode->Release();
		
        if (resultList) 
           resultList->Release();
		
        if (security) 
           delete security;
		
        if (xmlDoc) 
           xmlDoc->Release(); xmlDoc = NULL;

		try
		{
			throw e;
		}

		catch (...)
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, e.GetMessage().c_str());
		}
	}
	
	if ( (strcmp(Scust,"CUST") == 0) || (strcmp(Ccust,"CUST") == 0) )
	{
		if (security)
			delete security;
		security = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Customized CRAs not yet implemented");
	}

	return security;
}


// return FX_OPT_STRIP doc
MSXML2::IXMLDOMDocument* searchDocForExoticFxOptStrip(vector<string> listeTrade)
{
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMDocument* XMLDocForFxOptStrip = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	CCString tradeType;
	int trouveFXOPTSTRIP = 0;
	int deleteXMLDoc;

	int nbExotic = 0;

	HRESULT hr;
	VARIANT_BOOL bOK;
	long nbNodes;

	for (int i = 0; i < listeTrade.size (); i++)
	{
		try
		{
			hr = CoInitialize(NULL); 

			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = (listeTrade[i]).c_str();
			_bstr_t tmpChaine; 
			VariantTools::convert(listeTrade[i],tmpChaine); 


			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating XML document for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		deleteXMLDoc = 1;

		XMLDoc->selectSingleNode(_bstr_t((const char *)"*/*/TradeType"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			tradeType = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		if (strcmp(tradeType,"EXOTIC") == 0)
		{
			nbExotic ++;

			if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
					hr = S_FALSE;

					CCString msg((CCString)"Invalid XML string for an Exotic in getting PRCS (0 Asset)");
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}
				else
				{
					for (int j = 0; j < nbNodes; j++)
					{
						hr=resultList->get_item(j, &listItem);
						listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);
						if (theNode!=NULL)
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							if (strncmp((const char*)ff1,"FX_OPT_STRIP",12) == 0)
							{
								XMLDocForFxOptStrip = XMLDoc;
								trouveFXOPTSTRIP = 1;
								deleteXMLDoc = 0;
								i = listeTrade.size ();
							}

							theNode->Release();
							theNode=NULL;
							if (resultat) SysFreeString(resultat);
						}

						listItem->Release();
						listItem=NULL;
					}
				}
				resultList->Release();
				resultList=NULL;
			}
		}
		if (deleteXMLDoc == 1)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}
	}

	return XMLDocForFxOptStrip;
}



// return FX_OPT_STRIP doc
MSXML2::IXMLDOMDocument* searchDocForExoticFunding(vector<string> listeTrade)
{
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	vector<MSXML2::IXMLDOMDocument*> vectFunding;

	CCString tradeType;

	int trouveFUNDINGtmp = 0;

	int nbExotic = 0;

	HRESULT hr;
	VARIANT_BOOL bOK;
	long nbNodes;

	for (int i = 0; i < listeTrade.size (); i++)
	{
		trouveFUNDINGtmp = 0;

		try
		{
			hr = CoInitialize(NULL); 

			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = (listeTrade[i]).c_str();
			_bstr_t tmpChaine; 
			VariantTools::convert(listeTrade[i],tmpChaine); 


			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating XML document for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		XMLDoc->selectSingleNode(_bstr_t((const char *)"*/*/TradeType"), &theNode);
		if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			tradeType = ff1;

			theNode->Release();
			theNode = NULL;
			
            if (resultat) 
               SysFreeString(resultat);
		}

		if (strcmp(tradeType,"EXOTIC") == 0)
		{
			nbExotic ++;

			if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
				{
					hr = S_FALSE;

					CCString msg((CCString)"Invalid XML string for an Exotic in getting PRCS (0 Asset)");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}
				else
				{
					for (int j = 0; j < nbNodes; j++)
					{
						hr=resultList->get_item(j, &listItem);
						listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);
						if (theNode!=NULL)
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							if (strncmp((const char*)ff1,"FX_OPT_STRIP",12) == 0)
							{
								trouveFUNDINGtmp = 0;
								j = nbNodes;
							}
							else
							{
								trouveFUNDINGtmp = 1;
							}

							theNode->Release();
							theNode=NULL;
							if (resultat) SysFreeString(resultat);
						}
						else
						{
							trouveFUNDINGtmp = 1;
						}

						listItem->Release();
						listItem=NULL;
					}
				}
				resultList->Release();
				resultList=NULL;
			}
		}
		if (trouveFUNDINGtmp == 0)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}
		else
		{
			vectFunding.push_back(XMLDoc);
		}
	}

	int indexPay;
	int indexRecSameCcyAsPay;
	CCString ccy;

	if ( vectFunding.size() == 1 )
	{
	   return vectFunding[0];
	}
	else if (vectFunding.size() == 3)
	{
		// cas x-ccy
		for (int i = 0; i < vectFunding.size(); i++)
		{
			if (vectFunding[i]->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				hr=resultList->get_item(0, &listItem);
				listItem->selectSingleNode(_bstr_t((const char *)"PorS"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					if (strcmp((const char*)ff1,"S") == 0)
					{
						indexPay = i;
					}

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				if (indexPay == i)
				{
					listItem->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						ccy = ff1;

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					i = vectFunding.size();
				}

				listItem->Release();
				listItem = NULL;
			}

			resultList->Release();
			resultList = NULL;
		}
		for (i = 0; i < vectFunding.size(); i++)
		{
			if (i != indexPay)
			{
				if (vectFunding[i]->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
				{
					resultList->get_length(&nbNodes);

					hr=resultList->get_item(0, &listItem);
					listItem->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp((const char*)ff1,ccy) != 0)
						{
							indexRecSameCcyAsPay = i;
							i = vectFunding.size();
						}

						theNode->Release();
						theNode = NULL;
						
                        if (resultat) 
                           SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;
				}

                resultList->Release();
				resultList = NULL;
			}
		}
		for (i = 0; i < vectFunding.size(); i++)
		{
			if (i != indexRecSameCcyAsPay)
			{
				vectFunding[i]->Release();
				vectFunding[i] = NULL;
			}
		}
		return vectFunding[indexRecSameCcyAsPay];
	}
	else
	{
		for (int i = 0; i < vectFunding.size(); i++)
		{
			vectFunding[i]->Release();
			vectFunding[i] = NULL;
		}

		CCString msg((CCString)"funding : not yet implemented (just for 1 or 3 exotics)");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



MSXML2::IXMLDOMDocument* searchDocForSwaption(vector<string> listeTrade)
{
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	CCString tradeType;
	int trouve = 0;
	HRESULT hr;
	VARIANT_BOOL bOK;
	long nbNodes;

	VECTOR<MSXML2::IXMLDOMDocument*> listDoc;

	for (int i = 0; i < listeTrade.size (); i++)
	{
		try
		{
			hr = CoInitialize(NULL); 

			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = (listeTrade[i]).c_str();
			_bstr_t tmpChaine; 
			VariantTools::convert(listeTrade[i],tmpChaine); 


			XMLDoc->loadXML(tmpChaine, &bOK);
		}

		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;

			CCString msg((CCString)"Pb in creating XML document for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		XMLDoc->selectSingleNode(_bstr_t((const char *)"*/*/TradeType"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			tradeType = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		if (strcmp(tradeType,"SWOPT") == 0)
		{
			listDoc.push_back(XMLDoc);

			if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Env/ENV")), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if ( nbNodes != 1 )
				{
					//if (XMLDoc) XMLDoc->Release();
					//if (resultList) resultList->Release(); //Non release fait dans le catch
					hr = S_FALSE;

					CCString msg((CCString)"Invalid XML string for the dummy Swaption in getting PRCS (not 1 option)");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (const char *) msg);
				}
				else
				{
					hr = resultList->get_item(0, &listItem);
					listItem->selectSingleNode(_bstr_t((const char *)"ProductType"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp((const char*)ff1,"DSWAPTION") == 0)
						{
							trouve = 1;
						}

						theNode->Release();
						theNode = NULL;

                        if (resultat) 
                           SysFreeString(resultat);
					}

					listItem->Release();
					listItem = NULL;
				}

				resultList->Release();
				resultList = NULL;
			}
		}
		else
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}


		if ( trouve == 1 )
		{
			for (i = 0; i < listDoc.size()-1; i++)
			{
				listDoc[i]->Release();
				listDoc[i] = NULL;
			}

			return(XMLDoc);
		}
	}

	if ( listDoc.size() != 1 )
	{
		CCString msg((CCString)"Invalid XML string for the dummy Swaption in getting PRCS (more than 2 swaptions in the structure)");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	return(listDoc[0]);
}




ARM_Vector* getNotAndCancDates(MSXML2::IXMLDOMDocument *XMLDoc, ARM_Vector** vCancelDates, CCString& custId, CCString& dealId)
{
	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	CCString tradeType;
	int trouve = 0;
	HRESULT hr;
	long nbNodes;
	long nbNoticeDates;
	ARM_Vector* vNoticeDates;

	// Récupération du customer id
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Env/ENV")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes != 1 )
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting cancellation dates");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		hr = resultList->get_item(0, &listItem);

		listItem->selectSingleNode(_bstr_t((const char *)"Cust"), &theNode);
		if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			custId = ff1;

			theNode->Release();
			theNode = NULL;

			if (resultat) 
               SysFreeString(resultat);
		}
		else
			custId = "";


		listItem->selectSingleNode(_bstr_t((const char *)"DealId"), &theNode);
		if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			dealId = ff1;

			theNode->Release();
			theNode = NULL;
			if (resultat) 
               SysFreeString(resultat);
		}
		else
			dealId = "";
	}

	if (theNode)
	{
		theNode->Release();
		theNode = NULL;
	}

	if (resultList)
	{
		resultList->Release();
		resultList = NULL;
	}

	if (listItem)
	{
		listItem->Release();
		listItem = NULL;
	}

	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Option/OPTION")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes != 1 )
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting cancellation dates");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		hr=resultList->get_item(0, &listItem);

		// Cancellation Dates number
		listItem->selectSingleNode(_bstr_t((const char *)"NumOpEve"), &theNode);
		if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			nbNoticeDates = atol((const char*) ff1);

			theNode->Release();
			theNode = NULL;
			
            if (resultat) 
               SysFreeString(resultat);
		}

		*vCancelDates = new ARM_Vector(nbNoticeDates);
		vNoticeDates = new ARM_Vector(nbNoticeDates);

		listItem->selectNodes(_bstr_t((const char *)"OpEvents/OPEVENT"), &resultList2);

		for (int i=0; i<nbNoticeDates;i++)
		{
			resultList2->get_item(i,&listItem);

			(*vCancelDates)->Elt(i) = GetDateFromXMLNode(listItem,"Date").GetJulian();

			vNoticeDates->Elt(i) = GetDateFromXMLNode(listItem,"ADate").GetJulian();

			listItem->Release();
			listItem = NULL;
		}

		resultList2->Release();
		resultList2 = NULL;
	}

	resultList->Release();
	resultList = NULL;

	return(vNoticeDates);
}



ARM_SwapLeg* getInitPeriodLeg(MSXML2::IXMLDOMDocument* XMLDoc,
							  int typeDeal)  // si 0 : structure, si 1 monoTrade
{
	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * listItem2 = NULL, * theNode = NULL, * theNode2 = NULL;

	HRESULT hr;

	ARM_SwapLeg* initPeriodLeg = NULL;

	ARM_Date startDate;
	ARM_Date endDate;
	
	ARM_Date fstCpnEffDate;

	double fixCoupon;
	int fixFrequency;
	int daycountBasis;
	ARM_Currency* ccy = NULL;
	ARM_ReferenceValue* amount = NULL;
	double notional;
	char payCal[7];
	char myRefDate[11];
	long NxId;
	long intRuleId;

	VECTOR<double> ntlVal;
	VECTOR<double> ntlDate;

	long nbNodes;

	if (typeDeal == 0)
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes > 2)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Dont know how to use the exotic with more than 2 assets");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			if (nbNodes == 1)
				return NULL;

			for (int i = 0; i < nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);
				listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);

				if ( theNode == NULL )
				{
					// start of initPeriodLeg
					startDate = GetStartDate(listItem);
					fstCpnEffDate = startDate;
					ntlDate.push_back(startDate.GetJulian());

					// end of initPeriodLeg
					endDate = GetEndDate(listItem);

					// fix coupon of initPeriodLeg
					fixCoupon = GetDoubleFromXMLNode(listItem,"INTEREST_Rate") * 100.;

					// fix frequency of initPeriodLeg
					fixFrequency = GetPayFreq(listItem);

					// daycount Basis of initPeriodLeg
					listItem->selectSingleNode(_bstr_t((const char *)"INT_ACC_Basis"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						daycountBasis = FromSummitDaycountToARMDaycount((const char*) ff1);

						theNode->Release();
						theNode = NULL;
						
                        if (resultat) 
                           SysFreeString(resultat);
					}

					// int rule of initPeriodLeg
					listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_IntRule"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp((const char*) ff1,"N") == 0)
							intRuleId = K_UNADJUSTED;
						else
							intRuleId = K_ADJUSTED;

						theNode->Release();
						theNode = NULL;

						if (resultat) 
                           SysFreeString(resultat);
					}

					// currency of initPeriodLeg
					listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_Ccy"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						ccy = new ARM_Currency(ff1);

						theNode->Release();
						theNode = NULL;
						if (resultat) 
                           SysFreeString(resultat);
					}

					// notional of initPeriodLeg
					listItem->selectSingleNode(_bstr_t((const char *)"Notional"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						notional = atof((const char*) ff1);
						ntlVal.push_back(notional);

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}

					// paycal of initPeriodLeg
					GetPayCalendar(listItem,payCal);

					// reference Date
					listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp((const char*)ff1,"") == 0)
							strcpy(myRefDate,"NULL");
						else
						{
							ARM_Date tmpDate (ff1,"YYYYMMDD");
							tmpDate.JulianToStrDate(myRefDate);
						}

						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);
					}
					else
						strcpy(myRefDate,"NULL");

					// Echange de Notional
					NxId = GetNotionalExchangeFlag(listItem);

					// First Coupon Effective Date
					listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (strcmp((const char*)ff1,"") != 0)
						{
							ARM_Date tmpDate (ff1,"YYYYMMDD");

							fstCpnEffDate = tmpDate;
						}

						theNode->Release();
						theNode = NULL;
						if (resultat) 
                           SysFreeString(resultat);
					}

					// Recuperation des evenements
					listItem->selectSingleNode(_bstr_t((const char *)"Events"), &theNode);
					if (theNode!=NULL)
					{
						if (theNode->selectNodes(_bstr_t((const char *)("EVENT")), &resultList2) == S_OK)
						{
							long nbEvents;

							resultList2->get_length(&nbEvents);

							CCString type;

							for (int i = 0; i < nbEvents; i++)
							{
								hr=resultList2->get_item(i, &listItem2);
								
								listItem2->selectSingleNode(_bstr_t((const char *)"Type"), &theNode2);
								if (theNode2!=NULL)
								{
									BSTR resultat = NULL;
									theNode2->get_text(&resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									type = (const char*) ff1;

									theNode2->Release();
									theNode2=NULL;
									if (resultat) SysFreeString(resultat);
								}

								if (strcmp(type,"NTL") == 0)
								{
									int trouveNewNTL = 0;

									ARM_Date tmpDate = GetDateFromXMLNode(listItem2,"ADate");
									if (tmpDate != endDate)
									{
										ntlDate.push_back(tmpDate.GetJulian());
										trouveNewNTL = 1;
									}

									listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
									if (theNode2!=NULL)
									{
										BSTR resultat = NULL;
										theNode2->get_text(&resultat);

										_bstr_t ff(resultat,false);
										char * ff1=(char *)ff;

										double ntlTmp = atof(ff1);
										if (trouveNewNTL == 1)
											ntlVal.push_back(ntlTmp);

										theNode2->Release();
										theNode2 = NULL;
										
                                        if (resultat) 
                                           SysFreeString(resultat);
									}
								}

								listItem2->Release();
								listItem2 = NULL;
							}

							resultList2->Release();
							resultList2 = NULL;
						}
					}
				}
				else
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					theNode->Release();
					theNode=NULL;

					if (strncmp((const char*)ff1,"FX_OPT_STRIP",12) != 0)
					{
						if (resultat) SysFreeString(resultat);

						// start of initPeriodLeg
						startDate = GetStartDate(listItem);
						fstCpnEffDate = startDate;
						ntlDate.push_back(startDate.GetJulian());

						// end of initPeriodLeg
						endDate = GetEndDate(listItem);

						// fix coupon of initPeriodLeg
						fixCoupon = GetDoubleFromXMLNode(listItem,"INTEREST_Rate") * 100.;

						// fix frequency of initPeriodLeg
						fixFrequency = GetPayFreq(listItem);

						// daycount Basis of initPeriodLeg
						listItem->selectSingleNode(_bstr_t((const char *)"INT_ACC_Basis"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							daycountBasis = FromSummitDaycountToARMDaycount((const char*) ff1);

							theNode->Release();
							theNode = NULL;
							
                            if (resultat) 
                               SysFreeString(resultat);
						}

						// int rule of initPeriodLeg
						listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_IntRule"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							if (strcmp((const char*) ff1,"N") == 0)
								intRuleId = K_UNADJUSTED;
							else
								intRuleId = K_ADJUSTED;

							theNode->Release();
							theNode = NULL;
							
                            if (resultat) 
                               SysFreeString(resultat);
						}

						// currency of initPeriodLeg
						listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_Ccy"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							ccy = new ARM_Currency(ff1);

							theNode->Release();
							theNode = NULL;

							if (resultat) 
                               SysFreeString(resultat);
						}

						// notional of initPeriodLeg
						listItem->selectSingleNode(_bstr_t((const char *)"Notional"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							notional = atof((const char*) ff1);
							ntlVal.push_back(notional);

							theNode->Release();
							theNode = NULL;

							if (resultat) 
                               SysFreeString(resultat);
						}

						// paycal of initPeriodLeg
						GetPayCalendar(listItem,payCal);

						// reference Date
						listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							if (strcmp((const char*)ff1,"") == 0)
								strcpy(myRefDate,"NULL");
							else
							{
								ARM_Date tmpDate (ff1,"YYYYMMDD");
								tmpDate.JulianToStrDate(myRefDate);
							}

							theNode->Release();
							theNode = NULL;
							
                            if (resultat) 
                               SysFreeString(resultat);
						}
						else
							strcpy(myRefDate,"NULL");

						// Echange de Notional
						NxId = GetNotionalExchangeFlag(listItem);

						// First Coupon Effective Date
						listItem->selectSingleNode(_bstr_t((const char *)"FstCpEffect"), &theNode);
						if ( theNode != NULL )
						{
							BSTR resultat = NULL;
							theNode->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							if (strcmp((const char*)ff1,"") != 0)
							{
								ARM_Date tmpDate (ff1,"YYYYMMDD");

								fstCpnEffDate = tmpDate;
							}

							theNode->Release();
							theNode = NULL;
							
                            if (resultat) 
                               SysFreeString(resultat);
						}

						// Recuperation des evenements
						listItem->selectSingleNode(_bstr_t((const char *)"Events"), &theNode);
						if (theNode!=NULL)
						{
							if (theNode->selectNodes(_bstr_t((const char *)("EVENT")), &resultList2) == S_OK)
							{
								long nbEvents;

								resultList2->get_length(&nbEvents);

								CCString type;

								for (int i = 0; i < nbEvents; i++)
								{
									hr=resultList2->get_item(i, &listItem2);
									
									listItem2->selectSingleNode(_bstr_t((const char *)"Type"), &theNode2);
									if (theNode2!=NULL)
									{
										BSTR resultat = NULL;
										theNode2->get_text(&resultat);

										_bstr_t ff(resultat,false);
										char * ff1=(char *)ff;

										type = (const char*) ff1;

										theNode2->Release();
										theNode2 = NULL;
										if (resultat) 
                                           SysFreeString(resultat);
									}

									if (strcmp(type,"NTL") == 0)
									{
										int trouveNewNTL = 0;

										ARM_Date tmpDate = GetDateFromXMLNode(listItem2,"ADate");
										if (tmpDate != endDate)
										{
											ntlDate.push_back(tmpDate.GetJulian());
											trouveNewNTL = 1;
										}

										listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
										if (theNode2!=NULL)
										{
											BSTR resultat = NULL;
											theNode2->get_text(&resultat);

											_bstr_t ff(resultat,false);
											char * ff1=(char *)ff;

											double ntlTmp = atof(ff1);
											if (trouveNewNTL == 1)
												ntlVal.push_back(ntlTmp);

											theNode2->Release();
											theNode2 = NULL;
											
                                            if (resultat) 
                                               SysFreeString(resultat);
										}
									}

									listItem2->Release();
									listItem2 = NULL;
								}
				
                                resultList2->Release();
								resultList2 = NULL;
							}					
						}
					}
				}

				listItem->Release();
				listItem=NULL;
			}
		}

		resultList->Release();
		resultList = NULL;
	}
	else
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Option/OPTION/ProdData/cBLOB_I")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if ( nbNodes != 1 )
			{
				resultList->Release();
				resultList = NULL;

				return(NULL);
			}

			hr = resultList->get_item(0, &listItem);

			startDate = GetDateFromXMLNode(listItem,"Date1");
			fstCpnEffDate = startDate;
			ntlDate.push_back(startDate.GetJulian());

			endDate = GetDateFromXMLNode(listItem,"Date2");

			if (startDate == endDate)
			{
				listItem->Release();
				listItem=NULL;

				resultList->Release();
				resultList=NULL;

				return NULL;
			}

			// fix coupon of initPeriodLeg
			fixCoupon = GetDoubleFromXMLNode(listItem,"Rate1") * 100.;

			// daycount Basis of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"Basis1"), &theNode);
			if ( theNode != NULL )
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				daycountBasis = FromSummitDaycountToARMDaycount((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) 
                   SysFreeString(resultat);
			}

			// fix frequency of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"Num3"), &theNode);
			if ( theNode != NULL )
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				fixFrequency = FromSummitFreqToARMFreq((const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// int rule of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"Num6"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*) ff1,"2") == 0)
					intRuleId = K_UNADJUSTED;
				else
					intRuleId = K_ADJUSTED;

				theNode->Release();
				theNode = NULL;
				if (resultat) 
                   SysFreeString(resultat);
			}

			// currency of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"ResetCal1"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				ccy = new ARM_Currency(ff1);

				theNode->Release();
				theNode = NULL;
				if (resultat) 
                   SysFreeString(resultat);
			}

			// notional of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"Amount3"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				notional = atof((const char*) ff1);
				ntlVal.push_back(notional);

				theNode->Release();
				theNode = NULL;
				if (resultat) 
                   SysFreeString(resultat);
			}

			// paycal of initPeriodLeg
			listItem->selectSingleNode(_bstr_t((const char *)"ResetCal2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				strcpy(payCal,(const char*) ff1);

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// reference Date
			listItem->selectSingleNode(_bstr_t((const char *)"Date3"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") == 0)
					strcpy(myRefDate,"NULL");
				else
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");
					tmpDate.JulianToStrDate(myRefDate);
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}
			else
				strcpy(myRefDate,"NULL");

			// Echange de Notional
			listItem->selectSingleNode(_bstr_t((const char *)"Num2"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*) ff1,"3") == 0)
					NxId = K_NX_BOTH;
				else if (strcmp((const char*) ff1,"1") == 0)
					NxId = K_NX_START;
				else if (strcmp((const char*) ff1,"2") == 0)
					NxId = K_NX_END;
				else
					NxId = K_NX_NONE;

				theNode->Release();
				theNode = NULL;
				
                if (resultat) 
                   SysFreeString(resultat);
			}
			else
				NxId = K_NX_NONE;

			// First Coupon Effective Date
			listItem->selectSingleNode(_bstr_t((const char *)"Date4"), &theNode);
			if ( theNode != NULL )
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				if (strcmp((const char*)ff1,"") != 0)
				{
					ARM_Date tmpDate (ff1,"YYYYMMDD");

					fstCpnEffDate = tmpDate;
				}

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			// Recuperation des evenements
			listItem->selectSingleNode(_bstr_t((const char *)"Schedule1"), &theNode);
			if (theNode!=NULL)
			{
				if (theNode->selectNodes(_bstr_t((const char *)("cBLBLST_I")), &resultList2) == S_OK)
				{
					long nbEvents;

					resultList2->get_length(&nbEvents);

					if (nbEvents > 0)
					{
						ntlDate.clear();
						ntlVal.clear();
					}

					for (int i = 0; i < nbEvents; i++)
					{
						hr=resultList2->get_item(i, &listItem2);

						ARM_Date tmpDate = GetDateFromXMLNode(listItem2,"StartDate");
						ntlDate.push_back(tmpDate.GetJulian());

						listItem2->selectSingleNode(_bstr_t((const char *)"Amount1"), &theNode2);
						if ( theNode2 != NULL )
						{
							BSTR resultat = NULL;
							theNode2->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							double ntlTmp = atof(ff1);
							ntlVal.push_back(ntlTmp);

							theNode2->Release();
							theNode2 = NULL;

							if (resultat)
                               SysFreeString(resultat);
						}

						listItem2->Release();
						listItem2 = NULL;
					}
				}					
				resultList2->Release();
				resultList2 = NULL;
			}
			theNode->Release();
			theNode=NULL;
		}
		listItem->Release();
		listItem=NULL;

		resultList->Release();
		resultList=NULL;
	}

	long stubRuleId = K_SHORTSTART;

	initPeriodLeg = new ARM_SwapLeg(startDate, endDate, fixCoupon, K_RCV,
									fixFrequency, daycountBasis,
                                    K_COMP_PROP, K_ARREARS,
                                    intRuleId,
									stubRuleId, 
                                    ccy,
                                    payCal,
                                    NxId,
                                    myRefDate);

	ARM_Vector* vNtlDate;
	ARM_Vector* vNtlVal;
	long size = ntlDate.size();

	if ( size == 0 )
	{
		vNtlDate = NULL;
		vNtlVal = NULL;
	}
	else
	{
		vNtlDate = new ARM_Vector(size);
		vNtlVal = new ARM_Vector(size);
		for (int i = 0; i < size; i++)
		{
			vNtlDate->Elt(i) = ntlDate[i];
			vNtlVal->Elt(i) = ntlVal[i];
		}
	}
	
	vNtlDate->Elt(0) = initPeriodLeg->GetPaymentDates()->Elt(0);
	
	int j = 1;
	for (int i = 1; i < vNtlDate->GetSize(); i++)
	{
		while (( j < initPeriodLeg->GetFlowStartDates()->GetSize() ) 
               && 
               ( vNtlDate->Elt(i) < initPeriodLeg->GetFlowStartDates()->Elt(j) ))
			j++;

		if (j < initPeriodLeg->GetFlowStartDates()->GetSize())
		{
			if (vNtlDate->Elt(i) == initPeriodLeg->GetFlowStartDates()->Elt(j))
				vNtlDate->Elt(i) = initPeriodLeg->GetPaymentDates()->Elt(j);
		}
	}

	amount = new ARM_ReferenceValue(vNtlDate,vNtlVal);
	amount->SetCalcMethod(K_STEPUP_RIGHT);

	initPeriodLeg->SetAmount(amount);

	if ( fstCpnEffDate != startDate )
	   initPeriodLeg->CustomizeFirstPeriod(fstCpnEffDate);

	delete amount;
	delete ccy;

	return(initPeriodLeg);
}



// Return Fx Numeraire Leg
ARM_SwapLeg* getStructuredLegs(MSXML2::IXMLDOMDocument* XMLDoc,
							   ARM_SwapLeg** FxUnderLeg,
							   ARM_ReferenceValue** rCap,
							   ARM_ReferenceValue** rFloor,
							   ARM_ReferenceValue** rFx,
							   ARM_ReferenceValue** FxFixing,
							   int* dualFlag,
							   double* dualStrike,
							   CCString& bookName,
							   int typeDeal)  // si 0 : structure, si 1 monoTrade
{
	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL, * resultListRefValue = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * listItem2 = NULL, * theNode = NULL, * theNode2 = NULL;
	HRESULT hr;

	ARM_SwapLeg* FxNumLeg = NULL;

	ARM_Vector* vTmpCap = NULL;
	ARM_Vector* vTmpFloor = NULL;
	ARM_Vector* vTmpFx = NULL;
	ARM_Vector* vTmpCpnRate1 = NULL;
	ARM_Vector* vTmpCpnRate2 = NULL;

	ARM_Date startDate;
	ARM_Date endDate;
	int daycountBasis;
	int fixFrequency;
	int resetTiming;
	int payTiming;
	int resetGap;
	int payGap;
	double notional;
	char payCal[4];
	char resetCal[4];
	int payAnnDay;
	long NxId;
	int intrule;

	char stubDate[11];
	int stubType;

	ARM_ReferenceValue* amount = NULL;
	ARM_ReferenceValue* rFxUnder = NULL;

	ARM_Currency* ccy1 = NULL;
	ARM_Currency* ccy2 = NULL;

	ARM_IRIndex* FxNumIndex = NULL;
	ARM_IRIndex* FxUnderIndex = NULL;

	long nbNodes;
	long nbNodes2;
	long nbEvents;
	long nbStrip;

	VECTOR<double> fxFixing;
	VECTOR<double> resetFixing;

	CCString tmpChaine;

	if ( typeDeal == 0 )
		tmpChaine = "Response/EXOTIC/Env/ENV";
	else
		tmpChaine = "Response/SWAPTION/Env/ENV";

	if (XMLDoc->selectNodes(_bstr_t((const char *)tmpChaine), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes != 1 )
		{
			hr = S_FALSE;

			CCString msg((CCString)"No ENV tag in PRDC asset");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							(const char *) msg);
		}

		hr = resultList->get_item(0, &listItem);
		listItem->selectSingleNode(_bstr_t((const char *)"Book"), &theNode);
		if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			bookName = ff1;

			theNode->Release();
			theNode = NULL;

			if (resultat) 
               SysFreeString(resultat);
		}
	}

	resultList->Release();
	resultList = NULL;

	if ( typeDeal == 0 )
		tmpChaine = "Response/EXOTIC/Assets/ASSET";
	else
		tmpChaine = "Response/SWAPTION/Assets/ASSET";

	if (XMLDoc->selectNodes(_bstr_t((const char *)tmpChaine), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (nbNodes > 2)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Dont know how to use the exotic with more than 2 assets");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							(const char *) msg);
		}

		for (int i = 0; i < nbNodes; i++)
		{
			hr=resultList->get_item(i, &listItem);
			listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);
			if ( theNode != NULL )
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				theNode->Release();
				theNode = NULL;

				if ( strncmp((const char*)ff1,"FX_OPT_STRIP",12) == 0 )
				{
					if (strcmp((const char*)ff1,"FX_OPT_STRIP_DUALE_SPOT") == 0)
						*dualFlag = 2;
					else if (strcmp((const char*)ff1,"FX_OPT_STRIP_DUALE_MANDAT") == 0)
						*dualFlag = 1;

					if (resultat) SysFreeString(resultat);

					// start of legs
					startDate = GetStartDate(listItem);

					// end of legs
					endDate = GetEndDate(listItem);

					GetStubDate1(listItem,stubDate);

					stubType = GetStubRule(listItem);

					// Pay annDay
					listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_AnnDay"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						payAnnDay = atoi((const char*) ff1);

						theNode->Release();
						theNode = NULL;

						if (resultat) 
                           SysFreeString(resultat);
					}

					// notional of num leg
					listItem->selectSingleNode(_bstr_t((const char *)"Notional"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						notional = atof((const char*) ff1);

						theNode->Release();
						theNode = NULL;

						if (resultat) 
                           SysFreeString(resultat);
					}

					// daycount Basis of legs
					listItem->selectSingleNode(_bstr_t((const char *)"INT_ACC_Basis"), &theNode);
					if ( theNode != NULL )
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						daycountBasis = FromSummitDaycountToARMDaycount((const char*) ff1);

						theNode->Release();
						theNode = NULL;

						if (resultat) 
                           SysFreeString(resultat);
					}

					// frequency of Legs
					fixFrequency = GetPayFreq(listItem);

					// reset Timing
					resetTiming = GetResetTiming(listItem);

					// pay Timing
					payTiming = GetPayTiming(listItem);

					// reset FX Gap
					resetGap = GetResetGap(listItem);

					// pay Gap
					payGap = GetPayGap(listItem);

					// paycal of structured legs
					GetPayCalendar(listItem,payCal);

					// resetcal of structured legs
					GetResetCalendar(listItem,resetCal);

					intrule = GetIntRule(listItem);
					
					// Echange de Notional
					NxId = GetNotionalExchangeFlag(listItem);

					// Recuperation des evenements pour les fixings FX
					listItem->selectSingleNode(_bstr_t((const char *)"Events"), &theNode);
					if ( theNode != NULL )
					{
						if (theNode->selectNodes(_bstr_t((const char *)("EVENT")), &resultList2) == S_OK)
						{
							resultList2->get_length(&nbEvents);

							CCString type;

							for (int i = 0; i < nbEvents; i++)
							{
								hr=resultList2->get_item(i, &listItem2);
								
								listItem2->selectSingleNode(_bstr_t((const char *)"Type"), &theNode2);
								if (theNode2!=NULL)
								{
									BSTR resultat = NULL;
									theNode2->get_text(&resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									type = (const char*) ff1;

									theNode2->Release();
									theNode2=NULL;
									if (resultat) SysFreeString(resultat);
								}

								if ( strcmp(type,"FX") == 0 )
								{
									listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
									if (theNode2!=NULL)
									{
										BSTR resultat = NULL;
										theNode2->get_text(&resultat);

										_bstr_t ff(resultat,false);
										char * ff1=(char *)ff;

										double dfxfixing = atof(ff1);
										fxFixing.push_back(dfxfixing);

										theNode2->Release();
										theNode2=NULL;
										if (resultat) SysFreeString(resultat);
									}

									resetFixing.push_back(GetDateFromXMLNode(listItem2,"ADate").GetJulian());
								
								}
								else if ( strcmp(type,"SPR") == 0 )
								{
									listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
									if ( theNode2 != NULL )
									{
										BSTR resultat = NULL;
										theNode2->get_text(&resultat);

										_bstr_t ff(resultat,false);
										char * ff1=(char *)ff;

										*dualStrike = atof(ff1) * 10000.;

										theNode2->Release();
										theNode2=NULL;
										if (resultat) SysFreeString(resultat);
									}
								}

								listItem2->Release();
								listItem2 = NULL;
							}

							resultList2->Release();
							resultList2 = NULL;
						}
					}

					if (listItem->selectNodes(_bstr_t((const char *)("ProdData/cFXSTRP/Optlets/cFXSTRC")), &resultListRefValue) == S_OK)
					{
						resultListRefValue->get_length(&nbNodes2);

						if (nbNodes2 == 0)
						{
							hr = S_FALSE;

							CCString msg((CCString)"Invalid XML string for getting PRCS (RefValue Cap, Floor et FX");

							throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
											(char*) msg);
						}

						vTmpFx = new ARM_Vector(nbNodes2);
						vTmpCap = new ARM_Vector(nbNodes2);
						vTmpFloor = new ARM_Vector(nbNodes2);
						vTmpCpnRate1 = new ARM_Vector(nbNodes2);
						vTmpCpnRate2 = new ARM_Vector(nbNodes2);

						nbStrip = nbNodes2;
						long indexStrip = 0;
	
						for (long indexNode=0 ; indexNode<nbNodes2 ; indexNode++)
						{
							hr = resultListRefValue->get_item(indexNode, &listItem2);

							if ( ( indexNode == 0 ) && ( typeDeal == 1 ) )
							   startDate = GetDateFromXMLNode(listItem2,"FromDate");

							// FX
							listItem2->selectSingleNode(_bstr_t((const char *)"FXSpotRate"), &theNode);
							if ( theNode != NULL )
							{
								BSTR resultat = NULL;
								theNode->get_text(&resultat);

								_bstr_t ff(resultat,false);
								char * ff1=(char *)ff;

								vTmpFx->Elt(indexStrip) = atof((const char*) ff1);

								theNode->Release();
								theNode=NULL;
								if (resultat) 
                                   SysFreeString(resultat);

								listItem2->selectSingleNode(_bstr_t((const char *)"CoupRate1"), &theNode);
								if (theNode!=NULL)
								{
									resultat = NULL;
									theNode->get_text(&resultat);

									_bstr_t ff(resultat,false);
									ff1=(char *)ff;

									vTmpCpnRate1->Elt(indexStrip) = atof((const char*)ff1) * 100.;

									theNode->Release();
									theNode=NULL;
									if (resultat) 
                                       SysFreeString(resultat);
								}

								listItem2->selectSingleNode(_bstr_t((const char *)"CoupRate2"), &theNode);
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									vTmpCpnRate2->Elt(indexStrip) = atof((const char*)ff1) * 100.;

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}

								// Cap
								listItem2->selectSingleNode(_bstr_t((const char *)"CapRate"), &theNode);
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									double tmpCap = atof((const char*) ff1);
									if (tmpCap == 0.0)
										vTmpCap->Elt(indexStrip) = 10000.;
									else
										vTmpCap->Elt(indexStrip) = tmpCap * 100.;

									theNode->Release();
									theNode=NULL;
									if (resultat) SysFreeString(resultat);
								}
								else
									vTmpCap->Elt(indexStrip) = 10000.;

								// Floor
								listItem2->selectSingleNode(_bstr_t((const char *)"FloorRate"), &theNode);
								if (theNode!=NULL)
								{
									BSTR resultat = NULL;
									theNode->get_text(&resultat);

									_bstr_t ff(resultat,false);
									char * ff1=(char *)ff;

									vTmpFloor->Elt(indexStrip) = atof((const char*) ff1) * 100.;

									theNode->Release();
									theNode=NULL;

									if (resultat) 
                                       SysFreeString(resultat);
								}
								else
									vTmpFloor->Elt(indexStrip) = 0.;

								indexStrip ++;
							}
							else
							{
								nbStrip --;
							}
						}
					}
				}
			}
		}
	}

	resultList->Release();
	resultList = NULL;

	if ( typeDeal == 0 )
	   tmpChaine = "Response/EXOTIC/Assets/ASSET/ProdData/cFXSTRP/FXOptionData/FXOPTION";
	else
	   tmpChaine = "Response/SWAPTION/Assets/ASSET/ProdData/cFXSTRP/FXOptionData/FXOPTION";

	if ( XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK )
	{
		resultList->get_length(&nbNodes);

		if ( nbNodes != 1 )
		{
			hr = S_FALSE;

			CCString msg((CCString)"No FXOPTION tag in FX_OPT_STRIP asset");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							(const char *) msg);
		}

		hr = resultList->get_item(0, &listItem);

		// Main Currency of FX_OPT_STRIP asset
		listItem->selectSingleNode(_bstr_t((const char *)"MainCcy"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			ccy1 = new ARM_Currency(ff1);

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// Money Currency of FX_OPT_STRIP asser
		listItem->selectSingleNode(_bstr_t((const char *)"MoneyCcy"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			ccy2 = new ARM_Currency(ff1);

			theNode->Release();
			theNode=NULL;
	
            if (resultat) 
               SysFreeString(resultat);
		}
	}

	resultList->Release();
	resultList = NULL;

	ARM_Vector* vCap = new ARM_Vector(vTmpCap,0,nbStrip-1);
	ARM_Vector* vFloor = new ARM_Vector(vTmpFloor,0,nbStrip-1);
	ARM_Vector* vFx = new ARM_Vector(vTmpFx,0,nbStrip-1);
	ARM_Vector* vCpnRate1 = new ARM_Vector(vTmpCpnRate1, 0, nbStrip-1);
	ARM_Vector* vCpnRate2 = new ARM_Vector(vTmpCpnRate2, 0, nbStrip-1);

	delete vTmpCap;
	delete vTmpFloor;
	delete vTmpFx;
	delete vTmpCpnRate1;
	delete vTmpCpnRate2;

	ARM_Date endDateNa(endDate);
	endDateNa.SetDay(payAnnDay);

	FxUnderIndex = new ARM_IRIndex(ccy1->GetCcyName(),daycountBasis);

	FxUnderIndex->SetResetFrequency(fixFrequency);
	FxUnderIndex->SetPayFrequency(fixFrequency);

	FxUnderIndex->SetResetTiming(resetTiming);
	FxUnderIndex->SetPayTiming(payTiming);

	FxUnderIndex->SetCompMeth(K_COMP_PROP);
	FxUnderIndex->SetFwdRule(K_MOD_FOLLOWING);

	FxUnderIndex->SetResetGap(resetGap);
	FxUnderIndex->SetPayGap(payGap);

	FxUnderIndex->SetYearTerm(1./(double)FxUnderIndex->GetResetFrequency());
	FxUnderIndex->SetTerm(FxUnderIndex->GetResetFrequency());

	FxUnderIndex->SetIntRule(intrule);

	ARM_SwapLeg FxUnderTmpLeg(startDate,
							  endDateNa,
							  FxUnderIndex,						 
							  K_PAY,
							  0.0,
							  stubType,
							  K_COMP_PROP,
							  ccy1,
							  daycountBasis,
							  resetGap,
							  resetCal,
							  payCal,
							  1,
							  K_NX_NONE,
							  stubDate);

	ARM_Vector* vAmountFxUnder = new ARM_Vector(vFx->GetSize());
	for (int i = 0; i < vAmountFxUnder->GetSize(); i++)
	{
		if (vFx->Elt(i) != 0.0)
			vAmountFxUnder->Elt(i) = notional/vFx->Elt(i);
		else
			vAmountFxUnder->Elt(i) = 0.0;
	}

	ARM_Vector* resetUndDates = FxUnderTmpLeg.GetResetDates();
	ARM_ReferenceValue* rCpnRate2 = new ARM_ReferenceValue((ARM_Vector *) resetUndDates->Clone(),
                                                           (ARM_Vector *) vCpnRate2->Clone());
	delete vCpnRate2;

	ARM_FixLeg* fixUndLeg = new ARM_FixLeg();

	fixUndLeg->SetVarCoupons(rCpnRate2);
    if (rCpnRate2)
       delete rCpnRate2;

	*FxUnderLeg = (ARM_SwapLeg *) fixUndLeg;

	**FxUnderLeg = FxUnderTmpLeg;

	rFxUnder = new ARM_ReferenceValue((ARM_Vector *)((*FxUnderLeg)->GetPaymentDates())->Clone(),
                                      vAmountFxUnder);

	(*FxUnderLeg)->SetAmount(rFxUnder);
	delete rFxUnder;

	FxNumIndex = new ARM_IRIndex(ccy2->GetCcyName(),daycountBasis);

	FxNumIndex->SetResetFrequency(fixFrequency);
	FxNumIndex->SetPayFrequency(fixFrequency);

	FxNumIndex->SetCompMeth(K_COMP_PROP);
	FxNumIndex->SetFwdRule(K_MOD_FOLLOWING);

	FxNumIndex->SetYearTerm(1./(double)FxUnderIndex->GetResetFrequency());
	FxNumIndex->SetTerm(FxUnderIndex->GetResetFrequency());

	FxNumIndex->SetIntRule(intrule);

	ARM_SwapLeg FxNumTmpLeg(startDate,
							endDateNa,
							FxNumIndex,						 
							K_RCV,
							0.0,
							stubType,
							K_COMP_PROP,
							ccy2,
							daycountBasis,
							10000,
							resetCal,
							payCal,
							1,
							K_NX_END,
							stubDate);

    if (FxNumIndex)
    {
        delete FxNumIndex;

        FxNumIndex = NULL;
    }

    if (FxUnderIndex)
    {
       delete FxUnderIndex;

       FxUnderIndex = NULL;
    }

	if (ccy2)
       delete ccy2;

	ARM_Vector* resetNumDates = FxNumTmpLeg.GetResetDates();
	ARM_ReferenceValue* rCpnRate1 = new ARM_ReferenceValue((ARM_Vector *) resetNumDates->Clone(),
                                                           (ARM_Vector *) vCpnRate1->Clone());

	delete vCpnRate1;

	ARM_FixLeg* fixNumLeg = new ARM_FixLeg();

	fixNumLeg->SetVarCoupons(rCpnRate1);
    if (rCpnRate1)
       delete rCpnRate1;

	FxNumLeg = (ARM_SwapLeg *) fixNumLeg;

	*FxNumLeg = FxNumTmpLeg;

	amount = new ARM_ReferenceValue(notional);
	FxNumLeg->SetAmount(amount);
	delete amount;

	*rFx = new ARM_ReferenceValue((ARM_Vector *) resetUndDates->Clone(),(ARM_Vector*)vFx->Clone());
	*rCap = new ARM_ReferenceValue((ARM_Vector *) resetUndDates->Clone(),(ARM_Vector*)vCap->Clone());
	*rFloor = new ARM_ReferenceValue((ARM_Vector *) resetUndDates->Clone(),(ARM_Vector*)vFloor->Clone());

	delete vFx;
	delete vCap;
	delete vFloor;

	// Récupération des fixings FX
	*FxFixing = NULL;
	if ( resetFixing.size() > 0 )
	{
		ARM_Vector* vResetFixing = new ARM_Vector(resetFixing.size());
		ARM_Vector* vFxFixing = new ARM_Vector(resetFixing.size());

		for (int i = 0; i < resetFixing.size(); i++)
		{
			vResetFixing->Elt(i) = resetFixing[i];
			vFxFixing->Elt(i) = fxFixing[i];
		}

		*FxFixing = new ARM_ReferenceValue(vResetFixing,vFxFixing);
	}

    if (ccy1)
       delete ccy1;

	return(FxNumLeg);
}



ARM_SwapLeg* getFundingLeg(MSXML2::IXMLDOMDocument* XMLDoc,
						   ARM_ReferenceValue** LiborFixing,
						   int typeDeal) // si 0 : structure, si 1 : Monotrade
{
	ARM_SwapLeg* realFundLeg = NULL;

	MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * listItem2 = NULL, * theNode = NULL, * theNode2 = NULL;
	HRESULT hr;
	long nbNodes;
	long nbEvents;

	ARM_Date startDate;
	ARM_Date endDate;
	ARM_Currency ccy;

	int daycountBasis;
	int resetFreq;
	int payFreq;
	int resetGap;
	int payGap;
	int resetTiming;
	int payTiming;
	int refDay;
	int compMode;

	CCString Index;
	CCString Term;
	CCString Source;

	long NxId = K_NX_BOTH;

	char myRefDate[11];
	char payCal[4];
	char resetCal[4];

	ARM_Date startStubDate;
	ARM_Date endStubDate;

	VECTOR<double> spreadDate;
	VECTOR<double> spreadVal;

	VECTOR<double> liborFixing;
	VECTOR<double> resetFixing;
	VECTOR<double> startDateBefFix;

	double dNotional;

	VECTOR<double> FEEDate;
	VECTOR<double> FEEVal;

	VECTOR<double> DateOFF;
	VECTOR<double> DateROF;
	VECTOR<double> ADateROF;

	int index;

	CCString tmpChaine;

	if ( typeDeal == 0 )
	   tmpChaine = "Response/EXOTIC/Assets/ASSET";
	else
	   tmpChaine = "Response/SWAPTION/Assets/ASSET";

	// Recuperation de la funding leg
	if ( XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (!(( nbNodes == 1 ) || ( nbNodes == 2 )))
		{
			hr = S_FALSE;

			BSTR tmpChaine;

			hr = XMLDoc->get_xml(&tmpChaine);

			CCString msg((CCString)"Invalid String for getting Funding (just 1 asset)");

			if (tmpChaine)
				SysFreeString(tmpChaine);

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							(const char *) msg);
		}

		if ((typeDeal == 0) && (nbNodes == 2))
		{
			vector<ARM_Date> startDates, endDates;

			for (int i = 0; i < nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);

				// start of fundingLeg
				startDates.push_back(GetStartDate(listItem));

				// Maturity Date
				endDates.push_back(GetEndDate(listItem));

				listItem->Release();
				listItem = NULL;
			}

			if ( startDates[0] < startDates[1] )
			{
				startDate = startDates[0];
				endDate = endDates[1];

				index = 1;
			}
			else
			{
				startDate = startDates[1];
				endDate = endDates[0];
				
                index = 0;
			}
		}
		else if (( typeDeal == 0 ) && ( nbNodes == 1 ))
		{
			hr = resultList->get_item(0, &listItem);

			// start of fundingLeg
			startDate = GetStartDate(listItem);

			// Maturity Date
			endDate = GetEndDate(listItem);

			listItem->Release();
			listItem=NULL;

			index = 0;
		}
		// type MonoTrade, on cherche la jambe de funding (celle qui n'a pas un productName like FX_OPT_STRIP)
		else
		{
			for (int i = 0; i < nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);
				listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);
			
                if ( theNode != NULL )
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					theNode->Release();
					theNode = NULL;

					if ( strncmp((const char *) ff1, "FX_OPT_STRIP", 12) != 0 )
					{
						if (resultat) 
                           SysFreeString(resultat);

						// start of legs
						startDate = GetStartDate(listItem);

						// end of legs
						endDate = GetEndDate(listItem);

						index = i;
					}

					listItem->Release();
					listItem=NULL;
				}
			}
		}

		hr = resultList->get_item(index, &listItem);

		// ccy1
		listItem->selectSingleNode(_bstr_t((const char *)"Ccy"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			ARM_Currency tmpccy(ff1);
			ccy = tmpccy;

			theNode->Release();
			theNode=NULL;
		
            if (resultat) 
               SysFreeString(resultat);
		}

		// basis
		daycountBasis = GetAssetDayCount(listItem);

		// Notional
		listItem->selectSingleNode(_bstr_t((const char *)"Notional"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			dNotional = atof(ff1);

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// Echange de Notional
		NxId = GetNotionalExchangeFlag(listItem);
		
		// Index
		listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_dmIndex"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			Index = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// Term
		listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_Term"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			Term = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// Spread
		listItem->selectSingleNode(_bstr_t((const char *)"INTEREST_Spread"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			double dspread = atof(ff1) * 100.;
			spreadVal.push_back(dspread);
			spreadDate.push_back(startDate.GetJulian());

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// reset Frequency
		resetFreq = GetResetFreq(listItem);

		// reset Gap
		resetGap = GetResetGap(listItem);

		// reset Timing
		resetTiming = GetResetTiming(listItem);

		// pay Frequency
		payFreq = GetPayFreq(listItem);

		// pay Gap
		payGap = GetPayGap(listItem);

		// pay Timing
		payTiming = GetPayTiming(listItem);

		// compounding method
		compMode = GetCompMode(listItem);

		// reference Date
		listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date1"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			if (strcmp((const char*)ff1,"") == 0)
				strcpy(myRefDate,"NULL");
			else
			{
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				startStubDate = tmpDate;
			}

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// reference Date
		listItem->selectSingleNode(_bstr_t((const char *)"STUB_Date2"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			if (strcmp((const char*)ff1,"") != 0)
			{
				ARM_Date tmpDate (ff1,"YYYYMMDD");
				endStubDate = tmpDate;
			}

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// reference Day1
		listItem->selectSingleNode(_bstr_t((const char *)"SCHED_Pay_AnnDay"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			refDay = atoi((const char*) ff1);

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		// paycal
		GetPayCalendar(listItem,payCal);

		// resetcal
		GetResetCalendar(listItem,resetCal);

		int nbXNL = 0;

		// Recuperation des evenements
		listItem->selectSingleNode(_bstr_t((const char *)"Events"), &theNode);
		if (theNode!=NULL)
		{
			if (theNode->selectNodes(_bstr_t((const char *)("EVENT")), &resultList2) == S_OK)
			{
				resultList2->get_length(&nbEvents);

				CCString type;

				for (int i = 0; i < nbEvents; i++)
				{
					hr=resultList2->get_item(i, &listItem2);
					
					listItem2->selectSingleNode(_bstr_t((const char *)"Type"), &theNode2);
					if (theNode2!=NULL)
					{
						BSTR resultat = NULL;
						theNode2->get_text(&resultat);

						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						type = (const char*) ff1;

						theNode2->Release();
						theNode2=NULL;
						if (resultat) SysFreeString(resultat);
					}

					// Spread variable
					if (strcmp(type,"SPR") == 0)
					{
						listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
						if (theNode2!=NULL)
						{
							BSTR resultat = NULL;
							theNode2->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							double dspread = atof(ff1) * 100.;
							spreadVal.push_back(dspread);

							theNode2->Release();
							theNode2=NULL;
							
                            if (resultat) 
                               SysFreeString(resultat);
						}

						spreadDate.push_back(GetDateFromXMLNode(listItem2,"Date").GetJulian());
					}

					// Fixing
					else if (strcmp(type,"FRC") == 0)
					{
						listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
						if (theNode2!=NULL)
						{
							BSTR resultat = NULL;
							theNode2->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							double dfixing = atof(ff1)*100.;
							liborFixing.push_back(dfixing);

							theNode2->Release();
							theNode2 = NULL;
					
                            if (resultat) 
                               SysFreeString(resultat);
						}

						startDateBefFix.push_back(GetDateFromXMLNode(listItem2,"Date").GetJulian());

						resetFixing.push_back(GetDateFromXMLNode(listItem2,"ADate").GetJulian());
					}
					// Nominal variable
					else if (strcmp(type,"XNL") == 0)
					{
						listItem2->selectSingleNode(_bstr_t((const char *)"Amount"), &theNode2);
						if ( theNode2 != NULL )
						{
							BSTR resultat = NULL;
							theNode2->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1 = (char *)ff;

							// dans le cas MonoTrade, la jambe est stockée en payeuse
							// alors qu'elle est receveuse : on inverse le nominal
							if ( typeDeal == 1 )
								dNotional = -atof(ff1);
							else
								dNotional = atof(ff1);

							theNode2->Release();
							theNode2=NULL;
							if (resultat) SysFreeString(resultat);
						}
					}

					// Echeancier customizé
					else if ( strcmp(type, "OFF") == 0 )
					{
						DateOFF.push_back(GetDateFromXMLNode(listItem2,"Date").GetJulian());
					}
					else if ( strcmp(type, "ROF") == 0 )
					{
						DateROF.push_back(GetDateFromXMLNode(listItem2,"Date").GetJulian());
						
						ADateROF.push_back(GetDateFromXMLNode(listItem2,"ADate").GetJulian());
					}

					listItem2->Release();
					listItem2 = NULL;
				}

				resultList2->Release();
				resultList2 = NULL;
			}
		}


		listItem->Release();
		listItem=NULL;
	}

	resultList->Release();
	resultList=NULL;

	ARM_Date maturityDateNA (endDate);

	ARM_Date dateToday;

	if ( startStubDate != dateToday )
	   startStubDate.JulianToStrDate(myRefDate);
	else if ( endStubDate != dateToday )
	{
		ARM_Date tmpDate (startDate);
		tmpDate.AddPeriodMult(payFreq,1,ccy.GetCcyName());
		tmpDate.SetDay(endStubDate.GetDay());
		tmpDate.JulianToStrDate(myRefDate);
	}
	else 
		strcpy(myRefDate,"NULL");

	if ( endStubDate == dateToday )
	   maturityDateNA.SetDay(refDay);

	int indexType = FromStringToIndexType(Term,Index);

	if ( DateROF.size() == 0 )
	{
		long stubRuleId = K_SHORTSTART;

		if (strcmp((const char*) myRefDate, "NULL") != 0)
		{
			ARM_Date tmpDate (endStubDate);
			tmpDate.AddPeriodMult(payFreq, 1, ccy.GetCcyName());
			
            if ( tmpDate > maturityDateNA )
			   stubRuleId = K_SHORTEND;
			else
			   stubRuleId = K_LONGEND;

			ARM_Date refDate(myRefDate);

			if ( refDate.GetDay() != refDay )
			{
				refDate.SetDay(refDay);
				refDate.JulianToStrDate(myRefDate);
			}
		}

		realFundLeg = new ARM_SwapLeg(startDate,
									  maturityDateNA,
									  (ARM_INDEX_TYPE)indexType,
									  K_RCV,
									  0.0,
									  resetFreq,
									  payFreq,
									  resetTiming,
									  payTiming,
									  &ccy,
									  K_ADJUSTED,
									  resetGap,
									  resetCal,
									  payCal,
									  1,
									  NxId,
									  stubRuleId,
									  myRefDate,
									  0);

		ARM_ReferenceValue* rNotional = NULL;

		rNotional = new ARM_ReferenceValue(dNotional);

		realFundLeg->SetAmount(rNotional);
        delete rNotional;

	}
	else
    { 
		ARM_Vector* vStartDates = new ARM_Vector(DateROF.size()+DateOFF.size()-1);
		ARM_Vector* vEndDates = new ARM_Vector(DateROF.size()+DateOFF.size()-1);
		ARM_Vector* vPaymentDates = new ARM_Vector(DateROF.size()+DateOFF.size()-1);
		ARM_Vector* vResetDates = new ARM_Vector(DateROF.size()+DateOFF.size()-1);
		ARM_Vector* vIntDays = new ARM_Vector(DateROF.size()+DateOFF.size()-1);
		ARM_Vector* vNotional = new ARM_Vector(DateROF.size()+DateOFF.size()-1);

		ARM_Date tmpDate;

		for (int i = 0; i < startDateBefFix.size(); i++)
		{
			vStartDates->Elt(i) = startDateBefFix[i];
			vResetDates->Elt(i) = resetFixing[i];
			
			if ( i != startDateBefFix.size()-1 )
			{
				vEndDates->Elt(i) = startDateBefFix[i+1];
				tmpDate = startDateBefFix[i+1];
				vPaymentDates->Elt(i) = tmpDate.PreviousBusinessDay(-payGap,payCal).GetJulian();

		        vIntDays->Elt(i) = DaysBetweenDates(daycountBasis, vStartDates->Elt(i), vEndDates->Elt(i));
			}
		}

		vEndDates->Elt(startDateBefFix.size()-1) = DateROF[0];
		tmpDate = vEndDates->Elt(startDateBefFix.size()-1);
		vPaymentDates->Elt(startDateBefFix.size()-1) = tmpDate.PreviousBusinessDay(-payGap,payCal).GetJulian();

        vIntDays->Elt(startDateBefFix.size()-1) = DaysBetweenDates(daycountBasis, vStartDates->Elt(startDateBefFix.size()-1), vEndDates->Elt(startDateBefFix.size()-1));

		for (i = 0; i < DateROF.size(); i++)
		{
			vStartDates->Elt(startDateBefFix.size()+i) = DateROF[i];
			vResetDates->Elt(startDateBefFix.size()+i) = ADateROF[i];
			
			if ( i != DateROF.size()-1 )
			{
				vEndDates->Elt(startDateBefFix.size()+i) = DateROF[i+1];
				tmpDate = DateROF[i+1];
				vPaymentDates->Elt(startDateBefFix.size()+i) = tmpDate.PreviousBusinessDay(-payGap,payCal).GetJulian();

		        vIntDays->Elt(startDateBefFix.size()+i) = DaysBetweenDates(daycountBasis, vStartDates->Elt(startDateBefFix.size()+i), vEndDates->Elt(startDateBefFix.size()+i));
			}
		}

		vEndDates->Elt(DateROF.size()+startDateBefFix.size()-1) = DateOFF[DateOFF.size()-1];
		tmpDate = vEndDates->Elt(DateROF.size()+startDateBefFix.size()-1);
		vPaymentDates->Elt(DateROF.size()+startDateBefFix.size()-1) = tmpDate.PreviousBusinessDay(-payGap,payCal).GetJulian();

        vIntDays->Elt(DateROF.size()+startDateBefFix.size()-1) = DaysBetweenDates(daycountBasis, vStartDates->Elt(DateROF.size()+startDateBefFix.size()-1), vEndDates->Elt(DateROF.size()+startDateBefFix.size()-1));

		ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)indexType, resetFreq,payFreq, &ccy);

		ARM_ReferenceValue* amount = new ARM_ReferenceValue(dNotional);

		realFundLeg = new ARM_SwapLeg(vStartDates, vEndDates,
									  vPaymentDates, vResetDates,
									  vIntDays, NULL,
									  amount, newIndex,
									  K_RCV, 0.0,
									  K_FLOAT_RATE, 
									  &ccy);

		realFundLeg->SetResetCalName(resetCal);
		realFundLeg->SetPayCalName(payCal);

		realFundLeg->SetStartDateNA(startDate);
		realFundLeg->SetEndDateNA(endDate);

		realFundLeg->SetNxFlag(NxId);

		realFundLeg->SetAmount(amount);

		if (vStartDates)
			delete vStartDates;

		if (vEndDates)
			delete vEndDates;

		if (vPaymentDates)
			delete vPaymentDates;

		if (vResetDates)
			delete vResetDates;

		if (vIntDays)
			delete vIntDays;

        if (amount)
           delete amount;

        if (newIndex)
		   delete newIndex;
	}

	realFundLeg->SetCompoundingType(compMode);

	//construction du spread
	ARM_ReferenceValue* rSpread = NULL;
	if ( spreadDate.size() == 0 )
	   rSpread = new ARM_ReferenceValue(0.0);
	else if ( spreadDate.size() == 1 )
	   rSpread = new ARM_ReferenceValue(spreadVal[0]);
	else
	{
		ARM_Vector* dateSpread = new ARM_Vector(spreadDate.size());
		ARM_Vector* valSpread  = new ARM_Vector(spreadDate.size());

		if ( spreadDate.size() > 0 )
		{
			for (int i = 0; i < spreadDate.size(); i++)
			{
				valSpread->Elt(i) = spreadVal[i];
				ARM_Vector* tmpResetDates = realFundLeg->GetResetDates();
				ARM_Vector* tmpStartDates = realFundLeg->GetFlowStartDates();
				
                int j = 0;
				while (tmpStartDates->Elt(j) < spreadDate[i])
					j++;

				dateSpread->Elt(i) = tmpResetDates->Elt(j);
			}
		}

		rSpread = new ARM_ReferenceValue(dateSpread,valSpread);

		rSpread->SetCalcMethod(K_STEPUP_LEFT);
	}

	realFundLeg->SetVariableSpread(rSpread);
	delete rSpread;

	// Construction des fixings
	*LiborFixing = NULL;
	if (resetFixing.size() > 0)
	{
		ARM_Vector* vResetFixing = new ARM_Vector(resetFixing.size());
		ARM_Vector* vLiborFixing = new ARM_Vector(resetFixing.size());

		for (int i = 0; i < resetFixing.size(); i++)
		{
			vResetFixing->Elt(i) = resetFixing[i];
			vLiborFixing->Elt(i) = liborFixing[i];
		}

		*LiborFixing = new ARM_ReferenceValue(vResetFixing, vLiborFixing);
	}

	return(realFundLeg);
}



ARM_PowerReverse* ARMLOCAL_ParseGenPRCS(const char* chaineXML,
										const ARM_Date& date,
										CCString& bookName,
										CCString& structureId,
										CCString& custId,
										CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode  * listItem = NULL, * theNode = NULL;
	long nbNodes;

	CCString productName = "";

	// Récupération de la structureId
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting PRCS");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	// Recuperation du book et de la structure
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Option/OPTION")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (nbNodes != 1)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		hr = resultList->get_item(0, &listItem);

		listItem->selectSingleNode(_bstr_t((const char *)"ProductName"), &theNode);
		if ( theNode != NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			productName = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		listItem->Release();
		listItem=NULL;
	}

	resultList->Release();
	resultList = NULL;

	XMLDoc->Release();
	XMLDoc = NULL;

	if (strncmp((const char*)productName,"OPT_PRDC",8) == 0)
	{
		structureId = "";
		
		return ARMLOCAL_ParseMONOPRCS(chaineXML,
									  date,
									  bookName,
									  structureId,
									  custId,
									  dealId);
	}
	else
	{
		return ARMLOCAL_ParsePRCS(chaineXML,
								  date,
								  bookName,
								  structureId,
								  custId,
								  dealId);
	}
}



ARM_PowerReverse* ARMLOCAL_ParsePRCS(const char* chaineXML,
									 const ARM_Date& date,
									 CCString& bookName,
									 CCString& structureId,
									 CCString& custId,
									 CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	ARM_ReferenceValue* LiborFixing = NULL;
	ARM_ReferenceValue* FxFixing = NULL;

	ARM_SwapLeg* initPeriodLeg = NULL;
	ARM_SwapLeg* realFundLeg = NULL;
	ARM_SwapLeg* fxUnderLeg = NULL;
	ARM_SwapLeg* fxNumLeg = NULL;

	int isFirstPeriodLeg = 0;

	CCString desk;

    ARM_Date swapStartDate;
	ARM_Date fixedSwapEndDate;

	ARM_Date firstNoticeDate;
	ARM_Date noticeDate;
	ARM_Vector* vNoticeDates = NULL;
	ARM_Vector* vCancelDates = NULL;
	ARM_Vector* vPaymentDates = NULL;
	ARM_Vector* vFx = NULL;
	ARM_Vector* vCap = NULL;
	ARM_Vector* vFloor = NULL;

	ARM_ReferenceValue* rCap = NULL;
	ARM_ReferenceValue* rFloor = NULL;
	ARM_ReferenceValue* rFx = NULL;

	ARM_PowerReverse* newPowerReverse = NULL;

	CCString xmlResponse;
	CCString messageList;

	vector<string> listXMLResponse;


	// Récupération de la structureId
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting PRCS");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	long nbNodes;

	// Recuperation du book et de la structure
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Env/ENV")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (nbNodes != 1)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		hr = resultList->get_item(0, &listItem);

		listItem->selectSingleNode(_bstr_t((const char *)"StructureId"), &theNode);

        if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			structureId = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		listItem->selectSingleNode(_bstr_t((const char *)"Desk"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			desk = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}


		listItem->Release();
		listItem=NULL;
	}

	XMLDoc->Release();
	XMLDoc = NULL;

	resultList->Release();
	resultList=NULL;

	listXMLResponse = etoolkit_getTradeListByStructId(structureId,desk);

	// Recuperation de l'Exotic qui contient le FX_OPT_STRIP
	MSXML2::IXMLDOMDocument* XMLDocForFxOptStrip = NULL;
	MSXML2::IXMLDOMDocument* XMLDocForFunding = NULL;

	XMLDocForFxOptStrip = searchDocForExoticFxOptStrip(listXMLResponse);
	XMLDocForFunding = searchDocForExoticFunding(listXMLResponse);

	/*--------------------------------------*/

	// Recuperation de la swaption DUMMY
	MSXML2::IXMLDOMDocument* XMLDocSwaption = NULL;

	XMLDocSwaption = searchDocForSwaption(listXMLResponse);

	/**********************************************************/
	
	vNoticeDates = getNotAndCancDates(XMLDocSwaption,&vCancelDates,custId,dealId);
	/**********************************************************/

	initPeriodLeg = getInitPeriodLeg(XMLDocForFxOptStrip,0);
	/**********************************************************/

	int dualFlag = 0;
	double dualStrike = 30.0;

	fxNumLeg = getStructuredLegs(XMLDocForFxOptStrip,
								 &fxUnderLeg,
								 &rCap,
								 &rFloor,
								 &rFx,
								 &FxFixing,
								 &dualFlag,
								 &dualStrike,
								 bookName,
								 0);

	noticeDate = fxUnderLeg->GetResetDates()->Elt(fxUnderLeg->GetResetDates()->GetSize() - 1);

	realFundLeg = getFundingLeg(XMLDocForFunding,&LiborFixing,0);

	newPowerReverse = new ARM_PowerReverse(initPeriodLeg,
										   realFundLeg,
										   fxUnderLeg,
										   fxNumLeg,
										   vNoticeDates,
										   vCancelDates,
										   rFx,
										   rCap,
										   rFloor,
										   dualFlag,
										   dualStrike,
										   noticeDate);

    double  dNotional = fxNumLeg->GetAmount()->CptReferenceValue((ARM_Date) fxNumLeg->GetPaymentDates()->Elt(0));
	
    ARM_ReferenceValue* rNotional = new ARM_ReferenceValue(dNotional);

	newPowerReverse->SetAmount(rNotional);
	delete rNotional;

	XMLDocSwaption->Release();
	XMLDocSwaption = NULL;

	XMLDocForFxOptStrip->Release();
	XMLDocForFxOptStrip = NULL;

	XMLDocForFunding->Release();
	XMLDocForFunding = NULL;

	delete vNoticeDates;
	delete vCancelDates;

	delete fxNumLeg;
	delete fxUnderLeg;
	delete initPeriodLeg;
	delete realFundLeg;

    if (rFx)
       delete rFx;
    
    if (rCap)
       delete rCap;

    if (rFloor)
       delete rFloor;
    


	// Gestion des fixings
	// Valeur par défaut MinDate : Pas de gestion de fixing
	if ( date.isMinDate() == 0 )
	{
		double today = date.GetJulian();

		// d'abord les fixings LIBOR
		if ( LiborFixing != NULL ) /// otherwise no fixing in summit
		{
			/// reset dates for the funding leg
			ARM_Vector* fundingLegResetDates = newPowerReverse->GetItsRealFundLeg()->GetResetDates();
			/// reset dates from summit until summit current date
			ARM_Vector* summitResetDates = LiborFixing->GetDiscreteDates();
			/// computation date
			if ( fundingLegResetDates->Elt(0) <= today && summitResetDates->Elt(0) <= today ) 
			{
				int size = 0;
				for (int i = 0; i < MIN(fundingLegResetDates->size(), summitResetDates->size()); i++)
				{
					if ((fundingLegResetDates->Elt(i) == summitResetDates->Elt(i)) && 
						(fundingLegResetDates->Elt(i) <= today)) 
					{
						size++;
					}
					else 
                    {
						if (fundingLegResetDates->Elt(i) > today) 
						{
							break;
						}
						else 
						{
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								"Mismatch reset date"); 
						}
					}
				}

				ARM_Vector* summitResetValues = LiborFixing->GetDiscreteValues();
				ARM_Vector fundingLegResetValues(size);

				for (i = 0; i < size; i++) 
				{
					fundingLegResetValues.InitElt(i, summitResetValues->Elt(i));
				}

				newPowerReverse->GetItsRealFundLeg()->SetFixRates(&fundingLegResetValues);
			}
		}

		// Puis les fixings FX
		int rg = 0;

		ARM_Vector* fxResetDates = newPowerReverse->GetItsFxUnderLeg()->GetResetDates();

        if ( fxResetDates->Elt(0) <= today ) 
		{			
			for (int i = 1; i < fxResetDates->size(); i++)
			{
				if (( fxResetDates->Elt(i-1) <= today ) && ( fxResetDates->Elt(i) > today )) 
				{
					rg = i-1; /// find the reset date for the fx side
					break; //i = fxResetDates->size();
				}
			}			
			
			ARM_Date FxResetDate(fxResetDates->Elt(rg));			
			bool isFound(false);

			if (FxFixing)
			{
				for ( i = 0; i < FxFixing->GetSize(); i++)
				{
					if ( ARM_Date(FxFixing->GetDiscreteDates()->Elt(i)) == FxResetDate ) 
					{
						rg = i; /// find the associated reset date for the fx side in the table 
						isFound = true;
						break;
					}
				}
			}

			if ((!isFound) && ( FxResetDate != date )) 
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                                "No reset date associated to the FX reset date");
			}

			double fixing(0.0);

			if (isFound) 
			{
				fixing = FxFixing->GetDiscreteValues()->Elt(rg);
			}

			newPowerReverse->SetLastFxFixing(fixing);
		}
	}

	if (LiborFixing)
		delete LiborFixing;
	LiborFixing = NULL;

	if (FxFixing)
		delete FxFixing;
	FxFixing = NULL;

	return(newPowerReverse);
}



ARM_PowerReverse* ARMLOCAL_ParseMONOPRCS(const char* chaineXML,
										 const ARM_Date& date,
										 CCString& bookName,
										 CCString& structureId,
										 CCString& custId,
										 CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_ReferenceValue* LiborFixing = NULL;
	ARM_ReferenceValue* FxFixing = NULL;

	ARM_SwapLeg* initPeriodLeg = NULL;
	ARM_SwapLeg* realFundLeg = NULL;
	ARM_SwapLeg* fxUnderLeg = NULL;
	ARM_SwapLeg* fxNumLeg = NULL;

	int isFirstPeriodLeg = 0;

    ARM_Date swapStartDate;
	ARM_Date fixedSwapEndDate;

	ARM_Date firstNoticeDate;
	ARM_Date noticeDate;
	ARM_Vector* vNoticeDates = NULL;
	ARM_Vector* vCancelDates = NULL;
	ARM_Vector* vPaymentDates = NULL;
	ARM_Vector* vFx = NULL;
	ARM_Vector* vCap = NULL;
	ARM_Vector* vFloor = NULL;

	ARM_ReferenceValue* rCap = NULL;
	ARM_ReferenceValue* rFloor = NULL;
	ARM_ReferenceValue* rFx = NULL;

	ARM_PowerReverse* newPowerReverse = NULL;

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	CCString xmlResponse;
	CCString messageList;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting PRCS");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	long nbNodes;

	XMLDoc->selectSingleNode(_bstr_t((const char *)("Response/SWAPTION/Back/BACK")), &theNode);
	string status = GetStringFromXMLNode(theNode, "TermAssignStatus");
	theNode->Release();

	XMLDoc->selectSingleNode(_bstr_t((const char *)("Response/SWAPTION/Option/OPTION")), &theNode);
	string statusExercised = GetStringFromXMLNode(theNode, "Exercised");
	theNode->Release();

	// Recuperation du book et de la structure
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/SWAPTION/Env/ENV")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (nbNodes != 1)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		hr = resultList->get_item(0, &listItem);

		listItem->selectSingleNode(_bstr_t((const char *)"StructureId"), &theNode);

        if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			structureId = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		if (listItem)
			listItem->Release();

		if (resultList)
			resultList->Release();
	}

	/**********************************************************/
	
	vNoticeDates = getNotAndCancDates(XMLDoc,&vCancelDates,custId,dealId);
	/**********************************************************/

	initPeriodLeg = getInitPeriodLeg(XMLDoc,1);
	/**********************************************************/

	int dualFlag = 0;
	double dualStrike = 30.0;

	fxNumLeg = getStructuredLegs(XMLDoc,
								 &fxUnderLeg,
								 &rCap,
								 &rFloor,
								 &rFx,
								 &FxFixing,
								 &dualFlag,
								 &dualStrike,
								 bookName,
								 1);

	noticeDate = fxUnderLeg->GetResetDates()->Elt(fxUnderLeg->GetResetDates()->GetSize() - 1);

	realFundLeg = getFundingLeg(XMLDoc,&LiborFixing,1);

	newPowerReverse = new ARM_PowerReverse(initPeriodLeg,
										   realFundLeg,
										   fxUnderLeg,
										   fxNumLeg,
										   vNoticeDates,
										   vCancelDates,
										   rFx,
										   rCap,
										   rFloor,
										   dualFlag,
										   dualStrike,
										   noticeDate);

    double  dNotional = fxNumLeg->GetAmount()->CptReferenceValue((ARM_Date) fxNumLeg->GetPaymentDates()->Elt(0));
	ARM_ReferenceValue* rNotional = new ARM_ReferenceValue(dNotional);

	newPowerReverse->SetAmount(rNotional);
	delete rNotional;

	if ( (status == "TERM") || (statusExercised == "EXER") )
		newPowerReverse->SetFullterminated(true);

	XMLDoc->Release();
	XMLDoc = NULL;

	delete vNoticeDates;
	delete vCancelDates;

	delete fxNumLeg;
	delete fxUnderLeg;
	delete initPeriodLeg;
	delete realFundLeg;

	delete rFx;
	delete rCap;
	delete rFloor;

	// Gestion des fixings
	// Valeur par défaut MinDate : Pas de gestion de fixing
	if (date.isMinDate() == 0)
	{
		double today = date.GetJulian();

		// d'abord les fixings LIBOR
		if (LiborFixing != NULL) /// otherwise no fixing in summit
		{
			/// reset dates for the funding leg
			ARM_Vector* fundingLegResetDates = newPowerReverse->GetItsRealFundLeg()->GetResetDates();
			/// reset dates from summit until summit current date
			ARM_Vector* summitResetDates = LiborFixing->GetDiscreteDates();
			/// computation date
			if (fundingLegResetDates->Elt(0) <= today && summitResetDates->Elt(0) <= today) 
			{
				int size = 0;
				for (int i=0; i<MIN(fundingLegResetDates->size(), summitResetDates->size()); i++)
				{
					if ((fundingLegResetDates->Elt(i) == summitResetDates->Elt(i)) && 
						(fundingLegResetDates->Elt(i) <= today)) 
					{
						size++;
					}
					else {
						if (fundingLegResetDates->Elt(i) > today) 
						{
							break;
						}
						else 
						{
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								"Mismatch reset date"); 
						}
					}
				}

				ARM_Vector* summitResetValues = LiborFixing->GetDiscreteValues();
				ARM_Vector fundingLegResetValues(size);

				for (i=0; i<size; i++) 
					fundingLegResetValues.InitElt(i, summitResetValues->Elt(i));

				newPowerReverse->GetItsRealFundLeg()->SetFixRates(&fundingLegResetValues);
			}
		}

		// Puis les fixings FX
		int rg = 0;

		ARM_Vector* fxResetDates = newPowerReverse->GetItsFxUnderLeg()->GetResetDates();
		if (fxResetDates->Elt(0) <= today) 
		{			
			for(int i=1; i<fxResetDates->size(); i++)
			{
				if ((fxResetDates->Elt(i-1) <= today) && (fxResetDates->Elt(i) > today)) 
				{
					rg = i-1; /// find the reset date for the fx side
					break; //i = fxResetDates->size();
				}
			}			
			
			ARM_Date FxResetDate(fxResetDates->Elt(rg));			
			bool isFound(false);

			if (FxFixing)
			{
				for(int i=0; i<FxFixing->GetSize(); i++)
				{
					if (ARM_Date(FxFixing->GetDiscreteDates()->Elt(i)) == FxResetDate) 
					{
						rg = i; /// find the associated reset date for the fx side in the table 
						isFound = true;
						break;
					}
				}
			}

			if ((!isFound) && (FxResetDate != date)) 
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"No reset date associated to the FX reset date");
			}

			double fixing(0.0);

			if (isFound) 
			{
				fixing = FxFixing->GetDiscreteValues()->Elt(rg);
			}

			newPowerReverse->SetLastFxFixing(fixing);
		}
	}

	if (LiborFixing)
	   delete LiborFixing;
	LiborFixing = NULL;

	if (FxFixing)
	   delete FxFixing;
	FxFixing = NULL;

	return(newPowerReverse);
}

///////////////////////////////////////////////////////////////////////
// routine : ARMLOCAL_ParseUnderlyingPRCS
// action : parse un power reverse ayant uniquement le strip de FX options
//          on crée aussi une funding fictive (fixleg a taux 0)
//          il n'y a pas d'option d'annulation
///////////////////////////////////////////////////////////////////////
ARM_PowerReverse* ARMLOCAL_ParseUnderlyingPRCS ( const char* chaineXML,
												 const ARM_Date& date,
												 CCString& bookName,
												 CCString& structureId,
												 CCString& custId,
												 CCString& dealId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;

	ARM_ReferenceValue* LiborFixing = NULL;
	ARM_ReferenceValue* FxFixing = NULL;

	ARM_SwapLeg* initPeriodLeg = NULL;
	ARM_SwapLeg* realFundLeg = NULL;
	ARM_SwapLeg* fxUnderLeg = NULL;
	ARM_SwapLeg* fxNumLeg = NULL;

	int isFirstPeriodLeg = 0;

	CCString desk;

    ARM_Date startDate;
	ARM_Date endDate;
	ARM_Currency* currency;

	ARM_Date firstNoticeDate;
	ARM_Date noticeDate;
	ARM_Vector* vNoticeDates = NULL;
	ARM_Vector* vCancelDates = NULL;
	ARM_Vector* vPaymentDates = NULL;
	ARM_Vector* vFx = NULL;
	ARM_Vector* vCap = NULL;
	ARM_Vector* vFloor = NULL;

	ARM_ReferenceValue* rCap = NULL;
	ARM_ReferenceValue* rFloor = NULL;
	ARM_ReferenceValue* rFx = NULL;

	ARM_PowerReverse* newPowerReverse = NULL;

	CCString xmlResponse;
	CCString messageList;

	vector<string> listXMLResponse;


	// Récupération de la structureId
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting PRCS");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	long nbNodes;

	// Recuperation du book et de la structure
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Env/ENV")), &resultList) == S_OK)
	{
		resultList->get_length(&nbNodes);

		if (nbNodes != 1)
		{
			hr = S_FALSE;

			CCString msg((CCString)"Invalid XML string for getting PRCS");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (char*) msg);
		}

		hr = resultList->get_item(0, &listItem);

		listItem->selectSingleNode(_bstr_t((const char *)"StructureId"), &theNode);

        if ( theNode != NULL )
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			structureId = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}

		listItem->selectSingleNode(_bstr_t((const char *)"Desk"), &theNode);
		if (theNode!=NULL)
		{
			BSTR resultat = NULL;
			theNode->get_text(&resultat);

			_bstr_t ff(resultat,false);
			char * ff1=(char *)ff;

			desk = ff1;

			theNode->Release();
			theNode=NULL;
			if (resultat) SysFreeString(resultat);
		}


		listItem->Release();
		listItem=NULL;
	}

	// Recuperation de données pour la funding leg
	if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
	{
		hr = resultList->get_item(0, &listItem);
		startDate = GetStartDate(listItem);
		endDate = GetEndDate(listItem);
		currency = GetCcy(listItem);
	}
	
	XMLDoc->Release();
	XMLDoc = NULL;

	resultList->Release();
	resultList=NULL;

	listXMLResponse = etoolkit_getTradeListByStructId(structureId,desk);

	// Recuperation de l'Exotic qui contient le FX_OPT_STRIP
	MSXML2::IXMLDOMDocument* XMLDocForFxOptStrip = searchDocForExoticFxOptStrip(listXMLResponse);

	initPeriodLeg = getInitPeriodLeg(XMLDocForFxOptStrip,0);

	// Création d'une funding : on n'a besoin que de Start, End et Currency
	realFundLeg = new ARM_FixLeg(startDate, endDate, 0.0);
	realFundLeg->SetCurrencyUnit(currency);

	int dualFlag = 0;
	double dualStrike = 30.0;

	fxNumLeg = getStructuredLegs(XMLDocForFxOptStrip,
								 &fxUnderLeg,
								 &rCap,
								 &rFloor,
								 &rFx,
								 &FxFixing,
								 &dualFlag,
								 &dualStrike,
								 bookName,
								 0);

	noticeDate = fxUnderLeg->GetResetDates()->Elt(fxUnderLeg->GetResetDates()->GetSize() - 1);

	newPowerReverse = new ARM_PowerReverse(initPeriodLeg,
										   realFundLeg,
										   fxUnderLeg,
										   fxNumLeg,
										   vNoticeDates,
										   vCancelDates,
										   rFx,
										   rCap,
										   rFloor,
										   dualFlag,
										   dualStrike,
										   noticeDate);

    double  dNotional = fxNumLeg->GetAmount()->CptReferenceValue((ARM_Date) fxNumLeg->GetPaymentDates()->Elt(0));
	
    ARM_ReferenceValue* rNotional = new ARM_ReferenceValue(dNotional);

	newPowerReverse->SetAmount(rNotional);
	delete rNotional;

	XMLDocForFxOptStrip->Release();
	XMLDocForFxOptStrip = NULL;

	delete vNoticeDates;
	delete vCancelDates;

	delete fxNumLeg;
	delete fxUnderLeg;
	delete initPeriodLeg;
	delete realFundLeg;

    if (rFx)
       delete rFx;
    
    if (rCap)
       delete rCap;

    if (rFloor)
       delete rFloor;

	// Gestion des fixings
	// Valeur par défaut MinDate : Pas de gestion de fixing
	if ( date.isMinDate() == 0 )
	{
		double today = date.GetJulian();

		// d'abord les fixings LIBOR
		if ( LiborFixing != NULL ) /// otherwise no fixing in summit
		{
			/// reset dates for the funding leg
			ARM_Vector* fundingLegResetDates = newPowerReverse->GetItsRealFundLeg()->GetResetDates();
			/// reset dates from summit until summit current date
			ARM_Vector* summitResetDates = LiborFixing->GetDiscreteDates();
			/// computation date
			if ( fundingLegResetDates->Elt(0) <= today && summitResetDates->Elt(0) <= today ) 
			{
				int size = 0;
				for (int i = 0; i < MIN(fundingLegResetDates->size(), summitResetDates->size()); i++)
				{
					if ((fundingLegResetDates->Elt(i) == summitResetDates->Elt(i)) && 
						(fundingLegResetDates->Elt(i) <= today)) 
					{
						size++;
					}
					else 
                    {
						if (fundingLegResetDates->Elt(i) > today) 
						{
							break;
						}
						else 
						{
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								"Mismatch reset date"); 
						}
					}
				}

				ARM_Vector* summitResetValues = LiborFixing->GetDiscreteValues();
				ARM_Vector fundingLegResetValues(size);

				for (i = 0; i < size; i++) 
				{
					fundingLegResetValues.InitElt(i, summitResetValues->Elt(i));
				}

				newPowerReverse->GetItsRealFundLeg()->SetFixRates(&fundingLegResetValues);
			}
		}

		// Puis les fixings FX
		int rg = 0;

		ARM_Vector* fxResetDates = newPowerReverse->GetItsFxUnderLeg()->GetResetDates();

        if ( fxResetDates->Elt(0) <= today ) 
		{			
			for (int i = 1; i < fxResetDates->size(); i++)
			{
				if (( fxResetDates->Elt(i-1) <= today ) && ( fxResetDates->Elt(i) > today )) 
				{
					rg = i-1; /// find the reset date for the fx side
					break; //i = fxResetDates->size();
				}
			}			
			
			ARM_Date FxResetDate(fxResetDates->Elt(rg));			
			bool isFound(false);

			if (FxFixing)
			{
				for ( i = 0; i < FxFixing->GetSize(); i++)
				{
					if ( ARM_Date(FxFixing->GetDiscreteDates()->Elt(i)) == FxResetDate ) 
					{
						rg = i; /// find the associated reset date for the fx side in the table 
						isFound = true;
						break;
					}
				}
			}

			if ((!isFound) && ( FxResetDate != date )) 
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                                "No reset date associated to the FX reset date");
			}

			double fixing(0.0);

			if (isFound) 
			{
				fixing = FxFixing->GetDiscreteValues()->Elt(rg);
			}

			newPowerReverse->SetLastFxFixing(fixing);
		}
	}

	if (LiborFixing)
		delete LiborFixing;
	LiborFixing = NULL;

	if (FxFixing)
		delete FxFixing;
	FxFixing = NULL;

	return(newPowerReverse);
}


ARM_VolLInterpol* etoolkit_GetVolATMFromSummit(const CCString& index,
											   const CCString& currency,
											   const CCString& cvName,
											   ARM_Date date,
											   const CCString& vtype,
											   const CCString& impOrHist)
{
	ARM_result C_result;

	CCString xmlResponse;
	CCString msgList;

	long retCode = ARM_OK;

	ARM_Currency sCCY((const char*)currency);

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetVOLTenorsList((const char*)impOrHist,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		retCode = etoolkit_getcommsetname(index,
											   currency,
											   cvName,
											   vtype,
											   impOrHist,
											   date,
											   xmlResponse,
											   msgList);
	}

	if (( retCode == ARM_KO ) 
        || 
        ( strcmp((const char*)xmlResponse,"<Response></Response>") == 0 )
       )
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "etoolkit_GetVolATMFromSummit: no ATM Vol found!");

        return(NULL);
	}

	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listTenors;

	if ( vtype == "ISSUER" ) 
	   listTenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,&listRequetes, "P4");
	else 
	   listTenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,&listRequetes, "P5");

	ARM_Vector* vStrikes = new ARM_Vector(listTenors.size());

	int i;
	int j;

	for (i = 0; i < listTenors.size(); i++)
	{
	    double resMatu = 0.0;
	    int inMatu;

		char cMatu = ' ';
		long freqId;


		sscanf(listTenors[i], "%d%c", &inMatu, &cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
		   freqId = K_DAILY;
		else if ( cMatu == 'W' )  
		   freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
		   freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
		   freqId = K_ANNUAL;

        // Just compute year terms
        // Rq: asOfDate is used only if contract futurs!

        vStrikes->Elt(i) = FromStrMatuToDouble((const char *) listTenors[i], &date);
    }

	VECTOR<CCString> sortListTenors;
	VECTOR<CCString> sortRequetes;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	for (i = 0; i < vStrikes->GetSize(); i++)
	{
		j = 0;
		while ( j < vStrikes->GetSize() )
		{
			if ( vSortStrikes->Elt(i) != vStrikes->Elt(j) )
			{
			   j++;
			}
			else
			{
				sortListTenors.push_back(listTenors[j]);
				sortRequetes.push_back(listRequetes[j]);
			
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	CCString typeForGVC = "FWVL";

	if ( impOrHist == "HISTVOL" )
	   typeForGVC = "HISTVOL";

	// ** GIGASPACE
	ARM_XGigaToolKit::doGetVOLTenor((const char*)impOrHist,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)sortListTenors[0],date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		//SUMMIT
		retCode = etoolkit_getVolCurveByTenor(sortRequetes[0],
											  cvName,
											  date,
											  typeForGVC,
											  impOrHist,
											  xmlResponse,
											  msgList);
	}

	VECTOR<CCString> listYearTerms;

	retCode = ARMLOCAL_GetDatesAndDatasForVol(xmlResponse,
											  date,
											  0,
											  mVolatility,
											  vYearTerms,
											  (const char*)currency,
											  sortListTenors.size(),
											  &listYearTerms);

	for (i = 1; i < sortListTenors.size(); i++)
	{
		// ** GIGASPACE
		ARM_XGigaToolKit::doGetVOLTenor((const char*)impOrHist,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)sortListTenors[i],date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getVolCurveByTenor(sortRequetes[i],
												  cvName,
												  date,
												  typeForGVC,
												  impOrHist,
												  xmlResponse,
												  msgList);
		}
		retCode = ARMLOCAL_GetDatesAndDatasForVol(xmlResponse,
												  date,
												  i,
												  mVolatility,
												  vYearTerms,
												  (const char*)currency);
	}
	
	bool	isSorted = true;	
	size_t	sizeYearTerms = vYearTerms->GetSize();	
	
	VECTOR<CCString> NewlistYearTerms;
	
	for( i = 0; i < sizeYearTerms - 1; i++)	
	{
		if (vYearTerms->Elt(i) > vYearTerms->Elt(i+1) )			
		{
			isSorted = false;
			break;
		}
	}
	
	ARM_VolLInterpol* newVolCrv = NULL;
		
	if(!isSorted)		
	{		
		ARM_Vector* vSortYearTerms = new ARM_Vector(vYearTerms->Sort());		
			
		int nbCol = mVolatility->GetNumCols();
		int nbRow = mVolatility->GetNumLines();
				
		ARM_Matrix* mNewVolatility = new ARM_Matrix(nbRow, nbCol);
				
		if(nbRow != sizeYearTerms)
		{
			delete vSortYearTerms;
			delete mNewVolatility;	
			delete vStrikes;
			delete vYearTerms;
			delete mVolatility;
			
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							"etoolkit_GetVolATMFromSummit: nb rows in Vol matrix should be equal to year term vector size !");	
		}
		
		int oldIdx;
		for(i=0; i<sizeYearTerms; i++)		
		{			
			oldIdx = i;
			if( vSortYearTerms->Elt(i) != vYearTerms->Elt(i))
				oldIdx = vYearTerms->find(vSortYearTerms->Elt(i));			
			
			NewlistYearTerms.push_back(listYearTerms[oldIdx]);
			
			for(j=0; j<nbCol; j++)
				mNewVolatility->Elt(i,j) = mVolatility->Elt(oldIdx,j);
			
		}

		delete mVolatility;
		delete vYearTerms;
		
		newVolCrv = new ARM_VolLInterpol(date, vSortYearTerms, vSortStrikes, mNewVolatility, 1, K_ATMF_VOL);		
	}
	
	else		
	{
		newVolCrv = new ARM_VolLInterpol(date, vYearTerms,vSortStrikes, mVolatility, 1, K_ATMF_VOL);		
	}

	newVolCrv->SetInterpType(1);
	newVolCrv->SetCurrencyUnit(&sCCY);
	newVolCrv->SetIndexName((char*)(const char*)index);

	// Update Mkt data characteristics
	string	vVolType = vtype;
	if(vVolType == "IRG")
		vVolType = "CAP";
    
    string	vType = string("VOL ") + vVolType;
    string	vIndex((const char*) index);
	string	vCurrency((const char*) currency);
	string	vCrvId((const char*) cvName);

    if( ( index == "ROLIB" ) || ( index == "ROEUR" ) )
    {
		vIndex = "RO";
    }
    else if( ( index == "NULIB" ) || ( index == "NUEUR" ) )
    {
		vIndex = "NU";
    }
	else if( vIndex.substr(0, 2) == "CO" )
	{
		vType = "IRIR CORR";
		vIndex = "";
	}

    newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	for (i = 0; i < ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}
	
	if(isSorted)	
	{
		for (i = 0; i < listYearTerms.size(); i++)		
		{           
			strcpy(newVolCrv->itsYearTermsX[i], (const char *) listYearTerms[i]);
		}
	}
	else
	{
		for (i = 0; i < listYearTerms.size(); i++)			
		{           
			strcpy(newVolCrv->itsYearTermsX[i], (const char *) NewlistYearTerms[i]);
		}
	}

	for (i = 0; i < sortListTenors.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsY[i], (const char *) sortListTenors[i]);
	}

    // clean up
	if (vStrikes)
	   delete vStrikes;
	vStrikes = NULL;

	return(newVolCrv);
}

/*
ARM_VolLInterpol* etoolkit_GetVolATMFromSummit(const CCString& index,
											   const CCString& currency,
											   const CCString& cvName,
											   ARM_Date date,
											   const CCString& vtype,
											   const CCString& impOrHist)
{
	ARM_result C_result;

	ARM_Currency sCCY((const char*)currency);

	CCString xmlResponse;
	CCString msgList;

	long retCode;

	if (vtype =="ISSUER")
	{
		retCode = etoolkit_getcommsetname(index,
										  currency,
										  cvName,
										  vtype,
										  impOrHist,
										  date,
										  xmlResponse,
										  msgList);
	}
	else
	{
		retCode = etoolkit_newgetcommsetname(index,
											 currency,
											 cvName,
											 vtype,
											 impOrHist,
											 date,
											 xmlResponse,
											 msgList);
	}

	if (( retCode == ARM_KO ) 
        || 
        ( strcmp((const char*)xmlResponse,"<Response></Response>") == 0 )
       )
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "etoolkit_GetVolATMFromSummit: no ATM Vol found!");

        return(NULL);
	}

	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listTenors;

	if ( vtype == "ISSUER" ) 
	   listTenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,&listRequetes, "P4");
	else 
	   listTenors = ARMLOCAL_GetListTenorsForVol(xmlResponse,false);

	ARM_Vector* vStrikes = new ARM_Vector(listTenors.size());

	int i;

	for (i = 0; i < listTenors.size(); i++)
	{
	    double resMatu = 0.0;
	    int inMatu;

		char cMatu = ' ';
		long freqId;


		sscanf(listTenors[i], "%d%c", &inMatu, &cMatu);

		if ( cMatu == 'D' ) // Ex : "1D"
		   freqId = K_DAILY;
		else if ( cMatu == 'W' )  
		   freqId = K_WEEKLY;
		else if ( cMatu == 'M' ) 
		   freqId = K_MONTHLY;
		else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
		   freqId = K_ANNUAL;

        // Just compute year terms
        // Rq: asOfDate is used only if contract futurs!

        vStrikes->Elt(i) = FromStrMatuToDouble((const char *) listTenors[i], &date);
    }

	VECTOR<CCString> sortListTenors;
	VECTOR<CCString> sortRequetes;
	VECTOR<int> sortListTenorsId;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	for (i = 0; i < vStrikes->GetSize(); i++)
	{
		int j = 0;
		while ( j < vStrikes->GetSize() )
		{
			if ( vSortStrikes->Elt(i) != vStrikes->Elt(j) )
			{
			   j++;
			}
			else
			{
				sortListTenors.push_back(listTenors[j]);
				sortListTenorsId.push_back(j);
				if ( vtype == "ISSUER" ) 
					sortRequetes.push_back(listRequetes[j]);
			
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	CCString typeForGVC = "FWVL";

	if ( impOrHist == "HISTVOL" )
	   typeForGVC = "HISTVOL";

	ARM_VolLInterpol* newVolCrv = NULL;
	VECTOR<CCString> listYearTerms;
	VECTOR<CCString> sortYearTerms;
	VECTOR<int> sortYearTermsId;

	if ( vtype == "ISSUER" )
	{
		retCode = etoolkit_getVolCurveByTenor(sortRequetes[0],
											  cvName,
											  date,
											  typeForGVC,
											  impOrHist,
											  xmlResponse,
											  msgList);

		retCode = ARMLOCAL_GetDatesAndDatasForVol(xmlResponse,
												  date,
												  0,
												  mVolatility,
												  vYearTerms,
												  (const char*)currency,
												  sortListTenors.size(),
												  &listYearTerms);

		for (i = 1; i < sortListTenors.size(); i++)
		{
			retCode = etoolkit_getVolCurveByTenor(sortRequetes[i],
												  cvName,
												  date,
												  typeForGVC,
												  impOrHist,
												  xmlResponse,
												  msgList);

			retCode = ARMLOCAL_GetDatesAndDatasForVol(xmlResponse,
													  date,
													  i,
													  mVolatility,
													  vYearTerms,
													  (const char*)currency);
		}

		newVolCrv = new ARM_VolLInterpol(date, vYearTerms,
										 vSortStrikes, mVolatility, 1, K_ATMF_VOL);
	}
	else
	{
		// 1- Construction de la matrice
		retCode = ARMLOCAL_NewGetDatesAndDatasForVol(xmlResponse,
													 date,
													 mVolatility,
													 vYearTerms,
													 (const char*)currency,
													 &listYearTerms);
		// 2- Tri sur les expiries (par ligne)
		ARM_Vector* vSortYearTerms = new ARM_Vector(vYearTerms->Sort());

		for (i = 0; i < vYearTerms->GetSize(); i++)
		{
			int j = 0;
			while ( j < vYearTerms->GetSize() )
			{
				if ( vSortYearTerms->Elt(i) != vYearTerms->Elt(j) )
				{
				   j++;
				}
				else
				{
					sortYearTerms.push_back(listYearTerms[j]);
					sortYearTermsId.push_back(j);

					break;
				}
			}
		}

		// 3- reconstruction de la matrice de vol
		ARM_Matrix* mSortVol = new ARM_Matrix(vSortYearTerms->GetSize(),vSortStrikes->GetSize());
		for (i=0;i<vYearTerms->GetSize();i++)
			for(int j=0;j<vStrikes->GetSize();j++)
				mSortVol->Elt(i,j)=mVolatility->Elt(sortYearTermsId[i],sortListTenorsId[j]);

		// 4- Reconstruction de l'objet
		newVolCrv = new ARM_VolLInterpol(date, vSortYearTerms,
										 vSortStrikes, mSortVol, 1, K_ATMF_VOL);

		if (mVolatility)
			delete mVolatility;
		mVolatility = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;
	}

	newVolCrv->SetInterpType(1);
	newVolCrv->SetCurrencyUnit(&sCCY);
	newVolCrv->SetIndexName((char*)(const char*)index);

	// Update Mkt data characteristics
	string	vVolType = vtype;
	if(vVolType == "IRG")
		vVolType = "CAP";
    
    string	vType = string("VOL ") + vVolType;
    string	vIndex((const char*) index);
	string	vCurrency((const char*) currency);
	string	vCrvId((const char*) cvName);

    if( ( index == "ROLIB" ) || ( index == "ROEUR" ) )
    {
		vIndex = "RO";
    }
    else if( ( index == "NULIB" ) || ( index == "NUEUR" ) )
    {
		vIndex = "NU";
    }
	else if( vIndex.substr(0, 2) == "CO" )
	{
		vType = "IRIR CORR";
		vIndex = "";
	}

    newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	for (i = 0; i < ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}

	for (i = 0; i < listYearTerms.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsX[i], (const char *) listYearTerms[i]);
	}

	for (i = 0; i < sortListTenors.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsY[i], (const char *) sortListTenors[i]);
	}

    // clean up
	if (vStrikes)
	   delete vStrikes;
	vStrikes = NULL;

	return(newVolCrv);
}
*/


/*
ARM_VolLInterpol* etoolkit_GetSmileFromSummit(const CCString& index,
											  const CCString& currency,
											  const CCString& cvName,
											  ARM_Date date,
											  const CCString& vtype,
											  const CCString& matuIndex,
											  const CCString& impOrHist)
{
	ARM_VolLInterpol* newVolCrv = NULL;

    ARM_Currency sCCY((const char*)currency);
	
	ARM_result C_result;

	CCString xmlResponse;
	CCString msgList;

	long retCode = etoolkit_newgetcommsetname(index,
											  currency,
											  cvName,
											  vtype,
											  "SMILE",
											  date,
											  xmlResponse,
											  msgList,
											  matuIndex);

	if (retCode == ARM_KO)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "no Smile");
	}
	
	VECTOR<CCString> listStrikes;
	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listYearTerms;

	try
	{
		listStrikes = ARMLOCAL_GetListTenorsForVol(xmlResponse,true);
	}

    catch (Exception& e)
	{
		e.DebugPrint();

		throw e;
	}
	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "Pb in getting List of Strikes for smile");
	}
	
	ARM_Vector* vStrikes = new ARM_Vector(listStrikes.size());
	
	int i;

	for (i=0;i<vStrikes->GetSize();i++)
		vStrikes->Elt(i) = atof(listStrikes[i]) / 100.;

	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortListStrikes;
	VECTOR<int> sortListStrikesId;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	for (i=0;i<vStrikes->GetSize();i++)
	{
		int j=0;
		while (j<vStrikes->GetSize())
		{
			if (vSortStrikes->Elt(i) != vStrikes->Elt(j))
			{
				j++;
			}
			else
			{
				sortListStrikesId.push_back(j);
				sortRes.push_back(listStrikes[j]);

				double strike = atof(listStrikes[j]) / 100.;
				char s[128];
				sprintf(s, "%f", strike);
				sortListStrikes.push_back(s);
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	VECTOR<CCString> sortYearTerms;
	VECTOR<int> sortYearTermsId;

	// 1- Construction de la matrice
	retCode = ARMLOCAL_NewGetDatesAndDatasForVol(xmlResponse,
												 date,
												 mVolatility,
												 vYearTerms,
												 (const char*)currency,
												 &listYearTerms);
	// 2- Tri sur les expiries (par ligne)
	ARM_Vector* vSortYearTerms = new ARM_Vector(vYearTerms->Sort());

	for (i = 0; i < vYearTerms->GetSize(); i++)
	{
		int j = 0;
		while ( j < vYearTerms->GetSize() )
		{
			if ( vSortYearTerms->Elt(i) != vYearTerms->Elt(j) )
			{
			   j++;
			}
			else
			{
				sortYearTerms.push_back(listYearTerms[j]);
				sortYearTermsId.push_back(j);

				break;
			}
		}
	}

	// 3- reconstruction de la matrice de vol
	ARM_Matrix* mSortVol = new ARM_Matrix(vSortYearTerms->GetSize(),vSortStrikes->GetSize());
	for (i=0;i<vYearTerms->GetSize();i++)
		for(int j=0;j<vStrikes->GetSize();j++)
			mSortVol->Elt(i,j)=mVolatility->Elt(sortYearTermsId[i],sortListStrikesId[j]);

	// 4- Construction de l'objet
    newVolCrv = new ARM_VolLInterpol( date, vSortYearTerms,
								vSortStrikes, mSortVol, 1, K_SMILE_VOL);

	if (mVolatility)
		delete mVolatility;
	mVolatility = NULL;

	newVolCrv->SetCurrencyUnit(&sCCY);
	newVolCrv->SetIndexName((char*)(const char*)index);

	// Update Mkt data characteristics
	string	vVolType = vtype;
	if(vVolType == "IRG")
		vVolType = "CAP";
    
    string	vType = string("VOL ") + vVolType;
    string	vIndex((const char*) index);
	string	vCurrency((const char*) currency);
	string	vCrvId((const char*) cvName);

    if( ( index == "ROLIB" ) || ( index == "ROEUR" ) )
    {
		vIndex = "RO";
    }
    else if( ( index == "NULIB" ) || ( index == "NUEUR" ) )
    {
		vIndex = "NU";
    }
	else if( vIndex.substr(0, 2) == "CO" )
	{
		vType = "IRIR CORR";
		vIndex = "";
	}

    newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	for (i=0; i<ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}

	for (i=0; i<listYearTerms.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsX[i], (const char *) listYearTerms[i]);
	}

	for (i=0; i<sortListStrikes.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsY[i], (const char*) sortListStrikes[i]);
	}

	if (vStrikes)
		delete vStrikes;
	vStrikes = NULL;

	return newVolCrv;
}
*/

ARM_VolLInterpol* etoolkit_GetSmileFromSummit(const CCString& index,
											  const CCString& currency,
											  const CCString& cvName,
											  ARM_Date date,
											  const CCString& vtype,
											  const CCString& matuIndex,
											  int smileType,
											  const CCString& impOrHist)
{
	ARM_VolLInterpol* newVolCrv = NULL;

    ARM_Currency sCCY((const char*)currency);
	
	ARM_result C_result;

	CCString xmlResponse;
	CCString msgList;

	long retCode = ARM_OK;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetSMILETenorStrikeList((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)matuIndex,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		retCode = etoolkit_getcommsetname(index,
										  currency,
										  cvName,
										  vtype,
										  "SMILE",
										  date,
										  xmlResponse,
										  msgList,
										  matuIndex);
	}

	if (retCode == ARM_KO)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "no Smile");
	}
	
	VECTOR<CCString> listStrikes;
	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listYearTerms;

	try
	{
		listStrikes = ARMLOCAL_GetListStrikesFromXML(xmlResponse,&listRequetes);
	}

    catch (Exception& e)
	{
		e.DebugPrint();

		throw e;
	}
	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "Pb in getting List of Strikes for smile");
	}
	
	ARM_Vector* vStrikes = new ARM_Vector(listStrikes.size());
	
	int i;

	for (i=0;i<vStrikes->GetSize();i++)
		vStrikes->Elt(i) = atof(listStrikes[i]) / 100.;

	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortListStrikes;
	VECTOR<CCString> sortGSListStrikes;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	int j;

	for (i=0; i<vStrikes->GetSize(); i++)
	{
		j=0;
		while (j<vStrikes->GetSize())
		{
			if (vSortStrikes->Elt(i) != vStrikes->Elt(j))
			{
				j++;
			}
			else
			{
				sortRes.push_back(listStrikes[j]);

				double strike = atof(listStrikes[j]) / 100.;
				char s[128];
				sprintf(s, "%f", strike);
				sortListStrikes.push_back(s);

				// Pour Gigaspaces
				double strikeGiga = atof(listStrikes[j]);
				char sGiga[128];
				sprintf(sGiga, "%.0lf", strikeGiga);
				sortGSListStrikes.push_back(sGiga);
				break;
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	try
	{
		ARM_XGigaToolKit::doGetSMILETenorStrike((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)matuIndex,(const char*)sortGSListStrikes[0],date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getSmileCurveByStrike(sortRes[0],
													 listRequetes[0],
													 cvName,
													 date,
													 "SMILE",
													 xmlResponse,
													 msgList);
		}

		retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
													currency,
													date,
													0,
													mVolatility,
													vYearTerms,
													sortRes.size(),
													1,
													&listYearTerms,
													impOrHist);
		
		for (i=1;i<sortRes.size();i++)
		{
			ARM_XGigaToolKit::doGetSMILETenorStrike((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)matuIndex,(const char*)sortGSListStrikes[i],date,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				retCode = etoolkit_getSmileCurveByStrike(sortRes[i],
														 listRequetes[i],
														 cvName,
														 date,
														 "SMILE",
														 xmlResponse,
														 msgList);
			}

			retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
														currency,
														date,
														i,
														mVolatility,
														vYearTerms);
		}
	}
	catch (Exception& e)
	{
		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		if (vSortStrikes)
			delete vSortStrikes;
		vSortStrikes = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		if (mVolatility)
			delete mVolatility;
		mVolatility = NULL;

		throw e;
	}
	catch (...)
	{
		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		if (vSortStrikes)
			delete vSortStrikes;
		vSortStrikes = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		if (mVolatility)
			delete mVolatility;
		mVolatility = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "Pb in getting Datas for smile");
	}
	
	bool	isSorted = true;	
	size_t	sizeYearTerms = vYearTerms->GetSize();	
	
	VECTOR<CCString> NewlistYearTerms;
	
	for( i = 0; i < sizeYearTerms - 1; i++)	
	{
		if (vYearTerms->Elt(i) > vYearTerms->Elt(i+1) )			
		{
			isSorted = false;
			break;
		}
	}
	
	if(!isSorted)		
	{		
		ARM_Vector* vSortYearTerms = new ARM_Vector(vYearTerms->Sort());		
			
		int nbCol = mVolatility->GetNumCols();
		int nbRow = mVolatility->GetNumLines();
				
		ARM_Matrix* mNewVolatility = new ARM_Matrix(nbRow, nbCol);
				
		if(nbRow != sizeYearTerms)
		{
			delete vSortYearTerms;
			delete mNewVolatility;	
			delete vStrikes;
			delete vYearTerms;
			delete mVolatility;
			
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							"etoolkit_GetSmileFromSummit: nb rows in Vol matrix should be equal to year term vector size !");
		}
		
		int oldIdx;
		for(i=0; i<sizeYearTerms; i++)		
		{			
			oldIdx = i;
			if( vSortYearTerms->Elt(i) != vYearTerms->Elt(i))
				oldIdx = vYearTerms->find(vSortYearTerms->Elt(i));			
			
			NewlistYearTerms.push_back(listYearTerms[oldIdx]);
			
			for(j=0; j<nbCol; j++)
				mNewVolatility->Elt(i,j) = mVolatility->Elt(oldIdx,j);
			
		}

		delete mVolatility;
		delete vYearTerms;
		
		if(smileType == 0)
			newVolCrv = new ARM_VolLInterpol(date, vSortYearTerms, vSortStrikes, mNewVolatility, 1, K_SMILE_VOL);		
		else
			newVolCrv = new ARM_VolSplineInterpol(date, vSortYearTerms, vSortStrikes, mNewVolatility, 1, K_SMILE_VOL);		
	}
	
	else		
	{
		if(smileType == 0)
			newVolCrv = new ARM_VolLInterpol(date, vYearTerms,vSortStrikes, mVolatility, 1, K_SMILE_VOL);		
		else
			newVolCrv = new ARM_VolSplineInterpol(date, vYearTerms,vSortStrikes, mVolatility, 1, K_SMILE_VOL);		
	}

	newVolCrv->SetCurrencyUnit(&sCCY);
	newVolCrv->SetIndexName((char*)(const char*)index);

	// Update Mkt data characteristics
	string	vVolType = vtype;
	if(vVolType == "IRG")
		vVolType = "CAP";
    
    string	vType = string("VOL ") + vVolType;
    string	vIndex((const char*) index);
	string	vCurrency((const char*) currency);
	string	vCrvId((const char*) cvName);

    if( ( index == "ROLIB" ) || ( index == "ROEUR" ) )
    {
		vIndex = "RO";
    }
    else if( ( index == "NULIB" ) || ( index == "NUEUR" ) )
    {
		vIndex = "NU";
    }
	else if( vIndex.substr(0, 2) == "CO" )
	{
		vType = "IRIR CORR";
		vIndex = "";
	}

    newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	for (i=0; i<ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}

	if(isSorted)	
	{
		for (i = 0; i < listYearTerms.size(); i++)		
		{           
			strcpy(newVolCrv->itsYearTermsX[i], (const char *) listYearTerms[i]);
		}
	}
	else
	{
		for (i = 0; i < listYearTerms.size(); i++)			
		{           
			strcpy(newVolCrv->itsYearTermsX[i], (const char *) NewlistYearTerms[i]);
		}
	}

//	for (i=0; i<listYearTerms.size(); i++)
//	{	
//		strcpy(newVolCrv->itsYearTermsX[i], (const char *) listYearTerms[i]);
//	}

	for (i=0; i<sortListStrikes.size(); i++)
	{	
		strcpy(newVolCrv->itsYearTermsY[i], (const char*) sortListStrikes[i]);
	}

	if (vStrikes)
		delete vStrikes;
	vStrikes = NULL;

	return newVolCrv;
}


long etoolkit_GetInitialVolFromSummit (const CCString& index,
									   const CCString& currency,
									   const CCString& cvName,
									   double date,
									   const CCString& vtype,
									   const CCString& matuIndex,
									   VECTOR<CCString> *maturities,
									   VECTOR<CCString> *tenors,
									   VECTOR<double> *vol)
{
	ARM_result C_result;

	char strDate[11];
	Local_XLDATE2ARMDATE(date,strDate);
	ARM_Date myDate(strDate);

	CCString xmlResponse;
	CCString msgList;

	CCString typeData;
	char smileMat[20];

	if ( strcmp(matuIndex,"ATM") == 0 )
	{
	   typeData = "IRFWDVOL";
	}
	else if (( matuIndex[0] == 'S' ) || ( matuIndex[0] == 's' ))
	{
		char onechar;
		char dummy[20];

		sscanf(matuIndex,"%[^:]%c%s",dummy,&onechar,smileMat);

		typeData = "SMILE";
	}
	else
		return ARM_KO;

	long retCode = ARM_OK;
	std::string xmlOutput ;

	if (strcmp((const char*)matuIndex,"ATM") == 0)
	{
		// ** GIGASPACE
		ARM_XGigaToolKit::doGetVOLTenorsList((const char*)typeData,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,myDate,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getcommsetname(index,
											  currency,
											  cvName,
											  vtype,
											  typeData,
											  myDate,
											  xmlResponse,
											  msgList);
		}
	}
	else
	{
		// ** GIGASPACE
		ARM_XGigaToolKit::doGetSMILETenorStrikeList((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)smileMat,myDate,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getcommsetname(index,
											  currency,
											  cvName,
											  vtype,
											  typeData,
											  myDate,
											  xmlResponse,
											  msgList,
											  smileMat);
		}
	}
	if ( (retCode == ARM_KO) || (strcmp(xmlResponse,"<Response></Response>") == 0) )
	{
		return ARM_KO;
	}

	VECTOR<CCString> listRequetes;

	try
	{
		if ( strcmp((const char*)matuIndex,"ATM") == 0 )
		   *tenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,&listRequetes,"P5");
		else
		   *tenors = ARMLOCAL_GetListStrikesFromXML(xmlResponse,&listRequetes);
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		CCString msg((CCString) "Error unknown in etoolkit_GetInitialVolFromSummit");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	ARM_Vector* vStrikes = new ARM_Vector(tenors->size());

	int i;

	for (i = 0; i < tenors->size(); i++)
	{	
		if ( strcmp((const char*)matuIndex,"ATM") == 0)
		{
           vStrikes->Elt(i) = FromStrMatuToDouble((const char *) (*tenors)[i], &myDate);
		}
		else
		   vStrikes->Elt(i) = atof((*tenors)[i]);
	}

	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortRequetes;
	VECTOR<CCString> sortGSListStrikes;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	for (i = 0; i < vStrikes->GetSize(); i++)
	{
		int j = 0;
		
        while ( j < vStrikes->GetSize() )
		{
			if ( vSortStrikes->Elt(i) != vStrikes->Elt(j) )
			{
			   j++;
			}
			else
			{
				sortRes.push_back((*tenors)[j]);
			
                sortRequetes.push_back(listRequetes[j]);
				
				if (strcmp((const char*)matuIndex,"ATM") != 0)
				{
					// Pour Gigaspaces
					double strikeGiga = atof((*tenors)[j]);
					char sGiga[128];
					sprintf(sGiga, "%.0lf", strikeGiga);
					sortGSListStrikes.push_back(sGiga);
				}
                break;
			}
		}
	}

	*tenors = sortRes;

	ARM_Matrix* mVolatility = NULL;

	if (strcmp((const char*)matuIndex,"ATM") == 0)
	{
		ARM_XGigaToolKit::doGetVOLTenor((const char*)typeData,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)sortRes[0],myDate,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getSmileCurveByStrike(sortRes[0],
													 sortRequetes[0],
													 cvName,
													 myDate,
													 typeData,
													 xmlResponse,
													 msgList);
		}
	}
	else
	{
		ARM_XGigaToolKit::doGetSMILETenorStrike((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)smileMat,(const char*)sortGSListStrikes[0],myDate,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getSmileCurveByStrike(sortRes[0],
													 sortRequetes[0],
													 cvName,
													 myDate,
													 typeData,
													 xmlResponse,
													 msgList);
		}
	}

	retCode = ARMLOCAL_GetStringAndDatasForVol(xmlResponse,
											   myDate,
											   0,
											   mVolatility,
											   *maturities,
											   sortRes.size());

	ARM_Currency	vCurrency((const char*)currency);
	long	 vSpotDays = vCurrency.GetSpotDays();
    char*	 vPaymentCalendar = vCurrency.GetPayCalName(vCurrency.GetVanillaIndexType());
	ARM_Date vSettlementDate = myDate;
	vSettlementDate.NextBusinessDay(vSpotDays, vPaymentCalendar);

	bool	isSorted = true;	
	size_t	sizeMaturities = maturities->size();	
	ARM_Vector	vMaturities(sizeMaturities);

	for( i=0; i<sizeMaturities; i++)
	{
		vMaturities.Elt(i) = convPlotInYearTerm(maturities->at(i), myDate, vSettlementDate, vPaymentCalendar);
	}

	delete vPaymentCalendar;
		
	for( i=0; i<sizeMaturities-1; i++)	
	{
		if (vMaturities.Elt(i) > vMaturities.Elt(i+1) )			
		{
			isSorted = false;
			break;
		}
	}

	for (i=1;i<sortRes.size();i++)
	{
		if (strcmp((const char*)matuIndex,"ATM") == 0)
		{
			ARM_XGigaToolKit::doGetVOLTenor((const char*)typeData,(const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)sortRes[i],myDate,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				retCode = etoolkit_getSmileCurveByStrike(sortRes[i],
														 sortRequetes[i],
														 cvName,
														 myDate,
														 typeData,
														 xmlResponse,
														 msgList);
			}
		}
		else
		{
			ARM_XGigaToolKit::doGetSMILETenorStrike((const char*)index,(const char*)currency,(const char*)cvName,(const char*)vtype,(const char*)smileMat,(const char*)sortGSListStrikes[i],myDate,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				retCode = etoolkit_getSmileCurveByStrike(sortRes[i],
														 sortRequetes[i],
														 cvName,
														 myDate,
														 typeData,
														 xmlResponse,
														 msgList);
			}
		}

		retCode = ARMLOCAL_GetStringAndDatasForVol(xmlResponse,
												   myDate,
												   i,
												   mVolatility,
												   *maturities);
	}

	VECTOR<CCString>	vNewMaturitiesList;
	if(!isSorted)		
	{		
		ARM_Vector	vSortedMaturities(vMaturities.Sort());

		int nbCol = mVolatility->GetNumCols();
		int nbRow = mVolatility->GetNumLines();
				
		ARM_Matrix* mNewVolatility = new ARM_Matrix(nbRow, nbCol);
		
		int oldIdx;
		for(i=0; i<sizeMaturities; i++)		
		{			
			oldIdx = i;
			if( vSortedMaturities.Elt(i) != vMaturities.Elt(i))
				oldIdx = vMaturities.find(vSortedMaturities.Elt(i));			
			
			vNewMaturitiesList.push_back((*maturities)[oldIdx]);
			
			for(int j=0; j<nbCol; j++)
				mNewVolatility->Elt(i,j) = mVolatility->Elt(oldIdx,j);
		}

		for (i = 0; i < mVolatility->GetNumLines(); i++)
			for (int j = 0; j < mVolatility->GetNumCols(); j++)
				vol->push_back(mNewVolatility->Elt(i,j));

		*maturities = vNewMaturitiesList;

		delete mNewVolatility;
	}
	else
	{
		for (i = 0; i < mVolatility->GetNumLines(); i++)
			for (int j = 0; j < mVolatility->GetNumCols(); j++)
				vol->push_back(mVolatility->Elt(i,j));

	}

	delete mVolatility;

	if (vStrikes)
		delete vStrikes;
	vStrikes = NULL;

	if (vSortStrikes)
		delete vSortStrikes;
	vSortStrikes = NULL;

	return ARM_OK;
}



ARM_VolCube* etoolkit_GetVolCubeFromSummit(const CCString& index,
										   const CCString& currency,
										   const CCString& cvName,
										   ARM_Date date,
										   const CCString& vtype,
										   VECTOR<CCString>& tenors,
										   bool vCorrelCube)
{
    VECTOR<CCString> smileFallBacks;

	ARM_VolCube*  createdVolCube = NULL;

	ARM_VolCurve* ATMVolCrv = NULL;
	ARM_VolCurve* volCv     = NULL;
	ARM_VolCurve* inVols[200];
	double dTenors[200];

	char* ccy = (char*) currency;
	ARM_Currency* sCCY = new ARM_Currency(ccy);
	delete ccy;

	int nbCrv;
	char matIndex[20];


    // Temporary code for treating MO40 FALL BACKS
    // fallbacks should be a an array parameter

    if ( cvName == "MO40" )
    {
       smileFallBacks.push_back("MO");
    }

	ATMVolCrv = etoolkit_GetVolATMFromSummit(index, currency, cvName, date, vtype);

	nbCrv = 0;

	for (int i = 0; i < tenors.size(); i++)
	{
		strcpy(matIndex, tenors[i]);

		try
		{
			volCv = etoolkit_GetSmileFromSummit(index, currency, cvName, date,
                                                vtype, matIndex,
												0 // interpol = LINEAR
												);
		}
		
        catch(...)
		{
			volCv = NULL;

            // Try fall backs

            int size = smileFallBacks.size();
                
            int found = 0;

            for (int h = 0; (( h < size ) && (!(found))); h++)
            {
                try
                {     
                    volCv = etoolkit_GetSmileFromSummit(index, currency,
                                                        smileFallBacks[h],
                                                        date,
                                                        vtype,
                                                        matIndex,
														0 // interpol = LINEAR
														);

                    if ( volCv != NULL )
                    {
                       found = 1;
                    }
                }

                catch(...)
                {
                    volCv = NULL;
                }
            }
		}

		if ( volCv != NULL )
		{
			dTenors[nbCrv] = FromStrMatuToDouble(matIndex);			
			inVols[nbCrv] = volCv;
			nbCrv++;
		}
	}

	if ( nbCrv == 0 )
	{
		if(vCorrelCube)
		{
			if(ATMVolCrv)
			{
				ARM_Currency	vCurrency(currency.c_str());
				volCv = new ARM_VolFlat(date, 0., &vCurrency);
				if(volCv)
				{
					inVols[0] = volCv;
					dTenors[0] = FromStrMatuToDouble(matIndex);
					nbCrv = 1;
				}
				else
				{			
					delete	ATMVolCrv;
					CCString msg((CCString) "No Smile");
					throw	Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
				}
			}
			else
			{
				CCString msg((CCString) "No Correl");
				throw	Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
			}
		}
		else
		{
			delete	ATMVolCrv;
			CCString msg((CCString) "No Smile");
			throw	Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
		}
	}

	ARM_Vector existingTenors(nbCrv,dTenors);

	createdVolCube = new ARM_VolCube(ATMVolCrv,inVols,nbCrv,&existingTenors);

	createdVolCube->SetCurrencyUnit(sCCY);

	for (i = 0; i < nbCrv; i++)
	{
		if (inVols[i])
		   delete inVols[i];
		inVols[i] = NULL;
	}

	if (ATMVolCrv)
		delete ATMVolCrv;
	ATMVolCrv = NULL;

	if (sCCY)
	   delete sCCY;
	sCCY = NULL;

	return(createdVolCube);
}




ARM_VolLInterpol* etoolkit_GetFXVolATMFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist)
{
	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolCube* newVolCub = NULL;

	ARM_Matrix* mVolATM = NULL;
	ARM_Vector* vYearTerms = NULL;

	VECTOR<CCString>* matu = new VECTOR<CCString>;

	ARM_result result;
		
	CCString xmlResponse;
	CCString msgList;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetFXVOL((const char*)ccy1,(const char*)ccy2,(const char*)cvName,date,xmlOutput);

	long retCode = ARM_OK;

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		ARM_XGigaToolKit::doGetFXVOL((const char*)ccy2,(const char*)ccy1,(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getfxvol(ccy1,
											 ccy2,
											 cvName,
											 date,
											 impOrHist,
											 xmlResponse,
											 msgList);


			// on inverse l'ordre des devises
			if ( strcmp((const char*)xmlResponse,"") == 0 )
			{
				retCode = etoolkit_getfxvol(ccy2,
											ccy1,
											cvName,
											date,
											impOrHist,
											xmlResponse,
											msgList);
			}

			if ( strcmp((const char*)xmlResponse,"") == 0 )
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 "Pb in getting ATM FxVol");
 				return(NULL);
			}
		}
	}

	retCode = ARMLOCAL_GetDatesAndDatasForFxVol(xmlResponse,
												date,
												ccy1,
												mVolATM,
												vYearTerms,
												matu);

	ARM_Vector* vStrikes = new ARM_Vector(1,0.0);

	// construction de la matrice de vol forex (atm)
	
	newVolCrv = new ARM_VolLInterpol(date, vYearTerms,
								vStrikes, mVolATM, K_STK_TYPE_PRICE, K_ATMF_VOL);

	for (int i = 0; i < ARM_NB_TERMS; i++) 
	{
		sprintf(newVolCrv->itsYearTermsX[i], "X");
	}

	if ( matu->size() > ARM_NB_TERMS )
	{
	   throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "etoolkit_GetFXVolATMFromSummit: (matu->size() > ARM_NB_TERMS)");
	}

	for (i = 0; i < matu->size(); i++)
	{
		sprintf(newVolCrv->itsYearTermsX[i], "%s", (const char*) (*matu)[i]);
	}

	if (matu)
	   delete matu;
	matu = NULL;

	return newVolCrv;
}


long ARMLOCAL_GetDatesAndDatasForFxCorrel(const char* chaineXML,
										  ARM_Date asOfDate,
										  const char* ccy,
										  ARM_Matrix* &Volatility,
										  ARM_Vector* &vYearTerms)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	ARM_Date settleDate(asOfDate);
	ARM_Currency sCCY((char*)ccy);

	long spotDays = sCCY.GetSpotDays();
	settleDate.NextBusinessDay(spotDays,(char*)ccy);

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	
    catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting FX Correl");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response//Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting FX Correl");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			vYearTerms = new ARM_Vector(nbNodes);
			Volatility = new ARM_Matrix(nbNodes,1);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					char strYT[4];
					strcpy(strYT,ff1);

					long Nb;
					char cMatu;
					long freqId;

					sscanf(strYT, "%d%c", &Nb, &cMatu);

					cMatu = toupper(cMatu);

					if ( cMatu == 'D' ) // Ex : "1D"
						freqId = K_DAILY;
					else if ( cMatu == 'W' )
						freqId = K_WEEKLY;
					else if ( cMatu == 'M' ) 
						freqId = K_MONTHLY;
					else
						freqId = K_ANNUAL;

					ARM_Date tmpDate(settleDate);

					tmpDate.AddPeriodMult(freqId, Nb, (char*)ccy);

					vYearTerms->Elt(indexNode) = (tmpDate - asOfDate) / 365.;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}

				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					Volatility->Elt(indexNode,0) = atof(ff1) * 100.;

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return ARM_OK;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting FX Correl");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	return ARM_KO;
}



ARM_VolLInterpol* etoolkit_GetFXCorrelFromSummit(const CCString& ccy1,
												 const CCString& index,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const VECTOR<CCString>& tenors)
{

	ARM_VolLInterpol* newVolCrv = NULL;

	CCString xmlResponse;
	CCString msgList;

	ARM_Vector** vYearTerms = new ARM_Vector* [tenors.size()];
	ARM_Matrix** mVols = new ARM_Matrix* [tenors.size()];

	for (int i = 0; i < tenors.size(); i++)
	{
		vYearTerms[i] = NULL;
		mVols[i] = NULL;
	}

	long retCode;
	CCString msg ("");

	long nbTenors = 0;
	VECTOR<CCString> realTenors;

	for (i = 0; i < tenors.size(); i++)
	{
		retCode = etoolkit_getfxcorrel(ccy1,
									   index,
									   tenors[i],
									   ccy2,
									   cvName,
									   date,
									   xmlResponse,
									   msgList);

		if (retCode == ARM_OK)
		{
			realTenors.push_back(tenors[i]);

			retCode = ARMLOCAL_GetDatesAndDatasForFxCorrel(xmlResponse,
														   date,
														   ccy2,
														   mVols[nbTenors],
														   vYearTerms[nbTenors]);

			if (nbTenors != 0)
			{
				if (vYearTerms[nbTenors]->GetSize() != vYearTerms[nbTenors-1]->GetSize())
				{
					for (int j = 0; j < tenors.size(); j++)
					{
						if (vYearTerms[j])
							delete vYearTerms[j];
						vYearTerms[j] = NULL;

						if (mVols[j])
							delete mVols[j];
						mVols[j] = NULL;
					}

					delete [] vYearTerms;
					delete [] mVols;

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "not the same size of maturities between 2 tenors");
				}
				else
				{
					for (int j = 0; j < vYearTerms[nbTenors]->GetSize(); j++)
					{
						if (vYearTerms[nbTenors]->Elt(j) != vYearTerms[nbTenors-1]->Elt(j))
						{
							for (int k = 0; k < tenors.size(); k++)
							{
								if (vYearTerms[k])
									delete vYearTerms[k];
								vYearTerms[k] = NULL;

								if (mVols[k])
									delete mVols[k];
								mVols[k] = NULL;
							}

							delete [] vYearTerms;
							delete [] mVols;

							throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "not the same maturities between 2 tenors");
						}
					}
				}
			}

			nbTenors ++;
		}
	}

	if (nbTenors == 0)
	{
		for (int k = 0; k < tenors.size(); k++)
		{
			if (vYearTerms[k])
				delete vYearTerms[k];
			vYearTerms[k] = NULL;

			if (mVols[k])
				delete mVols[k];
			mVols[k] = NULL;
		}

		delete [] vYearTerms;
		delete [] mVols;

		CCString msg((CCString)"No data corresponding to parameters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	ARM_Vector* vRealYearTerms = new ARM_Vector (vYearTerms[0]->GetSize());
	ARM_Matrix* mRealVols = new ARM_Matrix (vYearTerms[0]->GetSize(),nbTenors);

	ARM_Vector* vStrikes = new ARM_Vector(nbTenors);
	char sTerm[5];

	for (i = 0; i < nbTenors; i++)
	{
		strcpy(sTerm,realTenors[i]);
		vStrikes->Elt(i) = GenMatuToYearTerm(StringToGenMatu(sTerm));
	}

	for (i = 0; i < vRealYearTerms->GetSize(); i++)
	{
		vRealYearTerms->Elt(i) = vYearTerms[0]->Elt(i);
		for (int j = 0; j < nbTenors; j++)
		{
			mRealVols->Elt(i,j) = mVols[j]->Elt(i,0);
		}
	}

	newVolCrv = new ARM_VolLInterpol(date, vRealYearTerms,
								vStrikes, mRealVols, K_STK_TYPE_YIELD, K_ATMF_VOL);

	newVolCrv->SetIndexName((char*)(const char*)index);

	for (i = 0; i < tenors.size(); i++)
	{
		if (vYearTerms[i])
			delete vYearTerms[i];
		vYearTerms[i] = NULL;

		if (mVols[i])
			delete mVols[i];
		mVols[i] = NULL;
	}

	delete [] vYearTerms;
	delete [] mVols;

	return newVolCrv;
}



ARM_VolLInterpol* etoolkit_GetCorrelFromSummit(const CCString& ccy1,
											   const CCString& index1,
											   const CCString& ccy2,
											   const CCString& index2,
											   ARM_Date date,
											   const CCString& cvName)
{
	ARM_VolLInterpol* newVolCrv = NULL;

	CCString xmlResponse;
	CCString msgList;

	CCString tmpCcy1 = ccy1;
	CCString tmpCcy2 = ccy2;
	CCString tmpIndex2 = index2;
	int isTranspose = 0;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetCORRELTenorsList((const char*)index1,(const char*)ccy1,(const char*)index2,(const char*)ccy2,(const char*)cvName,date,xmlOutput);

	long retCode = ARM_OK;

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		ARM_XGigaToolKit::doGetCORRELTenorsList((const char*)index2,(const char*)ccy2,(const char*)index1,(const char*)ccy1,(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		// Modif des var tmp
			tmpCcy1 = ccy2;
			tmpCcy2 = ccy1;
			tmpIndex2 = index1;
			isTranspose=1;
		}
		else
		{
			retCode = etoolkit_getcorrel(ccy1,
										  index1,
										  ccy2,
										  index2,
										  cvName,
										  date,
										  xmlResponse,
										  msgList);

			if (retCode == ARM_KO)
			{
				return NULL;
			}

			if (strcmp((const char*)xmlResponse,"<Response></Response>") == 0)
			{
				tmpCcy1 = ccy2;
				tmpCcy2 = ccy1;
				tmpIndex2 = index1;

				retCode = etoolkit_getcorrel(ccy2,
											 index2,
											 ccy1,
											 index1,
											 cvName,
											 date,
											 xmlResponse,
											 msgList);

				if (retCode == ARM_KO)
				{
					return NULL;
				}

				isTranspose = 1;
			}

			if (strcmp((const char*)xmlResponse,"<Response></Response>") == 0)
			{
				CCString msg((CCString)"Pb for getting Correl");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char*) msg);
			}
		}
	}

	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listTenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,
                                                                &listRequetes,
                                                                "P5");

	ARM_Vector* vStrikes = new ARM_Vector(listTenors.size());

	int i;

	ARM_Currency cCcy1((const char*)tmpCcy1);
	long spotDays = cCcy1.GetSpotDays();
    char* payCalTmp = cCcy1.GetPayCalName(cCcy1.GetVanillaIndexType());
    ARM_Date settleDate(date);
	settleDate.NextBusinessDay(spotDays, payCalTmp);

    char payCal[30];
    strcpy(payCal, payCalTmp);
    delete [] payCalTmp;

	for (i = 0; i < listTenors.size(); i++)
	{
		vStrikes->Elt(i) = convPlotInYearTerm((const char *) listTenors[i], date, 
                                                settleDate,
                                                payCal);

//        vStrikes->Elt(i) = FromStrMatuToDouble((const char *) listTenors[i], &date);
    }

	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortRequetes;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());

	for (i = 0; i < vStrikes->GetSize(); i++)
	{
		int j=0;
		while ( j < vStrikes->GetSize() )
		{
			if ( vSortStrikes->Elt(i) != vStrikes->Elt(j) )
			{
				j++;
			}
			else
			{
				sortRes.push_back(listTenors[j]);
				sortRequetes.push_back(listRequetes[j]);
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	if(!isTranspose) 
		ARM_XGigaToolKit::doGetCORRELTenor((const char*)index1,(const char*)ccy1,(const char*)index2,(const char*)ccy2,(const char*)sortRes[0],(const char*)cvName,date,xmlOutput);
	else
		ARM_XGigaToolKit::doGetCORRELTenor((const char*)index2,(const char*)ccy2,(const char*)index1,(const char*)ccy1,(const char*)sortRes[0],(const char*)cvName,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		retCode = etoolkit_getCorrelCurveByMatu(sortRequetes[0],
												tmpCcy2,
												tmpIndex2,
												cvName,
												date,
												xmlResponse,
												msgList);
	}

	VECTOR<CCString> listYearTerms;
	
	retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
												tmpCcy2,
												date,
												0,
												mVolatility,
												vYearTerms,
												sortRes.size(),
												1,
												&listYearTerms);

	for (i = 1; i < sortRes.size(); i++)
	{
		if(!isTranspose) 
			ARM_XGigaToolKit::doGetCORRELTenor((const char*)index1,(const char*)ccy1,(const char*)index2,(const char*)ccy2,(const char*)sortRes[i],(const char*)cvName,date,xmlOutput);
		else
			ARM_XGigaToolKit::doGetCORRELTenor((const char*)index2,(const char*)ccy2,(const char*)index1,(const char*)ccy1,(const char*)sortRes[i],(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getCorrelCurveByMatu(sortRequetes[i],
													tmpCcy2,
													tmpIndex2,
													cvName,
													date,
													xmlResponse,
													msgList);
		}
		retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
													tmpCcy2,
													date,
													i,
													mVolatility,
													vYearTerms,
													0,
													1);
	}



	if (isTranspose == 1)
	{
		mVolatility->Transpose();

		newVolCrv = new ARM_VolLInterpol(date,vSortStrikes,vYearTerms,
										 mVolatility, 1, K_ATMF_VOL);
	}
	else
	{
		newVolCrv = new ARM_VolLInterpol(date,vYearTerms,
										 vSortStrikes, mVolatility, 1, K_ATMF_VOL);

	}

	for (i = 0; i < ARM_NB_TERMS; i++) 
	{
		strcpy(newVolCrv->itsYearTermsX[i], "X");
		strcpy(newVolCrv->itsYearTermsY[i], "X");
	}
	
	if(ccy1==ccy2)
	{
		char* tmpccy = (char*)ccy1;
		ARM_Currency* myCcy = new ARM_Currency(tmpccy);
		delete [] tmpccy;

		if(myCcy != ARM_DEFAULT_CURRENCY)
			newVolCrv->SetCurrencyUnit(myCcy);

		if (myCcy)
			delete myCcy;
		myCcy = NULL;
	}

	// Correl case
	// we need to fill itsYearTerms
	for (i = 0; i < newVolCrv->GetExpiryTerms()->GetSize(); i++)
	{
		sprintf(newVolCrv->itsYearTermsX[i], "%s", (const char*) listYearTerms[i]);
	}

	for (i = 0; i < newVolCrv->GetStrikes()->GetSize(); i++)
	{
		sprintf(newVolCrv->itsYearTermsY[i], "%s", (const char*) sortRes[i]);
	}

	string	vIndex((const char *) ccy2);
	string	vCurrency((const char *) ccy1);
	string	vCrvId((const char *) cvName);
	string	vType("IRIR CORR");
	string	vNbFact((const char *) index1);
	if(vNbFact.substr(0, 3) == "COR")
	{
		vNbFact = vNbFact.substr(3, 2);
		vType += " " + vNbFact;
	}

	newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);


	if (vStrikes)
		delete vStrikes;
	vStrikes = NULL;

	return newVolCrv;
}


double etoolkit_getXMLMeanRevFromSummit(const CCString& C_ccy,
										const CCString& C_index,
										const CCString& C_cvname,
										const ARM_Date& Date,
										const CCString& C_NumFactor)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getMeanRev(C_ccy,
									 C_index,
									 C_cvname,
									 C_NumFactor,
									 Date,
									 xmlResponse,
									 messageList);

		if (retCode == ARM_OK)
		{
			double res = ARMLOCAL_ParseXMLForMeanRev(xmlResponse);

			return res;
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Pb for getting Mean Reversion" );

	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Pb for getting Mean Reversion" );
	}
}


double etoolkit_getXMLCutOffFromSummit(const CCString& C_ccy,
									   const CCString& C_index,
									   const CCString& C_cvname,
									   const CCString& C_NumFactor,
									   const ARM_Date& Date)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getCutOff(C_ccy,
									 C_index,
									 C_cvname,
									 C_NumFactor,
									 Date,
									 xmlResponse,
									 messageList);

		double res = ARMLOCAL_ParseXMLForCutOff(xmlResponse);

		return res;
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Pb for getting CutOff" );
	}
}


double ARMLOCAL_ParseXMLForMeanRev(const CCString& xmlResponse)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double res;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for getting Mean Reversion");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for getting Mean Reversion");
		
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				res = atof(ff1) * 100.;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->Release();
			listItem=NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for getting Mean Reversion");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char*) msg);
	}
}


double ARMLOCAL_ParseXMLForCutOff(const CCString& xmlResponse)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double res;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg ((CCString)"Pb in creating XML document for getting CutOff");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;

				CCString msg ((CCString)"Invalid XML string for getting CutOff");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				res = atof(ff1) * 100.;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->Release();
			listItem=NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg ((CCString)"Error in XML parsing for getting CutOff");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



ARM_Vector* etoolkit_getXMLQmodParamFromSummit(const CCString& C_ccy,
										  const CCString& C_index,
										  const CCString& C_cvname,
										  const ARM_Date& Date)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getQmodParam(C_ccy,
										C_index,
										C_cvname,
										Date,
										xmlResponse,
										messageList);

		ARM_Vector* res = (ARM_Vector*) ARMLOCAL_ParseXMLForQmodParam(xmlResponse);

		return res;
	}
	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		CCString msg((CCString)"Pb for getting Q model paramaters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



ARM_Vector* ARMLOCAL_ParseXMLForQmodParam(const CCString& xmlResponse)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Vector* res;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		// JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Q model paramaters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			res= new ARM_Vector(nbNodes);
			if (nbNodes != 2)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Q model paramaters");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{	
				hr=resultList->get_item(indexNode, &listItem);
				
				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					res->Elt(indexNode) = 100.0*atof(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}	

			
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Q model paramaters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


ARM_VolLInterpol* etoolkit_getXMLQFXParamFromSummit(const CCString& C_ccy,
												    const CCString& C_index,
												    const CCString& C_cvname,
												    const ARM_Date& Date)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getQFXParam(C_ccy,
										C_index,
										C_cvname,
										Date,
										xmlResponse,
										messageList);

		ARM_VolLInterpol* res = (ARM_VolLInterpol*) ARMLOCAL_ParseXMLForQFXParam(xmlResponse, Date);
        for (int i = 0; i < res->GetExpiryTerms()->GetSize(); i++)
		{
			sprintf(res->itsYearTermsX[i], "%s", (const char*) YearTermToStringMatu(res->GetExpiryTerms()->Elt(i)).c_str());
			// convert expiries from year terms to time lags (365 days per year ??)
			res->GetExpiryTerms()->Elt(i) *= 365.0;
		}

		return res;
	}
	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		CCString msg((CCString)"Pb for getting Q FX model paramaters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}

ARM_VolLInterpol* ARMLOCAL_ParseXMLForQFXParam(const CCString& xmlResponse, const ARM_Date& Date)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Vector* terms;
	ARM_Vector* values;
	ARM_VolLInterpol* res;

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;
		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		CCString msg((CCString)"Pb in creating XML document for getting QFX paramaters");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,(const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No QFX data");

			terms = new ARM_Vector(nbNodes);
			values = new ARM_Vector(nbNodes);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{	
				hr=resultList->get_item(indexNode, &listItem);
				
				listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					terms->Elt(indexNode) = StringMatuToYearTerm(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					values->Elt(indexNode) = 100.0*atof(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}	
		}

		res = new ARM_VolLInterpol(Date, terms, values);

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting QFX paramaters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}



void ARMLOCAL_GetListFxStrikesFromXML(const char* chaineXML, VECTOR<CCString>* listRequetes)
{

	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double* pdRate = NULL;
	double* pdMatu = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Strikes for Forex Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Commdata/Name2")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Strikes for Forex Smile");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{
				hr=resultList->get_item(indexNode, &listItem);

				listItem->selectSingleNode(_bstr_t((const char *)"Val"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					listRequetes->push_back(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}
		
		hr = S_OK;
	}

	catch(...)
	{		
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (theNode) theNode->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Strikes for Forex Smile");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}


}


ARM_VolLInterpol* etoolkit_GetXMLFxSmileFromSummit(const CCString& ccy1,
												  const CCString& ccy2,
												  const CCString& cvName,
												  ARM_Date date,
												  const char* CallOrPut,
												  int IncOrDec)
{
	ARM_VolLInterpol *newVolCrv = NULL;

	ARM_result C_result;

	CCString xmlResponse;
	CCString msgList;
	CCString request;

	CCString tmpCcy1(ccy1);
	CCString tmpCcy2(ccy2);

	long retCode;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy1,(const char*)ccy2,CallOrPut,(const char*)cvName,date,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy2,(const char*)ccy1,CallOrPut,(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			tmpCcy1 = ccy2;
			tmpCcy2 = ccy1;
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			request = (CCString)"SMILE/"  + ccy1 + (CCString)"/" + ccy2 + (CCString)"/FXOPT/" + CallOrPut + "/";
			retCode = etoolkit_getFxSmile(request,
												cvName,
												date,
												xmlResponse,
												msgList);

			
			// on inverse l'ordre des devises
			if (strcmp((const char*)xmlResponse,"<Response></Response>") == 0)
			{
				request=(CCString) "SMILE/" + ccy2 + "/" + ccy1 + "/FXOPT/C/";
				retCode = etoolkit_getFxSmile(request,
											cvName,
											date,
											xmlResponse,
											msgList);
			}
		}
	}

	if ( (retCode == ARM_KO) || (strcmp(xmlResponse,"<Response></Response>") == 0) )
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Error in getting forex smile");
	}

	VECTOR<CCString> listStrikes;

	try
	{
		ARMLOCAL_GetListFxStrikesFromXML(xmlResponse,&listStrikes);
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Error in getting forex smile");
//		return NULL;
	}

	ARM_Vector* vStrikes = new ARM_Vector(listStrikes.size());
	
	int i;

	for (i=0;i<vStrikes->GetSize();i++)
		vStrikes->Elt(i) = atof(listStrikes[i]) / 100.;

	// on trie les strikes
	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortGSListStrikes;

	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort(IncOrDec));

	
	for (i=0;i<vStrikes->GetSize();i++)
	{
		int j=0;
		while (j<vStrikes->GetSize())
		{
			if (vSortStrikes->Elt(i) != vStrikes->Elt(j))
			{
				j++;
			}
			else
			{
				sortRes.push_back(listStrikes[j]);

				// Pour Gigaspaces
				double strikeGiga = atof(listStrikes[j]);
				char sGiga[128];
				sprintf(sGiga, "%.0lf", strikeGiga);
				sortGSListStrikes.push_back(sGiga);

				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolatility = NULL;

	try
	{
		for (i=0;i<sortRes.size();i++)
		{
			ARM_XGigaToolKit::doGetSmileFXTenor((const char*)tmpCcy1,(const char*)tmpCcy2,CallOrPut,(const char*)sortGSListStrikes[i],(const char*)cvName,date,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				retCode = etoolkit_getFxSmileByStrike(sortRes[i],
														 request,
														 cvName,
														 date,
														 "SMILE",
														 xmlResponse,
														 msgList);
			}

			retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
														ccy1,
														date,
														i,
														mVolatility,
														vYearTerms,
														sortRes.size());
		}
	}
	catch (Exception& e)
	{
		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		if (vSortStrikes)
			delete vSortStrikes;
		vSortStrikes = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		if (mVolatility)
			delete mVolatility;
		mVolatility = NULL;

		throw e;
	}
	catch (...)
	{		
		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		if (vSortStrikes)
			delete vSortStrikes;
		vSortStrikes = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		if (mVolatility)
			delete mVolatility;
		mVolatility = NULL;

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Error in getting forex smile");
	}

	newVolCrv = new ARM_VolLInterpol( date, vYearTerms,
								vSortStrikes, mVolatility, 1, K_SMILE_VOL);

	if (vStrikes)
		delete vStrikes;
	vStrikes = NULL;

	return newVolCrv;
}



long etoolkit_GetInitialFXVolFromSummit (const CCString& ccy1,
										 const CCString& ccy2,
										 const CCString& cvName,
										 const ARM_Date& date,
										 const CCString& impOrHist,
										 VECTOR<CCString>* vMaturities,
										 ARM_Vector* vVolATM)
{
	CCString xmlResponse;
	CCString msgList;
	ARM_Matrix* mVolATM = NULL;
	ARM_Vector* vMatu = NULL;

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetFXVOL((const char*)ccy1,(const char*)ccy2,(const char*)cvName,date,xmlOutput);

	long retCode = ARM_OK;

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		ARM_XGigaToolKit::doGetFXVOL((const char*)ccy2,(const char*)ccy1,(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			long retCode = etoolkit_getfxvol(ccy1,
											 ccy2,
											 cvName,
											 date,
											 impOrHist,
											 xmlResponse,
											 msgList);


			// on inverse l'ordre des devises
			if (strcmp((const char*)xmlResponse,"") == 0)
			{
				retCode = etoolkit_getfxvol(ccy2,
											ccy1,
											cvName,
											date,
											impOrHist,
											xmlResponse,
											msgList);
			}

			if (strcmp((const char*)xmlResponse,"") == 0)
			{
				return ARM_KO;
			}
		}
	}

	retCode = ARMLOCAL_GetDatesAndDatasForFxVol(xmlResponse,
												date,
												ccy1,
												mVolATM,
												vMatu,
												vMaturities);

	if (retCode == ARM_OK)
	{
		vVolATM->Resize(vMatu->GetSize());

		for (int i = 0; i < vMatu->GetSize(); i++)
		{
			vVolATM->Elt(i) = mVolATM->Elt(i,0);
		}
	}

	if (vMatu)
		delete vMatu;
	vMatu = NULL;

	return retCode;
}




long etoolkit_GetInitialFXSmileFromSummit (const CCString& ccy1,
										   const CCString& ccy2,
										   const CCString& cvName,
										   const ARM_Date& date,
										   const CCString& impOrHist,
										   ARM_Matrix* mSmile,
										   ARM_Vector* vTenors)
{
	CCString xmlResponse("");
	CCString xmlResponseCall("");
	CCString xmlResponsePut("");
	CCString msgList;
	CCString request;
	CCString tmpCcy1 = ccy1;
	CCString tmpCcy2 = ccy2;

	long retCode;

	// ** GIGASPACE
	std::string xmlOutput;
	std::string xmlOutputCall;
	std::string xmlOutputPut;

	ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy1,(const char*)ccy2,"C",(const char*)cvName,date,xmlOutputCall);
	ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy1,(const char*)ccy2,"P",(const char*)cvName,date,xmlOutputPut);

	if ( (xmlOutputCall != "") || (xmlOutputPut != "") )
	{
		xmlResponseCall = (CCString)(xmlOutputCall.c_str());
		xmlResponsePut = (CCString)(xmlOutputPut.c_str());
	}
	else
	{
		ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy2,(const char*)ccy1,"C",(const char*)cvName,date,xmlOutputCall);
		ARM_XGigaToolKit::doGetSmileFXTenorList((const char*)ccy2,(const char*)ccy1,"P",(const char*)cvName,date,xmlOutputPut);
	
		if ( (xmlOutputCall != "") || (xmlOutputPut != "") )
		{
			xmlResponseCall = (CCString)(xmlOutputCall.c_str());
			xmlResponsePut = (CCString)(xmlOutputPut.c_str());
			tmpCcy1 = ccy2;
			tmpCcy2 = ccy1;
		}
		else
		{
			request = (CCString)"SMILE/"  + ccy1 + (CCString)"/" + ccy2 + (CCString)"/FXOPT/%/";

			retCode = etoolkit_getFxSmile(request,
											   cvName,
											   date,
											   xmlResponse,
											   msgList);

			// on inverse l'ordre des devises
			if (strcmp((const char*)xmlResponse,"<Response></Response>") == 0)
			{
				request=(CCString) "SMILE/" + ccy2 + "/" + ccy1 + "/FXOPT/%/";
			
				retCode = etoolkit_getFxSmile(request,
											  cvName,
											  date,
											  xmlResponse,
											  msgList);
			}

			if ( (retCode == ARM_KO) || (strcmp(xmlResponse,"<Response></Response>") == 0) )
			{
				return ARM_KO;
			}
		}
	}

	VECTOR<CCString> listStrikes;
	VECTOR<CCString> listStrikesETK;
	VECTOR<CCString> listStrikesCall;
	VECTOR<CCString> listStrikesPut;

	try
	{
		if (strcmp((const char*)xmlResponse,"") != 0)
		{
			ARMLOCAL_GetListFxStrikesFromXML(xmlResponse,&listStrikesETK);
			listStrikes = listStrikesETK;
		}
		else
		{
			// les Strikes PUT sont négatifs, les strikes CALL sont positifs
			ARMLOCAL_GetListFxStrikesFromXML(xmlResponsePut,&listStrikesPut);
			for (int i = 0; i < listStrikesPut.size(); i++)
				listStrikes.push_back(listStrikesPut[i]);

			ARMLOCAL_GetListFxStrikesFromXML(xmlResponseCall,&listStrikesCall);
			for (i = 0; i < listStrikesCall.size(); i++)
				listStrikes.push_back(listStrikesCall[i]);
		}
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Error in getting forex smile");
	}

	ARM_Vector* vStrikes = new ARM_Vector(listStrikes.size());

	int i;

	for (i=0;i<vStrikes->GetSize();i++)
		vStrikes->Elt(i) = atof(listStrikes[i]) / 100.;

	// on trie les strikes
	VECTOR<CCString> sortRes;
	VECTOR<CCString> sortGSListStrikes;
	
	vTenors->Resize(vStrikes->GetSize());
	ARM_Vector* vSortStrikes = new ARM_Vector(vStrikes->Sort());
	for (i = 0; i < vStrikes->GetSize(); i++)
		vTenors->Elt(i) = vSortStrikes->Elt(i);

	if (vSortStrikes)
		delete vSortStrikes;
	vSortStrikes = NULL;

	for (i=0;i<vStrikes->GetSize();i++)
	{
		int j=0;
		while (j<vStrikes->GetSize())
		{
			if (vTenors->Elt(i) != vStrikes->Elt(j))
			{
				j++;
			}
			else
			{
				sortRes.push_back(listStrikes[j]);

				// Pour Gigaspaces
				double strikeGiga = atof(listStrikes[j]);
				char sGiga[128];
				sprintf(sGiga, "%.0lf", strikeGiga);
				sortGSListStrikes.push_back(sGiga);
				break;
			}
		}
	}

	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mTmpSmiles = NULL;

	try
	{
		
		for (i=0;i<sortRes.size();i++)
		{
			if ( vStrikes->Elt(i) > 0 )
				ARM_XGigaToolKit::doGetSmileFXTenor((const char*)tmpCcy1,(const char*)tmpCcy2,"C",(const char*)sortGSListStrikes[i],(const char*)cvName,date,xmlOutput);
			else
				ARM_XGigaToolKit::doGetSmileFXTenor((const char*)tmpCcy1,(const char*)tmpCcy2,"P",(const char*)sortGSListStrikes[i],(const char*)cvName,date,xmlOutput);

			if (xmlOutput != "")
			{
				xmlResponse = (CCString)(xmlOutput.c_str());
			}
			else
			{
				retCode = etoolkit_getFxSmileByStrike(sortRes[i],
													  request,
													  cvName,
													  date,
													  "SMILE",
													  xmlResponse,
													  msgList);
			}

			retCode = ARMLOCAL_GetDatesAndDatasForCurve(xmlResponse,
														ccy1,
														date,
														i,
														mTmpSmiles,
														vYearTerms,
														sortRes.size());
		}

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		mSmile->Resize(mTmpSmiles->GetNumLines(),mTmpSmiles->GetNumCols());
		
		for (i = 0; i < mTmpSmiles->GetNumLines(); i++)
			for (int j = 0; j < mTmpSmiles->GetNumCols(); j++)
				mSmile->Elt(i,j) = mTmpSmiles->Elt(i,j);
	
		if (mTmpSmiles)
			delete mTmpSmiles;
		mTmpSmiles = NULL;

		return ARM_OK;
	}

	catch (...)
	{
		if (vStrikes)
			delete vStrikes;
		vStrikes = NULL;

		if (vYearTerms)
			delete vYearTerms;
		vYearTerms = NULL;

		return ARM_KO;
	}
}



ARM_VolLInterpol* etoolkit_CreateFXVolFromSummit(const CCString& ccy1,
												 const CCString& ccy2,
												 ARM_Date date,
												 const CCString& cvName,
												 const CCString& impOrHist,
												 const CCString& VolType)
{
	ARM_VolLInterpol* newVolCrv = NULL;
	ARM_VolLInterpol* FXVolATM = NULL;
	ARM_VolCube* newVolCub = NULL;

	/// voltype is only used when ETK is used
	/// hence included here
	int volType;
	ARM_Vector* vStrikes=NULL;
	ARM_Vector* vYearTerms = NULL;
	ARM_Matrix* mVolATM = NULL;
	ARM_Matrix* mVols=NULL;
	ARM_VolLInterpol* smile = NULL;


	FXVolATM = etoolkit_GetFXVolATMFromSummit(ccy1,ccy2,date,cvName,impOrHist);

	if (FXVolATM == NULL)
		return NULL;

	if (strcmp(VolType,"ATM")==0)	
	{
		smile=NULL;
		volType=K_ATMF_VOL;
	}
	else if ( (strcmp(VolType,"SMILE")==0)
			|| (strcmp(VolType,"FXSPI")==0)
			|| (strcmp(VolType,"CSMILE")==0)
			|| (strcmp(VolType,"CFXSPI")==0) )
	{
		smile = etoolkit_GetXMLFxSmileFromSummit(ccy1,ccy2,cvName,date);
		volType = K_FX_VOL_SP_INTERP;
	}
	else 
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Invalid VolType Argument for getting FX Vol Smile");

	// on recupere les strikes de l'objet smile;
	if ( (smile) && (volType==K_FX_VOL_SP_INTERP) && ((strcmp(VolType,"SMILE")==0) || (strcmp(VolType,"FXSPI")==0)) )
	{
		vStrikes= new ARM_Vector(smile->GetStrikes()->GetSize()+1);
		vStrikes->Elt(0)=0.0;
		for (int i=1; i<vStrikes->GetSize();i++)
		{
			vStrikes->Elt(i)=((ARM_Vector*) smile->GetStrikes())->Elt(i-1);
		}
		mVols= MergeMatrix(FXVolATM->GetVolatilities(), smile->GetVolatilities());
		
	}
	else if ( !(smile) && (volType==K_FX_VOL_SP_INTERP) )
	{
		if (FXVolATM)
			delete FXVolATM;
		FXVolATM = NULL;

		return NULL;
	}
	else
	{
		vStrikes = new ARM_Vector(1,0.0);
		mVols= (ARM_Matrix*) FXVolATM->GetVolatilities()->Clone();
	}

	if ((strcmp(VolType,"CSMILE")==0) || (strcmp(VolType,"CFXSPI")==0))
	{
		ARM_Vector* noStrikes = new ARM_Vector(1,0.0);
		
		ARM_VolCurve* inVols[1];
		inVols[0]=smile;
		
		ARM_Vector* underlying = new ARM_Vector(1,0.0);

		newVolCub = new ARM_VolCube((ARM_VolCurve*) FXVolATM,inVols,1,underlying,volType);

		delete FXVolATM;
		delete underlying;

		newVolCrv = newVolCub;
	}
	else
	{
		newVolCrv = new ARM_VolLInterpol(date, (ARM_Vector*) FXVolATM->GetExpiryTerms()->Clone(),
									vStrikes, mVols, K_STK_TYPE_PRICE, volType);
		for (int i=0; i<ARM_NB_TERMS; i++) 
		{
			strcpy(newVolCrv->itsYearTermsX[i], FXVolATM->itsYearTermsX[i]);
		}

		delete FXVolATM;
	}
	if (smile)
		delete smile;
	smile = NULL;

	// Update Mkt data characteristics    
    string	vType = string("FXVOL SMILE IMPROVED");
    string	vIndex((const char*) ccy2);
	string	vCurrency((const char*) ccy1);
	string	vCrvId((const char*) cvName);

   newVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	return newVolCrv;
}




ARM_CapFloor* ARMLOCAL_ParseCap(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
 
 
	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Cap from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	CCString structureId("");
 
 	string sBookName, sStructureId,sCustId,sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_CapFloor* sec = NULL; 
	VECTOR<string> listAssetId;

 	sec = (ARM_CapFloor*) ParserManager::BuildInstrument("IRG",
														 XMLDoc, date, "",
														 sBookName, sStructureId,
														 sCustId, sDealId, listAssetId);
	bookName = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId = sCustId.c_str();
 	dealId = sDealId.c_str();
 	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(sec);
}




double etoolkit_getXMLFixingFromSummit(const CCString& C_source,
									   const CCString& C_index,
									   const CCString& C_tenor,
									   const CCString& C_ccy,
									   const ARM_Date& Date)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getFixing(C_source,
									 C_index,
									 C_tenor,
									 C_ccy,
									 Date,
									 xmlResponse,
									 messageList);

		double res = ARMLOCAL_ParseXMLForFixing(xmlResponse);

		return res;
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch (...)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "Pb for getting Fixing" );
	}
}


double ARMLOCAL_ParseXMLForFixing(const CCString& xmlResponse)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	double res;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Pb in creating XML document for getting Fixing");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/REFRATE")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes != 1)
			{
				hr = S_FALSE;
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Invalid XML string for getting Fixing");
			}

			hr=resultList->get_item(0, &listItem);

			listItem->selectSingleNode(_bstr_t((const char *)"Rate"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);

				_bstr_t ff(resultat,false);
				char * ff1=(char *)ff;

				res = atof(ff1) * 100.;

				theNode->Release();
				theNode=NULL;
				if (resultat) SysFreeString(resultat);
			}

			listItem->Release();
			listItem=NULL;
		}

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
         "Error in XML parsing for getting Fx Volatility");
	}
}


ARM_SpreadOption* ARMLOCAL_ParseSpreadOption(const char* chaineXML, 
                                             const ARM_Date& date, CCString& bookName, CCString& custId, 
											 CCString& dealId, ARM_ReferenceValue& Premium, long isEtk,
											 int aResetFreqForCorridorOptim)
{
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;


	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

        wchar_t* wcharStr = constchar2wchar(chaineXML);

		_bstr_t tmpChaine(wcharStr); 

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) 
           XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	CCString structureId("");

	string sBookName, sStructureId,sCustId,sDealId;

	ParserManager::Init();

	VECTOR<string> listAssetId;
	ARM_SpreadOption* sec = (ARM_SpreadOption*) ParserManager::BuildInstrument("SPDOPT", XMLDoc, date, "",
													 sBookName, sStructureId,
													 sCustId, sDealId, listAssetId, 
													 aResetFreqForCorridorOptim);

	bookName = sBookName.c_str();
	structureId = sStructureId.c_str();
	custId = sCustId.c_str();
	dealId = sDealId.c_str();
	
	if (XMLDoc)
		XMLDoc->Release();


	return(sec);
}


ARM_MaturityCapCalculator* ARMLOCAL_ParseMaturityCap(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk)
{
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;


	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);

		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	CCString structureId("");

	string sBookName, sStructureId,sCustId,sDealId;

	ParserManager::Init(isEtk);

	ARM_MaturityCapCalculator* sec = NULL;

	VECTOR<string> listAssetId;
	sec = (ARM_MaturityCapCalculator*) ParserManager::BuildInstrument("MATURITYCAP", 
																	  XMLDoc, date, "",
																	  sBookName, sStructureId,
																	  sCustId, sDealId, listAssetId);

	bookName = sBookName.c_str();
	structureId = sStructureId.c_str();
	custId = sCustId.c_str();
	dealId = sDealId.c_str();
	
	if (XMLDoc)
		XMLDoc->Release();

	return(sec);
}



ARM_TARNCalculator* ARMLOCAL_ParseTarn(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
 
 
	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	CCString structureId("");
 
 	string sBookName, sStructureId,sCustId,sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_TARNCalculator* sec = NULL; 
	VECTOR<string> listAssetId;

 	sec = (ARM_TARNCalculator*) ParserManager::BuildInstrument(	"RFTARN", 
 																XMLDoc, date, "",
 																sBookName, sStructureId,
 																sCustId, sDealId, listAssetId);
	bookName = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId = sCustId.c_str();
 	dealId = sDealId.c_str();
 	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(sec);
}


ARM_BermudaSwaptionCalculator* ARMLOCAL_ParseBermudaSwaption(const char* chaineXML, const ARM_Date& asOfDate, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk )
{
	try
	{
		ARM_Swaption* swaption = ARMLOCAL_ParseSwaption(chaineXML, bookName, CCString(""), custId, dealId, Premium );

		int fixFreqIfZC = K_DEF_FREQ;
		// si la fixleg est ZC, on récupère en plus la frequence COMP_CompFreq
		if (swaption->GetFixedLeg()->GetPaymentFreq() == K_ZEROCOUPON)
		{
			VARIANT_BOOL bOK;
			HRESULT hr;
			MSXML2::IXMLDOMDocument *XMLDoc = NULL;
			MSXML2::IXMLDOMNodeList *resultList = NULL;
			MSXML2::IXMLDOMNode *listItem = NULL, *theNode = NULL, *item = NULL;
			long nbNodes;

			try
			{
				hr = CoInitialize(NULL); 
				hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
				SUCCEEDED(hr) ? 0 : throw hr;
				_bstr_t tmpChaine = chaineXML;
				XMLDoc->loadXML(tmpChaine, &bOK);
			}
			catch(...)
			{
				if (XMLDoc) XMLDoc->Release();
				hr = S_FALSE;
				CCString msg((CCString)"Pb in creating XML document for getting Swaption");
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
			}

			CCString tmpChaine;
			if (isEtk == 1)
				tmpChaine = (CCString)ETKSTRING + (CCString)"SWAPTION/Assets/ASSET";
		
			if (XMLDoc->selectNodes(_bstr_t((const char *)(tmpChaine)), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes != 2)
				{
					hr = S_FALSE;
					CCString msg((CCString)"Invalid XML string for getting Swaption");
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
				}

				for (int i = 0 ; i < nbNodes ; i++)
				{
					hr=resultList->get_item(i, &item);

					item->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
					if (theNode != NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;
						theNode->Release();
						theNode=NULL;
						if (resultat) SysFreeString(resultat);

						//UNDERLYING FIX LEG
						if (strcmp(ff1,"FIX") == 0)
						{
							fixFreqIfZC = GetCompFreq(item);
						}
					}
				}
			}

			if (item)
			{
				item->Release();
				item = NULL;
			}
			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}
			if (theNode)
			{
				theNode->Release();
				theNode = NULL;
			}
			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}
		}

		ARM_BermudaSwaptionCalculator* bermuda = new ARM_BermudaSwaptionCalculator(asOfDate,
																				  swaption,
																				  fixFreqIfZC);
		delete swaption;

		return bermuda;
	}
	catch(Exception& e)
	{
		throw e;
	}
	catch(...)
	{		
		CCString msg((CCString)"ParseBermudaSwaption : Error in XML parsing for getting BermudaSwaption");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
	}
}


ARM_CaptionCalculator* ARMLOCAL_ParseCaption(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	MSXML2::IXMLDOMDocument* XMLDoc = NULL;
 
 
	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	CCString structureId("");
 
 	string sBookName, sStructureId,sCustId,sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_CaptionCalculator* sec = NULL; 
	VECTOR<string> listAssetId;

 	sec = (ARM_CaptionCalculator*) ParserManager::BuildInstrument(	"ALMCAPTION", 
 																	XMLDoc, date, "",
 																	sBookName, sStructureId,
 																	sCustId, sDealId, listAssetId);

	bookName = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId = sCustId.c_str();
 	dealId = sDealId.c_str();
 	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(sec);
}


ARM_FxSpreadStripOption* ARMLOCAL_ParseFxStrip(	const char* chaineXML, 
												const ARM_Date& date, 	
												CCString& bookName,
												CCString& structureId,
												CCString& custId, 
												CCString& dealId,
												long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	IXMLDOMDocument* XMLDoc = NULL;
 
	structureId = "";

	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	string sBookName, sStructureId, sCustId, sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_FxSpreadStripOption* sec = NULL; 
	VECTOR<string> listAssetId;

 	sec = (ARM_FxSpreadStripOption*) ParserManager::BuildInstrument("FXSTRIP", 
 																	XMLDoc,
																	date,
																	"",
 																	sBookName,
																	sStructureId,
 																	sCustId,
																	sDealId,
																	listAssetId);

	bookName    = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId		= sCustId.c_str();
 	dealId		= sDealId.c_str();
	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(sec);

}

ARM_FXVanillaCalculator* ARMLOCAL_ParseFxStripCalculator(	const char* chaineXML, 
															const ARM_Date& date, 	
															CCString& bookName,
															CCString& structureId,
															CCString& custId, 
															CCString& dealId,
															long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	IXMLDOMDocument* XMLDoc = NULL;
 
	structureId = "";

	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);

        _bstr_t tmpChaine(wcharStr);
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	string sBookName, sStructureId, sCustId, sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_FXVanillaCalculator* calc = NULL; 
	VECTOR<string> listAssetId;

 	calc = (ARM_FXVanillaCalculator*) ParserManager::BuildInstrument("CFXSTRIP", 
 																	 XMLDoc,
																	 date,
																	 "",
 																	 sBookName,
																	 sStructureId,
 																	 sCustId,
																	 sDealId,
																	 listAssetId);

	bookName    = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId		= sCustId.c_str();
 	dealId		= sDealId.c_str();
	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(calc);

}

ARM_StdPortfolio* ARMLOCAL_ParseCorridorDblCondition(const char* chaineXML, 
												const ARM_Date& date, 	
												CCString& bookName,
												CCString& custId, 
												CCString& dealId,
												long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	IXMLDOMDocument* XMLDoc = NULL;
	try
 	{
 		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		wchar_t* wcharStr = constchar2wchar(chaineXML);
        _bstr_t tmpChaine(wcharStr);

		//on change <SWAP> en <EXOTIC>
		CCString myCCXml = tmpChaine;
		myCCXml.Replace("<SWAP>","<EXOTIC>");
		myCCXml.Replace("</SWAP>","</EXOTIC>");
		VariantTools::convert(CCSTringToSTLString(myCCXml),tmpChaine); 
		//Done
 
 		XMLDoc->loadXML(tmpChaine, &bOK);

        if (wcharStr)
           delete wcharStr;
 	}

 	catch(...)
 	{
 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
	CCString structureId("");
 	string sBookName, sStructureId,sCustId,sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_StdPortfolio* sec = NULL; 
	VECTOR<string> listAssetId;

	// On change les notations SWAP en EXOTIC dans la balise TradeType
	MSXML2::IXMLDOMNode* tmpNode = NULL;
	BSTR resultat = NULL;
	string exotic="EXOTIC";
	string resString = "";
	string query = "*/*/TradeType";
	XMLDoc->selectSingleNode((_bstr_t)query.c_str(), &tmpNode);
	if (tmpNode != NULL)
	{
		tmpNode->put_text((_bstr_t)exotic.c_str());
		tmpNode->Release();
		tmpNode = NULL;
	}
	//


 	sec = (ARM_StdPortfolio*) ParserManager::BuildInstrument("RNG_DOUBLE", 
 																	XMLDoc,
																	date,
																	"",
 																	sBookName,
																	sStructureId,
 																	sCustId,
																	sDealId,
																	listAssetId);

	bookName    = sBookName.c_str();
	structureId = sStructureId.c_str();
 	custId		= sCustId.c_str();
 	dealId		= sDealId.c_str();
	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	return(sec);

}


/*ARM_CallableSnowBallCalculator* ARMLOCAL_ParseCallableSnowBall(const char* chaineXML, const ARM_Date& date, CCString& bookName, CCString& custId, CCString& dealId, ARM_ReferenceValue& Premium, long isEtk)
{
 	HRESULT hr;
 	VARIANT_BOOL bOK;
 	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
 
 	wchar_t* xmlWCharText = NULL;
 
	try
 	{
 		hr = CoInitialize(NULL); 
 		SUCCEEDED(hr) ? 0 : throw hr;
 		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
 		SUCCEEDED(hr) ? 0 : throw hr;
 
 		xmlWCharText = constchar2wchar(chaineXML);
 
 		XMLDoc->loadXML((_bstr_t)xmlWCharText, &bOK);
 	}
 	catch(...)
 	{
 		if (xmlWCharText) delete xmlWCharText;
 		xmlWCharText = NULL;

 		if (XMLDoc) XMLDoc->Release();
 		hr = S_FALSE;
 
 		CCString msg((CCString)"Pb in creating XML document for getting Instrument from Summit");
 
 		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
 						 (const char *) msg);
 	}
 
 	CCString structureId("");
 
 	string sBookName, sStructureId,sCustId,sDealId;
 
 	ParserManager::Init(isEtk);
 
 	ARM_CallableSnowBallCalculator* sec = NULL; 
	VECTOR<string> listAssetId;

 	sec = (ARM_CallableSnowBallCalculator*) ParserManager::BuildInstrument(	"CSB", 
 																XMLDoc, date, "",
 																sBookName, sStructureId,
 																sCustId, sDealId, listAssetId);
	bookName = sBookName.c_str();
 	structureId = sStructureId.c_str();
 	custId = sCustId.c_str();
 	dealId = sDealId.c_str();
 	
 	if (XMLDoc)
 		XMLDoc->Release();
 
 	if (xmlWCharText) delete xmlWCharText;
 	xmlWCharText = NULL;
 
 	return sec;
}*/


ARM_VolCurve* ARMLOCAL_CreateNewFXVolFromSummitWithCurves(const CCString& ccy1,
															const CCString& ccy2,
															ARM_Date& Date,
															const CCString& cvName,
															ARM_ZeroCurve* domZc,
															ARM_ZeroCurve* forZc,
															double fxSpot,
															long WhatIsInterpolated,
															long correctSplineWithLinear,
															long isATM)
{
	ARM_VolLInterpol* PivotsAndInterp = NULL;


	PivotsAndInterp = etoolkit_GetPivotsAndInterpForFxVolCurve(ccy1,
			                                                   (CCString)"C"+ccy2,
															   cvName,
															   Date);

	if ( PivotsAndInterp == NULL )
    {
	   return etoolkit_CreateFXVolFromSummit(ccy1,
											 ccy2,
											 Date,
											 cvName,
											 "FXVOL",
											 "FXSPI");
    }

	ARM_FXVolCurve* newFxVolCrv = NULL;

	ARM_VolLInterpol* volATM = etoolkit_GetFXVolATMFromSummit(ccy1, ccy2, Date, 
                                                              cvName, "FXVOL");

	ARM_VolLInterpol* volCall = etoolkit_GetXMLFxSmileFromSummit(ccy1,
																 ccy2,
																 cvName,
																 Date,
																 "C",
																 K_DECREASING);

	ARM_VolLInterpol* volPut = etoolkit_GetXMLFxSmileFromSummit(ccy1,
																ccy2,
																cvName,
																Date,
																"P",
																K_DECREASING);

	int volSize = volATM->GetExpiryTerms()->GetSize();
	int pivotSize = PivotsAndInterp->GetExpiryTerms()->GetSize();

	if ( volSize != pivotSize )
	{
	   if (volATM)
		  delete volATM;
	   volATM = NULL;

	   if (volCall)
		  delete volCall;
	   volCall = NULL;

	   if (volPut)
		  delete volPut;
	   volPut = NULL;

	   if (PivotsAndInterp)
		  delete PivotsAndInterp;
	   PivotsAndInterp = NULL;

	   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					   "Vol and pivots don't have the same size");
	}

	for (int i = 0; i < volSize; i++)
	{
		if ( strcmp(volATM->itsYearTermsX[i],PivotsAndInterp->itsYearTermsX[i]) != 0 )
		{
			if (volATM)
				delete volATM;
			volATM = NULL;

			if (volCall)
				delete volCall;
			volCall = NULL;

			if (volPut)
				delete volPut;
			volPut = NULL;

			if (PivotsAndInterp)
				delete PivotsAndInterp;
			PivotsAndInterp = NULL;

			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							"Vol and pivots have not the same YearTerms");
		}
	}

	ARM_Vector* volFX =	volATM->GetVolatilities()->GetColumn(0);

//	FromRRSTRtoVol(volPut->GetVolatilities(),volCall->GetVolatilities(),volFX);

	ARM_Vector* tmpVector = PivotsAndInterp->GetVolatilities()->GetColumn(0);

	int size = tmpVector->GetSize();
	ARM_Vector* pivotType = new ARM_Vector(size);
	ARM_Vector* interpolTypesVect = new ARM_Vector(size);

	for (i = 0; i < volCall->GetStrikes()->GetSize(); i++)
	{
		volCall->GetStrikes()->Elt(i) /= 100.0;
	}

	for (i = 0; i < volPut->GetStrikes()->GetSize(); i++)
	{
		volPut->GetStrikes()->Elt(i) /= 100.0;
	}

	for (i = 0; i < size; i++)
	{
		pivotType->Elt(i) = (int) ((tmpVector->Elt(i) / 10 - 1));

		if ( (int) (tmpVector->Elt(i) - (pivotType->Elt(i) + 1)* 10) == 1)
			interpolTypesVect->Elt(i) = K_SPLINE;
		else if ( (int) (tmpVector->Elt(i) - (pivotType->Elt(i) + 1)* 10) == 2)
			interpolTypesVect->Elt(i) = K_LINEAR;
		else
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 "unexpected Interpol Method in new Fx Vol");
	}

	ARM_Vector inFxFwds(0);

	ARM_Currency currency((const char*) ccy1);
	forZc->SetCurrencyUnit( &currency );//new ARM_Currency( (char*) ccy1) );
	
	newFxVolCrv = new ARM_FXVolCurve(Date,
									 *(volATM->GetExpiryTerms()),
									 *volFX,
									 *pivotType,
									 *(volCall->GetStrikes()),
									 *(volCall->GetVolatilities()),
									 *(volPut->GetStrikes()),
									 *(volPut->GetVolatilities()),
									 inFxFwds,
									 *(interpolTypesVect),
									 WhatIsInterpolated,
									 correctSplineWithLinear,
									 fxSpot,
									 domZc,
									 forZc,
									 1,
									 isATM);

	for (i = 0; i < ARM_NB_TERMS; i++) 
	{
		sprintf(newFxVolCrv->itsYearTermsX[i], "%s", "X");

        sprintf(newFxVolCrv->itsYearTermsY[i], "%s", "X");
	}

	if ( volSize > ARM_NB_TERMS )
	{
		//TODO
		throw;
	}

	for (i = 0; i < volSize; i++)
	{
		sprintf(newFxVolCrv->itsYearTermsX[i], "%s", volATM->itsYearTermsX[i]);
	}

    // Update itsYearTermsY

    int idx = 0;
    sprintf(newFxVolCrv->itsYearTermsY[idx++], "%s", "ATM");

    for (i = 0; i < volPut->GetStrikes()->GetSize(); i++)
	{
		sprintf(newFxVolCrv->itsYearTermsY[idx], "%d",
                int(floor(volPut->GetStrikes()->Elt(i)*100.0)));
        idx++;
	}

    for (i = 0; i < volCall->GetStrikes()->GetSize(); i++)
	{
		sprintf(newFxVolCrv->itsYearTermsY[idx], "%d",
                int(floor(volCall->GetStrikes()->Elt(i)*100.0)));
        idx++;
	}

	if (volATM)
	   delete volATM;
	volATM = NULL;

	if (volCall)
	   delete volCall;
	volCall = NULL;

	if (volPut)
	   delete volPut;
	volPut = NULL;

	if (PivotsAndInterp)
	   delete PivotsAndInterp;
	PivotsAndInterp = NULL;

	if (interpolTypesVect)
	   delete interpolTypesVect;
	interpolTypesVect = NULL;

	if (pivotType)
	   delete pivotType;
	pivotType = NULL;

	if (tmpVector)
	   delete tmpVector;
	tmpVector = NULL;

	if (volFX)
	   delete volFX;
	volFX = NULL;

	// Update Mkt data characteristics    
    string	vType = string("FXVOL SMILE ");
	if(isATM)
		vType += "ATM";
	else
		vType += "IMPROVED";

	string	vIndex((const char*) ccy2);
	string	vCurrency((const char*) ccy1);
	string	vCrvId((const char*) cvName);

   newFxVolCrv->SetMktExternalCharacteristics(vIndex, vCurrency, vCrvId, vType);

	return	(newFxVolCrv);
}



ARM_VolCurve* ARMLOCAL_CreateNewFXVolFromSummitWithForwards(const CCString& ccy1,
															const CCString& ccy2,
															ARM_Date& Date,
															const CCString& cvName,
															const VECTOR<double>& Forwards,
															long WhatIsInterpolated,
															long correctSplineWithLinear,
															long isATM)
{
	return(NULL);	
}



ARM_VolLInterpol* etoolkit_GetPivotsAndInterpForFxVolCurve(const CCString& ccy,
														   const CCString& index,
														   const CCString& cvName,
														   ARM_Date date)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	ARM_Matrix* mVolATM = NULL;
	ARM_Vector* vYearTerms = NULL;

	ARM_VolLInterpol* newVolCrv = NULL;
	
	VECTOR<CCString>* matu = new VECTOR<CCString>;

	try
	{
		// ** GIGASPACE
		std::string xmlOutput ;
		ARM_XGigaToolKit::doGetMODELFACTOR((const char*)ccy,(const char*)index,"FXHD","SWOPT","FACTOR2",(const char*)cvName,date,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getFactor(ccy,
										 index,
										 cvName,
										 "SWOPT",
										 "FACTOR2",
										 date,
										 xmlResponse,
										 messageList);
		}

		retCode = ARMLOCAL_GetDatesAndDatasForFxVol(xmlResponse,
													date,
													ccy,
													mVolATM,
													vYearTerms,
													matu);

		ARM_Vector* vStrikes = new ARM_Vector(1,0.0);

		// construction de la matrice de vol forex (atm)
		
		newVolCrv = new ARM_VolLInterpol(date, vYearTerms,
									     vStrikes, mVolATM, 
                                         K_STK_TYPE_PRICE, K_ATMF_VOL);

		newVolCrv->SetIndexName((char *)(const char* ) index);

		for (int i = 0; i < ARM_NB_TERMS; i++) 
		{
			sprintf(newVolCrv->itsYearTermsX[i], "X");
		}

		if ( matu->size() > ARM_NB_TERMS )
		{
			return NULL;
		}
		
        for (i = 0; i < matu->size(); i++)
		{
			sprintf(newVolCrv->itsYearTermsX[i], "%s", (const char*) (*matu)[i]);
		}

		if (matu)
		   delete matu;
		matu = NULL;
	}

	catch(Exception&)
	{
		if (matu)
			delete matu;
		matu = NULL;

		return NULL;
	}

	catch (...)
	{
		if (matu)
			delete matu;
		matu = NULL;

		return NULL;
	}

	return newVolCrv;
}


ARM_Swap* ARMLOCAL_ParseSwap(const char* chaineXML,
							 const ARM_Date& date, 
							 CCString& bookName,
							 CCString& custId,
							 CCString& dealId,
							 VECTOR <string>& listAssetId,
							 long isEtk)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument *XMLDoc = NULL;
	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, *item = NULL;

	// Récupération de la structureId
	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Swap");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (char*) msg);
	}

	CCString structureId("");

	string sBookName, sStructureId,sCustId,sDealId;

	ParserManager::Init();
	//VECTOR<string> listAssetId;

	ARM_Swap* newSwap = (ARM_Swap*) ParserManager::BuildInstrument("SWAP", XMLDoc, date, "",
											 sBookName, sStructureId,
											 sCustId, sDealId, listAssetId);

	bookName = sBookName.c_str();
	structureId = sStructureId.c_str();
	custId = sCustId.c_str();
	dealId = sDealId.c_str();
	
	if (XMLDoc)
		XMLDoc->Release();

	return newSwap;
}


ARM_Object* ARMLOCAL_ParseCorridorSpreadOption(const char* chaineXML,
												 const ARM_Date& date, 
												 CCString& bookName,
												 const char* filter,
												 CCString& custId,
												 CCString& dealId,
												 VECTOR <string>& listAssetId,
												 long isEtk,
												 long isSwap)
{
	ARM_StdPortfolio* vPtf = NULL;
	vector<ARM_Security*> ListAssets;
//	CCString dealId, customerId, bookId;
//	vector<string>	vAssetIdList;
	vector<double> vWeights;
	vector<double> vMktPrices;
	vector<double> vVegas;
	double PoS = 1.0;

	if(isSwap)
	{
		ARM_Swap* newSwap = (ARM_Swap*) ARMLOCAL_ParseSwap(chaineXML,date,bookName,custId,dealId,listAssetId,isEtk);

		ARM_SwapLeg* firstSwapLeg = (ARM_SwapLeg*)newSwap->Get1stLeg()->Clone();
		ARM_SwapLeg* secondSwapLeg = (ARM_SwapLeg*)newSwap->Get2ndLeg()->Clone();
		delete newSwap;
		
		if(ARM_CorridorLeg* corridorLeg = dynamic_cast<ARM_CorridorLeg*>(firstSwapLeg))
		{
			MakeCorridorSpreadOption(corridorLeg,ListAssets,vWeights,vMktPrices,vVegas,dealId.c_str());
			delete corridorLeg;
		}
		else
		{
			PoS = firstSwapLeg->GetRcvOrPay();
			ListAssets.push_back((ARM_Security*)firstSwapLeg);
			vWeights.push_back(PoS);
			vVegas.push_back(1.0);
			vMktPrices.push_back(1.0);
		}

		PoS = secondSwapLeg->GetRcvOrPay();
		ListAssets.push_back((ARM_Security*)secondSwapLeg);
		vWeights.push_back(PoS);
		vVegas.push_back(1.0);
		vMktPrices.push_back(1.0);
	}
	else //GetType = EXOTIC
	{	
		CCString structureId("");

		ARM_StdPortfolio*	avPtf = (ARM_StdPortfolio*)ARMLOCAL_ParseExotic(chaineXML, date, filter, structureId,
																			bookName, custId, dealId, listAssetId);
		int size_Ptf = avPtf->GetSize();
		for(int i = 0; i < size_Ptf; i++)
		{
			

			if(ARM_CorridorLeg* corridorLeg = dynamic_cast<ARM_CorridorLeg*>(avPtf->GetAsset(i)))
			{
				MakeCorridorSpreadOption(corridorLeg,ListAssets,vWeights,vMktPrices,vVegas,dealId.c_str());
			}
			else
			{
				PoS = avPtf->GetAsset(i)->GetPorS();
				ListAssets.push_back((ARM_Security*)avPtf->GetAsset(i)->Clone());
				vWeights.push_back(PoS);
				vVegas.push_back(1.0);
				vMktPrices.push_back(1.0);
			}
		}
		
		delete avPtf;
	}

	vPtf = new ARM_StdPortfolio(ListAssets, vWeights, vMktPrices, vVegas);

	return vPtf;

	

}


ARM_SeasonalityManager* ARMLOCAL_ParseXMLForSeasonMgr(const CCString& Index,
													  const CCString& Ccy,
													  const CCString& CvName,
													  const ARM_Date& asof,
													  long modeId)
{
	CCString xmlResponse;
	CCString msgList;

	long retCode = ARM_OK;

	CCString postfix(4, Index);	
	CCString newIndex("S");
	newIndex += postfix;
	const char* type_vol = "IRFWDVOL";

	// ** GIGASPACE
	std::string xmlOutput ;
	ARM_XGigaToolKit::doGetVOLTenorsList(type_vol,(const char*)newIndex,(const char*)Ccy,(const char*)CvName,"IRG",asof,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		retCode = etoolkit_getcommsetname(newIndex,
										   Ccy,
										   CvName,
										   "IRG",
										   type_vol,
										   asof,
										   xmlResponse,
										   msgList);
	}

	if (( retCode == ARM_KO ) 
        || 
        ( strcmp((const char*)xmlResponse,"<Response></Response>") == 0 )
       )
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 "etoolkit_GetSeasonMgrFromSummit: no season Mgr found!");

        return(NULL);
	}

	VECTOR<CCString> listRequetes;
	VECTOR<CCString> listTenors;

	listTenors = ARMLOCAL_GetListTenorsFromXML(xmlResponse,&listRequetes, "P5");

	ARM_Vector vStrikes(listTenors.size());

	int i;
	int j;

	for (i = 0; i < listTenors.size(); i++)
	{
        vStrikes.Elt(i) = FromStrMatuToDouble((const char *) listTenors[i],(ARM_Date*) &asof);
    }

	VECTOR<CCString> sortListTenors;
	VECTOR<CCString> sortRequetes;

	ARM_Vector vSortStrikes(vStrikes.Sort());

	for (i = 0; i < vStrikes.GetSize(); i++)
	{
		j = 0;
		while ( j < vStrikes.GetSize() )
		{
			if ( vSortStrikes.Elt(i) != vStrikes.Elt(j) )
			{
			   j++;
			}
			else
			{
				sortListTenors.push_back(listTenors[j]);
				sortRequetes.push_back(listRequetes[j]);

				break;
			}
		}
	}

	CCString typeForGVC = type_vol;

	// ** GIGASPACE
	ARM_XGigaToolKit::doGetVOLTenor(type_vol,(const char*)newIndex,(const char*)Ccy,(const char*)CvName,"IRG",(const char*)sortListTenors[0],asof,xmlOutput);

	if (xmlOutput != "")
	{
		xmlResponse = (CCString)(xmlOutput.c_str());
	}
	else
	{
		retCode = etoolkit_getVolCurveByTenor(sortRequetes[0],
											  CvName,
											  asof,
											  typeForGVC,
											  type_vol,
											  xmlResponse,
											  msgList);
	}

	ARM_GP_Vector seasonSpreadList;
	vector<string> monthList;
	vector<double> horizonList;

	if (vSortStrikes.GetSize() > 1)
	{
		for (i = 0; i < vSortStrikes.GetSize(); i++)
			horizonList.push_back(vSortStrikes.Elt(i));
	}

	retCode = ARMLOCAL_GetMonthsAndDatasForSeasonMgr(xmlResponse,
													 seasonSpreadList,
													 monthList);

	for (i = 1; i < sortListTenors.size(); i++)
	{
		// ** GIGASPACE
		ARM_XGigaToolKit::doGetVOLTenor(type_vol,(const char*)newIndex,(const char*)Ccy,(const char*)CvName,"IRG",(const char*)sortListTenors[i],asof,xmlOutput);

		if (xmlOutput != "")
		{
			xmlResponse = (CCString)(xmlOutput.c_str());
		}
		else
		{
			retCode = etoolkit_getVolCurveByTenor(sortRequetes[i],
												  CvName,
												  asof,
												  typeForGVC,
												  type_vol,
												  xmlResponse,
												  msgList);
		}

		retCode = ARMLOCAL_GetMonthsAndDatasForSeasonMgr(xmlResponse,
														 seasonSpreadList,
														 monthList);
	}

	ARM_GP_Vector transposeSeasonSpreadList(seasonSpreadList);

	if (vSortStrikes.GetSize() > 1)
	{
		for (i = 0; i < monthList.size(); i++)
		{
			for (int j = 0; j < horizonList.size(); j++)
			{
				transposeSeasonSpreadList.Elt(i*horizonList.size()+j) = seasonSpreadList.Elt(j*monthList.size()+i);
			}
		}
	}

	ARM_SeasonalityManager* seasonMgr = NULL;

	if (vSortStrikes.GetSize() > 1)
	{
		seasonMgr = new ARM_SeasonalityManager(monthList,
											   transposeSeasonSpreadList,
											   horizonList,
											   (ARM_SeasonalityManager::CorrectionMode) modeId );
	}
	else
	{
		seasonMgr = new ARM_SeasonalityManager(monthList,
											   transposeSeasonSpreadList,
											   (ARM_SeasonalityManager::CorrectionMode) modeId );
	}

	return seasonMgr;
}


ARM_ResetManager* ParseReset(const ARM_Date& asof, const CCString& Source, const CCString& Ccy)
{
	if (Ccy.GetLen() != 7)
	{
		CCString msg((CCString)"Couple of currencies expected for FX (like 'EUR_USD')");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
	}

	string ccy1 = CCSTringToSTLString(Ccy).substr(0, 3);
	string ccy2 = CCSTringToSTLString(Ccy).substr(4, 3);

	if (ccy1 == ccy2)
	{
		CCString msg((CCString)"currencies should be different");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
	}

	vector<string> res;
	res.push_back("MKT");
	res.push_back("FX");
	res.push_back("CURRENCIES");
	res.push_back(CCSTringToSTLString(Ccy));
	
	string calendarStr = ccy1+ccy2;
	char* calendar = const_cast<char*>(calendarStr.c_str());

	ARM_Date date = asof;
	date.AddYears(-2); //temporairement 2 dernieres années
	while (asof > date.NextBusinessDay(calendar))
	{
		char dateStr[11];
		sprintf(dateStr, "%.0f", date.GetJulian());

		double value = ARMLOCAL_GetFxCurve(CCString(ccy1.c_str()), CCString(ccy2.c_str()), date, 1.0, Source);
		char valueStr[15];
		sprintf(valueStr, "%f", value);

		res.push_back(dateStr);
		res.push_back(valueStr);
	}

	return new ARM_ResetManager(res, res.size()/2, 2);
}

ARM_ResetManager* ParseReset(const CCString& xmlResponse, const ARM_Date& asof, const CCString& Index, const CCString& Ccy, const CCString& Term, long isInflatIndexId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument* XMLDoc = NULL;

	vector<string> vDates;
	vector<string> vValues;
	
	ARM_Date resDate;

	try
	{
		hr = CoInitialize(NULL);

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, 
                              CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();

        hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Reset Manager");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	MSXML2::IXMLDOMNodeList* resultList = NULL;
	MSXML2::IXMLDOMNode *listItem = NULL, *theNode = NULL;
	MSXML2::IXMLDOMNamedNodeMap *attributeMap=NULL;
	
	long nbNodes;

	int decimal, sign;
	
	try
	{
		if (XMLDoc->selectNodes((_bstr_t)((const char *)"Response/SummitQuery/xml/rs:data/*"), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
			{
				hr = S_FALSE;

				CCString msg((CCString)"Invalid XML string for getting Reset Manager");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (const char *) msg);
			}

			for (long indexNode = 0; indexNode<nbNodes; indexNode++)
			{
				hr = resultList->get_item(indexNode, &listItem);

				if (listItem->get_attributes(&attributeMap) == S_OK)
				{
					if (attributeMap->getNamedItem((_bstr_t)("c1"), &theNode) == S_OK)
					{
						VARIANT resultat;
						theNode->get_nodeValue(&resultat);

						_bstr_t ff(resultat.bstrVal,false);
						char * ff1=(char *)ff;

						if ((strcmp(ff1,"") != 0) && (strcmp(ff1,"SUN") != 0))
							resDate = ARM_Date(ff1,"YYYYMMDD");
					}
					theNode->Release();
					theNode = NULL;

					if ( ((resDate.GetDay() == 1) || (isInflatIndexId == K_NO)) && (resDate <= asof) )
					{
						char sDate[11];
						resDate.JulianToStrDate(sDate);
						ARM_Date aDate(resDate);
						vDates.push_back(ecvt(aDate.GetJulian(),7, &decimal, &sign));

						if (attributeMap->getNamedItem((_bstr_t)("c2"), &theNode) == S_OK)
						{
							VARIANT resultat;
							theNode->get_nodeValue(&resultat);

							_bstr_t ff(resultat.bstrVal,false);
							char * ff1=(char *)ff;

							vValues.push_back(ff1);
						}
						theNode->Release();
						theNode = NULL;
					}
					attributeMap->Release();
					attributeMap = NULL;
				}
				listItem->Release();
				listItem=NULL;
			}
		}
		if (resultList) resultList->Release();
	}
	catch(...)
	{
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Reset");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	if (XMLDoc) XMLDoc->Release();

	vector<string> res;

	res.push_back("MKT");

	if (isInflatIndexId == K_YES)
		res.push_back("INF");
	else
		res.push_back("IR");

	res.push_back("INDEX");

	CCString tmpIndex(Index);
	if (isInflatIndexId == K_NO)
	{
		if (tmpIndex == "EURIB")
			tmpIndex = "EURIBOR";
		else
			tmpIndex = "LIBOR";

		tmpIndex += Term;
		tmpIndex = Ccy + (const CCString&)"_" + tmpIndex;
	}

	res.push_back((const char*)tmpIndex);

	for (int i = 0; i < vDates.size(); i++)
	{
		res.push_back(vDates[i]);
		res.push_back(vValues[i]);
	}

	return new ARM_ResetManager(res,vDates.size()+2,2);
}


ARM_ResetManager* ARMLOCAL_ParseXMLForResetMgr(const ARM_Date& asof,
											   const CCString& Index,
											   const CCString& Source,
											   const CCString& Ccy,
											   long isInflatIndexId,
											   const CCString& Term)
{
	ARM_ResetManager* resetMgr = NULL;

	CCString xmlResponse;

	if (strcmp(Index, "FX") == 0)
	{
		resetMgr = ParseReset(asof, Source, Ccy);
	}
	else
	{
		xmlResponse = etoolkit_getXMLResetMgrFromSummit(Index, Source, Ccy, Term);
		if (strcmp(xmlResponse,"") != 0)
		{
			resetMgr = ParseReset(xmlResponse,asof,Index,Ccy,Term,isInflatIndexId);
		}
	}

	return resetMgr;
}


ARM_ReferenceValue* ARMLOCAL_ParseXMLCRFMeanRevParam(const CCString& xmlResponse, 
													 const ARM_Date& Date, 
													 const char* ccy, 
													 const char* index)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Vector* terms;
	ARM_Vector* values;
	ARM_ReferenceValue* res;

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		// JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting MRS parameter for CRF");
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * theDate = NULL, * theRate = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("ARMNT_CALL/CVLLIST/CURVEHEAD")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No MRS for CRF data");

			char key[100] = "";
			sprintf(key, "MODELFAC/MC/IRG/%s/%s/REVERSION", ccy, index);

			terms = new ARM_Vector();
			values = new ARM_Vector();

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{	
				hr=resultList->get_item(indexNode, &listItem);
				
				listItem->selectSingleNode(_bstr_t((const char *)"COM1"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					if (strcmp(ff1, key)==0)
					{
						listItem->selectSingleNode(_bstr_t((const char *)"ENTITYLIST/CURVE/Date"), &theDate);
						if (theDate!=NULL)
						{						
							BSTR resultat = NULL;
							theDate->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;
							ARM_Date& resDate = ARM_Date(ff1,"YYYYMMDD");

							terms->push_back((resDate.GetJulian()-Date.GetJulian())/365.);
						}
						listItem->selectSingleNode(_bstr_t((const char *)"ENTITYLIST/CURVE/Rate"), &theRate);
						if (theRate!=NULL)
						{
							BSTR resultat = NULL;
							theRate->get_text(&resultat);

							_bstr_t ff(resultat,false);
							char * ff1=(char *)ff;

							values->push_back(atof(ff1));
						}
					}
				}
			}
			if (terms->GetSize()==0)
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No MRS for CRF data");
			}
		}

		res = new ARM_ReferenceValue(terms, values);

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (theDate)
		{
			theDate->Release();
			theDate = NULL;
		}

		if (theRate)
		{
			theRate->Release();
			theRate = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();
		if (theDate) theDate->Release();
		if (theRate) theRate->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting MRS for CRF");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


VECTOR<CCString> ARMLOCAL_ParseListTenors(const CCString& index,
										  const CCString& ccy,
										  const CCString& cvname,
										  ARM_Date date,
										  const CCString& cvtype)
{
	VECTOR<CCString> listTenors;
	VECTOR<CCString> sortListTenors;

	CCString xmlResponse;

	xmlResponse = etoolkit_getListTenors(index, ccy, cvname, date, cvtype);

	if (strcmp(xmlResponse,"") != 0)
	{
		VARIANT_BOOL bOK;
		HRESULT hr;
		MSXML2::IXMLDOMDocument *XMLDoc = NULL;
		try
		{
			hr = CoInitialize(NULL); 
			hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
			SUCCEEDED(hr) ? 0 : throw hr;

			//JLA _bstr_t tmpChaine = xmlResponse;
			_bstr_t tmpChaine; 
			VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 


			XMLDoc->loadXML(tmpChaine, &bOK);
		}
		catch(...)
		{
			if (XMLDoc) XMLDoc->Release();
			hr = S_FALSE;
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Pb in creating XML document for getting Tenors List");
		}

		MSXML2::IXMLDOMNodeList * resultList = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * theDate = NULL, * theRate = NULL;
		long nbNodes;

		try
		{
			if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/TENOR_LIST/TENOR")), &resultList) == S_OK)
			{
				resultList->get_length(&nbNodes);

				if (nbNodes == 0)
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No SMILE data");

				for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
				{	
					hr=resultList->get_item(indexNode, &listItem);
					
					BSTR resultat = NULL;
					listItem->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;
					listTenors.push_back(ff1);

					if (listItem)
					{
						listItem->Release();
						listItem = NULL;
					}
				}
			}

			if (resultList)
			{
				resultList->Release();
				resultList = NULL;
			}

			if (XMLDoc)
			{
				XMLDoc->Release();
				XMLDoc = NULL;
			}

			hr = S_OK;
		}
		catch(...)
		{		
			if (listItem) listItem->Release();
			if (resultList) resultList->Release();
			if (XMLDoc) XMLDoc->Release();

			hr = S_FALSE;

			CCString msg((CCString)"Error in XML parsing for getting Tenor List");

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
							 (const char *) msg);
		}

		ARM_Vector vListTenor(nbNodes);
		for (int i = 0; i < nbNodes; i++)
			vListTenor.Elt(i) = FromStrMatuToDouble(listTenors[i]);

		ARM_Vector vSortTenor (vListTenor.Sort());

		for (i = 0; i < nbNodes; i++)
		{
			int j = 0;
			while ( j < vListTenor.GetSize() )
			{
				if ( vSortTenor.Elt(i) != vListTenor.Elt(j) )
				{
				   j++;
				}
				else
				{
					sortListTenors.push_back(listTenors[j]);				
					break;
				}
			}
		}
	}

	return sortListTenors;
}

ARM_ReferenceValue* etoolkit_getXMLModelParamFromSummit(const ARM_Date& date,
														const CCString& model,
														const CCString& type,
														const CCString& factorName,
														const CCString& ccy,
														const CCString& index,
														const CCString& cvName,
														long calcModId)
{
	CCString xmlResponse;
	CCString messageList;
	
	long retCode;

	try
	{
		retCode = etoolkit_getModelParam(date,
										 model,
										 type,
										 factorName,
										 ccy,
										 index,
										 cvName,
										 xmlResponse,
										 messageList);

		ARM_ReferenceValue* res = (ARM_ReferenceValue*) ARMLOCAL_ParseXMLForModelParam(xmlResponse,calcModId);

		// Recompute dates :
		// Add maturity and add spot days, then adjust
		ARM_Currency currency(CCSTringToSTLString(ccy).c_str());
		int spot = currency.GetSpotDays();

		ARM_Date settleDate = date;
		settleDate.NextBusinessDay(spot,CCSTringToSTLString(ccy).c_str());

		for (int i=0; i<res->GetDiscreteDates()->size(); i++)
		{
			ARM_Date tmpDate = settleDate;
			tmpDate.AddPeriod(YearTermToStringMatu(res->GetDiscreteDates()->Elt(i)).c_str(), ccy.c_str());
			tmpDate.AdjustToBusDate(ccy.c_str(), K_MOD_FOLLOWING);
			res->GetDiscreteDates()->Elt(i) = (tmpDate - date)/365.;
		}

		return res;
	}
	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		CCString msg((CCString)"Pb for getting model parameters");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}

ARM_ReferenceValue* ARMLOCAL_ParseXMLForModelParam(const CCString& xmlResponse, long calcModId)
{
	VARIANT_BOOL bOK;
	HRESULT hr;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	ARM_Vector* terms;
	ARM_Vector* values;
	ARM_ReferenceValue* res;

	try
	{
		hr = CoInitialize(NULL); 
		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;
		//JLA _bstr_t tmpChaine = xmlResponse;
		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;
		CCString msg((CCString)"Pb in creating XML document for getting model parameters");
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,(const char *) msg);
	}

	MSXML2::IXMLDOMNodeList * resultList = NULL;
	MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL;
	long nbNodes;

	try
	{
		if (XMLDoc->selectNodes(_bstr_t((const char *)("Response/Entity/COMMSET/CommData/COMMDATA")), &resultList) == S_OK)
		{
			resultList->get_length(&nbNodes);

			if (nbNodes == 0)
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No model param");

			terms = new ARM_Vector(nbNodes);
			values = new ARM_Vector(nbNodes);

			for (long indexNode=0 ; indexNode<nbNodes ; indexNode++)
			{	
				hr=resultList->get_item(indexNode, &listItem);
				
				listItem->selectSingleNode(_bstr_t((const char *)"Date"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					terms->Elt(indexNode) = StringMatuToYearTerm(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->selectSingleNode(_bstr_t((const char *)"Value"), &theNode);
				if (theNode!=NULL)
				{
					BSTR resultat = NULL;
					theNode->get_text(&resultat);

					_bstr_t ff(resultat,false);
					char * ff1=(char *)ff;

					values->Elt(indexNode) = atof(ff1);

					theNode->Release();
					theNode=NULL;
					if (resultat) SysFreeString(resultat);
				}
				listItem->Release();
				listItem=NULL;
			}	
		}

		res = new ARM_ReferenceValue(terms, values,0,1,calcModId);

		if (theNode)
		{
			theNode->Release();
			theNode = NULL;
		}

		if (resultList)
		{
			resultList->Release();
			resultList = NULL;
		}

		if (listItem)
		{
			listItem->Release();
			listItem = NULL;
		}

		if (XMLDoc)
		{
			XMLDoc->Release();
			XMLDoc = NULL;
		}

		hr = S_OK;

		return res;
	}

	catch(...)
	{		
		if (theNode) theNode->Release();
		if (listItem) listItem->Release();
		if (resultList) resultList->Release();
		if (XMLDoc) XMLDoc->Release();

		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting model parameter");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}


ARM_StripOption* ARMLOCAL_ParseFxOptionStrip(const char* chaineXML, const ARM_Date& date, 
											 CCString& bookName, CCString& structureId,
											 CCString& custId, CCString& dealId, long isEtk)
{
	VARIANT_BOOL bOK;
	HRESULT hr;

	MSXML2::IXMLDOMDocument*	XMLDoc = NULL;

	ARM_StripOption* FxOptionStrip = NULL;
	ARM_StripDigitalOption*	FxDigitalOptionStrip = NULL;

	structureId = "";

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		//JLA _bstr_t tmpChaine = chaineXML;
		_bstr_t tmpChaine; 
		VariantTools::convert(std::string(chaineXML),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in creating XML document for getting Fx Option Strip");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	try
	{
 		string sBookName, sStructureId,sCustId,sDealId;
 
 		ParserManager::Init(isEtk);
 
		VECTOR<string> listAssetId;

 		FxOptionStrip = (ARM_StripOption*) ParserManager::BuildInstrument ( "FXOPTSTRIP", 
 																			XMLDoc, date, "",
 																			sBookName, sStructureId,
 																			sCustId, sDealId, listAssetId );

		bookName = sBookName.c_str();
 		structureId = sStructureId.c_str();
 		custId = sCustId.c_str();
 		dealId = sDealId.c_str();
 		
 		if (XMLDoc)
 			XMLDoc->Release();

		return FxOptionStrip;
	}

	catch(...)
	{		
		hr = S_FALSE;

		CCString msg((CCString)"Error in XML parsing for getting Fx Option Strip");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}
}

void MakeCorridorSpreadOption(ARM_CorridorLeg* corridorLeg,
							 vector<ARM_Security*>& ListAssets,
							 vector<double>& vWeights,
							 vector<double>& vMktPrices,
							 vector<double>& vVegas,
							 const char* id)
{
	ARM_Date startDate = corridorLeg->GetStartDateNA();
	ARM_Date endDate   = corridorLeg->GetEndDateNA();
	
	ARM_IRIndex* PayIndex = corridorLeg->GetPaymentIndex();
	
	ARM_INDEX_TYPE liborType2 = corridorLeg->GetRefIndex()->GetIndexType();
	
	ARM_ReferenceValue* strikesUp      = corridorLeg->GetUpBarriers();
	ARM_ReferenceValue* strikesDown    = corridorLeg->GetDownBarriers();
	ARM_ReferenceValue* Fixing2        = corridorLeg->GetRefPastFixings();
	ARM_ReferenceValue* payIndexMargin = corridorLeg->GetSpreads();
	ARM_ReferenceValue* FixingPay      = corridorLeg->GetPayPastFixings();
	ARM_ReferenceValue* fee            = corridorLeg->GetFee();
	ARM_ReferenceValue* notional       = corridorLeg->GetAmount(); 

	ARM_ReferenceValue payFixedRate(0.0);
	ARM_ReferenceValue weight1(0.0);
	ARM_ReferenceValue weight2(1.0);
	
	ARM_Currency* discountCcy = corridorLeg->GetCurrencyUnit();

	int dayCount  = corridorLeg->GetDayCount();
	if(dayCount < 0)
		dayCount = corridorLeg->GetPaymentIndex()->GetDayCount();
	int resetFreq = corridorLeg->GetRefIndex()->GetResetFrequency();
	int payFreq   = corridorLeg->GetPaymentIndex()->GetPayFrequency();

	int resetTiming = corridorLeg->GetRefIndex()->GetResetTiming(); 
	int payTiming   = corridorLeg->GetPaymentIndex()->GetPayTiming();
	int resetGap    = corridorLeg->GetRefIndexResetGap();
	int intRule     = corridorLeg->GetRefIndex()->GetIntRule();
	int stubRule    = corridorLeg->GetStubMeth();
	
	char* resetCal = corridorLeg->GetResetCalName();
	char* payCal   = corridorLeg->GetPayCalName();
	
	 
	int PoS = corridorLeg->GetRcvOrPay();

	ARM_INDEX_TYPE liborType1 = ARM_INDEX_TYPE(K_CMS2);
	if(liborType2 == liborType1)
		liborType1 = ARM_INDEX_TYPE(K_CMS5);

	bool isCap = false, isFloor = false;
	int sizeCap = strikesDown->GetSize();
	int sizeFloor = strikesUp->GetSize();
	int i = 0, j = 0;
	while( (i < sizeCap) && (strikesDown->GetDiscreteValues()->Elt(i) < 0.05))
		i++;
	if( i!=sizeCap )
		isCap = true;

	while((j < sizeFloor) && (strikesUp->GetDiscreteValues()->Elt(j) >= 100.) )
		j++;
	if( j!=sizeFloor )
		isFloor = true;

	if(isCap && isFloor)
	{
		ARM_SpreadOption* spreadOptionUp = new ARM_SpreadOption(startDate,endDate, K_FLOOR, strikesUp,
															   PayIndex, &payFixedRate, liborType1, liborType2,  
															   &weight1,&weight2, dayCount,resetFreq,payFreq,
															   resetTiming, payTiming, discountCcy, resetGap,0.01,
															   -0.01, NULL, Fixing2, intRule, stubRule, resetCal,
															   payCal,1,1, payIndexMargin,1.0,
															   FixingPay);
		ARM_SwapLeg* Leg = spreadOptionUp->GetPayIndexLeg();
		if(strikesUp->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = strikesUp->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newStrikeUp = new ARM_ReferenceValue(resetDates,refValues);
			if(newStrikeUp)
				spreadOptionUp->SetStrikes(newStrikeUp);
			delete newStrikeUp;			
		}
		if(payIndexMargin->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = payIndexMargin->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newpayIndexMargin = new ARM_ReferenceValue(resetDates,refValues);
			if(newpayIndexMargin)
				spreadOptionUp->SetPayIndexMargins(newpayIndexMargin);
			delete newpayIndexMargin;		
		}

		spreadOptionUp->SetAmount(notional);
		spreadOptionUp->SetPorS(PoS);
		spreadOptionUp->SetFee(fee);
		spreadOptionUp->SetAssetId(id);

		ARM_SpreadOption* spreadOptionDown = new ARM_SpreadOption(startDate,endDate, K_FLOOR, strikesDown,
															   PayIndex, &payFixedRate, liborType1, liborType2,  
															   &weight1,&weight2, dayCount,resetFreq,payFreq,
															   resetTiming, payTiming, discountCcy, resetGap,0.01,
															   -0.01, NULL, Fixing2, intRule, stubRule, resetCal,
															   payCal,1,1, payIndexMargin,1.0,
															   FixingPay);

		
		Leg = spreadOptionUp->GetPayIndexLeg();
		if(strikesDown->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = strikesDown->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newStrikeDown = new ARM_ReferenceValue(resetDates,refValues);
			if(newStrikeDown)
				spreadOptionDown->SetStrikes(newStrikeDown);
			delete newStrikeDown;		
		}
		if(payIndexMargin->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = payIndexMargin->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newpayIndexMargin = new ARM_ReferenceValue(resetDates,refValues);
			if(newpayIndexMargin)
				spreadOptionDown->SetPayIndexMargins(newpayIndexMargin);
			delete newpayIndexMargin;		
		}
		
		spreadOptionDown->SetAmount(notional);
		spreadOptionDown->SetPorS(-PoS);
	//	spreadOptionDown->SetFee(fee);
		spreadOptionDown->SetAssetId(id);

		ListAssets.push_back((ARM_Security*)spreadOptionUp);
		vWeights.push_back(1);
		vVegas.push_back(1.0);
		vMktPrices.push_back(1.0);

		ListAssets.push_back((ARM_Security*)spreadOptionDown);
		vWeights.push_back(1);
		vVegas.push_back(1.0);
		vMktPrices.push_back(1.0);
	}
	else if(isCap)
	{
		ARM_SpreadOption* spreadOptionDown = new ARM_SpreadOption(startDate,endDate, K_CAP, strikesDown,
															   PayIndex, &payFixedRate, liborType1, liborType2,  
															   &weight1,&weight2, dayCount,resetFreq,payFreq,
															   resetTiming, payTiming, discountCcy, resetGap,0.01,
															   -0.01, NULL, Fixing2, intRule, stubRule, resetCal,
															   payCal,1,1, payIndexMargin,1.0,
															   FixingPay);

		ARM_SwapLeg* Leg = spreadOptionDown->GetPayIndexLeg();
		if(strikesDown->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = strikesDown->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newStrikeDown = new ARM_ReferenceValue(resetDates,refValues);
			if(newStrikeDown)
				spreadOptionDown->SetStrikes(newStrikeDown);
			delete newStrikeDown;		
		}
		if(payIndexMargin->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = payIndexMargin->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newpayIndexMargin = new ARM_ReferenceValue(resetDates,refValues);
			if(newpayIndexMargin)
				spreadOptionDown->SetPayIndexMargins(newpayIndexMargin);
			delete newpayIndexMargin;	
		}

		spreadOptionDown->SetAmount(notional);
		spreadOptionDown->SetPorS(PoS);
		spreadOptionDown->SetFee(fee);
	//	spreadOptionDown->SetAssetId(id);

		ListAssets.push_back((ARM_Security*)spreadOptionDown);
		vWeights.push_back(1);
		vVegas.push_back(1.0);
		vMktPrices.push_back(1.0);
	}
	else
	{
		ARM_SpreadOption* spreadOptionUp = new ARM_SpreadOption(startDate,endDate, K_FLOOR, strikesUp,
															   PayIndex, &payFixedRate, liborType1, liborType2,  
															   &weight1,&weight2, dayCount,resetFreq,payFreq,
															   resetTiming, payTiming, discountCcy, resetGap,0.01,
															   -0.01, NULL, Fixing2, intRule, stubRule, resetCal,
															   payCal,1,1, payIndexMargin,1.0,
															   FixingPay);

		ARM_SwapLeg* Leg = spreadOptionUp->GetPayIndexLeg();
		if(strikesUp->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = strikesUp->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newStrikeUp = new ARM_ReferenceValue(resetDates,refValues);
			if(newStrikeUp)
				spreadOptionUp->SetStrikes(newStrikeUp);
			delete newStrikeUp;			
		}
		if(payIndexMargin->GetSize())
		{
			ARM_Vector* startDates = Leg->GetFlowStartDates();
			ARM_Vector* resetDates = (ARM_Vector*)Leg->GetResetDates()->Clone();
			int vSize = resetDates->GetSize();
			ARM_Vector* refValues = new ARM_Vector(vSize);
			for(int k = 0; k < vSize; k++)
				refValues->Elt(k) = payIndexMargin->CptReferenceValue(startDates->Elt(k));
			ARM_ReferenceValue* newpayIndexMargin = new ARM_ReferenceValue(resetDates,refValues);
			if(newpayIndexMargin)
				spreadOptionUp->SetPayIndexMargins(newpayIndexMargin);
			delete newpayIndexMargin;		
		}

		spreadOptionUp->SetAmount(notional);
		spreadOptionUp->SetPorS(PoS);
		spreadOptionUp->SetFee(fee);
		spreadOptionUp->SetAssetId(id);
		
		ListAssets.push_back((ARM_Security*)spreadOptionUp);
		vWeights.push_back(1);
		vVegas.push_back(1.0);
		vMktPrices.push_back(1.0);
	}
}





ARM_Object* ARMLOCAL_ParseCalypsoObject(string xmlContent, const string type,const string modelType, const ARM_Date& date)
{

	try
	{
		string bookName, custId, dealId;
		
		if (!type.compare("CRF"))
			return ARMLOCAL_ParseCalypsoCRF(xmlContent, date,modelType,bookName, custId, dealId);
		
		else if (!type.compare("CDS"))
			return ARMLOCAL_ParseCalypsoCDS(xmlContent,
											date,
											modelType,
											bookName,
											custId);
		else if (!type.compare("NTD"))
			return ARMLOCAL_ParseCalypsoNTD(xmlContent,
											date,
											modelType,
											bookName,
											custId);
		
		else if (!type.compare("BERM"))
			return ARMLOCAL_ParseCalypsoBermudaSwaption(xmlContent, date,modelType,bookName, custId, dealId);
		else if (!type.compare("CDO"))
			return ARMLOCAL_ParseCalypsoCDO(xmlContent,
											date,
											modelType,
											bookName,
											custId);
		else
			return NULL;
		}
	catch (Exception& e)
	{
		throw e;
	}

	catch (...)
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			 "Error Unknown in XML parsing for getting Object");
	}

}



ARM_CRFCalculator* ARMLOCAL_ParseCalypsoCRF(string chaineXML,
									 const ARM_Date& asOfDate,
									 const string modeType,
									 string& bookName,
									 string& custId,
									 string& dealId)

{
    //vModelType
	ARM_PricingModelType::ModelType	modelType = ((modeType=="CRF_HW")?ARM_PricingModelType::HWM1F:ARM_PricingModelType::QGM1F);
	return ARMLOCAL_ParseCalypsoCRF(chaineXML, asOfDate, modelType,	bookName, custId, dealId);

}


ARM_CRFCalculator* ARMLOCAL_ParseCalypsoCRF(string chaineXML,
									 const ARM_Date& asOfDate,
									 ARM_PricingModelType::ModelType	modelType,
									 string& bookName,
									 string& custId,
									 string& dealId)

{
	ARM_ReferenceValue* exerFee = NULL;


	MSXML2::IXMLDOMDocumentPtr xmlDoc; 
	xmlDoc = XMLTools::LoadXML(chaineXML); 

    MSXML2::IXMLDOMNodePtr xmlNode;
	
    
    ARM_Date tempDate,startDate,endDate;
	xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows");
	
    XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"payLeg/startDate"),startDate,"YYYYMMDD");
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"receiveLeg/startDate"),tempDate,"YYYYMMDD");
    startDate=(startDate>tempDate)?tempDate:startDate;

	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"payLeg/endDate"),endDate,"YYYYMMDD");
	XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"receiveLeg/endDate"),tempDate,"YYYYMMDD");
	endDate=(endDate>tempDate)?endDate:tempDate;

    int buySell =1;
    double quantity; XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/quantity"),quantity);;
    if (quantity >=0 ) 
    	buySell = -1;
    buySell= calypso2ARMPayOrRec(buySell);


	string payLegType;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"payLeg/fixFloat"),payLegType);
	
    string fundingLeg= "payLeg";
    string structuredLeg= "receiveLeg";
	long payRec=  K_RCV;
    if((payLegType == "Fix")||(payLegType == "Form")){
    	payRec = K_PAY;
	   	fundingLeg = "receiveLeg";
       	structuredLeg = "payLeg";
    }

    string sProductDefinition;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/ProductDefinition"),sProductDefinition);
	//if(sProductDefinition =="CancellableSwap")
	if(sProductDefinition !="Swaption")
        payRec*=-buySell;


    xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/"+fundingLeg);
    CalypsoLeg fundLeg = getSwapLeg(xmlNode);
	xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/"+structuredLeg);
    CalypsoLeg cpnLeg = getCRFLeg(xmlNode);

	// Structured ARM Reference values initilization
	ARM_Vector* armStrikeDates = new ARM_Vector(cpnLeg.flowStartDatesVector);
	ARM_Vector* armStrikeValues = new ARM_Vector(cpnLeg.strikeVector);
	ARM_ReferenceValue* strike = new ARM_ReferenceValue(armStrikeDates, armStrikeValues);
	strike->SetCalcMethod(K_STEPUP_LEFT);

	ARM_ReferenceValue* leverage = new ARM_ReferenceValue(0.);
	if (cpnLeg.leverageVector.size() > 0)
	{
		ARM_Vector* armLeverageDates = new ARM_Vector(cpnLeg.flowStartDatesVector);
		ARM_Vector* armLeverageValues = new ARM_Vector(cpnLeg.leverageVector);
		leverage = new ARM_ReferenceValue(armLeverageDates, armLeverageValues);
	}

	ARM_ReferenceValue* cpnMin = new ARM_ReferenceValue(0.);
	if (cpnLeg.cpnMinVector.size() > 0)
	{
		ARM_Vector* armCpnMinDates = new ARM_Vector(cpnLeg.flowStartDatesVector);
		ARM_Vector* armCpnMinValues = new ARM_Vector(cpnLeg.cpnMinVector);
		cpnMin = new ARM_ReferenceValue(armCpnMinDates, armCpnMinValues);
	}

	ARM_ReferenceValue* cpnMax = new ARM_ReferenceValue(0.);
	if (cpnLeg.cpnMaxVector.size() > 0)
	{
		ARM_Vector* armCpnMaxDates = new ARM_Vector(cpnLeg.flowStartDatesVector);
		ARM_Vector* armCpnMaxValues = new ARM_Vector(cpnLeg.cpnMaxVector);
		cpnMax = new ARM_ReferenceValue(armCpnMaxDates, armCpnMaxValues);
	}

	ARM_Vector* armNominalDates = new ARM_Vector(cpnLeg.flowStartDatesVector);
	ARM_Vector* armNominalValues = new ARM_Vector(cpnLeg.flowNominalVector);
	ARM_ReferenceValue* nominal = new ARM_ReferenceValue(armNominalDates, armNominalValues);

	// Funding ARM Reference values initilization
	ARM_Vector* armFundNominalDates = new ARM_Vector(fundLeg.flowStartDatesVector);
	ARM_Vector* armFundNominalValues = new ARM_Vector(fundLeg.flowNominalVector);
	ARM_ReferenceValue* fundNominal	= new ARM_ReferenceValue(armFundNominalDates, armFundNominalValues);

	ARM_Vector* armFundSpreadDates = new ARM_Vector(fundLeg.flowStartDatesVector);
	ARM_Vector* armFundSpreadValues = new ARM_Vector(fundLeg.flowSpreadVector);
	ARM_ReferenceValue* fundSpread = new ARM_ReferenceValue(armFundSpreadDates, armFundSpreadValues);

	// Fees ARMReference_Value initialization
    xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/ListeputCalldate");
    exerFee = convertCalypsoOptionSchedule(xmlNode);

    double meanReversion = 0;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/IXIS_MeanRev"),meanReversion);
	meanReversion/=100.;

    //vModelType
//	ARM_PricingModelType::ModelType	modelType = ((modeType=="CRF_HW")?ARM_PricingModelType::HWM1F:ARM_PricingModelType::QGM1F);
	

	ARM_CRFCalculator* crfCalculator = new ARM_CRFCalculator(	asOfDate,
																startDate,
																endDate,
																*strike,
																payRec,
																cpnLeg.fixEndDate,
																cpnLeg.dayCount,	// fix daycount
																cpnLeg.dayCount,
																cpnLeg.freq,
																cpnLeg.resetTiming,
																cpnLeg.indexTerm,
																cpnLeg.indexDayCount,
																cpnLeg.resetCal,
																cpnLeg.payCal,
																cpnLeg.stubRule,
																cpnLeg.resetGap,
																*leverage,
																*cpnMin,
																*cpnMax,
																*fundSpread,
																fundLeg.freq,
																fundLeg.dayCount,
																*nominal,
																*exerFee,
																cpnLeg.currency,
																fundLeg.currency,
																cpnLeg.currency,
																cpnLeg.currency,
																fundLeg.currency,
																*fundNominal,
																modelType);

	crfCalculator->SetPorS(buySell);
	// Hedge management 
	/*mid = jni->GetMethodID(cls, "getPricingParamBooleanValue", "(Ljava/lang/String;)Z");
	jboolean strikeEquiv = jni->CallBooleanMethod(trade, mid, jni->NewStringUTF("CRF.SWAPTION_CALIB"));
	if(strikeEquiv) {
		crfCalculator->SetOSWCalibFlag(ARM::ARM_SigmaCalibrationType::strikeEquivalent);
	} else {
		crfCalculator->SetOSWCalibFlag(ARM::ARM_SigmaCalibrationType::Unknown);
	}	
	crfCalculator->SetCapCalibFlag((jni->CallBooleanMethod(trade, mid, jni->NewStringUTF("CRF.CAP_CALIB")))?true:false);
	crfCalculator->SetFloorCalibFlag((jni->CallBooleanMethod(trade, mid, jni->NewStringUTF("CRF.FLOOR_CALIB")))?true:false);
	crfCalculator->SetSkewReCalibFlag((jni->CallBooleanMethod(trade, mid, jni->NewStringUTF("CRF.SKEW_RECALIB")))?true:false);
	*/

    	if(true) {
			crfCalculator->SetOSWCalibFlag(ARM::ARM_SigmaCalibrationType::strikeEquivalent);
		} else {
			crfCalculator->SetOSWCalibFlag(ARM::ARM_SigmaCalibrationType::Unknown);
		}	
		crfCalculator->SetCapCalibFlag(true);
		crfCalculator->SetFloorCalibFlag(false);
		crfCalculator->SetSkewReCalibFlag(false);


	if (modelType == ARM::ARM_PricingModelType::HWM1F && meanReversion != 0)
	{
		ARM::ARM_CurveModelParam MRParam(ARM::ARM_ModelParamType::MeanReversion, meanReversion, "MRS");
		crfCalculator->SetMRS(&MRParam);
	}

	if (strike) delete strike;
	if (leverage) delete leverage;
	if (cpnMin) delete cpnMin;
	if (cpnMax) delete cpnMax;
	if (fundSpread) delete fundSpread;
	if (nominal) delete nominal;
	if (exerFee) delete exerFee;
	if (fundNominal) delete fundNominal;

	return crfCalculator;

}

//new version
ARM_Object* ARMLOCAL_ParseCalypsoBermudaSwaption(string chaineXML,
									 const ARM_Date& asOfDate,
									 const string modeType,
									 string& bookName,
									 string& custId,
									 string& dealId)
{
   //ModelType
    ARM_PricingModelType::ModelType	modelType = ((modeType=="BERM_HW")?ARM_PricingModelType::HWM1F:ARM_PricingModelType::QGM1F);
    if(modelType == ARM_PricingModelType::HWM1F )
         return ARMLOCAL_ParseCalypsoBermudaSwaption(chaineXML, asOfDate, modelType, bookName, custId, dealId);
    else
         return ARMLOCAL_ParseCalypsoCRF(chaineXML, asOfDate, modelType,	bookName, custId, dealId);
}
    
ARM_BermudaSwaptionCalculator* ARMLOCAL_ParseCalypsoBermudaSwaption(string chaineXML,
									 const ARM_Date& asOfDate,
									 ARM_PricingModelType::ModelType	modelType,
									 string& bookName,
									 string& custId,
									 string& dealId)
{
    int rcvPay;

	ARM_SwapLeg* fixLeg = NULL;
	ARM_SwapLeg* varLeg = NULL;
	ARM_Swap* swap = NULL;
	ARM_ExerciseStyle* exerStyle = NULL;
	ARM_ReferenceValue* spreads = NULL;
	ARM_ReferenceValue* strikes = NULL;
	ARM_ReferenceValue* fees = NULL;
	ARM_ReferenceValue* fundNotional = NULL;
	ARM_Swaption* newSwaption = NULL;

	ARM_BermudaSwaptionCalculator* newCalculator = NULL;



    MSXML2::IXMLDOMDocumentPtr xmlDoc; 
    MSXML2::IXMLDOMNodePtr xmlNode;
    xmlDoc = XMLTools::LoadXML(chaineXML); 
	
    
	// Start Date
	ARM_Date startDate; XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/StartDate"),startDate,"YYYYMMDD");
    
    // End Date
	ARM_Date endDate; XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/EndDate"),endDate,"YYYYMMDD");
        
	// As Of Date
	
	// Buy / Sell
	int buySell;
	double quantity; XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/quantity"),quantity);;
    if (quantity >=0 ) 
   		buySell = -1;
    else 
  		buySell = 1;
	buySell= calypso2ARMPayOrRec(buySell);


	// Receive or Pay
    string couponLeg;
	string fundingLeg;

	string payLegType;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/payLeg/fixFloat"),payLegType);


	if (payLegType == "Fix") {
		rcvPay = K_PAY;
		couponLeg = "payLeg";
	    fundingLeg = "receiveLeg";
        	
    }else {
		rcvPay = K_RCV;
		couponLeg = "receiveLeg";
    	fundingLeg = "payLeg";
    }
    //To modify update 1 TODO to finalize
    string sProductDefinition;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Description/ProductDefinition"),sProductDefinition);
	if(sProductDefinition =="CancellableSwap")
	rcvPay*=-buySell;

	/************* FIXED LEG ***********************/
	xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/"+couponLeg);
    fixLeg = convertFixLeg(xmlNode, startDate, endDate, rcvPay, strikes);
    
	/************* FLOAT LEG ***********************/
	xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/"+fundingLeg);
    varLeg = convertVarLeg(xmlNode, startDate, endDate, rcvPay, spreads,  fundNotional);
    
	if (spreads)
		varLeg->SetVariableSpread(spreads);

	swap = new ARM_Swap(fixLeg, varLeg);

	/************* EXERCICES & FEES ***********************/
	xmlNode = XMLTools::selectSingleNode(xmlDoc,"/Extraction/ListTrade/Trade/Flows/ListeputCalldate");
    convertCalypsoOptionSchedule(xmlNode, fees, exerStyle);

	newSwaption = new ARM_Swaption(swap, 
								   rcvPay,
								   exerStyle,
								   strikes,
								   0);

	newSwaption->SetFee(fees);

	newSwaption->SetAmount(fundNotional);
	newSwaption->SetPorS(buySell);

	int fixFreqIfZC = K_DEF_FREQ;
	// si la fixleg est ZC, on récupère en plus la frequence COMP_CompFreq
	if (newSwaption->GetFixedLeg()->GetPaymentFreq() == K_ZEROCOUPON) {
		fixFreqIfZC = fixLeg->GetDecompFreq();
	}

	newCalculator = new ARM_BermudaSwaptionCalculator(asOfDate, newSwaption, fixFreqIfZC);

	if(fundNotional) {
		delete fundNotional;
		fundNotional = NULL;
	}
	if(strikes) {
		delete strikes;
		strikes = NULL;
	}
	if(spreads) {
		delete spreads;
		spreads = NULL;
	}
	if(varLeg) {
		delete varLeg;
		varLeg = NULL;
	}
	if(fixLeg) {
		delete fixLeg;
		fixLeg = NULL;
	}
	if(swap) {
		delete swap;
		swap = NULL;
	}
	if(exerStyle) {
		delete exerStyle;
		exerStyle = NULL;
	}
	if(fees) {
		delete fees;
		fees = NULL;
	}
	if(newSwaption) {
		delete newSwaption;
		newSwaption = NULL;
	}

	return newCalculator;
}




// new function to get vol from Calypso
/**
Get data and size from xml
  **/

  //TODO deprecated
  double FromCalypsoStrMatuToDouble(string matu)
{
    if((matu.size()==8)&&!(isalpha(matu.c_str()[matu.size()-1]))){
        return (new ARM_Date(matu.c_str(),"YYYYMMDD"))->GetJulian();

    }else{
        return FromStrMatuToDouble(matu.c_str());

    }
    //return convPlotInYearTerm((char *)ff1,  ARM_Date asof,ARM_Date settleDate,char[] payCal);
}

  // Fonction de conversion des chaines (2D, 1M, Contrats Futur) en fraction d'année
double convCalypsoPlotInYearTerm(const CCString& inPlot, ARM_Date AsOf, ARM_Date Settle, 
                          char* calenDar)
{
	long isDate = 0;
	long Nb;
	char cMatu;
	long freqId;

    char plot[50];

    char* currency = calenDar;

    strcpy(plot, (const char *) inPlot);

	// Contrat
	if ( strlen((const char *) plot) == 5 )
	{
		int month, year;
		ARM_Date matDate;

		GetMonthYearFromExpiryDate(plot, &month, &year);
		matDate.ChangeDate(1, month, year);

		matDate.PiborDelivery();
		return ( (matDate.GetJulian() - AsOf.GetJulian()) /365.);
	}
	else
	{
		if ( strlen((const char *) plot) >= 8 )
		{
			isDate = 1;
		}
		else
		{
			sscanf(plot, "%ld%c", &Nb, &cMatu);

			cMatu = toupper(cMatu);

			if ( cMatu == 'D' ) // Ex : "1D"
				freqId = K_DAILY;
			else if ( cMatu == 'W' )
				freqId = K_WEEKLY;
			else if ( cMatu == 'M' ) 
				freqId = K_MONTHLY;
			else if ( cMatu == 'Y')  // ->implicitement ce sont des taux de swap
				freqId = K_ANNUAL;
			else
				isDate = 1;
		}
		
		if ( isDate == 1 )
		{
				if(plot[2] == '/')
				{
					ARM_Date tmpDate(plot, "DD/MM/YYYY");

					return ((tmpDate-AsOf)/365.);
				}
				else
				{
					ARM_Date tmpDate(plot, "YYYYMMDD");

					return ((tmpDate-AsOf)/365.);
				}
			
		}
		else
		{
			if ( freqId == K_DAILY )
			{
				Settle.AddDays(Nb);
				Settle.AdjustToBusDate(currency, K_FOLLOWING);
			}
			else
			{
				Settle.AddPeriodMult(freqId, (long) Nb, currency);
				Settle.AdjustToBusDate(currency, K_FOLLOWING);
			}

			return ( (Settle - AsOf) /365.);
		}
	}	
}

void loadData(MSXML2::IXMLDOMDocumentPtr xmlDoc, vector<string> *expiryVector, vector<string> *tenorVector, 
						vector<double> *strikeVector, vector<double> *valueVector ,
						set<string> *expirySet, set<string> *tenorSet,set<double> *strikeSet,
                        vector<string> *expirySetVector, vector<string> *tenorSetVector,vector<double> *strikeSetVector){
	
	MSXML2::IXMLDOMNodeListPtr xmlNodeList = XMLTools::selectNodes(xmlDoc,"/ListeVolSurface/VolSurface/Points/Point"); 
	long nbPoint=0;	
	
	xmlNodeList->get_length(&nbPoint); 

	string sexpiry;
	string stenor;
	double dstrike;
	double dvalue;

	for(int i=0;i<nbPoint;i++) 
	{
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Expiry"),sexpiry);
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Tenor"),stenor);
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Strike"),dstrike);
		XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Value"),dvalue);
			
		expiryVector->push_back(sexpiry);
		tenorVector->push_back(stenor);
		strikeVector->push_back(dstrike);
		valueVector->push_back(dvalue);

        if(expirySet->find(sexpiry)== expirySet->end( )){
			expirySet->insert(sexpiry);
            expirySetVector->push_back(sexpiry);
		}
			
        if(tenorSet->find(stenor)== tenorSet->end( )){
			tenorSet->insert(stenor);
            tenorSetVector->push_back(stenor);
		}

        if(strikeSet->find(dstrike)== strikeSet->end( )){
			strikeSet->insert(dstrike);
            strikeSetVector->push_back(dstrike);
		}
	}

}



/**
GetVolSmileFromCalypso
**/

	
ARM_VolCurve* GetVolSmileFromCalypso(ARM_Date asOfDate,ARM_Currency* armCurrency, vector<double> valueVector ,
						vector<string> expirySetVector,vector<double> strikeSetVector, double tenorCount, double tenorRank){


	/*jclass cls = jni->GetObjectClass(surface);
	jmethodID mid = jni->GetMethodID(cls, "getDate", "()J");
	jlong juliandate = jni->CallLongMethod(surface, mid);
	ARM_Date asOfDate = ARM_Date(juliandate);*/

	
	long juliandate =  asOfDate.GetJulian();

    long spotDays = armCurrency->GetSpotDays();
    char* payCalTmp = armCurrency->GetPayCalName(armCurrency->GetVanillaIndexType());
    char payCal[30];
    strcpy(payCal, payCalTmp);
    delete payCalTmp;

	ARM_Date settleDate = asOfDate;
	settleDate.NextBusinessDay(spotDays, payCal);

	//mid = jni->GetMethodID(cls, "getExpiryDateCount", "()I");
	//jint expDateCount = jni->CallIntMethod(surface, mid);
	int expDateCount = expirySetVector.size();

	//mid = jni->GetMethodID(cls, "getStrikeCount", "()I");
	//jint strikeCount = jni->CallIntMethod(surface, mid);
	int strikeCount = strikeSetVector.size();


	vector<double> expDates;
	//mid = jni->GetMethodID(cls, "getExpiryDate", "(I)J");
	//set <string> :: iterator expiry_Iter;
	//expiry_Iter= expirySet.begin();
	int i;
	for (i = 0; i < expDateCount; i++) {
		//jlong date = jni->CallLongMethod(surface, mid, i);
		//expDates.push_back((date-juliandate)/365.);
		//expDates.push_back((FromStrMatuToDouble((*expiry_Iter++).c_str())-juliandate)/365.);
        //expDates.push_back((FromCalypsoStrMatuToDouble((*expiry_Iter++))-juliandate)/365.);
        //expDates.push_back((FromCalypsoStrMatuToDouble((expirySetVector.at(i)))-juliandate)/365.);
        expDates.push_back(convCalypsoPlotInYearTerm((expirySetVector.at(i)).c_str(),  asOfDate,settleDate,payCal));
	}

	vector<double> strikes;
	//mid = jni->GetMethodID(cls, "getStrike", "(I)D");
	//set <double> :: iterator strike_Iter;
	//strike_Iter= strikeSet.begin();
	for (i = 0; i < strikeCount; i++) {
		//jdouble strike = jni->CallDoubleMethod(surface, mid, i);
		//strikes.push_back(strike*100.);
		//strikes.push_back((*strike_Iter++)*100.);
        strikes.push_back((strikeSetVector.at(i))*100.);
	}

	ARM_Matrix* volatilities = new ARM_Matrix(expDateCount, strikeCount);
	//mid = jni->GetMethodID(cls, "getVol", "(I)D");
	for (i = 0; i < expDateCount; i++)
	{
		for (int j = 0; j < strikeCount; j++)
		{
			//jdouble vol = jni->CallDoubleMethod(surface, mid, j+i*strikeCount);
			//volatilities->Elt(i, j) = vol*100.;
			volatilities->Elt(i, j) = valueVector.at(j+i*strikeCount*tenorCount+tenorRank*strikeCount)*100.;
		}
	}

	ARM_Vector* armExpDates = new ARM_Vector(expDates);
	ARM_Vector* armStrikes = new ARM_Vector(strikes);

	ARM_VolCurve* armVolCurve = new ARM_VolLInterpol(	asOfDate,
														armExpDates,
														armStrikes,
														volatilities,
														K_STK_TYPE_PRICE,
														K_SMILE_VOL);

	armVolCurve->SetCurrencyUnit(armCurrency);

    for (i = 0; i < expDateCount; i++) {
        sprintf(armVolCurve->itsYearTermsX[i], "%s", (expirySetVector.at(i)).c_str());
    }
	
    
    for (i = 0; i < tenorCount; i++) {
      //  sprintf(armVolCurve->itsYearTermsY[i], "%s", (strikeSetVector.at(i)).c_str());
        
	}


	/*if(armCurrency) {
		delete armCurrency;
		armCurrency = NULL;
	}*/

	return armVolCurve;
}

ARM_VolCurve* GetVolSmileFromCalypso(ARM_Date asOfDate, string xmlInput){
	
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 	
	xmlDoc = XMLTools::LoadXML(xmlInput); 

	vector<string> expiryVector;
	vector<string> tenorVector;
	vector<double> strikeVector;
	vector<double> valueVector;
	set<string> expirySet;
	set<string> tenorSet;
	set<double> strikeSet;
    vector<string> expirySetVector;
	vector<string> tenorSetVector;
	vector<double> strikeSetVector;

	loadData(xmlDoc, &expiryVector, &tenorVector, &strikeVector, &valueVector,&expirySet, &tenorSet, &strikeSet, &expirySetVector, &tenorSetVector, &strikeSetVector);

	/*mid = jni->GetMethodID(cls, "getCurrency","()Ljava/lang/String;");	
	jobject ccy = jni->CallObjectMethod(surface, mid);
	const char* ccyStr = jni->GetStringUTFChars((jstring)ccy, 0); 	    
	ARM_Currency* armCurrency = new ARM_Currency(ccyStr);
	jni->ReleaseStringUTFChars((jstring)ccy, ccyStr); */
	
	string ccy;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Currency"),ccy);
	ARM_Currency* armCurrency = new ARM_Currency(ccy.c_str());
	
	//ARM_VolCurve* armVolCurve =GetVolSmileFromCalypso(asOfDate, armCurrency,valueVector ,expirySet,strikeSet,tenorSet.size());
	ARM_VolCurve* armVolCurve =GetVolSmileFromCalypso(asOfDate, armCurrency,valueVector ,expirySetVector,strikeSetVector,tenorSetVector.size());
	
	if(armCurrency) {
		delete armCurrency;
		armCurrency = NULL;
	}

	return armVolCurve;
	
}	

/**
GetVolCubeFromCalypso
**/
ARM_VolCube* GetSmiledVolCubeFromCalypso(ARM_Date asOfDate, string xmlInput){
	
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 	
	xmlDoc = XMLTools::LoadXML(xmlInput); 
	
	vector<string> expiryVector;
	vector<string> tenorVector;
	vector<double> strikeVector;
	vector<double> valueVector;
	set<string> expirySet;
	set<string> tenorSet;
	set<double> strikeSet;
    vector<string> expirySetVector;
	vector<string> tenorSetVector;
	vector<double> strikeSetVector;

	loadData(xmlDoc, &expiryVector, &tenorVector, &strikeVector, &valueVector,&expirySet, &tenorSet, &strikeSet, &expirySetVector, &tenorSetVector, &strikeSetVector);

	string ccy;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Currency"),ccy);
	ARM_Currency* armCurrency = new ARM_Currency(ccy.c_str());

    string sindex;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Index"),sindex);


	// Builds smiles and gets their ID
	/*jclass cls = jni->GetObjectClass(cube); 
	jmethodID mid = jni->GetMethodID(cls, "getDate", "()J");
	jlong juliandate = jni->CallLongMethod(cube, mid);
	ARM_Date asOfDate = ARM_Date(juliandate);*/

	long juliandate =  asOfDate.GetJulian();

	//mid = jni->GetMethodID(cls, "getSmileCount", "()I");
	//jint smileCount = jni->CallIntMethod(cube, mid);
	int smileCount = tenorSetVector.size();

	//set <string> :: iterator tenor_Iter;
	//tenor_Iter= tenorSet.begin();

	vector<double> tenors;
	ARM_VolCurve** smileTab = new ARM_VolCurve*[smileCount];
	//mid = jni->GetMethodID(cls, "getSmile", "(I)Lcom/ixis/arm/common/wrapper/marketdata/ARMVolSurfaceSmile;");
	for (int smileIndex=0; smileIndex<smileCount; smileIndex++)
	{
		// Gets smile object
		//jobject smileSurface = jni->CallObjectMethod(cube, mid, smileIndex);
		//smileTab[smileIndex] = convertCalypsoVolatilitySmile(jni, smileSurface);
		//smileTab[smileIndex]=GetVolSmileFromCalypso(asOfDate,armCurrency, valueVector ,expirySet,strikeSet,tenorSet.size(),smileIndex);
        smileTab[smileIndex]=GetVolSmileFromCalypso(asOfDate,armCurrency, valueVector ,expirySetVector,strikeSetVector,tenorSetVector.size(),smileIndex);

		// Gets smile tenor
		//jclass smileCls = jni->GetObjectClass(smileSurface); 
		//jmethodID gettenorMid = jni->GetMethodID(smileCls, "getTenor", "()Ljava/lang/String;");
		//jobject tenorJStr = jni->CallObjectMethod(smileSurface, gettenorMid);
		//const char* tenorStr = jni->GetStringUTFChars((jstring)tenorJStr, 0); 	    
		//double tenor = CalypsoUtil::transcodeFromStrMatuToDouble(tenorStr);
		//tenors.push_back(tenor);
		//tenors.push_back(FromStrMatuToDouble((*tenor_Iter++).c_str()));
		//tenors.push_back(FromCalypsoStrMatuToDouble((*tenor_Iter++)));
        tenors.push_back(FromCalypsoStrMatuToDouble(tenorSetVector.at(smileIndex)));

		//jni->ReleaseStringUTFChars((jstring)tenorJStr, tenorStr);
	}

	ARM_Vector* armTenors = new ARM_Vector(tenors);
	
	// Builds Vol Cube (no ATM vol is used because the ATM volatility is added to the smile values in Calypso)
	ARM_VolCube* armVolCube = new ARM_VolCube(smileTab, smileCount, armTenors);

	if(smileTab) {
		for (smileIndex=0; smileIndex<smileCount; smileIndex++)	{
			if(smileTab[smileIndex]) {
				delete smileTab[smileIndex];
				smileTab[smileIndex] = NULL;
			}
		}
		delete[] smileTab;
	}
	if(armTenors) {
		delete armTenors;
		armTenors = NULL;
	}

	if(armCurrency) {
		delete armCurrency;
		armCurrency = NULL;
	}

	return armVolCube;
}

//TODO Change with a correct name
ARM_VolCube* GetVolCubeFromCalypso(ARM_Date asOfDate,string cvName, string xmlInput, string suffix){

    IXMLDOMDocumentPtr xmlDoc; 	
	xmlDoc = XMLTools::LoadXML(xmlInput); 
    string svolType;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Generator"),svolType);

    //Generator Smiled 
    if(svolType =="Smiled"){
    	IXMLDOMNodeListPtr xmlNodeList = XMLTools::selectNodes(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Parameters/Parameter"); 
	    long nbPoint=0;	
	
    	xmlNodeList->get_length(&nbPoint); 

	    string sname;
        string svalue;
        string ATMVolName;
        string smileVolName;
	
        ARM_VolCube* volSmiled;
        ARM_VolCurve* volATM;
    	for(int i=0;i<nbPoint;i++) 
    	{
	    	XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Name"),sname);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Value"),svalue);
            if(sname=="Vol Surface") ATMVolName=svalue;
            if(sname=="Smile Cube") smileVolName=svalue;
        }
        	//ARM_CalypsoToolkit::GetVolatilitySurface(volSurfName, cvName, date,	xmlInput);
        ARM_CalypsoToolkit::GetVolatilitySurface(smileVolName, cvName, asOfDate,	xmlInput);
	    volSmiled =GetSmiledVolCubeFromCalypso( asOfDate,xmlInput);
	    
        ARM_CalypsoToolkit::GetVolatilitySurface(ATMVolName, cvName, asOfDate,	xmlInput);
        volATM = GetVolSurfaceFromCalypso( asOfDate, xmlInput);
	    
        volSmiled->SetATMVol(volATM);
	    volSmiled->SetVolType(K_ATMF_VOL);
	    volSmiled->SetATMref(true);
	    //volSmiled = auto_ptr<ARM_VolLInterpol>((ARM_VolLInterpol*)volSmiled->Clone());
        return volSmiled; 
    }
        //Generator SABR
       else if(svolType =="SABR"){
    	IXMLDOMNodeListPtr xmlNodeList = XMLTools::selectNodes(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Parameters/Parameter"); 
	    long nbPoint=0;	
	
    	xmlNodeList->get_length(&nbPoint); 

	    string sname;
        string svalue;
	
        string betaVolName;
        string rhoVolName;
        string nuVolName;
        string ATMVolName;

    	for(int i=0;i<nbPoint;i++) 
    	{
	    	XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Name"),sname);
            XMLTools::convert(XMLTools::selectSingleNode(XMLTools::get_item(xmlNodeList,i),"Value"),svalue);
            if(sname=="Beta Surface") betaVolName=svalue;
            if(sname=="Rho Surface") rhoVolName=svalue;
            if(sname=="Nu Surface") nuVolName=svalue;
            if(sname=="Vol Surface") ATMVolName=svalue;
        
        }
        //ARM_CalypsoToolkit::GetVolatilitySurface(ATMVolName, "MO", asOfDate,	xmlInput);
        return GetSmiledVolCubeFromCalypso( asOfDate,xmlInput);
    }

    //Other Case
    else{
        return GetSmiledVolCubeFromCalypso( asOfDate,xmlInput);
    }

}

/**
GetVolSurfaceFromCalypso
**/
ARM_VolLInterpol* GetVolSurfaceFromCalypso(ARM_Date asOfDate, string xmlInput){

	
	MSXML2::IXMLDOMDocumentPtr xmlDoc; 	
	xmlDoc = XMLTools::LoadXML(xmlInput); 

	vector<string> expiryVector;
	vector<string> tenorVector;
	vector<double> strikeVector;
	vector<double> valueVector;
	set<string> expirySet;
	set<string> tenorSet;
	set<double> strikeSet;
    vector<string> expirySetVector;
	vector<string> tenorSetVector;
	vector<double> strikeSetVector;

	loadData(xmlDoc, &expiryVector, &tenorVector, &strikeVector, &valueVector,&expirySet, &tenorSet, &strikeSet, &expirySetVector, &tenorSetVector, &strikeSetVector);
	
	/*jclass cls = jni->GetObjectClass(surface); 
	jmethodID mid = jni->GetMethodID(cls, "getDate", "()J");
	jlong juliandate = jni->CallLongMethod(surface, mid);
	ARM_Date asOfDate = ARM_Date(juliandate);*/

	long juliandate =  asOfDate.GetJulian();

	/*mid = jni->GetMethodID(cls, "getCurrency","()Ljava/lang/String;");	
	jobject ccy = jni->CallObjectMethod(surface, mid);
	const char* ccyStr = jni->GetStringUTFChars((jstring)ccy, 0); 	    
	ARM_Currency* armCurrency = new ARM_Currency(ccyStr);
	jni->ReleaseStringUTFChars((jstring)ccy, ccyStr);*/

	string ccy;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Currency"),ccy);
	ARM_Currency* armCurrency = new ARM_Currency(ccy.c_str());

    string sindex;XMLTools::convert(XMLTools::selectSingleNode(xmlDoc,"/ListeVolSurface/VolSurface/SurfaceInfo/Index"),sindex);

    //string calendar;XMLTools::convert(XMLTools::selectSingleNode(xmlNode,"/ListeVolSurface/VolSurface/SurfaceInfo/Holidays"),calendar);
	//calendar = calypso2ARMCalendar(calendar);

    ARM_Currency sCCY(ccy.c_str());
    long spotDays = sCCY.GetSpotDays();
    char* payCalTmp = sCCY.GetPayCalName(sCCY.GetVanillaIndexType());
    char payCal[30];
    strcpy(payCal, payCalTmp);
    delete payCalTmp;

    ARM_Date settleDate = asOfDate;
	settleDate.NextBusinessDay(spotDays, payCal);

	//mid = jni->GetMethodID(cls, "getExpiryDateCount", "()I");
	//jint expDateCount = jni->CallIntMethod(surface, mid);
	int expDateCount = expirySet.size();

	//mid = jni->GetMethodID(cls, "getTenorCount", "()I");
	//jint tenorCount = jni->CallIntMethod(surface, mid);
	int tenorCount = tenorSet.size();

	vector<double> expDates;
	//set <string> :: iterator expiry_Iter;
	//expiry_Iter= expirySet.begin();
	//mid = jni->GetMethodID(cls, "getExpiryDate", "(I)J");
	int i;
	for (i = 0; i < expDateCount; i++) {
		//jlong date = jni->CallLongMethod(surface, mid, i);
		//expDates.push_back((date-juliandate)/365.);
		//expDates.push_back((FromStrMatuToDouble((*expiry_Iter++).c_str())-juliandate)/365.);
        //expDates.push_back(convPlotInYearTerm((*expiry_Iter++).c_str(),  asOfDate,settleDate,payCal));
        expDates.push_back(convCalypsoPlotInYearTerm((expirySetVector.at(i)).c_str(),  asOfDate,settleDate,payCal));
    }
	vector<double> tenors;
	//mid = jni->GetMethodID(cls, "getTenor", "(I)Ljava/lang/String;");
	//set <string> :: iterator tenor_Iter;
	//tenor_Iter= tenorSet.begin();
	for (i = 0; i < tenorCount; i++) {
		//jobject tenorJStr = jni->CallObjectMethod(surface, mid, i);
		//const char* tenorStr = jni->GetStringUTFChars((jstring)tenorJStr, 0); 	    
		//double yearFraction = CalypsoUtil::transcodeFromStrMatuToDouble(tenorStr);
		//jni->ReleaseStringUTFChars((jstring)tenorJStr, tenorStr);
		//tenors.push_back(yearFraction);
		//tenors.push_back((FromStrMatuToDouble((*tenor_Iter++).c_str())));
        tenors.push_back((FromStrMatuToDouble((tenorSetVector.at(i)).c_str())));
	}
 
	ARM_Matrix* volatilities = new ARM_Matrix(expDateCount, tenorCount);
	//mid = jni->GetMethodID(cls, "getVol", "(I)D");
	for (i = 0; i < expDateCount; i++)
	{
		for (int j=0; j < tenorCount; j++)
		{
			//jdouble vol = jni->CallDoubleMethod(surface, mid, i*tenorCount + j);
			//volatilities->Elt(i, j) = vol*100.;
			volatilities->Elt(i, j) = valueVector.at(i*tenorCount + j)*100.;
		}
	}

	ARM_Vector* armExpDates = new ARM_Vector(expDates);
	ARM_Vector* armTenors = new ARM_Vector(tenors);

	ARM_VolLInterpol* armVolCurve = new ARM_VolLInterpol(	asOfDate,
														armExpDates,
														armTenors,
														volatilities,
														K_STK_TYPE_PRICE,
														K_ATMF_VOL);
    armVolCurve->SetIndexName((char *)sindex.c_str());

	armVolCurve->SetCurrencyUnit(armCurrency);

	delete armCurrency;

    //set <string> :: iterator expiry_Iter;
	//expiry_Iter= expirySet.begin();
	for (i = 0; i < expDateCount; i++) {
        sprintf(armVolCurve->itsYearTermsX[i], "%s", (expirySetVector.at(i)).c_str());
    }
	//set <string> :: iterator tenor_Iter;
	//tenor_Iter= tenorSet.begin();
	for (i = 0; i < tenorCount; i++) {
        sprintf(armVolCurve->itsYearTermsY[i], "%s", (tenorSetVector.at(i)).c_str());
        
	}
   
	return armVolCurve;
}


ARM_ReferenceValue* ARMLOCAL_ParseFixingFromInstrument(	const CCString& xmlResponse )
{
	HRESULT hr;
	VARIANT_BOOL bOK;
	MSXML2::IXMLDOMDocument *XMLDoc = NULL;

	try
	{
		hr = CoInitialize(NULL); 

		hr = CoCreateInstance(__uuidof(MSXML2::DOMDocument30), NULL, CLSCTX_INPROC_SERVER, __uuidof(MSXML2::IXMLDOMDocument), (void**)&XMLDoc);
		SUCCEEDED(hr) ? 0 : throw hr;

		_bstr_t tmpChaine; 
		VariantTools::convert(CCSTringToSTLString(xmlResponse),tmpChaine); 

		XMLDoc->loadXML(tmpChaine, &bOK);
	}

	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in ARMLOCAL_ParseFixingFromInstrument");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
	}

	ARM_ReferenceValue* refValue = NULL;
	vector<double> dates;
	vector<double> values;

	try
	{
		MSXML2::IXMLDOMNodeList * resultList = NULL, * resultList2 = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, * item = NULL;

		if (XMLDoc->selectNodes(_bstr_t((const char *)((CCString)"Response/EXOTIC/Assets/ASSET")), &resultList) == S_OK)
		{
			long nbNodes = 0;
			resultList->get_length(&nbNodes);

			for (int i=0; i<nbNodes; i++)
			{
				hr=resultList->get_item(i, &listItem);
				listItem->selectNodes(_bstr_t((const char *)"Events/EVENT"), &resultList2);
				long nbOpEvents;
				resultList2->get_length(&nbOpEvents);

				for (int i = 0; i < nbOpEvents; i++)
				{
					hr=resultList2->get_item(i, &item);

					item->selectSingleNode(_bstr_t((const char *)"Type"), &theNode);
					if (theNode!=NULL)
					{
						BSTR resultat = NULL;
						theNode->get_text(&resultat);
						_bstr_t ff(resultat,false);
						char * ff1=(char *)ff;

						if (resultat) SysFreeString(resultat);
						theNode->Release();
						theNode=NULL;

						if ( strcmp((const char*) ff1, "FRC") == 0)
						{
							ARM_Date tmpDate = GetDateFromXMLNode(item, "ADate");
							if (tmpDate != ARM_DEFAULT_DATE)
							{
								dates.push_back(GetDateFromXMLNode(item, "Date").GetJulian());
								values.push_back(100.0*GetDoubleFromXMLNode(item,"Amount"));
							}
						}
					}
				}

				if (dates.size())
				{
					refValue = new ARM_ReferenceValue(CreateARMVectorFromVECTOR(dates), CreateARMVectorFromVECTOR(values));
					refValue->SetCalcMethod(K_STEPUP_LEFT);
					break;
				}
			}
		}
	}
	catch (Exception& e)
	{
		throw e;
	}
	catch(...)
	{
		if (XMLDoc) XMLDoc->Release();
		hr = S_FALSE;

		CCString msg((CCString)"Pb in getting past fixings from instrument");

		throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
	}

	if (XMLDoc)
		XMLDoc->Release();

	return refValue;
}


/*------------------------------------------------------------------------------------*/
/*---- End Of File ----*/
