/*
 *
 * Copyright (c) CDC IXIS CM October 2004 Paris
 *
 * $Log: $
 *
 */

/*! \file AssetAssetParserManager.cpp
 *
 *  \brief
 *	\author  M. Campet
 *	\version 1.0
 *	\date October 2004
 */

#include "firsttoinc.h"

#include <glob\expt.h>
#include <inst\security.h>
#include <inst\fixleg.h>
#include <inst\spreadoption.h>
#include <inst\armdigital.h>
#include <util\fromto.h>
#include <inst\barrier.h>
#include <inst\stripoption.h>
#include <inst\stripdigitaloption.h>
#include <inst\forex.h>
#include <inst\globalcap.h>
#include <inst\corridorDblCondition.h>

#include "AssetParserManager.h"
#include "PaserManagerUtilities.h"
// GP Base
#include <gpbase\curvetypedef.h>
#include "gpbase\gpvector.h"
#include "gpbase\gpmatrix.h"
#include "gpbase\datestrip.h"
#include "gpbase\gptrigomatrix.h"
#include "gpbase\singleton.h"
#include "gpbase\curve.h"
#include "gpbase\curveconvert.h"
// GP Infra
#include "gpinfra\gensecurity.h"
#include "gpinfra\pricingadviser.h"
#include "gpinfra\pricingmodel.h"
#include "gpinfra\curvemodelparam.h"
#include "gpinfra\correlmatparam.h"
// GP Calib
#include "gpcalib\calibmethod.h"
#include "gpcalib\modelparamsfactory.h"
// GP Calculator
#include "gpcalculators\maturitycapcalculator.h"
#include "gpcalculators\tarncalculator.h"
#include "gpcalculators\captioncalculator.h"
#include "gpcalculators\fxvanillacalculator.h"

#include "gpbase/timer.h"

#include "gpmodels/argconvdefault.h"

using ARM::ARM_PricingModelType;
using ARM::ARM_BasketType;
using ARM::ARM_VanillaType;
using ARM::ARM_DigitType;
using ARM::ARM_FXVanillaCalculator;

#include <ARMKernel\crv\volflat.h>



/// to fix the redefinition warning
/// be aware that the order is important as
/// this has to be before #include "arm_local_parsexml_util.h"
#if defined( va_start)
	#undef va_start
#endif

#if defined( va_end)
	#undef va_end
#endif

#include "arm_local_parsexml_util.h"
#include "ARM_local_parsexml.h"
#include <vector>

using namespace etoolkit;

using ARM::ARM_MaturityCapCalculator;
using ARM::ARM_TARNCalculator;
using ARM::ARM_CaptionCalculator;
using ARM::ARM_GenSecurity;
using ARM::ARM_PricingModel;
using ARM::ARM_CalibMethod;
using ARM::ARM_Curve;
using ARM::std::vector<double>;
using ARM::ARM_GP_Matrix;
using ARM::ARM_ModelParamType;
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_StringVector;
using ARM::ARM_CurveModelParam;
using ARM::ARM_CorrelMatParam;
using ARM::ARM_DateStrip;
using ARM::TrigoMatrix;
using ARM::ARM_ModelParamFactory;
using ARM::RefValueToCurve;

using ARM::ARM_Timer;
// Global Declarations

bool AssetParserManager::initDone = false;
int AssetParserManager::isEtkOrSummit = 1;
map<string, AbstractAssetParser*> AssetParserManager::itsAssetParsers;


void AssetParserManager::Release(void) 
{
	if ( !AssetParserManager::initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "AssetParserManager::Release, AssetParserManager::Init not done");
	}

	map<string, AbstractAssetParser*>::iterator iter = itsAssetParsers.begin();
	for (; iter != itsAssetParsers.end(); iter++) 
	{
		delete iter->second;
	}
}



AbstractAssetParser* AssetParserManager::GetAssetParser(string type)
{
	if ( !AssetParserManager::initDone ) 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "AssetParserManager::GetAssetParser, AssetParserManager::Init not done");
	}

	map<string, AbstractAssetParser*>::iterator iter = AssetParserManager::itsAssetParsers.find(type);	
	if ( iter != AssetParserManager::itsAssetParsers.end() )
	{
		return iter->second;
	}
	else 
	{
		CCString msg("AssetParserManager::GetAssetParser, no loader for Instrument type ");
		msg += type.c_str();
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION, msg);
	}	
}



bool AssetParserManager::RegisterParser(string& key, AbstractAssetParser* loader)
{
	pair<string, AbstractAssetParser*> p(key, loader);
	pair<map<string, AbstractAssetParser*>::iterator, bool> r = itsAssetParsers.insert(p);

	return r.second;
}



ARM_Object* AssetParserManager::BuildAsset(string type,
										   MSXML2::IXMLDOMDocument* XMLDoc,
										   string& assetId,
										   const ARM_Date& date)
{
	ARM_Object* objSec = NULL;

	AbstractAssetParser* p = GetAssetParser(type);

	if ( p ) 
	{
		objSec = p->BuildAsset(XMLDoc, date, assetId, IsEtkOrSummit());
	}
	else 
	{
		throw Exception (__LINE__, __FILE__, ERR_INVALID_OPERATION,
			 "AssetParserManager::BuildInstrument,no loader for Instrument type %s", type);
	}

	return objSec;
}


AbstractAssetParser::~AbstractAssetParser() 
{
}



//! \brief loader for Spread Option
class SpreadOptionAssetLoader : public AbstractAssetParser
{
public:

	SpreadOptionAssetLoader()
	{
	}

   ~SpreadOptionAssetLoader()
	{
	}

	void CorrectNotional(ARM_ReferenceValue* initNotional,
						 ARM_SpreadOption* object)
	{
		ARM_Vector* flowEndDates = object->GetFlowEndDates();
		ARM_Vector* paymentDates = object->GetSwapLeg()->GetTheoPayDates();

		if ( initNotional->GetDiscreteValues()->GetSize() != 1 )
		{
			//Correction sur le dernier : on lui met la date ajustée
			initNotional->GetDiscreteDates()->Elt(initNotional->GetDiscreteDates()->GetSize()-1) = 
				flowEndDates->Elt(flowEndDates->GetSize()-1);

			int j = 0;

			for (int i = 0; i < initNotional->GetDiscreteDates()->GetSize(); i++)
			{
				while ( (j < flowEndDates->GetSize())
					&& (initNotional->GetDiscreteDates()->Elt(i) > flowEndDates->Elt(j)) )
					j++;
				
				initNotional->GetDiscreteDates()->Elt(i) = paymentDates->Elt(j);
			}
		}
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		// Paramètres du produit

		char ProductName[30];

		int capFloor;
		double PorS = 1.0;

		ARM_Date start;
		ARM_Date end;

		ARM_ReferenceValue* notional = NULL;
		ARM_Currency* discountCcy = NULL;
		int basis;
		
		ARM_ReferenceValue* strikesGet = NULL;
		ARM_ReferenceValue* strikes = NULL;

		int fwdRule = 0;

		int decompPricingFlag;
		int intRule;

		int stubRule;

		int freezeFixing = 0;
		
		ARM_SpreadOption* newSpreadOpt = NULL;

		int payFreq = K_DEF_FREQ;
		char payCal[4];
		int payTiming;
		int paygap;

		int isQuantoOrNot = 0;

		// Paramètres du spread

		ARM_Currency* Idx1Ccy = NULL;
		ARM_Currency* Idx2Ccy = NULL;
		ARM_Currency* Idx1CcyTmp = NULL;
		ARM_Currency* Idx2CcyTmp = NULL;

		int liborType1;
		int liborType2;

		ARM_IRIndex* Idx1 = NULL;
		ARM_IRIndex* Idx2 = NULL;

		double weight1;
		double weight2;

		double spread1 =  0.01;
		double spread2 = -0.01;

		int resetFreq = K_DEF_FREQ;
		char resetCal[4];
		int resetTiming;
		int resetgap;
		
		ARM_ReferenceValue* spreadFirstFixing = NULL;
		ARM_ReferenceValue* spreadSecondFixing = NULL;
		
		// Paramètres de la paiement leg dans le cas des spread options corridor

		ARM_Currency* Idx3Ccy = NULL;
		ARM_Currency* Idx3CcyTmp = NULL;
		
		int liborType3;

		ARM_IRIndex* Idx3 = NULL;

		double weight3;

		double payIdxSpread = 0;

		int payIdxResetTiming;
		char payIdxResetCal[4];
		int payIdxResetFreq;
		int payIdxResetGap;

		int adjStartDate = 1;

		ARM_ReferenceValue* payIdxPastFixings = NULL;	// When pay index is variable
		ARM_ReferenceValue* payFixedRate = NULL; // When pay index is fixed

		int payBasis;

		// Paramètres de pricing

		int slopeFlag=1;
		int cptStrikeMethod=1;				//cptStrikeMethod is not available in Summit. 
											//Still, we must set it and pass it to the leg constructor
											//Set value as 1 for new method

		// Paramètres utilitaires
	
		char cust[20];	
		
		char myRefDate[11];
		char myRefDate2[11];

		ARM_SwapLeg* payLeg = NULL; 
		
		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			
			GetFormula(node,ProductName);

			capFloor = GetCapFloor(node);
			PorS = GetPorS(node);

			start = GetStartDate(node);
			end = GetEndDate(node);
			
			discountCcy = GetCcy(node);
			basis = GetAssetDayCount(node);

			strikesGet = GetSpreadVariable(node);// strikes correspond with the spreads in the Summit for the Spreadoption

			fwdRule = GetFwdRule(node);
			decompPricingFlag = GetDecompFlag(node);
			intRule = GetIntRule(node);
			
			stubRule = GetStubRule(node);

			resetFreq = GetResetFreq(node);
			GetResetCalendar(node,resetCal);
			resetgap = GetResetGap(node);
			
			payFreq = GetPayFreq(node);
			payTiming = GetPayTiming(node);
			paygap = GetPayGap(node);
			GetPayCalendar(node,payCal);

			GetStubDate1(node,myRefDate);		// Date de départ de la 1ere période non rompue si rompu début
			GetStubDate2(node,myRefDate2);		// Date de départ de la dernière période si rompu fin

//			strcpy(myRefDate,"NULL");
//			strcpy(myRefDate2,"NULL");

			char resetRolldate[4];
			GetPayRollDate(node,resetRolldate);

			ARM_Date startStubDate;
			ARM_Date endStubDate;

			if (strcmp(myRefDate, "NULL") != 0)
				startStubDate = (ARM_Date)myRefDate;
			if (strcmp(myRefDate2, "NULL") != 0)
				endStubDate = (ARM_Date)myRefDate2;

			int tmpFreq = payFreq;				// Cas des index Zero Coupon
			if (payFreq == K_ZEROCOUPON)
				tmpFreq = resetFreq;

			DeduceRefDateAndStub(start,end,startStubDate,endStubDate,tmpFreq,atoi(resetRolldate),discountCcy->GetCcyName(),(char*)myRefDate,stubRule);
			
			spreadFirstFixing = GetRefValSpreadFirstFixing(node,date);
			spreadSecondFixing = GetRefValSpreadSecondFixing(node,date);

			GetCustom(node, cust);

			// Initialisation des currency : initialement on considère un cas non quanto
			Idx1CcyTmp = (ARM_Currency*) discountCcy->Clone();
			Idx2CcyTmp = (ARM_Currency*) discountCcy->Clone();
			Idx3CcyTmp = (ARM_Currency*) discountCcy->Clone();

			if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR") == 0)
			{
				// Certaines données à parser ne se situent pas au même endroit selon qu'on
				// considère une spread option corridor ou les autres types de spreadoption

				// On récupère les index et on écrase les currency initiales dans le cas d'un index quanto
				liborType1 = GetBLOBIdx(node,Idx1CcyTmp,1);
				liborType2 = GetBLOBIdx(node,Idx2CcyTmp,2);
				liborType3 = GetBLOBIdx(node,Idx3CcyTmp,3);

				// On recrée les currency correspondant à la devise de l'index
				Idx1Ccy = new ARM_Currency(Idx1CcyTmp->GetCcyName());
				Idx2Ccy = new ARM_Currency(Idx2CcyTmp->GetCcyName());
				Idx3Ccy = new ARM_Currency(Idx3CcyTmp->GetCcyName());
				
				weight1 = GetBLOBAmount(node,1);
				if (weight1==0.0)
					weight1 = 1.0; 
				weight2 = GetBLOBAmount(node,2);
				if (weight2==0.0)
					weight2 = 1.0;
				weight3 = GetBLOBAmount(node,3);
				if (weight3==0.0)
					weight3 = 1.0;

				spread1 = GetBLOBAmount(node,4)/100.;		// le epsilon dans le menu spread option du deal dans Summit
				if (spread1 == 0.0)
					spread1 = 0.01;

				if (spread1 < 0.0)
					spread1 = -spread1;
				spread2 = -spread1;

				slopeFlag = GetBLOBNum(node,1); // NONE 0, NO 1, YES 2
				if (slopeFlag > 0)
					slopeFlag -= 1;

				resetTiming = GetBLOBTiming(node,1);		// reset timing du spread 
				payIdxResetTiming = GetBLOBTiming(node,2);	// reset timing de la pay leg
				
				if (GetBLOBNum(node,3) == 1)				// Bond basis
					payBasis = K30_360;
				else
					payBasis = basis;

				freezeFixing = GetBLOBNum(node,4);

				GetPayIndexResetCal(node, payIdxResetCal);
				payIdxResetGap = GetBLOBResetGap(node, 1);
				
				payIdxSpread = GetBLOBAmount(node, 5);
				
				payIdxResetFreq = GetBLOBNum(node,2);

				switch (payIdxResetFreq)
				{
					case 3:
					{
						payIdxResetFreq = 4;
						break;
					}
					case 4:
					{
						payIdxResetFreq = 12;
						break;
					}
					case 5 :
					{
						payIdxResetFreq = 52;
						break;
					}
					case 6 :
					{
						payIdxResetFreq = 365;
						break;
					}
					case 7 :
					{
						payIdxResetFreq = 0;
						break;
					}
				}

				if ( IsFixedIndex( (ARM_INDEX_TYPE)liborType3 ) )
				{
                    if ( payFixedRate == NULL )
                       payFixedRate = new ARM_ReferenceValue();

					GetPayIdxPastFixings(node, date, liborType3, payFixedRate);
				}
				else
				{
                    if ( payIdxPastFixings == NULL )
                       payIdxPastFixings = new ARM_ReferenceValue();

					GetPayIdxPastFixings(node, date, liborType3, payIdxPastFixings);
				}

				payLeg = GetPayLeg(node,liborType3,Idx3Ccy,payBasis,payFreq,end,date);
			}
			else
			{
				liborType1 = GetSpreadIndexType1(node,Idx1CcyTmp);
				liborType2 = GetSpreadIndexType2(node,Idx2CcyTmp);

				Idx1Ccy = new ARM_Currency(Idx1CcyTmp->GetCcyName());
				Idx2Ccy = new ARM_Currency(Idx2CcyTmp->GetCcyName());

				weight1 = GetWeight1(node);
				if (weight1==0.0)
					weight1 = 1.0; 
				weight2 = GetWeight2(node);
				if (weight2==0.0)
					weight2 = 1.0;

				// NB : Dans le cas des spread options non corridor le epsilon n'est pas renseigné dans Summit
				//		Dans le cas des spread options non corridor le slopeflag n'est pas renseigné dans Summit
				
				resetTiming = GetResetTiming(node);
			}

			ARM_ReferenceValue vWeight1(weight1);
			ARM_ReferenceValue vWeight2(weight2);

			//# IDENTIFICATION DES CAS QUANTO //
			if (strcmp(Idx1Ccy->GetCcyName(), discountCcy->GetCcyName()) || strcmp(Idx2Ccy->GetCcyName(), discountCcy->GetCcyName()))
				isQuantoOrNot = 1;

			if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR")==0)
			{
				if (strcmp(Idx3Ccy->GetCcyName(), discountCcy->GetCcyName()))
					isQuantoOrNot = 1;
			}
			//#//

			if ( (isQuantoOrNot) && payLeg)
			{
				delete payLeg;
				payLeg = NULL;
			}

			if (!isQuantoOrNot)			// Cas non quanto
			{
				if (strcmp(cust,"CUST") != 0)
				{
					notional = GetNotional(node);

					if (strcmp(ProductName,"SPREADOPTIONLOG")==0)
					{					
						newSpreadOpt = new ARM_SpreadOption(start, end,
															 capFloor, (ARM_ReferenceValue*)NULL,
															 (ARM_INDEX_TYPE) liborType1,  
															 (ARM_INDEX_TYPE) liborType2,  
															 &vWeight1, &vWeight2,
															 basis, 
															 resetFreq, 
															 payFreq, 
															 resetTiming, 
															 payTiming, 
															 discountCcy,
															 resetgap,
															 spreadFirstFixing,
															 spreadSecondFixing,
															 intRule,
															 stubRule,
															 cptStrikeMethod,
															 resetCal,
															 payCal,
															 K_SPREADOPTION_FORMULA,
															 fwdRule,
															 myRefDate,
															 paygap);
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGDIGITAL")==0)
					{
						ARM_ReferenceValue* payoff = GetVariableSpreadPayOff(node);
						slopeFlag = GetSpreadSlopeFlag(node);

						newSpreadOpt= new ARM_SpreadOption(start, end,
															capFloor, (ARM_ReferenceValue*)NULL,
															payoff,
															(ARM_INDEX_TYPE) liborType1,  
															(ARM_INDEX_TYPE) liborType2,  
															&vWeight1, &vWeight2, 
															basis, 
															resetFreq, 
															payFreq, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spread1,
															spread2,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															resetCal,
															payCal,
															slopeFlag,
															cptStrikeMethod,
															K_SPREADOPTION_FORMULA,
															fwdRule,
															myRefDate,
															paygap);
						if (payoff)
							delete payoff;
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGFLTDIGITAL")==0)
					{
						int PayOffliborType = GetSpreadPayoffLibor(node,discountCcy);
						slopeFlag = GetSpreadSlopeFlag(node);
						newSpreadOpt= new ARM_SpreadOption(start, end,
														capFloor, (ARM_ReferenceValue*)NULL,
														(ARM_INDEX_TYPE) PayOffliborType, 
														(ARM_INDEX_TYPE) liborType1,
														(ARM_INDEX_TYPE) liborType2,
														&vWeight1, &vWeight2,
														basis, 
														resetFreq, 
														payFreq, 
														resetTiming, 
														payTiming, 
														discountCcy,
														resetgap,
														spread1,
														spread2,
														spreadFirstFixing,
														spreadSecondFixing,
														intRule,
														stubRule,
														resetCal,
														payCal,
														slopeFlag,
														cptStrikeMethod,
														NULL,
														NULL,
														K_SPREADOPTION_FORMULA,
														fwdRule,
														myRefDate,
														paygap);
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR")==0)
					{
						ARM_ReferenceValue vWeight3(weight3);

						// Optim pour Mercure
						if(resetFreq > GetResetFreqForCorridorOptim())
							resetFreq = GetResetFreqForCorridorOptim();

						// CONSTRUCTION DE LA PAY LEG

/*						ARM_SwapLeg* payLeg = NULL; 

						if (IsCMSIndex((ARM_INDEX_TYPE) liborType3))
						{    
						   payLeg = (ARM_SwapLeg *) new ARM_CMSLeg(	start,
																	end,
																	(ARM_INDEX_TYPE) liborType3,
																	PorS,
																	payIdxSpread,
																	decompPricingFlag,
																	payBasis,
																	intRule,
																	payIdxResetFreq,
																	payFreq, 
																	payIdxResetGap,
																	discountCcy,
																	payIdxResetTiming,
																	stubRule,
																	payIdxResetCal,
																	payCal,
																	adjStartDate,			// defaulted : AdjStartDate
																	fwdRule,
																	myRefDate,
																	paygap);
						   payLeg->SetIndexType((ARM_INDEX_TYPE) liborType3);
						}
						else
						{
						   if (IsLiborIndex((ARM_INDEX_TYPE) liborType3))
						   {       
							  payLeg = (ARM_SwapLeg *) new ARM_SwapLeg( start, end,
																		(ARM_INDEX_TYPE) liborType3,
																		PorS,
																		payIdxSpread, 
																		payIdxResetFreq, 
																		payFreq,
																		payIdxResetTiming, payTiming,
																		discountCcy,
																		intRule,
																		payIdxResetGap,
																		payIdxResetCal,
																		payCal,
																		decompPricingFlag,
																		K_NX_NONE,
																		stubRule,
																		myRefDate,
																		adjStartDate,		// defaulted : AdjStartDate
																		basis, 
																		fwdRule,
																		paygap);
						   }       
						   else 
						   {
							   if (IsFixedIndex((ARM_INDEX_TYPE) liborType3)) // si FIXED, recuperation du spread, considere comme le taux fixe
							   {
								  payLeg = (ARM_SwapLeg *) new ARM_SwapLeg( start, end,
																			payIdxSpread,
																			PorS,
																			payIdxResetFreq, 
																			basis,
																			payIdxResetFreq, // ?
																			payTiming, intRule,
																			stubRule, 
																			discountCcy);
								  payLeg->SetPayCalName(payCal);
								  payLeg->SetResetCalName(payIdxResetCal);

                                  if (payFixedRate)
                                  {
                                     payLeg->SetVariableSpread(payFixedRate);

                                     delete payFixedRate;

                                     payFixedRate = NULL;
                                  }
                                  
								  payLeg->CptCashFlowDates();
 
							   }
							   else if (IsCMTIndex((ARM_INDEX_TYPE) liborType3))
							   {
								   int resetGapVal = payIdxResetGap <= 0 ? payIdxResetGap : (-discountCcy->GetSpotDays());
               
								   payLeg = (ARM_SwapLeg *) new ARM_CMTLeg(	start,
																			end,
																			(ARM_INDEX_TYPE) liborType3,
																			discountCcy->GetFixedPayFreq(),
																			discountCcy->GetFixedDayCount(), 
																			PorS,
																			payIdxSpread,
																			payIdxResetFreq, // ?,
																			basis,
																			intRule,
																			resetGapVal,
																			payIdxResetFreq,
																			100.0,			// ??
																			discountCcy,
																			payIdxResetTiming,
																			stubRule);

								   payLeg->SetIndexType((ARM_INDEX_TYPE) liborType3);
							   }
						   }
						}
*/
						// CONSTRUCTION DE LA SPREAD OPTION
/*						newSpreadOpt = new ARM_SpreadOption(start,
															end,
															capFloor,
															(ARM_ReferenceValue*)NULL,
															payLeg,
															(ARM_INDEX_TYPE) liborType1,
															(ARM_INDEX_TYPE) liborType2,
															&vWeight1, &vWeight2,
															payLeg->GetDayCount(),
															resetFreq,
															payFreq,
															resetTiming,
															payTiming,
															discountCcy,
															resetgap,
															spread1,
															spread2,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															resetCal,
															payCal,
															slopeFlag,
															cptStrikeMethod,
															weight3,
															freezeFixing,
															K_SPREADOPTION_FORMULA,
															fwdRule,
															myRefDate,
															paygap);
*/

							newSpreadOpt = new ARM_SpreadOption(start,
																end,
																capFloor,
																(ARM_ReferenceValue*)NULL,
																payLeg,//(ARM_INDEX_TYPE)payoffliborType1Id, 
																(ARM_INDEX_TYPE) liborType1,
																(ARM_INDEX_TYPE) liborType2,
																&vWeight1, &vWeight2,
																payLeg->GetDayCount(),
																resetFreq,
																payFreq,
																resetTiming,
																payTiming,
																discountCcy,
																resetgap,
																spread1,
																spread2,
																spreadFirstFixing,
																spreadSecondFixing,
																intRule,
																stubRule,
																resetCal,
																payCal,
																slopeFlag,
																cptStrikeMethod,
																weight3,
																freezeFixing,
																K_SPREADOPTION_FORMULA,
																fwdRule,
																myRefDate,
																paygap);

							if (payLeg)
								delete payLeg;
							payLeg = NULL;
					}
					else
					{
						CCString msg((CCString)"Invalid Product Formula for getting SpreadOption \n" + ProductName);
						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
					}
				}
				else
				{
					ARM_SwapLeg* firstLeg = NULL;
					ARM_SwapLeg* secondLeg = NULL;
					ARM_Vector* flowStartDates = NULL;
					ARM_Vector* flowEndDates = NULL;
					ARM_Vector* resetDates = NULL;
					ARM_Vector* paymentDates = NULL;
					ARM_Vector* interestDays = NULL;

					GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

					interestDays = new ARM_Vector(flowStartDates->GetSize());
					for (int i = 0; i < flowStartDates->GetSize(); i++)
						interestDays->Elt(i) = DaysBetweenDates(basis, flowStartDates->Elt(i), flowEndDates->Elt(i));

					ARM_SecurityFlows* custEcheancier = new ARM_SecurityFlows();
					custEcheancier->SetValues(string("STARTDATES"),flowStartDates);
					custEcheancier->SetValues(string("ENDDATES"),flowEndDates);
					custEcheancier->SetValues(string("RESETDATES"),resetDates);
					custEcheancier->SetValues(string("PAYMENTDATES"),paymentDates);
					custEcheancier->SetValues(string("INTERESTDAYS"),interestDays);

					notional = GetNotionalCust(node);

					if (IsCMSIndex((ARM_INDEX_TYPE)liborType1))
					{
						firstLeg = (ARM_SwapLeg *) new ARM_CMSLeg(flowStartDates,
																  flowEndDates,
																  resetDates,
																  paymentDates,
																  interestDays,
																  (ARM_INDEX_TYPE)liborType1,
																  notional,
																  K_RCV,
																  0.0,
																  K_COMP_PROP,
																  basis,
																  intRule,
																  resetFreq,
																  payFreq,
																  10000,
																  discountCcy,
																  resetTiming,
																  stubRule,
																  resetCal,
																  payCal);
						firstLeg->SetIndexType((ARM_INDEX_TYPE)liborType1);
					}
					else
					{
						ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)liborType1,resetFreq,payFreq,discountCcy);
						
						firstLeg = new ARM_SwapLeg(flowStartDates,
												   flowEndDates,
												   paymentDates,
												   resetDates,
												   interestDays,
												   NULL,
												   notional,
												   newIndex,
												   K_RCV,
												   0.0,
												   K_FLOAT_RATE,
												   discountCcy,
												   K_NX_NONE,
												   decompPricingFlag,
												   basis);

						firstLeg->SetStartDate(start);
						firstLeg->SetEndDate(end);

	//					firstLeg->CptRealCashFlowDates();

						delete newIndex;
					}    

					if (IsCMSIndex((ARM_INDEX_TYPE)liborType2))
					{
						secondLeg = (ARM_SwapLeg *) new ARM_CMSLeg(flowStartDates,
																   flowEndDates,
																   resetDates,
																   paymentDates,
																   interestDays,
																   (ARM_INDEX_TYPE)liborType2,
																   notional,
																   K_RCV,
																   0.0,
																   K_COMP_PROP,
																   basis,
																   intRule,
																   resetFreq,
																   payFreq,
																   10000,
																   discountCcy,
																   resetTiming,
																   stubRule,
																   resetCal,
																   payCal);
						secondLeg->SetIndexType((ARM_INDEX_TYPE)liborType2);
					}
					else
					{
						ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)liborType2,resetFreq,payFreq,discountCcy);
						
						secondLeg = new ARM_SwapLeg(flowStartDates,
													flowEndDates,
													paymentDates,
													resetDates,
													interestDays,
													NULL,
													notional,
													newIndex,
													K_RCV,
													0.0,
													K_FLOAT_RATE,
													discountCcy,
													K_NX_NONE,
													decompPricingFlag,
													basis);

						secondLeg->SetStartDate(start);
						secondLeg->SetEndDate(end);

						delete newIndex;
					}    

					if (strcmp(ProductName,"SPREADOPTIONLOG")==0)
					{
	/*					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   &vWeight1, &vWeight2,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   cptStrikeMethod);
	*/
						newSpreadOpt= new ARM_SpreadOption(start,end,custEcheancier,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   (ARM_INDEX_TYPE)liborType1,(ARM_INDEX_TYPE)liborType2,
														   &vWeight1, &vWeight2, notional,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   basis,resetFreq,payFreq,resetTiming,discountCcy,
														   intRule,stubRule,
														   cptStrikeMethod,resetCal,payCal);
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGDIGITAL")==0)
					{
						
						ARM_ReferenceValue* payoff = GetVariableSpreadPayOff(node);
						slopeFlag = GetSpreadSlopeFlag(node);
	/*					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   payoff,
														   &vWeight1, &vWeight2,
														   spread1,
														   spread2,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   slopeFlag,
														   cptStrikeMethod);
	*/
						newSpreadOpt= new ARM_SpreadOption(start,end,custEcheancier,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   payoff,
														   (ARM_INDEX_TYPE)liborType1,(ARM_INDEX_TYPE)liborType2,
														   &vWeight1, &vWeight2,
														   spread1, spread2, notional,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   basis,resetFreq,payFreq,resetTiming,discountCcy,
														   intRule,stubRule,
														   cptStrikeMethod,resetCal,payCal,
														   slopeFlag);

						if (payoff)
							delete payoff;
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGFLTDIGITAL")==0)
					{
						int PayOffliborType = GetSpreadPayoffLibor(node,discountCcy);
						slopeFlag = GetSpreadSlopeFlag(node);

	/*					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   (ARM_INDEX_TYPE) PayOffliborType,
														   &vWeight1, &vWeight2,
														   spread1,
														   spread2,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   slopeFlag,
														   cptStrikeMethod);
	*/
					
						newSpreadOpt= new ARM_SpreadOption(start,end,custEcheancier,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   (ARM_INDEX_TYPE)PayOffliborType,
														   (ARM_INDEX_TYPE)liborType1,(ARM_INDEX_TYPE)liborType2,
														   &vWeight1, &vWeight2,
														   spread1, spread2, notional,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   NULL, NULL,
														   basis,resetFreq,payFreq,resetTiming,discountCcy,
														   intRule,stubRule,
														   cptStrikeMethod,resetCal,payCal,
														   slopeFlag);
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR")==0)
					{
						// POURQUOI IL N'EXISTE PAS DE CONSTRUCTEURS AVEC ECHEANCIER POUR LE CORRIDOR ???

						// CONSTRUCTION DE LA PAY LEG

						ARM_SwapLeg* payLeg = NULL; 

						if (IsCMSIndex((ARM_INDEX_TYPE) liborType3))
						{    
						   payLeg = (ARM_SwapLeg *) new ARM_CMSLeg(	start,
																	end,
																	(ARM_INDEX_TYPE) liborType3,
																	PorS,
																	payIdxSpread,
																	decompPricingFlag,
																	payBasis,
																	intRule,
																	payIdxResetFreq,
																	payFreq, 
																	payIdxResetGap,
																	discountCcy,
																	payIdxResetTiming,
																	stubRule,
																	payIdxResetCal,
																	payCal,
																	adjStartDate,			// defaulted : AdjStartDate
																	fwdRule,
																	myRefDate,
																	paygap);
						   payLeg->SetIndexType((ARM_INDEX_TYPE) liborType3);
						}
						else
						{
						   if (IsLiborIndex((ARM_INDEX_TYPE) liborType3))
						   {       
							  payLeg = (ARM_SwapLeg *) new ARM_SwapLeg( start, end,
																		(ARM_INDEX_TYPE) liborType3,
																		PorS,
																		payIdxSpread, 
																		payIdxResetFreq, 
																		payFreq,
																		payIdxResetTiming, payTiming,
																		discountCcy,
																		intRule,
																		payIdxResetGap,
																		payIdxResetCal,
																		payCal,
																		decompPricingFlag,
																		K_NX_NONE,
																		stubRule,
																		myRefDate,
																		adjStartDate,		// defaulted : AdjStartDate
																		basis, 
																		fwdRule,
																		paygap);
						   }       
						   else 
						   {
							   if (IsFixedIndex((ARM_INDEX_TYPE) liborType3)) // si FIXED, recuperation du spread, considere comme le taux fixe
							   {
								  payLeg = (ARM_SwapLeg *) new ARM_SwapLeg( start, end,
																			payIdxSpread,
																			PorS,
																			payIdxResetFreq, 
																			basis,
																			payIdxResetFreq, // ?
																			payTiming, intRule,
																			stubRule, 
																			discountCcy);
								  payLeg->SetPayCalName(payCal);
								  payLeg->SetResetCalName(payIdxResetCal);

                                  if (payFixedRate)
                                  {
                                     payLeg->SetVariableSpread(payFixedRate);

                                     delete payFixedRate;

                                     payFixedRate = NULL;
                                  }

								  payLeg->CptCashFlowDates();
 
							   }
							   else if (IsCMTIndex((ARM_INDEX_TYPE) liborType3))
							   {
								   int resetGapVal = payIdxResetGap <= 0 ? payIdxResetGap : (-discountCcy->GetSpotDays());
               
								   payLeg = (ARM_SwapLeg *) new ARM_CMTLeg(	start,
																			end,
																			(ARM_INDEX_TYPE) liborType3,
																			discountCcy->GetFixedPayFreq(),
																			discountCcy->GetFixedDayCount(), 
																			PorS,
																			payIdxSpread,
																			payIdxResetFreq, // ?,
																			basis,
																			intRule,
																			resetGapVal,
																			payIdxResetFreq,
																			100.0,			// ??
																			discountCcy,
																			payIdxResetTiming,
																			stubRule);

								   payLeg->SetIndexType((ARM_INDEX_TYPE) liborType3);
							   }
						   }
						}

						// CONSTRUCTION DE LA SPREAD OPTION

						newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
														   capFloor, (ARM_ReferenceValue*)NULL,
														   payLeg,
														   &vWeight1, &vWeight2,
														   spread1,
														   spread2,
														   spreadFirstFixing,
														   spreadSecondFixing,
														   slopeFlag,
														   cptStrikeMethod,
														   weight3,
														   freezeFixing);

						if (payLeg)
							delete payLeg;
						payLeg = NULL;
					}
					else
					{
						CCString msg((CCString)"Invalid Product Formula for getting SpreadOption \n" + ProductName);
						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
					}

					newSpreadOpt->GetSwapLeg()->CptTheoPayDates();
					newSpreadOpt->SetPaymentDates((ARM_Vector*)newSpreadOpt->GetSwapLeg()->GetTheoPayDates()->Clone());

					if (flowStartDates)
						delete flowStartDates;
					flowStartDates = NULL;

					if (flowEndDates)
						delete flowEndDates;
					flowEndDates = NULL;

					if (resetDates)
						delete resetDates;
					resetDates = NULL;

					if (paymentDates)
						delete paymentDates;
					paymentDates = NULL;

					if (interestDays)
						delete interestDays;
					interestDays = NULL;

					if (firstLeg)
						delete firstLeg;
					firstLeg = NULL;

					if (secondLeg)
						delete secondLeg;
					secondLeg = NULL;

					if (custEcheancier)
						delete custEcheancier;
					custEcheancier= NULL;

				}
			}
			else		// Cas quanto
			{	
				/*// CONSTRUCTION DES INDEX DU SPREAD //*/

				Idx1 = new ARM_IRIndex((ARM_INDEX_TYPE) liborType1, resetFreq, payFreq, Idx1Ccy, Idx1Ccy->GetLiborIndexDayCount());
				Idx2 = new ARM_IRIndex((ARM_INDEX_TYPE) liborType2, resetFreq, payFreq, Idx2Ccy, Idx2Ccy->GetLiborIndexDayCount());

				if ( !(( liborType1 >= K_CMS1) && ( liborType1 <= K_CMS30 ))	// Idx1 est Libor
					&& Idx1->IsLiborIndex() )
				{

					Idx1->SetResetTiming(resetTiming);
					Idx1->SetPayTiming(payTiming);

					if (resetgap != 10000)
						Idx1->SetResetGap(resetgap);

					if (paygap != 10000)
						Idx1->SetPayGap(paygap);

					Idx1->SetIndexStyle(IN_ARREARS);

					if ( (Idx1->GetResetTiming() == K_ADVANCE)
						&& (Idx1->GetPayTiming() == K_ARREARS)
						&& (Idx1->GetResetFrequency() == Idx1->GetPayFrequency())		// Cas Index1 Vanille
					   )
					{
						Idx1->SetIndexStyle(VANILLA);
					}

					Idx1->SetIntRule(intRule);

				}
				else		// Cas Idx1 est CMS
				{
					ARM_INDEX_TYPE Libor_Type1 = Idx1Ccy->GetVanillaIndexType();

					Idx1->Set((ARM_INDEX_TYPE) liborType1, Libor_Type1, resetFreq, -1, Idx1Ccy);	// defaulted : cmCompMeth

					Idx1->Set(Idx1Ccy->GetLiborIndexDayCount(), resetFreq, payFreq,
								-1, -1, fwdRule, intRule,					// defaulted : term, compMethod
								resetTiming, resetgap, payTiming, paygap,
								Idx1Ccy, (ARM_INDEX_TYPE) liborType1, -1);	// defaulted : decompFreq
				}

				if ( !(( liborType2 >= K_CMS1) && ( liborType2 <= K_CMS30 ))		// Idx2 est Libor
					&& Idx2->IsLiborIndex() )
				{
					Idx2->SetResetTiming(resetTiming);
					Idx2->SetPayTiming(payTiming);

					if (resetgap != 10000)
						Idx2->SetResetGap(resetgap);

					if (paygap != 10000)
						Idx2->SetPayGap(paygap);

					Idx2->SetIndexStyle(IN_ARREARS);

					if ( (Idx2->GetResetTiming() == K_ADVANCE)						// Cas Index2 Vanille
						&& (Idx2->GetPayTiming() == K_ARREARS)
						&& (Idx2->GetResetFrequency() == Idx2->GetPayFrequency())
					   )
					{
						Idx2->SetIndexStyle(VANILLA);
					}

					Idx2->SetIntRule(intRule);

				}
				else
				{

					ARM_INDEX_TYPE Libor_Type2 = Idx2Ccy->GetVanillaIndexType();

					Idx2->Set((ARM_INDEX_TYPE) liborType2, Libor_Type2, resetFreq, -1, Idx2Ccy);

					Idx2->Set(Idx2Ccy->GetLiborIndexDayCount(), resetFreq, payFreq,
						-1, -1, fwdRule, intRule,
						resetTiming, resetgap, payTiming, paygap,
						Idx2Ccy, (ARM_INDEX_TYPE) liborType2, -1);
				}

				if (strcmp(cust,"CUST") != 0)
				{
					notional = GetNotional(node);

					/*//CONSTRUCTION DE LA SPREAD OPTION //*/

					if (strcmp(ProductName,"SPREADOPTIONLOG")==0)
					{					
						newSpreadOpt = new ARM_SpreadOption(start, end, 
															capFloor, strikesGet,
															Idx1, Idx2,
															&vWeight1, &vWeight2,
															basis, 
															resetFreq, 
															payFreq, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															cptStrikeMethod,
															resetCal,
															payCal,
															K_SPREADOPTION_FORMULA,			// defaulted : computedFormula
															fwdRule,
															myRefDate,
															paygap, (ARM_Vector*) NULL,		// defaulted : calibInfos
															notional);


					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGDIGITAL")==0)
					{
						ARM_ReferenceValue* payoff = GetVariableSpreadPayOff(node);
						slopeFlag = GetSpreadSlopeFlag(node);

						newSpreadOpt = new ARM_SpreadOption(start, end, 
															capFloor, strikesGet,
															payoff,
															Idx1, Idx2,
															&vWeight1, &vWeight2,
															basis, 
															resetFreq, 
															payFreq, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spread1,
															spread2,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															cptStrikeMethod,
															resetCal,
															payCal,
															slopeFlag,
															K_SPREADOPTION_FORMULA,			// defaulted : computedFormula
															fwdRule,
															myRefDate,
															paygap, (ARM_Vector*) NULL,		// defaulted : calibInfos
															notional);
						
						if (payoff)
							delete payoff;
						payoff = NULL;

					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGFLTDIGITAL")==0)
					{
						// CONSTRUCTION DE LA LEG DE PAIEMENT
						int PayOffliborType = GetSpreadPayoffLibor(node,Idx3CcyTmp);
						
						if ( (PayOffliborType == Idx1->GetIndexType()) && ( strcmp(Idx3CcyTmp->GetCcyName(), Idx1->GetCurrencyUnit()->GetCcyName() ) == 0 ) )  // La paiement leg corrrespond à la leg1 du spread
							Idx3 = (ARM_IRIndex *) Idx1->Clone();
						else if ( (PayOffliborType == Idx2->GetIndexType()) && ( strcmp(Idx3CcyTmp->GetCcyName(), Idx2->GetCurrencyUnit()->GetCcyName() ) == 0 ) )  // La paiement leg corrrespond à la leg1 du spread		// La paiement leg corrrespond à la leg1 du spread
							Idx3 = (ARM_IRIndex *) Idx2->Clone();
						else	// Problème...
						{
							CCString msg((CCString)"Pb on paiement Index when parsing Digitale Floater SpreadOption \n");
							throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
						}	

						slopeFlag = GetSpreadSlopeFlag(node);

						newSpreadOpt = new ARM_SpreadOption(start, end, 
															capFloor, strikesGet,
															Idx3,
															Idx1, Idx2,
															&vWeight1, &vWeight2,
															basis, 
															resetFreq, 
															payFreq, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spread1,
															spread2,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															resetCal,
															payCal,
															slopeFlag,
															cptStrikeMethod,
															(ARM_ReferenceValue*) NULL,		// defaulted : payMargins
															payIdxPastFixings,
															K_SPREADOPTION_FORMULA,			// defaulted : computedFormula
															fwdRule,
															myRefDate,
															paygap, (ARM_Vector*) NULL,		// defaulted : calibInfos
															notional);
					}
					else if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR")==0)
					{
						ARM_ReferenceValue vWeight3(weight3);
			
						if (IsFixedIndex((ARM_INDEX_TYPE) liborType3))
						{	
							Idx3 = GetPayIndex(node, (ARM_INDEX_TYPE) liborType3, Idx3Ccy, payBasis, payFreq);
						}
						else 
						{
							Idx3 = new ARM_IRIndex((ARM_INDEX_TYPE) liborType3, payIdxResetFreq, payFreq, Idx3Ccy, payBasis);
							
							if ( !(( liborType3 >= K_CMS1) && ( liborType3 <= K_CMS30 ))	// Idx3 est Libor
								&& Idx3->IsLiborIndex() )
							{

								Idx3->SetResetTiming(payIdxResetTiming);
								Idx3->SetPayTiming(payTiming);

								if (payIdxResetGap != 10000)
									Idx3->SetResetGap(payIdxResetGap);

								if (paygap != 10000)
									Idx3->SetPayGap(paygap);

								Idx3->SetIndexStyle(IN_ARREARS);

								if ( (Idx3->GetResetTiming() == K_ADVANCE)
									&& (Idx3->GetPayTiming() == K_ARREARS)
									&& (Idx3->GetResetFrequency() == Idx3->GetPayFrequency())		// Cas Index3 Vanille
								   )
								{
									Idx3->SetIndexStyle(VANILLA);
								}

								Idx3->SetIntRule(intRule);
							}
							else
							{

								ARM_INDEX_TYPE Libor_Type3 = Idx3Ccy->GetVanillaIndexType();

								Idx3->Set((ARM_INDEX_TYPE) liborType3, Libor_Type3, payIdxResetFreq, -1, Idx3Ccy);

								Idx3->Set(Idx3Ccy->GetLiborIndexDayCount(), payIdxResetFreq, payFreq,
									-1, -1, fwdRule, intRule,
									payIdxResetTiming, payIdxResetGap, payTiming, paygap,
									Idx3Ccy, (ARM_INDEX_TYPE) liborType3, -1);
							}
						}

						newSpreadOpt = new ARM_SpreadOption(start, end, 
															capFloor, (ARM_ReferenceValue*)NULL,
															Idx1, Idx2, Idx3, payFixedRate,
															&vWeight1, &vWeight2, &vWeight3,
															discountCcy, 
															spreadFirstFixing,
															spreadSecondFixing,
															payIdxPastFixings,
															spread1,
															spread2,
															slopeFlag,
															cptStrikeMethod,
															K_SPREADOPTION_FORMULA,
															basis, 
															resetFreq, 
															payFreq, 
															resetTiming, 
															payTiming, 
															intRule,
															stubRule,
															resetgap,
															resetCal,
															payCal,
															payIdxResetCal,
															fwdRule,
															myRefDate,
															notional,
															freezeFixing);
					}
					else
					{
						CCString msg((CCString)"Invalid Product Formula for getting SpreadOption \n" + ProductName);
						throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
									 (char*) msg);
					}

					if (Idx1)
						delete Idx1;
					Idx1 = NULL;

					if (Idx2)
						delete Idx2;
					Idx2 = NULL;

					if (Idx3)
						delete Idx3;
					Idx3 = NULL;

					if (Idx1Ccy)
						delete Idx1Ccy;
					Idx1Ccy = NULL;

					if (Idx2Ccy)
						delete Idx2Ccy;
					Idx2Ccy = NULL;

					if (Idx3Ccy)
						delete Idx3Ccy;
					Idx3Ccy = NULL;
				}
				else
				{
					// CUST NON IMPLEMENTE POUR LES CAS QUANTO
					CCString msg((CCString)"Cust not implemented for quanto sprea doptions \n");
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
				}
			}
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : SpreadOptionAssetLoader::BuildAsset()");
		} 

		// Gestion du strike
		if(strikesGet == NULL)
			strikes = new ARM_ReferenceValue(GetSpread(node)*100.0);
		else
		{
			ARM_SwapLeg* refLeg;

			if (newSpreadOpt->IsCorridorSpreadOption())
				refLeg = newSpreadOpt->GetPayIndexLeg();
			else
				refLeg = newSpreadOpt->GetSpreadLeg()->GetFirstLeg();

			strikes  = ARM_FromStartRefToResetRef(strikesGet,refLeg);
		}

		strikes->SetCalcMethod(K_STEPUP_LEFT);

		newSpreadOpt->SetStrikes(strikes);

		CorrectNotional(notional,newSpreadOpt);

		// Recup des Primes et Fees
		ARM_ReferenceValue* tmpPrimes = GetPrimes(node);
		if (tmpPrimes != NULL)
		{
			newSpreadOpt->SetFee(tmpPrimes);
			delete tmpPrimes;
		}

		newSpreadOpt->SetAmount(notional);
		newSpreadOpt->SetPorS(PorS);
		newSpreadOpt->SetAssetId(assetId);

		if (Idx1CcyTmp)
			delete Idx1CcyTmp;
		Idx1CcyTmp = NULL;

		if (Idx2CcyTmp)
			delete Idx2CcyTmp;
		Idx2CcyTmp = NULL;
		
		if (Idx3CcyTmp)
			delete Idx3CcyTmp;
		Idx3CcyTmp = NULL;
		
		if (Idx1)
			delete Idx1;
		Idx1 = NULL;

		if (Idx2)
			delete Idx2;
		Idx2 = NULL;

		if (Idx3)
			delete Idx3;
		Idx3 = NULL;

		if (Idx1Ccy)
			delete Idx1Ccy;
		Idx1Ccy = NULL;

		if (Idx2Ccy)
			delete Idx2Ccy;
		Idx2Ccy = NULL;

		if (Idx3Ccy)
			delete Idx3Ccy;
		Idx3Ccy = NULL;

		if (discountCcy)
			delete discountCcy;
		discountCcy = NULL;

		if (notional)
			delete notional;
		notional = NULL;

		if (spreadFirstFixing)
			delete spreadFirstFixing;
		spreadFirstFixing = NULL;

		if (spreadSecondFixing)
			delete spreadSecondFixing;
		spreadFirstFixing = NULL;

		if (payIdxPastFixings)
			delete payIdxPastFixings;
		payIdxPastFixings = NULL;

		if (payFixedRate)
			delete payFixedRate;
		payFixedRate = NULL;

		if (strikesGet)
			delete strikesGet;
		strikesGet = NULL;

		if (strikes)
			delete strikes;
		strikes = NULL;

		return newSpreadOpt;
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate,double notional,int isEtkOrNot = 1)
	{
		ARM_Date start;
		ARM_Date end;
		double spread=0.0;
		int capFloor;
		int resetTiming;
		int payTiming;
		int liborType1;
		int liborType2;
		int resetgap;
		int paygap;
		int resetFreq = K_DEF_FREQ;
		int basis;
		char ProductName[30];
		double weight1;
		double weight2;
		double spread1 =  0.01;
		double spread2 = -0.01;
		int intRule;
		double PorS = 1.0;
		int stubRule;
		char cust[20];
		char resetCal[4];
		char payCal[4];
		int slopeFlag=1;
		int freezeFixing = 0;
		int fwdRule = 0;
		int decompPricingFlag;

		char myRefDate[11];
		char myRefDate2[11];

		//cptStrikeMethod is not available in Summit. 
		//Still, we must set it and pass it to the leg constructor
		//Set value as 1 for new method
		int cptStrikeMethod=1;

		ARM_ReferenceValue* strikesGet = NULL;
		ARM_ReferenceValue* strikes = NULL;
		ARM_Currency* discountCcy = NULL;

		ARM_Currency* Idx1Ccy = NULL;
		ARM_Currency* Idx2Ccy = NULL;
		ARM_Currency* Idx3Ccy = NULL;
		ARM_Currency* Idx1CcyTmp = NULL;
		ARM_Currency* Idx2CcyTmp = NULL;
		ARM_Currency* Idx3CcyTmp = NULL;

		ARM_IRIndex* Idx1 =NULL;
		ARM_IRIndex* Idx2 = NULL;
		ARM_IRIndex* Idx3 = NULL;

		ARM_ReferenceValue* cstNotional = new ARM_ReferenceValue(notional);

		ARM_ReferenceValue* spreadFirstFixing = NULL;
		ARM_ReferenceValue* spreadSecondFixing = NULL;
		ARM_SwapLeg* payLeg = NULL;

		ARM_SpreadOption* newSpreadOpt = NULL;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			start = GetStartDate(node);
			end = GetEndDate(node);
			capFloor = GetCapFloor(node);
			PorS = GetPorS(node);

			discountCcy = GetCcy(node);
			basis = GetAssetDayCount(node);
			fwdRule = GetFwdRule(node);
			decompPricingFlag = GetDecompFlag(node);
			GetFormulaInMemorySO(node,ProductName);

			Idx1CcyTmp = (ARM_Currency*) discountCcy->Clone();
			Idx2CcyTmp = (ARM_Currency*) discountCcy->Clone();
			Idx3CcyTmp = (ARM_Currency*) discountCcy->Clone();

			// Pas de Corridor
			liborType1 = GetSpreadIndexType1(node,Idx1CcyTmp);
			liborType2 = GetSpreadIndexType2(node,Idx2CcyTmp);

			Idx1Ccy = new ARM_Currency(Idx1CcyTmp->GetCcyName());
			Idx2Ccy = new ARM_Currency(Idx2CcyTmp->GetCcyName());

			weight1 = GetWeight1(node);
			if (weight1==0.0)
				weight1 = 1.0; 
			weight2 = GetWeight2(node);
			if (weight2==0.0)
				weight2 = 1.0;

			resetTiming = GetResetTiming(node);

			resetFreq = GetResetFreq(node);
			payTiming = GetPayTiming(node);
			resetgap = GetResetGap(node);
			paygap = GetPayGap(node);
			intRule = GetIntRule(node);
			stubRule = GetStubRule(node);

			GetResetCalendar(node,resetCal);
			GetPayCalendar(node,payCal);
			// reference Date
			GetStubDate1(node,myRefDate);
			GetStubDate2(node,myRefDate2);
			char resetRolldate[4];
			GetPayRollDate(node,resetRolldate);

			ARM_Date startStubDate;
			ARM_Date endStubDate;

			if (strcmp(myRefDate, "NULL") != 0)
				startStubDate = (ARM_Date)myRefDate;
			if (strcmp(myRefDate2, "NULL") != 0)
				endStubDate = (ARM_Date)myRefDate2;

			int tmpFreq = resetFreq;
			DeduceRefDateAndStub(start,(ARM_Date)endDate,startStubDate,endStubDate,tmpFreq,atoi(resetRolldate),discountCcy->GetCcyName(),(char*)myRefDate,stubRule);

			spreadFirstFixing = GetRefValSpreadFirstFixing(node,date);
			spreadSecondFixing = GetRefValSpreadSecondFixing(node,date);

			strikesGet = GetSpreadVariable(node);// strikes correspond with the spreads in the Summit for the Spreadoption

			GetCustom(node, cust);

			int isQuantoOrNot = 0;

			if ((Idx1Ccy->GetCcyName() != discountCcy->GetCcyName()) || (Idx2Ccy->GetCcyName() != discountCcy->GetCcyName()))
				isQuantoOrNot = 1;

			if (strcmp(ProductName,"SPREADOPTIONLOGCORRIDOR")==0)
			{
				if (Idx3Ccy->GetCcyName() != discountCcy->GetCcyName())
					isQuantoOrNot = 1;
			}

			if (strcmp(cust,"CUST") != 0)
			{
				ARM_ReferenceValue vWeight1(weight1);
				ARM_ReferenceValue vWeight2(weight2);

				if (strcmp(ProductName,"SPREADOPTIONLOG")==0)
				{
					if (!isQuantoOrNot)
					{
						newSpreadOpt= new ARM_SpreadOption(start, (ARM_Date)endDate,
															 capFloor, (ARM_ReferenceValue*)NULL,
															 (ARM_INDEX_TYPE) liborType1,  
															 (ARM_INDEX_TYPE) liborType2,  
															 &vWeight1, &vWeight2,
															 basis, 
															 resetFreq, 
															 K_ZEROCOUPON, 
															 resetTiming, 
															 payTiming, 
															 discountCcy,
															 resetgap,
															 spreadFirstFixing,
															 spreadSecondFixing,
															 intRule,
															 stubRule,
															 cptStrikeMethod,
															 resetCal,
															 payCal,
															 K_SPREADOPTION_FORMULA,
															 fwdRule,
															 myRefDate,
															 paygap);
					}
					else
					{
						Idx1 = new ARM_IRIndex((ARM_INDEX_TYPE) liborType1, resetFreq, K_ZEROCOUPON, Idx1Ccy, basis);
						Idx2 = new ARM_IRIndex((ARM_INDEX_TYPE) liborType2, resetFreq, K_ZEROCOUPON, Idx2Ccy, basis);

						if ( !(( liborType1 >= K_CMS1) && ( liborType1 <= K_CMS30 ))
							&& Idx1->IsLiborIndex() )
						{
							Idx1->SetResetTiming(resetTiming);
							Idx1->SetPayTiming(payTiming);

							if (resetgap != 10000)
								Idx1->SetResetGap(resetgap);

							if (paygap != 10000)
								Idx1->SetPayGap(paygap);

							Idx1->SetIndexStyle(IN_ARREARS);

							if ( (Idx1->GetResetTiming() == K_ADVANCE)
								&& (Idx1->GetPayTiming() == K_ARREARS)
								&& (Idx1->GetResetFrequency() == Idx1->GetPayFrequency())
							   )
							{
								Idx1->SetIndexStyle(VANILLA);
							}

							Idx1->SetIntRule(intRule);

						}
						else
						{
							ARM_INDEX_TYPE Libor_Type1 = Idx1Ccy->GetVanillaIndexType();

							Idx1->Set((ARM_INDEX_TYPE) liborType1, Libor_Type1, resetFreq, -1, Idx1Ccy);

							Idx1->Set(basis, resetFreq, K_ZEROCOUPON,
												-1, -1, fwdRule, intRule,
												resetTiming, resetgap, payTiming, paygap,
												Idx1Ccy, (ARM_INDEX_TYPE) liborType1, -1);
						}

						if ( !(( liborType2 >= K_CMS1) && ( liborType2 <= K_CMS30 ))
							&& Idx2->IsLiborIndex() )
						{
							Idx2->SetResetTiming(resetTiming);
							Idx2->SetPayTiming(payTiming);

							if (resetgap != 10000)
								Idx2->SetResetGap(resetgap);

							if (paygap != 10000)
								Idx2->SetPayGap(paygap);

							Idx2->SetIndexStyle(IN_ARREARS);

							if ( (Idx2->GetResetTiming() == K_ADVANCE)
								&& (Idx2->GetPayTiming() == K_ARREARS)
								&& (Idx2->GetResetFrequency() == Idx2->GetPayFrequency())
							   )
							{
								Idx2->SetIndexStyle(VANILLA);
							}

							Idx2->SetIntRule(intRule);

						}
						else
						{

							ARM_INDEX_TYPE Libor_Type2 = Idx2Ccy->GetVanillaIndexType();

							Idx2->Set((ARM_INDEX_TYPE) liborType2, Libor_Type2, resetFreq, -1, Idx2Ccy);

							Idx2->Set(basis, resetFreq, K_ZEROCOUPON,
												-1, -1, fwdRule, intRule,
												resetTiming, resetgap, payTiming, paygap,
												Idx2Ccy, (ARM_INDEX_TYPE) liborType2, -1);
						}

						newSpreadOpt = new ARM_SpreadOption(start, (ARM_Date) endDate, 
															capFloor, (ARM_ReferenceValue*)NULL,
															Idx1, Idx2,
															&vWeight1, &vWeight2,
															basis, 
															resetFreq, 
															K_ZEROCOUPON, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															cptStrikeMethod,
															resetCal,
															payCal,
															K_SPREADOPTION_FORMULA,
															fwdRule,
															myRefDate,
															paygap);

					}
				}
				else if (strcmp(ProductName,"SPREADOPTIONLOGDIGITAL")==0)
				{
					if (!isQuantoOrNot)
					{
						ARM_ReferenceValue* payoff = GetVariableSpreadPayOff(node);
						slopeFlag = GetSpreadSlopeFlag(node);
						newSpreadOpt= new ARM_SpreadOption(start, (ARM_Date)endDate,
															capFloor, (ARM_ReferenceValue*)NULL,
															payoff,
															(ARM_INDEX_TYPE) liborType1,  
															(ARM_INDEX_TYPE) liborType2,  
															&vWeight1, &vWeight2, 
															basis, 
															resetFreq, 
															K_ZEROCOUPON, 
															resetTiming, 
															payTiming, 
															discountCcy,
															resetgap,
															spread1,
															spread2,
															spreadFirstFixing,
															spreadSecondFixing,
															intRule,
															stubRule,
															resetCal,
															payCal,
															slopeFlag,
															cptStrikeMethod,
															K_SPREADOPTION_FORMULA,
															fwdRule,
															myRefDate,
															paygap);
						if (payoff)
							delete payoff;
					}
					else
					{
					}
				}
				else if (strcmp(ProductName,"SPREADOPTIONLOGFLTDIGITAL")==0)
				{
					if (!isQuantoOrNot)
					{
						int PayOffliborType = GetSpreadPayoffLibor(node,discountCcy);
						slopeFlag = GetSpreadSlopeFlag(node);
						newSpreadOpt= new ARM_SpreadOption(start, (ARM_Date)endDate,
														capFloor, (ARM_ReferenceValue*)NULL,
														(ARM_INDEX_TYPE) PayOffliborType, 
														(ARM_INDEX_TYPE) liborType1,
														(ARM_INDEX_TYPE) liborType2,
														&vWeight1, &vWeight2,
														basis, 
														resetFreq, 
														K_ZEROCOUPON, 
														resetTiming, 
														payTiming, 
														discountCcy,
														resetgap,
														spread1,
														spread2,
														spreadFirstFixing,
														spreadSecondFixing,
														intRule,
														stubRule,
														resetCal,
														payCal,
														slopeFlag,
														cptStrikeMethod,
														NULL,
														NULL,
														K_SPREADOPTION_FORMULA,
														fwdRule,
														myRefDate,
														paygap);
					}
					{
					}
				}
				else
				{
					CCString msg((CCString)"Invalid Product Formula for getting SpreadOption \n" + ProductName);
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
				}

				if (Idx1)
					delete Idx1;
				Idx1 = NULL;

				if (Idx2)
					delete Idx2;
				Idx2 = NULL;

				if (Idx3)
					delete Idx3;
				Idx3 = NULL;

				if (Idx1Ccy)
					delete Idx1Ccy;
				Idx1Ccy = NULL;

				if (Idx2Ccy)
					delete Idx2Ccy;
				Idx2Ccy = NULL;

				if (Idx3Ccy)
					delete Idx3Ccy;
				Idx3Ccy = NULL;

			}
			else
			{
				int payFreq = GetPayFreq(node);

				ARM_SwapLeg* firstLeg = NULL;
				ARM_SwapLeg* secondLeg = NULL;
				ARM_Vector* flowStartDates = NULL;
				ARM_Vector* flowEndDates = NULL;
				ARM_Vector* resetDates = NULL;
				ARM_Vector* paymentDates = NULL;

				GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

				int compteur;

				// On coupe l'échéancier à la date de paiement qui nous intéresse
				for (int i = flowEndDates->GetSize() - 1; i>= 0; i--)
				{
					if ( fabs(flowEndDates->Elt(i)-endDate.GetJulian()) < 10. )
					{
						compteur = i;
						i = 0;
					}
				}

				ARM_Vector newflowStartDates (flowStartDates,0,compteur);
				ARM_Vector newflowEndDates (flowEndDates,0,compteur);
				ARM_Vector newresetDates (resetDates,0,compteur);
				ARM_Vector newpaymentDates (compteur+1,paymentDates->Elt(compteur));
				ARM_Vector interestDays(newflowStartDates.GetSize());

				for (i = 0; i < newflowStartDates.GetSize(); i++)
					interestDays.Elt(i) = DaysBetweenDates(basis, newflowStartDates.Elt(i), newflowEndDates.Elt(i));

				if (IsCMSIndex((ARM_INDEX_TYPE)liborType1))
				{
					firstLeg = (ARM_SwapLeg *) new ARM_CMSLeg(&newflowStartDates,
															  &newflowEndDates,
															  &newresetDates,
															  &newpaymentDates,
															  &interestDays,
															  (ARM_INDEX_TYPE)liborType1,
															  cstNotional,
															  K_RCV,
															  0.0,
															  K_COMP_PROP,
															  basis,
															  intRule,
															  resetFreq,
															  K_ZEROCOUPON,
															  10000,
															  discountCcy,
															  resetTiming,
															  stubRule,
															  resetCal,
															  payCal);
					firstLeg->SetIndexType((ARM_INDEX_TYPE)liborType1);
				}
				else
				{
					ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)liborType1,resetFreq,K_ZEROCOUPON,discountCcy);
					
					firstLeg = new ARM_SwapLeg(&newflowStartDates,
											   &newflowEndDates,
											   &newpaymentDates,
											   &newresetDates,
											   &interestDays,
											   NULL,
											   cstNotional,
											   newIndex,
											   K_RCV,
											   0.0,
											   K_FLOAT_RATE,
											   discountCcy,
											   K_NX_NONE,
											   decompPricingFlag,
											   basis);

					firstLeg->SetStartDate(start);
					firstLeg->SetEndDate((ARM_Date)endDate);

//					firstLeg->CptRealCashFlowDates();

					delete newIndex;
				}    

				if (IsCMSIndex((ARM_INDEX_TYPE)liborType2))
				{
					secondLeg = (ARM_SwapLeg *) new ARM_CMSLeg(&newflowStartDates,
															   &newflowEndDates,
															   &newresetDates,
															   &newpaymentDates,
															   &interestDays,
															   (ARM_INDEX_TYPE)liborType2,
															   cstNotional,
															   K_RCV,
															   0.0,
															   K_COMP_PROP,
															   basis,
															   intRule,
															   resetFreq,
															   K_ZEROCOUPON,
															   10000,
															   discountCcy,
															   resetTiming,
															   stubRule,
															   resetCal,
															   payCal);
					secondLeg->SetIndexType((ARM_INDEX_TYPE)liborType2);
				}
				else
				{
					ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)liborType2,resetFreq,K_ZEROCOUPON,discountCcy);
					
					secondLeg = new ARM_SwapLeg(&newflowStartDates,
												&newflowEndDates,
												&newpaymentDates,
												&newresetDates,
												&interestDays,
												NULL,
												cstNotional,
												newIndex,
												K_RCV,
												0.0,
												K_FLOAT_RATE,
												discountCcy,
												K_NX_NONE,
												decompPricingFlag,
												basis);

					secondLeg->SetStartDate(start);
					secondLeg->SetEndDate((ARM_Date)endDate);

					delete newIndex;
				}    

				ARM_ReferenceValue vWeight1(weight1);
				ARM_ReferenceValue vWeight2(weight2);


				if (strcmp(ProductName,"SPREADOPTIONLOG")==0)
				{
					
					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
													   capFloor, (ARM_ReferenceValue*)NULL,
													   &vWeight1, &vWeight2,
													   spreadFirstFixing,
													   spreadSecondFixing,
													   cptStrikeMethod);
				}
				else if (strcmp(ProductName,"SPREADOPTIONLOGDIGITAL")==0)
				{
					
					ARM_ReferenceValue* payoff = GetVariableSpreadPayOff(node);
					slopeFlag = GetSpreadSlopeFlag(node);
					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
													   capFloor, (ARM_ReferenceValue*)NULL,
													   payoff,
													   &vWeight1, &vWeight2,
													   spread1,
													   spread2,
													   spreadFirstFixing,
													   spreadSecondFixing,
													   slopeFlag,
													   cptStrikeMethod);
					if (payoff)
						delete payoff;
				}
				else if (strcmp(ProductName,"SPREADOPTIONLOGFLTDIGITAL")==0)
				{
					int PayOffliborType = GetSpreadPayoffLibor(node,discountCcy);
					slopeFlag = GetSpreadSlopeFlag(node);
					newSpreadOpt= new ARM_SpreadOption(firstLeg, secondLeg,
													   capFloor, (ARM_ReferenceValue*)NULL,
													   (ARM_INDEX_TYPE) PayOffliborType,
													   &vWeight1, &vWeight2,
													   spread1,
													   spread2,
													   spreadFirstFixing,
													   spreadSecondFixing,
													   slopeFlag,
													   cptStrikeMethod);
				}
				else
				{
					CCString msg((CCString)"Invalid Product Formula for getting SpreadOption \n" + ProductName);
					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 (char*) msg);
				}

				newSpreadOpt->GetSwapLeg()->CptTheoPayDates();
				newSpreadOpt->SetPaymentDates((ARM_Vector*)newSpreadOpt->GetSwapLeg()->GetTheoPayDates()->Clone());

				if (flowStartDates)
					delete flowStartDates;
				flowStartDates = NULL;

				if (flowEndDates)
					delete flowEndDates;
				flowEndDates = NULL;

				if (resetDates)
					delete resetDates;
				resetDates = NULL;

				if (paymentDates)
					delete paymentDates;
				paymentDates = NULL;

				if (firstLeg)
					delete firstLeg;
				firstLeg = NULL;

				if (secondLeg)
					delete secondLeg;
				secondLeg = NULL;

			}
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : SpreadOptionAssetLoader::BuildAsset()");
		} 

		// Gestion du strike
		if(strikesGet == NULL)
			strikes = new ARM_ReferenceValue(GetSpread(node)*100.0);
		else
		{
			ARM_SwapLeg* refLeg;

			if (newSpreadOpt->IsCorridorSpreadOption())
				refLeg = newSpreadOpt->GetPayIndexLeg();
			else
				refLeg = newSpreadOpt->GetSpreadLeg()->GetFirstLeg();

			strikes  = ARM_FromStartRefToResetRef(strikesGet,refLeg);
		}

		strikes->SetCalcMethod(K_STEPUP_LEFT);

		newSpreadOpt->SetStrikes(strikes);

		// Recup des Primes et Fees
		if (end == endDate)
		{
			ARM_ReferenceValue* tmpPrimes = GetPrimes(node);
			if (tmpPrimes != NULL)
			{
				newSpreadOpt->SetFee(tmpPrimes);
				delete tmpPrimes;
			}
		}

		newSpreadOpt->SetAmount(cstNotional);
		newSpreadOpt->SetPorS(PorS);
		newSpreadOpt->SetAssetId(assetId);

		if (discountCcy)
			delete discountCcy;
		discountCcy = NULL;

		if (cstNotional)
			delete cstNotional;
		cstNotional = NULL;

		if (spreadFirstFixing)
			delete spreadFirstFixing;
		spreadFirstFixing = NULL;

		if (spreadSecondFixing)
			delete spreadSecondFixing;
		spreadFirstFixing = NULL;

		if (strikesGet)
			delete strikesGet;
		strikesGet = NULL;

		if (strikes)
			delete strikes;
		strikes = NULL;

		return newSpreadOpt;
	}
};


//! \brief loader for Spread Option
class CMSAssetLoader : public AbstractAssetParser
{
public:
	CMSAssetLoader()
	{
	}

	~CMSAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_CMSLeg* cmsLeg = NULL;

		ARM_Date startDate, endDate;
		ARM_Currency* discountCcy = NULL;
		ARM_ReferenceValue* notional = NULL;
		ARM_ReferenceValue* spreads = NULL;
		int rcvOrPay, cmsIndexType, yieldDecompFreq, swapLegDayCount, intRule, resetFreq, payFreq, resetGap, resetTiming;
		int stubRule, adjStartDate, fwdRule, paygap;
		double spread;
		char RollDate[4];
		char myRefDate[11];
		char myRefDate2[11];

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			startDate = GetStartDate(node);
			endDate = GetEndDate(node);
			discountCcy = GetCcy(node);
			cmsIndexType = GetSpreadIndexType1(node,discountCcy);
			rcvOrPay = GetPorS(node);
			spread = GetSpread(node)*100;
			spreads = GetSpreadVariable(node);
			swapLegDayCount = GetAssetDayCount(node);
			intRule = GetIntRule(node);
			fwdRule = GetFwdRule(node);
			resetFreq = GetResetFreq(node);
			payFreq = GetPayFreq(node);
			resetGap = GetResetGap(node);
			resetTiming = GetResetTiming(node);
			stubRule = GetStubRule(node);
			GetResetRollDate(node,RollDate);
			yieldDecompFreq = K_COMP_PROP;
			paygap = GetPayGap(node);

			GetStubDate1(node,myRefDate);
			GetStubDate2(node,myRefDate2);
			char resetRolldate[4];
			GetPayRollDate(node,resetRolldate);

			ARM_Date startStubDate;
			ARM_Date endStubDate;

			if (strcmp(myRefDate, "NULL") != 0)
				startStubDate = (ARM_Date)myRefDate;
			if (strcmp(myRefDate2, "NULL") != 0)
				endStubDate = (ARM_Date)myRefDate2;

			int tmpFreq = payFreq;
			if (payFreq == K_ZEROCOUPON)
				tmpFreq = resetFreq;
			DeduceRefDateAndStub(startDate,endDate,startStubDate,endStubDate,tmpFreq,atoi(resetRolldate),discountCcy->GetCcyName(),(char*)myRefDate,stubRule);

			char ProductName[30];
			GetFormula(node,ProductName);
			if ( (!strcmp(ProductName, "CUO_SPLK"))
				|| (!strcmp(ProductName, "CUI_SPLK"))
				|| (!strcmp(ProductName, "PDO_SPLK"))
				|| (!strcmp(ProductName, "PDI_SPLK"))
				|| (!strcmp(ProductName, "DIGITAL_SPLK")) )
			{
				char* arg = "";
				GetNthArgInFormula(node, 5, arg);
				// SPLK(index,strike,barrier,round,decomp,...)
				if (atoi(arg) == 1)
					yieldDecompFreq = payFreq;
				else if (atoi(arg) == 2)
					yieldDecompFreq = resetFreq;
			}

/*			if (resetFreq != K_DAILY)
			{
				try
				{
					endDate.ChangeDate(atoi(RollDate),endDate.GetMonth(),endDate.GetYear());
				}
				catch (Exception& ) 
				{
					endDate.ChangeDate(DaysInMonth(endDate.GetMonth(),endDate.GetYear()),endDate.GetMonth(),endDate.GetYear());
				}
			}
*/
			adjStartDate = 0;
			char resetCal[4];
			char payCal[4];

			char cust[20];
			GetCustom(node, cust);
			if (strcmp(cust,"CUST") != 0)
			{
				notional = GetNotional(node);//amortizing conv à ARM_Curve

				ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

				cmsLeg = new ARM_CMSLeg(startDate,
										endDate,
										(ARM_INDEX_TYPE)cmsIndexType,
										rcvOrPay,
										spread,
										yieldDecompFreq,
										swapLegDayCount,
										intRule,
										resetFreq,
										payFreq,
										resetGap,
										discountCcy,
										resetTiming,
										stubRule,
										NULL, //resetCal,
										NULL, //payCal,
										adjStartDate,
										fwdRule,
										myRefDate,
										paygap);

				// notional dates must match payment dates, which are always adjusted
				// if intRule = 0, these notional dates will not be adjusted in leg construction
				// so we must adjust them by ourselves
				if ( (intRule == K_UNADJUSTED)
					&&
					 (notional->GetDiscreteDates() != NULL)
					 )
				{
					int size = notional->GetDiscreteDates()->GetSize();
					for (int i = 0; i < size; i++)
					{
						ARM_Date notionalCurDate (notional->GetDiscreteDates()->Elt(i));
						for (int j = i; j < cmsLeg->GetTheoPayDates()->GetSize(); j++)
						{
							double curPayDate = cmsLeg->GetTheoPayDates()->Elt(j);
							if ( notionalCurDate.AdjustToBusDate(discountCcy->GetCcyName(), K_ADJUSTED).GetJulian() == curPayDate )
							{
								notional->GetDiscreteDates()->Elt(i) = curPayDate;
								break;
							}
						}
					}
				}
				cmsLeg->SetAmount(notional);
				cmsLeg->SetFixRates(vFixing);

				delete discountCcy;
				delete notional;
				if (vFixing)
					delete vFixing;
				vFixing = NULL;
			}
			else
			{
				// Build leg from schedule
				ARM_Vector* flowStartDates = NULL;
				ARM_Vector* flowEndDates = NULL;
				ARM_Vector* resetDates = NULL;
				ARM_Vector* paymentDates = NULL;
				ARM_Vector* interestDays = NULL;

				// Parsing method differs whether the leg is part of a SWAP or of another asset
				CCString assetType = GetCCStringFromXMLNode(node,"Type");
				if (assetType == "SWAP")
					GetSwapCustomSchedule(node,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);
				else
					GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

				interestDays = new ARM_Vector(flowStartDates->GetSize());
				for (int i = 0; i < flowStartDates->GetSize(); i++)
					interestDays->Elt(i) = DaysBetweenDates(swapLegDayCount, flowStartDates->Elt(i), flowEndDates->Elt(i));

				notional = GetNotionalCust(node);
				GetResetCalendar(node,resetCal);
				GetPayCalendar(node,payCal);

				ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)cmsIndexType,resetFreq,payFreq,discountCcy);

				cmsLeg = new ARM_CMSLeg (	flowStartDates, flowEndDates, resetDates, paymentDates, interestDays,
											(ARM_INDEX_TYPE)cmsIndexType, notional, rcvOrPay, spread, yieldDecompFreq,
											swapLegDayCount, intRule, resetFreq, payFreq, resetGap,
											discountCcy, resetTiming, stubRule, resetCal, payCal );

				cmsLeg->SetStartDate(startDate);
				cmsLeg->SetEndDate(endDate);

				delete newIndex;
			}

			// Gestion du spread
			if(spreads != NULL)
			{
				ARM_ReferenceValue* newspreads  = ARM_FromStartRefToResetRef(spreads,cmsLeg);
				newspreads->SetCalcMethod(K_STEPUP_LEFT);

				cmsLeg->SetVariableSpread(newspreads);

				if (newspreads)
					delete newspreads;
				newspreads = NULL;

				if (spreads)
					delete spreads;
				spreads = NULL;
			}
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "CMSAssetLoader::BuildAsset : Failed in constructing CMSLeg");
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : CMSAssetLoader::BuildAsset()");
		} 

		cmsLeg->SetAssetId(assetId);

		return cmsLeg;
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate,double notional,int isEtkOrNot = 1)
	{
		ARM_CMSLeg* cmsLeg = NULL;

		ARM_Date startDate;
		ARM_Currency* discountCcy = NULL;
		ARM_ReferenceValue* spreads = NULL;
		int rcvOrPay, cmsIndexType, yieldDecompFreq, swapLegDayCount, intRule, resetFreq, payFreq, resetGap, resetTiming;
		int stubRule, adjStartDate, fwdRule;
		double spread;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			startDate = GetStartDate(node);
			discountCcy = GetCcy(node);
			cmsIndexType = GetSpreadIndexType1(node,discountCcy);
			rcvOrPay = GetPorS(node);
			spread = GetSpread(node)*100;
			spreads = GetSpreadVariable(node);
			yieldDecompFreq = K_COMP_PROP;
			swapLegDayCount = GetAssetDayCount(node);
			intRule = GetIntRule(node);
			fwdRule = GetFwdRule(node);
			resetFreq = GetResetFreq(node);
			resetGap = GetResetGap(node);
			resetTiming = GetResetTiming(node);
			stubRule = GetStubRule(node);
			adjStartDate = 0;
			char resetCal[4];
			char payCal[4];

			char cust[20];
			GetCustom(node, cust);
			if (strcmp(cust,"CUST") != 0)
			{
				ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

				cmsLeg = new ARM_CMSLeg(startDate,
										(ARM_Date)endDate,
										(ARM_INDEX_TYPE)cmsIndexType,
										rcvOrPay,
										spread,
										yieldDecompFreq,
										swapLegDayCount,
										intRule,
										resetFreq,
										K_ZEROCOUPON,
										resetGap,
										discountCcy,
										resetTiming,
										stubRule,
										NULL, //resetCal,
										NULL, //payCal,
										adjStartDate,
										fwdRule);

				ARM_ReferenceValue cstNotional(notional);
				cmsLeg->SetAmount(&cstNotional);
				cmsLeg->SetFixRates(vFixing);

				delete discountCcy;

				if (vFixing)
					delete vFixing;
				vFixing = NULL;
			}
			else
			{
				payFreq = GetPayFreq(node);

				// Build leg from schedule
				ARM_Vector* flowStartDates = NULL;
				ARM_Vector* flowEndDates = NULL;
				ARM_Vector* resetDates = NULL;
				ARM_Vector* paymentDates = NULL;

				// Parsing method differs whether the leg is part of a SWAP or of another asset
				CCString assetType = GetCCStringFromXMLNode(node,"Type");
				if (assetType == "SWAP")
					GetSwapCustomSchedule(node,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);
				else
					GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

				int compteur;

				// On coupe l'échéancier à la date de paiement qui nous intéresse
				for (int i = flowEndDates->GetSize() - 1; i>= 0; i--)
				{
					if ( fabs(flowEndDates->Elt(i)-endDate.GetJulian()) < 10. )
					{
						compteur = i;
						i = 0;
					}
				}

				ARM_Vector newflowStartDates (flowStartDates,0,compteur);
				ARM_Vector newflowEndDates (flowEndDates,0,compteur);
				ARM_Vector newresetDates (resetDates,0,compteur);
				ARM_Vector newpaymentDates (compteur+1,paymentDates->Elt(compteur));
				ARM_Vector interestDays(newflowStartDates.GetSize());

				for (i = 0; i < newflowStartDates.GetSize(); i++)
					interestDays.Elt(i) = DaysBetweenDates(swapLegDayCount, newflowStartDates.Elt(i), newflowEndDates.Elt(i));

				GetResetCalendar(node,resetCal);
				GetPayCalendar(node,payCal);

				ARM_ReferenceValue cstNotional(notional);
				ARM_IRIndex* newIndex = new ARM_IRIndex((ARM_INDEX_TYPE)cmsIndexType,resetFreq,K_ZEROCOUPON,discountCcy);

				cmsLeg = new ARM_CMSLeg (	&newflowStartDates, &newflowEndDates, &newresetDates, &newpaymentDates, &interestDays,
											(ARM_INDEX_TYPE)cmsIndexType, &cstNotional, rcvOrPay, spread, yieldDecompFreq,
											swapLegDayCount, intRule, resetFreq, payFreq, resetGap,
											discountCcy, resetTiming, stubRule, resetCal, payCal );

				cmsLeg->SetStartDate(startDate);
				cmsLeg->SetEndDate((ARM_Date)endDate);

				delete newIndex;

				// Gestion du spread
				if(spreads != NULL)
				{
					ARM_ReferenceValue* newspreads  = ARM_FromStartRefToResetRef(spreads,cmsLeg);
					newspreads->SetCalcMethod(K_STEPUP_LEFT);

					cmsLeg->SetVariableSpread(newspreads);

					if (newspreads)
						delete newspreads;
					newspreads = NULL;

					if (spreads)
						delete spreads;
					spreads = NULL;

					if (flowStartDates)
						delete flowStartDates;
					flowStartDates = NULL;

					if (flowEndDates)
						delete flowEndDates;
					flowEndDates = NULL;

					if (resetDates)
						delete resetDates;
					resetDates = NULL;

					if (paymentDates)
						delete paymentDates;
					paymentDates = NULL;
				}
			}
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : CMSAssetLoader::BuildAssetForMemory()");
		} 

		cmsLeg->SetAssetId(assetId);

		return cmsLeg;
	}
};


//! \brief loader for Spread Option
class CMTAssetLoader : public AbstractAssetParser
{
public:
	CMTAssetLoader()
	{
	}

	~CMTAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_CMTLeg* cmtLeg = NULL;

		ARM_Date startDate, endDate;
		ARM_Currency* discountCcy = NULL;
		ARM_ReferenceValue* notional = NULL;
		int rcvOrPay, cmtIndexType, yieldDecompFreq, swapLegDayCount, intRule, resetFreq, payFreq,
			resetGap, resetTiming, bondCouponFreq, bondDayCount, stubRule;
		double spread;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			discountCcy = GetCcy(node);
			cmtIndexType = GetSpreadIndexType1(node, discountCcy, 1);
			swapLegDayCount = GetAssetDayCount(node);
			resetFreq = GetResetFreq(node);
			payFreq = GetPayFreq(node);

			startDate = GetStartDate(node);
			endDate = GetEndDate(node);
			bondCouponFreq = GetYieldFreq(node); //K_ANNUAL; // ????
			bondDayCount = swapLegDayCount; //K30_360; // ????
			rcvOrPay = GetPorS(node);
			spread = GetSpread(node)*100;
			yieldDecompFreq = K_COMP_QUARTERLY; //K_COMP_PROP;
			intRule = GetIntRule(node);
			resetGap = GetResetGap(node);
			char cust[20];
			GetCustom(node, cust);
			notional = GetNotional(node);//amortizing conv à ARM_Curve
			resetTiming = GetResetTiming(node);
			stubRule = GetStubRule(node);
			
			ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

			cmtLeg = new ARM_CMTLeg(startDate, endDate, (ARM_INDEX_TYPE)cmtIndexType, 
										bondCouponFreq, bondDayCount, rcvOrPay, spread, 
										yieldDecompFreq, swapLegDayCount, intRule, resetGap,
										resetFreq, 0,
										discountCcy, resetTiming, stubRule);
			
			// notional dates must match payment dates, which are always adjusted
			// if intRule = 0, these notional dates will not be adjusted in leg construction
			// so we must adjust them by ourselves
			if ( (intRule == K_UNADJUSTED)
				&&
				 (notional->GetDiscreteDates() != NULL)
				 )
			{
				int size = notional->GetDiscreteDates()->GetSize();
				for (int i = 0; i < size; i++)
				{
					ARM_Date notionalCurDate (notional->GetDiscreteDates()->Elt(i));
					for (int j = i; j < cmtLeg->GetTheoPayDates()->GetSize(); j++)
					{
						double curPayDate = cmtLeg->GetTheoPayDates()->Elt(j);
						if ( notionalCurDate.AdjustToBusDate(discountCcy->GetCcyName(), K_ADJUSTED).GetJulian() == curPayDate )
						{
							notional->GetDiscreteDates()->Elt(i) = curPayDate;
							break;
						}
					}
				}
			}
			cmtLeg->SetAmount(notional);
			cmtLeg->SetFixRates(vFixing);

			delete discountCcy;
			delete notional;
			if (vFixing)
				delete vFixing;
			vFixing = NULL;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "CMSAssetLoader::BuildAsset : Failed in constructing CMSLeg");
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : CMSAssetLoader::BuildAsset()");
		} 

		cmtLeg->SetAssetId(assetId);

		return cmtLeg;
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate,double notional,int isEtkOrNot = 1)
	{
		ARM_CMTLeg* cmtLeg = NULL;

		ARM_Date startDate;
		ARM_Currency* discountCcy = NULL;
		int rcvOrPay, cmtIndexType, yieldDecompFreq, swapLegDayCount, intRule, resetFreq,
			resetGap, resetTiming, bondDayCount, stubRule;
		double spread;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			discountCcy = GetCcy(node);
			cmtIndexType = GetSpreadIndexType1(node, discountCcy, 1);
			swapLegDayCount = GetAssetDayCount(node);
			resetFreq = GetResetFreq(node);

			startDate = GetStartDate(node);
			bondDayCount = swapLegDayCount; //K30_360; // ????
			rcvOrPay = GetPorS(node);
			spread = GetSpread(node)*100;
			yieldDecompFreq = K_COMP_QUARTERLY; //K_COMP_PROP;
			intRule = GetIntRule(node);
			resetGap = GetResetGap(node);
			char cust[20];
			GetCustom(node, cust);
			resetTiming = GetResetTiming(node);
			stubRule = GetStubRule(node);
			
			ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

			cmtLeg = new ARM_CMTLeg(startDate, (ARM_Date)endDate, (ARM_INDEX_TYPE)cmtIndexType, 
										K_ZEROCOUPON, bondDayCount, rcvOrPay, spread, 
										yieldDecompFreq, swapLegDayCount, intRule, resetGap,
										resetFreq, 0,
										discountCcy, resetTiming, stubRule);

			ARM_ReferenceValue cstNotional(notional);
			cmtLeg->SetAmount(&cstNotional);
			cmtLeg->SetFixRates(vFixing);

			delete discountCcy;
			if (vFixing)
				delete vFixing;
			vFixing = NULL;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				 "Unrecognized Failure in : CMSAssetLoader::BuildAssetForMemory()");
		} 

		cmtLeg->SetAssetId(assetId);

		return cmtLeg;
	}
};

//! \brief loader for Fix Leg
class FixAssetLoader : public AbstractAssetParser
{
public:
	FixAssetLoader()
	{
	}

	~FixAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_FixLeg* fixLeg = NULL;
		ARM_Date startDate,fstCpnEffDate,endDate;
		ARM_ReferenceValue* fixCoupon = NULL;
		int fixFrequency;
		int daycountBasis;
		int decompFreq;
		int payTiming;
		ARM_Currency* discountCcy = NULL;
		ARM_ReferenceValue* notional = NULL;
		char payCal[7];
		char myRefDate[11];
		long NxId;
		long intRule;
		int rcvOrPay;
		char payRollDate[4];
		int rollDay;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			startDate = GetStartDate(node);
			endDate = GetEndDate(node);
			rcvOrPay = GetPorS(node);

			// fix coupon of initPeriodLeg
			fixCoupon = GetStepUpFixCoupon(node);

			// fix frequency of initPeriodLeg
			fixFrequency = GetPayFreq(node);

			// daycount Basis of initPeriodLeg
			daycountBasis = GetInterestDayCount(node);

			decompFreq = GetDecompFreq(node);
			payTiming = GetPayTiming(node);

			// int rule of initPeriodLeg
			intRule = GetPayIntRule(node);

			// currency of initPeriodLeg
			discountCcy = GetCcy(node);

			// notional of initPeriodLeg
			notional = GetNotional(node);

			// paycal of initPeriodLeg
			GetPayCalendar(node,payCal);

			// reference Date
			GetStubDate1(node,myRefDate);

			// Echange de Notional
			NxId = GetNotionalExchangeFlag(node);

			// First Coupon Effective Date
			fstCpnEffDate = GetFstEffDate(node,startDate);

			// Roll Date
			GetPayRollDate(node, payRollDate);
			rollDay	=	atoi(payRollDate);

			long stubRuleId = K_SHORTSTART;

			fixLeg = new ARM_FixLeg(startDate,endDate,0.0,rcvOrPay,
									fixFrequency,daycountBasis,decompFreq,payTiming,intRule,
									stubRuleId,discountCcy,payCal,NxId,myRefDate,1,rollDay);

			fixLeg->SetVarCoupons(fixCoupon);

			if (fixCoupon)
				delete fixCoupon;

			// notional dates must match payment dates, which are always adjusted
			// if intRule = 0, these notional dates will not be adjusted in leg construction
			// so we must adjust them by ourselves
			if ( (intRule == K_UNADJUSTED)
				&&
				 (notional->GetDiscreteDates() != NULL)
				 )
			{
				int size = notional->GetDiscreteDates()->GetSize();
				for (int i = 0; i < size; i++)
				{
					ARM_Date notionalCurDate (notional->GetDiscreteDates()->Elt(i));
					for (int j = i; j < fixLeg->GetTheoPayDates()->GetSize(); j++)
					{
						double curPayDate = fixLeg->GetTheoPayDates()->Elt(j);
						if ( notionalCurDate.AdjustToBusDate(payCal, K_ADJUSTED).GetJulian() == curPayDate )
						{
							notional->GetDiscreteDates()->Elt(i) = curPayDate;
							break;
						}
					}
				}
			}
			fixLeg->SetAmount(notional);

			if (fstCpnEffDate != startDate)
				fixLeg->CustomizeFirstPeriod(fstCpnEffDate);

			delete notional;
			delete discountCcy;

		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "FixAssetLoader::BuildAsset : Failed in constructing FixLeg");
		} 

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : FixAssetLoader::BuildAsset()");
		} 

		fixLeg->SetAssetId(assetId);

		return fixLeg;
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate, double notional,int isEtkOrNot = 1)
	{
		ARM_SwapLeg* fixLeg = NULL;
		ARM_Date startDate,fstCpnEffDate;
		double fixCoupon;
		int fixFrequency;
		int daycountBasis;
		int decompFreq;
		int payTiming;
		ARM_Currency* discountCcy = NULL;
		char payCal[7];
		char myRefDate[11];
		long NxId;
		long intRule;
		int rcvOrPay;

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			startDate = GetStartDate(node);
			rcvOrPay = GetPorS(node);

			// fix coupon of initPeriodLeg
			fixCoupon = GetDoubleFromXMLNode(node,"INTEREST_Rate") * 100.;

			// fix frequency of initPeriodLeg
			fixFrequency = GetPayFreq(node);

			// daycount Basis of initPeriodLeg
			daycountBasis = GetInterestDayCount(node);

			decompFreq = GetDecompFreq(node);
			payTiming = GetPayTiming(node);

			// int rule of initPeriodLeg
			intRule = GetPayIntRule(node);

			// currency of initPeriodLeg
			discountCcy = GetCcy(node);

			// paycal of initPeriodLeg
			GetPayCalendar(node,payCal);

			// reference Date
			GetStubDate1(node,myRefDate);

			// Echange de Notional
			NxId = GetNotionalExchangeFlag(node);

			// First Coupon Effective Date
			fstCpnEffDate = GetFstEffDate(node,startDate);

			long stubRuleId = K_SHORTSTART;

			fixLeg = new ARM_SwapLeg(startDate,(ARM_Date)endDate,fixCoupon,rcvOrPay,
									K_ZEROCOUPON,daycountBasis,decompFreq,payTiming,intRule,
									stubRuleId,discountCcy,payCal,NxId,myRefDate);

			ARM_ReferenceValue cstNotional(notional);
			fixLeg->SetAmount(&cstNotional);

			if (fstCpnEffDate != startDate)
				fixLeg->CustomizeFirstPeriod(fstCpnEffDate);

			delete discountCcy;

		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : FixAssetLoader::BuildAssetForMemory()");
		} 

		fixLeg->SetAssetId(assetId);

		return fixLeg;
	}
};


//! \brief loader for Floating Leg
class FloatingAssetLoader : public AbstractAssetParser
{
public:
	FloatingAssetLoader()
	{
	}

	~FloatingAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_SwapLeg* floatingLeg = NULL;

		ARM_Date startDate, endDate;
		double spread = 0.0;
		ARM_ReferenceValue* stepUpSpread = NULL;
		int rcvOrPay, dayCount, liborType, resetFreq, resetTiming, payFreq, 
			payTiming, intRule, stubRule, nxChange, resetGap, decompPricingFlag, adjStartDate, fwdRule,
			rollDay, compMode;
		ARM_Currency* discountCcy;
		ARM_Currency* resetCcy;
		ARM_ReferenceValue* notional = NULL;
		char payCalName[4];
		char resetCalName[4];
		char refDate[11];
		char RollDate[4];
		strcpy(refDate,"NULL");

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			discountCcy = GetCcy(node);
			// WARNING : we must clone CCY !
			// in case of quanto asset, GetFormulaIndexType() will set a new value for CCY
			resetCcy = (ARM_Currency*)discountCcy->Clone();

			bool isQuanto = false;
			liborType = GetFormulaIndexType(node, resetCcy);
			
			if (liborType != K_LIVRET_A)
			{
				if (liborType != -1)
				{
					isQuanto = true;
				}
				else
				{
					liborType = GetIndexType(node);
				}
			}

			char ProductName[30];
			GetFormula(node,ProductName);
			// Careful : formula can specify a FIX RATE !!
			// so we must create a fixed leg

			if ( (string(ProductName).find("FIX") != -1) && (GetDoubleFromXMLNode(node, "INTEREST_Rate") != 0.0) )
			{
				FixAssetLoader loader;
				floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId,isEtkOrNot));
			}

			if ( string(ProductName).find("CMT") != -1 )
			{
				CMTAssetLoader loader;
				floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date, assetId, isEtkOrNot));
			}
			else if ( (string(ProductName).find("CMS") != -1) || (string(ProductName).find("TEC") != -1) 
					|| (IsCMSIndex((ARM_INDEX_TYPE)liborType)) )
			{
				CMSAssetLoader loader;
				floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date, assetId, isEtkOrNot));
			}
			else 
			{
				startDate		= GetStartDate(node);
				endDate			= GetEndDate(node);
				rcvOrPay		= GetPorS(node);
				dayCount		= GetAssetDayCount(node);
				intRule			= GetPayIntRule(node); //same value as ResetIntRule returned by GetIntRule()
				resetGap		= GetResetGap(node);
				resetFreq		= GetResetFreq(node);
				payFreq			= GetPayFreq(node);
				resetTiming		= GetResetTiming(node);
				stubRule		= GetStubRule(node);
				fwdRule			= GetFwdRule(node);
				compMode		= GetCompMode(node);
				
				if(strcmp(ProductName, "GLOBALCAP") == 0)
					liborType = GetGlobalCapIndexType(node, discountCcy);

				GetResetRollDate(node,RollDate);
				try
				{
					rollDay	= atoi(RollDate);
					endDate.ChangeDate(atoi(RollDate),endDate.GetMonth(),endDate.GetYear());
				}
				catch (Exception& ) 
				{
					endDate.ChangeDate(DaysInMonth(endDate.GetMonth(),endDate.GetYear()),endDate.GetMonth(),endDate.GetYear());
				}
				
				//spread considered as an event.
				stepUpSpread	= GetSwaplegStepUpSpread(node, false);
				//spread considered as a static data.
				if (!stepUpSpread)
				{
					spread		= GetSpread(node) * 100;
				}

				decompPricingFlag		= GetDecompFlag(node);
				nxChange				= GetNotionalExchangeFlag(node);
				ARM_IRIndex* newIndex	= NULL;

				char cust[20];
				GetCustom(node, cust);
				if (strcmp(cust,"CUST") != 0)
				{
					payTiming		= GetPayTiming(node);
					GetResetCalendar(node, resetCalName);
					GetPayCalendar(node, payCalName);
					adjStartDate	= 0;

					ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

					newIndex = new ARM_IRIndex(resetCcy->GetLiborIndexDayCount(),//dayCount,
											   resetFreq,
											   payFreq,
											   1./(double)FromLiborTypeToFrequency(liborType),
											   K_COMP_PROP,
											   fwdRule,
											   intRule,
											   resetTiming,
											   resetGap,
											   payTiming,
											   0,
											   resetCcy,
											   (ARM_INDEX_TYPE)liborType);
					floatingLeg = new ARM_SwapLeg(startDate,
													  endDate,
													  newIndex,
													  rcvOrPay,
													  spread,
													  stubRule,
													  K_COMP_PROP,
													  discountCcy,
													  dayCount,
													  resetGap,
													  resetCalName,
													  payCalName,
													  decompPricingFlag,
													  nxChange,
													  refDate,
													  adjStartDate);

					if (stepUpSpread)
					{
						ARM_ReferenceValue* stepUpSpreadOnReset = ARM_FromStartRefToResetRef(stepUpSpread, floatingLeg);
						delete stepUpSpread;

						floatingLeg->SetVariableSpread(stepUpSpreadOnReset);
						delete stepUpSpreadOnReset;
					}

					notional = GetNotional(node);//amortizing conv à ARM_Curve
					// notional dates must match payment dates, which are always adjusted
					// if intRule = 0, these notional dates will not be adjusted in leg construction
					// so we must adjust them by ourselves
					if ( (intRule == K_UNADJUSTED) && (notional->GetDiscreteDates() != NULL) )
					{
						int size = notional->GetDiscreteDates()->GetSize();
						for (int i = 0; i < size; i++)
						{
							ARM_Date notionalCurDate (notional->GetDiscreteDates()->Elt(i));
							for (int j = i; j < floatingLeg->GetTheoPayDates()->GetSize(); j++)
							{
								double curPayDate = floatingLeg->GetTheoPayDates()->Elt(j);
								if ( notionalCurDate.AdjustToBusDate(payCalName, K_ADJUSTED).GetJulian() == curPayDate )
								{
									notional->GetDiscreteDates()->Elt(i) = curPayDate;
									break;
								}
							}
						}
					}
					floatingLeg->SetAmount(notional);
					floatingLeg->SetFixRates(vFixing);

					if (vFixing)
						delete vFixing;
					vFixing = NULL;
				}
				else
				{
					// Build leg from schedule
					ARM_Vector* flowStartDates	= NULL;
					ARM_Vector* flowEndDates	= NULL;
					ARM_Vector* resetDates		= NULL;
					ARM_Vector* paymentDates	= NULL;
					ARM_Vector* interestDays	= NULL;

					// Parsing method differs whether the leg is part of a SWAP or of another asset
					CCString assetType = GetCCStringFromXMLNode(node,"Type");
					if (assetType == "SWAP")
						GetSwapCustomSchedule(node,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);
					else
						GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

					interestDays = new ARM_Vector(flowStartDates->GetSize());
					for (int i = 0; i < flowStartDates->GetSize(); i++)
						interestDays->Elt(i) = DaysBetweenDates(dayCount, flowStartDates->Elt(i), flowEndDates->Elt(i));

					notional = GetNotionalCust(node);

					newIndex	= new ARM_IRIndex( (ARM_INDEX_TYPE)liborType,
													resetFreq,
													payFreq,
													discountCcy);

					floatingLeg = new ARM_SwapLeg (	flowStartDates, 
													flowEndDates, 
													paymentDates, 
													resetDates,
													interestDays, 
													NULL, 
													notional, 
													newIndex, 
													rcvOrPay, 
													spread,
													K_FLOAT_RATE, 
													discountCcy, 
													nxChange, 
													decompPricingFlag );

					if (stepUpSpread)
					{
						ARM_ReferenceValue* stepUpSpreadOnReset = ARM_FromStartRefToResetRef(stepUpSpread, floatingLeg);
						delete stepUpSpread;

						floatingLeg->SetVariableSpread(stepUpSpreadOnReset);
						delete stepUpSpreadOnReset;
					}

					ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

					floatingLeg->SetFixRates(vFixing);

					if (vFixing)
						delete vFixing;
					vFixing = NULL;

					floatingLeg->SetStartDate(startDate);
					floatingLeg->SetEndDate(endDate);

					if (newIndex)
					{
						delete newIndex;
						newIndex = NULL;
					}
					if (flowStartDates)
					{
						delete flowStartDates;
						flowStartDates = NULL;
					}
					if (flowEndDates)
					{
						delete flowEndDates;
						flowEndDates = NULL;
					}
					if (paymentDates)
					{
						delete paymentDates;
						paymentDates = NULL;
					}
					if (resetDates)
					{
						delete resetDates;
						resetDates = NULL;
					}
					if (interestDays)
					{
						delete interestDays;
						interestDays = NULL;
					}
				}
				delete newIndex;
				newIndex = NULL;

				delete notional;
				notional = NULL;
			}

			delete discountCcy;
			discountCcy = NULL;
			
			if (resetCcy)
			{
				delete resetCcy;
				resetCcy = NULL;
			}
			
			floatingLeg->SetCompoundingType(compMode);

			floatingLeg->SetAssetId(assetId);
			if (isQuanto)
			{
				floatingLeg->SetQuantoAsset(1);
			}

			return floatingLeg;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : FloatingAssetLoader::BuildAsset()");
		} 
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId, const ARM_Date& endDate, double notional, int isEtkOrNot = 1)
	{
		ARM_SwapLeg* floatingLeg = NULL;

		ARM_Date startDate;
		double spread = 0.0;
		ARM_ReferenceValue* stepUpSpread = NULL;
		int rcvOrPay, dayCount, liborType, resetFreq, resetTiming, payFreq, 
			payTiming, intRule, stubRule, nxChange, resetGap, decompPricingFlag, adjStartDate;
		ARM_Currency* discountCcy;
		ARM_Currency* resetCcy;
		char payCalName[4];
		char resetCalName[4];
		char refDate[11];
		strcpy(refDate,"NULL");

		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			discountCcy = GetCcy(node);
			// WARNING : we must clone CCY !
			// in case of quanto asset, GetQuantoIndexType() will set a new value for CCY
			resetCcy = (ARM_Currency*)discountCcy->Clone();

			bool isQuanto = false;
			liborType = GetFormulaIndexType(node, resetCcy);
			if (liborType != K_LIVRET_A)
			{
				if (liborType != -1)
				{
					isQuanto = true;
				}
				else
				{
					liborType = GetIndexType(node);
				}
			}
			
			char ProductName[30];
			GetFormula(node,ProductName);
			// Careful : formula can specify a FIX RATE !!
			// so we must create a fixed leg
			if ( (string(ProductName).find("FIX") != -1) && (GetDoubleFromXMLNode(node, "INTEREST_Rate") != 0.0) )
			{
				FixAssetLoader loader;
				floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));
			}

			if ( (string(ProductName).find("CMS") != -1) || (string(ProductName).find("TEC") != -1) || (IsCMSIndex((ARM_INDEX_TYPE)liborType)) )
			{
				CMSAssetLoader loader;
				floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));
			}
			else 
			{
				startDate		= GetStartDate(node);
				rcvOrPay		= GetPorS(node);
				dayCount		= GetAssetDayCount(node);
				intRule			= GetPayIntRule(node); //same value as ResetIntRule returned by GetIntRule()
				resetGap		= GetResetGap(node);
				resetFreq		= GetResetFreq(node);
				resetTiming		= GetResetTiming(node);
				stubRule		= GetStubRule(node);

				//spread considered as an event.
				stepUpSpread	= GetSwaplegStepUpSpread(node,false);
				//spread considered as a static data.
				if (!stepUpSpread)
				{
					spread		= GetSpread(node) * 100;
				}

				if (string(ProductName).find("CMT") != -1)
				{
					CMTAssetLoader loader;
					floatingLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));
				}
				else 
				{
					decompPricingFlag		= GetDecompFlag(node);
					nxChange				= GetNotionalExchangeFlag(node);
					ARM_IRIndex* newIndex	= NULL;

					char cust[20];
					GetCustom(node, cust);
					if (strcmp(cust,"CUST") != 0)
					{
						payTiming		= GetPayTiming(node);
						GetResetCalendar(node, resetCalName);
						GetPayCalendar(node, payCalName);
						adjStartDate	= 0;

						ARM_ReferenceValue* vFixing = GetRefValSpreadFirstFixing(node,date);

						newIndex = new ARM_IRIndex( dayCount, 
												    resetFreq, 
													K_ZEROCOUPON, 
													1./(double)FromLiborTypeToFrequency(liborType), 
													K_COMP_PROP, 
													K_MOD_FOLLOWING, 
													intRule, 
													resetTiming, 
													resetGap, 
													payTiming, 
													0, 
													resetCcy,
													(ARM_INDEX_TYPE)liborType);

						floatingLeg	= new ARM_SwapLeg(startDate, 
													(ARM_Date)endDate, 
													newIndex, 
													rcvOrPay, 
													spread, 
													stubRule, 
													K_COMP_PROP, 
													discountCcy, 
													dayCount, 
													resetGap, 
													resetCalName, 
													payCalName, 
													decompPricingFlag,
													nxChange, 
													refDate, 
													adjStartDate );

						if (stepUpSpread)
						{
							ARM_ReferenceValue* stepUpSpreadOnReset = ARM_FromStartRefToResetRef(stepUpSpread, floatingLeg);
							delete stepUpSpread;

							floatingLeg->SetVariableSpread(stepUpSpreadOnReset);
							delete stepUpSpreadOnReset;
						}

						ARM_ReferenceValue cstNotional(notional);
						floatingLeg->SetAmount(&cstNotional);
						floatingLeg->SetFixRates(vFixing);
					}
					else
					{
						payFreq = GetPayFreq(node);

						// Build leg from schedule
						ARM_Vector* flowStartDates	= NULL;
						ARM_Vector* flowEndDates	= NULL;
						ARM_Vector* resetDates		= NULL;
						ARM_Vector* paymentDates	= NULL;

						// Parsing method differs whether the leg is part of a SWAP or of another asset
						CCString assetType = GetCCStringFromXMLNode(node,"Type");
						if (assetType == "SWAP")
							GetSwapCustomSchedule(node,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);
						else
							GetEcheancier(node,resetFreq,payFreq,&flowStartDates,&flowEndDates,&resetDates,&paymentDates);

						int compteur;

						// On coupe l'échéancier à la date de paiement qui nous intéresse
						for (int i = paymentDates->GetSize() - 1; i>= 0; i--)
						{
							if ( fabs(paymentDates->Elt(i)-endDate.GetJulian()) < 10. )
							{
								compteur = i;
								i = 0;
							}
						}

						ARM_Vector newflowStartDates (flowStartDates,0,compteur);
						ARM_Vector newflowEndDates (flowEndDates,0,compteur);
						ARM_Vector newresetDates (resetDates,0,compteur);
						ARM_Vector newpaymentDates (compteur+1,paymentDates->Elt(compteur));
						ARM_Vector interestDays	= NULL;

						//JLA interestDays = new ARM_Vector(newflowStartDates.GetSize());
						interestDays = ARM_Vector(newflowStartDates.GetSize());
						for (i = 0; i < newflowStartDates.GetSize(); i++)
							interestDays.Elt(i) = DaysBetweenDates(dayCount, newflowStartDates.Elt(i), newflowEndDates.Elt(i));

						ARM_ReferenceValue cstNotional(notional);

						newIndex	= new ARM_IRIndex( (ARM_INDEX_TYPE)liborType,
														resetFreq,
														K_ZEROCOUPON,
														discountCcy);

						floatingLeg = new ARM_SwapLeg (	&newflowStartDates, 
														&newflowEndDates, 
														&newpaymentDates, 
														&newresetDates,
														&interestDays, 
														NULL, 
														&cstNotional, 
														newIndex, 
														rcvOrPay, 
														spread,
														K_FLOAT_RATE, 
														discountCcy, 
														nxChange,
														decompPricingFlag );

						if (stepUpSpread)
						{
							ARM_ReferenceValue* stepUpSpreadOnReset = ARM_FromStartRefToResetRef(stepUpSpread, floatingLeg);
							delete stepUpSpread;

							floatingLeg->SetVariableSpread(stepUpSpreadOnReset);
							delete stepUpSpreadOnReset;
						}

						floatingLeg->SetStartDate(startDate);
						floatingLeg->SetEndDate((ARM_Date)endDate);

						floatingLeg->CptRealCashFlowDates();

						if (newIndex)
						{
							delete newIndex;
							newIndex = NULL;
						}
						if (flowStartDates)
						{
							delete flowStartDates;
							flowStartDates = NULL;
						}
						if (flowEndDates)
						{
							delete flowEndDates;
							flowEndDates = NULL;
						}
						if (paymentDates)
						{
							delete paymentDates;
							paymentDates = NULL;
						}
						if (resetDates)
						{
							delete resetDates;
							resetDates = NULL;
						}
					}
					delete newIndex;
					newIndex = NULL;
				}
			}

			delete discountCcy;
			discountCcy = NULL;
			
			if (resetCcy)
			{
				delete resetCcy;
				resetCcy = NULL;
			}
			floatingLeg->SetAssetId(assetId);
			if (isQuanto)
			{
				floatingLeg->SetQuantoAsset(1);
			}

			return floatingLeg;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : FloatingAssetLoader::BuildAssetForMemory()");
		} 
	}
};


class MaturityCapAssetLoader : public AbstractAssetParser
{
public:
	MaturityCapAssetLoader()
	{
	}

	~MaturityCapAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_MaturityCapCalculator* newMapCap = NULL;

		ARM_Date start;
		ARM_Date realstart;
		ARM_Date end;
		ARM_Date underlyingDate;
		int dPorS = 1.0;
		int capFloor;
		int resetFreq = K_DEF_FREQ;
		int payFreq = K_DEF_FREQ;
		string indexTerm;
		int basis;
		int intRule;	
		double spread=0.0;
		double nominal;
		double initTri;
		double annuity ;
		double refNotional ;
		ARM_MaturityCapCalculator::ProductMode matCapMode;
		double coef ;
		ARM_ReferenceValue* notional = NULL;
		int resetgap;
		char payCalName[4];
		char resetCalName[4];
		ARM_MaturityCapCalculator::CalibrationMode calibMode;
		ARM_Currency* ccy = NULL;
		//int amortFreq;
		int stubMeth;

		ARM_ReferenceValue* spreads = NULL;

		try
		{
			ccy = GetCcy(node);
			resetgap = GetResetGap(node);

			// Get actual flow 
			// Result : if startdate < AsOf, we search for the first flow after the AsOf				
			//			if startdate > AsOf, we take the first flow
            if ( isEtkOrNot )
			{
			   ARM_SummitFlow theFlow;

			   CCString tradeId = CCString((const char *)GetTradeId().c_str());

			   CCString tradeType = CCString((const char *)GetTradeType().c_str());

			   CCString curveId("");
			
			   RetrieveSummitFlow( tradeId, tradeType, curveId, date, &theFlow);

			   start = theFlow.fStartDate; // actual start date
			   nominal = theFlow.fNotional;// actual nominal
			   initTri = theFlow.fStrike;  // actual initTri
			}
			else
			{
			   start = GetStartDateSpecific(node);
		 	   nominal = GetNotionalSpecific(node);
			   initTri = GetStrikeSpecific(node); 
            }

			realstart = GetStartDate(node); // start date of the input
			end = GetEndDate(node);			// end date of the input
			underlyingDate = GetBLOBDate(node,1);
			
			dPorS = (int) GetPorS(node);
			capFloor = GetCapFloor(node);
			resetFreq = GetResetFreq(node);
			payFreq = GetPayFreq(node);
			char* tmpIndexTerm = GetBLOBIdxTerm(node, ccy,1);
			indexTerm = string(tmpIndexTerm);
			delete [] tmpIndexTerm;
			basis = GetAssetDayCount(node);	
			intRule = GetIntRule(node);
			spread = GetSpread(node); // Let As is because of GP!!!!
			spreads = GetSpreadVariable(node);
			annuity = GetBLOBAmount(node,1);
			refNotional = GetBLOBAmount(node,3);
			matCapMode = (ARM_MaturityCapCalculator::ProductMode) GetBLOBNum(node,1);
			coef = GetBLOBRate(node,1)*100.0;			
			
			GetPayCalendar(node,payCalName);
			GetResetCalendar(node,resetCalName);

			stubMeth = GetStubRule(node);

			int calibModeFlag = GetBLOBNum(node,2);
			switch(calibModeFlag)
			{
			    case K_MCAP_EX_BOUNDARY :
					{
						calibMode = (ARM_MaturityCapCalculator::CalibrationMode) ARM_MaturityCapCalculator::EX_BOUNDARY;
					};
					break;

				case K_MCAP_ATM:
					{
						calibMode = (ARM_MaturityCapCalculator::CalibrationMode) ARM_MaturityCapCalculator::ATM;
					};
					break;

				case K_MCAP_FLAT:
					{
						calibMode = (ARM_MaturityCapCalculator::CalibrationMode) ARM_MaturityCapCalculator::FLAT;
					};
					break;

				default:
					{
						calibMode = (ARM_MaturityCapCalculator::CalibrationMode) ARM_MaturityCapCalculator::EX_BOUNDARY;
					};
					break;

			};

			int nbIterations = 10000;// par default

			// construct the amortisation curve
			char cust[20];
			GetCustom(node, cust);
			if (strcmp(cust,"CUST") != 0)
				notional = GetNotional(node);//amortizing conv à ARM_Curve
			else
				notional = GetNotionalCust(node);

			ARM_Vector* dates = notional->GetDiscreteDates(); // WARNING: These are endDates of the flows
			ARM_Vector* values = notional->GetDiscreteValues();

			// generation du sched 
			int indexType = FromFrequencyToXiborType(payFreq, ccy->GetCcyName()); // a changer en ARM_INDEX_TYPE, ex : EURIBOR6M
			int payTiming = GetPayTiming(node); // meme que reset timing
			ARM_SwapLeg sched(realstart, end, (ARM_INDEX_TYPE)indexType, K_RCV, 0.0, 
							  resetFreq, payFreq, payTiming, payTiming, ccy, intRule);

			ARM_Vector *schedEndDates   = sched.GetFlowEndDates();
			ARM_Vector *schedPayDates   = sched.GetPaymentDates();
			if ( !schedEndDates || !schedPayDates )
			{
				CCString msg((CCString)"Error generating payment schedule");

				throw Exception (__LINE__, __FILE__, ERR_OBJECT_NULL, (char*) msg);
			}		


			ARM_Curve* amortizCurve = NULL;

			if (dates)
			{
				// on construit le vecteur avec les dates de paiement
				ARM_Vector payDates(dates->GetSize());
				int idx;
				for (int i = 0; i < dates->GetSize(); i++)
				{
					idx = schedEndDates->find(dates->Elt(i));
					if (idx == -1)
					{ // enddate du sched sont ajustees mais pas la enddate du deal
						idx = schedPayDates->size() - 1;
					}
					payDates[i] = schedPayDates->Elt(idx);
				}

				std::vector<double> gpDates(payDates.GetSize(),payDates.GetElt());
				gpDates -= date.GetJulian();
				std::vector<double> gpValues(values->GetSize(),values->GetElt());
				amortizCurve = new ARM_Curve(gpDates,gpValues,new ARM::ARM_StepUpRightOpenCstExtrapolDble);
			}
			else
			{
				std::vector<double> gpDates(values->GetSize(),0.0);
				std::vector<double> gpValues(values->GetSize(),values->GetElt());
				amortizCurve = new ARM_Curve(gpDates,gpValues,new ARM::ARM_StepUpRightOpenCstExtrapolDble);
			}


			// search for the first notional because in Summit Schedule, there aren't completes flow values			
			// nominal donne par la startDate des specific data
			int firstFutIdx = schedEndDates->find(start.GetJulian()) + 1; // la start est la end de la periode precedente
			double payLag = schedPayDates->Elt(firstFutIdx) - date.GetJulian();
			double firstNotionalPrevu = amortizCurve->Interpolate(payLag);

			std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(4, false);
			productsToPrice[0] = true;
			
			double nominalrapport = nominal/firstNotionalPrevu; //init nominal par rapport à notional
			double annuityrapport = annuity/refNotional;

			string payCal = (string) payCalName;
			string resetCal = (string) resetCalName;
			
			newMapCap = new ARM_MaturityCapCalculator(date,
													  start,
													  end,
													  underlyingDate,
													  dPorS,
													  capFloor,
													  resetFreq,
													  payFreq, 
													  indexTerm, 
													  basis, 
													  intRule,
													  spread, 
													  nominalrapport,
													  initTri, 
													  annuityrapport, 
													  matCapMode,
													  coef,
													  *amortizCurve, 
													  resetgap, 
													  resetCal,
													  payCal,
													  calibMode, 
													  nbIterations,
													  productsToPrice, 
													  *ccy);

			if (ccy)
			   delete ccy;
			ccy = NULL;

			if (notional)
			   delete notional;
			notional = NULL;

			if (amortizCurve)
			   delete amortizCurve;
			amortizCurve = NULL;

			if (spreads)
				delete spreads;
			spreads = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "MaturityCapAssetLoader::BuildAsset : Failed in constructing Maturity Cap");
		}

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : MaturityCapAssetLoader::BuildAsset() ");
		}

		return newMapCap;
	}
};


// SBGM: In order to take in account the SBGM model we define a macro symbol _SBGM

#define _SBGM

class TarnAssetLoader : public AbstractAssetParser
{
public:
	TarnAssetLoader()
	{
	}

	~TarnAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		return BuildTARNAsset(node, date, assetId, isEtkOrNot, false);
	}

	ARM_Object* BuildTARNAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot, bool oswCalibFlag)
	{
		FloatingAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId, isEtkOrNot));

		int dPorS;
		int cpnDayCount;
		int cpnFreq;
		int cpnTiming;
		int cpnIdxDayCount;
		int cpnResetGap = 0; //For customised deals, field not used
		int intRule;
        int fundFreq;
        int fundDayCount;

		double lifeTimeCap;
		double lifeTimeFloor;
		double coef;
		double meanRev;
		double correl;

		char payCal[4];
		char* fundCcyName;
		char* calibName;

		bool oswFlag;
		bool controlVariableFlag(false);
		bool digSmoothing(false);

		string cpnIdxTerm;
		string cpnResetCal;

		ARM_Date start;
		ARM_Date end;

		ARM_ReferenceValue* strikes			= NULL;
		ARM_ReferenceValue* notionals		= NULL;
		ARM_ReferenceValue* fundNotionals	= NULL;
		
		ARM_ReferenceValue* leverage		= NULL;
		ARM_ReferenceValue* fundSpreads		= NULL;

		ARM_Currency* ccy			= NULL;
		ARM_Currency* fundCcy		= NULL;

		ARM_Curve* strikesCurve     = NULL;
		ARM_Curve* notionalsCurve   = NULL;
		ARM_Curve* leverageCurve    = NULL;
		ARM_Curve* cpnMinCurve		= NULL;
		ARM_Curve* cpnMaxCurve		= NULL;
		ARM_Curve* fundSpreadsCurve = NULL;
		ARM_Curve* fundNomCurve     = NULL;

		ARM_CurveModelParam* mr		= NULL;

		ARM_CorrelMatParam* correlModelParam = NULL;

		//std::vector<double> nbIterations(2, 10000);

		ARM_MarketData_ManagerRep mktDataMgr(date);
		ARM_StringVector mdmKeys;

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(10, false);
		productsToPrice[0] = true;

		ARM_TARNCalculator::CalibrationMode calib;
		
		ARM_TARNCalculator* newTarn = NULL;

		try
		{
			start    = GetStartDate(node);
			end      = GetEndDate(node);
			dPorS    = GetPorS(node);
			intRule  = GetIntRule(node);
			cpnFreq  = GetPayFreq(node);
			ccy      = GetCcy(node);

			GetPayCalendar(node, payCal);

			cpnDayCount    = GetBLOBBasis(node,3);
			cpnIdxDayCount = GetBLOBBasis(node);
			fundDayCount   = GetBLOBBasis(node,2);

			cpnTiming      = GetBLOBTiming(node);

			cpnIdxTerm     = GetBLOBIdxTerm(node,ccy);

			cpnResetCal    = GetBLOBResetCal(node);

			coef		   = GetBLOBRate(node);
			lifeTimeCap    = GetBLOBRate(node,4);
			lifeTimeFloor  = GetBLOBRate(node,5);

			meanRev        = GetBLOBAmount(node, 3);
			correl         = GetBLOBAmount(node, 7);

			fundCcyName    = GetBLOBString(node,3);
			fundCcy        = new ARM_Currency(fundCcyName);

			if (oswCalibFlag)
				oswFlag    = true;
			else
				oswFlag    = (GetBLOBNum(node,3) == 1);

			fundFreq       = GetBLOBNum(node);
			
            if ( fundFreq == 0 )
			{
			   fundFreq = cpnFreq;
			}
			else
            {  
                // convertir en ARM
				
                char frqStr[5];
				sprintf(frqStr, "%i", fundFreq);
				fundFreq = FromSummitFreqToARMFreq(frqStr);
			}
	
			// strikes
			strikes        = GetBLOBSchedAmount(node, 1, 1);
			ARM_ReferenceValue* newStrikes;

			if (strikes->GetDiscreteDates() != NULL)
				newStrikes = ARM_FromStartRefToResetRef(strikes,newSwapLeg);
			else
				newStrikes = (ARM_ReferenceValue*) strikes->Clone(); 

			strikesCurve   = FromRefValueToARM_Curve(newStrikes, date, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());
			strikesCurve->GetOrdinates() /= 100.0;

			// leverage
			leverage = GetBLOBSchedAmount(node, 1, 2); // RefValue with Schedule1/StartDate and schedule1/Amount2
			
            if (leverage == NULL)
				leverage = new ARM_ReferenceValue(coef);
			
            leverageCurve = FromRefValueToARM_Curve(leverage, date, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());

			std::vector<double> cpnMinAbs(1,0.0),cpnMinOrd(1,0.0), cpnMaxAbs(1,0.0), cpnMaxOrd(1,1000.0);

			cpnMinCurve = new ARM_Curve(cpnMinAbs, cpnMinOrd, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());
			cpnMaxCurve = new ARM_Curve(cpnMaxAbs, cpnMaxOrd, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());

			// fundSpreads
			fundSpreads = GetBLOBSchedAmount(node, 2, 1); // RefValue with Schedule2/StartDate and schedule2/Amount1
			if (fundSpreads == NULL)
				fundSpreads = new ARM_ReferenceValue( GetBLOBAmount(node, 6) );
			fundSpreadsCurve   = FromRefValueToARM_Curve(fundSpreads, date, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());
			fundSpreadsCurve->GetOrdinates() /= 100.0;

			// funding nominal
			fundNotionals = GetFundNotionalWithPayDates(node);
			fundNomCurve  = FromRefValueToARM_Curve(fundNotionals, date, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());

			// calibration
			calibName = GetBLOBString(node,4);
			if (!strcmp(calibName, "NO"))
				calib = ARM_TARNCalculator::NOCalib;
			else if (!strcmp(calibName, "RF"))
				calib = ARM_TARNCalculator::RFStrike;
			else if (!strcmp(calibName, "ATM"))
				calib = ARM_TARNCalculator::ATMStrike;
			else 
				calib = ARM_TARNCalculator::ExerStrike;

			//Nb Iterations
			double nbIter = GetNbIterations();
			std::vector<double> nbIterations(3, nbIter);

			//Coupon Reset Gap:
			string sname = GetTarnNodeName(node, "Custom");
			bool isCust = false;
			std::vector<double>& resetDates;
			if (sname == "CUST")
			{
				isCust = true;
				//Reset dates
				resetDates = GetTarnCustomDates(node);
			}
			else //Only one coupon reset gap
			{
				cpnResetGap = GetBLOBResetGap(node);
				resetDates = NULL;
			}
			
			//Nominaux
			notionals      = GetNotionalWithPayDates(node, isCust);
			notionalsCurve = FromRefValueToARM_Curve(notionals, date, new ARM::ARM_StepUpLeftOpenCstExtrapolDble());

			// At this moment the fees are not handled by SUMMIT.
			std::vector<double> dates(1,start.GetJulian());
			std::vector<double> values(1,0.0);
			ARM_Curve fees(dates, values, new ARM::ARM_StepUpLeftOpenCstExtrapolDble);

			//Started deal: AsoFDate > FirstResetDate
			std::vector<double>& pastFixings = GetPastFixings(node, date);
			std::vector<double>& pastStart   = GetPastStart(node, date);
			
            ARM_Date futureStart = start;

            if ( pastFixings != NULL ) //Past fixings
			{
				int size = pastFixings->size();

				double currCap = 0;
				double currFloor = 0;
				double currFixing;
				double currStrike;
				double currLeverage;
				
                currCap = lifeTimeCap;
				currFloor = lifeTimeFloor;

				ARM_Date currStart;
				ARM_Date currEnd;
				
                double cvg;

				for (int i = 0; i < size ; i++)
				{
					currStart = ARM_Date(pastStart->Elt(i));
					currEnd   = ARM_Date(pastStart->Elt(i));
					
                    

                    if ( cpnFreq == K_ANNUAL )
					{
						currEnd.AddMonths(12);
						currEnd.AdjustToBusDate(ccy->GetCcyName(), (long)intRule);	
					}

					if ( cpnFreq == K_SEMIANNUAL )
					{
						currEnd.AddMonths(6);
						currEnd.AdjustToBusDate(ccy->GetCcyName(), (long)intRule);	
					}

					if ( cpnFreq == K_QUARTERLY )
					{
						currEnd.AddMonths(3);
						currEnd.AdjustToBusDate(ccy->GetCcyName(), (long)intRule);	
					}

					if ( cpnFreq == K_MONTHLY )
					{
						currEnd.AddMonths(1);
						currEnd.AdjustToBusDate(ccy->GetCcyName(), (long)intRule);			
					}
					
					cvg = CountYears(cpnDayCount, currStart.GetJulian(), currEnd.GetJulian());

					currStrike   = newStrikes->GetDiscreteValues()->Elt(i)/100;
					currLeverage = leverage->GetDiscreteValues()->Elt(i);
					currFixing   = pastFixings->Elt(i);

                    double coupon = (currStrike-currLeverage*currFixing)*cvg;

                    if ( coupon < 0 ) // Coupons are floored to 0
                       coupon = 0.0;
                    
					currCap   = currCap-coupon;   // Update with cvg.
					currFloor = currFloor-coupon; // Update with cvg.
				
                    // Update the real start date

                    futureStart = currEnd;
                }
			
                lifeTimeCap	  = currCap;
				lifeTimeFloor = currFloor;
				
                if ( lifeTimeCap < 0 )
				   lifeTimeCap = 0.0;
				
                if ( lifeTimeFloor < 0 )
				   lifeTimeFloor = 0.0;

				delete pastFixings;
				pastFixings = NULL;

                delete pastStart;
				pastStart = NULL;
			}


#ifndef _SBGM // case of old TARN
			
            if(!isCust)
			{
				newTarn = new ARM_TARNCalculator(date,
												futureStart,
												end,
												*strikesCurve,
												dPorS,
												cpnDayCount,
												cpnFreq,
												cpnTiming,
												cpnIdxTerm,
												cpnIdxDayCount,
												cpnResetCal,
												payCal,
												cpnResetGap,
												intRule,
												*leverageCurve,
												cpnMinCurve,
												cpnMaxCurve,
												lifeTimeCap,
												lifeTimeFloor,
												*fundSpreadsCurve,
												fundFreq,
												fundDayCount,
												*notionalsCurve,
												fees,
												nbIterations,
												calib,
												calib,
												oswFlag,
												controlVariableFlag,
												digSmoothing,
												productsToPrice,
												*ccy,
												*fundCcy,
												*ccy,
												*ccy,
												*fundCcy,
												*fundNomCurve,
												&mdmKeys);
			}
			else
			{
				newTarn = new ARM_TARNCalculator(date,
												futureStart,
												end,
												*strikesCurve,
												dPorS,
												cpnDayCount,
												cpnFreq,
												cpnTiming,
												cpnIdxTerm,
												cpnIdxDayCount,
												cpnResetCal,
												payCal,
												cpnResetGap,
												intRule,
												*leverageCurve,
												lifeTimeCap,
												lifeTimeFloor,
												*fundSpreadsCurve,
												fundFreq,
												fundDayCount,
												*notionalsCurve,
												fees,
												nbIterations,
												calib,
												calib,
												oswFlag,
												controlVariableFlag,
												digSmoothing,
												productsToPrice,
												*ccy,
												*fundCcy,
												*ccy,
												*ccy,
												*fundCcy,
												*fundNomCurve,
												&mdmKeys,
												NULL,
												true,
												*resetDates);
			}

			// Set Mean Reversion param
			std::vector<double> mrsValues(1, meanRev);
			std::vector<double> mrsTimes (1, 0.0);
			mr = new ARM_CurveModelParam(ARM_ModelParamType::MeanReversion, 
										 &mrsValues, 
										 &mrsTimes, 
										 "MRS");
			newTarn->SetMRS(mr);
			

			// Set correl calib param
			
			// Buid a DateStrip to have business days in times
			ARM_DateStrip resetLegSched(futureStart, 
										end,
										cpnFreq,
										cpnDayCount,
										cpnResetCal.c_str(),
										K_MOD_FOLLOWING,
										intRule,
										K_SHORTSTART,
										cpnResetGap,
										cpnFreq,
										GETDEFAULTVALUE,
										"DEFAULT",
										cpnTiming);

			// Conversion from ARM_Vector To ARM_RealVector
			std::vector<double>& correlTimes = NULL;
			if (isCust)
				correlTimes = resetDates;		
			else
				correlTimes = resetLegSched.GetResetDates();
			
			*correlTimes -= date.GetJulian(); // times will be deleted in DateStrip destructor

			// Conversion from ARM_GP_Matrix To ARM_RealVector
			int n = correlTimes->size();
			ARM_GP_Matrix* trigoMatrix = TrigoMatrix(n,correl);
			std::vector<double>& correlValues = new std::vector<double>(n*n);

			size_t i = 0,j = 0;
			for (i = 0; i < n; ++i)
				for (j = 0; j < n; ++j)
					(*correlValues)[i*n+j] = (*trigoMatrix)(i,j);

			correlModelParam = static_cast<ARM_CorrelMatParam*>(ARM_ModelParamFactory.Instance()->CreateModelParam(
				ARM_ModelParamType::Correlation,
				correlValues,
				correlTimes));

			newTarn->SetCorrel(correlModelParam);

#else // new TARN

// -----> SBGM case: See SBGM constructor in tarncalculator.cpp

            // Initialize to defaults

            string genType1   = "MRGK5";
	        string genType2   = "Sobol";
	        int firstNbTimes  = 0;
	        int firstNbDims   = 0;
	        string pathScheme = "Incremental";
	        string pathOrder  = "BucketOrder";
	        bool antithetic   = (bool) 1;
			bool globalCap = GetBLOBNum(node,4) ? true : false;

			nbIterations[2] = 10000;

            if (!isCust)
			{
				newTarn = new ARM_TARNCalculator(date,
												futureStart,
												end,
												*strikesCurve,
												dPorS,
												cpnDayCount,
												cpnFreq,
												cpnTiming,
												cpnIdxTerm,
												cpnIdxDayCount,
												cpnResetCal,
												payCal,
												cpnResetGap,
												intRule,
												*leverageCurve,
												*cpnMinCurve,
												*cpnMaxCurve,
												lifeTimeCap,
												globalCap,
												lifeTimeFloor,
												*fundSpreadsCurve,
												fundFreq,
												fundDayCount,
												*notionalsCurve,
												fees,
												nbIterations,
                                                ARM_PricingModelType::SBGM,
												calib,
												calib,
												oswFlag,
												controlVariableFlag,
												digSmoothing,
												false,
                                                // new Fields
                                                genType1,
	                                            genType2,
	                                            firstNbTimes,
	                                            firstNbDims,
	                                            pathScheme,
	                                            pathOrder,
	                                            antithetic,
                                                // End new fields
												productsToPrice,
												*ccy,
												*fundCcy,
												*ccy,
												*ccy,
												*fundCcy,
												*fundNomCurve,
												&mdmKeys);
			}
			else
			{
				newTarn = new ARM_TARNCalculator(date,
												futureStart,
												end,
												*strikesCurve,
												dPorS,
												cpnDayCount,
												cpnFreq,
												cpnTiming,
												cpnIdxTerm,
												cpnIdxDayCount,
												cpnResetCal,
												payCal,
												cpnResetGap,
												intRule,
												*leverageCurve,
												*cpnMinCurve,
												*cpnMaxCurve,
												lifeTimeCap,
												globalCap,
												lifeTimeFloor,
												*fundSpreadsCurve,
												fundFreq,
												fundDayCount,
												*notionalsCurve,
												fees,
												nbIterations,
                                                ARM_PricingModelType::SBGM,
												calib,
												calib,
												oswFlag,
												controlVariableFlag,
												digSmoothing,
												false,
                                                // new Fields
                                                genType1,
	                                            genType2,
	                                            firstNbTimes,
	                                            firstNbDims,
	                                            pathScheme,
	                                            pathOrder,
	                                            antithetic,
                                                // End new fields
												productsToPrice,
												*ccy,
												*fundCcy,
												*ccy,
												*ccy,
												*fundCcy,
												*fundNomCurve,
												&mdmKeys,
												NULL,
												true,
												*resetDates);
			}

            // Set the betaCorrel and the hump
            ARM_CurveModelParam hump;
            ARM_CurveModelParam betaCorrel;

            // The hump is defaulted for now
            const double theHUMP_DEFAULT_VALUE       = 1.;
            std::vector<double>       hValues(1, theHUMP_DEFAULT_VALUE);
	        std::vector<double>       hTimes (1, 0.0);

		    hump = ARM_CurveModelParam(ARM_ModelParamType::Hump, 
								       &hValues, 
								       &hTimes, 
								       "HUMP");

            newTarn->SetHump(&hump);


            // Set the Beta Correl
            std::vector<double>       bcValues(1, meanRev);
	        std::vector<double>       bcTimes (1, 0.0);

		    betaCorrel = ARM_CurveModelParam(ARM_ModelParamType::BetaCorrelation, 
									         &bcValues, 
									         &bcTimes, 
									         "BETACORREL");



            newTarn->SetBetaCorrel(&betaCorrel);

#endif

			if (ccy)
			   delete ccy;
			ccy = NULL;

			if (fundCcy)
				delete fundCcy;
			fundCcy = NULL;

			if (fundCcyName)
				delete fundCcyName;
			fundCcyName = NULL;

			if (strikesCurve)
			   delete strikesCurve;
			strikesCurve = NULL;

			if (notionalsCurve)
			   delete notionalsCurve;
			notionalsCurve = NULL;

			if (leverageCurve)
				delete leverageCurve;
			leverageCurve = NULL;

			if (fundSpreadsCurve)
				delete fundSpreadsCurve;
			fundSpreadsCurve = NULL;

			if (fundNomCurve)
				delete fundNomCurve;
			fundNomCurve = NULL;

			if (strikes)
				delete strikes;
			strikes = NULL;

			if (newStrikes)
				delete newStrikes;
			newStrikes = NULL;

			if (newSwapLeg)
				delete newSwapLeg;
			newSwapLeg = NULL;

			if (notionals)
				delete notionals;
			notionals = NULL;

			if (leverage)
				delete leverage;
			leverage = NULL;

			if (fundSpreads)
				delete fundSpreads;
			fundSpreads = NULL;

#ifndef _SBGM
			
            if (trigoMatrix)
			   delete trigoMatrix;
			trigoMatrix = NULL;

			if (correlValues)
			   delete correlValues;
			correlValues = NULL;

#else

#endif

			if (calibName)
			   delete calibName;
			calibName = NULL;

			if (mr)
			   delete mr;
			mr = NULL;
		}

		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		}

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : TarnAssetLoader::BuildAsset() ");
		}

		return newTarn;
	}
};


class CaptionAssetLoader : public AbstractAssetParser
{
public:
	CaptionAssetLoader()
	{
	}

	~CaptionAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_Date startDate, endDate;
		string cpnIdxTerm;
		int payRec, capFloor;
		ARM_ReferenceValue* cpnProfile = NULL;
		string fundIdxTerm;
		int notifDays, nonCall;
		ARM_ReferenceValue* exerciseStyle = NULL;
		ARM_ReferenceValue* notional = NULL;
		int cpnDayCount, cpnResetTiming;
		char cpnPayCal[7];
		char cpnResetCal[7];
		ARM_ReferenceValue* cpnSpread = NULL;
		int fundDayCount;
		char fundPayCal[7];
		char fundResetCal[7];
		ARM_ReferenceValue* fundSpread = NULL;
		
		ARM_CaptionCalculator* captionCalculator = NULL;

		try
		{
			// CAP FLOOR DATA
			ARM_Currency* ccy = GetCcy(node);
			startDate = GetStartDate(node);
			endDate = GetEndDate(node);
			payRec = GetCaptionUnderPorS(node);
			capFloor = GetCapFloor(node);

			// CAPPED LEG DATA
			int index = GetCaptionIndex(node, ccy);
			cpnIdxTerm = YearTermToStringMatu(1.0/(double)FromLiborTypeToFrequency(index));
			cpnDayCount = GetCaptionCpnDayCount(node);
			cpnResetTiming = GetResetTiming(node);
			GetResetCalendar(node, cpnResetCal);
			GetPayCalendar(node, cpnPayCal);
			cpnSpread = GetCaptionCpnSpread(node);
			cpnProfile = GetCapStrikes(node);
			*cpnProfile /= 100.0;
			if (cpnProfile->GetDiscreteValues()->GetSize() == 1)
			{
				*cpnProfile += cpnSpread->GetDiscreteValues()->Elt(0);
			}
			else
			{
				int decompPricingFlag = GetDecompFlag(node);
				int nxChange = GetNotionalExchangeFlag(node);
				int stubRule = GetStubRule(node);
				int adjStartDate = 0;

				for (int i=0; i<cpnProfile->GetDiscreteValues()->GetSize(); i++)
				{
					cpnProfile->GetDiscreteValues()->Elt(i) += cpnSpread->Interpolate(cpnProfile->GetDiscreteDates()->Elt(i));
				}

				// strikes will be interpolated on reset dates
				// while our reference value is set on start dates
				// just update the reference value with corresponding reset dates
				ARM_SwapLeg tmpLeg( startDate, endDate, (ARM_INDEX_TYPE)index, payRec, 0.0,
									GetResetFreq(node), GetPayFreq(node), GetResetTiming(node),
									GetPayTiming(node), ccy, GetIntRule(node), GetResetGap(node),
									cpnResetCal, cpnPayCal, decompPricingFlag, nxChange,
									stubRule, "NULL", adjStartDate);

				cpnProfile = ARM_FromStartRefToResetRef(cpnProfile, &tmpLeg);
			}

			// FUNDING LEG DATA
			fundIdxTerm = cpnIdxTerm;
			fundDayCount = GetCaptionFundDayCount(node);
			GetResetCalendar(node, fundResetCal);
			GetPayCalendar(node, fundPayCal);
			// spread information is stored in the XML "caption" block (constant spread case only), or in "asset" block 
			fundSpread = GetCaptionFundSpread(node);
			if (fundSpread == NULL)
			{
				fundSpread = GetSpreadVariable(node);
				if (fundSpread == NULL)
				{
					fundSpread = new ARM_ReferenceValue(0.0);
				}
				else
				{
					// spreads will be interpolated on reset dates
					// while our reference value is set on start dates
					// just update the reference value with corresponding reset dates
					int decompPricingFlag = GetDecompFlag(node);
					int nxChange = GetNotionalExchangeFlag(node);
					int stubRule = GetStubRule(node);
					int adjStartDate = 0;

					ARM_SwapLeg fundingLeg( startDate, endDate, (ARM_INDEX_TYPE)index, -payRec, 0.0,
											GetResetFreq(node), GetPayFreq(node), GetResetTiming(node),
											GetPayTiming(node), ccy, GetIntRule(node), GetResetGap(node),
											fundResetCal, fundPayCal, decompPricingFlag, nxChange,
											stubRule, "NULL", adjStartDate);

					fundSpread = ARM_FromStartRefToResetRef(fundSpread, &fundingLeg);
				}
				*fundSpread /= 100.;
			}

			// OPTION DATA
			notifDays = GetCaptionNotifDays(node);
			ARM_Vector* exerciseDates = GetCaptionExerciseDates(node, date, (ARM_INDEX_TYPE)index,
																cpnResetTiming, K_ARREARS, ccy,
																notifDays, cpnResetCal, cpnPayCal);

			ARM_Vector* fees = NULL;
			// fees might be null, with no data in the XML file
			// in that case GetFees() returns an error !
			try
			{
				fees = GetFees(node);
			}
			catch(...)
			{
				delete fees;
				fees = new ARM_Vector(exerciseDates->GetSize(), 0.0);
			}
			*fees /= 100.;

			exerciseStyle = new ARM_ReferenceValue(exerciseDates, fees);
			notional = GetNotional(node);
			nonCall = 0; // TMP	
			
			captionCalculator = new ARM_CaptionCalculator ( date, startDate, endDate, cpnIdxTerm,
															payRec, capFloor, *RefValueToCurve(*cpnProfile, date.GetJulian()),
															fundIdxTerm, notifDays, nonCall,
															*RefValueToCurve(*exerciseStyle, date.GetJulian()),
															*RefValueToCurve(*notional, date.GetJulian()),
															cpnDayCount, cpnResetTiming,
															cpnResetCal, cpnPayCal,
															*RefValueToCurve(*cpnSpread, date.GetJulian()),
															fundDayCount, fundResetCal, fundPayCal,
															*RefValueToCurve(*fundSpread, date.GetJulian()) );

			int PorS = GetPorS(node);
			captionCalculator->SetPorS(PorS);

			delete exerciseStyle;
			delete cpnProfile;
			delete fundSpread;
			delete cpnSpread;
			delete notional;
		}
		catch(Exception& x)
		{
	        x.DebugPrint();

			throw x;
		}

		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : CaptionAssetLoader::BuildAsset() ");
		}

		return captionCalculator;
	}
};

//! \brief loader for a swap
// BuildAsset return one swapleg
// To get the whole swap, call BuildAsset for each swapleg
class SwapAssetLoader : public AbstractAssetParser
{
public:
	SwapAssetLoader()
	{
	}

	~SwapAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_SwapLeg* newSwapLeg = NULL;
		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, *item = NULL;
			node->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				theNode->Release();
				theNode = NULL;
				if (strcmp((const char*)(_bstr_t)resultat,"FLO") == 0) //Jambe flottante
				{
					FloatingAssetLoader loader;
					newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId, isEtkOrNot));
				}
				else //Jambe fixe
				{
					FixAssetLoader loader;
					newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId, isEtkOrNot));
				}
				
				ARM_ReferenceValue* primes = GetPrimes(node);
				if (primes)
				{
					newSwapLeg->SetFee(primes);
				}			
			}

			newSwapLeg->SetAssetId(assetId);

			return newSwapLeg;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : SwapAssetLoader::BuildAsset()");
		} 
	}

	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate,double notional,int isEtkOrNot = 1)
	{
		ARM_SwapLeg* newSwapLeg = NULL;
		try
		{
			assetId = GetCCStringFromXMLNode(node,"dmAssetId");
			MSXML2::IXMLDOMNode * listItem = NULL, * theNode = NULL, *item = NULL;
			node->selectSingleNode(_bstr_t((const char *)"INTEREST_FixFloat"), &theNode);
			if (theNode!=NULL)
			{
				BSTR resultat = NULL;
				theNode->get_text(&resultat);
				theNode->Release();
				theNode = NULL;
				if (strcmp((const char*)(_bstr_t)resultat,"FLO") == 0) //Jambe flottante
				{
					FloatingAssetLoader loader;
					newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));
				}
				else //Jambe fixe
				{
					FixAssetLoader loader;
					newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));
				}
				
				ARM_ReferenceValue* primes = GetPrimes(node);
				if (primes)
				{
					newSwapLeg->SetFee(primes);
				}			
			}

			newSwapLeg->SetAssetId(assetId);

			return newSwapLeg;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : SwapAssetLoader::BuildAsset()");
		} 
	}
};


//! \brief loader for a cap or a floor
class CapFloorAssetLoader : public AbstractAssetParser
{
public:
	CapFloorAssetLoader()
	{
	}

	~CapFloorAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		SwapAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date, assetId, isEtkOrNot));

		ARM_CapFloor* newCapFloor = NULL;

		ARM_Barrier* newCFBarrier = NULL;

		MSXML2::IXMLDOMNodeList * resultList = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL;

		int capFloor;
		int liborType;
		int resetFreq = K_DEF_FREQ;
		int payFreq = K_DEF_FREQ;
		char RollDate[5];
		double dPorS = 1.0;
		bool withBarrier = false;
		int pOrS = 1;

		ARM_ReferenceValue* strikes = NULL;
		ARM_ReferenceValue* barriers = NULL;
		
		ARM_Currency* discountCcy = NULL;
		
		// Global Cap
		ARM_ReferenceValue* globalCapSpreads = NULL;
		ARM_ReferenceValue* globalCapFixedRates = NULL;
		ARM_ReferenceValue* globalCapBarriers = NULL;
		ARM_ReferenceValue* globalCapPastFixings = NULL;
		double globalCapNumerator = 1.0;
		double globalCapDenominator = 1.0;
		double nbIter = 0.0;

		try
		{
			// CF a barriere : il y a en plus une formule
			// de type : CUO(Index,Strike,Barrier)
			char ProductName[30];
			GetFormula(node, ProductName);
			
			if ( (strncmp(ProductName, "CUO", 3) == 0)
				|| (strncmp(ProductName, "CUI", 3) == 0)
				|| (strncmp(ProductName, "PDO", 3) == 0)
				|| (strncmp(ProductName, "PDI", 3) == 0) )
			{
				withBarrier = true;
				strikes = GetRefValueFromFormula(node, 2);
				barriers = GetRefValueFromFormula(node, 3);
			}
			else
			{
				strikes = GetCapStrikes(node);
			}
			dPorS = GetPorS(node);
			capFloor = GetCapFloor(node);
			discountCcy = GetCcy(node);
			resetFreq = GetResetFreq(node);
			payFreq = GetPayFreq(node);
			GetResetRollDate(node,RollDate);
			char payRolldate[15];
			GetPayRollDate(node,payRolldate);

			bool isQuanto = false;
			
			liborType = GetFormulaIndexType(node, discountCcy);
			if (liborType != K_LIVRET_A)
			{
				if (liborType != -1)
				{
					isQuanto = true;
				}
				else
				{
					liborType = GetIndexType(node);
				}
			}

			if ( (resetFreq != K_DAILY) && (strcmp(payRolldate,RollDate)!=0))
			{
				CCString msg((CCString)"ResetRollDate must be equal to PayRollDate in XML parsing for getting Cap \n");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
			}

			if (strikes->GetDiscreteDates() != NULL)
				strikes = ARM_FromStartRefToResetRef(strikes,newSwapLeg);
			
			// Recup des Primes et Fees
			ARM_ReferenceValue* Premium = GetPrimes(node);

			pOrS = GetPorS(node);

			if(strcmp(ProductName, "GLOBALCAP") == 0)
			{
				GetGlobalCapInfo(node, globalCapSpreads, globalCapFixedRates, globalCapBarriers, newSwapLeg, globalCapNumerator, globalCapDenominator, nbIter, date);
				liborType = GetGlobalCapIndexType(node, discountCcy);

				globalCapPastFixings = (ARM_ReferenceValue*)GetRefValSpreadFirstFixing(node,date);
				if (globalCapPastFixings)
					*globalCapPastFixings /= 100.;

				newCapFloor = (ARM_GlobalCap*) new ARM_GlobalCap(newSwapLeg, capFloor, strikes->GetDiscreteValues()->Elt(0), globalCapSpreads, globalCapFixedRates, globalCapBarriers, globalCapNumerator/globalCapDenominator, nbIter, globalCapPastFixings, pOrS);
			}
			else
			{
				newCapFloor = new ARM_CapFloor(newSwapLeg, capFloor, strikes);
				newCapFloor->SetPorS(pOrS);
			}

			newCapFloor->SetPorS(dPorS);

			if (Premium)
				newCapFloor->SetFee(Premium);

			if (isQuanto)
				newCapFloor->SetQuantoAsset(1);

			if (strikes)
				delete strikes;
			strikes = NULL;

			// Cas Cap a barriere : on construit un objet ARM_Barrier
			if (withBarrier)
			{
				ARM_ExerciseStyle* barrierStyle = new ARM_ExerciseStyle(newCapFloor->GetSwapLeg()->GetResetDates(),
																		newCapFloor->GetSwapLeg()->GetFlowStartDates());
				
				ARM_ReferenceValue* newBarriers = NULL;
				
				if (barriers->GetDiscreteDates() != NULL)
					barriers = ARM_FromStartRefToResetRef(barriers, newSwapLeg);

				int upDown = K_UP;
				int inOut = K_IN;

				if ( (strncmp(ProductName, "PDO", 3) == 0)
					|| (strncmp(ProductName, "PDI", 3) == 0) )
				{
					upDown = K_DOWN;
				}

				if ( (strncmp(ProductName, "CUO", 3) == 0)
					|| (strncmp(ProductName, "PDO", 3) == 0) )
				{
					inOut = K_OUT;
				}

				newCFBarrier = new ARM_Barrier( newCapFloor, newCapFloor, 
												barrierStyle, barriers,
												upDown, inOut );

				if (newCapFloor)
					delete newCapFloor;
				newCapFloor = NULL;

				if (barrierStyle)
					delete barrierStyle;
				barrierStyle = NULL;

				if (barriers)
					delete barriers;
				barriers = NULL;

				return newCFBarrier;
			}

			if (newSwapLeg)
				delete newSwapLeg;
			newSwapLeg = NULL;
			
			if (globalCapSpreads)
				delete globalCapSpreads;
			globalCapSpreads = NULL;

			if (globalCapFixedRates)
				delete globalCapFixedRates;
			globalCapFixedRates = NULL;

			if (globalCapBarriers)
				delete globalCapBarriers;
			globalCapBarriers = NULL;

			if (globalCapPastFixings)
				delete globalCapPastFixings;
			globalCapPastFixings = NULL;

			return newCapFloor;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : CapFloorAssetLoader::BuildAsset()");
		} 
	}

	// Méthode spécifique pour contruire les memoryCap depuis le loader de Cap
	ARM_Object* BuildAssetForMemory(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,const ARM_Date& endDate,double notional,int isEtkOrNot = 1)
	{
		SwapAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAssetForMemory(node, date, assetId, endDate, notional, isEtkOrNot));

		ARM_CapFloor* newCapFloor = NULL;

		MSXML2::IXMLDOMNodeList * resultList = NULL;
		MSXML2::IXMLDOMNode * listItem = NULL;

		int capFloor;
		int liborType;
		int resetFreq = K_DEF_FREQ;
		char RollDate[5];
		double dPorS = 1.0;

		ARM_ReferenceValue* strikes = NULL;
		ARM_Currency* discountCcy = NULL;
		ARM_Date end;

		try
		{
			end = GetEndDate(node);

			strikes = GetCapStrikes(node);
			dPorS = GetPorS(node);
			capFloor = GetCapFloor(node);
			discountCcy = GetCcy(node);
			resetFreq = GetResetFreq(node);
			GetResetRollDate(node,RollDate);
			char payRolldate[15];
			GetPayRollDate(node,payRolldate);

			bool isQuanto = false;
			liborType = GetFormulaIndexType(node, discountCcy);
			if (liborType != K_LIVRET_A)
			{
				if (liborType != -1)
				{
					isQuanto = true;
				}
				else
				{
					liborType = GetIndexType(node);
				}
			}

			if ( (resetFreq != K_DAILY) && (strcmp(payRolldate,RollDate)!=0))
			{
				CCString msg((CCString)"ResetRollDate must be equal to PayRollDate in XML parsing for getting Cap \n");

				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						 (const char *) msg);
			}
				
			ARM_ReferenceValue* newStrikes;

			if (strikes->GetDiscreteDates() != NULL)
				newStrikes = ARM_FromStartRefToResetRef(strikes,newSwapLeg);
			else
				newStrikes = (ARM_ReferenceValue*) strikes->Clone(); 


			newCapFloor = new ARM_CapFloor(newSwapLeg, capFloor, strikes);

			newCapFloor->SetPorS(dPorS);

			// Recup des Primes et Fees
			if (end == endDate)
			{
				ARM_ReferenceValue* Premium = GetPrimes(node);
				if (Premium)
				{
					newCapFloor->SetFee(Premium);
					delete Premium;
				}
			}

			if (isQuanto)
				newCapFloor->SetQuantoAsset(1);

			if (strikes)
				delete strikes;
			strikes = NULL;

			if (newStrikes)
				delete newStrikes;
			newStrikes = NULL;

			return newCapFloor;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : CapFloorAssetLoader::BuildAsset()");
		} 
	}
};

class DigitalAssetLoader : public AbstractAssetParser
{
public:
	DigitalAssetLoader()
	{
	}

	~DigitalAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node, const ARM_Date& date, string& assetId, int isEtkOrNot = 1)
	{
		ARM_ReferenceValue* strikes = NULL;
		ARM_ReferenceValue* payOffs = NULL;
		MSXML2::IXMLDOMNode* theNode=NULL;
		MSXML2::IXMLDOMNodeList* nodeList = NULL;
		long nb=0;
		ARM_Vector* vDates   = NULL;
		ARM_Vector* vStrikes = NULL;
		ARM_Vector* vPayOffs = NULL;
		int capFloor;
		double dPorS = 1.0;

		ARM_Digital* newDigital = NULL;
		ARM_Currency* discountCcy = NULL;
		ARM_ReferenceValue* Premium = NULL;

		SwapAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date, assetId, isEtkOrNot));

		capFloor = GetCapFloor(node);
		try
		{

			node->selectNodes(_bstr_t("Formula/FORMULA"), &nodeList);
			if(nodeList != NULL)
			{
				nodeList->get_length(&nb);
				if(nb == 0)
				{
					CCString msg((CCString)"Invalid XML string for getting Digital Strikes \n");

					throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (char*) msg);
				}
				discountCcy = GetCcy(node);
				Premium = GetPrimes(node);
				dPorS = GetPorS(node);
				vDates = new ARM_Vector(nb);
				vStrikes = new ARM_Vector(nb);
				vPayOffs = new ARM_Vector(nb);
				for(int i = 0; i < nb; i++)
				{
					nodeList->get_item(i,&theNode);
					vDates->Elt(i) = (GetDateFromXMLNode(theNode,"Date")).GetJulian();
					char tmpStrike[20];
					char tmpPayOff[20];

					GetNthArgInFormula(theNode,2,tmpStrike,1);
					vStrikes->Elt(i) = atof(tmpStrike);
					GetNthArgInFormula(theNode,3,tmpPayOff,1);
					vPayOffs->Elt(i) = atof(tmpPayOff);

					if (theNode)
					theNode->Release();
				}
				nodeList->Release();

				if(nb == 1)
				{
					strikes = new ARM_ReferenceValue(vStrikes->Elt(0));
					payOffs = new ARM_ReferenceValue(vPayOffs->Elt(0));
				}
				else
				{
					strikes = new ARM_ReferenceValue(vDates,vStrikes);
					strikes->SetCalcMethod(K_STEPUP_LEFT);
					payOffs = new ARM_ReferenceValue((ARM_Vector*)vDates->Clone(),vPayOffs);
					payOffs->SetCalcMethod(K_STEPUP_LEFT);
				}

				if (strikes->GetDiscreteDates() != NULL)
					strikes = ARM_FromStartRefToResetRef(strikes, newSwapLeg);

				if (payOffs->GetDiscreteDates() != NULL)
					payOffs = ARM_FromStartRefToResetRef(payOffs, newSwapLeg);

				newDigital = new ARM_Digital(newSwapLeg, capFloor, strikes, 0.01, -0.01, payOffs);

				newDigital->SetPorS(dPorS);

				int libortype = GetFormulaIndexType(node, discountCcy);

				if ( (libortype != -1) && (libortype != K_LIVRET_A) )
					newDigital->SetQuantoAsset(1);

				if(Premium)
					newDigital->SetFee(Premium);

				delete newSwapLeg;
				delete strikes;
				delete payOffs;
				delete Premium;
			}

			return newDigital;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw x;
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : DigitalAssetLoader::BuildAsset()");
		} 

	}
};



//! \brief loader for a Memory cap or a floor
class MemoryCapFloorAssetLoader : public AbstractAssetParser
{
public:
	MemoryCapFloorAssetLoader()
	{
	}

	~MemoryCapFloorAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_StdPortfolio* pfCap = NULL;
		vector<ARM_Security*> assets;
		vector<double> weights;
		vector<double> prices;
		vector<double> vegas;

		SwapAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId, isEtkOrNot));

		ARM_Date startDate = newSwapLeg->GetEndDateNA();
		ARM_Date endDate = newSwapLeg->GetEndDateNA();
		int resetFreq = newSwapLeg->GetIRIndex()->GetResetFrequency();
		ARM_ReferenceValue* varNotional = newSwapLeg->GetAmount();
		ARM_Vector* paymentDates = newSwapLeg->GetPaymentDates();

		CapFloorAssetLoader cfloader;
		double notional;
		ARM_Date tmpEndDate(endDate);

		for (int i=paymentDates->GetSize()-1; i>=0; i--)
		{
			notional = varNotional->CptReferenceValue(paymentDates->Elt(i));
			ARM_CapFloor* newCap = dynamic_cast<ARM_CapFloor*>(cfloader.BuildAssetForMemory(node, date,assetId, endDate, notional, isEtkOrNot));
			assets.push_back(dynamic_cast<ARM_Security*>(newCap));
			weights.push_back(1.0);
			prices.push_back(1.0);
			vegas.push_back(0.0);
			endDate = tmpEndDate;
			endDate.AddPeriodMult(-resetFreq,paymentDates->GetSize()-i); // pour reculer
		}

		pfCap= new ARM_StdPortfolio( assets, weights, prices, vegas);

        // Free Cloned Assets
             
        int sz = assets.size();

        for (i = 0; i < sz; i++)
        {
            delete assets[i];

            assets[i] = NULL;
        }

		return pfCap;
	}
};


//! \brief loader for a Memory SO
class MemorySpreadOptionAssetLoader : public AbstractAssetParser
{
public:
	MemorySpreadOptionAssetLoader()
	{
	}

	~MemorySpreadOptionAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_StdPortfolio* pfSO = NULL;
		vector<ARM_Security*> assets;
		vector<double> weights;
		vector<double> prices;
		vector<double> vegas;

		SwapAssetLoader loader;
		ARM_SwapLeg* newSwapLeg = dynamic_cast<ARM_SwapLeg*>(loader.BuildAsset(node, date,assetId, isEtkOrNot));

		ARM_Date startDate = newSwapLeg->GetEndDateNA();
		ARM_Date endDate = newSwapLeg->GetEndDateNA();
		int resetFreq = newSwapLeg->GetIRIndex()->GetResetFrequency();
		ARM_ReferenceValue* varNotional = newSwapLeg->GetAmount();
		ARM_Vector* paymentDates = newSwapLeg->GetPaymentDates();
		ARM_Vector* endDates = newSwapLeg->GetFlowEndDates();

		SpreadOptionAssetLoader soloader;
		double notional;
		ARM_Date tmpEndDate(endDate);

		for (int i = paymentDates->GetSize()-1; i >= 0; i--)
		{
			notional = varNotional->CptReferenceValue(paymentDates->Elt(i));
			ARM_SpreadOption* newSO = dynamic_cast<ARM_SpreadOption*>(soloader.BuildAssetForMemory(node, date,assetId, endDate, notional, isEtkOrNot));
			assets.push_back(dynamic_cast<ARM_Security*>(newSO));
			weights.push_back(1.0);
			prices.push_back(1.0);
			vegas.push_back(0.0);
			endDate = tmpEndDate;
			endDate.AddPeriodMult(-resetFreq,paymentDates->GetSize()-i); // pour reculer
		}

		pfSO = new ARM_StdPortfolio(assets, weights, prices, vegas);

        // Free Cloned Assets

               
        int sz = assets.size();

        for (i = 0; i < sz; i++)
        {
            delete assets[i];

            assets[i] = NULL;
        }

		return(pfSO);
	}
};


// Corridor summit loader (EXOTIC / SWAP)
// Keywords recognized: RNGLIBORFIX / RNGLIBORLIBOR.
class CorridorAssetLoader : public AbstractAssetParser
{
public:
	CorridorAssetLoader()
	{
	}

	~CorridorAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_CorridorLeg* corridor			= NULL;

		ARM_Date startDate;
		ARM_Date endDate;
		double porS = 1.0;
		int payFreq;
		int resetFreq;
		ARM_Currency* currency				= NULL;
		ARM_ReferenceValue* notional		= NULL;
		ARM_ReferenceValue* fees			= NULL;
		int dayCount;
		int stub;
		int refFwdRule;
		int payFwdRule;
		char payCal[7];
		char resetCal[7];
		char cust[20];
		ARM_IRIndex* paymentIndex			= NULL; 
		ARM_IRIndex* refIndex				= NULL;
		ARM_ReferenceValue* fixRateOrSpread	= NULL;
		ARM_ReferenceValue* barDown			= NULL;
		ARM_ReferenceValue* barUp			= NULL;
		int decompPricingFlag				= 0;
		int resetPayTiming;
		int resetRefTiming;
		int payIndexType					= K_FIXED;
		int refIndexType					= K_LIBOR;
		int refIndexTerm;
		int paidIndexTerm;
		int refIndexResetGap;
		int payIndexResetGap;
		int payIndexPayGap;
		double term = 0.5;
		double alpha = 0.5;
		double beta = 0.5;
		double range = 0.;
		ARM_ReferenceValue* rangeWeight = NULL;

		ARM_ReferenceValue*	pastRefFixings  = NULL;
		ARM_ReferenceValue*	pastPayFixings  = NULL;

		try
		{
			startDate			= GetStartDate(node);
			endDate				= GetEndDate(node);
			char tmpDay[15];
			GetPayRollDate(node,tmpDay);
			int payDay			= atoi(tmpDay);
			porS				= GetPorS(node);
			payFreq				= GetPayFreq(node);
			resetFreq			= GetResetFreq(node);
			currency			= GetCcy(node);
			notional			= GetNotional(node);
			dayCount			= GetAssetDayCount(node);
			stub				= GetStubRule(node);
			refFwdRule			= GetFwdRule(node);
			payFwdRule			= GetPayIndexFwdRule(node);
			decompPricingFlag   = GetDecompFlag(node);
			fees				= GetPrimes(node);

			GetPayCalendar(node, payCal);
			GetResetCalendar(node, resetCal);
			GetCustom(node, cust);
			
			//Reference Index Timing: ADVANCE or ARREARS.
			resetRefTiming		= GetResetTiming(node);
			
			//Payment Index Timing: ADVANCE or ARREARS.
			resetPayTiming		= K_ADVANCE; //GetPayTiming(node);

			refIndexResetGap	= GetResetGap(node);
			payIndexResetGap	= GetPayResetGap(node);
			payIndexPayGap		= GetPayGap(node);

			//Get if Restrikable
			int IsRestrikable = GetIfResrtikableType(node,payIndexType,refIndexType);

			//Started deal.
			if (startDate.GetJulian() < date.GetJulian())
			{
				pastRefFixings	= GetCorridorPastFixings(node,0, date);
				pastPayFixings	= GetCorridorPastFixings(node,1, date); 
			}
			if(IsRestrikable)
			{
				refIndexTerm = GetRestrikableRefIndex(node,currency);
				resetPayTiming = GetPayTiming(node);
				alpha		= GetRestrikableAlpha(node);
				beta		= GetRestrikableBeta(node);
				range		= 0.01 * GetRestrikableRange(node);
				rangeWeight = new ARM_ReferenceValue(range);
			}
			else
			{
				//NGLIBORLIBOR or RNGLIBORFIX 
				payIndexType		= GetRngProductType(node);
		
				//Low Barrier
				barDown				= GetLowBarrier(node, refIndexResetGap, resetCal); 
				
				//Up Barrier
				barUp				= GetUpBarrier(node, refIndexResetGap, resetCal);
				
				//Reference Index
				refIndexTerm		= GetCorridorIndexType(node, currency,  1);
			}

			//FIXED PAYMENT INDEX.
			if (payIndexType == K_FIXED)
			{
				fixRateOrSpread = GetBoostedRate(node);
				if ((!payFreq == K_ZEROCOUPON))
				{
					term = 1.0 / ((double)payFreq);
				}

				 paymentIndex   = new ARM_IRIndex(	dayCount,						//DayCount
													payFreq,						//ResetFreq
													payFreq,						//PayFreq
													term,							//Term
													K_COMP_PROP,					//CompoundMeth
													payFwdRule,						//FwdRule
													K_ADJUSTED,						//IntRule
													K_ADVANCE,						//ResetTiming
													payIndexResetGap,				//ResetGap
													K_ARREARS,						//PayTiming
													payIndexPayGap,					//PayGap
													currency,						//Currency
													(ARM_INDEX_TYPE) IDXFIXED,		//IndexTerm
													K_COMP_PROP);					//DecompFreq
			}
			//VARIABLE PAYMENT INDEX.
			else
			{
				//spread considered as an event.
				fixRateOrSpread = GetCorridorStepUpSpread(node);
				//spread considered as a static data.
				if (!fixRateOrSpread)
				{
					double spread	= GetSpread(node);
					fixRateOrSpread = new ARM_ReferenceValue(100 * spread);
				}
				if(IsRestrikable)
					paidIndexTerm   = GetRestrikablePayIndex(node,currency);
				else
					paidIndexTerm	= GetCorridorIndexType(node, currency, 4);

				term			= 1.0/((double)FromLiborTypeToFrequency(paidIndexTerm));
				
				paymentIndex    = new ARM_IRIndex(	dayCount,						//DayCount
													payFreq,						//ResetFreq
													payFreq,						//PayFreq
													term,							//Term
													K_COMP_PROP,					//CompoundMeth
													payFwdRule,						//FwdRule
													K_ADJUSTED,						//IntRule
													K_ADVANCE,						//ResetTiming
													payIndexResetGap,				//ResetGap
													K_ARREARS,						//PayTiming
													payIndexPayGap,					//PayGap
													currency,						//Currency
													(ARM_INDEX_TYPE) paidIndexTerm,	//IndexTerm
													K_COMP_PROP);					//DecompFreq
				if(payIndexType == K_CMS)
				{
					ARM_INDEX_TYPE Libor_Type = currency->GetVanillaIndexType();
					paymentIndex->Set((ARM_INDEX_TYPE) paidIndexTerm, Libor_Type, payFreq, K_COMP_PROP, currency);
					paymentIndex->Set(dayCount, payFreq, payFreq,
									term, K_COMP_PROP, payFwdRule, K_ADJUSTED,
									K_ADVANCE, payIndexResetGap, K_ARREARS, payIndexPayGap,
									currency, (ARM_INDEX_TYPE) paidIndexTerm, K_COMP_PROP);
				}
			}

			//FIXED REFERENCE INDEX.
			if (refIndexType == K_FIXED)
			{
				term = 1.0 / (double)payFreq;
				refIndex    = new ARM_IRIndex(	dayCount,						//DayCount
												resetFreq,						//ResetFreq
												resetFreq,						//PayFreq
												term,							//Term
												K_COMP_PROP,					//CompoundMeth
												refFwdRule,						//FwdRule
												K_ADJUSTED,						//IntRule
												K_ADVANCE,						//ResetTiming
												refIndexResetGap,				//ResetGap
												K_ARREARS,						//PayTiming
												0,								//PayGap
												currency,						//Currency
												(ARM_INDEX_TYPE) IDXFIXED,		//IndexTerm
												K_COMP_PROP);					//DecompFreq
			}
			//VARIABLE REFERENCE INDEX.
			else
			{
				term = 1.0/((double)FromLiborTypeToFrequency(refIndexTerm));
				refIndex    = new ARM_IRIndex(	dayCount,						//DayCount
												resetFreq,						//ResetFreq
												resetFreq,						//PayFreq
												term,							//Term
												K_COMP_PROP,					//CompoundMeth
												refFwdRule,						//FwdRule
												K_ADJUSTED,						//IntRule
												resetRefTiming,					//ResetTiming
												refIndexResetGap,				//ResetGap
												K_ARREARS,						//PayTiming
												0,								//PayGap
												currency,						//Currency
												(ARM_INDEX_TYPE) refIndexTerm,	//IndexTerm
												K_COMP_PROP);					//DecompFreq
				if(refIndexType == K_CMS)
				{
					ARM_INDEX_TYPE Libor_Type = currency->GetVanillaIndexType();
					refIndex->Set((ARM_INDEX_TYPE) refIndexTerm, Libor_Type, resetFreq, K_COMP_PROP, currency);
					refIndex->Set(dayCount, resetFreq, resetFreq,
									term, K_COMP_PROP, refFwdRule, K_ADJUSTED,
									resetRefTiming, refIndexResetGap, K_ARREARS, 0,
									currency, (ARM_INDEX_TYPE) refIndexTerm, K_COMP_PROP);
				}

			}
			if(IsRestrikable)
			{
				refIndex->SetResetTiming(resetRefTiming);
				refIndex->SetPayTiming(resetPayTiming);
				paymentIndex->SetResetTiming(resetRefTiming);
				paymentIndex->SetPayTiming(resetPayTiming);

				corridor = new ARM_CorridorLeg( startDate, 
												endDate,
												porS, 
												paymentIndex,
												fixRateOrSpread,
												refIndex, 
												stub,
												rangeWeight,
												currency,
												alpha, 
												beta,
												decompPricingFlag,
												1,
												resetCal, 
												payCal,
												pastRefFixings,
												pastPayFixings);
			/*	ARM_Date&			startDate, 
					ARM_Date&			endDate,
                    int					rcvOrPay = K_RCV, 
                    ARM_IRIndex*		PaymentIndex = NULL,
                    ARM_ReferenceValue*	Spread = NULL,
                    ARM_IRIndex*		RefIndex = NULL, 
                    int					stubRule = K_SHORTSTART,
                    ARM_ReferenceValue*	RangeWidth = NULL,
                    ARM_Currency*		discountCcy = ARM_DEFAULT_CURRENCY,
					double				AlphaRestrikable = 0.5, 
					double				BetaRestrikable = 0.5,
                    int					decompPricingFlag=1,
					int					AdjStartDateFlag=1,
					char*				resetCalName = NULL, 
					char*				payCalName = NULL,
					ARM_ReferenceValue* refPastFixings = NULL,
					ARM_ReferenceValue* payPastFixings = NULL);*/
			}
			else
			{
				corridor = new ARM_CorridorLeg(	startDate, 
												endDate,
												porS, 
												paymentIndex,
												payFreq,
												fixRateOrSpread,
												refIndex, 
												resetFreq,
												resetPayTiming, 
												resetRefTiming,
												stub,
												barDown, 
												K_STD,
												barUp, 
												K_STD,
												currency,
												K_DEF_FREQ, 
												K_LINEAR,
												K_DIGITALE,
												decompPricingFlag,
												resetCal,
												payCal, 
												refIndexResetGap,
												pastRefFixings,
												pastPayFixings);
			}

			corridor->SetAmount(notional);

			if(fees)
			{
				corridor->SetFee(fees);
				delete fees;
			}

			if (notional)
			{
				delete notional;
				notional = NULL;
			}

			if (paymentIndex)
			{
				delete paymentIndex;
				paymentIndex = NULL;
			}
			if (refIndex)
			{
				delete refIndex;
				refIndex = NULL;
			}
			if (fixRateOrSpread)
			{
				delete fixRateOrSpread;
				fixRateOrSpread = NULL;
			}
			if (barDown)
			{
				delete barDown;
				barDown = NULL;
			}
			if (barUp)
			{
				delete barUp;
				barUp = NULL;
			}
			if (currency)
			{
				delete currency;
				currency = NULL;
			}
			if (pastRefFixings)
			{
				delete pastRefFixings;
				pastRefFixings = NULL;
			}
			if (pastPayFixings)
			{
				delete pastPayFixings;
				pastPayFixings = NULL;
			}
		
			return corridor;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "CorridorAssetLoader::BuildAsset : Failed in constructing CorridorLeg");
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : CorridorAssetLoader::BuildAsset()");
		} 
	}
};


//! \brief loader for a strip of FX options
class FXOptionStripAssetLoader : public AbstractAssetParser
{
public:
	FXOptionStripAssetLoader()
	{
	}

	~FXOptionStripAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		int isCallOrPut;
		ARM_Date startDate;
		ARM_Date endDate;
		int	resetFreq;
		int	dayCount;
		char resetCalendar[4];
		int fwdRule;
		int intRule;
		int stubRule;
		int resetGap;
		int payFreq;
		int payGap;
		char payCalendar[4];
		int resetTiming;
		int payTiming;
		char* paymentCcy;
		double dPorS = 1.0;
		bool isDigital = false;
		int callSpreadFlag = 1;
		double epsilon = 0.0001;
		ARM_Curve* strikesCurve = NULL;
		ARM_ReferenceValue*	notionalRefValue = NULL;
		ARM_Curve* notionalCurve = NULL;
		ARM_ReferenceValue*	fxFixingsRefValue = NULL;
		ARM_Curve* fxFixingsCurve = NULL;
		ARM_Curve* payOffCurve = NULL;
		ARM_Curve* leverageCurve = NULL;
		ARM_Forex underlying;

		char formula[35];

		ARM_Currency*	Ccy1 = NULL;
		ARM_Currency*	Ccy2 = NULL;

		ARM_StripOption* FxOptionStrip = NULL;
		ARM_StripDigitalOption*	FxDigitalOptionStrip = NULL;

		try
		{
			isCallOrPut = GetCapFloor(node);
			startDate = GetStartDate(node);
			endDate = GetEndDate(node);
			resetFreq = GetResetFreq(node);
			dayCount = GetAssetDayCount(node);
			GetResetCalendar(node, resetCalendar);
			fwdRule = GetFwdRule(node);
			intRule = GetIntRule(node);
			stubRule = GetStubRule(node);
			resetGap = GetResetGap(node);
			payFreq = GetPayFreq(node);
			payGap = GetPayGap(node);
			GetPayCalendar(node, payCalendar);
			resetTiming = GetResetTiming(node);
			payTiming = GetPayTiming(node);
			paymentCcy = GetCcy(node)->GetCcyName();
			dPorS = GetPorS(node);

			GetFormula(node, formula);
			string vFormula(formula);
			if (vFormula.find("FX_CALL_STRIP") != string::npos)
			{
				if( vFormula.find("DIGITAL") != string::npos )
					isDigital = true;

				char subFormula[8];
				GetNthArgInFormula(node, 1, subFormula);

				string	vCcy1 = string(subFormula).substr(1, 3);
				string::size_type vPos = string(subFormula).find("/");
				string	vCcy2 = string(subFormula).substr(vPos+1, 3);
				Ccy1 = new ARM_Currency(vCcy1.c_str());
				Ccy2 = new ARM_Currency(vCcy2.c_str());
			}
			else
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "No 'FX_CALL_STRIP' formula found in XML");
			}

			underlying = ARM_Forex(Ccy1, Ccy2);
			delete	Ccy1;
			delete	Ccy2;

			ARM_ReferenceValue*	strikesRefValue = GetCapStrikes(node);
			strikesRefValue->SetCalcMethod(K_STEPUP_LEFT);
			strikesCurve = RefValueToCurve(*strikesRefValue, date.GetJulian());
			delete strikesRefValue;

			notionalRefValue = GetNotional(node); //STEP UP RIGHT !
			// WARNING : adjust payment dates to business days !
			if (notionalRefValue->size())
			{
				ARM_Vector* newDates = new ARM_Vector(notionalRefValue->size());
				for (int i=0; i<notionalRefValue->size(); i++)
				{
					ARM_Date date(notionalRefValue->GetDiscreteDates()->Elt(i));
					date.AdjustToBusDate(payCalendar, fwdRule);
					newDates->Elt(i) = date.GetJulian();
				}
				notionalRefValue->SetDiscreteDates(newDates);
				delete newDates;
			}
			notionalCurve = RefValueToCurve(*notionalRefValue, date.GetJulian());

			ARM_ReferenceValue*	fxFixingsRefValue = GetFxFixings(node);
			if (fxFixingsRefValue)
			{
				fxFixingsRefValue->SetCalcMethod(K_STEPUP_LEFT);
				fxFixingsCurve = RefValueToCurve(*fxFixingsRefValue, date.GetJulian());
				delete fxFixingsRefValue;
			}

			ARM_ReferenceValue* leverageRefValue = GetRefValueFromFormula(node, 2); //STEP UP LEFT
			if (leverageRefValue)
			{
				leverageCurve = RefValueToCurve(*leverageRefValue, date.GetJulian());
				delete leverageRefValue;
			}

			if(isDigital)
			{
				// payoff is stored in spread !
				ARM_ReferenceValue* payOffRefValue = GetSpreadVariable(node);
				if (!payOffRefValue)
				{
					payOffRefValue = new ARM_ReferenceValue(GetSpread(node)*100.0);
				}
				payOffRefValue->SetCalcMethod(K_STEPUP_LEFT);
				payOffCurve = RefValueToCurve(*payOffRefValue, date.GetJulian());
				delete payOffRefValue;

				char epsilonStr[10];
				GetNthArgInFormula(node, 3, epsilonStr);
				if (strcmp(epsilonStr, "") != 0)
					epsilon = atof(epsilonStr);

				FxDigitalOptionStrip = new ARM_StripDigitalOption(date, &underlying, *strikesCurve, isCallOrPut, startDate, endDate, 
																  resetFreq, dayCount, *notionalCurve, resetCalendar, 
																  fwdRule, intRule, stubRule, resetGap, payFreq, payGap, 
																  payCalendar, resetTiming, payTiming, paymentCcy, dPorS,
																  callSpreadFlag, epsilon, fxFixingsCurve, payOffCurve, leverageCurve);
				FxDigitalOptionStrip->SetPorS(dPorS);
				FxDigitalOptionStrip->SetAmount(notionalRefValue);

				delete notionalRefValue;
				delete strikesCurve;
				delete notionalCurve;
				delete fxFixingsCurve;
				delete leverageCurve;
				delete payOffCurve;

				return FxDigitalOptionStrip;
			}
			else
			{
				FxOptionStrip = new ARM_StripOption(date, &underlying, *strikesCurve, isCallOrPut, startDate, endDate, 
													resetFreq, dayCount, *notionalCurve, resetCalendar, 
													fwdRule, intRule, stubRule, resetGap, payFreq, payGap, 
													payCalendar, resetTiming, payTiming, paymentCcy, dPorS, 
													fxFixingsCurve, leverageCurve);
				FxOptionStrip->SetPorS(dPorS);
				FxOptionStrip->SetAmount(notionalRefValue);

				delete notionalRefValue;
				delete strikesCurve;
				delete notionalCurve;
				delete fxFixingsCurve;
				delete leverageCurve;

				return FxOptionStrip;
			}
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();
			throw x;
		} 
		catch(...)
		{		
			if (strikesCurve)
				delete strikesCurve;

			if (notionalRefValue)
				delete notionalRefValue;

			if (notionalCurve)
				delete notionalCurve;

			if (fxFixingsCurve)
				delete fxFixingsCurve;
			
			if (payOffCurve)
				delete payOffCurve;

			if (leverageCurve)
				delete leverageCurve;

			CCString msg((CCString)"Error in XML parsing for getting Fx Option Strip");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
		}
	}
};

//! \brief loader for a strip of FX Quanto options
class FXStripQuantoAssetLoader : public AbstractAssetParser
{
public:
	FXStripQuantoAssetLoader()
	{
	}

	~FXStripQuantoAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_Date		startDate, endDate;
		int				expiryGap, setlmentGap, payGap, resetGap, payFreq, resetFreq;
		int				dayCount, callPut, callPut2;
		int				fwdRule, intRule, stubRule, resetTiming, payTiming;
		int				leverageIndex, epsilonIndex;
		string			resetCal, payCal, fx1Name, fx2Name;
		string			currency1, currency2, currency3, currency4, payCurrency;
		double			epsilon(0.001), leverageDbl(1.0);
		double			alphaDbl(1.0), betaDbl(-1.0);
		double			spread(0.0);
		bool			isDigital(false);
		bool			isQuanto(false), isSpread(false);
		int				i(0), adjFirstDate(0), perf(1);
		int				pOrS(0);
		
		ARM_ReferenceValue* notionalRef	= NULL;
		ARM_ReferenceValue* leverageRef	= NULL;
		ARM_ReferenceValue* strikeRef	= NULL;
		ARM_ReferenceValue* strike2Ref  = NULL;
		ARM_ReferenceValue* alphaRef	= NULL;
		ARM_ReferenceValue* betaRef		= NULL;
		ARM_ReferenceValue* spreadRef	= NULL;
		ARM_ReferenceValue*	fxFixingsRefValue = NULL;
		ARM_ReferenceValue*	fxFixingsRefValue2 = NULL;

		ARM_Curve*			notionalCrv	= NULL;
		ARM_Curve*			leverageCrv	= NULL;
		ARM_Curve*			strikeCrv	= NULL;
		ARM_Curve*			strike2Crv  = NULL;
		ARM_Curve*			alphaCrv	= NULL;
		ARM_Curve*			betaCrv		= NULL;
		ARM_Curve*			spreadCrv	= NULL;
		ARM_Curve*			fxFixingsCurve = NULL;
		ARM_Curve*			fxFixingsCurve2 = NULL;

		ARM_Currency*		payCcy		= NULL;

		ARM_Forex*			Fx1			= NULL;
		ARM_Forex*			Fx2			= NULL;

		char resetCalendar[4];
		char payCalendar[4];
		char formula[35];

		ARM_FxSpreadStripOption* FxSpreadStripOption = NULL;

		try
		{
			// Parsing of all standard parameters
			//---------------------------------------
			startDate	= GetStartDate(node);
			endDate		= GetEndDate(node);
			resetGap	= GetResetGap(node);
			payGap		= GetPayGap(node);
			resetFreq	= GetResetFreq(node);
			payFreq		= GetPayFreq(node);
			dayCount	= GetAssetDayCount(node);
			callPut		= GetCapFloor(node);

			// Default
			callPut2	= callPut;
			expiryGap	= resetGap;
			setlmentGap	= ARM_DFLTSETTLEMENTGAP;
			
			GetResetCalendar(node, resetCalendar);
			GetPayCalendar  (node, payCalendar  );

			resetCal	= resetCalendar;
			payCal		= payCalendar;

			fwdRule		= GetFwdRule(node);
			intRule		= GetIntRule(node);
			stubRule	= GetStubRule(node);
			resetTiming	= GetResetTiming(node);
			payTiming	= GetPayTiming(node);
			payCcy		= GetCcy(node);
			pOrS		= GetPorS(node);

			
			// dateStrip construction
			//--------------------------------------
			ARM_DateStrip dateStrip(startDate,
									endDate, 
									resetFreq,
									dayCount,
									resetCalendar,
									fwdRule,
									intRule,
									stubRule,
									resetGap,
									payFreq,
									payGap,
									payCalendar,
									resetTiming,
									payTiming,
									adjFirstDate);

			// Get the Summit formula
			//--------------------------------------
			GetFormula(node, formula);
			string vFormula(formula);

			
			// Check the formula
			//--------------------------------------
			if (vFormula.find("FX_QUANTO_STRIP") != string::npos)
			{
				isQuanto	  = true;
				leverageIndex = 2;
				epsilonIndex  = 4;
			}
			else if (vFormula.find("FX_SPREAD_STRIP") != string::npos)
			{
				isSpread	  = true;
				leverageIndex = 5;
				epsilonIndex  = 6;	
			}
			else
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 "No 'FX_QUANTO_STRIP' or 'FX_SPREAD_STRIP' formula found in XML");
			}

			
			// Get the FX names (name1 and name2)
			//--------------------------------------
			if( vFormula.find("DIGITAL") != string::npos )
			{
				isDigital = true;
			}

			char subFormula[9];

			GetNthArgInFormula(node, 1, subFormula);
			
			currency1	= string(subFormula).substr(1, 3);
			currency2	= string(subFormula).substr(5, 3);
			payCurrency	= string(payCcy->GetCcyName());

			Fx1 = new ARM_Forex(&ARM_Currency(currency1.c_str()), &ARM_Currency(currency2.c_str()));
			
			if ( isQuanto )
			{
				if ( (currency1 == payCurrency) || (currency2 == payCurrency) )
				{
					Fx2 = static_cast<ARM_Forex*>(Fx1->Clone());
				}
				else
				{
					Fx2 = new ARM_Forex(&ARM_Currency(currency1.c_str()), &ARM_Currency(payCurrency.c_str()));
				}
			}
			else
			{
				GetNthArgInFormula(node, 2, subFormula);

				currency3 = string(subFormula).substr(1, 3);
				currency4 = string(subFormula).substr(5, 3);

				Fx2 = new ARM_Forex(&ARM_Currency(currency3.c_str()), &ARM_Currency(currency4.c_str()));
			}

			// Get the strikes (strike1 & strike2) 
			//--------------------------------------
			strikeRef = GetCapStrikes(node);
			strikeRef->SetCalcMethod(K_STEPUP_LEFT);
			ARM_ReferenceValue* newStrikeRef = ARM_FromStartDSToResetDS(strikeRef,dateStrip);
			strikeCrv = RefValueToCurve(*newStrikeRef, date.GetJulian());

			delete newStrikeRef;
			newStrikeRef = NULL;
			delete strikeRef;
			strikeRef = NULL;

			// Default
			strike2Ref  = new ARM_ReferenceValue(0.0);
			strike2Crv  = RefValueToCurve(*strike2Ref, date.GetJulian());
			
			delete strike2Ref;
			strike2Ref  = NULL;

			// Get The Spreads
			//----------------
			spread		= GetSpread(node);
			spreadRef	= GetSpreadVariable(node);

			if ( spreadRef )
			{
				ARM_ReferenceValue* newSpreadRef = ARM_FromStartDSToResetDS(spreadRef,dateStrip);
				spreadCrv = RefValueToCurve(*newSpreadRef, date.GetJulian());

				delete spreadRef;
				spreadRef = NULL;
				delete newSpreadRef;
				newSpreadRef = NULL;
			}
			else
			{
				if ( !(isDigital) )
				{
					spread = 1;
				}

				ARM_ReferenceValue* newSpreadRef = new ARM_ReferenceValue((100.0 * spread));
				spreadCrv = RefValueToCurve(*newSpreadRef, date.GetJulian());

				delete newSpreadRef;
				newSpreadRef = NULL;
			}

			// Get the notional 
			//--------------------------------------
			notionalRef = GetNotional(node, false);
			ARM_ReferenceValue* newNotional = ARM_FromStartDSToResetDS(notionalRef,dateStrip);
			notionalCrv = RefValueToCurve(*newNotional, date.GetJulian());
			

			if ( spreadCrv )
			{
				// Special treatment when the spread is not constant
				// --> notionalCrv receives all spreadCrv abscisses that it originally misses

				int notionalCrvSize = notionalCrv->GetAbscisses().size();
				int spreadCrvSize	= spreadCrv->GetAbscisses().size();

				// temporary objects
				ARM_Vector newAbscisses(notionalCrvSize + spreadCrvSize);
				
				ARM_Vector* tmpVec;
				ARM_Curve*  tmpCurve = new ARM_Curve(*notionalCrv);

				for (i=0; i<notionalCrvSize; i++)
				{
					newAbscisses[i] = notionalCrv->GetAbscisse(i);
				}

				for (i=0; i<spreadCrvSize; i++)
				{
					newAbscisses[notionalCrvSize+i] = spreadCrv->GetAbscisse(i);
				}

				// Sort tmpVec and delete multiple values
				tmpVec = newAbscisses.Sort_Compact();

				delete notionalCrv;
				notionalCrv = new ARM_Curve;

				notionalCrv->SetInterpolator(tmpCurve->GetInterpolator()->Clone());

				for (i=0; i<tmpVec->size(); i++)
				{
					notionalCrv->insert((*tmpVec)[i], tmpCurve->Interpolate((*tmpVec)[i]));
				}

				// notionalCrv computation
				for (i=0; i<notionalCrv->size(); i++)
				{
					(*notionalCrv)[i] *= spreadCrv->Interpolate(notionalCrv->GetAbscisse(i)) / 100.;
				}

				delete tmpCurve;
				delete tmpVec;
			}

			delete spreadCrv;
			spreadCrv = NULL;

			// Get the leverage 
			//--------------------------------------
			leverageRef = GetRefValueFromFormula(node, leverageIndex);

			// Default case when leverageRef is NULL
			if ( !leverageRef )
			{
				leverageRef = new ARM_ReferenceValue(leverageDbl);
			}

			ARM_ReferenceValue* newLeverage = ARM_FromStartDSToResetDS(leverageRef,dateStrip);
			leverageCrv = RefValueToCurve(*newLeverage, date.GetJulian());

			delete newLeverage;
			newLeverage = NULL;
			delete leverageRef;
			leverageRef = NULL;


			// Get alpha and beta coefficients
			//--------------------------------------
			if ( isSpread)
			{
				// Alpha
				//--------------------------------------
				alphaRef = GetRefValueFromFormula(node, 3);


				// Beta
				//--------------------------------------
				betaRef = GetRefValueFromFormula(node, 4);
			}

			// Default case when alphaRef is NULL
			if ( !alphaRef )
			{
				alphaRef = new ARM_ReferenceValue(alphaDbl);
			}

			// Default case when betaRef is NULL
			if ( !betaRef )
			{
				betaRef = new ARM_ReferenceValue(betaDbl);
			}

			ARM_ReferenceValue* newAlpha = ARM_FromStartDSToResetDS(alphaRef,dateStrip);
			alphaCrv = RefValueToCurve(*newAlpha, date.GetJulian());

			delete newAlpha;
			newAlpha = NULL;
			delete alphaRef;
			alphaRef = NULL;

			ARM_ReferenceValue* newBeta = ARM_FromStartDSToResetDS(betaRef,dateStrip);
			betaCrv  = RefValueToCurve(*newBeta, date.GetJulian());
			(*betaCrv) *= -1.;

			delete newBeta;
			newBeta = NULL;
			delete betaRef;
			betaRef =  NULL;

			ARM_ReferenceValue*	fxFixingsRefValue = GetFxFixings(node);
			if (fxFixingsRefValue)
			{
				fxFixingsRefValue->SetCalcMethod(K_STEPUP_LEFT);
				fxFixingsCurve = RefValueToCurve(*fxFixingsRefValue, date.GetJulian());
				delete fxFixingsRefValue;
			}

			fxFixingsRefValue2 = GetRefValSpreadFirstFixing(node,date);
			if (fxFixingsRefValue2)
			{
				fxFixingsRefValue2->SetCalcMethod(K_STEPUP_LEFT);
				fxFixingsCurve2 = RefValueToCurve(*fxFixingsRefValue2, date.GetJulian());
				delete fxFixingsRefValue2;
			}

			// Get performance (FX_QUANTO_STRIP case)
			//--------------------------------------
			if ( isQuanto )
			{
				GetNthArgInFormula(node, 3, subFormula);
				
				if ( subFormula[0] != 0 )
				{
					perf = atoi(subFormula);
				}
			}

			// Get epsilon (DIGITAL case)
			//--------------------------------------
			if ( isDigital )
			{
				GetNthArgInFormula(node, epsilonIndex, subFormula);

				if ( subFormula[0] != 0 )
				{
					epsilon = atof(subFormula);
				}
			}

			// FxSpreadStripOption construction
			//--------------------------------------
			FxSpreadStripOption = new ARM_FxSpreadStripOption(	date,
																dateStrip,
																(*strikeCrv),
																Fx1,
																Fx2,
																callPut,
																NULL,
																0.5,
																setlmentGap,
																payCurrency.c_str(),
																pOrS,
																notionalCrv,
																leverageCrv,
																alphaCrv,
																betaCrv,
																false,
																isQuanto,
																isDigital,
																epsilon,
																perf);

			FxSpreadStripOption->SetPorS(pOrS);
			FxSpreadStripOption->SetAmount(newNotional);

			// Once amount is set, newNotional can be deleted
			delete newNotional;
			newNotional = NULL;
			delete notionalRef;
			notionalRef = NULL;

			if (fxFixingsCurve || fxFixingsCurve2)
			{
				ARM_FixingSched* fixingsSched = new ARM_FixingSched(date);

				if (fxFixingsCurve)
				{
					if (isSpread)
						fixingsSched->AddFxFixing(currency3 + "_" + currency4,fxFixingsCurve);
					else
						fixingsSched->AddFxFixing(currency1 + "_" + currency2,fxFixingsCurve);
				}

				if (fxFixingsCurve2 && isSpread)
				{
					fixingsSched->AddFxFixing(currency1 + "_" + currency2,fxFixingsCurve2);
				}

				FxSpreadStripOption->SetFixings(fixingsSched);

				if (fxFixingsCurve)
				{
					delete fxFixingsCurve;
					fxFixingsCurve = NULL;
				}

				if (fxFixingsCurve2)
				{
					delete fxFixingsCurve2;
					fxFixingsCurve2 = NULL;
				}

				if (fixingsSched)
					delete fixingsSched;
				fixingsSched = NULL;
			}

			// Delete Currencies
			//--------------------------------------
			if ( Fx1 )
			{
				delete Fx1;
				Fx1 = NULL;
			}

			if ( Fx2 )
			{
				delete Fx2;
				Fx2 = NULL;
			}

			// End of function
			//--------------------------------------
			return FxSpreadStripOption;
		}

		// Catch managed exceptions
		//--------------------------------------
		catch (Exception& x) 
		{
	        x.DebugPrint();
			throw x;
		}

		// Catch non managed exceptions
		// --> Delete everything
		//--------------------------------------
		catch(...)
		{
			if (payCcy)
			{
				delete payCcy;
				payCcy = NULL;
			}

			if (notionalRef)
			{
				delete notionalRef;
				notionalRef = NULL;
			}

			if (strikeRef)
			{
				delete strikeRef;
				strikeRef = NULL;
			}

			if (strike2Ref)
			{
				delete strike2Ref;
				strike2Ref = NULL;
			}

			if (leverageRef)
			{
				delete leverageRef;
				leverageRef = NULL;
			}

			if (alphaRef)
			{
				delete alphaRef;
				alphaRef = NULL;
			}

			if (betaRef)
			{
				delete betaRef;
				betaRef = NULL;
			}

			if (notionalCrv)
			{
				delete notionalCrv;
				notionalCrv = NULL;
			}

			if (strikeCrv)
			{
				delete strikeCrv;
				strikeCrv = NULL;
			}

			if (strike2Crv)
			{
				delete strike2Crv;
				strike2Crv = NULL;
			}

			if (leverageCrv)
			{
				delete leverageCrv;
				leverageCrv = NULL;
			}

			if (alphaCrv)
			{
				delete alphaCrv;
				alphaCrv = NULL;
			}

			if (betaCrv)
			{
				delete betaCrv;
				betaCrv = NULL;
			}

			if ( Fx1 )
			{
				delete Fx1;
				Fx1 = NULL;
			}

			if ( Fx2 )
			{
				delete Fx2;
				Fx2 = NULL;
			}

			if (fxFixingsCurve)
			{
				delete fxFixingsCurve;
				fxFixingsCurve = NULL;
			}

			if (fxFixingsCurve2)
			{
				delete fxFixingsCurve2;
				fxFixingsCurve2 = NULL;
			}

			CCString msg((CCString) "Error in XML parsing for getting Fx Quanto Option Strip");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
		}
	}
};

//! \brief loader for a strip of FX Quanto options
class FXStripQuantoCalculatorLoader : public AbstractAssetParser
{
public:
	FXStripQuantoCalculatorLoader()
	{
	}

	~FXStripQuantoCalculatorLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_Date		startDate, endDate;
		int				expiryGap, setlmentGap, payGap, resetGap, payFreq, resetFreq;
		int				dayCount, callPut, callPut2;
		int				fwdRule, intRule, stubRule, resetTiming, payTiming;
		int				leverageIndex, epsilonIndex;
		string			resetCal, payCal, fx1Name, fx2Name;
		string			currency1, currency2, currency3, currency4, payCurrency;
		double			epsilon(0.001), leverageDbl(1.0);
		double			alphaDbl(1.0), betaDbl(-1.0);
		double			spread(0.0);
		bool			isDigital(false);
		bool			isQuanto(false), isSpread(false);
		int				i(0), adjFirstDate(0), perf(1);
		int				pOrS(0);
		
		ARM_ReferenceValue* notionalRef	= NULL;
		ARM_ReferenceValue* leverageRef	= NULL;
		ARM_ReferenceValue* strikeRef	= NULL;
		ARM_ReferenceValue* strike2Ref  = NULL;
		ARM_ReferenceValue* alphaRef	= NULL;
		ARM_ReferenceValue* betaRef		= NULL;
		ARM_ReferenceValue* spreadRef	= NULL;
		ARM_ReferenceValue*	fxFixingsRefValue = NULL;
		ARM_ReferenceValue*	fxFixingsRefValue2 = NULL;

		ARM_Curve*			notionalCrv	= NULL;
		ARM_Curve*			leverageCrv	= NULL;
		ARM_Curve*			strikeCrv	= NULL;
		ARM_Curve*			strike2Crv  = NULL;
		ARM_Curve*			alphaCrv	= NULL;
		ARM_Curve*			betaCrv		= NULL;
		ARM_Curve*			spreadCrv	= NULL;
		ARM_Curve*			fxFixingsCurve = NULL;
		ARM_Curve*			fxFixingsCurve2 = NULL;

		// Calculator parameters
		ARM::ARM_BasketType		minMax;
		ARM::ARM_DigitType		digitType;
		ARM::ARM_VanillaType	vanillaType;

		ARM_Currency*		payCcy		= NULL;

		char resetCalendar[4];
		char payCalendar[4];
		char formula[35];

		ARM_FXVanillaCalculator* calculator = NULL;

		try
		{
			// Parsing of all standard parameters
			//---------------------------------------
			startDate	= GetStartDate(node);
			endDate		= GetEndDate(node);
			resetGap	= GetResetGap(node);
			payGap		= GetPayGap(node);
			resetFreq	= GetResetFreq(node);
			payFreq		= GetPayFreq(node);
			dayCount	= GetAssetDayCount(node);
			callPut		= GetCapFloor(node);

			// Default
			callPut2	= callPut;
			expiryGap	= resetGap;
			setlmentGap	= ARM_DFLTSETTLEMENTGAP;
			
			GetResetCalendar(node, resetCalendar);
			GetPayCalendar  (node, payCalendar  );

			resetCal	= resetCalendar;
			payCal		= payCalendar;

			fwdRule		= GetFwdRule(node);
			intRule		= GetIntRule(node);
			stubRule	= GetStubRule(node);
			resetTiming	= GetResetTiming(node);
			payTiming	= GetPayTiming(node);
			payCcy		= GetCcy(node);
			pOrS		= GetPorS(node);

			
			// dateStrip construction
			//--------------------------------------
			ARM_DateStrip dateStrip(startDate,
									endDate, 
									resetFreq,
									dayCount,
									resetCalendar,
									fwdRule,
									intRule,
									stubRule,
									resetGap,
									payFreq,
									payGap,
									payCalendar,
									resetTiming,
									payTiming,
									adjFirstDate);

			// Get the Summit formula
			//--------------------------------------
			GetFormula(node, formula);
			string vFormula(formula);

			
			// Check the formula
			//--------------------------------------
			if (vFormula.find("FX_QUANTO_STRIP") != string::npos)
			{
				isQuanto	  = true;
				leverageIndex = 2;
				epsilonIndex  = 4;
			}
			else if (vFormula.find("FX_SPREAD_STRIP") != string::npos)
			{
				isSpread	  = true;
				leverageIndex = 5;
				epsilonIndex  = 6;	
			}
			else
			{
				throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
								 "No 'FX_QUANTO_STRIP' or 'FX_SPREAD_STRIP' formula found in XML");
			}

			
			// Get the FX names (name1 and name2)
			//--------------------------------------
			if( vFormula.find("DIGITAL") != string::npos )
			{
				isDigital = true;
			}

			char subFormula[9];

			GetNthArgInFormula(node, 1, subFormula);
			
			currency1	= string(subFormula).substr(1, 3);
			currency2	= string(subFormula).substr(5, 3);
			payCurrency	= string(payCcy->GetCcyName());

			fx1Name = currency1 + currency2;
			
			if ( isQuanto )
			{
				if ( (currency1 == payCurrency) || (currency2 == payCurrency) )
				{
					fx2Name = fx1Name;
				}
				else
				{
					fx2Name = currency1 + payCurrency;
				}
			}
			else
			{
				GetNthArgInFormula(node, 2, subFormula);

				currency3 = string(subFormula).substr(1, 3);
				currency4 = string(subFormula).substr(5, 3);

				fx2Name = currency3 + currency4;
			}

			// Get the strikes (strike1 & strike2) 
			//--------------------------------------
			strikeRef = GetCapStrikes(node);
			strikeRef->SetCalcMethod(K_STEPUP_LEFT);
			ARM_ReferenceValue* newStrikeRef = ARM_FromStartDSToResetDS(strikeRef,dateStrip);
			strikeCrv = RefValueToCurve(*newStrikeRef, date.GetJulian());

			delete newStrikeRef;
			newStrikeRef = NULL;
			delete strikeRef;
			strikeRef = NULL;

			// Default
			strike2Ref  = new ARM_ReferenceValue(0.0);
			strike2Crv  = RefValueToCurve(*strike2Ref, date.GetJulian());
			
			delete strike2Ref;
			strike2Ref  = NULL;

			// Get The Spreads
			//----------------
			spread		= GetSpread(node);
			spreadRef	= GetSpreadVariable(node);

			if ( spreadRef )
			{
				ARM_ReferenceValue* newSpreadRef = ARM_FromStartDSToResetDS(spreadRef,dateStrip);
				spreadCrv = RefValueToCurve(*newSpreadRef, date.GetJulian());

				delete spreadRef;
				spreadRef = NULL;
				delete newSpreadRef;
				newSpreadRef = NULL;
			}
			else
			{
				if ( !(isDigital) )
				{
					spread = 1;
				}

				ARM_ReferenceValue* newSpreadRef = new ARM_ReferenceValue((100.0 * spread));
				spreadCrv = RefValueToCurve(*newSpreadRef, date.GetJulian());

				delete newSpreadRef;
				newSpreadRef = NULL;
			}

			// Get the notional 
			//--------------------------------------
			notionalRef = GetNotional(node, false);
			ARM_ReferenceValue* newNotional = ARM_FromStartDSToResetDS(notionalRef,dateStrip);
			notionalCrv = RefValueToCurve(*newNotional, date.GetJulian());
			

			if ( spreadCrv )
			{
				// Special treatment when the spread is not constant
				// --> notionalCrv receives all spreadCrv abscisses that it originally misses

				int notionalCrvSize = notionalCrv->GetAbscisses().size();
				int spreadCrvSize	= spreadCrv->GetAbscisses().size();

				// temporary objects
				ARM_Vector newAbscisses(notionalCrvSize + spreadCrvSize);
				
				ARM_Vector* tmpVec;
				ARM_Curve*  tmpCurve = new ARM_Curve(*notionalCrv);

				for (i=0; i<notionalCrvSize; i++)
				{
					newAbscisses[i] = notionalCrv->GetAbscisse(i);
				}

				for (i=0; i<spreadCrvSize; i++)
				{
					newAbscisses[notionalCrvSize+i] = spreadCrv->GetAbscisse(i);
				}

				// Sort tmpVec and delete multiple values
				tmpVec = newAbscisses.Sort_Compact();

				delete notionalCrv;
				notionalCrv = new ARM_Curve;

				notionalCrv->SetInterpolator(tmpCurve->GetInterpolator()->Clone());

				for (i=0; i<tmpVec->size(); i++)
				{
					notionalCrv->insert((*tmpVec)[i], tmpCurve->Interpolate((*tmpVec)[i]));
				}

				// notionalCrv computation
				for (i=0; i<notionalCrv->size(); i++)
				{
					(*notionalCrv)[i] *= spreadCrv->Interpolate(notionalCrv->GetAbscisse(i)) / 100.;
				}

				delete tmpCurve;
				delete tmpVec;
			}

			delete spreadCrv;
			spreadCrv = NULL;


			// Get the leverage 
			//--------------------------------------
			leverageRef = GetRefValueFromFormula(node, leverageIndex);

			// Default case when leverageRef is NULL
			if ( !leverageRef )
			{
				leverageRef = new ARM_ReferenceValue(leverageDbl);
			}

			ARM_ReferenceValue* newLeverage = ARM_FromStartDSToResetDS(leverageRef,dateStrip);
			leverageCrv = RefValueToCurve(*newLeverage, date.GetJulian());

			delete newLeverage;
			newLeverage = NULL;
			delete leverageRef;
			leverageRef = NULL;


			// Get alpha and beta coefficients
			//--------------------------------------
			if ( isSpread)
			{
				// Alpha
				//--------------------------------------
				alphaRef = GetRefValueFromFormula(node, 3);


				// Beta
				//--------------------------------------
				betaRef = GetRefValueFromFormula(node, 4);
			}

			// Default case when alphaRef is NULL
			if ( !alphaRef )
			{
				alphaRef = new ARM_ReferenceValue(alphaDbl);
			}

			// Default case when betaRef is NULL
			if ( !betaRef )
			{
				betaRef = new ARM_ReferenceValue(betaDbl);
			}

			ARM_ReferenceValue* newAlpha = ARM_FromStartDSToResetDS(alphaRef,dateStrip);
			alphaCrv = RefValueToCurve(*newAlpha, date.GetJulian());

			delete newAlpha;
			newAlpha = NULL;
			delete alphaRef;
			alphaRef = NULL;

			ARM_ReferenceValue* newBeta = ARM_FromStartDSToResetDS(betaRef,dateStrip);
			betaCrv  = RefValueToCurve(*newBeta, date.GetJulian());
			(*betaCrv) *= -1.;

			delete newBeta;
			newBeta = NULL;
			delete betaRef;
			betaRef =  NULL;

			ARM_ReferenceValue*	fxFixingsRefValue = GetFxFixings(node);
			if (fxFixingsRefValue)
			{
				fxFixingsRefValue->SetCalcMethod(K_STEPUP_LEFT);
				fxFixingsCurve = RefValueToCurve(*fxFixingsRefValue, date.GetJulian());
				delete fxFixingsRefValue;
			}

			fxFixingsRefValue2 = GetRefValSpreadFirstFixing(node,date);
			if (fxFixingsRefValue2)
			{
				fxFixingsRefValue2->SetCalcMethod(K_STEPUP_LEFT);
				fxFixingsCurve2 = RefValueToCurve(*fxFixingsRefValue2, date.GetJulian());
				delete fxFixingsRefValue2;
			}

			// Get performance (FX_QUANTO_STRIP case)
			//--------------------------------------
			if ( isQuanto )
			{
				GetNthArgInFormula(node, 3, subFormula);
				
				if ( subFormula[0] != 0 )
				{
					perf = atoi(subFormula);
				}
			}

			// Get epsilon (DIGITAL case)
			//--------------------------------------
			if ( isDigital )
			{
				GetNthArgInFormula(node, epsilonIndex, subFormula);

				if ( subFormula[0] != 0 )
				{
					epsilon = atof(subFormula);
				}
			}

			// Basket Type
			//--------------------------------------
			minMax		= ARM::ARM_FXBasketType::max;


			// Digit Type
			//--------------------------------------
			digitType	= ARM::ARM_FXDigitType::centred; 


			// Vanilla Type
			//--------------------------------------
			if ( isQuanto )
			{
				if ( isDigital )
				{
					vanillaType = ARM::ARM_FXVanillaType::digit;
				}
				else
				{
					if ( perf )
					{
						vanillaType = ARM::ARM_FXVanillaType::perf;
					}
					else
					{
						vanillaType	= ARM::ARM_FXVanillaType::vanilla;
					}
				}
			}
			else
			{
				if ( isDigital )
				{
					vanillaType = ARM::ARM_FXVanillaType::digitspread;
				}
				else
				{
					vanillaType	= ARM::ARM_FXVanillaType::spread;
				}
			}

			// FXVanillaCalculator construction
			//--------------------------------------
			calculator = new ARM_FXVanillaCalculator(date,
													 dateStrip,
													 startDate,
													 endDate,
													 (*payCcy),
													 expiryGap,
													 setlmentGap,
													 resetFreq,
													 dayCount,
													 resetCal,
													 payCal,
													 fx1Name,
													 fx2Name,
													 (*leverageCrv),
													 (*notionalCrv),
													 (*strikeCrv),
													 callPut,
													 vanillaType,
													 (*alphaCrv),
													 (*betaCrv),
													 digitType,
													 epsilon,
													 (*strikeCrv),						// check: strikeCrv2
													 callPut,							// check: callPut2
													 minMax);





			//FxSpreadStripOption->SetPorS(pOrS);
			//FxSpreadStripOption->SetAmount(newNotional);

			// Once amount is set, newNotional can be deleted
			delete newNotional;
			newNotional = NULL;
			delete notionalRef;
			notionalRef = NULL;

			if (fxFixingsCurve || fxFixingsCurve2)
			{
				ARM_FixingSched* fixingsSched = new ARM_FixingSched(date);

				if (fxFixingsCurve)
				{
					if (isSpread)
						fixingsSched->AddFxFixing(currency3 + "_" + currency4,fxFixingsCurve);
					else
						fixingsSched->AddFxFixing(currency1 + "_" + currency2,fxFixingsCurve);
				}

				if (fxFixingsCurve2 && isSpread)
				{
					fixingsSched->AddFxFixing(currency1 + "_" + currency2,fxFixingsCurve2);
				}

				calculator->SetFixings(fixingsSched);

				if (fxFixingsCurve)
				{
					delete fxFixingsCurve;
					fxFixingsCurve = NULL;
				}

				if (fxFixingsCurve2)
				{
					delete fxFixingsCurve2;
					fxFixingsCurve2 = NULL;
				}

				if (fixingsSched)
					delete fixingsSched;
				fixingsSched = NULL;
			}

			// End of function
			//--------------------------------------
			return calculator;
		}

		// Catch managed exceptions
		//--------------------------------------
		catch (Exception& x) 
		{
	        x.DebugPrint();
			throw x;
		}

		// Catch non managed exceptions
		// --> Delete everything
		//--------------------------------------
		catch(...)
		{
			if (payCcy)
			{
				delete payCcy;
				payCcy = NULL;
			}

			if (notionalRef)
			{
				delete notionalRef;
				notionalRef = NULL;
			}

			if (strikeRef)
			{
				delete strikeRef;
				strikeRef = NULL;
			}

			if (strike2Ref)
			{
				delete strike2Ref;
				strike2Ref = NULL;
			}

			if (leverageRef)
			{
				delete leverageRef;
				leverageRef = NULL;
			}

			if (alphaRef)
			{
				delete alphaRef;
				alphaRef = NULL;
			}

			if (betaRef)
			{
				delete betaRef;
				betaRef = NULL;
			}

			if (notionalCrv)
			{
				delete notionalCrv;
				notionalCrv = NULL;
			}

			if (strikeCrv)
			{
				delete strikeCrv;
				strikeCrv = NULL;
			}

			if (strike2Crv)
			{
				delete strike2Crv;
				strike2Crv = NULL;
			}

			if (leverageCrv)
			{
				delete leverageCrv;
				leverageCrv = NULL;
			}

			if (alphaCrv)
			{
				delete alphaCrv;
				alphaCrv = NULL;
			}

			if (betaCrv)
			{
				delete betaCrv;
				betaCrv = NULL;
			}

			if (fxFixingsCurve)
			{
				delete fxFixingsCurve;
				fxFixingsCurve = NULL;
			}

			if (fxFixingsCurve2)
			{
				delete fxFixingsCurve2;
				fxFixingsCurve2 = NULL;
			}

			CCString msg((CCString) "Error in XML parsing for getting Fx Quanto Option Strip");
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT, (const char *) msg);
		}
	}
};

// Corridor Double Condition summit loader (SWAP)
// Keywords recognized: RNG_DOUBLE
class CorridorDoubleConditionAssetLoader : public AbstractAssetParser
{
public:
	CorridorDoubleConditionAssetLoader()
	{
	}

	~CorridorDoubleConditionAssetLoader()
	{
	}

	ARM_Object* BuildAsset(MSXML2::IXMLDOMNode *node,const ARM_Date& date,string& assetId,int isEtkOrNot = 1)
	{
		ARM_CorridorDblCondition* corridor			= NULL;

		ARM_Date startDate;
		ARM_Date endDate;
		int payFreq;
		ARM_Currency* currency				= NULL;
		ARM_Currency* Idx1CcyTmp			= NULL;
		ARM_Currency* Idx2CcyTmp				= NULL;
		ARM_Currency* Idx3CcyTmp				= NULL;
		ARM_Currency* Idx4CcyTmp				= NULL;
		ARM_Currency* Idx5CcyTmp				= NULL;

		ARM_ReferenceValue* refFixingCond1STIndex = NULL;
		ARM_ReferenceValue* refFixingCond1LTIndex = NULL;
		ARM_ReferenceValue* refFixingCond2STIndex = NULL;
		ARM_ReferenceValue* refFixingCond2LTIndex = NULL;
		ARM_ReferenceValue* refFixingPay = NULL;

		ARM_ReferenceValue* digitalBarrier	= NULL;
		ARM_ReferenceValue* spreadBarrier	= NULL;
		ARM_ReferenceValue* payIndexMargin	= NULL;
		
		ARM_ReferenceValue* notional		= NULL;

		int dayCount;
		int stubRule						= K_SHORTSTART;
		int intRule							= K_ADJUSTED;

		ARM_IRIndex* paymentIndex			= NULL; 
		
		int payIndexType					= K_FIXED;
		
		double term = 0.5;
		int fixCondTiming					=K_ADVANCE;
		int payTiming						=K_ARREARS;
		int fixCondResetGap;
		int capOrFloorCondition2;
		int capOrFloorCondition1;
		int resetFreq;

		double condition1STIndexWeight=0;
		double condition1LTIndexWeight=0;
		double condition2STIndexWeight=0;
		double condition2LTIndexWeight=0;
		double payIndexCoeff;

		int freezeFixing					= 0;

		int Condition1ShortTermIndexType = K_LIBOR3M;
		int Condition1LongTermIndexType = K_LIBOR3M;
		int Condition2ShortTermIndexType = K_LIBOR3M;
		int Condition2LongTermIndexType = K_LIBOR3M;
		
		ARM_ReferenceValue* refPayFixedRate   = NULL;

		ARM_VolLInterpol* shiftCorrelVol = NULL;

		ARM_Date tmpDate;
		tmpDate=date;
		int usedCorrel;
		double correlValue;

		try
		{
			startDate			= GetStartDate(node);
			endDate				= GetEndDate(node);
			payFreq				= GetBLOBNum(node,4);

			// convertir en ARM
            char payfrqStr[5];
			sprintf(payfrqStr, "%i", payFreq);
			payFreq = FromSummitFreqToARMFreq(payfrqStr);

			fixCondTiming		= GetBLOBTiming(node,1);
			fixCondResetGap		= GetBLOBResetGap(node,1);
			currency			= GetCcy(node);
			dayCount			= GetAssetDayCount(node);
			resetFreq			= GetBLOBNum(node,3);

			// convertir en ARM
            char resetfrqStr[5];
			sprintf(resetfrqStr, "%i", resetFreq);
			resetFreq = FromSummitFreqToARMFreq(resetfrqStr);

			notional			= GetNotional(node);

			usedCorrel			= GetBLOBNum(node,5);

			// Initialisation des currency : initialement on considère un cas non quanto
			Idx1CcyTmp = (ARM_Currency*) currency->Clone();
			Idx2CcyTmp = (ARM_Currency*) currency->Clone();
			Idx3CcyTmp = (ARM_Currency*) currency->Clone();
			Idx4CcyTmp = (ARM_Currency*) currency->Clone();
			Idx5CcyTmp = (ARM_Currency*) currency->Clone();


			payIndexType = GetBLOBIdx(node,Idx1CcyTmp,5);
			Condition1ShortTermIndexType = GetBLOBIdx(node,Idx2CcyTmp,1);
			Condition1LongTermIndexType = GetBLOBIdx(node,Idx3CcyTmp,2);
			Condition2ShortTermIndexType = GetBLOBIdx(node,Idx4CcyTmp,3);
			Condition2LongTermIndexType = GetBLOBIdx(node,Idx5CcyTmp,4);

			condition1STIndexWeight=GetBLOBAmount(node,1);
			condition1LTIndexWeight=GetBLOBAmount(node,2);
			condition2STIndexWeight=GetBLOBAmount(node,3);
			condition2LTIndexWeight=GetBLOBAmount(node,4);

			

			payIndexCoeff=GetBLOBAmount(node,7);

			ARM_ReferenceValue aWeight1(condition1STIndexWeight);
			ARM_ReferenceValue aWeight2(condition1LTIndexWeight);
			ARM_ReferenceValue aWeight3(condition2STIndexWeight);
			ARM_ReferenceValue aWeight4(condition2LTIndexWeight);

			capOrFloorCondition2   = GetBLOBNum(node,2);
			capOrFloorCondition1   = GetBLOBNum(node,1);

			//Convert from Summit C/F to ARM C/F
			if (capOrFloorCondition2==1)
			{
				capOrFloorCondition2=K_CAP;
			}
			else if (capOrFloorCondition2==2)
			{
				capOrFloorCondition2=K_FLOOR;
			}

			if (capOrFloorCondition1==1)
			{
				capOrFloorCondition1=K_CAP;
			}
			else if (capOrFloorCondition1==2)
			{
				capOrFloorCondition1=K_FLOOR;
			}


			digitalBarrier = GetBLOBSchedRateMulti(node,"StartDate",1,2,100);
			spreadBarrier = GetBLOBSchedRateMulti(node,"StartDate",1,1,100);
			payIndexMargin = GetBLOBSchedRateMulti(node,"StartDate",1,3,100);

			refFixingCond1STIndex= GetBLOBSchedRateMultiWithLastDate(node,"FixingDate",date,2,1,100);
			refFixingCond1LTIndex= GetBLOBSchedRateMultiWithLastDate(node,"FixingDate",date,2,2,100);
			refFixingCond2STIndex= GetBLOBSchedRateMultiWithLastDate(node,"FixingDate",date,2,3,100);
			refFixingCond2LTIndex= GetBLOBSchedRateMultiWithLastDate(node,"FixingDate",date,2,4,100);
			refFixingPay = GetRefValSpreadFirstFixing(node,date);


			refPayFixedRate = new ARM_ReferenceValue(0.);



			//FIXED PAYMENT INDEX.
			if (payIndexType != K_FIXED)
			{
		
				if ((!payFreq == K_ZEROCOUPON))
				{
					term = 1.0 / ((double)payFreq);
				}

				paymentIndex   = new ARM_IRIndex(	dayCount,						//DayCount
													payFreq,						//ResetFreq
													payFreq,						//PayFreq
													term,							//Term
													K_COMP_PROP,					//CompoundMeth
													K_MOD_FOLLOWING,				//FwdRule
													K_ADJUSTED,						//IntRule
													fixCondTiming,					//ResetTiming
													10000.,							//ResetGap
													K_ARREARS,						//PayTiming
													10000.,							//PayGap
													currency,						//Currency
													(ARM_INDEX_TYPE) payIndexType,	//IndexType
													K_COMP_PROP);					//DecompFreq

				if(payIndexType == K_CMS)
				{
					ARM_INDEX_TYPE Libor_Type = currency->GetVanillaIndexType();
					paymentIndex->Set((ARM_INDEX_TYPE) payIndexType, Libor_Type, payFreq, K_COMP_PROP, currency);
					paymentIndex->Set(dayCount, payFreq, payFreq,
									term, K_COMP_PROP, K_MOD_FOLLOWING, K_ADJUSTED,
									fixCondTiming, 10000., K_ARREARS, 10000.,
									currency, (ARM_INDEX_TYPE) payIndexType, K_COMP_PROP);
				}
			}
			else
			{
				paymentIndex = new ARM_IRIndex(LIBOR3M, payFreq, payFreq, currency, dayCount);
				paymentIndex->SetResetTiming(K_ADVANCE);
				paymentIndex->SetIntRule(K_ADJUSTED);
			}	

			shiftCorrelVol = new ARM_VolFlat( tmpDate,0.0, currency);

			if (usedCorrel==0)
			{
				corridor = new ARM_CorridorDblCondition(startDate, 
																endDate,
																capOrFloorCondition2,
																capOrFloorCondition1,													 
																digitalBarrier,
																spreadBarrier,
																paymentIndex,							 
																refPayFixedRate,
																payIndexMargin,
																(ARM_INDEX_TYPE) Condition2ShortTermIndexType,
    															(ARM_INDEX_TYPE) Condition2LongTermIndexType,  
																(ARM_INDEX_TYPE) Condition1ShortTermIndexType,
																(ARM_INDEX_TYPE) Condition1LongTermIndexType,
																&aWeight3, 
																&aWeight4,
																&aWeight1,
																&aWeight2,
																refFixingCond2STIndex,
																refFixingCond2LTIndex,
																refFixingCond1STIndex,
																refFixingCond1LTIndex,
																refFixingPay,
  															    dayCount, 
																resetFreq, 
																payFreq, 
																fixCondTiming,
																payTiming, 
 																currency,
																fixCondResetGap,
																fabs(-0.01),
																fabs(0.01),
																intRule,
																stubRule,
																NULL,
																NULL,	
																payIndexCoeff,
																freezeFixing,
																10000,
																NULL,
																shiftCorrelVol,
																false);
			}

			if (usedCorrel==1)
			{
				correlValue=GetBLOBAmount(node,5);

				//On construit un objet market data : c'est mal
				ARM_VolFlat* createdVolCurveFlat = NULL;
	
				ARM_Date ADate; 
				ADate=date;
				createdVolCurveFlat = new ARM_VolFlat(ADate,correlValue, currency);

				corridor = new ARM_CorridorDblCondition(startDate, 
																endDate,
																capOrFloorCondition2,
																capOrFloorCondition1,													 
																digitalBarrier,
																spreadBarrier,
																paymentIndex,							 
																refPayFixedRate,
																payIndexMargin,
																(ARM_INDEX_TYPE) Condition2ShortTermIndexType,
    															(ARM_INDEX_TYPE) Condition2LongTermIndexType,  
																(ARM_INDEX_TYPE) Condition1ShortTermIndexType,
																(ARM_INDEX_TYPE) Condition1LongTermIndexType,
																&aWeight3, 
																&aWeight4,
																&aWeight1,
																&aWeight2,
																refFixingCond2STIndex,
																refFixingCond2LTIndex,
																refFixingCond1STIndex,
																refFixingCond1LTIndex,
																refFixingPay,
  															    dayCount, 
																resetFreq, 
																payFreq, 
																fixCondTiming,
																payTiming, 
 																currency,
																fixCondResetGap,
																fabs(-0.01),
																fabs(0.01),
																intRule,
																stubRule,
																NULL,
																NULL,	
																payIndexCoeff,
																freezeFixing,
																10000,
																(ARM_VolLInterpol*) createdVolCurveFlat,															
																shiftCorrelVol,
																false);

				if (createdVolCurveFlat)
				delete createdVolCurveFlat;
				createdVolCurveFlat = NULL;
			}

			corridor->SetAmount(notional);

			if (refFixingCond2STIndex)
				delete refFixingCond2STIndex;
				refFixingCond2STIndex = NULL;

			if (refFixingCond2LTIndex)
				delete refFixingCond2LTIndex;
				refFixingCond2LTIndex = NULL;

			if (refFixingCond1STIndex)
				delete refFixingCond1STIndex;
				refFixingCond1STIndex = NULL;

			if (refFixingCond1LTIndex)
				delete refFixingCond1LTIndex;
				refFixingCond1LTIndex = NULL;

			if (refFixingPay)
				delete refFixingPay;
				refFixingPay = NULL;

			if (paymentIndex)
				delete paymentIndex;
				paymentIndex = NULL;

			if (shiftCorrelVol)
				delete shiftCorrelVol;
				shiftCorrelVol = NULL;

			if (refPayFixedRate)
				delete refPayFixedRate;
				refPayFixedRate = NULL;

			if (digitalBarrier)
				delete digitalBarrier;
				digitalBarrier = NULL;

			if (spreadBarrier)
				delete spreadBarrier;
				spreadBarrier = NULL;

			if (payIndexMargin)
				delete payIndexMargin;
				payIndexMargin = NULL;

			if (notional)
				delete notional;
				notional = NULL;
				
		
			return corridor;;
		}
		catch (Exception& x) 
		{
	        x.DebugPrint();

			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "CorridorDoubleConditionAssetLoader::BuildAsset : Failed in constructing CorridorDblConditionLeg");
		} 
		catch (...) 
		{
			throw Exception (__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					 "Unrecognized Failure in : CorridorDoubleConditionAssetLoader::BuildAsset()");
		} 
	}
};

void AssetParserManager::Init(int isEtkOrNot) 
{
	if ( !initDone ) 
	{
		isEtkOrSummit = isEtkOrNot;

		//KEY is like "TYPE.PRODUCTNAME.INDEX"
		RegisterParser(string("SPDOPT"), new SpreadOptionAssetLoader());
	
		// Exotic screen:
		RegisterParser(string("EXOTIC..FIXED"),							new SwapAssetLoader());			// FixedLeg
		RegisterParser(string("EXOTIC.TEC_SPLK.FORM"),					new SwapAssetLoader());			// CMSLeg
		RegisterParser(string("EXOTIC..EUR1M"),							new SwapAssetLoader());			// LiborLeg
		RegisterParser(string("EXOTIC..EUR3M"),							new SwapAssetLoader());			// LiborLeg
		RegisterParser(string("EXOTIC..EURIB"),							new SwapAssetLoader());			// LiborLeg
		RegisterParser(string("EXOTIC..EUR12"),							new SwapAssetLoader());			// LiborLeg
		RegisterParser(string("EXOTIC..FORM"),							new SwapAssetLoader());			// LiborLeg
		RegisterParser(string("EXOTIC..LIBOR"),							new SwapAssetLoader());			// LiborLeg USD
		RegisterParser(string("EXOTIC.RNGLIBORFIX.FORM"),				new CorridorAssetLoader());		// CorridorLeg
		RegisterParser(string("EXOTIC.RNGLIBORFIX"),					new CorridorAssetLoader());		// CorridorLeg
		RegisterParser(string("EXOTIC.RNGLIBORLIBOR.FORM"),				new CorridorAssetLoader());		// CorridorLeg
		RegisterParser(string("EXOTIC.RNGLIBORLIBOR"),					new CorridorAssetLoader());		// CorridorLeg
		RegisterParser(string("EXOTIC.RESTRIKABLECMSLIBOR.FORM"),		new CorridorAssetLoader());		// RestrikableLeg
		RegisterParser(string("EXOTIC.RESTRIKABLECMSCMS.FORM"),			new CorridorAssetLoader());		// RestrikableLeg
		RegisterParser(string("EXOTIC.RESTRIKABLECMSFIX.FORM"),			new CorridorAssetLoader());		// RestrikableLeg
		RegisterParser(string("EXOTIC.RNG_DOUBLE"),						new CorridorDoubleConditionAssetLoader());		// RestrikableLeg

		// Cap/Floor screen:
		RegisterParser(string("IRG.SPREADOPTIONLOG.FORM"),				new SpreadOptionAssetLoader());	// Spreadoption
		RegisterParser(string("IRG.SPREADOPTIONLOGDIGITAL.FORM"),		new SpreadOptionAssetLoader());	// Spreadoption
		RegisterParser(string("IRG.SPREADOPTIONLOGFLTDIGITAL.FORM"),	new SpreadOptionAssetLoader());	// Spreadoption
		RegisterParser(string("IRG.SPREADOPTIONLOGCORRIDOR.FORM"),		new SpreadOptionAssetLoader());	// Spreadoption


		RegisterParser(string("IRG.MATURITYCAP.FORM"),					new MaturityCapAssetLoader());	// Maturity Cap
		RegisterParser(string("IRG.RFTARN.FORM"),						new TarnAssetLoader());			// Tarn
		RegisterParser(string("IRG.ALMCAPTION.FORM"),					new CaptionAssetLoader());		// Caption

		RegisterParser(string("IRG..EURIB"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG..EUR3M"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG..EUR12"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG..LIBOR"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG..STIBO"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG..FORM"),								new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.CUO_SMILE.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.CUO_SPLK.FORM"),						new CapFloorAssetLoader());		//CapFloor CMS
		RegisterParser(string("IRG.CUI_SMILE.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.CUI_SPLK.FORM"),						new CapFloorAssetLoader());		//CapFloor CMS
		RegisterParser(string("IRG.PDO_SMILE.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.PDO_SPLK.FORM"),						new CapFloorAssetLoader());		//CapFloor CMS
		RegisterParser(string("IRG.PDI_SMILE.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.PDI_SPLK.FORM"),						new CapFloorAssetLoader());		//CapFloor CMS
		RegisterParser(string("IRG.TEC_SPLK_CF.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.RCMS.FORM"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.RCMSD.FORM"),						new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.RCMT.FORM"),							new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.RCMTD.FORM"),						new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.LIVRET_A.FORM"),						new CapFloorAssetLoader());		//CapFloor LivretA
		RegisterParser(string("IRG.LIVRET_A_ARRONDI.FORM"),				new CapFloorAssetLoader());		//CapFloor LivretA
		RegisterParser(string("IRG.GLOBALCAP.FORM"),					new CapFloorAssetLoader());		//CapFloor
		RegisterParser(string("IRG.FX_CALL_STRIP.FORM"),				new FXOptionStripAssetLoader());//StripOption
		RegisterParser(string("IRG.FX_CALL_STRIP_DIGITAL.FORM"),		new FXOptionStripAssetLoader());//StripDigitalOption	
		
		RegisterParser(string("IRG.DIGITAL_SMILE.FORM"),				new DigitalAssetLoader());		//Digital
		RegisterParser(string("IRG.DIGITAL_SPLK.FORM"),					new DigitalAssetLoader());		//Digital CMS

		RegisterParser(string("IRG.FX_QUANTO_STRIP.FORM"),				new FXStripQuantoAssetLoader());// FX_Quanto_Strip
		RegisterParser(string("IRG.FX_QUANTO_STRIP_DIGITAL.FORM"),		new FXStripQuantoAssetLoader());// FX_Quanto_Strip_Digital
		RegisterParser(string("IRG.FX_SPREAD_STRIP.FORM"),				new FXStripQuantoAssetLoader());// FX_Spread_Strip
		RegisterParser(string("IRG.FX_SPREAD_STRIP_DIGITAL.FORM"),		new FXStripQuantoAssetLoader());// FX_Spread_Strip_Digital

		RegisterParser(string("C.IRG.FX_QUANTO_STRIP.FORM"),			new FXStripQuantoCalculatorLoader()); // FX_Quanto_Strip
		RegisterParser(string("C.IRG.FX_QUANTO_STRIP_DIGITAL.FORM"),	new FXStripQuantoCalculatorLoader()); // FX_Quanto_Strip_Digital
		RegisterParser(string("C.IRG.FX_SPREAD_STRIP.FORM"),			new FXStripQuantoCalculatorLoader()); // FX_Spread_Strip
		RegisterParser(string("C.IRG.FX_SPREAD_STRIP_DIGITAL.FORM"),	new FXStripQuantoCalculatorLoader()); // FX_Spread_Strip_Digital

		// Memory
		RegisterParser(string("IRG.MEMORYCAP.FORM"),					new MemoryCapFloorAssetLoader());		//MemoryCapFloor
		RegisterParser(string("IRG.MEMORYSO.FORM"),						new MemorySpreadOptionAssetLoader());	//MemorySO

		// Swap screen: 
		RegisterParser(string("SWAP..FIXED"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..EURIB"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..LIBOR"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..EUR3M"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..EUR1M"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..EUR12"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..STIBO"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP..FORM"),							new SwapAssetLoader());			// Swap
		RegisterParser(string("SWAP.RCMS.FORM"),						new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("SWAP.RCMSD.FORM"),						new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("SWAP.RCMT.FORM"),						new SwapAssetLoader());			// Swap CMT
		RegisterParser(string("SWAP.RCMTD.FORM"),						new SwapAssetLoader());			// Swap CMT
		RegisterParser(string("SWAP.TEC_SPLK.FORM"),					new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("SWAP.LIVRET_A.FORM"),					new SwapAssetLoader());			// Swap LivretA
		RegisterParser(string("SWAP.LIVRET_A_ARRONDI.FORM"),			new SwapAssetLoader());			// Swap LivretA
		RegisterParser(string("EXOTIC.RCMS.FORM"),						new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("EXOTIC.RCMSD.FORM"),						new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("EXOTIC.RCMT.FORM"),						new SwapAssetLoader());			// Swap CMT
		RegisterParser(string("EXOTIC.RCMTD.FORM"),						new SwapAssetLoader());			// Swap CMT
		RegisterParser(string("EXOTIC.TEC_SPLK.FORM"),					new SwapAssetLoader());			// Swap CMS
		RegisterParser(string("EXOTIC.LIVRET_A.FORM"),					new SwapAssetLoader());			// Swap LivretA
		RegisterParser(string("EXOTIC.LIVRET_A_ARRONDI.FORM"),			new SwapAssetLoader());			// Swap LivretA

		RegisterParser(string("IRG.RNGLIBORFIX.FORM"),					new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("IRG.RNGLIBORFIX"),						new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("IRG.RNGLIBORLIBOR"),						new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("IRG.RNGLIBORLIBOR.FORM"),				new CorridorAssetLoader());		// Corridor Leg

		RegisterParser(string("SWAP.RNGLIBORFIX.FORM"),					new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("SWAP.RNGLIBORFIX"),						new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("SWAP.RNGLIBORLIBOR"),					new CorridorAssetLoader());		// Corridor Leg
		RegisterParser(string("SWAP.RNGLIBORLIBOR.FORM"),				new CorridorAssetLoader());		// Corridor Leg

		RegisterParser(string("SWAP.RNG_DOUBLE.FORM"),					new CorridorDoubleConditionAssetLoader());	// Corridor Double Condition Leg

		initDone = true;
	}
	isEtkOrSummit = isEtkOrNot;
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

