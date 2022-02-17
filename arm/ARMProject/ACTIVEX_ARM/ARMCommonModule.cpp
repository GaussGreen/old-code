//	
//		CAUTION 
//
//		This file is the COM DELEGATION implementation
//
//		It is compiled by the ARMFILib and ARMCCSLib projects
// 

#include "stdafx.h"
#ifdef _ACTIVEX_FI_
	#include "..\ARMFILib\ARMFILib.h"
	#include "..\ARMFILib\ARMFIModule.h"
#endif 

#ifdef _ACTIVEX_CCS_
	#include "..\ARMCCSLib\ARMCCSLib.h"
	#include "..\ARMCCSLib\ARMCCSModule.h"
#endif 

#ifdef _ACTIVEX_OLD_
	#include "ARMModule.h"
#endif 

#include "..\ActiveX_ARM\ActiveXModule.h"

static ActiveXModule theInstance; 


//	----------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_Credit_CDS(double	pEffectiveDate,
									   double	pEndDate, 
									   double	pSpread, 
									   BSTR		pFixingFreq, 
									   BSTR		pDayCount, 
									   double	pFirst_period_refdate, 
									   double	pFixedPayerAmount, 
									   double	pFloatingPayerAmount, 
									   BSTR		StubRule,
									   BSTR     pCurrency,
									   BSTR		Adjusted,
									   int		l_CreditLag,
									   BSTR	    pIncludeMaturity,	   
									   double	pProtectionStartDate,
									   double	pProtectionEndDate, 
									   BSTR		IssuerName,
									   double   Binary,
									   double	pFstCpnEffDate,
									   BSTR		StartAdj,
									   VARIANT *pRet)
{
	return theInstance.ARM_Credit_CDS(	pEffectiveDate,
									   	pEndDate, 
									   	pSpread, 
									   		pFixingFreq, 
									   		pDayCount, 
									   	pFirst_period_refdate, 
									   	pFixedPayerAmount, 
									   	pFloatingPayerAmount, 
									   		StubRule,
									        pCurrency,
									   		Adjusted,
									   		l_CreditLag,
									   	    pIncludeMaturity,
									   	pProtectionStartDate,
									   	pProtectionEndDate, 
									   		IssuerName,
									      Binary,
									   	pFstCpnEffDate,
											StartAdj,
									   pRet); 
}

//	----------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_Credit_FTD(VARIANT *pEffectiveDate,
									   VARIANT *pEndDate, 
									   VARIANT *pSpread, 
									   VARIANT *pLabels, 
									   VARIANT *pFixingFreq, 
									   VARIANT *pDayCountFrq, 
									   VARIANT *pFirst_period_refdate, 
									   VARIANT *pIssuerNotional, 
									   VARIANT *pAccruedOnDefault, 	
									   VARIANT *pCurrency,
									   VARIANT *pPayCreditLag,
									   VARIANT *pStub,	
									   double	pFstCpnEffDate,
									   VARIANT *pintRule,
									   VARIANT *pstartAdj,
									   VARIANT *pRet)
{
	return theInstance.ARM_Credit_FTD(pEffectiveDate,
									   pEndDate, 
									   pSpread, 
									   pLabels, 
									   pFixingFreq, 
									   pDayCountFrq, 
									   pFirst_period_refdate, 
									   pIssuerNotional, 
									   pAccruedOnDefault, 	
									   pCurrency,
									   pPayCreditLag,
									   pStub,	
									   pFstCpnEffDate,
									   pintRule,
									   pstartAdj,
									   pRet); 
	
}

//	----------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_Credit_NTD(double pEffectiveDate,
									   double pEndDate, 
									   double pSpread, 
									   int	  pFirstNumDefault, 	
									   int	  pLastNumDefault, 	
									   VARIANT *pLabels, 	
									   BSTR		pFixingFreq, 
									   BSTR		pDayCountFrq, 
									   double	pFirst_period_refdate, 
									   double	pIssuerNotional, 
									   BSTR	    pAccruedOnDefault, 	
									   BSTR		pCurrency,
									   double	pPayCreditLag,	
									   BSTR		pStub,
									   BSTR		pFrequencyDefLeg,
									   double	pBinary,
									   BSTR		pPayCal,	
									   BSTR		pRcvFeeLeg,
									   double   TradedNotional,
									   BSTR		InclMaturity,
									   double	pFstCpnEffDate,
									   BSTR		intRule,
									   BSTR		startAdj,
									   VARIANT *pRet)
{
return theInstance.ARM_Credit_NTD(  pEffectiveDate,
								    pEndDate, 
								    pSpread, 
									pFirstNumDefault, 	
									pLastNumDefault, 	
									pLabels, 	
									pFixingFreq, 
									pDayCountFrq, 
									pFirst_period_refdate, 
									pIssuerNotional, 
									pAccruedOnDefault, 	
									pCurrency,
									pPayCreditLag,	
									pStub,
									pFrequencyDefLeg,
									pBinary,
									pPayCal,	
									pRcvFeeLeg,
									TradedNotional,
									InclMaturity,
									pFstCpnEffDate,
									intRule,
									startAdj,
									pRet) ;
}

//	----------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_Credit_Mezzanine(double pEffectiveDate,
											 double pEndDate, 
											 double pSpread, 
											 double pMezzAmount, 
											 double pSubAmount,
											 VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 BSTR	pFreqFeeLeg, 
											 BSTR	pDayCountFrq, 
											 double pFirst_period_refdate, 
											 BSTR	pAccruedOnDefault, 	
											 BSTR	pCurrency,
											 double pPayCreditLag,	
											 BSTR	pStub,	
											 BSTR	pFreqDefLeg, 
											 double	pBinary, 
											 BSTR	pPayCal,	
											 BSTR	pRcvFeeLeg,
											 double TradedNotional,
											 BSTR	InclMaturity,
											 double	pFstCpnEffDate,
											 BSTR intRule,
											 BSTR adjStartDate,
											 VARIANT *pRet)
{
return 
theInstance.ARM_Credit_Mezzanine( pEffectiveDate,
											  pEndDate, 
											  pSpread, 
											  pMezzAmount, 
											  pSubAmount,
											 pLabels, 
											 pNotionals, 
											 	pFreqFeeLeg, 
											 	pDayCountFrq, 
											  pFirst_period_refdate, 
											 	pAccruedOnDefault, 	
											 	pCurrency,
											  pPayCreditLag,	
											 	pStub,	
											 	pFreqDefLeg, 
											 	pBinary, 
											 	pPayCal,	
											 	pRcvFeeLeg,
											  TradedNotional,
											 	InclMaturity,
											 	pFstCpnEffDate,
												intRule,
											 adjStartDate,
											 pRet) ;
}

//	----------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_Credit_CDO2(double pEffectiveDate,
										double pEndDate,
										BSTR   pPortfolio,
										double pSpread,
										double pMezzAmount,
										double pSubAmount,
										BSTR pFreqFeeLeg,
										BSTR pFreqDefLeg,
										BSTR pDayCountFrq,
										double pFirst_period_refdate,
										BSTR	pAccruedOnDefault,
										BSTR	pCurrency,
										double	pPayCreditLag,
										BSTR	pStub,
										double	pBinary,
										BSTR	pPayCal,	
										BSTR	pRcvFeeLeg,
										double  TradedNotional,
										BSTR	CrossSubordination,
										BSTR	InclMaturity,
										double	pFstCpnEffDate,
										VARIANT *pRet)
{
	return theInstance.ARM_Credit_CDO2( pEffectiveDate,
										 pEndDate,
										   pPortfolio,
										 pSpread,
										 pMezzAmount,
										 pSubAmount,
										 pFreqFeeLeg,
										 pFreqDefLeg,
										 pDayCountFrq,
										 pFirst_period_refdate,
											pAccruedOnDefault,
											pCurrency,
											pPayCreditLag,
											pStub,
											pBinary,
											pPayCal,	
											pRcvFeeLeg,
										  TradedNotional,
											CrossSubordination,
											InclMaturity,
											pFstCpnEffDate,
										pRet) ;
}
//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule:: ARM_Credit_Portfolio(VARIANT *pSecuritiesID, 
											  BSTR pParameters,
										   	  VARIANT *pRet)
{
return theInstance. ARM_Credit_Portfolio(pSecuritiesID, 
											   pParameters,
										   	  pRet) ;
}


//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule::ARM_Credit_CMTranche(VARIANT *pEffectiveDate,
											 VARIANT *pEndDate, 
											 VARIANT *pSpread, 
											 VARIANT *pMezzAmount, 
											 VARIANT *pSubAmount,
											 VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 VARIANT *pIndex,
											 VARIANT *pFreqFeeLeg, 
											 VARIANT *pDayCountFrq, 
											 VARIANT *pFirst_period_refdate, 
											 VARIANT *pAccruedOnDefault, 	
											 VARIANT *pCurrency,
											 VARIANT *pPayCreditLag,	
											 VARIANT *pStub,	
											 VARIANT *pFreqDefLeg, 
											 VARIANT *pBinary, 
											 VARIANT *pPayCal,	
											 BSTR	 pRcvFeeLeg,
											 double  TradedNotional,
											 double  FwdFixedDate,
											 BSTR	 InclMaturity,
											 double	 pFstCpnEffDate,
											 VARIANT *pRet)
{
	return  theInstance.ARM_Credit_CMTranche(pEffectiveDate,
											 pEndDate, 
											 pSpread, 
											 pMezzAmount, 
											 pSubAmount,
											 pLabels, 
											 pNotionals, 
											 pIndex,
											 pFreqFeeLeg, 
											 pDayCountFrq, 
											 pFirst_period_refdate, 
											 pAccruedOnDefault, 	
											 pCurrency,
											 pPayCreditLag,	
											 pStub,	
											 pFreqDefLeg, 
											 pBinary, 
											 pPayCal,	
											 pRcvFeeLeg,
											 TradedNotional,
											 FwdFixedDate,
											 InclMaturity,
											 pFstCpnEffDate,
											 pRet) ;

}
//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule::ARM_Credit_Index(VARIANT *pLabels,
										 double YearFrac,
										 double Spread, 
										 BSTR Method, 
										 BSTR Basis,
										 BSTR ResetFreq,
										 BSTR PayFreq, 
										 BSTR ccy,
										 BSTR DefCurve, 
										 BSTR fwdRule,
										 BSTR resetTiming,
										 int resetGap,
										 BSTR payTiming,
										 int payGap,
										 BSTR intRule,
										 BSTR AdjCalType,
										 int cm_resetWeekDay,
										 int cm_resetOccur, 
										 BSTR *pRet)
{
return theInstance.ARM_Credit_Index(pLabels,
										 YearFrac,
										 Spread, 
										 Method,
										 Basis,
										 ResetFreq,
										 PayFreq, 
										 ccy,
										 DefCurve,  
										 fwdRule,
										 resetTiming,
										 resetGap,
										 payTiming,
										 payGap,
										 intRule,
										 AdjCalType,
										 cm_resetWeekDay,
										 cm_resetOccur, 
										 pRet) ;
}

STDMETHODIMP ARMCommonModule::ARM_Credit_IndexCompo(BSTR IndexName,
										VARIANT *pLabels,
										 VARIANT * YearFrac,		 			 
										 VARIANT * Spread, 
										 BSTR Method, 			 
										 BSTR Basis,
										 BSTR ResetFreq,
										 BSTR PayFreq,
										 BSTR ccy, 
										 BSTR fwdRule,
										 BSTR resetTiming,
										 int resetGap,
										 BSTR payTiming,
										 int payGap,
										 BSTR intRule,
										 BSTR AdjCalType,
										 int cm_resetWeekDay,
										 int cm_resetOccur, 
										 BSTR *pRet)
{
return theInstance.ARM_Credit_IndexCompo(IndexName,
										 pLabels,
										 YearFrac,										 
										 Spread, 
										 Method, 
										 Basis,
										 ResetFreq,
										 PayFreq, 
										 ccy,
										 fwdRule,
										 resetTiming,
										 resetGap,
										 payTiming,
										 payGap,
										 intRule,
										 AdjCalType,
										 cm_resetWeekDay,
										 cm_resetOccur, 
										 pRet) ;
}

//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule::ARM_Credit_CDSIndex(VARIANT *pEffectiveDate,
											VARIANT *pEndDate, 
											VARIANT *pSpread,
											VARIANT *pIndex,
											VARIANT *pFixingFreq, 
											VARIANT *pDayCountFrq, 
											VARIANT *pFirst_period_refdate, 
											VARIANT *pFixedPayerAmount, 
											VARIANT *pFloatingPayerAmount, 
											BSTR	StubRule,
											VARIANT *pCurrency,
											BSTR	Adjusted,
											int		pCreditLag,											
											BSTR	pIncludeMaturity,
											VARIANT *pProtectionStartDate,
											VARIANT *pProtectionEndDate, 
											VARIANT *pRet)
{
	return theInstance.ARM_Credit_CDSIndex(pEffectiveDate,
											pEndDate, 
											pSpread,
											pIndex,
											pFixingFreq, 
											pDayCountFrq, 
											pFirst_period_refdate, 
											pFixedPayerAmount, 
											pFloatingPayerAmount, 
											StubRule,
											pCurrency,
											Adjusted,
											pCreditLag,											
											pIncludeMaturity,
											pProtectionStartDate,
											pProtectionEndDate, 
											pRet) ;

}
//	----------------------------------------------------------------------------------
//


STDMETHODIMP ARMCommonModule::ARM_Credit_Option(VARIANT *  pUnderlyingMaturity,
												double OptionExpiry,
												BSTR Currency,
												BSTR CdsAdj,
												BSTR EndAdj,
												double pStrike,
											    VARIANT *pOptionType,
												VARIANT *pKoType,
												double notional,

												VARIANT *pRet
										  )
{


	return theInstance.ARM_Credit_Option(pUnderlyingMaturity,
										OptionExpiry,
										Currency,
										CdsAdj,
										EndAdj,
										pStrike,
										pOptionType,
										pKoType,
										notional,
										pRet) ;
}

//  ---------------------------------------------------------------------------------
//
STDMETHODIMP ARMCommonModule::ARM_GETINSTRUMENTFROMCALYPSO(BSTR CalypsoId,
													BSTR Type,
													double AsOf,
													BSTR ModelType,
													VARIANT* pRet)
{
	return theInstance.ARM_GETINSTRUMENTFROMCALYPSO(CalypsoId,
													Type,
													 AsOf,
													 ModelType,
													 pRet);
}

//  ---------------------------------------------------------------------------------
//

//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule::ARM_GETINSTRUMENTFROMSUMMIT(BSTR SummitId,
													BSTR Type,
													double AsOf,
													BSTR Exoticfilter,
													VARIANT* pRet)
{
	return theInstance.ARM_GETINSTRUMENTFROMSUMMIT(SummitId,
													Type,
													AsOf,
													Exoticfilter,
													pRet) ;
	
}

//	----------------------------------------------------------------------------------
//

STDMETHODIMP ARMCommonModule::ARM_Credit_EmptyLeg(VARIANT *pRet)
{
	return theInstance.ARM_Credit_EmptyLeg(pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_IRLEGTOCREDITLEG(BSTR SwapLegId,
													BSTR LegType,
													BSTR creditindexId,
													BSTR PricerId,
													BSTR *pRet)
{
	return theInstance.ARM_Credit_IRLEGTOCREDITLEG(SwapLegId,
													LegType,
													creditindexId,
													PricerId,
													pRet);
}



STDMETHODIMP ARMCommonModule::ARM_Credit_Collateral(VARIANT *pLabels, 
											 VARIANT *pNotionals, 
											 VARIANT *pRet)
{
return theInstance.ARM_Credit_Collateral(pLabels, 
											 pNotionals, 
											 pRet) ;

}

STDMETHODIMP ARMCommonModule::ARM_Credit_VariableCollateral(VARIANT pLabels, 
											 VARIANT pNotionals, 
											 VARIANT *pRet)
{
return theInstance.ARM_Credit_VariableCollateral(pLabels, 
										 pNotionals, 
											 pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_CDSGEN(BSTR FeeLegId,
										  BSTR DefLegId,
										  double RcvFee,
										  double TradedNot,
										  BSTR *pRet)
{
	return theInstance.ARM_Credit_CDSGEN(FeeLegId,
										  DefLegId,
										  RcvFee,
										  TradedNot,
										  pRet) ;

}

STDMETHODIMP ARMCommonModule::ARM_Credit_NTDGEN(BSTR CdsId,
										int firstnumdef,
										int lastnumdef,
										BSTR CollateralId,
										double binary,
										double rcvfee,
										BSTR *pRet)
{
	return theInstance.ARM_Credit_NTDGEN(CdsId,
										firstnumdef,
										lastnumdef,
										CollateralId,
										binary,
										rcvfee,
										pRet) ; 


}

STDMETHODIMP ARMCommonModule::ARM_Credit_CDOGEN(BSTR CdsId,
										double subamount,
										BSTR CollateralId,
										double binary,
										double rcvfee,
										BSTR *pRet)
{
	return theInstance.ARM_Credit_CDOGEN(CdsId,
										subamount,
										CollateralId,
										binary,
										rcvfee,
										pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_CDO_SQUARE_GEN(BSTR CdsId,
												  double subamount,
												  BSTR PortfolioId,
												  double binary,
												  double rcvfee,
												  BSTR *pRet)
{
	return theInstance.ARM_Credit_CDO_SQUARE_GEN(CdsId,
												 subamount,
												 PortfolioId,
												 binary,
												 rcvfee,
												 pRet);

}



STDMETHODIMP ARMCommonModule::ARM_Credit_GenLeg(double	StartDate,
										  double	EndDate, 
										  double	FixedRate, 
										  double	FixedNotional, 
										  BSTR		RefValNotional, 
										  BSTR		RefValRate, 
										  BSTR		XChangeNotional, 
										  BSTR		Frequency, 
									      BSTR		Basis,
									      BSTR      payTiming,
									      BSTR		intrule,
									      BSTR		stubrule,
									      BSTR	    ccyid,
									      BSTR		paycalname,
									      double	refdate, 
									      BSTR		includematurity,
									      BSTR	    adjstartdate,
										  BSTR	    legtype,
										  BSTR		indexobj,
										  int		creditlag,
										  double    binary,
										  BSTR      Name,
										  BSTR      Nxchange,
										  BSTR		baccruedOnDef,
									      VARIANT	*pRet)
{
	return theInstance.ARM_Credit_GenLeg(	StartDate,
						EndDate, 
						FixedRate, 
						FixedNotional, 
						RefValNotional, 
						RefValRate, 
						XChangeNotional, 
						Frequency, 
						Basis,
						payTiming,
						intrule,
						stubrule,
						ccyid,
						paycalname,
						refdate, 
						includematurity,
						adjstartdate,
						legtype,
						indexobj,
						creditlag,
						binary,
						Name,
						Nxchange,
						baccruedOnDef,
						pRet) ;
}


//	------------------------------------------------------------------------------------
STDMETHODIMP ARMCommonModule::ARM_Credit_GetLastFixingDate(
										BSTR instId_,
										VARIANT xlAsofDate,
										VARIANT*pRet)
{
	return theInstance.ARM_Credit_GetLastFixingDate(
										 instId_,
										 xlAsofDate,
										pRet) ;
}

//	------------------------------------------------------------------------------------
STDMETHODIMP ARMCommonModule::ARM_Credit_SetPastFixing(
									BSTR instId_,
									VARIANT resetDate,
									VARIANT fixingValue,
									VARIANT*pRet)
{
	return theInstance.ARM_Credit_SetPastFixing(
									 instId_,
									 resetDate,
									 fixingValue,
									pRet) ;

}

STDMETHODIMP ARMCommonModule::ARM_Credit_GetBounds(
									BSTR instId_,
									double* low,
									double* up)
{
	return theInstance.ARM_Credit_GetBounds(
									 instId_,
									 low,
									 up) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_Customized_CDO(
									VARIANT* pLabels, 
									VARIANT *pNotionals,
									BSTR pCurrency,
									BSTR pDefaultLeg, 
									BSTR pPremiumLeg, 
									BSTR pParameters, 
									BSTR *pRet)
{

	return theInstance.ARM_Credit_Customized_CDO(
									 pLabels, 
									 pNotionals,
									 pCurrency,
									 pDefaultLeg, 
									 pPremiumLeg, 
									 pParameters, 
									pRet) ;


}


STDMETHODIMP ARMCommonModule::ARM_Credit_CLN(double	pEffectiveDate,
									   double	pEndDate, 
									   double	pSpread, 
									   BSTR		pIndexId, 
									   double	pFirst_period_refdate, 
									   double	pFstCpnEffDate,
									   double	pNotional, 
									   BSTR		AccOnDef,
									   BSTR		pDayCount, 
									   BSTR		pDecompFreq, 
									   BSTR		StubRule,
									   double	resetgap,
									   BSTR     pCurrency,
									   BSTR     ResetCal,
									   BSTR     PayCal,
									   BSTR     Nxchange,
									   BSTR	    pIncludeMaturity,
									   BSTR		AdjustedStartDate,
									   double   Binary,
									   BSTR		*pRet)
{
	return theInstance.ARM_Credit_CLN(	pEffectiveDate,
									   	pEndDate, 
									   	pSpread, 
									   		pIndexId, 
									   	pFirst_period_refdate, 
									   	pFstCpnEffDate,
									   	pNotional, 
									   		AccOnDef,
									   		pDayCount, 
									   		pDecompFreq, 
									   		StubRule,
									   	resetgap,
									        pCurrency,
									        ResetCal,
									        PayCal,
									        Nxchange,
									   	    pIncludeMaturity,
									   		AdjustedStartDate,
									      Binary,
									   pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetInitialCurveFromSummit(
															   VARIANT *pIndex, 
															   VARIANT *pCurrency, 
															   VARIANT *pCvName, 
															   VARIANT *pDate, 
															   VARIANT* pTypeValue, 
															   VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetInitialCurveFromSummit(pIndex, 
		pCurrency, 
		pCvName, 
		pDate, 
		pTypeValue, 
		pRet);
}



STDMETHODIMP ARMCommonModule::ARM_Credit_GetDPFromSummit(double AsOfDate,
												BSTR Issuer, 
												BSTR CurveName, 
												BSTR PWCcurveId, 
												BSTR label, 
												VARIANT *pRet)
{
return theInstance.ARM_Credit_GetDPFromSummit( AsOfDate,
												 Issuer, 
												 CurveName, 
												 PWCcurveId, 
												 label, 
												 pRet);
}


STDMETHODIMP ARMCommonModule::ARM_Credit_DPMktDataFromSummit(double AsOfDate,
												BSTR Issuer, 
												BSTR CurveName, 
												BSTR Parameter, 
												VARIANT *pRet)
{
return theInstance.ARM_Credit_DPMktDataFromSummit( AsOfDate,
												 Issuer, 
												 CurveName, 
												 Parameter, 
												 pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_DPMktDataFromCalypso(double AsOfDate,
												BSTR pricingEnv,
												BSTR issuer, 
												BSTR seniority, 
												BSTR ccy, 
												BSTR forceCurveName,
												BSTR xmlFile,
												VARIANT *Parameter, 
												VARIANT *pRet)
{
return theInstance.ARM_Credit_DPMktDataFromCalypso( AsOfDate,
												 pricingEnv, 
												 issuer, 
												 seniority,
												 ccy,
												 forceCurveName,
												 xmlFile,
												 Parameter,
												 pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_GetDPFromCalypso(double pDate,
										  BSTR pricingEnv, 
										  BSTR issuer,
										  BSTR seniority,
										  BSTR ccy,
										  BSTR forceCurveName,
										  BSTR xmlFile,
										  BSTR irCurveId,
										  BSTR label,
										  VARIANT *ret)
{
	return theInstance.ARM_Credit_GetDPFromCalypso(pDate,
											pricingEnv,
											issuer,
											seniority,
											ccy,
											forceCurveName,
											xmlFile,
											irCurveId,
											label,
											ret);

}



STDMETHODIMP ARMCommonModule::ARM_Credit_ConstantDefaultCurve(DATE AsOfDate,
														VARIANT pTenors, 
														VARIANT pRates, 
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	bAdjCalType, 
														BSTR	bIsSummit, 
														BSTR	calibrationData,
														int lag,
														BSTR	calibrationAlgo,
														BSTR	paramId,
														BSTR	intRule,
														BSTR	adjStartRule,
														VARIANT *pRet)
{
	return theInstance.ARM_Credit_ConstantDefaultCurve( AsOfDate,
														pTenors, 
														pRates, 
														Recovery, 
														IRCurveId,
														bCurrency, 
														bLabel, 
														bAdjCalType, 
														bIsSummit, 
														calibrationData,
														lag,
														calibrationAlgo,
														paramId,
														intRule,
														adjStartRule,
														pRet) ;
}


STDMETHODIMP ARMCommonModule::ARM_Credit_ZeroCouponDefaultCurveFromSummit(double AsOfDate,
																	BSTR   bIssuer, 
																	BSTR   bCurrency, 
																	BSTR   bCvName,	
															        BSTR   IRCurveId,
																	BSTR   bLabel,
																	VARIANT *pRet)
{
return 
theInstance.ARM_Credit_ZeroCouponDefaultCurveFromSummit(AsOfDate,
														bIssuer, 
														bCurrency, 
														bCvName,	
														IRCurveId,
														bLabel,
														pRet); 
}
STDMETHODIMP ARMCommonModule::ARM_Credit_ZeroCouponDefaultCurveFromCalypso(double AsOfDate,
												BSTR pricingEnv,
												BSTR issuer, 
												BSTR seniority, 
												BSTR ccy, 
												BSTR forceCurveName,
												BSTR xmlFile,
												BSTR irCurveId, 
												BSTR label, 
												VARIANT *pRet)
{
return theInstance.ARM_Credit_ZeroCouponDefaultCurveFromCalypso( AsOfDate,
												 pricingEnv, 
												 issuer, 
												 seniority,
												 ccy,
												 forceCurveName,
												 xmlFile,
												 irCurveId,
												 label,
												 pRet);
}







STDMETHODIMP ARMCommonModule::ARM_Credit_InputDefaultCurve(VARIANT AsOfDate,
														VARIANT  pDates, 
														VARIANT  pRates, 
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														// BSTR	bInterpolType, 
														VARIANT *pRet)
{
	return theInstance.ARM_Credit_InputDefaultCurve(AsOfDate,
													pDates, 
													pRates, 
													Recovery, 
													IRCurveId,
													bCurrency, 
													bLabel, 
													// bInterpolType, 
													pRet);
}

STDMETHODIMP ARMCommonModule::ARM_FxConvert(VARIANT *pccy1, 
									  VARIANT *pccy2, 
									  VARIANT *pAsOfDate,
									  VARIANT *pCvName,
									  VARIANT *pRet)
{
	return theInstance.ARM_FxConvert(pccy1, 
									 pccy2, 
									 pAsOfDate,
									 pCvName,
									 pRet); 
	
}


STDMETHODIMP ARMCommonModule::ARM_FxConvertFromCalypso(BSTR ccy1, 
									  BSTR ccy2, 
									  BSTR cvName,
									  DATE pAsOfDate,
									  VARIANT *pRet)
{
	return theInstance.ARM_FxConvertFromCalypso(ccy1, 
									 ccy2, 
									 cvName,
									 pAsOfDate,
									 pRet); 
	
}

STDMETHODIMP ARMCommonModule::ARM_DiscountPrice(VARIANT *pCurve, VARIANT *pMatu, VARIANT *pRet)
{
	return theInstance.ARM_DiscountPrice(pCurve, pMatu, pRet);
}


// -----------------------------------------------------------------------------------------------
// Fonctionalités Credit
// -----------------------------------------------------------------------------------------------



STDMETHODIMP ARMCommonModule::ARM_Credit_Delivery(VARIANT *pAsOfDate, 
											VARIANT *pTenorContract, 
											VARIANT *pRet)
{
	return theInstance.ARM_Credit_Delivery(pAsOfDate, 
											pTenorContract, 
											pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_CptInterpolDefCurve(BSTR pCurve,VARIANT pTENOR,VARIANT *pRet)
{
	return theInstance.ARM_Credit_CptInterpolDefCurve( pCurve, pTENOR,pRet); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_createFlatCurve(BSTR pCurve,VARIANT pTENOR,VARIANT *pRet)
{
	return theInstance.ARM_Credit_createFlatCurve( pCurve, pTENOR,pRet); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_createDefCurveFromBase(BSTR pCurveCDS,BSTR pCurveIndex,VARIANT vBase,VARIANT *pRet)
{
	return theInstance.ARM_Credit_createDefCurveFromBase( pCurveCDS,pCurveIndex, vBase,pRet); 
}
STDMETHODIMP ARMCommonModule::ARM_Credit_DefaultProba(VARIANT *pCurve, 
												VARIANT *pMatu, 
												VARIANT *pRet)
{
	return theInstance.ARM_Credit_DefaultProba(pCurve, 
												pMatu, 
												pRet); 
	}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetBeta (VARIANT *pPricer, 
											VARIANT *pLabel,
											VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetBeta (pPricer, 
											pLabel,
											pRet); 
	}

STDMETHODIMP ARMCommonModule::ARM_Credit_Price (VARIANT *pPricer , 
										  VARIANT *pAsOfDate,
										  VARIANT *pRet)
{
	return theInstance.ARM_Credit_Price (pPricer , 
										  pAsOfDate,
										  pRet); 
	
}


STDMETHODIMP ARMCommonModule::ARM_Credit_Spread (VARIANT *pPricer , 
										  VARIANT *pMTM,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_Spread (pPricer , 
										  pMTM,
										  pRet); 
}



STDMETHODIMP ARMCommonModule::ARM_Credit_RiskyDuration (BSTR pDefCurve, 
												  VARIANT Date,
												  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_RiskyDuration (pDefCurve, 
												  Date,
												  pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetDefProbTranche (BSTR pPricer, 
													  double yearterm,
													  VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetDefProbTranche ( pPricer, 
													  yearterm,
													  pRet); 
	
}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetCleanSpreadTranche (VARIANT *pPricer, 
											VARIANT *pPlot,
											VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetCleanSpreadTranche (pPricer, 
											pPlot,
											pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_CDONPV (VARIANT *pPricer , 
										  VARIANT *pCPTTYPE,
										  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_CDONPV (pPricer , 
										  pCPTTYPE,
										  pRet) ;
}


STDMETHODIMP ARMCommonModule::ARM_Credit_CorrMatrix(VARIANT *pLabels,
											  VARIANT *pCoefs,
											  double AsOf,
											  BSTR	 Name,
											  VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_CorrMatrix(pLabels,
											  pCoefs,
											  AsOf,
											  	 Name,
											  pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_ExtractCorrMatrix(VARIANT *pCorrMatrixId,
													 VARIANT *pLabels,
													 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_ExtractCorrMatrix(pCorrMatrixId,
													 pLabels,
													 pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetLabel (VARIANT *pCurveId , 
										  	 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_GetLabel (pCurveId , 
										  	 pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_SetLabel (VARIANT *pCurveId , 
											 VARIANT *pLabel, 
										  	 VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_SetLabel (pCurveId , 
											 pLabel, 
										  	 pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_Sensitivity(VARIANT *pPricer , 
											   VARIANT *pType,
											   VARIANT *pPlot,
											   VARIANT *pLabel,
											   VARIANT *pEpsilon,
											   double epsilonGamma,
											   VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_Sensitivity(pPricer , 
											   pType,
											   pPlot,
											   pLabel,
											   pEpsilon,
											   epsilonGamma,
											   pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_GenSchedule(VARIANT *pAccStartDate,
											   VARIANT *pAccEndDate,
											   VARIANT *pFixingFreq,
											   VARIANT *pDayCountFrq,
											   VARIANT *prefdate,
											   VARIANT *pCurrency,
											   VARIANT *ptypeDates,
											   VARIANT *pModFol,
											   VARIANT *pGapCredit,
											   VARIANT *pRet)
{
	return theInstance.ARM_Credit_GenSchedule(pAccStartDate,
											   pAccEndDate,
											   pFixingFreq,
											   pDayCountFrq,
											   prefdate,
											   pCurrency,
											   ptypeDates,
											   pModFol,
											   pGapCredit,
											   pRet) ;
}


STDMETHODIMP ARMCommonModule::ARM_Credit_CashFlows(VARIANT *pCoefs, 
											  VARIANT *pRet)
{
	return theInstance.ARM_Credit_CashFlows(pCoefs, 
											  pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_Parameters(VARIANT *pCoefs, 
											  long	   nbcols,
											  VARIANT *pRet)
{
	return theInstance.ARM_Credit_Parameters(pCoefs, 
											  nbcols,
											  pRet); 

}
 


STDMETHODIMP ARMCommonModule::ARM_View (VARIANT *pObjet,
								  VARIANT *pRet)
{
	return theInstance.ARM_View (pObjet,
								  pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_Version(VARIANT *pRet)
{
	return theInstance.ARM_Credit_Version(pRet); 

}


/**
STDMETHODIMP ARMCommonModule::ARM_Credit_CorrelationSmile(VARIANT *pPricerId,
													double pSmileType,
													double pMktPrice,
													double pSeed,
													double pUpfrontPay,
													int datatype,
													VARIANT *pRet)	
{
	return theInstance.ARM_Credit_CorrelationSmile(pPricerId,
													pSmileType,
													pMktPrice,
													pSeed,
													pUpfrontPay,
													datatype,
													pRet)	 ;
	
}
**/ 

/*STDMETHODIMP ARMCommonModule::ARM_Credit_CpteDefLegTrancheEquivalente(VARIANT *pPricerId,
																double pStrikeUp1,
																double pStrikeUp2,
																double pBC1,
																double pBC2,
																double pSmileType,
																VARIANT *pRet)
{
	return theInstance.ARM_Credit_CpteDefLegTrancheEquivalente(pPricerId,
																pStrikeUp1,
																pStrikeUp2,
																pBC1,
																pBC2,
																pSmileType,
																pRet); 
	
}
*/
/*STDMETHODIMP ARMCommonModule::ARM_Credit_CpteStrikesEquivalents(VARIANT *pPricerId,
														  double pStrikeDown,
														  double pStrikeUp,
														  double pIndexNotional,
														  double pELTranche,
														  VARIANT *pRet)
{
	return theInstance.ARM_Credit_CpteStrikesEquivalents(pPricerId,
														  pStrikeDown,
														  pStrikeUp,
														  pIndexNotional,
														  pELTranche,
														  pRet); 

}
*/
STDMETHODIMP ARMCommonModule::ARM_Credit_CptBaseCorrelation(double AsOf,
															BSTR name,
															BSTR CalMethod,
															BSTR IndexId,
															VARIANT *pStrikeLow,
															VARIANT *pStrikeHigh,
															VARIANT *pVMktBid,
															VARIANT *pVMktAsk,
															VARIANT *pVUpfBid,
															VARIANT *pVUpfAsk,
															VARIANT *pVInitialCorrel,
															VARIANT *pVDeltaLevrage,
															BSTR ModelId,
															double integrationStep,
															double lagStartDate,
															double creditLag,
															VARIANT *pVectorPrevIndexId,
															VARIANT *pMatrixPrevBC,
															double step,
															BSTR CalMeth,
															VARIANT *pRet)
													  
{
	return theInstance.ARM_Credit_CptBaseCorrelation(   AsOf,
													    name,
														CalMethod,
														IndexId,
														pStrikeLow,
														pStrikeHigh,
														pVMktBid,
														pVMktAsk,
														pVUpfBid,
														pVUpfAsk,
														pVInitialCorrel,
														pVDeltaLevrage,
														ModelId,
														integrationStep,
														lagStartDate,
														creditLag,
														pVectorPrevIndexId,
														pMatrixPrevBC,
														step,
														CalMeth,
														pRet); 
}


/*STDMETHODIMP ARMCommonModule::ARM_Credit_GetImpliedCorrelation(VARIANT *pPricerId,
													  double pRange,
													  VARIANT *pMktSpreads,
													  double pSmileType,
													  VARIANT *pSeeds,
													  VARIANT *pUpfronts,
													  VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetImpliedCorrelation(pPricerId,
													  pRange,
													  pMktSpreads,
													  pSmileType,
													  pSeeds,
													  pUpfronts,
													  pRet); 
}
*/
STDMETHODIMP ARMCommonModule::ARM_Credit_GetDuration (VARIANT *pPricer,VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetDuration (pPricer,pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetCorrelationMatrix (BSTR pModelMultiCurves, 
														 BSTR pCorrMatrix, 
														 BSTR *pRet)
{
	return theInstance.ARM_Credit_SetCorrelationMatrix ( pModelMultiCurves, 
														  pCorrMatrix, 
														  pRet); 
	
}


STDMETHODIMP ARMCommonModule:: ARM_Credit_CloneCorrMatrixBary(BSTR pCorrMatrix, 
														double Beta,
														int UpOrDown,
														BSTR *pRet)
{
	return theInstance. ARM_Credit_CloneCorrMatrixBary(pCorrMatrix, 
														Beta,
														UpOrDown,
														pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_FwdSpreadPricer(VARIANT *pPricer , 
												   double Mty1,
												   double Mty2,
												   VARIANT *pRet)
{
	return theInstance.ARM_Credit_FwdSpreadPricer(pPricer , 
												   Mty1,
												   Mty2,
												   pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_ImpliedVol (VARIANT *pPricer , 
											   double pMktPrice,
											   VARIANT *pRet)
{
	return theInstance.ARM_Credit_ImpliedVol (pPricer , 
											   pMktPrice,
											   pRet); 
}



STDMETHODIMP ARMCommonModule::ARM_Credit_VirtualCdsSpread (VARIANT *pPricer,
													 VARIANT *pMaturity,
													 VARIANT *pRet)
{
	return theInstance.ARM_Credit_VirtualCdsSpread (pPricer,
													 pMaturity,
													 pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_BSGreeks(VARIANT *pPricer,
											VARIANT *pGreekType,
											VARIANT *pRet)
{
	return theInstance.ARM_Credit_BSGreeks(pPricer,
											pGreekType,
											pRet); 

}



STDMETHODIMP ARMCommonModule::ARM_Credit_SetCorrelation (BSTR pModelMultiCurvesId , 
												   BSTR pCorrelationId,
												   VARIANT *pRet)
{
	return theInstance.ARM_Credit_SetCorrelation (pModelMultiCurvesId , 
												   pCorrelationId,
												   pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_CorrelationStrike(VARIANT* pLabels, 
													 VARIANT* pVolCurves, 
													 VARIANT* pProportions,
													 VARIANT* pSmileStrikeLow,
													 VARIANT* pSmileStrikeHigh,
													 VARIANT* pVIndex,
													 double AsOf,
													 BSTR Name,
													 BSTR *pRet)
{
	return theInstance.ARM_Credit_CorrelationStrike( pLabels, 
													 pVolCurves, 
													 pProportions,
													 pSmileStrikeLow,
													 pSmileStrikeHigh,
													 pVIndex,
													 AsOf,
													 Name,
													 pRet); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_CorrelationSmileStrike(VARIANT* pLabels, 
															  VARIANT* pVolCurves, 
															  VARIANT* pProportions,
															  double AsOf,
															  VARIANT* pSmileStrikeLow,
													          VARIANT* pSmileStrikeHigh,
													          VARIANT* pVIndex,
														      BSTR Name,
															  VARIANT*  pFullStrikeLow,
															  VARIANT*  pFullStrikeUp,
															  BSTR *pRet) 
{
	return theInstance.ARM_Credit_CorrelationSmileStrike( pLabels, 
													 pVolCurves, 
													 pProportions,
													 AsOf,
													 pSmileStrikeLow,
													 pSmileStrikeHigh,
													 pVIndex,
													 Name,
													 pFullStrikeLow,
													 pFullStrikeUp,
													 pRet); 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_Beta_Correlation(VARIANT *pLabels,
													VARIANT *pCoefs, 
													double AsOf,
													BSTR Name,
													BSTR idIndex1,
													BSTR idIndex2,
													VARIANT *pRet)
{
	return theInstance.ARM_Credit_Beta_Correlation(pLabels,
													pCoefs, 
													AsOf,
													Name,
													idIndex1,
													idIndex2,
													pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetEqStrikeDown(BSTR correlId,
												   BSTR indexname, 
													double *pRet)
{
	return theInstance.ARM_Credit_GetEqStrikeDown(correlId,
												   indexname, 
													pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_GetEqStrikeUp(BSTR correlId,
												 BSTR indexname, 
												 double *pRet)
{
	return theInstance.ARM_Credit_GetEqStrikeUp(correlId,
												 indexname, 
												 pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_GetCorrelStrikeDown(BSTR correlId,
													   double yfmaturity,
													   double *pRet)
{
	return theInstance.ARM_Credit_GetCorrelStrikeDown(correlId,
													   yfmaturity,
													   pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetCorrelStrikeUp(BSTR correlId,
													 double yfmaturity,
													 double *pRet)
{
	return theInstance.ARM_Credit_GetCorrelStrikeUp(correlId,
													 yfmaturity,
													 pRet); 


}


STDMETHODIMP ARMCommonModule::ARM_Credit_GetCorrelation(BSTR ModelId,
												  BSTR *pRet)
{
	return theInstance.ARM_Credit_GetCorrelation(ModelId,
												  pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_NextCpnDate(double AsOfDate,
										double maturity,
										BSTR frequency,
										BSTR rule,
										BSTR currency,
										BSTR intrule,
 									    double *pRet)
{
	return theInstance.ARM_NextCpnDate(AsOfDate,
										maturity,
										frequency,
										rule,
										currency,
										intrule,
 									    pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_SetProportionsInfos(BSTR correlId,
													   BSTR IndexName,
													   double proportion,
													   double forcedstrikelow,
													   double forcedstrikehigh)
{
	return theInstance.ARM_Credit_SetProportionsInfos(correlId,
													   IndexName,
													   proportion,
													   forcedstrikelow,
													   forcedstrikehigh); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_CptImplCvForCDO2(BSTR pricerId,
													BSTR Name,
													BSTR Tenor,
													double *pRet)
{
	return theInstance.ARM_Credit_CptImplCvForCDO2(pricerId,
													Name,
													Tenor,
													pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_Credit_AddPeriod(double pAsOf,
										BSTR Maturity,
										BSTR pCcy,
										BSTR AdjRule,
										BSTR AdjCDS,
 									    VARIANT *pRet)
{
	return theInstance.ARM_Credit_AddPeriod(pAsOf,
										Maturity,
										pCcy,
										AdjRule,
										AdjCDS,
 									    pRet); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_SetCoupons (BSTR CdsorCdoId, 
											   BSTR CouponsId, 
											   BSTR TypesId, 
											   BSTR PartId, 
											   BSTR *pRet)
{
	return theInstance.ARM_Credit_SetCoupons (CdsorCdoId, 
											   CouponsId, 
											   TypesId, 
											   PartId, 
											   pRet) ;
	
}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetRiskyProfile (BSTR CdsorCdoId, 
											   BSTR CouponsId, 
											   BSTR TypesId, 
											   BSTR *pRet)
{
	return theInstance.ARM_Credit_SetRiskyProfile (CdsorCdoId, 
											   CouponsId, 
											   TypesId, 
											   pRet); 
	
}


STDMETHODIMP ARMCommonModule::ARM_Credit_DataFromLabel (BSTR	 pPricer,
												  BSTR	 pLabel,
												  double *pRet)
{
	return theInstance.ARM_Credit_DataFromLabel (	 pPricer,
												  	 pLabel,
												  pRet); 
	
}



STDMETHODIMP ARMCommonModule::ARM_Credit_GetEqStrike(BSTR correlId,
												BSTR indexname, 
												BSTR UpOrLow, 
												VARIANT *pRetMatu,
												VARIANT *pRetStrikes)
{
	return theInstance.ARM_Credit_GetEqStrike(correlId,
												indexname, 
												UpOrLow, 
												pRetMatu,
												pRetStrikes); 

}


STDMETHODIMP ARMCommonModule::ARM_Credit_DefaultIntensity(BSTR pricerId,
													double Maturity,
													double *pRet)
{
return theInstance.ARM_Credit_DefaultIntensity(pricerId,
													Maturity,
													pRet);
}


STDMETHODIMP ARMCommonModule::ARM_Credit_SetMatuLabel (BSTR pCurveId , 
											 VARIANT *pMatuLabels, 
										  	 double *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_SetMatuLabel (pCurveId , 
											 pMatuLabels, 
										  	 pRet) ;

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetPricerForRatesComputation(BSTR legId,
															  BSTR pricerId,
															  BSTR *pRet)
{
	return theInstance.ARM_Credit_SetPricerForRatesComputation(legId,
															  pricerId,
															  pRet) ;

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetRecovCoef (BSTR pSecId , 
											 double RecovCoef, 
										  	 double *pRet)
{
	return theInstance.ARM_Credit_SetRecovCoef (pSecId , 
											 RecovCoef, 
										  	 pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetFees(BSTR securityId,
										   BSTR RefvalueId)
{
	return theInstance.ARM_Credit_SetFees(securityId,
										   RefvalueId); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetInterpolationType(VARIANT *pVolCurve,
														VARIANT *pInterpolType,
														VARIANT *pRet)
{
	return theInstance.ARM_Credit_SetInterpolationType(pVolCurve,
														pInterpolType,
														pRet); 

}




STDMETHODIMP ARMCommonModule::ARM_Credit_DefProbModelNew(BSTR pDefCurve,
												   BSTR pIRcurve,
												   BSTR pVolcurve,
												   VARIANT *pRet)
{
	return  theInstance.ARM_Credit_DefProbModelNew(pDefCurve,
												   pIRcurve,
												   pVolcurve,
												   pRet) ; 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_ModelMultiCurves(VARIANT *pIRcurve,
													VARIANT *pDefCurves, 
													VARIANT *pRecoveryRates, 
													BSTR	CorrelId,
													BSTR	pVolcurve,
													BSTR	pCpnInfcurve,
													BSTR	pCpnIRcurve,
													VARIANT *pRet)
{
	return theInstance.ARM_Credit_ModelMultiCurves(pIRcurve,
													pDefCurves, 
													pRecoveryRates, 
														CorrelId,
														pVolcurve,
														pCpnInfcurve,
														pCpnIRcurve,
													pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_ModelMultiCvMktDataMng(VARIANT *pIRcurve,
													VARIANT *pDefCurves, 
													VARIANT *pRecoveryRates, 
													BSTR	CorrelId,
													BSTR	MktDataMngId,
													BSTR	pVolcurve,
													BSTR	cloneOrNot,
													VARIANT *pRet)
{
	return theInstance.ARM_Credit_ModelMultiCvMktDataMng(pIRcurve,
														pDefCurves, 
														pRecoveryRates, 
														CorrelId,
														MktDataMngId,
														pVolcurve,
														cloneOrNot,
														pRet); 

}



STDMETHODIMP ARMCommonModule::ARM_Credit_Pricer(BSTR pSecurity, 
										  BSTR pModel,
										  BSTR pPricerType,
										  int		l_nbpaths,
										  BSTR     pParameters,
										  double   AsOfDate,
										  VARIANT *pRet)
{
	return theInstance.ARM_Credit_Pricer(pSecurity, 
										  pModel,
										  pPricerType,
										  l_nbpaths,
										      pParameters,
										    AsOfDate,
										  pRet) ; 


}


/**
STDMETHODIMP ARMCommonModule::ARM_Credit_GetModelFromSummit(BSTR IRcurve,
													  BSTR IDSummit, 
													  BSTR type,	
													  VARIANT *pRet)
{
	return theInstance.ARM_Credit_GetModelFromSummit(IRcurve,
													  IDSummit, 
													  type,	
													  pRet); 

}

**/ 
STDMETHODIMP ARMCommonModule::ARM_Credit_SetVolatility(VARIANT *pPricer , 
												 BSTR pVolcurve,
												 VARIANT *pRet)
{
	return theInstance.ARM_Credit_SetVolatility(pPricer , 
												 pVolcurve,
												 pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_SetVolCurve(BSTR Model , 
												 BSTR pVolcurve,
												 VARIANT *pRet)
{
	return theInstance.ARM_Credit_SetVolCurve(Model , 
												 pVolcurve,
												 pRet); 

}

STDMETHODIMP ARMCommonModule::ARM_Credit_GenerateImpliedCurve (BSTR pricerId, 
														 BSTR *pRet)
{
	return theInstance.ARM_Credit_GenerateImpliedCurve (pricerId, 
														 pRet); 

}




STDMETHODIMP ARMCommonModule::ARMSetEtoolkit(BSTR pUserName,
									   BSTR pPassWord,
									   BSTR pDatabaseContext,
									   BSTR pItConfigDomainDir,
									   BSTR pItDomainName,
									   long *pRet)
{
	return theInstance.ARMSetEtoolkit(pUserName,
									   pPassWord,
									   pDatabaseContext,
									   pItConfigDomainDir,
									   pItDomainName,
									    pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMConnectionEtoolkit(VARIANT *pRet)
{
	return theInstance.ARMConnectionEtoolkit(pRet); 

}

STDMETHODIMP ARMCommonModule::ARMDeconnectionEtoolkit(VARIANT *pRet)
{
	return theInstance.ARMDeconnectionEtoolkit(pRet); 

}






STDMETHODIMP ARMCommonModule::ARMcptBonibor(double pDate,
									  BSTR pDepart,
									  BSTR pMatStruct,
									  BSTR pMatTot,
									  BSTR pAmort,
									  BSTR pFreq,
									  BSTR pSjUSD,
									  BSTR pTiming,
									  double pBarriere,
									  double pSpdPostBar,
									  double pMarge,
									  double pFunding,
									  BSTR pFundingFreq,
									  double pSpd2phase,
									  double pSoulte,
									  BSTR pYcModId,
									  BSTR pBsModId,
									  BSTR pBsModVolUSDId,
									  BSTR pBsModCorrPlusId,
									  BSTR pBsModCorrMoinsId,
									  BSTR pCrossModId,
									  BSTR pProbaMarge,
									  double pInt,
									  VARIANT *pRet)
{
	return theInstance.ARMcptBonibor(pDate,
									  pDepart,
									  pMatStruct,
									  pMatTot,
									  pAmort,
									  pFreq,
									  pSjUSD,
									  pTiming,
									  pBarriere,
									  pSpdPostBar,
									  pMarge,
									  pFunding,
									  pFundingFreq,
									  pSpd2phase,
									  pSoulte,
									  pYcModId,
									  pBsModId,
									  pBsModVolUSDId,
									  pBsModCorrPlusId,
									  pBsModCorrMoinsId,
									  pCrossModId,
									  pProbaMarge,
									  pInt,
									  pRet); 
}




STDMETHODIMP ARMCommonModule::ARMcptDigital(double pAsOf,
								   double pStartDate,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pBsmod,
								   BSTR pBsmodDelta,
								   BSTR pBsmodVega,
								   BSTR pFreqP,
								   BSTR pResetTiming,
								   BSTR pPorR,
								   BSTR pCcy,
								   BSTR pCcyIdx,
								   BSTR pDayCount,
								   BSTR pCapOrFloor,
								   BSTR pAmort,
								   BSTR pStrike,
								   BSTR pPayOff,
								   BSTR pSpd,
								   double pResetGap,
								   double pSpreadBelow,
								   double pSpreadAbove,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,
								   double pRefDate,
								   VARIANT *pRet)

{
	return theInstance.ARMcptDigital(pAsOf,
								   pStartDate,
								   pEndDate,
								   pNtl,
								   pIndex,
								   pBsmod,
								   pBsmodDelta,
								   pBsmodVega,
								   pFreqP,
								   pResetTiming,
								   pPorR,
								   pCcy,
								   pCcyIdx,
								   pDayCount,
								   pCapOrFloor,
								   pAmort,
								   pStrike,
								   pPayOff,
								   pSpd,
								   pResetGap,
								   pSpreadBelow,
								   pSpreadAbove,
								   pFwdRule,
								   pIntRule,
								   pStubRule,
								   pFreqAmort,
								   pTxAmort,
								   pAmountAmort,
								   pRefDate,
								   pRet); 
}



STDMETHODIMP ARMCommonModule::ARMcptSPTQTF(double pAsOf,
										   double pStartDate,
										   double pStartDatePhase2,
										   double pStartDatePhase3, 
										   double pEndDate, 
										   double pNtl, 
										   BSTR pIndexPhase1, 
										   BSTR pIndexPhase2, 
										   BSTR pIndexFund,
										   BSTR pIndexPhase3,  
										   BSTR pFreqPPhase1,
										   BSTR pFreqPPhase2, 
										   BSTR pFreqPFund, 
										   BSTR pFreqPPhase3, 
										   BSTR pFreqR, 
										   BSTR pResetTimingPhase1,
										   BSTR pResetTimingPhase2,
										   BSTR pResetTimingPhase3, 
										   BSTR pCcy, 
										   BSTR pCcyIdx, 
										   BSTR pDayCount,
										   double pFee, 
										   BSTR pIsRateFixedPhase2, 
										   double pFixedRatePhase2,
										   BSTR pBarrier, 
										   BSTR pSpdPhase1, 
										   BSTR  pSpdPhase1fund, 
										   BSTR pSpdPhase2Tf,
										   BSTR pSpdPhase2fund,
										   BSTR pSpdPhase3, 
										   BSTR pSpdPhase3fund, 
										   double pResetGapPhase1, 
										   double pResetGapPhase2, 
										   double pResetGapPhase3,
										   BSTR pAmort,
										   BSTR pBsmod,
										   BSTR pBsmodDeltaCcy1,
										   BSTR pBsmodVegaCcy1,
										   BSTR pBsmodDeltaCcy2,
										   BSTR pBsmodVegaCcy2,
										   BSTR pBsmodFxCorrel,
										   BSTR pFwdRule, 
										   BSTR pIntRule,
										   BSTR pStubRule,
										   BSTR pFreqAmort, 
										   double pTxAmort, 
										   double pAmountAmort, 
										   double pRefDate, 
										   VARIANT *pRet)

{
	return theInstance.ARMcptSPTQTF(pAsOf,
									pStartDate,
									pStartDatePhase2,
									pStartDatePhase3, 
									pEndDate, 
									pNtl, 
									pIndexPhase1, 
									pIndexPhase2, 
									pIndexFund,
									pIndexPhase3,  
									pFreqPPhase1,
									pFreqPPhase2, 
									pFreqPFund, 
									pFreqPPhase3, 
									pFreqR, 
									pResetTimingPhase1,
									pResetTimingPhase2,
									pResetTimingPhase3, 
									pCcy, 
									pCcyIdx, 
									pDayCount,
									pFee, 
									pIsRateFixedPhase2, 
									pFixedRatePhase2,
									pBarrier, 
									pSpdPhase1, 
									pSpdPhase1fund, 
									pSpdPhase2Tf,
									pSpdPhase2fund,
									pSpdPhase3, 
									pSpdPhase3fund, 
									pResetGapPhase1, 
									pResetGapPhase2, 
									pResetGapPhase3,
									pAmort,
									pBsmod,
									pBsmodDeltaCcy1,
									pBsmodVegaCcy1,
									pBsmodDeltaCcy2,
									pBsmodVegaCcy2,
									pBsmodFxCorrel,
									pFwdRule, 
									pIntRule,
									pStubRule,
									pFreqAmort, 
									pTxAmort, 
									pAmountAmort, 
									pRefDate, 
									pRet);
}




//Ajout TD le 04/07/2005
STDMETHODIMP ARMCommonModule::ARMcomputeBilibor(double pAsOf,
								   double pStartDate,
								   double pDateSecondPhase,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pIndexFund,
								   BSTR pIndexSF,
								   BSTR pBsmod,
								   BSTR pBsmodFund,
								   BSTR pBsmodDeltaCcy1,
								   BSTR pBsmodDeltaFund,
								   BSTR pBsmodDeltaCcy2,
								   BSTR pBsmodFxCorrel,
								   BSTR pFreqP,
								   BSTR pFreqR,
								   BSTR pFreqPFund,
								   BSTR pFreqRFund, 
								   BSTR pFreqPSF,
								   BSTR pFreqRSF,
								   BSTR pResetTiming,
								   BSTR pResetTimingSF,
								   BSTR pCcy1,
								   BSTR pCcy2,
								   BSTR pDayCount,
								   BSTR pDayCountSF,
								   double pSpdPF,
								   double pSpdSF,
								   double pSpdfund,
								   double pSpdfund2, 
								   double pResetGap,
								   double pResetGapSF,
								   BSTR pAmort,
								   double pRefDate,
								   double pFee,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,						   
								   VARIANT *pRet)

{
	return theInstance.ARMcomputeBilibor(pAsOf,
								   pStartDate,
								   pDateSecondPhase,
								   pEndDate,
								   pNtl,
								   pIndex,
								   pIndexFund,
								   pIndexSF,
								   pBsmod,
								   pBsmodFund,
								   pBsmodDeltaCcy1,
								   pBsmodDeltaFund,
								   pBsmodDeltaCcy2,
								   pBsmodFxCorrel,
								   pFreqP,
								   pFreqR,
								   pFreqPFund,
								   pFreqRFund, 
								   pFreqPSF,
								   pFreqRSF,
								   pResetTiming,
								   pResetTimingSF,
								   pCcy1,
								   pCcy2,
								   pDayCount,
								   pDayCountSF,
								   pSpdPF,
								   pSpdSF,
								   pSpdfund,
								   pSpdfund2, 
								   pResetGap,
								   pResetGapSF,
								   pAmort,
								   pRefDate,
								   pFee,
								   pFwdRule,
								   pIntRule,
								   pStubRule,
								   pFreqAmort,
								   pTxAmort,
								   pAmountAmort,								   
								   pRet) ; 

	}




//Ajout TD le 16/08/2005
STDMETHODIMP ARMCommonModule::ARMcomputeOptilix(double pAsOf,
								   double pStartDate,
								   double pDateSecondPhase,
								   double pEndDate,
								   double pNtl,
								   BSTR pIndex,
								   BSTR pIndexFund,
								   BSTR pIndexSF,
								   BSTR pBsmod,
								   BSTR pBsmodFund,
								   BSTR pBsmodDeltaCcy,
								   BSTR pBsmodDeltaFund,
								   BSTR pFreqP,
								   BSTR pFreqR,
								   BSTR pFreqPFund,
								   BSTR pFreqRFund, 
								   BSTR pFreqPSF,
								   BSTR pFreqRSF,
								   BSTR pResetTiming,
								   BSTR pResetTimingSF,
								   BSTR pCcy,
								   BSTR pDayCount,
								   BSTR pDayCountSF,
								   double pSpdSF,
								   VARIANT pSpdfund,
								   double pResetGap,
								   double pResetGapSF,
								   BSTR pAmort,
								   double pRefDate,
								   double pFee,
								   BSTR pFwdRule,
								   BSTR pIntRule,
								   BSTR pStubRule,
								   BSTR pFreqAmort,
								   double pTxAmort,
								   double pAmountAmort,
								   VARIANT *pRet)

{
	return theInstance.ARMcomputeOptilix(pAsOf,
								   pStartDate,
								   pDateSecondPhase,
								   pEndDate,
								   pNtl,
								   pIndex,
								   pIndexFund,
								   pIndexSF,
								   pBsmod,
								   pBsmodFund,
								   pBsmodDeltaCcy,
								   pBsmodDeltaFund,
								   pFreqP,
								   pFreqR,
								   pFreqPFund,
								   pFreqRFund, 
								   pFreqPSF,
								   pFreqRSF,
								   pResetTiming,
								   pResetTimingSF,
								   pCcy,
								   pDayCount,
								   pDayCountSF,
								   pSpdSF,
								   pSpdfund,
								   pResetGap,
								   pResetGapSF,
								   pAmort,
								   pRefDate,
								   pFee,
								   pFwdRule,
								   pIntRule,
								   pStubRule,
								   pFreqAmort,
								   pTxAmort,
								   pAmountAmort,
								   pRet) ; 
}




STDMETHODIMP ARMCommonModule::ARMcomputePentifix(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund,  
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pFloorOrCap, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub, 
										   BSTR pCorrManager, 
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   VARIANT *pBumpBsGenMod,
										   VARIANT *pBumpVolBsGenMod,
										   VARIANT *pRet
										   )

{	
	return theInstance.ARMcomputePentifix(pNtl, 
										   pStartDatePhase1,  
										   pCcy,
										   pIndexPhase1, 										   
										   pSpreadPhase1,
										   pDayCountPhase1,
										   pPayFreqPhase1,
										   pResetFreqPhase1,
										   pResetTimingPhase1,
										   pRoll, 
										   pAdjPhase1, 
										   pStub,
										   pIndexPhase2DIG,
										   pIndexLongPhase2DIG,
										   pStrikePhase2DIG,
										   pResetTimingPhase2DIG,
										   pAdjPhase2DIG, 
										   pStartDatePhase2,
										   pSpreadPhase2,										   
										   pDayCountPhase2, 
										   pPayFreqPhase2, 
										   pResetFreqPhase2,
										   pAdjPhase2,
										   pStartDatePhase3,
										   pEndDatePhase3,
										   pIndexPhase3,
										   pSpreadPhase3, 
										   pDayCountPhase3, 
										   pPayFreqPhase3, 
										   pResetFreqPhase3,
										   pResetTimingPhase3,
										   pAdjPhase3,
										   pIndexFund,
										    pSpreadFund,
										   pDayCountFund, 
										   pPayFreqFund, 
										   pResetFreqFund,
										   pResetTimingFund,
										   pAdjFund,  
										   pEndDateAmort,
										   pDayCountAmort,
										   pIntRuleAmort, 
										   pTxAmort,
										   pFreqAmort,
										   pAmountAmort, 
										   pTypeAmort, 
										   pFloorOrCap, 
										   pFee, 
										   pVolCurvFromMatriceShift,
										   pVol, 
										   pVolCub, 
										   pCorrManager, 
										   pConvexityManager, 
										   pZc, 
										   pSmiledMod, 
										   pSmiledModBump, 					
										   pHyperCubeCorrel,
										   pBumpBsGenMod,
										   pBumpVolBsGenMod, 
										   pRet
										   ) ; 

}



STDMETHODIMP ARMCommonModule::ARMcomputePentibor(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPay, 					
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub,  
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   BSTR pIndexIndexCorrelCube, 
										   BSTR pCorrEUR, 
										   BSTR pInterCorr, 
										   VARIANT *pRet
										   )

{	
	return theInstance.ARMcomputePentibor(pNtl, 
										   pStartDatePhase1,  
										   pCcy,
										   pIndexPay, 					
										   pIndexPhase1, 										   
										   pSpreadPhase1,
										   pDayCountPhase1,
										   pPayFreqPhase1,
										   pResetFreqPhase1,
										   pResetTimingPhase1,
										   pRoll, 
										   pAdjPhase1, 
										   pStub,
										   pIndexPhase2DIG,
										   pIndexLongPhase2DIG,
										   pStrikePhase2DIG,
										   pResetTimingPhase2DIG,
										   pAdjPhase2DIG, 
										   pStartDatePhase2,
										   pSpreadPhase2,										   
										   pDayCountPhase2, 
										   pPayFreqPhase2, 
										   pResetFreqPhase2,
										   pAdjPhase2,
										   pStartDatePhase3,
										   pEndDatePhase3,
										   pIndexPhase3,
										   pSpreadPhase3, 
										   pDayCountPhase3, 
										   pPayFreqPhase3, 
										   pResetFreqPhase3,
										   pResetTimingPhase3,
										   pAdjPhase3,
										   pIndexFund,
										    pSpreadFund,
										   pDayCountFund, 
										   pPayFreqFund, 
										   pResetFreqFund,
										   pResetTimingFund,
										   pAdjFund, 
										   pEndDateAmort,
										   pDayCountAmort,
										   pIntRuleAmort, 
										   pTxAmort,
										   pFreqAmort,
										   pAmountAmort, 
										   pTypeAmort, 
										   pFee, 
										   pVolCurvFromMatriceShift,
										   pVol, 
										   pVolCub,  
										   pConvexityManager, 
										   pZc, 
										   pSmiledMod, 
										   pSmiledModBump, 					
										   pHyperCubeCorrel,
										   pIndexIndexCorrelCube, 
										   pCorrEUR, 
										   pInterCorr, 
										   pRet
										   ) ;
}



STDMETHODIMP ARMCommonModule::ARMNextBusinessDay(double pDate, BSTR pCalendrier, long pNbDays, double *pDate2)
{
	return theInstance.ARMNextBusinessDay(pDate, pCalendrier, pNbDays, pDate2); 

}

STDMETHODIMP ARMCommonModule::ARMAdjustToBusDate(double pDate, BSTR pCalendrier, BSTR pRule, double *pDate2)
{
	return theInstance.ARMAdjustToBusDate(pDate, pCalendrier, pRule, pDate2) ;

}


STDMETHODIMP ARMCommonModule::ARMFreeObject(BSTR pId, long *pRet)
{
	return theInstance.ARMFreeObject(pId, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMIsBusinessDay(double pDate, BSTR pCalendrier, long *pRes)
{
	return theInstance.ARMIsBusinessDay(pDate, pCalendrier, pRes) ; 

}


STDMETHODIMP ARMCommonModule::ARMGetZCFromSummit(BSTR pIndex, 
												 BSTR pCurrency, 
												 BSTR pCvName, 
												 double pDate, 
												 BSTR pInterpMethod, 
												 BSTR *pRet)
{
	return theInstance.ARMGetZCFromSummit(pIndex, pCurrency, pCvName, pDate, pInterpMethod, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMCreateZCFromSummit(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate,BSTR pAdj,BSTR pRaw, BSTR *pRet)
{
	return theInstance.ARMCreateZCFromSummit( pIndex,  pCurrency,  pCvName, pDate, pAdj, pRaw,  pRet) ; 
	
}


STDMETHODIMP ARMCommonModule::ARMBumpCurve(BSTR pZc, double pEpsilon, long pMethod, BSTR pPLot, BSTR *pRet)
{
	return theInstance.ARMBumpCurve( pZc, pEpsilon, pMethod, pPLot, pRet); 

}


STDMETHODIMP ARMCommonModule::ARMFreeAllObjects(long *pRet)
{
	return theInstance.ARMFreeAllObjects(pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMYcMod(BSTR pZc, BSTR pZcDiscount, BSTR *pRet)
{
	return theInstance.ARMYcMod( pZc,  pZcDiscount,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMForwardYield(BSTR pZc, double pMatu1, double pMatu2, BSTR pMeth, BSTR pAdj, double *pRet)
{
	return theInstance.ARMForwardYield( pZc, pMatu1, pMatu2,  pMeth,  pAdj, pRet); 

}


STDMETHODIMP ARMCommonModule::ARMDiscountYield(VARIANT *pZc, VARIANT *pMatu, VARIANT *pMeth, VARIANT *pRet)
{
	return theInstance.ARMDiscountYield(pZc, pMatu, pMeth, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMLiborSwap(VARIANT *pStartDate, VARIANT *pEndDate, VARIANT *pLiborType, VARIANT *pRecOrPay, VARIANT *pFixedRate, VARIANT *pSpread, VARIANT *pCcy, BSTR pDaycount, BSTR pFloatingDaycount, VARIANT *pRet)
{
	return theInstance.ARMLiborSwap(pStartDate, pEndDate, pLiborType, pRecOrPay, pFixedRate, pSpread, pCcy, pDaycount, pFloatingDaycount, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMSwapPriceToRate(VARIANT *pSwap, VARIANT *pDate, VARIANT *pPrice, VARIANT *pModel, VARIANT *pRet)
{
	return theInstance.ARMSwapPriceToRate(pSwap, pDate, pPrice, pModel, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMPrice(VARIANT *pSec, VARIANT *pModel, VARIANT *pRet)
{
	return theInstance.ARMPrice(pSec, pModel, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMBetweenDates(VARIANT *pDate1, VARIANT *pDate2, VARIANT *pDaycount, VARIANT *pIsYearFrac, VARIANT *pRet)
{
	return theInstance.ARMBetweenDates(pDate1, pDate2, pDaycount, pIsYearFrac, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMAddPeriod(VARIANT *pDate, VARIANT *pFreq, VARIANT *pCcy, VARIANT *pNbPeriod, VARIANT *pAdjRule, VARIANT *pRet)
{
	return theInstance.ARMAddPeriod(pDate, pFreq, pCcy, pNbPeriod, pAdjRule, pRet) ; 
	
}


STDMETHODIMP ARMCommonModule::ARMIsoCcy(BSTR pCcy, BSTR pRefObj, BSTR *pRet)
{
	return theInstance.ARMIsoCcy( pCcy,  pRefObj,  pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMGetSpotDays(VARIANT *pCcy, VARIANT *pRet)
{
	return theInstance.ARMGetSpotDays(pCcy, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMGetLiborIndexDaycount(VARIANT *pCcy, VARIANT *pRet)
{
	return theInstance.ARMGetLiborIndexDaycount(pCcy, pRet) ; 

}

STDMETHODIMP ARMCommonModule::ARMGetLiborTerm(VARIANT *pCcy, VARIANT *pRet)
{
	return theInstance.ARMGetLiborTerm(pCcy, pRet) ; 

}

STDMETHODIMP ARMCommonModule::ARMGetFixedDayCount(VARIANT *pCcy, VARIANT *pRet)
{
	return theInstance.ARMGetFixedDayCount(pCcy, pRet) ; 

}

STDMETHODIMP ARMCommonModule::ARMGetFixedPayFreq(VARIANT *pCcy, VARIANT *pRet)
{
	return theInstance.ARMGetFixedPayFreq(pCcy, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMComputeVolatility(VARIANT *pVol, VARIANT* pMatu, VARIANT* pStrike, VARIANT* pTenor, VARIANT *pRet)
{
	return theInstance.ARMComputeVolatility(pVol,  pMatu,  pStrike,  pTenor, pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMVolCurv(VARIANT *pMatu, VARIANT* pStrikes, VARIANT* pVols, double pAsOf, BSTR pStkType, BSTR pVolType, BSTR pCcy, BSTR pIndexId,BSTR *pRet)
{
	return theInstance.ARMVolCurv(pMatu,  pStrikes,  pVols, pAsOf,  pStkType,  pVolType,  pCcy,  pIndexId,pRet) ;

}



STDMETHODIMP ARMCommonModule::ARMGetVolCubeFromSummit(BSTR pIndex, BSTR pCcy, BSTR pCvName, double pAsOf, BSTR pType, VARIANT* pSmiles, BSTR pTypeCube, BSTR indexId,BSTR *pRet)
{
	return theInstance.ARMGetVolCubeFromSummit(pIndex,  pCcy,  pCvName,  pAsOf,  pType,  pSmiles,  pTypeCube, indexId,pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMParallelShift(BSTR pZc, double pBump, BSTR *pRet)
{
	return theInstance.ARMParallelShift( pZc, pBump,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMBumpVolatility(BSTR pVol, double pValue, long pNthLine, long pNthCol, BSTR pIsCumul, BSTR *pRet)
{
	return theInstance.ARMBumpVolatility( pVol, pValue, pNthLine, pNthCol,  pIsCumul,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMGlobDFBS(BSTR pDomBSId, BSTR pDomCurrId, BSTR pFrgBSId, BSTR pFrgCurrId, BSTR pFxVolCrvId, BSTR pFFxCorrId, BSTR pRatesCorrId, BSTR pFxVolModelId, BSTR *pRet)
{
	return theInstance.ARMGlobDFBS( pDomBSId,  pDomCurrId,  pFrgBSId,  pFrgCurrId,  pFxVolCrvId,  pFFxCorrId,  pRatesCorrId,  pFxVolModelId, pRet) ;
}

STDMETHODIMP ARMCommonModule::ARMDFFXBS(BSTR pDVolId,BSTR pFVolId,BSTR pDZcId,BSTR pFZcId,BSTR pDFxCorrId,BSTR pFFxCorrId,BSTR pFxVolId,double pRatesCorr,BSTR *pRet)
{
return theInstance.ARMDFFXBS( pDVolId, pFVolId, pDZcId, pFZcId, pDFxCorrId, pFFxCorrId, pFxVolId,pRatesCorr, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMTRIBSMODEL(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,BSTR *pRet)
{
	return theInstance.ARMTRIBSMODEL( pModel1, pModel2, pDiscModel, pFX1DiscVol, pFX2DiscVol, pIdx1Idx2Corr, pIdx1DiscIdxCorr, pIdx2DiscIdxCorr, pIdx1FxCorr, pIdx2FxCorr,pQuantoFlag, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMTRIBSDUAL(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,double pCorrelForAdj,int pWithslopeflag,BSTR *pRet)
{
return theInstance.ARMTRIBSDUAL( pModel1, pModel2, pDiscModel, pFX1DiscVol, pFX2DiscVol, pIdx1Idx2Corr, pIdx1DiscIdxCorr, pIdx2DiscIdxCorr, pIdx1FxCorr, pIdx2FxCorr,pQuantoFlag,pCorrelForAdj,pWithslopeflag, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMBsSmiledModel(double pDate, double pSpot, BSTR pDividend, BSTR pDiscrate, BSTR pVolATM, BSTR pRo, BSTR pNu, BSTR pIsSABR, BSTR pBeta, BSTR *pRet)
{
	return theInstance.ARMBsSmiledModel(pDate, pSpot,  pDividend,  pDiscrate,  pVolATM,  pRo,  pNu,  pIsSABR,  pBeta,  pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMGetVolFromSummit(VARIANT *pIndex, VARIANT* pCcy, VARIANT* pCvName, VARIANT* pAsOf, VARIANT* pType, VARIANT* pMatIndex, VARIANT* pImpOrHist, BSTR indexId,VARIANT *pRet)
{
	return theInstance.ARMGetVolFromSummit(pIndex,  pCcy,  pCvName,  pAsOf,  pType,  pMatIndex,  pImpOrHist, indexId,pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARMGetFXVolFromSummit(BSTR pCcy1, BSTR pCcy2, double pDate, BSTR pCvName, BSTR pType, BSTR *pRet)
{
	return theInstance.ARMGetFXVolFromSummit( pCcy1,  pCcy2, pDate,  pCvName,  pType,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMGetFXCorrelFromSummit(BSTR pCcy1, BSTR pIndex, BSTR pCcy2, double pDate, BSTR pCvName, VARIANT* pTenors, BSTR *pRet)
{
	return theInstance.ARMGetFXCorrelFromSummit( pCcy1,  pIndex,  pCcy2, pDate,  pCvName,  pTenors,  pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARMVolFlat(double pVol, double pDate, BSTR pCcy, BSTR* pRet)
{
	return theInstance.ARMVolFlat(pVol, pDate,  pCcy,  pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMVolCube(BSTR pATMVol,VARIANT *pSmileCurveIds,VARIANT *pTenors,BSTR pVolType,BSTR pRefObj, BSTR *pRet)
{
	return theInstance.ARMVolCube( pATMVol,pSmileCurveIds,pTenors, pVolType, pRefObj, pRet) ;
}


STDMETHODIMP ARMCommonModule::ARMZcFlat(double pZc, double pDate, BSTR pCcy, BSTR* pRet)
{
	return theInstance.ARMZcFlat(pZc, pDate,  pCcy,  pRet); 

}


STDMETHODIMP ARMCommonModule::ARMBsModel(double pDate, double pSpot, BSTR pDividend, BSTR pDiscrate, BSTR pVol, BSTR pTypeStk, BSTR *pRet)
{
	return theInstance.ARMBsModel(pDate, pSpot,  pDividend,  pDiscrate,  pVol,  pTypeStk,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMBsSlModel(double pDate, BSTR pZc, BSTR pVolSpreadLock, BSTR pCvCapVol, BSTR pCvIndexVol, BSTR *pRet)
{
	return theInstance.ARMBsSlModel(pDate,  pZc,  pVolSpreadLock,  pCvCapVol,  pCvIndexVol,  pRet) ; 


}


STDMETHODIMP ARMCommonModule::ARMSwitchToETK()
{
	return theInstance.ARMSwitchToETK() ;
}

STDMETHODIMP ARMCommonModule::ARMSwitchToWSETK()
{
	return theInstance.ARMSwitchToWSETK(); 
}

STDMETHODIMP ARMCommonModule::ARMShutdownETK()
{
	return theInstance.ARMShutdownETK(); 
}


STDMETHODIMP ARMCommonModule::ARMSwitchToFLATFILE()
{
	return theInstance.ARMSwitchToFLATFILE() ;
}


STDMETHODIMP ARMCommonModule::ARMInfocentreConnect()
{
	return theInstance.ARMInfocentreConnect() ; 
}


STDMETHODIMP ARMCommonModule::ARMBaseReplicationConnect()
{
	return theInstance.ARMBaseReplicationConnect() ;
}

STDMETHODIMP ARMCommonModule::ARMZCLINT(VARIANT* pMatu, VARIANT* pRate, BSTR pMeth, double pDate, BSTR pCurrency, BSTR pInterpMeth,double pMatuAreDoubles,VARIANT *pSTerms, BSTR *pRet)
{
	return theInstance.ARMZCLINT( pMatu,  pRate,  pMeth, pDate,  pCurrency,  pInterpMeth, pMatuAreDoubles,pSTerms, pRet) ; 
	
}


STDMETHODIMP ARMCommonModule::ARM_zcspreaded(BSTR zcSprId,
										BSTR zcInitId,
										double date,
										BSTR MMFreq,
										BSTR SwapFreq,
										BSTR ccyId,
										BSTR *pRet)
{
	return theInstance.ARM_zcspreaded( zcSprId,
										 zcInitId,
										date,
										 MMFreq,
										 SwapFreq,
										 ccyId,
										 pRet) ; 
	

}


STDMETHODIMP ARMCommonModule::ARMCreateZCSwapInt(double pDate,VARIANT *pMatu,VARIANT *pRate,BSTR pMMVsFut,
										   BSTR pSwapVsFut,BSTR pRaw,BSTR pInterp,BSTR pCcy,
										   BSTR pRefObj,BSTR *pRet)
{
	return theInstance.ARMCreateZCSwapInt(pDate,pMatu,pRate, pMMVsFut,
										    pSwapVsFut, pRaw, pInterp, pCcy,
										    pRefObj, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMGetInitialCurveFromSummit( BSTR pIndex, BSTR pCurrency,  
													  BSTR pCvName,  double pDate, BSTR pAdjOrNot,  
													  VARIANT *pRetMat, VARIANT *pRetRate)
{
	return theInstance.ARMGetInitialCurveFromSummit(  pIndex,  pCurrency,  
													   pCvName,  pDate,  pAdjOrNot,  
													  pRetMat, pRetRate) ;

}


STDMETHODIMP ARMCommonModule::ARMGetInitialCurveFromCalypso(BSTR pIndex,
												  BSTR pCurrency,
												  BSTR pTerm,
												  BSTR pricingEnv,
												  double pDate,
												  BSTR forceCurveName,
												  BSTR xmlFile, 
												  BSTR pDoAdj, 
												  VARIANT *pRetMat,
												  VARIANT *pRetRate) 
{
	return theInstance.ARMGetInitialCurveFromCalypso(  pIndex,  
														pCurrency,  
														pTerm, 
														pricingEnv,
														pDate,
														forceCurveName,
														xmlFile,
														pDoAdj,
														pRetMat,
														pRetRate) ;

}


STDMETHODIMP ARMCommonModule::ARMGetInitialVolFromSummit(BSTR pIndex, BSTR pCurrency,
												   BSTR pCvName,  double pDate, BSTR pType,
												   BSTR pMatIndex, VARIANT *pRetMat, VARIANT *pRetTenor,
												   VARIANT* pRetVol)
{
	return theInstance.ARMGetInitialVolFromSummit( pIndex,  pCurrency,
												    pCvName,  pDate,  pType,
												    pMatIndex, pRetMat, pRetTenor,
												   pRetVol); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_CreateBasketCorrelMkDataFromCalypso(BSTR pricingEnv,
																			 double date,
																			 BSTR forceCurveName,
																			 BSTR Ccy,
																			 BSTR xmlFileName,
																			 BSTR indexId,
																			 BSTR *pRet)
{

	return theInstance.ARM_Credit_CreateBasketCorrelMkDataFromCalypso(pricingEnv,
														date, 
														forceCurveName, 
														Ccy,
														xmlFileName, 
														indexId,
														pRet) ; 
}
STDMETHODIMP ARMCommonModule::ARM_Credit_GetBasketCorrelMkDataFromCalypso(BSTR pricingEnv,
																  double date,
																  BSTR forceCurveName, 
																  BSTR xmlFileName,
																  VARIANT *pRetMat, 
																  VARIANT *pRetTenor,
																  VARIANT* pRetVol)
{
	return theInstance.ARM_Credit_GetBasketCorrelMkDataFromCalypso(pricingEnv,
														date, 
														forceCurveName,  
														xmlFileName, 
														pRetMat,
														pRetTenor,
														pRetVol) ; 
}



STDMETHODIMP ARMCommonModule::ARMGetInitialFXVolFromSummit(BSTR pCcy1,BSTR pCcy2,double pDate,BSTR pCvName,
													 BSTR pImpOrHist,BSTR pVolType,VARIANT *pRetMat,
													 VARIANT *pRetTenor,VARIANT *pRetVol)
{
	return theInstance.ARMGetInitialFXVolFromSummit( pCcy1, pCcy2,pDate, pCvName,
													  pImpOrHist, pVolType,pRetMat,
													 pRetTenor,pRetVol) ;

}


STDMETHODIMP ARMCommonModule::ARMTHREEMONTHFUT(BSTR pDelivery,
										 long pMarket,
										 BSTR pCcy,
										 BSTR pRefObj,
										 BSTR *pRet)
{
	return theInstance.ARMTHREEMONTHFUT( pDelivery,
										 pMarket,
										  pCcy,
										  pRefObj,
										  pRet) ;
}


STDMETHODIMP ARMCommonModule::ARMFutPibor(BSTR pDelivery,BSTR *pRet)
{
	return theInstance.ARMFutPibor( pDelivery, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMIRFUT(double pDelivery,BSTR pIdUnderlying,BSTR pRefObj,BSTR *pRet)
{
	return theInstance.ARMIRFUT(pDelivery, pIdUnderlying, pRefObj, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMLibor(BSTR pLiborTypeId,BSTR pCcyId, BSTR pResetFreqId, 
								 BSTR pPayFreqId,BSTR pRefObj,BSTR pBasis,BSTR pIntRule,BSTR *pRet)
{
	return theInstance.ARMLibor( pLiborTypeId, pCcyId,  pResetFreqId, 
								  pPayFreqId, pRefObj, pBasis, pIntRule, pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARMLiborSwaption(double pStartDate,
										 double pEndDate,
										 BSTR pReceiveOrPay,
										 double pStrike, 
										 double pMaturity,
										 BSTR pLiborType,
										 double pSpread,
										 BSTR pExerciseType,
										 BSTR pResetFreq,
										 BSTR pPayFreq,
										 BSTR pCcyId,
										 BSTR pRefObj,
										 BSTR *pRet)
{
	return theInstance.ARMLiborSwaption(pStartDate,
										 pEndDate,
										  pReceiveOrPay,
										 pStrike, 
										 pMaturity,
										  pLiborType,
										 pSpread,
										  pExerciseType,
										  pResetFreq,
										  pPayFreq,
										  pCcyId,
										  pRefObj,
										  pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMFixedLeg(double pStartDate, double pEndDate, BSTR pReceiveOrPay, double pFixRate, BSTR pDayCount, BSTR pFreq, BSTR pDecompFreq, BSTR pPayTiming, BSTR pIntRule, BSTR pStubRule, BSTR pCcyId, BSTR pPayCalName, BSTR pNxChange, double pRefDate, BSTR pAdjStartDate,BSTR *pRet)
{
	return  theInstance.ARMFixedLeg(pStartDate,pEndDate, pReceiveOrPay, 
									pFixRate, pDayCount, 
									 pFreq, pDecompFreq, pPayTiming,  
									 pIntRule, pStubRule, pCcyId, 
									 pPayCalName, pNxChange, pRefDate, pAdjStartDate, pRet) ; 
}



STDMETHODIMP ARMCommonModule::ARMAswPrice(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pMargin, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet)
{
	return theInstance.ARMAswPrice(pMaturity, pCpn, pFreq,  pBase, pMargin, pRedemptionPrice, pAsOf, pDelivery,  pFixDecompfreq,  pCcy1,  pIndex1,  pFwdCurve1,  pDiscCurve1,  pCcy2,  pIndex2,  pFwdCurve2,  pDiscCurve2,  pAmortizationId, pSolve, pMinValue, pMaxValue, pRet) ;
}


STDMETHODIMP ARMCommonModule::ARMAswMargin(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pPrice, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet)
{
	return theInstance.ARMAswMargin(pMaturity, pCpn,  pFreq,  pBase, pPrice, pRedemptionPrice, pAsOf, pDelivery,  pFixDecompfreq,  pCcy1,  pIndex1,  pFwdCurve1,  pDiscCurve1,  pCcy2,  pIndex2,  pFwdCurve2,  pDiscCurve2,  pAmortizationId, pSolve, pMinValue, pMaxValue, pRet); 

}

STDMETHODIMP ARMCommonModule::ARMFrnMargin(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pPrice, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet)
{
	return theInstance.ARMFrnMargin(pAsOf, pDelivery, pMaturity,  pCcy1,  pIndex1,  pFwdCurve1,  pDiscCurve1, pFacialMargin, pPrice,  pCcy2,  pIndex2,  pFwdCurve2,  pDiscCurve2, pFixing, pSpread, pOutMode, pSolve,  pAmortizationId, pRet); 

}


STDMETHODIMP ARMCommonModule::ARMFrnPrice(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pValoMargin, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet)
{
	return theInstance.ARMFrnPrice(pAsOf, pDelivery, pMaturity,  pCcy1,  pIndex1,  pFwdCurve1,  pDiscCurve1, pFacialMargin, pValoMargin,  pCcy2,  pIndex2,  pFwdCurve2,  pDiscCurve2, pFixing, pSpread, pOutMode, pSolve,  pAmortizationId, pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMRefValue(VARIANT* pdates,
									VARIANT* pvalues,
									VARIANT* pvalues2,
									long valueType,
									long conversion,
									BSTR calcMethod,
									BSTR *pRefVal)
{
	return theInstance.ARMRefValue( pdates,
									 pvalues,
									 pvalues2,
									valueType,
									conversion,
									 calcMethod,
									 pRefVal); 

}

STDMETHODIMP ARMCommonModule::ARMCreateGenCorrelManager(VARIANT* pMktTags,
												  VARIANT* pIntraMktTags,
												  VARIANT* pCorrelCurveIds,
												  BSTR *pRet)
{
	return theInstance.ARMCreateGenCorrelManager( pMktTags,
												   pIntraMktTags,
												   pCorrelCurveIds,
												   pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMBSConvAdjust(BSTR pSUMMITFormulaeUsed,
											  BSTR pUseSABRCMS,
										 BSTR *pRet)
{
	return theInstance.ARMBSConvAdjust( pSUMMITFormulaeUsed,
										pUseSABRCMS,
										pRet) ;


}


STDMETHODIMP ARMCommonModule::ARMBsModelGen(BSTR pYieldCurve,
									  BSTR pVolatility,
									  BSTR pCorrMgr,
									  BSTR pCnvxManager,
									  BSTR pCapletVol,
									  BSTR pSpreadLock,
									  BSTR pDiscCurve,
									  BSTR pCorrel,
									  BSTR pCashVol,
									  BSTR pSpreadVol,
									  BSTR pModelType,
									  BSTR pSpreadVolType,
									  BSTR pSabrMod,
									  BSTR pLnOrNorVol,
									  long pNumSteps,
									  BSTR *pRet)
{
	return theInstance.ARMBsModelGen( pYieldCurve,
									   pVolatility,
									   pCorrMgr,
									   pCnvxManager,
									   pCapletVol,
									   pSpreadLock,
									   pDiscCurve,
									   pCorrel,
									   pCashVol,
									   pSpreadVol,
									   pModelType,
									   pSpreadVolType,
									   pSabrMod,
									   pLnOrNorVol,
									   pNumSteps,
									   pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMIrIndex(BSTR pDaycount,BSTR pPayFreq,double pMaturity,BSTR pCompMethod,BSTR pFwdRule,BSTR pResetTiming,double pResetGap,BSTR pPayTiming,double pPayGap,BSTR pCcy,BSTR pIndexType,double pDecompFreq,BSTR pIntRule,BSTR pResetFreq, BSTR *pRet)
{
	return theInstance.ARMIrIndex( pDaycount, pPayFreq,pMaturity, pCompMethod, pFwdRule, pResetTiming,pResetGap, pPayTiming,pPayGap, pCcy, pIndexType,pDecompFreq, pIntRule, pResetFreq,  pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMSwapleg(BSTR pIndexId, double pStartDate, double pEndDate, BSTR pRecOrPay, VARIANT pSpread, BSTR pCcy, BSTR pDayCount, double pResetGap, BSTR pResetCal, BSTR pPayCal, double pDecompPricingFlag, BSTR pNxChange, BSTR pStubRule, double pRefDate, BSTR pAdjStartDate, BSTR *pRet)
{
	return theInstance.ARMSwapleg( pIndexId, pStartDate, pEndDate,  pRecOrPay, pSpread,  pCcy,  pDayCount, pResetGap,  pResetCal,  pPayCal, pDecompPricingFlag,  pNxChange,  pStubRule, pRefDate,  pAdjStartDate,  pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMConstRefvalue(double pValue,  BSTR *pRet)
{
	return theInstance.ARMConstRefvalue(pValue,  pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMBond(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedempPrice,double pPeriodicity,VARIANT pDaycount,double pSettleGap,double pCpnDateFlag,BSTR pCcy,BSTR *pRet)
{
	return theInstance.ARMBond(pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedempPrice,pPeriodicity, pDaycount,pSettleGap,pCpnDateFlag, pCcy, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMINFCreateOATLeg(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,double pCoMultiple,BSTR *pRet)
{
	return theInstance.ARMINFCreateOATLeg(pStartDate,pEndDate, pInfIdx, pRcvOrPay, pInterpType,pLeverage,pSpread, pResetFreq, pDaycount, pResetCal, pFwdRule, pIntRule, pStubRule,pResetNumGap,pResetDenomGap, pPayFreq,pPayGap, pPayCal, pFinalNotionalType,pFirstReset,pCoMultiple, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMSwap(BSTR pSwapleg1,BSTR pSwapleg2,double pMinPay,BSTR *pRet)
{
	return theInstance.ARMSwap( pSwapleg1, pSwapleg2,pMinPay, pRet) ; 


}


STDMETHODIMP ARMCommonModule::ARMPToYield(BSTR pBond,double pSettleDate,double pPrice,VARIANT *pRet)
{
	return theInstance.ARMPToYield( pBond,pSettleDate,pPrice,pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMYToPrice(BSTR pBond,double pSettleDate,double pYield,VARIANT *pRet)
{
	return theInstance.ARMYToPrice( pBond,pSettleDate,pYield,pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMYToDuration(BSTR pBond,double pSettleDate,double pActuRate,double pFlagCpn,VARIANT *pRet)
{
	return theInstance.ARMYToDuration( pBond,pSettleDate,pActuRate,pFlagCpn,pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMLiborleg(double pStartDate,double pEndDate,BSTR pLiborType,BSTR pRecOrPay,VARIANT pSpread,BSTR pResetFreq,BSTR pPayFreq,BSTR pResetTiming,BSTR pPayTiming,BSTR pCcy,BSTR pIntRule,double pResetGap,BSTR pResetCal,BSTR pPayCal,double pDecompPricingFlag,BSTR pNxChange,BSTR pStubRule,double pRefDate,BSTR pAdjStartDate,BSTR pCpnDaycount,BSTR *pRet)
{
	return theInstance.ARMLiborleg(pStartDate,pEndDate, pLiborType, pRecOrPay, pSpread, pResetFreq, pPayFreq, pResetTiming, pPayTiming, pCcy, pIntRule,pResetGap, pResetCal, pPayCal,pDecompPricingFlag, pNxChange, pStubRule,pRefDate, pAdjStartDate, pCpnDaycount, pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARMImpliedSpread(BSTR pSwap,BSTR pModel,double pPrice,double pLeg1or2,double *pRet)
{
	return theInstance.ARMImpliedSpread( pSwap, pModel,pPrice,pLeg1or2,pRet); 
}


STDMETHODIMP ARMCommonModule::ARMDiscountPrice(BSTR pZeroCurve,double pMatu,double *pRet)
{
	return theInstance.ARMDiscountPrice( pZeroCurve,pMatu,pRet); 

}


STDMETHODIMP ARMCommonModule::ARMINFCreateCurve(double pAsOf,BSTR pIndexName,double pCPIIndexValue,double pCPIIndexDate,VARIANT* pMatu,VARIANT* pRate,BSTR pMonthlyInterpType,BSTR pDailyInterpType,BSTR pDCFMonthly,BSTR pDCFDaily,BSTR pExtrapolType,BSTR pResetManager,BSTR pSeasonManager,BSTR *pRet)
{
	return theInstance.ARMINFCreateCurve(pAsOf, pIndexName,pCPIIndexValue,pCPIIndexDate, pMatu, pRate, pMonthlyInterpType, pDailyInterpType, pDCFMonthly, pDCFDaily, pExtrapolType, pResetManager, pSeasonManager, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMINFInterpCPI(BSTR pZc,double pCPIDate,BSTR pDCFlag,BSTR pDailyInterpType,BSTR pCPIlag,double pWeight,double *pRet)
{
	return theInstance.ARMINFInterpCPI( pZc,pCPIDate, pDCFlag, pDailyInterpType, pCPIlag,pWeight,pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMINFSeasonManager(VARIANT* pMonthList,VARIANT* pValues,BSTR pSeasonAdjMode,BSTR *pRet)
{
	return theInstance.ARMINFSeasonManager( pMonthList, pValues, pSeasonAdjMode, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMINFResetManager(VARIANT* pDatas,double pNbIndex,BSTR *pRet)
{
	return theInstance.ARMINFResetManager( pDatas,pNbIndex, pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMINFYcMod(BSTR pYieldCurve,BSTR pInfCurve,BSTR *pRet)
{
	return theInstance.ARMINFYcMod( pYieldCurve, pInfCurve, pRet) ; 

}



STDMETHODIMP ARMCommonModule::ARMAccrued(BSTR pSec,double pDate,BSTR pModel,double *pRet)
{
	return theInstance.ARMAccrued( pSec,pDate, pModel,pRet); 

}


STDMETHODIMP ARMCommonModule::ARMClonedAndSetNotional(BSTR bLegId,BSTR bAmortId, BSTR *pRet)
{
	return theInstance.ARMClonedAndSetNotional( bLegId, bAmortId, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARM_INF_GetZcFromSummit(BSTR Index,
											    BSTR Ccy,
												BSTR cvname,
												double date,
												BSTR seasonAdj,
												BSTR seasonAdjMode,
												BSTR *pRet)
{
	return theInstance.ARM_INF_GetZcFromSummit( Index,
											     Ccy,
												 cvname,
												date,
												 seasonAdj,
												 seasonAdjMode,
												 pRet) ; 

}

/*
STDMETHODIMP ARMCommonModule::ARM_INF_CreateGenericLeg(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,BSTR *pRet)
{
	return theInstance.ARM_INF_CreateGenericLeg(pStartDate,pEndDate, pInfIdx, pRcvOrPay, pInterpType,pLeverage,pSpread, pResetFreq, pDaycount, pResetCal, pFwdRule, pIntRule, pStubRule,pResetNumGap,pResetDenomGap, pPayFreq,pPayGap, pPayCal, pFinalNotionalType,pFirstReset, pRet) ; 

}
*/


STDMETHODIMP ARMCommonModule::ARMRiskyBond(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet)
{
	return theInstance.ARMRiskyBond(pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedemptionPrice,pPeriodicity, pDaycount,pSettleGap,pCpnDateFlag, pCcyId,pRepo,pSsl,pRecoveryRate, pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMRiskyBondWithCF(double pAsOfDate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,VARIANT *pYearTerms,VARIANT *pCashFlows,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet)
{
	return theInstance.ARMRiskyBondWithCF(pAsOfDate,pRedemptionPrice,pPeriodicity, pDaycount,pYearTerms,pCashFlows,pSettleGap,pCpnDateFlag, pCcyId,pRepo,pSsl,pRecoveryRate, pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMGetFixing(BSTR source,
									 BSTR index,
									 BSTR term,
									 BSTR ccy,
									 double date,
									 double *pRet)
{
	return theInstance.ARMGetFixing( source,
									  index,
									  term,
									  ccy,
									 date,
									 pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMGetFixingFromCalypso(BSTR source,
									 BSTR index,
									 BSTR term,
									 BSTR ccy,
									 BSTR curveName,
									 DATE  date,
									 VARIANT *pRet)
{
	return theInstance.ARMGetFixingFromCalypso( source,
									  index,
									  term,
									  ccy,
									  curveName,
									 date,
									 pRet) ;

}

STDMETHODIMP ARMCommonModule::ARMIndexIndexCorrelCube(VARIANT* pVolCurveId,
									 VARIANT* pTenors1List,
									 VARIANT* pTenors2List, 
									 BSTR pInterSurfInterp,
									 BSTR *pRet)
{
	return theInstance.ARMIndexIndexCorrelCube( pVolCurveId,
									  pTenors1List,
									  pTenors2List, 
									  pInterSurfInterp,
									  pRet) ; 

}




STDMETHODIMP ARMCommonModule::ARMCmsLeg(double startDate,
					  double endDate,
					  BSTR cmsTypeId,
					  BSTR receiveOrPay,
					  BSTR yieldDecompFreq,
					  BSTR swapLegDayCount,
					  BSTR resetFreq,
					  BSTR payFreq,
					  long resetGap,		 
					  BSTR intRule,
					  BSTR ccyName,
					  BSTR resetTiming,
					  BSTR stubRule,
					  BSTR adjStartDate,
					  BSTR *pRet)
{
	return theInstance.ARMCmsLeg(startDate,
					  endDate,
					   cmsTypeId,
					   receiveOrPay,
					   yieldDecompFreq,
					   swapLegDayCount,
					   resetFreq,
					   payFreq,
					  resetGap,		 
					   intRule,
					   ccyName,
					   resetTiming,
					   stubRule,
					   adjStartDate,
					   pRet) ;
}



STDMETHODIMP ARMCommonModule::ARM_ReplicConvAdjust_Create(
										BSTR Payoff_ReplicMode,
										double Payoff_StepOrReplicPrecision,
										BSTR Payoff_StopMode,
										double Payoff_StopThreshold,
										BSTR Sensi_ReplicMode,
										double Sensi_StepOrReplicPrecision,
										BSTR Sensi_StopMode,
										double Sensi_StopThreshold,
										BSTR UsedModelId,
										double StrikeMinReplic,
										double StrikeMaxReplic,
										BSTR *pRet)
{
	return theInstance.ARM_ReplicConvAdjust_Create(
										 Payoff_ReplicMode,
										Payoff_StepOrReplicPrecision,
										 Payoff_StopMode,
										Payoff_StopThreshold,
										 Sensi_ReplicMode,
										Sensi_StepOrReplicPrecision,
										 Sensi_StopMode,
										Sensi_StopThreshold,
										 UsedModelId,
										StrikeMinReplic,
										StrikeMaxReplic,
										 pRet) ; 

}

STDMETHODIMP ARMCommonModule::ARM_MapConvAdjust_Create(
										BSTR LiborArrearAdj,
										BSTR NaturalCMSAdj,
										BSTR PaymentLagAdj,
										BSTR *pRet)
{
	return theInstance.ARM_MapConvAdjust_Create(
										 LiborArrearAdj,
										 NaturalCMSAdj,
										 PaymentLagAdj,
										 pRet) ; 

}


STDMETHODIMP ARMCommonModule::ARMCreateGenCorrelatorManager(VARIANT *pMktTags,
									 VARIANT* pHyperDiagVol,
									 VARIANT *pIndexIndexVol,
									 VARIANT* pCorrelVol,
									 VARIANT* pIndexVol,
									 BSTR *pRet)
{
	return theInstance.ARMCreateGenCorrelatorManager(pMktTags,
									  pHyperDiagVol,
									 pIndexIndexVol,
									  pCorrelVol,
									  pIndexVol,
									  pRet) ;

}


STDMETHODIMP ARMCommonModule::ARMHyperCube(VARIANT *pVolCurvId,
									 VARIANT* pKeys,
									 BSTR *pRet)
{
	return theInstance.ARMHyperCube(pVolCurvId,
									pKeys,
									pRet) ;
}

STDMETHODIMP ARMCommonModule::ARMDisplaySchedule(BSTR pLegId,
										   BSTR pDataType,
										   VARIANT *pRet)
{
	return theInstance.ARMDisplaySchedule( pLegId,
										    pDataType,
										    pRet) ;
}

STDMETHODIMP ARMCommonModule::ARMGetCorrelFromSummit(BSTR pCcy1, BSTR pIndex1, BSTR pCcy2, BSTR pIndex2, double pDate, BSTR pCvName, BSTR *pRet)
{
	return theInstance.ARMGetCorrelFromSummit( pCcy1,  pIndex1,  pCcy2,  pIndex2,  pDate,  pCvName, pRet);
}
STDMETHODIMP ARMCommonModule::ARMcomputePentilix(double pNtl, 
										   double pStartDatePhase1,  
										   BSTR pCcy,
										   BSTR pIndexPay, 					
										   BSTR pIndexPhase1, 										   
										   double pSpreadPhase1,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG, 
										   double pStartDatePhase2,
										   double pSpreadPhase2,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   double pStartDatePhase3,
										   double pEndDatePhase3,
										   BSTR pIndexPhase3,
										   double pSpreadPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub,  
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 					
										   BSTR pHyperCubeCorrel,
										   BSTR pIndexIndexCorrelCube, 
										   BSTR pCorrEUR, 
										   BSTR pInterCorr, 
										   VARIANT *pRet
										   )
{
	return theInstance.ARMcomputePentilix(pNtl, 
						pStartDatePhase1,  
						pCcy,
						pIndexPay, 					
						pIndexPhase1, 										   
						pSpreadPhase1,
						pDayCountPhase1,
						pPayFreqPhase1,
						pResetFreqPhase1,
						pResetTimingPhase1,
						pRoll, 
						pAdjPhase1, 
						pStub,
						pIndexPhase2DIG,
						pIndexLongPhase2DIG,
						pStrikePhase2DIG,
						pResetTimingPhase2DIG,
						pAdjPhase2DIG, 
						pStartDatePhase2,
						pSpreadPhase2,										   
						pDayCountPhase2, 
						pPayFreqPhase2, 
						pResetFreqPhase2,
						pAdjPhase2,
						pStartDatePhase3,
						pEndDatePhase3,
						pIndexPhase3,
						pSpreadPhase3, 
						pDayCountPhase3, 
						pPayFreqPhase3, 
						pResetFreqPhase3,
						pResetTimingPhase3,
						pAdjPhase3,
						pIndexFund,
						pSpreadFund,
						pDayCountFund, 
						pPayFreqFund, 
						pResetFreqFund,
						pResetTimingFund,
						pAdjFund, 
						pEndDateAmort,
						pDayCountAmort,
						pIntRuleAmort, 
						pTxAmort,
						pFreqAmort,
						pAmountAmort, 
						pTypeAmort, 
						pFee, 
						pVolCurvFromMatriceShift,
						pVol, 
						pVolCub,  
						pConvexityManager, 
						pZc, 
						pSmiledMod, 
						pSmiledModBump, 					
						pHyperCubeCorrel,
						pIndexIndexCorrelCube, 
						pCorrEUR, 
						pInterCorr, 
						pRet) ;
}


STDMETHODIMP ARMCommonModule::ARM_Credit_RiskyPV01 (BSTR pDefCurve,
													VARIANT Date1,
													VARIANT Date2,
													VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_RiskyPV01 (pDefCurve,Date1,Date2,pRet) ; 
}


STDMETHODIMP ARMCommonModule::ARM_Credit_RiskyPV01AsSensitivity (BSTR pDefCurve,
													BSTR tenor,
													VARIANT *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARM_Credit_RiskyPV01AsSensitivity (pDefCurve,tenor,pRet) ; 
}

STDMETHODIMP ARMCommonModule::ARMGetResetMgrFromSummit (double pAsOf,
														 BSTR pIndex,
														 BSTR pSource,
														 BSTR pCcy,
														 BSTR pIsInflationIndex,
														 BSTR pTerm,
														 BSTR *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARMGetResetMgrFromSummit (pAsOf,
												  pIndex,
												  pSource,
												  pCcy,
												  pIsInflationIndex,
												  pTerm,
												  pRet);
}


STDMETHODIMP ARMCommonModule::ARMGetReset (BSTR pResetMgr,
											double pDate,
											double *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARMGetReset (pResetMgr,
									 pDate,
									 pRet);
}


STDMETHODIMP ARMCommonModule::ARMSetLastFixing (BSTR pSecurityId,
											double pRate,
											double pAsOf,
											double pBeforeLastFixingDate, 
											double pResetDate,
											BSTR *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARMSetLastFixing (pSecurityId,
									 pRate,
									 pAsOf,
									 pBeforeLastFixingDate,
									 pResetDate,
									 pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_Sectorial_Correlation(DATE AsOf,
															   BSTR StructName,
															   BSTR correlation_Type,
															VARIANT * vLabels,
															VARIANT * vector_Membership,
															double intra_Sector_Correlation,
															double inter_Sector_Correlation,
															VARIANT * vBetas,
															VARIANT * vLambdas,
															VARIANT * vBetas_Down,
															VARIANT * vLambdas_Down,
															VARIANT *pRet)
{
	return theInstance.ARM_Credit_Sectorial_Correlation(AsOf,
														StructName,
														correlation_Type,
														vLabels,
														vector_Membership,
														intra_Sector_Correlation,
														inter_Sector_Correlation,
														vBetas,
														vLambdas,
														vBetas_Down,
														vLambdas_Down,
														pRet);
}


STDMETHODIMP ARMCommonModule::ARMcomputeReviPentix(double pNtl, 
										   VARIANT pDate,  
										   BSTR pCcy,
										   BSTR pIndexPhase1, 										   
										   VARIANT pSpread,
										   BSTR pDayCountPhase1,
										   BSTR pPayFreqPhase1,
										   BSTR pResetFreqPhase1,
										   BSTR pResetTimingPhase1,
										   BSTR pRoll, 
										   BSTR pAdjPhase1, 
										   BSTR pStub,
										   BSTR pIndexPhase2DIG,
										   BSTR pIndexLongPhase2DIG,
										   BSTR pStrikePhase2DIG,
										   BSTR pResetTimingPhase2DIG,
										   BSTR pAdjPhase2DIG,										   
										   BSTR pDayCountPhase2, 
										   BSTR pPayFreqPhase2, 
										   BSTR pResetFreqPhase2,
										   BSTR pAdjPhase2,
										   BSTR pIndexPhase3, 
										   BSTR pDayCountPhase3, 
										   BSTR pPayFreqPhase3, 
										   BSTR pResetFreqPhase3,
										   BSTR pResetTimingPhase3,
										   BSTR pAdjPhase3,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund,  
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pFloorOrCap, 
										   double pFee, 
										   double pLevier, 
										   double pTxFixeMax, 
										   BSTR pIsCapped, 
										   double pTxCap,
										   BSTR pVolCurvFromMatriceShift,
										   BSTR pVol, 
										   BSTR pVolCub, 
										   BSTR pCorrManager, 
										   BSTR pConvexityManager, 
										   BSTR pZc, 
										   BSTR pSmiledMod, 					
										   BSTR pHyperCubeCorrel,
										   VARIANT *pBumpBsGenMod,
										   VARIANT *pBumpVolBsGenMod, 
										   VARIANT *pRet
										  )

{	return theInstance.ARMcomputeReviPentix(pNtl, 
										   pDate,  
										   pCcy,
										   pIndexPhase1, 										   
										   pSpread,
										   pDayCountPhase1,
										   pPayFreqPhase1,
										   pResetFreqPhase1,
										   pResetTimingPhase1,
										   pRoll, 
										   pAdjPhase1, 
										   pStub,
										   pIndexPhase2DIG,
										   pIndexLongPhase2DIG,
										   pStrikePhase2DIG,
										   pResetTimingPhase2DIG,
										   pAdjPhase2DIG, 										   
										   pDayCountPhase2, 
										   pPayFreqPhase2, 
										   pResetFreqPhase2,
										   pAdjPhase2,
										   pIndexPhase3, 
										   pDayCountPhase3, 
										   pPayFreqPhase3, 
										   pResetFreqPhase3,
										   pResetTimingPhase3,
										   pAdjPhase3,
										   pIndexFund,
										   pSpreadFund,
										   pDayCountFund, 
										   pPayFreqFund, 
										   pResetFreqFund,
										   pResetTimingFund,
										   pAdjFund,  
										   pDayCountAmort,
										   pIntRuleAmort, 
										   pTxAmort,
										   pFreqAmort,
										   pAmountAmort, 
										   pTypeAmort, 
										   pFloorOrCap, 
										   pFee, 
										   pLevier, 
										   pTxFixeMax, 
										   pIsCapped,
										   pTxCap,
										   pVolCurvFromMatriceShift,
										   pVol,
										   pVolCub, 
										   pCorrManager, 
										   pConvexityManager, 
										   pZc, 
										   pSmiledMod, 					
										   pHyperCubeCorrel,
										   pBumpBsGenMod,
										   pBumpVolBsGenMod, 
										   pRet
										   );
										   }

STDMETHODIMP ARMCommonModule::ARMGenAmortization (BSTR pSwaplegId,
												  BSTR pAmortMethod,
												  BSTR pAmortFreq,
												  double pAmortAmount,
												  BSTR pDaycount,
												  double pLegNotional,
												  double pAmortRate,
												  double pReducedMaturity,
												  BSTR pModelId,
												  double pCleanUp,
												  BSTR *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARMGenAmortization (pSwaplegId,
										   pAmortMethod,
										   pAmortFreq,
										   pAmortAmount,
										   pDaycount,
										   pLegNotional,
										   pAmortRate,
										   pReducedMaturity,
										   pModelId,
										   pCleanUp,
										   pRet);
}



STDMETHODIMP ARMCommonModule::ARMCptRefvalue (BSTR pRefvalId,
											double pDate,
											double *pRet)
{
	// TODO: Add your implementation code here
	return theInstance.ARMCptRefvalue (pRefvalId,
									 pDate,
									 pRet);
}

STDMETHODIMP ARMCommonModule::ARMCalypsoDevConnect()
{
	return theInstance.ARMCalypsoDevConnect() ; 
}

STDMETHODIMP ARMCommonModule::ARMCalypsoProdConnect()
{
	return theInstance.ARMCalypsoProdConnect() ; 
}

STDMETHODIMP ARMCommonModule::ARMCalypsoRecConnect()
{
	return theInstance.ARMCalypsoRecConnect() ; 
}

STDMETHODIMP ARMCommonModule::Local_ARM_ProdConnect()
{
	return theInstance.Local_ARM_ProdConnect() ; 
}

STDMETHODIMP ARMCommonModule::ARM_SetDefaultCurrency(BSTR isoCCy, VARIANT *pRet)
{
	return theInstance.ARM_SetDefaultCurrency(isoCCy,pRet ) ; 
}

STDMETHODIMP ARMCommonModule::ARMGetZCFromCalypso(BSTR pIndex,
												  BSTR pCurrency,
												  BSTR pTerm,
												  BSTR pricingEnv,
												  double pDate,
												  BSTR pInterpMethod, 
												  BSTR forceCurveName,
												  BSTR xmlFile, 
												  BSTR *pRet)
{
	return  theInstance.ARMGetZCFromCalypso( pIndex,
												   pCurrency,
												   pTerm,
												   pricingEnv,
												   pDate,
												   pInterpMethod, 
												   forceCurveName,
												   xmlFile, 
												  pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_DefProbInverse(BSTR pCurveId,
													double dDefProba,
													VARIANT *pRet)
{
	return theInstance.ARM_Credit_DefProbInverse(pCurveId,
												dDefProba,
												pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_QMatrix(VARIANT* pQMatrix, 
											   VARIANT *pRet) 
{
	return theInstance.ARM_Credit_QMatrix(pQMatrix,pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_MarketDataMng(VARIANT *pstrVect, 
									  VARIANT *pRet)
{
	return theInstance.ARM_Credit_MarketDataMng(pstrVect, pRet);
}


STDMETHODIMP ARMCommonModule::ARMLivretALeg(double pStartDate,
											double pEndDate,
											BSTR pRcvOrPay,
											VARIANT pSpread,
											BSTR pResetFreq,
											BSTR pPayFreq,
											BSTR pResetTiming,
											BSTR pPayTiming,
											BSTR pCcy,
											BSTR pIntRule,
											double pResetGap,
											BSTR pResetCal,
											BSTR pPayCal,
											double pDecompPricingFlag,
											BSTR pNxChange,
											BSTR pStubRule,
											double pRefDate,
											BSTR pAdjStartDate,
											BSTR pDayCount,
											BSTR *pRet)
{
	return theInstance.ARMLivretALeg(pStartDate,pEndDate,pRcvOrPay,pSpread,pResetFreq,pPayFreq,pResetTiming,pPayTiming,
									 pCcy,pIntRule,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,
									 pRefDate,pAdjStartDate,pDayCount,pRet);
}

STDMETHODIMP ARMCommonModule::ARMLivretACurve(double pAsOf,
											  BSTR pInfCurvId,
											  BSTR pEuribCurvId,
											  double pFlagRouding,
											  BSTR pInfResetMgrId,
											  BSTR pFixingLivretAId,
											  BSTR pFixingEuribId,
											  BSTR pMonthForAugust,
											  BSTR pMonthForFebruary,
											  BSTR *pRet)
{
	return theInstance.ARMLivretACurve(pAsOf,pInfCurvId,pEuribCurvId,pFlagRouding,pInfResetMgrId,pFixingLivretAId,pFixingEuribId,pMonthForAugust,pMonthForFebruary,pRet);
}

STDMETHODIMP ARMCommonModule::ARMcomputeLivretA(double pNtl, 
										   double pStartDateLeg1,  
										   double pEndDateLeg1,
										   BSTR pCcy,
										   BSTR pIndexLeg1, 										   
										   VARIANT pSpreadLeg1,
										   BSTR pDayCountLeg1,
										   BSTR pPayFreqLeg1,
										   BSTR pResetFreqLeg1,
										   BSTR pResetTimingLeg1,
										   BSTR pAdjLeg1,
										   BSTR pRoll, 
										   BSTR pStub,
										   double pEndDateLA,
										   double pSpreadLeg2,
										   BSTR pDayCountLA,
										   BSTR pPayFreqLA,
										   BSTR pResetFreqLA,
										   BSTR pResetTimingLA, 									   
										   BSTR pAdjLA, 
										   BSTR pIndexLeg2, 
										   BSTR pDayCountLeg2, 
										   BSTR pPayFreqLeg2, 
										   BSTR pResetFreqLeg2,
										   BSTR pResetTimingLeg2,
										   BSTR pAdjLeg2, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 
										   BSTR pLAMod, 
										   BSTR pLAModBump,
										   BSTR pLAModBumpInflation,
										   VARIANT *pResetMgrIds, 
										   VARIANT *pRet)
{
	return theInstance.ARMcomputeLivretA( pNtl,  pStartDateLeg1,  pEndDateLeg1,  pCcy,  pIndexLeg1, pSpreadLeg1, pDayCountLeg1,  pPayFreqLeg1, pResetFreqLeg1, pResetTimingLeg1, pAdjLeg1, pRoll, pStub, pEndDateLA, pSpreadLeg2, pDayCountLA, pPayFreqLA, pResetFreqLA, pResetTimingLA, pAdjLA,  pIndexLeg2, pDayCountLeg2, pPayFreqLeg2, pResetFreqLeg2, pResetTimingLeg2, pAdjLeg2, pEndDateAmort,  pDayCountAmort, pIntRuleAmort, pTxAmort, pFreqAmort, pAmountAmort, pTypeAmort, pFee, pSmiledMod, pSmiledModBump, pLAMod, pLAModBump, pLAModBumpInflation,pResetMgrIds, pRet);
}

STDMETHODIMP ARMCommonModule::ARMFutDelivery(BSTR pFut, BSTR pCcy, double *pRet)
{
	return theInstance.ARMFutDelivery(pFut, pCcy, pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_FwdSpread(BSTR defcurveId, 
												  double Maturity1,
												  double Maturity2,
												  double FwdStartDate,
												  double FwdEndDate,
												  BSTR VolId, 
												  double *pRet)
{
	return theInstance.ARM_Credit_FwdSpread(defcurveId, Maturity1, Maturity2,FwdStartDate,FwdEndDate, VolId,pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_FwdSpreadAsIndex(BSTR defcurveId, 
												  double Maturity1,
												  double Maturity2, 
												  double *pRet)
{
	return theInstance.ARM_Credit_FwdSpreadAsIndex(defcurveId, Maturity1, Maturity2,pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_Flat_Correlation(DATE AsOf, 
														  BSTR structName,
														  double correlValue,
														  BSTR idIndex1,
														  BSTR idIndex2,
														  VARIANT *pRet)
{
	return theInstance.ARM_Credit_Flat_Correlation(AsOf,structName,correlValue,idIndex1,idIndex2,pRet); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_CorridorLeg_Sche(
														double Notional,
														BSTR RecieveOrPay,
														BSTR RefValueSpreads,
									                    BSTR floatingIdx,
														double leverageFloatIdx,
														BSTR creditIdx,
														BSTR refvalueKUP,
														BSTR refvalueKDW,
														BSTR ScheduleInfoId,
														VARIANT* accondef,
														BSTR disc_ccy,
														BSTR Name,
														VARIANT		*pRet)
{
	return theInstance.ARM_Credit_CorridorLeg_Sche(Notional,RecieveOrPay, RefValueSpreads, floatingIdx, leverageFloatIdx, creditIdx,
													refvalueKUP, refvalueKDW, ScheduleInfoId, accondef,disc_ccy,Name, pRet);
}


 STDMETHODIMP ARMCommonModule::ARM_Credit_Schedule_Info (double 	EffectiveDate,
																double 	MaturityDate,			
																BSTR 	payFrequency,
																BSTR 	ResetFreq ,
																BSTR 	DayCount,
																BSTR 	Stubrule,
																BSTR 	intRule,
																BSTR 	payCalName,
																BSTR 	PayTiming,
																BSTR 	ResetTiming,
																BSTR 	fwdRule,
																BSTR 	IncludeMaturity,
																BSTR 	adj,
																BSTR 	intStartAdj,
																BSTR 	AccDayCount,
																double 	ReferenceDate,
																double 	FirstCpnEffDate,
																BSTR 	AdjCal,
																VARIANT		*pRet)
 {
		return theInstance.ARM_Credit_Schedule_Info (EffectiveDate,	MaturityDate,payFrequency,	ResetFreq ,DayCount,Stubrule,
													intRule,payCalName,PayTiming,ResetTiming,fwdRule,IncludeMaturity,adj,
													intStartAdj,AccDayCount,ReferenceDate,	FirstCpnEffDate,AdjCal,pRet);
 }
STDMETHODIMP ARMCommonModule::ARMInfCurveSetResetMgr(BSTR pInfCurve, BSTR pResetMgr,BSTR *pRet)
{
	return theInstance.ARMInfCurveSetResetMgr(pInfCurve, pResetMgr, pRet); 
}


STDMETHODIMP ARMCommonModule::ARMcomputeTxFixed(double pNtl, 
										   double pStartDateLeg1,  
										   double pEndDateLeg1,
										   BSTR pCcy,
										   BSTR pIndexLeg1, 										   
										   VARIANT pSpreadLeg1,
										   BSTR pDayCountLeg1,
										   BSTR pPayFreqLeg1,
										   BSTR pResetFreqLeg1,
										   BSTR pResetTimingLeg1,
										   BSTR pAdjLeg1,
										   BSTR pRoll, 
										   BSTR pStub,
										   double pEndDateFixed,
										   double pSpreadLeg2,
										   BSTR pDayCountFixed,
										   BSTR pPayFreqFixed,
										   BSTR pResetFreqFixed,
										   BSTR pResetTimingFixed, 									   
										   BSTR pAdjFixed, 
										   BSTR pIndexLeg2, 
										   BSTR pDayCountLeg2, 
										   BSTR pPayFreqLeg2, 
										   BSTR pResetFreqLeg2,
										   BSTR pResetTimingLeg2,
										   BSTR pAdjLeg2, 
										   double pEndDateAmort,
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   double pFee, 
										   BSTR pSmiledMod, 
										   BSTR pSmiledModBump, 
										   VARIANT *pRet)
{
	return theInstance.ARMcomputeTxFixed( pNtl,  pStartDateLeg1,  pEndDateLeg1,  pCcy,  pIndexLeg1, pSpreadLeg1, pDayCountLeg1,  pPayFreqLeg1, pResetFreqLeg1, pResetTimingLeg1, pAdjLeg1, pRoll, pStub, pEndDateFixed, pSpreadLeg2, pDayCountFixed, pPayFreqFixed, pResetFreqFixed, pResetTimingFixed, pAdjFixed,  pIndexLeg2, pDayCountLeg2, pPayFreqLeg2, pResetFreqLeg2, pResetTimingLeg2, pAdjLeg2, pEndDateAmort,  pDayCountAmort, pIntRuleAmort, pTxAmort, pFreqAmort, pAmountAmort, pTypeAmort, pFee, pSmiledMod, pSmiledModBump, pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_CPDO(BSTR pRiskyLeg,
											 BSTR pRollLeg, 
											 BSTR pNoRiskyLeg, 
											 double pInitialValo, 
											 double pTarget,
											 double pMaturity, 
											 BSTR pCpnType, 
											 double	pUFFees, 
											 double	pRunningFees, 
											 double pVExpo, 
											 double	pV0Expo, 	
											 double	pAlpha,
											 double pBeta,	
											 double	pDesactivation,	
											 int	pNbAssets,
											 VARIANT *pRet)
{
return 
theInstance.ARM_Credit_CPDO(pRiskyLeg,
						pRollLeg, 
						pNoRiskyLeg, 
						pInitialValo, 
						pTarget,
						pMaturity, 
						pCpnType, 
						pUFFees, 
						pRunningFees, 
						pVExpo, 
						pV0Expo, 	
						pAlpha,
						pBeta,	
						pDesactivation,	
						pNbAssets,
						pRet) ;
}

STDMETHODIMP ARMCommonModule::ARM_Credit_PriceVector (VARIANT *pPricer , 
													BSTR pCPTTYPE,
													VARIANT *pRetVectorValos)
{
return 
theInstance.ARM_Credit_PriceVector(pPricer , 
							pCPTTYPE,
							pRetVectorValos);

}

STDMETHODIMP ARMCommonModule::ARM_Credit_GenPrice (VARIANT *pPricer , 
										  BSTR pCPTTYPE,
										  VARIANT *pRet)
{
return 
theInstance.ARM_Credit_GenPrice(pPricer , 
							pCPTTYPE,
							pRet);

}

STDMETHODIMP ARMCommonModule::ARM_Credit_FixingCurve(VARIANT *pDates, 
														VARIANT *pValues, 
														double  AsOfDate,
														BSTR	B_IndexName,
														BSTR	B_IndexID, 
														VARIANT *pRet) 
{

return  theInstance.ARM_Credit_FixingCurve(pDates, pValues, AsOfDate,
										   B_IndexName, B_IndexID, pRet);

}

STDMETHODIMP ARMCommonModule::ARM_Credit_CptInterpolDefCurveOLD(BSTR pCurve,BSTR pTENOR,double pSlope,double pDate,double pInterpDate, VARIANT *pRet)
{
	return theInstance.ARM_Credit_CptInterpolDefCurveOLD( pCurve, pTENOR,pSlope,pDate,pInterpDate, pRet); 
}


STDMETHODIMP ARMCommonModule::ARMSetCalendar(BSTR pFileName,double *pRet)
{
	return theInstance.ARMSetCalendar( pFileName, pRet); 
}


STDMETHODIMP ARMCommonModule::ARMInitGigaSpaces(BSTR pUrl,double *pRet)
{
	return theInstance.ARMInitGigaSpaces( pUrl, pRet); 
}

STDMETHODIMP ARMCommonModule::ARM_Credit_GetExpectedLoss (BSTR 	pricerId,
														double 	YearTerm, 
														double		*pRet)
{
	return theInstance.ARM_Credit_GetExpectedLoss (pricerId, YearTerm, pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_FunctionRegister (long address)
{
	return theInstance.ARM_Credit_FunctionRegister (address);
}

STDMETHODIMP ARMCommonModule::ARMBermudanXStyle (VARIANT *pxDates,VARIANT *pexpiryDates,BSTR *pRet)
{
	return theInstance.ARMBermudanXStyle (pxDates,pexpiryDates,pRet);
}

STDMETHODIMP ARMCommonModule::ARMcomputeCRA(BSTR pFixorFloat,
											double pFee,
										  double pAsOf, 
											double pStartDate,
											double pEndDate,
											BSTR pCcy,
											double pLevelUp,
											BSTR pUpSpec,
											double pLevelDown,
											BSTR pDownSpec,
											BSTR pRefIndex,
											BSTR pDayCount,
											BSTR pPayFreqPayIndex,
											BSTR pResetFreqRefIndex,
											BSTR pPaidRstTiming,
											BSTR pRefRstTiming,
											BSTR pStubRule,
											BSTR pPOrR,
											double pStartCallDate,
											BSTR pXStyle,
											BSTR pFundingIndex,
											BSTR pResetFreqFunding,
											BSTR pPayFreqFunding,
											VARIANT pSpreadFunding,
											BSTR pPOrRFunding,
											double pDecompPricingFlag,
											double pdiscMarginFactor,
											BSTR pPreInitFlag,
											double pMeanReversion,
											VARIANT pCalibParams,
											VARIANT pCalibParamsPF,
											double pKernelToGP,
											VARIANT pMarkovTreeParams,
											double pMarkovTreePathNumber,
					                        BSTR pBsmodId,
								            BSTR pBsmodSwoptId,
											BSTR pBsmodSwoptBumpId,
											BSTR pzcId,
											VARIANT *pRet)
{
	return theInstance.ARMcomputeCRA(pFixorFloat, pFee, pAsOf, pStartDate,	pEndDate,pCcy,pLevelUp,pUpSpec,pLevelDown,pDownSpec,pRefIndex,pDayCount,pPayFreqPayIndex,pResetFreqRefIndex,pPaidRstTiming, pRefRstTiming, pStubRule, pPOrR, pStartCallDate,	 pXStyle, pFundingIndex, pResetFreqFunding, pPayFreqFunding, pSpreadFunding,pPOrRFunding, pDecompPricingFlag,pdiscMarginFactor,pPreInitFlag,pMeanReversion, pCalibParams,pCalibParamsPF,pKernelToGP,pMarkovTreeParams,pMarkovTreePathNumber, pBsmodId,pBsmodSwoptId,pBsmodSwoptBumpId,pzcId,pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_Math_BivNormale (double x, double y, double rho, double *pRet)
{
	return theInstance.ARM_Credit_Math_BivNormale (x,y,rho,pRet);
}

STDMETHODIMP ARMCommonModule::ARMSetDiscountPricingMode(BSTR pModelId,int pDiscountPricingMode,BSTR *pRet)
{
	return theInstance.ARMSetDiscountPricingMode(pModelId, pDiscountPricingMode, pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_Math_RandUniform (double seed, double *pRet)
{
	return theInstance.ARM_Credit_Math_RandUniform (seed,pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_Math_Interpol(VARIANT *X ,
													 VARIANT *Y,
													 double value,
													 double type,
													 double smooth,
													 VARIANT *Weights,
													 double modeSpline,
													 double withC1condition,
													 double leftSlope,
													 double rightSlope,
													 double *pRet)
{
	return theInstance.ARM_Credit_Math_Interpol(X ,
												Y,
												value,
												type,
												smooth,
												Weights,
												modeSpline,
												withC1condition,
												leftSlope,
												rightSlope,
												pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_Random_Generator (BSTR RandomType, 
														 BSTR ParamId, 
														 VARIANT *pRet)
{

	return theInstance.ARM_Credit_Random_Generator(RandomType ,
												ParamId,
												pRet);

}
STDMETHODIMP ARMCommonModule::ARM_Credit_GenerateOneRandom(BSTR RandomId, 
													double *pRet)
{

	return theInstance.ARM_Credit_GenerateOneRandom(RandomId,
												pRet);

}
STDMETHODIMP ARMCommonModule::ARM_Credit_GenerateRandoms(BSTR RandomId, 
														 int DimVector,
														 VARIANT *pRet)
{

	return theInstance.ARM_Credit_GenerateRandoms(RandomId,
												DimVector,
												pRet);

}


STDMETHODIMP ARMCommonModule::ARM_Credit_ResetRandom(BSTR RandomId)
{

	return theInstance.ARM_Credit_ResetRandom(RandomId);

}

STDMETHODIMP ARMCommonModule::ARMPF (VARIANT *pinsts,
								  VARIANT *pcoeffs,
								  VARIANT *pmarketPrices,
								  VARIANT *pprecisions,
								  BSTR *pRet)
{
	return theInstance.ARMPF (pinsts,pcoeffs,pmarketPrices,pprecisions,pRet);
}

STDMETHODIMP ARMCommonModule::ARMBondTEC(double pIssueDate,  double pMaturityDate,  double pFirstCpnDate,  double pCpnRate,  double pRedempPrice, long pPeriodicity,VARIANT pDaycount,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double ptec,BSTR pPFTecId,BSTR pModTecId,BSTR *pRet)
{
		return theInstance.ARMBondTEC( pIssueDate,   pMaturityDate,   pFirstCpnDate,   pCpnRate,   pRedempPrice,  pPeriodicity, pDaycount, pSettleGap, pCpnDateFlag, pCcyId, ptec, pPFTecId, pModTecId, pRet);
}

STDMETHODIMP ARMCommonModule::ARMPFModFit(BSTR pmodName,BSTR ppf,double psettlement,BSTR pzc,VARIANT *pvList,VARIANT *pfList,long nag_algo,long pstep,double phorizon,BSTR *pRet)
{
		return theInstance.ARMPFModFit( pmodName, ppf, psettlement, pzc, pvList, pfList, nag_algo, pstep, phorizon,pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_Restrikable_CDO(double TriggerStartDate,
													   double Expiry,
													   double Strike,
													   double InitSpread,
													   int	OptionType,
													   BSTR pUnderlying,
													   double Rehauss,
													   BSTR TriggerFreq,
													   int DiffCDO,
													   int IsCMSpread,
													   double CMSpreadMatu,
													   VARIANT *pRet)
{
	return theInstance.ARM_Credit_Restrikable_CDO(TriggerStartDate,
												Expiry, 
												Strike,
												InitSpread,
												OptionType,
												pUnderlying,
												Rehauss,
												TriggerFreq,
												DiffCDO,
												IsCMSpread,
												CMSpreadMatu,
												pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_PropertyList(VARIANT attrNames,
													   VARIANT attrValues,
													   VARIANT attrTypes,
													   VARIANT *pRet)
{
	return theInstance.ARM_Credit_PropertyList(attrNames, 
												attrValues, 
												attrTypes,
												pRet);
}

STDMETHODIMP ARMCommonModule::ARMTMLeg(BSTR ptmIxType, 
									   double pstartDate, 
									   double pendDate, 
									   BSTR pPorR, 
									   double pspread,
									   BSTR ppayFrequency, 
									   BSTR presetFrequency, 
									   BSTR pinterestRule, 
									   BSTR pfwdRule, 
									   BSTR pstubRule, 
									   BSTR pccy,
									   BSTR *pRet)
{
		return theInstance.ARMTMLeg( ptmIxType, pstartDate, pendDate, pPorR, pspread, ppayFrequency, presetFrequency, pinterestRule, pfwdRule,pstubRule,pccy,pRet);
}
STDMETHODIMP ARMCommonModule::ARM_Credit_DefCurveIntensityPWC( double AsOfDate,
														VARIANT *pMatuRates, 
														VARIANT *pInputs, 
														double Type,
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	VolCurveId, 
														BSTR	calibrationAlgo,
														int lag,
														BSTR paramId,
														VARIANT *pRet)
{
	return theInstance.ARM_Credit_DefCurveIntensityPWC(AsOfDate,
													pMatuRates,
													pInputs,
													Type,
													Recovery,
													IRCurveId,
													bCurrency,
													bLabel,
													VolCurveId,
													calibrationAlgo,
													lag,
													paramId,
													pRet);

}

STDMETHODIMP ARMCommonModule::ARMGetNewFXVolFromSummit(BSTR pCcy1, 
													 BSTR pCcy2, 
													 double pDate, 
													 BSTR pCvName, 
													 BSTR pdomZcId, 
													 BSTR pforZcId,
													 double pfxSpot,
													 VARIANT *pForwards, 
													 BSTR pwhatIsInterpolated, 
													 double pcorrectSplineWithLinear,
													 BSTR pisATM,
													 BSTR *pRet)
{
		return theInstance.ARMGetNewFXVolFromSummit(pCcy1, pCcy2, pDate, pCvName, pdomZcId, pforZcId,pfxSpot,pForwards, pwhatIsInterpolated, pcorrectSplineWithLinear,	pisATM,pRet);
}



STDMETHODIMP ARMCommonModule::ARMGetInfoFromFxVolatility(BSTR pCurveId,	VARIANT *pRetMatData)
{
		return theInstance.ARMGetInfoFromFxVolatility(pCurveId, pRetMatData);
}


STDMETHODIMP ARMCommonModule::ARMGPModelParamCreate(BSTR pModelParamType, 
													 VARIANT *pParamTimes, 
													 VARIANT *pParamValues, 
													 BSTR pModelParamName, 
													 VARIANT *pLowerBoundaries, 
													 VARIANT *pUpperBoundaries, 
													 BSTR pInterpolMethod,
													 double pAdviseBreakPointTimes,
													 BSTR pCurrency,
													 BSTR *pRet)
{
		return theInstance.ARMGPModelParamCreate(pModelParamType, pParamTimes, pParamValues, pModelParamName, pLowerBoundaries,  pUpperBoundaries, pInterpolMethod,pAdviseBreakPointTimes,pCurrency, pRet);
}


STDMETHODIMP ARMCommonModule::ARMFXModelCreate(BSTR pZeroCurveId, 
											 VARIANT *pParamsId,
											 double pSpot,
											 BSTR pForCurveId,
											 long pModelType,
											 BSTR *pRet)
{
		return theInstance.ARMFXModelCreate(pZeroCurveId, pParamsId, pSpot, pForCurveId, pModelType,  pRet);
}


STDMETHODIMP ARMCommonModule::ARMcomputeHelvetix(double pNtl, 
										   double pAsOfDate, 
										   double pStartDate, 
										   double pEndDate, 
										   double pHelvetixType,
										   double pFee,
										   BSTR pCcy,
										   BSTR pIndexFund,
										   VARIANT pSpreadFund,
										   BSTR pDayCountFund, 
										   BSTR pPayFreqFund, 
										   BSTR pResetFreqFund,
										   BSTR pResetTimingFund,
										   BSTR pAdjFund, 
										   BSTR pDayCountAmort,
										   BSTR pIntRuleAmort, 
										   double pTxAmort,
										   BSTR pFreqAmort,
										   double pAmountAmort, 
										   BSTR pTypeAmort, 
										   BSTR pIndexPhasePre, 										   
										   double pSpreadPhasePre,
										   BSTR pDayCountPhasePre, 
										   BSTR pPayFreqPhasePre,
										   BSTR pResetFreqPhasePre,
										   BSTR pResetTimingPhasePre,
										   BSTR pAdjPhasePre, 
										   double pStartDatePhaseStruct,
										   double pEndDatePhaseStruct,
										   BSTR pResetTimingPhaseStruct,
										   BSTR pAdjPhaseStruct,										   
										   BSTR pDayCountPhaseStruct, 
										   BSTR pPayFreqPhaseStruct, 
										   BSTR pResetFreqPhaseStruct,
										   double pSpotPhaseStruct,
										   BSTR pOptionTypePhaseStruct,
										   BSTR pForeignCcyPhaseStruct,
										   BSTR pDomesticCcyPhaseStruct,
										   double pLevierPhaseStruct,
										   BSTR pStrikePhaseStruct, 
										   BSTR pPayOffDigitalPhaseStruct, 
										   BSTR pResetTimingOptionPhaseStruct,
										   BSTR pAdjOptionPhaseStruct,										   
										   BSTR pDayCountOptionPhaseStruct, 
										   BSTR pPayFreqOptionPhaseStruct, 
										   BSTR pResetFreqOptionPhaseStruct,
										   BSTR pPayTimingOptionPhaseStruct,
										   double pPayGapOptionPhaseStruct,
										   double pResetGapOptionPhaseStruct,								   
										   BSTR pIsCapped,
   										   double pTxCap, 
										   BSTR pIndexPhasePost,
										   double pSpreadPhasePost,
										   BSTR pDayCountPhasePost, 
										   BSTR pPayFreqPhasePost, 
										   BSTR pResetFreqPhasePost,
										   BSTR pResetTimingPhasePost,
										   BSTR pAdjPhasePost,
										   BSTR pRoll, 
										   BSTR pFwdRule,
										   BSTR pIntRule,
										   BSTR pStubRule,
										   VARIANT *pProba,
										   VARIANT *pMT,
										   BSTR pSmiledModEUR,
										   BSTR pSmiledModCHF,
										   BSTR pDeltaSmiledModEUR,
										   BSTR pMixtureModEURCHF,										   
										   BSTR pDeltaEURMixtureMod,
										   BSTR pBSEURMixtureMod,
										   BSTR pVolFxEURMixtureMod,
										   BSTR pDeltaCHFMixtureMod,
										   BSTR pBSCHFMixtureMod,										   
										   VARIANT *pRet)

										   
{
	return theInstance.ARMcomputeHelvetix(pNtl, 
										    pAsOfDate, 
											pStartDate, 
										    pEndDate, 
											pHelvetixType,
										    pFee,
										    pCcy,
										    pIndexFund,
										    pSpreadFund,
										    pDayCountFund, 
										    pPayFreqFund, 
										    pResetFreqFund,
										    pResetTimingFund,
										    pAdjFund, 
										    pDayCountAmort,
										    pIntRuleAmort, 
										    pTxAmort,
										    pFreqAmort,
										    pAmountAmort, 
										    pTypeAmort, 
										    pIndexPhasePre, 										   
										    pSpreadPhasePre,
										    pDayCountPhasePre, 
										    pPayFreqPhasePre,
										    pResetFreqPhasePre,
										    pResetTimingPhasePre,
										    pAdjPhasePre, 
										    pStartDatePhaseStruct,
										    pEndDatePhaseStruct,
										    pResetTimingPhaseStruct,
										    pAdjPhaseStruct,										   
										    pDayCountPhaseStruct, 
										    pPayFreqPhaseStruct, 
										    pResetFreqPhaseStruct,
										    pSpotPhaseStruct,
										    pOptionTypePhaseStruct,
										    pForeignCcyPhaseStruct,
										    pDomesticCcyPhaseStruct,
										    pLevierPhaseStruct,
										    pStrikePhaseStruct, 
											pPayOffDigitalPhaseStruct,
										    pResetTimingOptionPhaseStruct,
										    pAdjOptionPhaseStruct,										   
										    pDayCountOptionPhaseStruct, 
										    pPayFreqOptionPhaseStruct, 
										    pResetFreqOptionPhaseStruct,
										    pPayTimingOptionPhaseStruct,
										    pPayGapOptionPhaseStruct,
										    pResetGapOptionPhaseStruct,								   
										    pIsCapped,
   										    pTxCap, 
										    pIndexPhasePost,
										    pSpreadPhasePost,
										    pDayCountPhasePost, 
										    pPayFreqPhasePost, 
										    pResetFreqPhasePost,
										    pResetTimingPhasePost,
										    pAdjPhasePost,
										    pRoll, 
										    pFwdRule,
										    pIntRule,
										    pStubRule,
											pProba,
											pMT,
										    pSmiledModEUR,
											pSmiledModCHF,
										    pDeltaSmiledModEUR,
										    pMixtureModEURCHF,										   
										    pDeltaEURMixtureMod,
										    pBSEURMixtureMod,
										    pVolFxEURMixtureMod,
										    pDeltaCHFMixtureMod,
										    pBSCHFMixtureMod,									   
										    pRet);

}

STDMETHODIMP ARMCommonModule::ARMCapFloor(BSTR pSwapLegId, 
										BSTR pIsItCapOrFloor,
										VARIANT pStrike,
										BSTR *pRet)
{
		return theInstance.ARMCapFloor(pSwapLegId, pIsItCapOrFloor,pStrike,pRet);
}

STDMETHODIMP ARMCommonModule::ARMSpreadDigital(double pStartDate,
											 double pEndDate,
											 BSTR pCapOrFloor,
											 VARIANT pStrike,
											 VARIANT *pSpread,
											 VARIANT pPayOff,
											 BSTR pLiborType1,
											 BSTR pLiborType2,
											 VARIANT *pWeight,
											 BSTR pDayCount,
											 BSTR pResetFreq,
											 BSTR pPayFreq,
											 BSTR pResetTiming,
											 BSTR pPayTiming,
											 BSTR pCurrency,
											 double pResetGap,
											 BSTR pIntRule,
											 BSTR pStubRule,
											 VARIANT *pFixing1,
											 VARIANT *pFixing2,
											 BSTR *pRet)
{
		return theInstance.ARMSpreadDigital(pStartDate,
											 pEndDate,
											 pCapOrFloor,
											 pStrike,
											 pSpread,
											 pPayOff,
											 pLiborType1,
											 pLiborType2,
											 pWeight,
											 pDayCount,
											 pResetFreq,
											 pPayFreq,
											 pResetTiming,
											 pPayTiming,
											 pCurrency,
											 pResetGap,
											 pIntRule,
											 pStubRule,
											 pFixing1,
											 pFixing2,
											 pRet);
}


STDMETHODIMP ARMCommonModule::ARMSpreadDigitalFlt(double pStartDate,
											 double pEndDate,
											 BSTR pCapOrFloor,
											 VARIANT pStrike,
											 VARIANT *pSpread,
											 BSTR pPayOffLiborType,
											 BSTR pLiborType1,
											 BSTR pLiborType2,
											 VARIANT *pWeight,
											 BSTR pDayCount,
											 BSTR pResetFreq,
											 BSTR pPayFreq,
											 BSTR pResetTiming,
											 BSTR pPayTiming,
											 BSTR pCurrency,
											 double pResetGap,
											 BSTR pIntRule,
											 BSTR pStubRule,
											 VARIANT *pFixing1,
											 VARIANT *pFixing2,
											 BSTR *pRet)
{
		return theInstance.ARMSpreadDigitalFlt(pStartDate,
											 pEndDate,
											 pCapOrFloor,
											 pStrike,
											 pSpread,
											 pPayOffLiborType,
											 pLiborType1,
											 pLiborType2,
											 pWeight,
											 pDayCount,
											 pResetFreq,
											 pPayFreq,
											 pResetTiming,
											 pPayTiming,
											 pCurrency,
											 pResetGap,
											 pIntRule,
											 pStubRule,
											 pFixing1,
											 pFixing2,
											 pRet);
}

STDMETHODIMP ARMCommonModule::ARM_Credit_DefCurvePWC_ABS(double AsOfDate,
														VARIANT *Tenors, 
														VARIANT *Rates,
														VARIANT	*UpFront,
														VARIANT	*RefValIds,
														double  Recovery, 
														BSTR	IRCurveId,
														BSTR	bCurrency, 
														BSTR	bLabel, 
														BSTR	adjCalType,
														BSTR	IsSummitCurve,
														BSTR	VolCurveId, 
														BSTR	Accrued, //????
														BSTR	calibrationAlgo,
														BSTR	calibrationData,
														int lag,
														VARIANT *pRet)
{
	return theInstance.ARM_Credit_DefCurvePWC_ABS(AsOfDate,
													Tenors, 
													Rates,
													UpFront,
													RefValIds, 
													Recovery, 
													IRCurveId,
													bCurrency, 
													bLabel, 
													adjCalType,
													IsSummitCurve,
													VolCurveId, 
													Accrued, //????
													calibrationAlgo,
													calibrationData,
													lag,
													pRet);
}


