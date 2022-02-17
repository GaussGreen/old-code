// ARMModule.h: Definition of the ARMModule class
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARMMODULE_H__F153E6FC_B3F3_4A73_9A2D_6A0FE0A84A8F__INCLUDED_)
#define AFX_ARMMODULE_H__F153E6FC_B3F3_4A73_9A2D_6A0FE0A84A8F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "resource.h"       // main symbols
#include <CCString.h>

/////////////////////////////////////////////////////////////////////////////
// ARMModule
#include "local_DLLARM.h"

//ERROR TYPE
#define ARM_ERROR_IRCURVE 1
#define ARM_ERROR_DFCURVE 2
#define ARM_ERROR_DIFSZCURVE 3
#define ARM_ERROR_SZMATRIX 4
#define ARM_ERROR_PRICER 5
#define ARM_ERROR_FREQ 6
#define ARM_ERROR_DAYCOUNT 7
#define ARM_ERROR_ACCONDEF 8
#define ARM_ERROR_AMONDEF 9
#define ARM_ERROR_INTONDEF 10
#define ARM_ERROR_CORRMATRIX 11
#define ARM_ERROR_SECURITY 12
#define ARM_ERROR_MODEL 13

//void ERROR_MSG(CCString Err_mess,VARIANT* pRet,long noerr = 0);

class ARMCommonModule : 
	public IDispatchImpl<IARMModule, &IID_IARMModule, &LIBID_LOCAL_DLLARMLib>, 
	public ISupportErrorInfo,
	public CComObjectRoot,
	public CComCoClass<ARMCommonModule ,&CLSID_ARMModule>
{
public:
	ARMCommonModule () {}
BEGIN_COM_MAP(ARMCommonModule )
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY(IARMModule)
	COM_INTERFACE_ENTRY(ISupportErrorInfo)
END_COM_MAP()
//DECLARE_NOT_AGGREGATABLE(ARMModule) 
// Remove the comment from the line above if you don't want your object to 
// support aggregation. 

DECLARE_REGISTRY_RESOURCEID(IDR_ARMModule)
// ISupportsErrorInfo
	STDMETHOD(InterfaceSupportsErrorInfo)(REFIID riid);

// IARMModule
public:
	#include "..\ActiveX_ARM\ActiveXLib.h" 
	/**STDMETHOD(ARMAdjustToBusDate)(double pDate, BSTR pCalendrier, BSTR pRule, double *pDate2);
	STDMETHOD(ARMNextBusinessDay)(double pDate, BSTR pCalendrier, long pNbDays, double *pDate2);
	STDMETHOD(ARMFreeObject)(BSTR pId, long *pRet);
	STDMETHOD(ARMIsBusinessDay)(double pDate, BSTR pCalendrier, long *pRet);
	STDMETHOD(ARMGetZCFromSummit)(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate, BSTR *pRet);
	STDMETHOD(ARMFreeAllObjects)(long *pRet);
	STDMETHOD(ARMYcMod)(BSTR pZc, BSTR pZcDiscount, BSTR *pRet);
	STDMETHOD(ARMForwardYield)(BSTR pZc, double pMatu1, double pMatu2, BSTR pMeth, BSTR pAdj, double *pRet);
	STDMETHOD(ARMDiscountYield)(VARIANT *pZc, VARIANT *pMatu, VARIANT *pMeth, VARIANT *pRet);
	STDMETHOD(ARMLiborSwap)(VARIANT *pStartDate, VARIANT *pEndDate, VARIANT *pLiborType, VARIANT *pRecOrPay, VARIANT *pFixedRate, VARIANT *pSpread, VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMSwapPriceToRate)(VARIANT *pSwap, VARIANT *pDate, VARIANT *pPrice, VARIANT *pModel, VARIANT *pRet);
	STDMETHOD(ARMPrice)(VARIANT *pSec, VARIANT *pModel, VARIANT *pRet);
	STDMETHOD(ARMBetweenDates)(VARIANT *pDate1, VARIANT *pDate2, VARIANT *pDaycount, VARIANT *pIsYearFrac, VARIANT *pRet);
	STDMETHOD(ARMAddPeriod)(VARIANT *pDate, VARIANT *pFreq, VARIANT *pCcy, VARIANT *pNbPeriod, VARIANT *pAdjRule, VARIANT *pRet);
	STDMETHOD(ARMIsoCcy)(BSTR pCcy, BSTR pRefObj, BSTR *pRet);
	STDMETHOD(ARMGetSpotDays)(VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMGetLiborIndexDaycount)(VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMGetLiborTerm)(VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMGetFixedDayCount)(VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMGetFixedPayFreq)(VARIANT *pCcy, VARIANT *pRet);
	STDMETHOD(ARMComputeVolatility)(VARIANT *pVol, VARIANT* pMatu, VARIANT* pStrike, VARIANT* pTenor, VARIANT *pRet);
	STDMETHOD(ARMVolCurv)(VARIANT *pMatu,VARIANT* pStrikes,VARIANT* pVols,double pAsOf,BSTR pStkType,BSTR pVolType,BSTR pCcy,BSTR *pRet);
	STDMETHOD(ARMGetVolCubeFromSummit)(VARIANT *pIndex, VARIANT* pCcy, VARIANT* pCvName, VARIANT* pAsOf, VARIANT* pType, VARIANT* pSmiles, VARIANT* pTypeCube, VARIANT *pRet);
	STDMETHOD(ARMParallelShift)(BSTR pZc, double pBump, BSTR *pRet);
    STDMETHOD(ARMBumpVolatility)(BSTR pVol, double pValue, long pNthLine, long pNthCol,BSTR pIsCumul, BSTR *pRet);
	STDMETHOD(ARMBsSmiledModel)(double pDate, double pSpot, BSTR pDividend, BSTR pDiscrate, BSTR pVolATM, BSTR pRo, BSTR pNu, BSTR pIsSABR, BSTR pBeta, BSTR *pRet);
	STDMETHOD(ARM_FxConvert)(VARIANT *pccy1,VARIANT *pccy2,VARIANT *pAsOfDate,VARIANT *pCvName,VARIANT *pRet);
	STDMETHOD(ARM_DiscountPrice)(VARIANT *pCurve,VARIANT *pMatu,VARIANT *pRet);
	STDMETHOD(ARM_View)(VARIANT *pObjet,VARIANT *pRet);
	STDMETHOD(ARMGetVolFromSummit)(VARIANT *pIndex, VARIANT* pCcy, VARIANT* pCvName, VARIANT* pAsOf, VARIANT* pType, VARIANT* pMatIndex, VARIANT* pImpOrHist, VARIANT *pRet);
	STDMETHOD(ARMGetFXVolFromSummit)(BSTR pCcy1, BSTR pCcy2, double pDate, BSTR pCvName, BSTR pType, BSTR *pRet);
	STDMETHOD(ARMGetFXCorrelFromSummit)(BSTR pCcy1, BSTR pIndex, BSTR pCcy2, double pDate, BSTR pCvName, VARIANT* pTenors, BSTR *pRet);
	STDMETHOD(ARMGetCorrelFromSummit)(BSTR pCcy1, BSTR pIndex1, BSTR pCcy2, BSTR pIndex2, double pDate, BSTR pCvName, BSTR *pRet);
	STDMETHOD(ARMcptBonifixEUR)(double pDate, BSTR pDepart, BSTR pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pTiming, double pBarriere, double pSpdPostBar, double pMarge, double pFunding, double pSpd2phase, double pSoulte, BSTR pBsModId1, BSTR pBsModId2, BSTR pBsModId3, BSTR pProbaMarge, double pInt, VARIANT *pRet);
	STDMETHOD(ARMcptBonifixUSD)(double pDate, BSTR pDepart, BSTR pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pTiming, double pBarriere, double pSpdPostBar, double pMarge, double pFunding, BSTR pFundingFreq, double pSpd2phase, double pSoulte, BSTR pYcModId, BSTR pBsModId1, BSTR pBsModId2, BSTR pProbaMarge, double pInt, VARIANT *pRet);
	STDMETHOD(ARMcptOptilix)(double pDate, BSTR pDepart, BSTR pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pSjOptibor, double pGapPostFixe, double pMarge, double pFunding, double pSpd2phase, double pSoulte, BSTR pBsModId1, BSTR pBsModId2, BSTR pBsModId3, double pCap, double pInt, VARIANT *pRet);
	STDMETHOD(ARMSetEtoolkit)(BSTR pUserName,BSTR pPassWord,BSTR pDatabaseContext,BSTR pItConfigDomainDir,BSTR pItDomainName,long *pRet);
	STDMETHOD(ARMConnectionEtoolkit)(VARIANT *pRet);
	STDMETHOD(ARMVolFlat)(double pVol,double pDate,BSTR pCcy,BSTR *pRet);
	STDMETHOD(ARMVolCube)(BSTR pATMVol,VARIANT *pSmileCurveIds,VARIANT *pTenors,BSTR pVolType,BSTR pRefObj,BSTR *pRet);
	STDMETHOD(ARMDeconnectionEtoolkit)(VARIANT *pRet);
	STDMETHOD(ARMZcFlat)(double pZc,double pDate,BSTR pCcy,BSTR *pRet);
	STDMETHOD(ARMBsModel)(double pDate,double pSpot,BSTR pDividend,BSTR pDiscrate,BSTR pVol,BSTR pTypeStk,BSTR *pRet);
	STDMETHOD(ARMSwitchToETK)();
	STDMETHOD(ARMSwitchToWSETK)();
	STDMETHOD(ARMInfocentreConnect)();
	STDMETHOD(ARMSwitchToFLATFILE)();
	STDMETHOD(ARMZCLINT)(VARIANT* pMatu, VARIANT* pRate, BSTR pMeth, double pDate, BSTR pCurrency, BSTR pInterpMeth, BSTR *pRet);
	STDMETHOD(ARM_zcspreaded)(BSTR zcSprId,BSTR zcInitId,double date,BSTR MMFreq,BSTR SwapFreq,BSTR ccyId,BSTR *pRet);
	STDMETHOD(ARMCreateZCSwapInt)(double pDate,VARIANT *pMatu,VARIANT *pRate,BSTR pMMVsFut,BSTR pSwapVsFut,BSTR pRaw,BSTR pInterp,BSTR pCcy,BSTR pRefObj,BSTR *pRet);
	STDMETHOD(ARMGetInitialCurveFromSummit)(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate, BSTR pAdjOrNot, VARIANT *pRetMat,VARIANT *pRetRate);
	STDMETHOD(ARMTHREEMONTHFUT)(BSTR pDelivery, long pMarket,BSTR pCcy,BSTR pRefObj,BSTR *pRet);
	STDMETHOD(ARMFutPibor)( BSTR pDelivery,BSTR *pRet);
	STDMETHOD(ARMIRFUT)(double pDelivery,BSTR pIdUnderlying,BSTR pRefObj,BSTR *pRet);
	STDMETHOD(ARMLibor)(BSTR pLiborTypeId,BSTR pCcyId, BSTR pResetFreqId,BSTR pPayFreqId,BSTR pRefObj,BSTR pBasis,BSTR pIntRule, BSTR *pRet);
	STDMETHOD(ARMLiborSwaption)(double pStartDate, double pEndDate, BSTR pReceiveOrPay, double pStrike, double pMaturity, BSTR pLiborType, double pSpread, BSTR pExerciseType, BSTR pResetFreq, BSTR pPayFreq, BSTR pCcyId, BSTR pRefObj, BSTR *pRet);
	STDMETHOD(ARMFixedLeg)(double pStartDate, double pEndDate, BSTR pReceiveOrPay, double pFixRate, BSTR pDayCount, BSTR pFreq, BSTR pDecompFreq, BSTR pPayTiming, BSTR pIntRule, BSTR pStubRule, BSTR pCcyId, BSTR pPayCalName, BSTR pNxChange,BSTR pRefObj,BSTR *pRet);
	STDMETHOD(ARMcptOptilixCMS)(double pDate, BSTR pDepart, BSTR pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pSjOptibor, BSTR pResetTiming, double pMarge, double pFunding, BSTR pFundingFreq, double pSpd2phase, double pSoulte, BSTR pBsModId1, double pInt, VARIANT *pRet);
	STDMETHOD(ARMcptBilibor)(double pDate, BSTR pDepart, BSTR pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pSjCHF, BSTR pResetTiming, double pMarge, double pFunding, BSTR pFundingFreq, double pSpd2phase, double pSoulte, BSTR pBsModId, BSTR pBsModShiftEURId, BSTR pBsModShiftCHFId, BSTR pBsModCorrId, double pInt, VARIANT *pRet);
	STDMETHOD(ARMBsSlModel)(double pDate,BSTR pZc,BSTR pVolSpreadLock,BSTR pCvCapVol,BSTR pCvIndexVol,BSTR *pRet);
	STDMETHOD(ARMGlobDFBS)(BSTR pDomBSId,BSTR pDomCurrId,BSTR pFrgBSId,BSTR pFrgCurrId,BSTR pFxVolCrvId,BSTR pFFxCorrId,BSTR pRatesCorrId,BSTR *pRet);
	STDMETHOD(ARMDFFXBS)(BSTR pDVolId,BSTR pFVolId,BSTR pDZcId,BSTR pFZcId,BSTR pDFxCorrId,BSTR pFFxCorrId,BSTR pFxVolId,double pRatesCorr,BSTR *pRet);
	STDMETHOD(ARMTRIBSMODEL)(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,BSTR *pRet);
	STDMETHOD(ARMTRIBSDUAL)(BSTR pModel1,BSTR pModel2,BSTR pDiscModel,BSTR pFX1DiscVol,BSTR pFX2DiscVol,BSTR pIdx1Idx2Corr,BSTR pIdx1DiscIdxCorr,BSTR pIdx2DiscIdxCorr,BSTR pIdx1FxCorr,BSTR pIdx2FxCorr,int pQuantoFlag,double pCorrelForAdj,int pWithslopeflag,BSTR *pRet);
	STDMETHOD(ARMcptBonibor)(double pDate,BSTR pDepart,BSTR pMatStruct,BSTR pMatTot,BSTR pAmort,BSTR pFreq,BSTR pSjUSD,BSTR pTiming,double pBarriere,double pSpdPostBar,double pMarge,double pFunding,BSTR pFundingFreq,double pSpd2phase,double pSoulte,BSTR pYcModId,BSTR pBsModId,BSTR pBsModVolUSDId,BSTR pBsModCorrPlusId,BSTR pBsModCorrMoinsId,BSTR pCrossModId,BSTR pProbaMarge,double pInt,VARIANT *pRet);
	STDMETHOD(ARMAswPrice)(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pMargin, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet);
	STDMETHOD(ARMAswMargin)(double pMaturity, double pCpn, BSTR pFreq, BSTR pBase, double pPrice, double pRedemptionPrice, double pAsOf, double pDelivery, BSTR pFixDecompfreq, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, BSTR pAmortizationId, long pSolve, double pMinValue, double pMaxValue, double *pRet);
	STDMETHOD(ARMFrnPrice)(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pValoMargin, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet);
	STDMETHOD(ARMFrnMargin)(double pAsOf, double pDelivery, double pMaturity, BSTR pCcy1, BSTR pIndex1, BSTR pFwdCurve1, BSTR pDiscCurve1, double pFacialMargin, double pPrice, BSTR pCcy2, BSTR pIndex2, BSTR pFwdCurve2, BSTR pDiscCurve2, double pFixing, double pSpread, double pOutMode, long pSolve, BSTR pAmortizationId, double *pRet);
	STDMETHOD(ARMcptBonifixProgEUR)(double pDate, BSTR pDepart, VARIANT* pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pTiming, VARIANT* pBarriere, double pSpdPostBar, double pMarge, double pFunding, BSTR pFundingFreq, double pSpd2phase, double pSoulte, BSTR pBsModId1, BSTR pBsModId2, BSTR pBsModId3, BSTR pProbaMarge, double pInt, VARIANT *pRet);
	STDMETHOD(ARMcptBonifixProgUSD)(double pDate, BSTR pDepart, VARIANT* pMatStruct, BSTR pMatTot, BSTR pAmort, BSTR pFreq, BSTR pTiming, VARIANT* pBarriere, double pSpdPostBar, double pMarge, double pFunding, BSTR pFundingFreq, double pSpd2phase, double pSoulte, BSTR pYcModId, BSTR pBsModId1, BSTR pBsModId2, BSTR pProbaMarge, double pInt, VARIANT *pRet);
	STDMETHOD(ARMcptDigital)(double pAsOf,double pStartDate, double pEndDate, double pNtl, BSTR pIndex, BSTR pBsmod,BSTR pBsmodDelta,BSTR pBsmodVega, BSTR pFreqP, BSTR pResetTiming, BSTR pPorR, BSTR pCcy,  BSTR pCcyIdx, BSTR pDayCount, BSTR pCapOrFloor, BSTR pAmort, BSTR pStrike, BSTR pPayOff, BSTR pSpd, double pResetGap, double pSpreadBelow, double pSpreadAbove, BSTR pFwdRule, BSTR pIntRule,BSTR pStubRule,BSTR pFreqAmort, double pTxAmort, double pAmountAmort, double pRefDate, VARIANT *pRet);
	STDMETHOD(ARMRefValue)(VARIANT *pdates,VARIANT* pvalues,VARIANT* pvalues2,long valueType,long conversion,BSTR calcMethod,BSTR *pRet);
	STDMETHOD(ARMCreateGenCorrelManager)(VARIANT *pMktTags,VARIANT *pIntraMktTags,VARIANT *pCorrelCurveIds,BSTR *pRet);
	STDMETHOD(ARMBSConvAdjust)(BSTR pSUMMITFormulaeUsed,BSTR *pRet);
	STDMETHOD(ARMBsModelGen)(BSTR pYieldCurve,BSTR pVolatility,BSTR pCorrMgr,BSTR pCnvxManager,BSTR pCapletVol,BSTR pSpreadLock,BSTR pDiscCurve,BSTR pCorrel,BSTR pCashVol,BSTR pSpreadVol,BSTR pModelType,BSTR pSpreadVolType,BSTR pSabrMod, BSTR *pRet);
	STDMETHOD(ARMcptQTF)(double pAsOf,double pStartDate, double pEndDate, double pNtl, BSTR pIndex, BSTR pBsmod, BSTR pFreqP, BSTR pResetTiming, BSTR pCcy, BSTR pDayCount, BSTR pAmort, BSTR pTauxBoni, BSTR pBarrier, BSTR pSpdTf, double pResetGap, double pSpreadBelow, double pSpreadAbove, BSTR pFwdRule, BSTR pIntRule,BSTR pStubRule,BSTR pFreqAmort, double pTxAmort, double pAmountAmort, double pRefDate, VARIANT *pRet);
	STDMETHOD(ARMcptSPTQTF)(double pAsOf,double pStartDate,double pDateSecondPhase, double pEndDate, double pNtl, BSTR pIndex, BSTR pIndexFund,BSTR pIndexSF, BSTR pBsmod,BSTR pBsmodDeltaCcy1,BSTR pBsmodVegaCcy1,BSTR pBsmodDeltaCcy2,BSTR pBsmodVegaCcy2,BSTR pBsmodFxCorrel, BSTR pFreqP,BSTR pFreqPFund,BSTR pFreqPSF,BSTR pFreqR, BSTR pResetTiming, BSTR pResetTimingSF, BSTR pCcy, BSTR pCcyIdx, BSTR pDayCount, BSTR pDayCountSF,double pFee, BSTR pBarrier, BSTR pSpdTf,BSTR pSpd,BSTR pSpdSF, BSTR pSpdSFfund, double pResetGap, double pResetGapSF, BSTR pFwdRule, BSTR pIntRule,BSTR pStubRule,BSTR pAmort,BSTR pFreqAmort, double pTxAmort, double pAmountAmort, double pRefDate, VARIANT *pRet);
	STDMETHOD(ARMDisplaySchedule)(BSTR pLegId,BSTR pDataType,VARIANT *pRet);
	STDMETHOD(ARMIrIndex)(BSTR pDaycount,BSTR pPayFreq,double pMaturity,BSTR pCompMethod,BSTR pFwdRule,BSTR pResetTiming,double pResetGap,BSTR pPayTiming,double pPayGap,BSTR pCcy,BSTR pIndexType,double pDecompFreq,BSTR pIntRule,BSTR pResetFreq, BSTR *pRet);
	STDMETHOD(ARMSwapleg)(BSTR pIndexId, double pStartDate, double pEndDate, BSTR pRecOrPay, VARIANT pSpread, BSTR pCcy, BSTR pDayCount, double pResetGap, BSTR pResetCal, BSTR pPayCal, double pDecompPricingFlag, BSTR pNxChange, BSTR pStubRule, double pRefDate, BSTR pAdjStartDate, BSTR *pRet);
	STDMETHOD(ARMConstRefvalue)(double pValue, BSTR *pRet);
	STDMETHOD(ARMShutdownETK)();
	STDMETHOD(ARMGetInitialVolFromSummit)(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate, BSTR pType, BSTR pMatIndex, VARIANT *pRetMat, VARIANT *pRetTenor, VARIANT *pRetVol);
	STDMETHOD(ARMBond)(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedempPrice,double pPeriodicity,VARIANT pDaycount,double pSettleGap,double pCpnDateFlag,BSTR pCcy,BSTR *pRet);
	STDMETHOD(ARMINFCreateOATLeg)(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,BSTR *pRet);
	STDMETHOD(ARMSwap)(BSTR pSswapleg1,BSTR pSswapleg2,double pMinPay,BSTR *pRet);
	STDMETHOD(ARMPToYield)(BSTR pBond,double pSettleDate,double pPrice,VARIANT *pRet);
	STDMETHOD(ARMYToPrice)(BSTR pBond,double pSettleDate,double pYield,VARIANT *pRet);
	STDMETHOD(ARMYToDuration)(BSTR pBond,double pSettleDate,double pActuRate,double pFlagCpn,VARIANT *pRet);
	STDMETHOD(ARMLiborleg)(double pStartDate,double pEndDate,BSTR pLiborType,BSTR pRecOrPay,VARIANT pSpread,BSTR pResetFReq,BSTR pPayFreq,BSTR pResetTiming,BSTR pPayTiming,BSTR pCcy,BSTR pIntRule,double pResetGap,BSTR pResetCal,BSTR pPayCal,double pDecompPricingFlag,BSTR pNxChange,BSTR pStubRule,double pRefDate,BSTR pAdjStartDate,BSTR *pRet);
	STDMETHOD(ARMImpliedSpread)(BSTR pSwap,BSTR pModel,double pPrice,double pLeg1or2,double *pRet);
	STDMETHOD(ARMDiscountPrice)(BSTR pZeroCurve,double pMatu,double *pRet);
	STDMETHOD(ARMINFCreateCurve)(double pAsOf,BSTR pIndexName,double pCPIIndexValue,double pCPIIndexDate,VARIANT* pMatu,VARIANT* pRate,BSTR pMonthlyInterpType,BSTR pDailyInterpType,BSTR pDCFMonthly,BSTR pDCFDaily,BSTR pExtrapolType,BSTR pResetManager,BSTR pSeasonManager,BSTR *pRet);
	STDMETHOD(ARMINFInterpCPI)(BSTR pZc,double pCPIDate,BSTR pDCFlag,BSTR pDailyInterpType,BSTR pCPIlag,double pWeight,double *pRet);
	STDMETHOD(ARMINFSeasonManager)(VARIANT* pMonthList,VARIANT* pValues,BSTR pSeasonAdjMode,BSTR *pRet);
	STDMETHOD(ARMINFResetManager)(VARIANT* pDatas,double pNbIndex,BSTR *pRet);
	STDMETHOD(ARMINFYcMod)(BSTR pYieldCurve,BSTR pInfCurve,BSTR *pRet);
	STDMETHOD(ARMBaseReplicationConnect)();
	STDMETHOD(ARMGetInitialFXVolFromSummit)(BSTR pCcy1,BSTR pCcy2,double pDate,BSTR pCvName,BSTR pImpOrHist,BSTR pVolType,VARIANT *pRetMat,VARIANT *pRetTenor,VARIANT *pRetVol);
	STDMETHOD(ARMCreateZCFromSummit)(BSTR pIndex, BSTR pCurrency, BSTR pCvName, double pDate,BSTR pAdj,BSTR pRaw,BSTR *pRet);
	STDMETHOD(ARMBumpCurve)(BSTR pZc,double pEpsilon,long pMethod,BSTR pPlot,BSTR *pRet);
	STDMETHOD(ARMAccrued)(BSTR pSec,double pDate,BSTR pModel,double *pRet);
	STDMETHOD(ARMcomputeBilibor)(double pAsOf,double pStartDate,double pDateSecondPhase,double pEndDate,double pNtl,BSTR pIndex,BSTR pIndexFund,BSTR pIndexSF,BSTR pBsmod,BSTR pBsmodFund,BSTR pBsmodDeltaCcy1, BSTR pBsmodDeltaFund,BSTR pBsmodDeltaCcy2,BSTR pBsmodFxCorrel,BSTR pFreqP,BSTR pFreqR,BSTR pFreqPFund,BSTR pFreqRFund,BSTR pFreqPSF,BSTR pFreqRSF, BSTR pResetTiming,BSTR pResetTimingSF,BSTR pCcy1,BSTR pCcy2,BSTR pDayCount,BSTR pDayCountSF,double pSpdPF,double pSpdSF,double pSpdfund, double pSpdfund2,double pResetGap,double pResetGapSF,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,BSTR pAmort,BSTR pFreqAmort,double pTxAmort,double pAmountAmort, double pRefDate, double pFee, VARIANT *pRet);
	STDMETHOD(ARMcomputeOptilix)(double pAsOf,double pStartDate,double pDateSecondPhase,double pEndDate,double pNtl,BSTR pIndex,BSTR pIndexFund,BSTR pIndexSF,BSTR pBsmod,BSTR pBsmodFund,BSTR pBsmodDeltaCcy, BSTR pBsmodDeltaFund,BSTR pFreqP,BSTR pFreqR,BSTR pFreqPFund,BSTR pFreqRFund,BSTR pFreqPSF,BSTR pFreqRSF, BSTR pResetTiming,BSTR pResetTimingSF,BSTR pCcy,BSTR pDayCount,BSTR pDayCountSF, double pSpdSF,VARIANT pSpdfund,double pResetGap,double pResetGapSF,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,BSTR pAmort,BSTR pFreqAmort,double pTxAmort,double pAmountAmort, double pRefDate, double pFee, VARIANT *pRet);
	STDMETHOD(ARMClonedAndSetNotional)(BSTR bLegId,BSTR bAmortId,BSTR *pRet);
	STDMETHOD(ARM_INF_GetZcFromSummit)(BSTR Index,BSTR Ccy,BSTR cvname,double date,BSTR seasonAdj,BSTR seasonAdjMode,BSTR *pRet);
	STDMETHOD(ARM_INF_CreateGenericLeg)(double pStartDate,double pEndDate,BSTR pInfIdx,BSTR pRcvOrPay,BSTR pInterpType,double pLeverage,double pSpread,BSTR pResetFreq,BSTR pDaycount,BSTR pResetCal,BSTR pFwdRule,BSTR pIntRule,BSTR pStubRule,double pResetNumGap,double pResetDenomGap,BSTR pPayFreq,double pPayGap,BSTR pPayCal,BSTR pFinalNotionalType,double pFirstReset,BSTR *pRet);
	STDMETHOD(ARMRiskyBond)(double pIssueDate,double pMaturityDate,double pFirstCpnDate,double pCpnRate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet);
	STDMETHOD(ARMRiskyBondWithCF)(double pAsOfDate,double pRedemptionPrice,long pPeriodicity,VARIANT pDaycount,VARIANT *pYearTerms,VARIANT *pCashFlows,long pSettleGap,long pCpnDateFlag,BSTR pCcyId,double pRepo,double pSsl,double pRecoveryRate,BSTR *pRet);
	STDMETHOD(ARMHyperCube)(VARIANT *pVolCurvId,VARIANT* pKeys,BSTR *pRet);
	STDMETHOD(ARMcomputePentifix)(double pNtl,double pStartdatePhase1,BSTR pCcy,BSTR pIndexPhase1,double pSpreadPhase1,BSTR pDayCountPhase1,BSTR pPayFreqPhase1,BSTR pResetFreqPhase1,BSTR pResetTimingPhase1,BSTR pRoll,BSTR pAdjPhase1,BSTR pStub,BSTR pIndexPhase2DIG,BSTR pIndexLongPhase2DIG,BSTR pStrikePhase2DIG,BSTR pResetTimingPhase2DIG,BSTR pAdjPhase2DIG,double pStartDatePhase2,double pSpreadPhase2,BSTR pDayCountPhase2,BSTR pPayFreqPhase2,BSTR pResetFreqPhase2,BSTR pAdjPhase2,double pStartDatePhase3,double pEndDatePhase3,BSTR pIndexPhase3,double pSpreadPhase3,BSTR pDayCountPhase3,BSTR pPayFreqPhase3,BSTR pResetFreqPhase3,BSTR pResetTimingPhase3,BSTR pAdjPhase3,BSTR pIndexFund,VARIANT pSpreadFund,BSTR pDayCountFund,BSTR pPayFreqFund,BSTR pResetFreqFund,BSTR pResetTimingFund,BSTR pAdjFund, double pEndDateAmort, BSTR pDayCountAmort,BSTR pIntRuleAmort,double pTxAmort,BSTR pFreqAmort,double pAmountAmort,BSTR pTypeAmort,BSTR pFloorOrcap,double pFee,BSTR pVolCurvFromMatriceShift,BSTR pVol,BSTR pVolCub,BSTR pCorrManager,BSTR pConvexityManager,BSTR pZc,BSTR pSmiledMod,BSTR pSmiledModBump, BSTR pHyperCubeCorrel,VARIANT *pRet);
	STDMETHOD(ARMCmsLeg)(double startDate,double endDate,BSTR cmsTypeId,BSTR receiveOrPay,BSTR yieldDecompFreq,BSTR swapLegDayCount,BSTR resetFreq,BSTR payFreq,long resetGap,BSTR intRule,BSTR ccyName,BSTR resetTiming,BSTR stubRule,BSTR adjStartDate,BSTR *pRet);
	STDMETHOD(ARM_ReplicConvAdjust_Create)(BSTR Payoff_ReplicMode, double Payoff_StepOrReplicPrecision, BSTR Payoff_StopMode, double Payoff_StopThreshold, BSTR Sensi_ReplicMode, double Sensi_StepOrReplicPrecision, BSTR Sensi_StopMode, double Sensi_StopThreshold, BSTR UsedModelId, double StrikeMinReplic, double StrikeMaxReplic, BSTR *pRet);
	STDMETHOD(ARM_MapConvAdjust_Create)(BSTR LiborArrearAdj,BSTR NaturalCMSAdj,BSTR PaymentLagAdj,BSTR *pRet);
	STDMETHOD(ARMcomputePentibor)(double pNtl,double pStartdatePhase1,BSTR pCcy,BSTR pIndexPay, BSTR pIndexPhase1,double pSpreadPhase1,BSTR pDayCountPhase1,BSTR pPayFreqPhase1,BSTR pResetFreqPhase1,BSTR pResetTimingPhase1,BSTR pRoll,BSTR pAdjPhase1,BSTR pStub,BSTR pIndexPhase2DIG,BSTR pIndexLongPhase2DIG,BSTR pStrikePhase2DIG,BSTR pResetTimingPhase2DIG,BSTR pAdjPhase2DIG,double pStartDatePhase2,double pSpreadPhase2,BSTR pDayCountPhase2,BSTR pPayFreqPhase2,BSTR pResetFreqPhase2,BSTR pAdjPhase2,double pStartDatePhase3,double pEndDatePhase3,BSTR pIndexPhase3,double pSpreadPhase3,BSTR pDayCountPhase3,BSTR pPayFreqPhase3,BSTR pResetFreqPhase3,BSTR pResetTimingPhase3,BSTR pAdjPhase3,BSTR pIndexFund,VARIANT pSpreadFund,BSTR pDayCountFund,BSTR pPayFreqFund,BSTR pResetFreqFund,BSTR pResetTimingFund,BSTR pAdjFund,double pEndDateAmort, BSTR pDayCountAmort,BSTR pIntRuleAmort,double pTxAmort,BSTR pFreqAmort,double pAmountAmort,BSTR pTypeAmort,double pFee,BSTR pVolCurvFromMatriceShift,BSTR pVol,BSTR pVolCub,BSTR pConvexityManager,BSTR pZc,BSTR pSmiledMod,BSTR pSmiledModBump, BSTR pHyperCubeCorrel,BSTR pIndexIndexCorrelCube, BSTR pCorrEUR, BSTR pInterCorr, VARIANT *pRet);
	STDMETHOD(ARMIndexIndexCorrelCube)(VARIANT *pVolCurvId, VARIANT* pTenors1List, VARIANT* pTenors2List, BSTR pInterSurfInterp, BSTR *pRet);
	STDMETHOD(ARMCreateGenCorrelatorManager)(VARIANT *pMktTags,VARIANT *pHyperDiagVol, VARIANT *pIndexIndexVol, VARIANT* pCorrelVol, VARIANT* pIndexVol,BSTR *pRet);
	STDMETHOD(ARMcomputePentilix)(double pNtl,double pStartdatePhase1,BSTR pCcy,BSTR pIndexPay, BSTR pIndexPhase1,double pSpreadPhase1,BSTR pDayCountPhase1,BSTR pPayFreqPhase1,BSTR pResetFreqPhase1,BSTR pResetTimingPhase1,BSTR pRoll,BSTR pAdjPhase1,BSTR pStub,BSTR pIndexPhase2DIG,BSTR pIndexLongPhase2DIG,BSTR pStrikePhase2DIG,BSTR pResetTimingPhase2DIG,BSTR pAdjPhase2DIG,double pStartDatePhase2,double pSpreadPhase2,BSTR pDayCountPhase2,BSTR pPayFreqPhase2,BSTR pResetFreqPhase2,BSTR pAdjPhase2,double pStartDatePhase3,double pEndDatePhase3,BSTR pIndexPhase3,double pSpreadPhase3,BSTR pDayCountPhase3,BSTR pPayFreqPhase3,BSTR pResetFreqPhase3,BSTR pResetTimingPhase3,BSTR pAdjPhase3,BSTR pIndexFund,VARIANT pSpreadFund,BSTR pDayCountFund,BSTR pPayFreqFund,BSTR pResetFreqFund,BSTR pResetTimingFund,BSTR pAdjFund,double pEndDateAmort, BSTR pDayCountAmort,BSTR pIntRuleAmort,double pTxAmort,BSTR pFreqAmort,double pAmountAmort,BSTR pTypeAmort,double pFee,BSTR pVolCurvFromMatriceShift,BSTR pVol,BSTR pVolCub,BSTR pConvexityManager,BSTR pZc,BSTR pSmiledMod,BSTR pSmiledModBump, BSTR pHyperCubeCorrel,BSTR pIndexIndexCorrelCube, BSTR pCorrEUR, BSTR pInterCorr, VARIANT *pRet);

	//Section Credit
	STDMETHOD(ARM_Credit_Delivery)(VARIANT *pAsOfDate,VARIANT *pTenorContract,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CptInterpolDefCurve)(BSTR pCurve, BSTR pTenor, double pSlope, double pDate, double pInterpDate, VARIANT *pRet);
	STDMETHOD(ARM_Credit_DefaultProba)(VARIANT *pCurve,VARIANT *pMatu ,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetBeta )(VARIANT *pPricer,VARIANT *pLabel,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetBaseCorrelation )(VARIANT *pPricer,double pEquityAmount, VARIANT *pMktSpreads, double pSmileType,VARIANT *pSeeds,VARIANT *pUpfronts, VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetImpliedCorrelation )(VARIANT *pPricer,double pRange, VARIANT *pMktSpreads, double pSmileType,VARIANT *pSeeds,VARIANT *pUpfronts, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Mezzanine) (double pEffectiveDate,double pEndDate,double pSpread,double pMezzAmount,double pSubAmount,VARIANT *pLabels,VARIANT *pNotionals,BSTR pFreqFeeLeg,BSTR pDayCountFrq,double pFirst_period_refdate,BSTR pAccruedOnDefault,BSTR Currency,double pPayCreditLag,BSTR pStub,BSTR pFreqDefLeg,double pBinary,BSTR pPayCal,BSTR LongOrShortRisk,double  TradedNotional,BSTR IncludeMatu,double pFstCpnEffDate,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CDO2)(double pEffectiveDate,double pEndDate,BSTR pPortfolio,double pSpread,double pSubAmount,double pMezzAmount,BSTR pFreqFeeLeg,BSTR pFreqDefLeg,BSTR pDayCountFrq,double pFirst_period_refdate,BSTR pAccruedOnDefault,BSTR Currency,double pPayCreditLag,BSTR pStub,double pBinary,BSTR pPayCal,BSTR LongOrShortRisk,double  TradedNotional,BSTR CrossSub,BSTR IncludeMatu,double pFstCpnEffDate,VARIANT *pRet);
	STDMETHOD(ARM_Credit_ModelMultiCurves) (VARIANT *pIRcurve,VARIANT *pDefCurves,VARIANT *pRecovery, BSTR CorrelationId,BSTR	pVolcurve,BSTR	pCpnInfcurve, BSTR	pCpnIRcurve,VARIANT *pRet);
	STDMETHOD(ARM_Credit_FTD) (VARIANT *pEffectiveDate,VARIANT *pEndDate,VARIANT *pSpread,VARIANT *pLabels ,VARIANT *pFixingFreq,VARIANT *pDayCountFrq,VARIANT *pFirst_period_refdate,VARIANT *pIssuerNotional,VARIANT *pAccruedOnDefault,VARIANT *pCurrency,VARIANT *pPayCreditLag,VARIANT *pStub,double pFstCpnEffDate,VARIANT *pRet);
	STDMETHOD(ARM_Credit_NTD) (double pEffectiveDate,double pEndDate,double pSpread,int pFirstNumDefault,int pLastNumDefault,VARIANT *pLabels,BSTR pFixingFreq,BSTR pDayCountFrq,double pFirst_period_refdate,double	pIssuerNotional,BSTR pAccruedOnDefault,BSTR pCurrency,double pPayCreditLag,BSTR	pStub,BSTR pFrequencyDefLeg,double pBinary,BSTR	pPayCal,BSTR pRcvFeeLeg,double TradedNotional,BSTR IncludeMatu,double pFstCpnEffDate, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Price) (VARIANT *pPricer ,VARIANT *pAsOfDate,VARIANT *pRet);
	STDMETHOD(ARM_Credit_RiskyDuration) (VARIANT *pDefCurve,double Date,BSTR Tenor,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CDONPV) (VARIANT *pPricer ,VARIANT *pCPTTYPE,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CorrMatrix) (VARIANT *pLabels,VARIANT *pCoefs,double AsOf,BSTR Name,VARIANT *pRet);
	STDMETHOD(ARM_Credit_ExtractCorrMatrix) (VARIANT *pCorrMatrixId,VARIANT *pLabels,VARIANT *pRet);
	STDMETHOD(ARM_Credit_Spread)(VARIANT *pPricer,VARIANT *pMTM,VARIANT *pRet);
	STDMETHOD(ARM_Credit_SetLabel)(VARIANT *pCurve,VARIANT *pLabel,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetLabel)(VARIANT *pCurve,VARIANT *pRet);
	STDMETHOD(ARM_Credit_Sensitivity)(VARIANT *pPricer,VARIANT *pType,VARIANT *pPlot,VARIANT *pLabel,VARIANT *pEpsilon,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetCleanSpreadTranche)(VARIANT *pPricer,VARIANT *pPlot,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetDefProbTranche)(BSTR PricerId,double pPlot,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetDuration)(VARIANT *pPricer,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GenSchedule)(VARIANT *pAccStartDate,VARIANT *pAccEndDate,VARIANT *pFixingFreq,VARIANT *pDayCountFrq,VARIANT *prefdate,VARIANT *pCurrency,VARIANT *ptypeDates,VARIANT *pModFol, VARIANT *pCreditGap, VARIANT *pRet);
	STDMETHOD(ARM_Credit_CashFlows)(VARIANT *pCoefs,VARIANT *pRet);
	STDMETHOD(ARM_Credit_Version)(VARIANT *pRet);
	STDMETHOD(ARM_Credit_CDS) (double pEffectiveDate,double pEndDate,double pSpread,BSTR pFixingFreq,BSTR pDayCountFrq,double pFirst_period_refdate,double pFixedPayerAmount,double pFloatingPayerAmount,BSTR StubRule,BSTR pCurrency,BSTR Adjusetd,int CreditLag,BSTR IncludeMatu,double StartProtection, double EndProtection,BSTR name,double Binary, double pFstCpnEffDate, VARIANT *pRet);
	STDMETHOD(ARM_Credit_CorrelationSmile)(VARIANT *pPricerId,double pSmileType,double pMktPrice ,double pSeed, double pUpfrontPay,int dataType,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CpteDefLegTrancheEquivalente)(VARIANT *pPricerId,double pStrikeUp1, double pStrikeUp2, double pBC1, double pBC2, double pSmileType, VARIANT *pRet);
	STDMETHOD(ARM_Credit_CpteStrikesEquivalents)(VARIANT *pPricerId,double pStrikeDpwn, double pStrikeUp, double pIndexNotional, double pELTranche, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Pricer)(BSTR pSecurity,BSTR pModel,BSTR pPricerType,int Nbpaths, BSTR pParameters,double AsOf, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Portfolio)(VARIANT *pSecuritiesID,BSTR parameters, VARIANT *pRet);
	STDMETHOD(ARM_Credit_SetCorrelationMatrix)(BSTR pModelMultiCurves,BSTR pCorrMatrix, BSTR *pRet);
	STDMETHOD(ARM_Credit_DPMktDataFromSummit)(double AsOfDate,BSTR Issuer,BSTR CurveName,BSTR Parameter,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetDPFromSummit)(double AsOfDate,BSTR Issuer,BSTR CurveName,BSTR ircurveId,BSTR label,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CloneCorrMatrixBary)(BSTR CorrMatrixId,double Beta,int UpOrDown,BSTR *pRet);
	STDMETHOD(ARM_Credit_SetCorrelation)(BSTR pModelMultiCurvesId,BSTR pCorrelationId,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CorrelationStrike)(VARIANT* pLabels,VARIANT* pVolCurves,VARIANT* pProportions,VARIANT* pSmileStrikeLow,VARIANT* pSmileStrikeHigh,VARIANT *pIndexVector,double AsOf,BSTR Name, BSTR *pRet);
	STDMETHOD(ARM_Credit_Beta_Correlation)(VARIANT *pLabels,VARIANT *pCoefs,double AsOf,BSTR Name,VARIANT *pRet);
	STDMETHOD(ARM_Credit_ConstantDefaultCurve)(double AsOfDate,VARIANT *pTenors,VARIANT *pRates,double  Recovery,BSTR IRCurveId,BSTR bCurrency,BSTR	bLabel,BSTR AdjCalType,BSTR IsSummit, VARIANT *pRet);
	STDMETHOD(ARM_Credit_DefProbModelNew) (BSTR pDefProb,BSTR pIrCurve,BSTR pVolCurve,VARIANT *pRet);
	STDMETHOD(ARM_Credit_ZeroCouponDefaultCurveFromSummit)(double AsOfDate,BSTR   bIssuer,BSTR   bCurrency,BSTR   bCvName,BSTR   IRCurveId,BSTR   bLabel,VARIANT *pRet);
	STDMETHOD(ARM_Credit_GetInitialCurveFromSummit)(VARIANT *pIndex, VARIANT *pCurrency, VARIANT *pCvName, VARIANT *pDate, VARIANT* pTypeValue, VARIANT *pRet); 
	STDMETHOD(ARM_Credit_CMTranche) (VARIANT *pEffectiveDate,VARIANT *pEndDate,VARIANT *pSpread,VARIANT *pMezzAmount,VARIANT *pSubAmount,VARIANT *pLabels,VARIANT *pNotionals,VARIANT* pIndex, VARIANT *pFreqFeeLeg,VARIANT *pDayCountFrq,VARIANT *pFirst_period_refdate,VARIANT *pAccruedOnDefault,VARIANT *Currency,VARIANT *pPayCreditLag,VARIANT *pStub,VARIANT *pFreqDefLeg,VARIANT *pBinary,VARIANT *pPayCal,BSTR LongOrShortRisk,double  TradedNotional,double  FwdFixedDate,BSTR IncludeMatu,double pFstCpnEffDate, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Index) (VARIANT *pLabels,double YearFrac,double MaturityDate,double Spread, BSTR Method, BSTR ccy,BSTR Basis,BSTR ResetFreq,BSTR PayFreq,BSTR DefCurve,BSTR VolCurve,BSTR fwdRule,BSTR resetTiming,int resetGap,BSTR payTiming,int payGap,BSTR intRule,BSTR AdjCalType, int cm_resetWeekDay,int cm_resetOccur,BSTR *pRet);
	STDMETHOD(ARM_Credit_Parameters)(VARIANT *pCoefs,long nbcols, VARIANT *pRet);
 	STDMETHOD(ARM_Credit_CDSIndex) (VARIANT *pEffectiveDate,VARIANT *pEndDate,VARIANT *pSpread,VARIANT *pIndex,VARIANT *pFixingFreq,VARIANT *pDayCountFrq,VARIANT *pFirst_period_refdate,VARIANT *pFixedPayerAmount,VARIANT *pFloatingPayerAmount,BSTR StubRule,VARIANT *pCurrency,BSTR Adjusetd,int pCreditLag,BSTR IncludeMatu,VARIANT* StartProtection, VARIANT* EndProtection,VARIANT *pRet);
	STDMETHOD(ARM_Credit_Option) (VARIANT *pUnderlyingPricer,double pStrike,VARIANT *pOptionType, VARIANT *pKoType, VARIANT *pRet);
	STDMETHOD(ARM_Credit_FwdSpreadPricer) (VARIANT *pPricer,VARIANT *pMaturity1, VARIANT *pMaturity2, VARIANT *pRet);
	STDMETHOD(ARM_Credit_ImpliedVol) (VARIANT *pPricer, double pMktPrice, VARIANT *pRet) ;
	STDMETHOD(ARM_Credit_VirtualCdsSpread) (VARIANT *pPricer, VARIANT *pMaturity, VARIANT *pRet);
	STDMETHOD(ARM_Credit_BSGreeks)(VARIANT *pPricer,VARIANT *pGreekType, VARIANT *pRet);
	STDMETHOD(ARM_GETINSTRUMENTFROMSUMMIT)(BSTR SummitId,BSTR Type,double AsOf,BSTR Exoticfilter,VARIANT* pRet);
	STDMETHOD(ARM_Credit_GetEqStrikeDown)(BSTR correlId,BSTR indexname, double *pRet);
	STDMETHOD(ARM_Credit_GetEqStrikeUp)(BSTR correlId,BSTR indexname,double *pRet);
	STDMETHOD(ARM_Credit_GetCorrelStrikeDown)(BSTR correlId,double yfmaturity, double *pRet);
	STDMETHOD(ARM_Credit_GetCorrelStrikeUp)(BSTR correlId,double yfmaturity,double *pRet);
	STDMETHOD(ARM_Credit_GetCorrelation)(BSTR ModelId,BSTR *pRet);
	STDMETHOD(ARM_Credit_GetModelFromSummit)(BSTR IRcurve,BSTR IDSummit,BSTR type,VARIANT *pRet);
	STDMETHOD(ARM_Credit_SetVolatility) (VARIANT *pPricer, BSTR pVolCurve, VARIANT *pRet) ;
	STDMETHOD(ARM_NextCpnDate) (double AsOfDate,double maturity,BSTR frequency,BSTR rule,BSTR currency,BSTR intrule,double *pRet) ;
	STDMETHOD(ARM_Credit_SetProportionsInfos)(BSTR correlId,BSTR IndexName,double proportion,double forcedstrikelow,double forcedstrikehigh) ;
	STDMETHOD(ARM_Credit_CptImplCvForCDO2)(BSTR pricerId,BSTR Name,BSTR Tenor,double *pRet);
	STDMETHOD(ARM_Credit_AddPeriod)(double pAsOf, BSTR Maturity,BSTR pCcy,BSTR AdjRule,BSTR AdjCDS, VARIANT *pRet);
	STDMETHOD(ARM_Credit_Str_CDS)(VARIANT* Plot, VARIANT* Spread, double Margin, double Maturity, double Pas, double Recov, double Nominal, double Taux, BSTR DateToShift, double ValToShift, BSTR Fonction, VARIANT *pRet);
	STDMETHOD(ARM_Credit_SetCoupons)(BSTR CdsorCdoId,BSTR CouponsId,BSTR TypesId,BSTR PartId,BSTR *pRet);
	STDMETHOD(ARM_Credit_SetRiskyProfile)(BSTR CdsorCdoId,BSTR CouponsId,BSTR TypesId,BSTR *pRet);
	STDMETHOD(ARM_Credit_InputDefaultCurve)(double AsOfDate,VARIANT *pDates,VARIANT *pRates,double  Recovery, BSTR IRCurveId,BSTR bCurrency,BSTR bLabel, BSTR bInterpolType, VARIANT *pRet);
	STDMETHOD(ARM_Credit_DataFromLabel)(BSTR pPricer,BSTR pLabel,double *pRet);
	STDMETHOD(ARM_Credit_GenerateImpliedCurve)(BSTR pricerId, BSTR *pRet);
	STDMETHOD(ARM_Credit_GetEqStrike)(BSTR CorrelId,BSTR IndexName,BSTR UpOrLow,  VARIANT *pRetMatu, VARIANT *pRetStrikes);
	STDMETHOD(ARM_Credit_DefaultIntensity)(BSTR pricerId,double Maturity,double *pRet);	
	STDMETHOD(ARM_Credit_EmptyLeg)(VARIANT *pRet);
	STDMETHOD(ARM_Credit_IRLEGTOCREDITLEG)(BSTR SwapLegId,BSTR LegType,BSTR creditindexId,BSTR PricerId,BSTR *pRet);
	STDMETHOD(ARM_Credit_Collateral)(VARIANT *pLabels,VARIANT *pNotionals,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CDSGEN)(BSTR FeeLegId,BSTR DefLegId,double RcvFee,double TradedNot,BSTR *pRet);
	STDMETHOD(ARM_Credit_NTDGEN)(BSTR CdsId,int firstnumdef,int lastnumdef,BSTR CollateralId,double binary,double rcvfee,BSTR *pRet);
	STDMETHOD(ARM_Credit_CDOGEN)(BSTR CdsId,double subamount,BSTR CollateralId,double binary,double rcvfee,BSTR *pRet);
	STDMETHOD(ARM_Credit_GenLeg)(double	StartDate,double EndDate,double	FixedRate,double	FixedNotional,BSTR RefValNotional,BSTR RefValRate,BSTR XChangeNotional, BSTR Frequency,BSTR Basis, BSTR payTiming, BSTR intrule,BSTR stubrule, BSTR ccyid, BSTR paycalname,double refdate,BSTR includematurity,BSTR adjstartdate,BSTR legtype,BSTR indexobj,int creditlag,double binary,BSTR Name,BSTR Nxchange,VARIANT *pRet);
	STDMETHOD(ARM_Credit_CDO_SQUARE_GEN)(BSTR CdsId,double subamount,BSTR portfolioId,double binary,double rcvfee,BSTR *pRet);
	STDMETHOD(ARM_Credit_GetLastFixingDate)(BSTR CdsId,VARIANT asofDate,VARIANT*pRet);
	STDMETHOD(ARM_Credit_SetPastFixing)(BSTR CdsId,VARIANT resetDate,VARIANT fixingValue,VARIANT*pRet);
	STDMETHOD(ARMGetFixing)(BSTR source,BSTR index,BSTR term,BSTR ccy,double date, double *pRet);
	STDMETHOD(ARM_Credit_SetPricerForRatesComputation)(BSTR legId,BSTR pricerId,BSTR *pRet);
	STDMETHOD(ARM_Credit_SetMatuLabel)(BSTR pCurveId, VARIANT *pMatuLabels, double *pRet);
	STDMETHOD(ARM_Credit_SetRecovCoef)(BSTR pSecId, double RecovCoef, double *pRet);
	STDMETHOD(ARM_Credit_SetFees)(BSTR securityId,BSTR RefvalueId);
	STDMETHOD(ARM_Credit_GetBounds)(BSTR instId_,double* low,double* up);
	STDMETHOD(ARM_Credit_SetInterpolationType)(VARIANT * pVolCurveId,VARIANT *pInterpMeth , VARIANT *pRet);
	STDMETHOD(ARM_Credit_Customized_CDO)(VARIANT* pLabels, BSTR pDefaultLeg, BSTR pPremiumLeg, BSTR pParameters, BSTR *pRet);
	STDMETHOD(ARM_Credit_CLN)(double pEffectiveDate,double	pEndDate,double	pSpread,BSTR pIndexId,double pFirst_period_refdate,double pFstCpnEffDate,double	pNotional,BSTR AccOnDef, BSTR pDayCount,BSTR pDecompFreq,BSTR StubRule,double resetgap,BSTR pCurrency,BSTR ResetCal,BSTR PayCal,BSTR Nxchange,BSTR pIncludeMaturity,BSTR AdjustedStartDate,double Binary,BSTR *pRet);
	//End Section Credit
	*/ 
};

#endif // !defined(AFX_ARMMODULE_H__F153E6FC_B3F3_4A73_9A2D_6A0FE0A84A8F__INCLUDED_)


