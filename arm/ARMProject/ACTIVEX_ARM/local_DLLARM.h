/* this ALWAYS GENERATED file contains the definitions for the interfaces */


/* File created by MIDL compiler version 5.01.0164 */
/* at Fri May 18 15:56:40 2007
 */
/* Compiler settings for .\local_DLLARM.idl:
    Oicf (OptLev=i2), W1, Zp8, env=Win32, ms_ext, c_ext
    error checks: allocation ref bounds_check enum stub_data 
*/
//@@MIDL_FILE_HEADING(  )


/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 440
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif // __RPCNDR_H_VERSION__

#ifndef COM_NO_WINDOWS_H
#include "windows.h"
#include "ole2.h"
#endif /*COM_NO_WINDOWS_H*/

#ifndef __local_DLLARM_h__
#define __local_DLLARM_h__

#ifdef __cplusplus
extern "C"{
#endif 

/* Forward Declarations */ 

#ifndef __IARMModule_FWD_DEFINED__
#define __IARMModule_FWD_DEFINED__
typedef interface IARMModule IARMModule;
#endif 	/* __IARMModule_FWD_DEFINED__ */


#ifndef __ARMModule_FWD_DEFINED__
#define __ARMModule_FWD_DEFINED__

#ifdef __cplusplus
typedef class ARMModule ARMModule;
#else
typedef struct ARMModule ARMModule;
#endif /* __cplusplus */

#endif 	/* __ARMModule_FWD_DEFINED__ */


/* header files for imported files */
#include "oaidl.h"
#include "ocidl.h"

void __RPC_FAR * __RPC_USER MIDL_user_allocate(size_t);
void __RPC_USER MIDL_user_free( void __RPC_FAR * ); 

#ifndef __IARMModule_INTERFACE_DEFINED__
#define __IARMModule_INTERFACE_DEFINED__

/* interface IARMModule */
/* [unique][helpstring][dual][uuid][object] */ 


EXTERN_C const IID IID_IARMModule;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("CA2D4B2C-934D-4643-BA08-CAEDAB80B9B8")
    IARMModule : public IDispatch
    {
    public:
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMNextBusinessDay( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [in] */ long pNbDays,
            /* [retval][out] */ double __RPC_FAR *pDate2) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMAdjustToBusDate( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [in] */ BSTR pRule,
            /* [retval][out] */ double __RPC_FAR *pDate2) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFreeObject( 
            /* [in] */ BSTR pId,
            /* [retval][out] */ long __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMIsBusinessDay( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [retval][out] */ long __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetZCFromSummit( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFreeAllObjects( 
            /* [retval][out] */ long __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMYcMod( 
            /* [in] */ BSTR pZc,
            /* [defaultvalue][in][optional] */ BSTR pZcDiscount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMForwardYield( 
            /* [in] */ BSTR pZc,
            /* [in] */ double pMatu1,
            /* [in] */ double pMatu2,
            /* [in] */ BSTR pMeth,
            /* [in] */ BSTR pAdj,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMDiscountYield( 
            /* [in] */ VARIANT __RPC_FAR *pZc,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pMeth,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLiborSwap( 
            /* [in] */ VARIANT __RPC_FAR *pStartDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pLiborType,
            /* [in] */ VARIANT __RPC_FAR *pRecOrPay,
            /* [in] */ VARIANT __RPC_FAR *pFixedRate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pFloatingDaycount,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwapPriceToRate( 
            /* [in] */ VARIANT __RPC_FAR *pSwap,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pPrice,
            /* [in] */ VARIANT __RPC_FAR *pModel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMPrice( 
            /* [in] */ VARIANT __RPC_FAR *pSec,
            /* [in] */ VARIANT __RPC_FAR *pModel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBetweenDates( 
            /* [in] */ VARIANT __RPC_FAR *pDate1,
            /* [in] */ VARIANT __RPC_FAR *pDate2,
            /* [in] */ VARIANT __RPC_FAR *pDaycount,
            /* [in] */ VARIANT __RPC_FAR *pIsYearFrac,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMAddPeriod( 
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pFreq,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [in] */ VARIANT __RPC_FAR *pNbPeriods,
            /* [in] */ VARIANT __RPC_FAR *pAdjRule,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMIsoCcy( 
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetSpotDays( 
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetLiborIndexDaycount( 
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetLiborTerm( 
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFixedDayCount( 
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFixedPayFreq( 
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMComputeVolatility( 
            /* [in] */ VARIANT __RPC_FAR *pVol,
            /* [in] */ VARIANT __RPC_FAR *pmatu,
            /* [in] */ VARIANT __RPC_FAR *pStrike,
            /* [in] */ VARIANT __RPC_FAR *pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMVolCurv( 
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pStrike,
            /* [in] */ VARIANT __RPC_FAR *pVol,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pStrikeType,
            /* [in] */ BSTR pVolType,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIndexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetVolCubeFromSummit( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [in] */ BSTR pType,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pSmiles,
            /* [defaultvalue][in][optional] */ BSTR pTypeCube,
            /* [defaultvalue][in][optional] */ BSTR indexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_FxConvert( 
            /* [in] */ VARIANT __RPC_FAR *pccy1,
            /* [in] */ VARIANT __RPC_FAR *pccy2,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pCvName,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_DiscountPrice( 
            /* [in] */ VARIANT __RPC_FAR *pcurve,
            /* [in] */ VARIANT __RPC_FAR *pmatu,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Delivery( 
            /* [in] */ VARIANT __RPC_FAR *pAsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pTenorContract,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CptInterpolDefCurve( 
            /* [in] */ BSTR pCurve,
            /* [in] */ VARIANT pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DefaultProba( 
            /* [in] */ VARIANT __RPC_FAR *pCurve,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetBeta( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Mezzanine( 
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [in] */ double pMezzAmount,
            /* [in] */ double pSubAmount,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [defaultvalue][in] */ BSTR pFreqFeeLeg,
            /* [defaultvalue][in] */ BSTR pDayCount,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in] */ BSTR Currency,
            /* [defaultvalue][in] */ double pPayCreditLag,
            /* [defaultvalue][in] */ BSTR pStub,
            /* [defaultvalue][in] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in] */ double pBinary,
            /* [defaultvalue][in] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR adjStartDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ModelMultiCurves( 
            /* [in] */ VARIANT __RPC_FAR *pIRcurve,
            /* [in] */ VARIANT __RPC_FAR *pDefCurves,
            /* [in] */ VARIANT __RPC_FAR *pRecoveryRates,
            /* [defaultvalue][in][optional] */ BSTR CorrelId,
            /* [defaultvalue][in][optional] */ BSTR VolCurve,
            /* [defaultvalue][in][optional] */ BSTR CpnInfCurve,
            /* [defaultvalue][in][optional] */ BSTR CpnIRCurve,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FTD( 
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pIssuerNotional,
            /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
            /* [in] */ VARIANT __RPC_FAR *pStub,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pintRule,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pstartAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_NTD( 
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [defaultvalue][in] */ int pFirstNumDefault,
            /* [defaultvalue][in] */ int pLastNumDefault,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ BSTR pFixingFreq,
            /* [defaultvalue][in] */ BSTR pDayCountFrq,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ double pIssuerNotional,
            /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in] */ BSTR pCurrency,
            /* [defaultvalue][in] */ double pPayCreditLag,
            /* [defaultvalue][in] */ BSTR pStub,
            /* [defaultvalue][in] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in] */ double pBinary,
            /* [defaultvalue][in] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR startAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Price( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pAsofDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_RiskyDuration( 
            /* [in] */ BSTR pDefCurve,
            /* [in] */ VARIANT date,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDONPV( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pCPTTYPE,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CorrMatrix( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [defaultvalue][in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ExtractCorrMatrix( 
            /* [in] */ VARIANT __RPC_FAR *pCorrMatrixId,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Spread( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pMTM,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetLabel( 
            /* [in] */ VARIANT __RPC_FAR *pCurveId,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetLabel( 
            /* [in] */ VARIANT __RPC_FAR *pCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Sensitivity( 
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [in] */ VARIANT __RPC_FAR *pType,
            /* [in] */ VARIANT __RPC_FAR *pPlot,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [in] */ VARIANT __RPC_FAR *pEpsilon,
            /* [defaultvalue][in][optional] */ double epsilonGamma,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetCleanSpreadTranche( 
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [in] */ VARIANT __RPC_FAR *pPlot,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetDefProbTranche( 
            /* [in] */ BSTR PricerId,
            /* [in] */ double Yearterm,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetDuration( 
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenSchedule( 
            /* [in] */ VARIANT __RPC_FAR *pAccStartDate,
            /* [in] */ VARIANT __RPC_FAR *pAccEndDate,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *prefDate,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *ptypeDates,
            /* [in] */ VARIANT __RPC_FAR *pModFol,
            /* [in] */ VARIANT __RPC_FAR *pCreditGap,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CashFlows( 
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_View( 
            /* [in] */ VARIANT __RPC_FAR *pObjet,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Version( 
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDS( 
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [defaultvalue][in] */ BSTR pFixingFreq,
            /* [defaultvalue][in] */ BSTR pDayCountFrq,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ double pFixedPayerAmount,
            /* [defaultvalue][in] */ double pFloatingPayerAmount,
            /* [defaultvalue][in] */ BSTR StubRule,
            /* [defaultvalue][in] */ BSTR pCurrency,
            /* [defaultvalue][in] */ BSTR Adjusted,
            /* [defaultvalue][in] */ int CreditDefLag,
            /* [defaultvalue][in] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double StartProtection,
            /* [defaultvalue][in] */ double EndProtection,
            /* [defaultvalue][in] */ BSTR name,
            /* [defaultvalue][in] */ double binary,
            /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR StartAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetVolFromSummit( 
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [in] */ VARIANT __RPC_FAR *pCvName,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pType,
            /* [in] */ VARIANT __RPC_FAR *pMatIndex,
            /* [in] */ VARIANT __RPC_FAR *pImpOrHist,
            /* [defaultvalue][in][optional] */ BSTR pindexId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMParallelShift( 
            /* [in] */ BSTR pZc,
            /* [in] */ double pBump,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBumpVolatility( 
            /* [in] */ BSTR pVolCurv,
            /* [in] */ double pBump,
            /* [in] */ long pNthLine,
            /* [in] */ long pNthCol,
            /* [in] */ BSTR pIsCumul,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBsSmiledModel( 
            /* [in] */ double pDate,
            /* [in] */ double pSpot,
            /* [in] */ BSTR pDividend,
            /* [in] */ BSTR pDiscrate,
            /* [in] */ BSTR pVolATM,
            /* [in] */ BSTR pRo,
            /* [in] */ BSTR pNu,
            /* [in] */ BSTR pIsSABR,
            /* [defaultvalue][in][optional] */ BSTR pBeta,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSetEtoolkit( 
            /* [in] */ BSTR pUserName,
            /* [in] */ BSTR pPassWord,
            /* [in] */ BSTR pDatabaseContext,
            /* [in] */ BSTR pItConfigDomainDir,
            /* [in] */ BSTR pItDomainName,
            /* [retval][out] */ long __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMConnectionEtoolkit( 
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMVolFlat( 
            /* [in] */ double pVol,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMVolCube( 
            /* [in] */ BSTR pATMVol,
            /* [in] */ VARIANT __RPC_FAR *pSmileCurveIds,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [defaultvalue][in] */ BSTR pVolType,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMDeconnectionEtoolkit( 
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMZcFlat( 
            /* [in] */ double pZc,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBsModel( 
            /* [in] */ double pDate,
            /* [in] */ double pSpot,
            /* [in] */ BSTR pDividend,
            /* [in] */ BSTR pDiscrate,
            /* [in] */ BSTR pVol,
            /* [defaultvalue][in][optional] */ BSTR pTypeStk,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwitchToETK( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwitchToFLATFILE( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMZCLINT( 
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [in] */ BSTR pMeth,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCurrency,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Pricer( 
            /* [in] */ BSTR pSecurity,
            /* [in] */ BSTR pModel,
            /* [in] */ BSTR pPricerType,
            /* [in] */ int l_nbpaths,
            /* [defaultvalue][in][optional] */ BSTR pParameters,
            /* [defaultvalue][in][optional] */ double valuationdate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_zcspreaded( 
            /* [in] */ BSTR zcSprId,
            /* [in] */ BSTR zcInitId,
            /* [in] */ double date,
            /* [in] */ BSTR MMFreq,
            /* [in] */ BSTR SwapFreq,
            /* [in] */ BSTR ccyId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CptBaseCorrelation( 
            /* [in] */ double AsOf,
            /* [in] */ BSTR name,
            /* [in] */ BSTR CalMethod,
            /* [in] */ BSTR IndexId,
            /* [in] */ VARIANT __RPC_FAR *pStrikeLow,
            /* [in] */ VARIANT __RPC_FAR *pStrikeHigh,
            /* [in] */ VARIANT __RPC_FAR *pVMktBid,
            /* [in] */ VARIANT __RPC_FAR *pVMktAsk,
            /* [in] */ VARIANT __RPC_FAR *pVUpfBid,
            /* [in] */ VARIANT __RPC_FAR *pVUpfAsk,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVInitialCorrel,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVDeltaLevrage,
            /* [defaultvalue][in][optional] */ BSTR ModelId,
            /* [defaultvalue][in][optional] */ double integrationStep,
            /* [defaultvalue][in][optional] */ double lagStartDate,
            /* [defaultvalue][in][optional] */ double creditLag,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVectorPrevIndexId,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pMatrixPrevBC,
            /* [defaultvalue][in][optional] */ double step,
            /* [defaultvalue][in][optional] */ BSTR CalMeth,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDO2( 
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pPortfolio,
            /* [in] */ double pSpread,
            /* [in] */ double pSubAmount,
            /* [in] */ double pMezzAmount,
            /* [defaultvalue][in][optional] */ BSTR pFreqFeeLeg,
            /* [defaultvalue][in][optional] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in][optional] */ BSTR pDayCountFrq,
            /* [defaultvalue][in][optional] */ double pFirst_period_refdate,
            /* [defaultvalue][in][optional] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in][optional] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ double pPayCreditLag,
            /* [defaultvalue][in][optional] */ BSTR pStub,
            /* [defaultvalue][in][optional] */ double pBinary,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR CrossSub,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Portfolio( 
            /* [in] */ VARIANT __RPC_FAR *pSecuritiesID,
            /* [defaultvalue][in][optional] */ BSTR Parameters,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCreateZCSwapInt( 
            /* [in] */ double pDate,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [in] */ BSTR pMMVsFut,
            /* [in] */ BSTR pSwapVsFut,
            /* [in] */ BSTR pRaw,
            /* [in] */ BSTR pInterp,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetInitialCurveFromSummit( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in] */ BSTR pAdjOrNot,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetRate) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMTHREEMONTHFUT( 
            /* [in] */ BSTR pDelivery,
            /* [in] */ long pMarket,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFutPibor( 
            /* [in] */ BSTR pDelivery,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMIRFUT( 
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pIdUnderlying,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLibor( 
            /* [in] */ BSTR pLiborTypeId,
            /* [in] */ BSTR pCcyId,
            /* [in] */ BSTR pResetFreqId,
            /* [in] */ BSTR pPayFreqId,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [defaultvalue][in] */ BSTR pBasis,
            /* [defaultvalue][in] */ BSTR pIntrule,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLiborSwaption( 
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pReceiveOrPay,
            /* [in] */ double pStrike,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pLiborType,
            /* [in] */ double pSpread,
            /* [in] */ BSTR pExerciseType,
            /* [in] */ BSTR pResetFreq,
            /* [in] */ BSTR pPayFreq,
            /* [in] */ BSTR pCcyId,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFixedLeg( 
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pReceiveOrPay,
            /* [in] */ double pFixRate,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [defaultvalue][in][optional] */ BSTR pFreq,
            /* [defaultvalue][in][optional] */ BSTR pDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ BSTR pPayCalName,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetCorrelationMatrix( 
            /* [in] */ BSTR pModelMultiCurvesId,
            /* [in] */ BSTR pCorrMatrixId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DPMktDataFromSummit( 
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR Issuer,
            /* [in] */ BSTR CurveName,
            /* [in] */ BSTR Parameter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetDPFromSummit( 
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR Issuer,
            /* [in] */ BSTR CurveName,
            /* [in] */ BSTR ircurveId,
            /* [in] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CloneCorrMatrixBary( 
            /* [in] */ BSTR CorrMatrixId,
            /* [in] */ double Beta,
            /* [in] */ int UpOrDown,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ConstantDefaultCurve( 
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [in] */ VARIANT __RPC_FAR *pRates,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR Ccy,
            /* [in] */ BSTR Label,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ BSTR IsSummit,
            /* [defaultvalue][optional][in] */ BSTR calibrationData,
            /* [defaultvalue][optional][in] */ int lag,
            /* [optional][in] */ BSTR calibrationAlgo,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DefProbModelNew( 
            /* [in] */ BSTR pDefCurve,
            /* [in] */ BSTR pIRcurve,
            /* [defaultvalue][in][optional] */ BSTR VolCurve,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ZeroCouponDefaultCurveFromSummit( 
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR bIssuer,
            /* [in] */ BSTR bCurrency,
            /* [in] */ BSTR bCvName,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR bLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBsSlModel( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pVolSpreadLock,
            /* [in] */ BSTR pCvCapVol,
            /* [in] */ BSTR pCvIndexVol,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFXVolFromSummit( 
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [defaultvalue][in][optional] */ BSTR pType,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFXCorrelFromSummit( 
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMInfocentreConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGlobDFBS( 
            /* [in] */ BSTR pDomBSId,
            /* [in] */ BSTR pDomCurrId,
            /* [in] */ BSTR pFrgBSId,
            /* [in] */ BSTR pFrgCurrId,
            /* [in] */ BSTR pFxVolCrvId,
            /* [in] */ BSTR pFFxCorrId,
            /* [in] */ BSTR pRatesCorrId,
            /* [defaultvalue][in][optional] */ BSTR pFxVolModelId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetInitialCurveFromSummit( 
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *pCvName,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *value,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetCorrelFromSummit( 
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ BSTR pIndex2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMDFFXBS( 
            /* [in] */ BSTR pDVolId,
            /* [in] */ BSTR pFVolId,
            /* [in] */ BSTR pDZcId,
            /* [in] */ BSTR pFZcId,
            /* [in] */ BSTR pDFxCorrId,
            /* [in] */ BSTR pFFxCorrId,
            /* [in] */ BSTR pFxVolId,
            /* [in] */ double pRatesCorr,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMTRIBSMODEL( 
            /* [in] */ BSTR pModel1,
            /* [in] */ BSTR pModel2,
            /* [in] */ BSTR pDiscModel,
            /* [in] */ BSTR pFX1DiscVol,
            /* [in] */ BSTR pFX2DiscVol,
            /* [in] */ BSTR pIdx1Idx2Corr,
            /* [in] */ BSTR pIdx1DiscIdxCorr,
            /* [in] */ BSTR pIdx2DiscIdxCorr,
            /* [in] */ BSTR pIdx1FxCorr,
            /* [in] */ BSTR pIdx2FxCorr,
            /* [in] */ int pQuantoFlag,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMTRIBSDUAL( 
            /* [in] */ BSTR pModel1,
            /* [in] */ BSTR pModel2,
            /* [in] */ BSTR pDiscModel,
            /* [in] */ BSTR pFX1DiscVol,
            /* [in] */ BSTR pFX2DiscVol,
            /* [in] */ BSTR pIdx1Idx2Corr,
            /* [in] */ BSTR pIdx1DiscIdxCorr,
            /* [in] */ BSTR pIdx2DiscIdxCorr,
            /* [in] */ BSTR pIdx1FxCorr,
            /* [in] */ BSTR pIdx2FxCorr,
            /* [in] */ int pQuantoFlag,
            /* [in] */ double pCorrelForAdj,
            /* [in] */ int pWithslopeflag,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcptBonibor( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pDepart,
            /* [in] */ BSTR pMatStruct,
            /* [in] */ BSTR pMatTot,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pSjUSD,
            /* [in] */ BSTR pTiming,
            /* [in] */ double pBarriere,
            /* [in] */ double pSpdPostBar,
            /* [in] */ double pMarge,
            /* [in] */ double pFunding,
            /* [in] */ BSTR pFundingFreq,
            /* [in] */ double pSpd2phase,
            /* [in] */ double pSoulte,
            /* [in] */ BSTR pYcModId,
            /* [in] */ BSTR pBsModId,
            /* [in] */ BSTR pBsModVolUSDId,
            /* [in] */ BSTR pBsModCorrPlusId,
            /* [in] */ BSTR pBsModCorrMoinsId,
            /* [in] */ BSTR pCrossModId,
            /* [in] */ BSTR pProbaMarge,
            /* [defaultvalue][in][optional] */ double pInt,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CMTranche( 
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pParticipationRate,
            /* [in] */ VARIANT __RPC_FAR *pMezzAmount,
            /* [in] */ VARIANT __RPC_FAR *pSubAmount,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pFreqFeeLeg,
            /* [in] */ VARIANT __RPC_FAR *pDayCount,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
            /* [in] */ VARIANT __RPC_FAR *Currency,
            /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
            /* [in] */ VARIANT __RPC_FAR *pStub,
            /* [in] */ VARIANT __RPC_FAR *pFreqDefLeg,
            /* [in] */ VARIANT __RPC_FAR *pBinary,
            /* [in] */ VARIANT __RPC_FAR *pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ double FwdFixedDate,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Index( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ double YearFrac,
            /* [defaultvalue][in] */ double pSpread,
            /* [defaultvalue][in] */ BSTR Method,
            /* [defaultvalue][in] */ BSTR Basis,
            /* [defaultvalue][in] */ BSTR ResetFreq,
            /* [defaultvalue][in] */ BSTR PayFreq,
            /* [defaultvalue][in] */ BSTR ccy,
            /* [defaultvalue][in] */ BSTR DefCurve,
            /* [defaultvalue][in] */ BSTR fwdRule,
            /* [defaultvalue][in] */ BSTR resetTiming,
            /* [defaultvalue][in] */ int resetGap,
            /* [defaultvalue][in] */ BSTR payTiming,
            /* [defaultvalue][in] */ int payGap,
            /* [defaultvalue][in] */ BSTR intRule,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ int cm_resetWeekDay,
            /* [defaultvalue][in] */ int cm_resetOccur,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Parameters( 
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [in] */ long nbcols,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMAswPrice( 
            /* [in] */ double pMaturity,
            /* [in] */ double pCpn,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pBase,
            /* [in] */ double pMargin,
            /* [defaultvalue][in][optional] */ double pRedemptionPrice,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pFixDecompfreq,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ double pMinValue,
            /* [defaultvalue][in][optional] */ double pMaxValue,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMAswMargin( 
            /* [in] */ double pMaturity,
            /* [in] */ double pCpn,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pBase,
            /* [in] */ double pPrice,
            /* [defaultvalue][in][optional] */ double pRedemptionPrice,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pFixDecompfreq,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ double pMinValue,
            /* [defaultvalue][in][optional] */ double pMaxValue,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDSIndex( 
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pFixedPayerAmount,
            /* [in] */ VARIANT __RPC_FAR *pFloatingPayerAmount,
            /* [in] */ BSTR StubRule,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [defaultvalue][in] */ BSTR Adjusted,
            /* [defaultvalue][in] */ int CreditLag,
            /* [defaultvalue][in] */ BSTR IncludeMaturity,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionStartDate,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionEndDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Option( 
            /* [in] */ VARIANT __RPC_FAR *UnderlyingMaturity,
            /* [in] */ double OptionExpiry,
            /* [in] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ BSTR CdsAdj,
            /* [defaultvalue][in][optional] */ BSTR EndAdj,
            /* [defaultvalue][in][optional] */ double pStrike,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pOptionType,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pKoType,
            /* [defaultvalue][in][optional] */ double Notional,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FwdSpreadPricer( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ double Maturity1,
            /* [in] */ double Maturity2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ImpliedVol( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ double pMktPrice,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_VirtualCdsSpread( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pMaturity,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_BSGreeks( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pGreekType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetCorrelation( 
            /* [in] */ BSTR pModelMultiCurvesId,
            /* [in] */ BSTR pCorrelationId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CorrelationStrike( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pVolCurves,
            /* [in] */ VARIANT __RPC_FAR *pProportions,
            /* [in] */ VARIANT __RPC_FAR *pSmileStrikeLow,
            /* [in] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
            /* [in] */ VARIANT __RPC_FAR *pIndexVector,
            /* [defaultvalue][in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Beta_Correlation( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [defaultvalue][in] */ BSTR idIndex1,
            /* [defaultvalue][in] */ BSTR idIndex2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_GETINSTRUMENTFROMSUMMIT( 
            /* [in] */ BSTR SummitId,
            /* [in] */ BSTR Type,
            /* [defaultvalue][in][optional] */ double AsOf,
            /* [defaultvalue][in][optional] */ BSTR ExoticFilter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFrnPrice( 
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [in] */ double pFacialMargin,
            /* [in] */ double pValoMargin,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ double pFixing,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ double pOutMode,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFrnMargin( 
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [in] */ double pFacialMargin,
            /* [in] */ double pPrice,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ double pFixing,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ double pOutMode,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetEqStrikeDown( 
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetEqStrikeUp( 
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetCorrelStrikeDown( 
            /* [in] */ BSTR CorrelId,
            /* [in] */ double yfmaturity,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetCorrelStrikeUp( 
            /* [in] */ BSTR CorrelId,
            /* [in] */ double yfmaturity,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetCorrelation( 
            /* [in] */ BSTR ModelId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetModelFromSummit( 
            /* [in] */ BSTR IRcurve,
            /* [in] */ BSTR IDSummit,
            /* [in] */ BSTR type,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetVolatility( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR VolCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_NextCpnDate( 
            /* [in] */ double AsOfDate,
            /* [in] */ double maturity,
            /* [in] */ BSTR frequency,
            /* [defaultvalue][in][optional] */ BSTR rule,
            /* [defaultvalue][in][optional] */ BSTR currency,
            /* [defaultvalue][in][optional] */ BSTR intrule,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetProportionsInfos( 
            /* [in] */ BSTR correlId,
            /* [in] */ BSTR IndexName,
            /* [in] */ double proportion,
            /* [defaultvalue][in][optional] */ double forcedstrikelow = -999,
            /* [defaultvalue][in][optional] */ double forcedstrikehigh = -999) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcptDigital( 
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodDelta,
            /* [in] */ BSTR pBsmodVega,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pPorR,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCcyIdx,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pCapOrFloor,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pStrike,
            /* [in] */ BSTR pPayOff,
            /* [in] */ BSTR pSpd,
            /* [in] */ double pResetGap,
            /* [defaultvalue][in][optional] */ double pSpreadBelow,
            /* [defaultvalue][in][optional] */ double pSpreadAbove,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMRefValue( 
            /* [in] */ VARIANT __RPC_FAR *pdates,
            /* [in] */ VARIANT __RPC_FAR *pvalues,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pvalues2,
            /* [defaultvalue][in][optional] */ long valueType,
            /* [defaultvalue][in][optional] */ long conversion,
            /* [defaultvalue][in][optional] */ BSTR calcMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCreateGenCorrelManager( 
            /* [in] */ VARIANT __RPC_FAR *pMktTags,
            /* [in] */ VARIANT __RPC_FAR *pIntraMktTags,
            /* [in] */ VARIANT __RPC_FAR *pCorrelCurveIds,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBSConvAdjust( 
            /* [defaultvalue][in][optional] */ BSTR pSUMMITFormulaeUsed,
            /* [defaultvalue][in][optional] */ BSTR pUseSABRCMS,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBsModelGen( 
            /* [in] */ BSTR pYieldCurve,
            /* [in] */ BSTR pVolatility,
            /* [defaultvalue][in][optional] */ BSTR pCorrMgr,
            /* [defaultvalue][in][optional] */ BSTR pCnvxManager,
            /* [defaultvalue][in][optional] */ BSTR pCapletVol,
            /* [defaultvalue][in][optional] */ BSTR pSpreadLock,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve,
            /* [defaultvalue][in][optional] */ BSTR pCorrel,
            /* [defaultvalue][in][optional] */ BSTR pCashVol,
            /* [defaultvalue][in][optional] */ BSTR pSpreadVol,
            /* [defaultvalue][in][optional] */ BSTR pModelType,
            /* [defaultvalue][in][optional] */ BSTR pSpreadVolType,
            /* [defaultvalue][in][optional] */ BSTR pSabrMod,
            /* [defaultvalue][in][optional] */ BSTR pLnorNorVol,
            /* [defaultvalue][in][optional] */ long pNumSteps,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcptSPTQTF( 
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ BSTR pIndexPhase2,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ BSTR pFreqPPhase1,
            /* [in] */ BSTR pFreqPPhase2,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqPPhase3,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pResetTimingPhase2,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCcyIdx,
            /* [in] */ BSTR pDayCount,
            /* [in] */ double pFee,
            /* [in] */ BSTR pIsRateFixedPhase2,
            /* [in] */ double pFixedRatePhase2,
            /* [in] */ BSTR pBarrier,
            /* [in] */ BSTR pSpdPhase1,
            /* [in] */ BSTR pSpdPhase1Fund,
            /* [in] */ BSTR pSpdPhase2Tf,
            /* [in] */ BSTR pSpdPhase2fund,
            /* [in] */ BSTR pSpdPhase3,
            /* [in] */ BSTR pSpdPhase3fund,
            /* [in] */ double pResetGapPhase1,
            /* [in] */ double pResetGapPhase2,
            /* [in] */ double pResetGapPhase3,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodDeltaCcy1,
            /* [in] */ BSTR pBsmodVegaCcy1,
            /* [defaultvalue][in][optional] */ BSTR pBsmodDeltaCcy2,
            /* [defaultvalue][in][optional] */ BSTR pBsmodVegaCcy2,
            /* [defaultvalue][in][optional] */ BSTR pBsmodFxCorrel,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMDisplaySchedule( 
            /* [in] */ BSTR legId,
            /* [in] */ BSTR dataType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMIrIndex( 
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ double pMaturity,
            /* [defaultvalue][in][optional] */ BSTR pCompMethod,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ double pPayGap,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIndexType,
            /* [defaultvalue][in][optional] */ double pDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwapleg( 
            /* [in] */ BSTR pIndexId,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pRecOrPay,
            /* [in] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMConstRefvalue( 
            /* [in] */ double pValue,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CptImplCvForCDO2( 
            /* [in] */ BSTR pricerId,
            /* [in] */ BSTR Name,
            /* [in] */ BSTR Tenor,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_AddPeriod( 
            /* [in] */ double pAsOf,
            /* [in] */ BSTR Maturity,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR AdjRule,
            /* [defaultvalue][in][optional] */ BSTR AdjCDS,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMShutdownETK( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetCoupons( 
            /* [in] */ BSTR CdsorCdoId,
            /* [in] */ BSTR CouponsId,
            /* [in] */ BSTR TypesId,
            /* [defaultvalue][in][optional] */ BSTR PartId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_InputDefaultCurve( 
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pDates,
            /* [in] */ VARIANT __RPC_FAR *pRates,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR bCurrency,
            /* [in] */ BSTR bLabel,
            /* [defaultvalue][in][optional] */ BSTR bInterpolType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetInitialVolFromSummit( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [in] */ BSTR pType,
            /* [defaultvalue][in][optional] */ BSTR pMatIndex,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBond( 
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedempPrice,
            /* [in] */ double pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ double pSettleGap,
            /* [defaultvalue][in][optional] */ double pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFCreateOATLeg( 
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pInfIdx,
            /* [in] */ BSTR pRcvOrPay,
            /* [defaultvalue][in][optional] */ BSTR pInterpType,
            /* [defaultvalue][in][optional] */ double pLeverage,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pResetNumGap,
            /* [defaultvalue][in][optional] */ double pResetDenomGap,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ double pPayGap,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR pFinalNotionalType,
            /* [defaultvalue][in][optional] */ double pFirstReset,
            /* [defaultvalue][in][optional] */ double pCoMultiple,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwap( 
            /* [in] */ BSTR pSwapleg1,
            /* [in] */ BSTR pSwapleg2,
            /* [defaultvalue][in][optional] */ double pMinPay,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMPToYield( 
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pPrice,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMYToPrice( 
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pYield,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMYToDuration( 
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pActuRate,
            /* [defaultvalue][in][optional] */ double pFlagCpn,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLiborleg( 
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pLiborType,
            /* [in] */ BSTR pRecOrPay,
            /* [defaultvalue][in][optional] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFReq,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [defaultvalue][in][optional] */ BSTR pCpnDaycount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMImpliedSpread( 
            /* [in] */ BSTR pSwap,
            /* [in] */ BSTR pModel,
            /* [in] */ double pPrice,
            /* [in] */ double pLeg1or2,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMDiscountPrice( 
            /* [in] */ BSTR pZeroCurve,
            /* [in] */ double pMatu,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFCreateCurve( 
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pIndexName,
            /* [in] */ double pCPIIndexValue,
            /* [in] */ double pCPIIndexDate,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [defaultvalue][in][optional] */ BSTR pMonthlyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pDCFMonthly,
            /* [defaultvalue][in][optional] */ BSTR pDCFDaily,
            /* [defaultvalue][in][optional] */ BSTR pExtrapolType,
            /* [defaultvalue][in][optional] */ BSTR pResetManager,
            /* [defaultvalue][in][optional] */ BSTR pSeasonManager,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFInterpCPI( 
            /* [in] */ BSTR pZc,
            /* [in] */ double pCPIDate,
            /* [defaultvalue][in][optional] */ BSTR pDCFlag,
            /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pCPIlag,
            /* [defaultvalue][in][optional] */ double pWeight,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFSeasonManager( 
            /* [in] */ VARIANT __RPC_FAR *pMonthList,
            /* [in] */ VARIANT __RPC_FAR *pValues,
            /* [defaultvalue][in][optional] */ BSTR pSeasonAdjMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFResetManager( 
            /* [in] */ VARIANT __RPC_FAR *pDatas,
            /* [in] */ double pNbIndex,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMINFYcMod( 
            /* [in] */ BSTR pYieldCurve,
            /* [in] */ BSTR pInfCurve,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBaseReplicationConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetInitialFXVolFromSummit( 
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCvName,
            /* [defaultvalue][in][optional] */ BSTR pImpOrHist,
            /* [defaultvalue][in][optional] */ BSTR pVolType,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCreateZCFromSummit( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pAdj,
            /* [defaultvalue][in][optional] */ BSTR pRaw,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBumpCurve( 
            /* [in] */ BSTR pZc,
            /* [in] */ double pEpsilon,
            /* [defaultvalue][in][optional] */ long pMethod,
            /* [defaultvalue][in][optional] */ BSTR pPlot,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMAccrued( 
            /* [in] */ BSTR pSec,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pModel,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSwitchToWSETK( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DataFromLabel( 
            /* [in] */ BSTR pricer,
            /* [in] */ BSTR label,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeBilibor( 
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pDateSecondPhase,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexSF,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodFund,
            /* [in] */ BSTR pBsmodDeltaCcy1,
            /* [in] */ BSTR pBsmodDeltaFund,
            /* [in] */ BSTR pBsmodDeltaCcy2,
            /* [in] */ BSTR pBsmodFxCorrel,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqRFund,
            /* [in] */ BSTR pFreqPSF,
            /* [in] */ BSTR pFreqRSF,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pResetTimingSF,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pDayCountSF,
            /* [in] */ double pSpdPF,
            /* [in] */ double pSpdSF,
            /* [in] */ double pSpdfund,
            /* [in] */ double pSpdfund2,
            /* [in] */ double pResetGap,
            /* [in] */ double pResetGapSF,
            /* [in] */ BSTR pAmort,
            /* [in] */ double pRefDate,
            /* [in] */ double pFee,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenerateImpliedCurve( 
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetEqStrike( 
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [in] */ BSTR UpOrLow,
            /* [out] */ VARIANT __RPC_FAR *pRetMatu,
            /* [out] */ VARIANT __RPC_FAR *pRetStrikes) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeOptilix( 
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pDateSecondPhase,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexSF,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodFund,
            /* [in] */ BSTR pBsmodDeltaCcy,
            /* [in] */ BSTR pBsmodDeltaFund,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqRFund,
            /* [in] */ BSTR pFreqPSF,
            /* [in] */ BSTR pFreqRSF,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pResetTimingSF,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pDayCountSF,
            /* [in] */ double pSpdSF,
            /* [in] */ VARIANT pSpdfund,
            /* [in] */ double pResetGap,
            /* [in] */ double pResetGapSF,
            /* [in] */ BSTR pAmort,
            /* [in] */ double pRefDate,
            /* [in] */ double pFee,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DefaultIntensity( 
            /* [in] */ BSTR pricerId,
            /* [in] */ double Maturity,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMRiskyBond( 
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedemptionPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double pRepo,
            /* [defaultvalue][in][optional] */ double pSsl,
            /* [defaultvalue][in][optional] */ double pRecoveryRate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMRiskyBondWithCF( 
            /* [in] */ double pAsOfDate,
            /* [in] */ double pRedemptionPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [in] */ VARIANT __RPC_FAR *pYearTerms,
            /* [in] */ VARIANT __RPC_FAR *pCashFlows,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double pRepo,
            /* [defaultvalue][in][optional] */ double pSsl,
            /* [defaultvalue][in][optional] */ double pRecoveryRate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_EmptyLeg( 
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMClonedAndSetNotional( 
            /* [in] */ BSTR bLegId,
            /* [defaultvalue][in][optional] */ BSTR bAmortId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_INF_GetZcFromSummit( 
            /* [in] */ BSTR Index,
            /* [in] */ BSTR Ccy,
            /* [in] */ BSTR cvname,
            /* [in] */ double date,
            /* [defaultvalue][in][optional] */ BSTR seasonAdj,
            /* [defaultvalue][in][optional] */ BSTR seasonAdjMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_IRLEGTOCREDITLEG( 
            /* [in] */ BSTR SwapLegId,
            /* [in] */ BSTR LegType,
            /* [in] */ BSTR creditindexId,
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Collateral( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDSGEN( 
            /* [in] */ BSTR FeeLegId,
            /* [in] */ BSTR DefLegId,
            /* [in] */ double RcvFee,
            /* [in] */ double TradedNot,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_NTDGEN( 
            /* [in] */ BSTR CdsId,
            /* [in] */ int firstnumdef,
            /* [in] */ int lastnumdef,
            /* [in] */ BSTR CollateralId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDOGEN( 
            /* [in] */ BSTR CdsId,
            /* [in] */ double subamount,
            /* [in] */ BSTR CollateralId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenLeg( 
            /* [in] */ double StartDate,
            /* [in] */ double EndDate,
            /* [in] */ double FixedRate,
            /* [in] */ double FixedNotional,
            /* [defaultvalue][in][optional] */ BSTR RefValNotional,
            /* [defaultvalue][in][optional] */ BSTR RefValRate,
            /* [defaultvalue][in][optional] */ BSTR XChangeNotional,
            /* [defaultvalue][in][optional] */ BSTR Frequency,
            /* [defaultvalue][in][optional] */ BSTR Basis,
            /* [defaultvalue][in][optional] */ BSTR payTiming,
            /* [defaultvalue][in][optional] */ BSTR intrule,
            /* [defaultvalue][in][optional] */ BSTR stubrule,
            /* [defaultvalue][in][optional] */ BSTR ccyid,
            /* [defaultvalue][in][optional] */ BSTR paycalname,
            /* [defaultvalue][in][optional] */ double refdate,
            /* [defaultvalue][in][optional] */ BSTR includematurity,
            /* [defaultvalue][in][optional] */ BSTR adjstartdate,
            /* [defaultvalue][in][optional] */ BSTR legtype,
            /* [defaultvalue][in][optional] */ BSTR indexobj,
            /* [defaultvalue][in][optional] */ int creditlag,
            /* [defaultvalue][in][optional] */ double binary,
            /* [defaultvalue][in][optional] */ BSTR name,
            /* [defaultvalue][in][optional] */ BSTR Nxchange,
            /* [defaultvalue][in][optional] */ BSTR baccruedOnDef,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CDO_SQUARE_GEN( 
            /* [in] */ BSTR CdsId,
            /* [in] */ double subamount,
            /* [in] */ BSTR portfolioId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetLastFixingDate( 
            /* [in] */ BSTR instId,
            /* [in] */ VARIANT asofDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetPastFixing( 
            /* [in] */ BSTR instId,
            /* [in] */ VARIANT resetDate,
            /* [in] */ VARIANT fixingValue,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFixing( 
            /* [in] */ BSTR source,
            /* [in] */ BSTR index,
            /* [in] */ BSTR term,
            /* [in] */ BSTR ccy,
            /* [in] */ double date,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetRiskyProfile( 
            /* [in] */ BSTR CdsorCdoId,
            /* [in] */ BSTR CouponsId,
            /* [in] */ BSTR TypesId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMHyperCube( 
            /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
            /* [in] */ VARIANT __RPC_FAR *pKeys,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetPricerForRatesComputation( 
            /* [in] */ BSTR legId,
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCmsLeg( 
            /* [in] */ double startDate,
            /* [in] */ double endDate,
            /* [in] */ BSTR cmsTypeId,
            /* [defaultvalue][in][optional] */ BSTR receiveOrPay,
            /* [defaultvalue][in][optional] */ BSTR yieldDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR swapLegDayCount,
            /* [defaultvalue][in][optional] */ BSTR resetFreq,
            /* [defaultvalue][in][optional] */ BSTR payFreq,
            /* [defaultvalue][in][optional] */ long resetGap,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR ccyName,
            /* [defaultvalue][in][optional] */ BSTR resetTiming,
            /* [defaultvalue][in][optional] */ BSTR stubRule,
            /* [defaultvalue][in][optional] */ BSTR adjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetMatuLabel( 
            /* [in] */ BSTR pCurveId,
            /* [in] */ VARIANT __RPC_FAR *pMatuLabels,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputePentifix( 
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ BSTR pFloorOrCap,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pCorrManager,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
            /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_ReplicConvAdjust_Create( 
            /* [in] */ BSTR Payoff_ReplicMode,
            /* [in] */ double Payoff_StepOrReplicPrecision,
            /* [in] */ BSTR Payoff_StopMode,
            /* [in] */ double Payoff_StopThreshold,
            /* [in] */ BSTR Sensi_ReplicMode,
            /* [in] */ double Sensi_StepOrReplicPrecision,
            /* [in] */ BSTR Sensi_StopMode,
            /* [in] */ double Sensi_StopThreshold,
            /* [in] */ BSTR UsedModelId,
            /* [in] */ double StrikeMinReplic,
            /* [in] */ double StrikeMaxReplic,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_MapConvAdjust_Create( 
            /* [in] */ BSTR LiborArrearAdj,
            /* [in] */ BSTR NaturalCMSAdj,
            /* [in] */ BSTR PaymentLagAdj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetFees( 
            /* [in] */ BSTR securityId,
            /* [in] */ BSTR RefvalueId) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetBounds( 
            /* [in] */ BSTR securityId,
            /* [out] */ double __RPC_FAR *down,
            /* [out] */ double __RPC_FAR *up) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputePentibor( 
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPay,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ BSTR pIndexIndexCorrelCube,
            /* [in] */ BSTR pCorrEUR,
            /* [in] */ BSTR pInterCorr,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetRecovCoef( 
            /* [in] */ BSTR pCurveId,
            /* [in] */ double RecovCoef,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMIndexIndexCorrelCube( 
            /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
            /* [in] */ VARIANT __RPC_FAR *pTenors1List,
            /* [in] */ VARIANT __RPC_FAR *pTenors2List,
            /* [defaultvalue][in][optional] */ BSTR pInterSurfInterp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetInterpolationType( 
            /* [in] */ VARIANT __RPC_FAR *pVolCurveId,
            /* [in] */ VARIANT __RPC_FAR *pInterpolType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Customized_CDO( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [defaultvalue][in] */ BSTR Currency,
            /* [in] */ BSTR pDefaultLeg,
            /* [in] */ BSTR pPremiumLeg,
            /* [in] */ BSTR pParameters,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCreateGenCorrelatorManager( 
            /* [in] */ VARIANT __RPC_FAR *pMktTags,
            /* [in] */ VARIANT __RPC_FAR *pHyperDiagVol,
            /* [in] */ VARIANT __RPC_FAR *pIndexIndexVol,
            /* [in] */ VARIANT __RPC_FAR *pCorrelVol,
            /* [in] */ VARIANT __RPC_FAR *pIndexVol,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CLN( 
            /* [in] */ double EffectiveDate,
            /* [in] */ double EndDate,
            /* [in] */ double Spread,
            /* [defaultvalue][in][optional] */ BSTR IndexId,
            /* [defaultvalue][in][optional] */ double Refdate,
            /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ double Notional,
            /* [defaultvalue][in][optional] */ BSTR AccOnDef,
            /* [defaultvalue][in][optional] */ BSTR DayCount,
            /* [defaultvalue][in][optional] */ BSTR DecompFreq,
            /* [defaultvalue][in][optional] */ BSTR StubRule,
            /* [defaultvalue][in][optional] */ double resetgap,
            /* [defaultvalue][in][optional] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ BSTR ResetCal,
            /* [defaultvalue][in][optional] */ BSTR PayCal,
            /* [defaultvalue][in][optional] */ BSTR Nxchange,
            /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
            /* [defaultvalue][in][optional] */ BSTR AdjustedStartDate,
            /* [defaultvalue][in][optional] */ double Binary,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputePentilix( 
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPay,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ BSTR pIndexIndexCorrelCube,
            /* [in] */ BSTR pCorrEUR,
            /* [in] */ BSTR pInterCorr,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_RiskyPV01( 
            /* [in] */ BSTR pDefCurve,
            /* [optional][in] */ VARIANT Date1,
            /* [in] */ VARIANT Date2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetResetMgrFromSummit( 
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pSource,
            /* [defaultvalue][in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pIsInflationIndex,
            /* [defaultvalue][in] */ BSTR pTerm,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetReset( 
            /* [in] */ BSTR pResetMgr,
            /* [in] */ double pDate,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSetLastFixing( 
            /* [in] */ BSTR pSecurityId,
            /* [in] */ double pRate,
            /* [in] */ double pAsOf,
            /* [defaultvalue][in][optional] */ double pBeforeLastFixingDate,
            /* [defaultvalue][in][optional] */ double pResetDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Sectorial_Correlation( 
            /* [in] */ DATE AsOf,
            /* [in] */ BSTR structName,
            /* [in] */ BSTR correlation_Type,
            /* [in] */ VARIANT __RPC_FAR *vLabels,
            /* [in] */ VARIANT __RPC_FAR *vector_Membership,
            /* [defaultvalue][in][optional] */ double intra_Sector_Correlation,
            /* [defaultvalue][in][optional] */ double inter_Sector_Correlation,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas_Down,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas_Down,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeReviPentix( 
            /* [in] */ double pNtl,
            /* [in] */ VARIANT pDate,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ VARIANT pSpread,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ BSTR pFloorOrCap,
            /* [in] */ double pFee,
            /* [in] */ double pLevier,
            /* [in] */ double pTxFixeMax,
            /* [in] */ BSTR pIsCapped,
            /* [in] */ double pTxCap,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pCorrManager,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
            /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGenAmortization( 
            /* [in] */ BSTR pSwaplegId,
            /* [in] */ BSTR pAmortMethod,
            /* [defaultvalue][in][optional] */ BSTR pAmortFreq,
            /* [defaultvalue][in][optional] */ double pAmortAmount,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ double pLegNotional,
            /* [defaultvalue][in][optional] */ double pAmortRate,
            /* [defaultvalue][in][optional] */ double pReducedMaturity,
            /* [defaultvalue][in][optional] */ BSTR pModelId,
            /* [defaultvalue][in][optional] */ double pCleanUp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCptRefvalue( 
            /* [in] */ BSTR pRefValId,
            /* [in] */ double pDate,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetDPFromCalypso( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR irCurveId,
            /* [defaultvalue][in][optional] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *ret) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCalypsoDevConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetZCFromCalypso( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pTerm,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DefProbInverse( 
            /* [in] */ BSTR pCurveId,
            /* [in] */ double dDefProba,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DPMktDataFromCalypso( 
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [in] */ VARIANT __RPC_FAR *Parameter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ZeroCouponDefaultCurveFromCalypso( 
            /* [in] */ double pDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR irCurveId,
            /* [defaultvalue][in][optional] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *ret) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetInitialCurveFromCalypso( 
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pTerm,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR pDoAdj,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetRate) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCalypsoProdConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMCalypsoRecConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetBasketCorrelMkDataFromCalypso( 
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double date,
            /* [in] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFileName,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_QMatrix( 
            /* [in] */ VARIANT __RPC_FAR *pQMatrix,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_MarketDataMng( 
            /* [in] */ VARIANT __RPC_FAR *pstrVect,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ModelMultiCvMktDataMng( 
            /* [in] */ VARIANT __RPC_FAR *pIRcurve,
            /* [in] */ VARIANT __RPC_FAR *pDefCurves,
            /* [in] */ VARIANT __RPC_FAR *pRecovery,
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR MktdataMngId,
            /* [defaultvalue][in][optional] */ BSTR pVolcurve,
            /* [defaultvalue][in][optional] */ BSTR cloneorNot,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE Local_ARM_ProdConnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_SetDefaultCurrency( 
            /* [in] */ BSTR isoCCy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLivretALeg( 
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pRcvOrPay,
            /* [defaultvalue][in][optional] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CorrelationSmileStrike( 
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pVolCurves,
            /* [in] */ VARIANT __RPC_FAR *pProportions,
            /* [in] */ double AsOf,
            /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeLow,
            /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
            /* [in][optional] */ VARIANT __RPC_FAR *pIndexVector,
            /* [defaultvalue][in][optional] */ BSTR Name,
            /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeLow,
            /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeUp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMFutDelivery( 
            /* [in] */ BSTR pFut,
            /* [in] */ BSTR pCcy,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMLivretACurve( 
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pInfCurvId,
            /* [in] */ BSTR pEuribCurvId,
            /* [defaultvalue][in][optional] */ double pFlagRouding,
            /* [defaultvalue][in][optional] */ BSTR pInfResetMgrId,
            /* [defaultvalue][in][optional] */ BSTR pFixingLivretAId,
            /* [defaultvalue][in][optional] */ BSTR pFixingEuribId,
            /* [defaultvalue][in][optional] */ BSTR pMonthForAugust,
            /* [defaultvalue][in][optional] */ BSTR pMonthForFebruary,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeLivretA( 
            /* [in] */ double pNtl,
            /* [in] */ double pStartDateLeg1,
            /* [in] */ double pEndDateLeg1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexLeg1,
            /* [in] */ VARIANT pSpreadLeg1,
            /* [in] */ BSTR pDayCountLeg1,
            /* [in] */ BSTR pPayFreqLeg1,
            /* [in] */ BSTR pResetFreqLeg1,
            /* [in] */ BSTR pResetTimingLeg1,
            /* [in] */ BSTR pAdjLeg1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pStub,
            /* [in] */ double pEndDateLA,
            /* [in] */ double pSpreadLeg2,
            /* [in] */ BSTR pDayCountLA,
            /* [in] */ BSTR pPayFreqLA,
            /* [in] */ BSTR pResetFreqLA,
            /* [in] */ BSTR pResetTimingLA,
            /* [in] */ BSTR pAdjLA,
            /* [in] */ BSTR pIndexLeg2,
            /* [in] */ BSTR pDayCountLeg2,
            /* [in] */ BSTR pPayFreqLeg2,
            /* [in] */ BSTR pResetFreqLeg2,
            /* [in] */ BSTR pResetTimingLeg2,
            /* [in] */ BSTR pAdjLeg2,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pLAMod,
            /* [in] */ BSTR pLAModBump,
            /* [in] */ BSTR pLAModBumpInflation,
            /* [in] */ VARIANT __RPC_FAR *pResetMgrIds,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMGetFixingFromCalypso( 
            /* [in] */ BSTR source,
            /* [in] */ BSTR index,
            /* [in] */ BSTR term,
            /* [in] */ BSTR ccy,
            /* [in] */ BSTR curveName,
            /* [in] */ DATE date,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_FxConvertFromCalypso( 
            /* [in] */ BSTR ccy1,
            /* [in] */ BSTR ccy2,
            /* [defaultvalue][in][optional] */ BSTR pCvName,
            /* [in] */ DATE pDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FwdSpread( 
            /* [in] */ BSTR defcurveId,
            /* [in] */ double Maturity1,
            /* [in] */ double Maturity2,
            /* [defaultvalue][in][optional] */ double FwdStartDate,
            /* [defaultvalue][in][optional] */ double FwdEndDate,
            /* [defaultvalue][in][optional] */ BSTR VolId,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Flat_Correlation( 
            /* [in] */ DATE AsOf,
            /* [in] */ BSTR structName,
            /* [in] */ double correlValue,
            /* [defaultvalue][in] */ BSTR idIndex1,
            /* [defaultvalue][in] */ BSTR idIndex2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CorridorLeg_Sche( 
            /* [in] */ double Notional,
            /* [in] */ BSTR RecieveOrPay,
            /* [in] */ BSTR RefValueSpreadsInBP,
            /* [in] */ BSTR floatingIdx,
            /* [in] */ double leverageFloatIdx,
            /* [in] */ BSTR creditIdx,
            /* [in] */ BSTR refvalueKUPinBP,
            /* [in] */ BSTR refvalueKDWinBP,
            /* [in] */ BSTR ScheduleInfoId,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *accondef,
            /* [defaultvalue][in][optional] */ BSTR disc_ccy,
            /* [defaultvalue][in][optional] */ BSTR Name,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Schedule_Info( 
            /* [in] */ double EffectiveDate,
            /* [in] */ double MaturityDate,
            /* [defaultvalue][in][optional] */ BSTR payFrequency,
            /* [defaultvalue][in][optional] */ BSTR ResetFreq,
            /* [defaultvalue][in][optional] */ BSTR DayCount,
            /* [defaultvalue][in][optional] */ BSTR Stubrule,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR payCalName,
            /* [defaultvalue][in][optional] */ BSTR PayTiming,
            /* [defaultvalue][in][optional] */ BSTR ResetTiming,
            /* [defaultvalue][in][optional] */ BSTR fwdRule,
            /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
            /* [defaultvalue][in][optional] */ BSTR adj,
            /* [defaultvalue][in][optional] */ BSTR intStartAdj,
            /* [defaultvalue][in][optional] */ BSTR AccDayCount,
            /* [defaultvalue][in][optional] */ double ReferenceDate,
            /* [defaultvalue][in][optional] */ double FirstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR AdjCal,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMInfCurveSetResetMgr( 
            /* [in] */ BSTR pInfCurve,
            /* [in] */ BSTR pResetMgr,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CPDO( 
            /* [in] */ BSTR pRiskyLeg,
            /* [in] */ BSTR pRollLeg,
            /* [in] */ BSTR pNoRiskyLeg,
            /* [in] */ double pInitialValo,
            /* [in] */ double pTarget,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCpnType,
            /* [in] */ double pUFFees,
            /* [in] */ double pRunningFees,
            /* [in] */ double pVExpo,
            /* [in] */ double pV0Expo,
            /* [in] */ double pAlpha,
            /* [in] */ double pBeta,
            /* [in] */ double pDesactivation,
            /* [in] */ int pNbAssets,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_PriceVector( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR pCPTTYPE,
            /* [out] */ VARIANT __RPC_FAR *pRetVectorValos) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenPrice( 
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR pCPTTYPE,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeTxFixed( 
            /* [in] */ double pNtl,
            /* [in] */ double pStartDateLeg1,
            /* [in] */ double pEndDateLeg1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexLeg1,
            /* [in] */ VARIANT pSpreadLeg1,
            /* [in] */ BSTR pDayCountLeg1,
            /* [in] */ BSTR pPayFreqLeg1,
            /* [in] */ BSTR pResetFreqLeg1,
            /* [in] */ BSTR pResetTimingLeg1,
            /* [in] */ BSTR pAdjLeg1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pStub,
            /* [in] */ double pEndDateFixed,
            /* [in] */ double pSpreadLeg2,
            /* [in] */ BSTR pDayCountFixed,
            /* [in] */ BSTR pPayFreqFixed,
            /* [in] */ BSTR pResetFreqFixed,
            /* [in] */ BSTR pResetTimingFixed,
            /* [in] */ BSTR pAdjFixed,
            /* [in] */ BSTR pIndexLeg2,
            /* [in] */ BSTR pDayCountLeg2,
            /* [in] */ BSTR pPayFreqLeg2,
            /* [in] */ BSTR pResetFreqLeg2,
            /* [in] */ BSTR pResetTimingLeg2,
            /* [in] */ BSTR pAdjLeg2,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FixingCurve( 
            /* [in] */ VARIANT __RPC_FAR *pDates,
            /* [in] */ VARIANT __RPC_FAR *pValues,
            /* [defaultvalue][in][optional] */ double AsOfDate,
            /* [defaultvalue][in][optional] */ BSTR B_IndexName,
            /* [defaultvalue][in][optional] */ BSTR B_IndexID,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CptInterpolDefCurveOLD( 
            /* [in] */ BSTR pCurve,
            /* [in] */ BSTR pTenor,
            /* [in] */ double pSlope,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ double pInterpDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_CreateBasketCorrelMkDataFromCalypso( 
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double date,
            /* [in] */ BSTR forceCurveName,
            /* [in] */ BSTR Ccy,
            /* [defaultvalue][in][optional] */ BSTR xmlFileName,
            /* [defaultvalue][in][optional] */ BSTR indexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSetCalendar( 
            /* [in] */ BSTR pFileName,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMInitGigaSpaces( 
            /* [in] */ BSTR pUrl,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GetExpectedLoss( 
            /* [in] */ BSTR pricerId,
            /* [in] */ double YearTerm,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_VariableCollateral( 
            /* [in] */ VARIANT pLabels,
            /* [in] */ VARIANT pNotionals,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_IndexCompo( 
            /* [in] */ BSTR IndexName,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *YearFrac,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *pSpread,
            /* [defaultvalue][in] */ BSTR Method,
            /* [defaultvalue][in] */ BSTR Basis,
            /* [defaultvalue][in] */ BSTR ResetFreq,
            /* [defaultvalue][in] */ BSTR PayFreq,
            /* [defaultvalue][in] */ BSTR ccy,
            /* [defaultvalue][in] */ BSTR fwdRule,
            /* [defaultvalue][in] */ BSTR resetTiming,
            /* [defaultvalue][in] */ int resetGap,
            /* [defaultvalue][in] */ BSTR payTiming,
            /* [defaultvalue][in] */ int payGap,
            /* [defaultvalue][in] */ BSTR intRule,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ int cm_resetWeekDay,
            /* [defaultvalue][in] */ int cm_resetOccur,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_createFlatCurve( 
            /* [in] */ BSTR pCurve,
            /* [in] */ VARIANT pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FunctionRegister( 
            /* [in] */ long address) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBermudanXStyle( 
            /* [in] */ VARIANT __RPC_FAR *pxDates,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pexpiryDates,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMcomputeCRA( 
            /* [in] */ BSTR pFixorFloat,
            /* [in] */ double pFee,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pCcy,
            /* [in] */ double pLevelUp,
            /* [in] */ BSTR pUpSpec,
            /* [in] */ double pLevelDown,
            /* [in] */ BSTR pDownSpec,
            /* [in] */ BSTR pRefIndex,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pPayFreqPayIndex,
            /* [in] */ BSTR pResetFreqRefIndex,
            /* [in] */ BSTR pPaidRstTiming,
            /* [in] */ BSTR pRefRstTiming,
            /* [in] */ BSTR pStubRule,
            /* [in] */ BSTR pPOrR,
            /* [in] */ double pStartCallDate,
            /* [in] */ BSTR pXStyle,
            /* [in] */ BSTR pFundingIndex,
            /* [in] */ BSTR pResetFreqFunding,
            /* [in] */ BSTR pPayFreqFunding,
            /* [in] */ VARIANT pSpreadFunding,
            /* [in] */ BSTR pPOrRFunding,
            /* [in] */ double pDecompPricingFlag,
            /* [in] */ double pdiscMarginFactor,
            /* [in] */ BSTR pPreInitFlag,
            /* [in] */ double pMeanReversion,
            /* [in] */ VARIANT pCalibParams,
            /* [in] */ VARIANT pCalibParamsPF,
            /* [in] */ double pKernelToGP,
            /* [in] */ VARIANT pMarkovTreeParams,
            /* [in] */ double pMarkovTreePathNumber,
            /* [in] */ BSTR pBsmodId,
            /* [in] */ BSTR pBsmodSwoptId,
            /* [in] */ BSTR pBsmodSwoptBumpId,
            /* [in] */ BSTR pzcId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Math_BivNormale( 
            /* [in] */ double x,
            /* [in] */ double y,
            /* [in] */ double rho,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_GETINSTRUMENTFROMCALYPSO( 
            /* [in] */ BSTR CalypsoId,
            /* [in] */ BSTR Type,
            /* [defaultvalue][optional][in] */ double AsOf,
            /* [defaultvalue][optional][in] */ BSTR ModelType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMSetDiscountPricingMode( 
            /* [in] */ BSTR pModelId,
            /* [in] */ int pDiscountPricingMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Math_RandUniform( 
            /* [defaultvalue][optional][in] */ double seed,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Math_Interpol( 
            /* [in] */ VARIANT __RPC_FAR *X,
            /* [in] */ VARIANT __RPC_FAR *Y,
            /* [in] */ double value,
            /* [defaultvalue][in][optional] */ double type,
            /* [defaultvalue][in][optional] */ double smooth,
            /* [in][optional] */ VARIANT __RPC_FAR *Weights,
            /* [defaultvalue][in][optional] */ double modeSpline,
            /* [defaultvalue][in][optional] */ double withC1condition,
            /* [defaultvalue][in][optional] */ double leftSlope,
            /* [defaultvalue][in][optional] */ double rightSlope,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Random_Generator( 
            /* [in] */ BSTR RandomType,
            /* [defaultvalue][optional][in] */ BSTR ParamId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenerateOneRandom( 
            /* [in] */ BSTR RandomId,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_GenerateRandoms( 
            /* [in] */ BSTR RandomId,
            /* [in] */ int DimVector,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_ResetRandom( 
            /* [in] */ BSTR RandomId) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_FwdSpreadAsIndex( 
            /* [in] */ BSTR DefCurveId,
            /* [in] */ double matu1,
            /* [in] */ double Matu2,
            /* [retval][out] */ double __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_createDefCurveFromBase( 
            /* [in] */ BSTR pCurveCDS,
            /* [in] */ BSTR pCurveIndex,
            /* [in] */ VARIANT vBase,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_RiskyPV01AsSensitivity( 
            /* [in] */ BSTR pDefCurve,
            /* [in] */ BSTR Tenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_SetVolCurve( 
            /* [in] */ BSTR Model,
            /* [in] */ BSTR VolCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMPF( 
            /* [in] */ VARIANT __RPC_FAR *pinsts,
            /* [in] */ VARIANT __RPC_FAR *pcoeffs,
            /* [in] */ VARIANT __RPC_FAR *pmarketPrices,
            VARIANT __RPC_FAR *pprecisions,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMBondTEC( 
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedempPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double ptec,
            /* [defaultvalue][in][optional] */ BSTR pPFTecId,
            /* [defaultvalue][in][optional] */ BSTR pModTecId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMPFModFit( 
            /* [in] */ BSTR pmodName,
            /* [in] */ BSTR ppf,
            /* [in] */ double psettlement,
            /* [in] */ BSTR pzc,
            /* [in] */ VARIANT __RPC_FAR *pvList,
            /* [in] */ VARIANT __RPC_FAR *pfList,
            /* [defaultvalue][in][optional] */ long nag_algo,
            /* [defaultvalue][in][optional] */ long pstep,
            /* [defaultvalue][in][optional] */ double phorizon,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_Restrikable_CDO( 
            /* [in] */ double UnderlyingMatu,
            /* [in] */ double Expiry,
            /* [in] */ double Strike,
            /* [in] */ int OptionType,
            /* [in] */ BSTR pUnderlying,
            /* [in] */ double Rehauss,
            /* [in] */ double TriggerFreq,
            /* [in] */ int DiffCDO,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_PropertyList( 
            /* [in] */ VARIANT attrNames,
            /* [in] */ VARIANT attrValues,
            /* [in][optional] */ VARIANT attrTypes,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARM_Credit_DefCurveIntensityPWC( 
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pMatuRates,
            /* [in] */ VARIANT __RPC_FAR *pInputs,
            /* [in] */ double Type,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [defaultvalue][in][optional] */ BSTR bCurrency,
            /* [defaultvalue][in][optional] */ BSTR bLabel,
            /* [defaultvalue][in][optional] */ BSTR VolCurveId,
            /* [defaultvalue][in][optional] */ BSTR calibrationAlgo,
            /* [defaultvalue][in][optional] */ int lag,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ARMTMLeg( 
            /* [in] */ BSTR ptmIxType,
            /* [in] */ double pstartDate,
            /* [in] */ double pendDate,
            /* [in] */ BSTR pPorR,
            /* [defaultvalue][in][optional] */ double pspread,
            /* [defaultvalue][in][optional] */ BSTR ppayFrequency,
            /* [defaultvalue][in][optional] */ BSTR presetFrequency,
            /* [defaultvalue][in][optional] */ BSTR pinterestRule,
            /* [defaultvalue][in][optional] */ BSTR pfwdRule,
            /* [defaultvalue][in][optional] */ BSTR pstubRule,
            /* [in] */ BSTR pccy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IARMModuleVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IARMModule __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IARMModule __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTypeInfoCount )( 
            IARMModule __RPC_FAR * This,
            /* [out] */ UINT __RPC_FAR *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTypeInfo )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo __RPC_FAR *__RPC_FAR *ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetIDsOfNames )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR __RPC_FAR *rgszNames,
            /* [in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID __RPC_FAR *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Invoke )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ DISPID dispIdMember,
            /* [in] */ REFIID riid,
            /* [in] */ LCID lcid,
            /* [in] */ WORD wFlags,
            /* [out][in] */ DISPPARAMS __RPC_FAR *pDispParams,
            /* [out] */ VARIANT __RPC_FAR *pVarResult,
            /* [out] */ EXCEPINFO __RPC_FAR *pExcepInfo,
            /* [out] */ UINT __RPC_FAR *puArgErr);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMNextBusinessDay )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [in] */ long pNbDays,
            /* [retval][out] */ double __RPC_FAR *pDate2);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMAdjustToBusDate )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [in] */ BSTR pRule,
            /* [retval][out] */ double __RPC_FAR *pDate2);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFreeObject )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pId,
            /* [retval][out] */ long __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMIsBusinessDay )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCalendrier,
            /* [retval][out] */ long __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetZCFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFreeAllObjects )( 
            IARMModule __RPC_FAR * This,
            /* [retval][out] */ long __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMYcMod )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZc,
            /* [defaultvalue][in][optional] */ BSTR pZcDiscount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMForwardYield )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZc,
            /* [in] */ double pMatu1,
            /* [in] */ double pMatu2,
            /* [in] */ BSTR pMeth,
            /* [in] */ BSTR pAdj,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMDiscountYield )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pZc,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pMeth,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLiborSwap )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pStartDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pLiborType,
            /* [in] */ VARIANT __RPC_FAR *pRecOrPay,
            /* [in] */ VARIANT __RPC_FAR *pFixedRate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pFloatingDaycount,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwapPriceToRate )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pSwap,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pPrice,
            /* [in] */ VARIANT __RPC_FAR *pModel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pSec,
            /* [in] */ VARIANT __RPC_FAR *pModel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBetweenDates )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pDate1,
            /* [in] */ VARIANT __RPC_FAR *pDate2,
            /* [in] */ VARIANT __RPC_FAR *pDaycount,
            /* [in] */ VARIANT __RPC_FAR *pIsYearFrac,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMAddPeriod )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pFreq,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [in] */ VARIANT __RPC_FAR *pNbPeriods,
            /* [in] */ VARIANT __RPC_FAR *pAdjRule,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMIsoCcy )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetSpotDays )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetLiborIndexDaycount )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetLiborTerm )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFixedDayCount )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFixedPayFreq )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMComputeVolatility )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pVol,
            /* [in] */ VARIANT __RPC_FAR *pmatu,
            /* [in] */ VARIANT __RPC_FAR *pStrike,
            /* [in] */ VARIANT __RPC_FAR *pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMVolCurv )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pStrike,
            /* [in] */ VARIANT __RPC_FAR *pVol,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pStrikeType,
            /* [in] */ BSTR pVolType,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIndexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetVolCubeFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [in] */ BSTR pType,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pSmiles,
            /* [defaultvalue][in][optional] */ BSTR pTypeCube,
            /* [defaultvalue][in][optional] */ BSTR indexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_FxConvert )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pccy1,
            /* [in] */ VARIANT __RPC_FAR *pccy2,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pCvName,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_DiscountPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pcurve,
            /* [in] */ VARIANT __RPC_FAR *pmatu,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Delivery )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pAsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pTenorContract,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CptInterpolDefCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurve,
            /* [in] */ VARIANT pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DefaultProba )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCurve,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetBeta )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Mezzanine )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [in] */ double pMezzAmount,
            /* [in] */ double pSubAmount,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [defaultvalue][in] */ BSTR pFreqFeeLeg,
            /* [defaultvalue][in] */ BSTR pDayCount,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in] */ BSTR Currency,
            /* [defaultvalue][in] */ double pPayCreditLag,
            /* [defaultvalue][in] */ BSTR pStub,
            /* [defaultvalue][in] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in] */ double pBinary,
            /* [defaultvalue][in] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR adjStartDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ModelMultiCurves )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pIRcurve,
            /* [in] */ VARIANT __RPC_FAR *pDefCurves,
            /* [in] */ VARIANT __RPC_FAR *pRecoveryRates,
            /* [defaultvalue][in][optional] */ BSTR CorrelId,
            /* [defaultvalue][in][optional] */ BSTR VolCurve,
            /* [defaultvalue][in][optional] */ BSTR CpnInfCurve,
            /* [defaultvalue][in][optional] */ BSTR CpnIRCurve,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FTD )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pIssuerNotional,
            /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
            /* [in] */ VARIANT __RPC_FAR *pStub,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pintRule,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pstartAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_NTD )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [defaultvalue][in] */ int pFirstNumDefault,
            /* [defaultvalue][in] */ int pLastNumDefault,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ BSTR pFixingFreq,
            /* [defaultvalue][in] */ BSTR pDayCountFrq,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ double pIssuerNotional,
            /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in] */ BSTR pCurrency,
            /* [defaultvalue][in] */ double pPayCreditLag,
            /* [defaultvalue][in] */ BSTR pStub,
            /* [defaultvalue][in] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in] */ double pBinary,
            /* [defaultvalue][in] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR startAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Price )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pAsofDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_RiskyDuration )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDefCurve,
            /* [in] */ VARIANT date,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDONPV )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pCPTTYPE,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CorrMatrix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [defaultvalue][in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ExtractCorrMatrix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCorrMatrixId,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Spread )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pMTM,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetLabel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCurveId,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetLabel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Sensitivity )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [in] */ VARIANT __RPC_FAR *pType,
            /* [in] */ VARIANT __RPC_FAR *pPlot,
            /* [in] */ VARIANT __RPC_FAR *pLabel,
            /* [in] */ VARIANT __RPC_FAR *pEpsilon,
            /* [defaultvalue][in][optional] */ double epsilonGamma,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetCleanSpreadTranche )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [in] */ VARIANT __RPC_FAR *pPlot,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetDefProbTranche )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR PricerId,
            /* [in] */ double Yearterm,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetDuration )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricerId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenSchedule )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pAccStartDate,
            /* [in] */ VARIANT __RPC_FAR *pAccEndDate,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *prefDate,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *ptypeDates,
            /* [in] */ VARIANT __RPC_FAR *pModFol,
            /* [in] */ VARIANT __RPC_FAR *pCreditGap,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CashFlows )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_View )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pObjet,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Version )( 
            IARMModule __RPC_FAR * This,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDS )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pSpread,
            /* [defaultvalue][in] */ BSTR pFixingFreq,
            /* [defaultvalue][in] */ BSTR pDayCountFrq,
            /* [defaultvalue][in] */ double pFirst_period_refdate,
            /* [defaultvalue][in] */ double pFixedPayerAmount,
            /* [defaultvalue][in] */ double pFloatingPayerAmount,
            /* [defaultvalue][in] */ BSTR StubRule,
            /* [defaultvalue][in] */ BSTR pCurrency,
            /* [defaultvalue][in] */ BSTR Adjusted,
            /* [defaultvalue][in] */ int CreditDefLag,
            /* [defaultvalue][in] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double StartProtection,
            /* [defaultvalue][in] */ double EndProtection,
            /* [defaultvalue][in] */ BSTR name,
            /* [defaultvalue][in] */ double binary,
            /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR StartAdj,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetVolFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pCcy,
            /* [in] */ VARIANT __RPC_FAR *pCvName,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *pType,
            /* [in] */ VARIANT __RPC_FAR *pMatIndex,
            /* [in] */ VARIANT __RPC_FAR *pImpOrHist,
            /* [defaultvalue][in][optional] */ BSTR pindexId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMParallelShift )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZc,
            /* [in] */ double pBump,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBumpVolatility )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pVolCurv,
            /* [in] */ double pBump,
            /* [in] */ long pNthLine,
            /* [in] */ long pNthCol,
            /* [in] */ BSTR pIsCumul,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBsSmiledModel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ double pSpot,
            /* [in] */ BSTR pDividend,
            /* [in] */ BSTR pDiscrate,
            /* [in] */ BSTR pVolATM,
            /* [in] */ BSTR pRo,
            /* [in] */ BSTR pNu,
            /* [in] */ BSTR pIsSABR,
            /* [defaultvalue][in][optional] */ BSTR pBeta,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSetEtoolkit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pUserName,
            /* [in] */ BSTR pPassWord,
            /* [in] */ BSTR pDatabaseContext,
            /* [in] */ BSTR pItConfigDomainDir,
            /* [in] */ BSTR pItDomainName,
            /* [retval][out] */ long __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMConnectionEtoolkit )( 
            IARMModule __RPC_FAR * This,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMVolFlat )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pVol,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMVolCube )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pATMVol,
            /* [in] */ VARIANT __RPC_FAR *pSmileCurveIds,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [defaultvalue][in] */ BSTR pVolType,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMDeconnectionEtoolkit )( 
            IARMModule __RPC_FAR * This,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMZcFlat )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pZc,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBsModel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ double pSpot,
            /* [in] */ BSTR pDividend,
            /* [in] */ BSTR pDiscrate,
            /* [in] */ BSTR pVol,
            /* [defaultvalue][in][optional] */ BSTR pTypeStk,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwitchToETK )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwitchToFLATFILE )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMZCLINT )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [in] */ BSTR pMeth,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCurrency,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Pricer )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSecurity,
            /* [in] */ BSTR pModel,
            /* [in] */ BSTR pPricerType,
            /* [in] */ int l_nbpaths,
            /* [defaultvalue][in][optional] */ BSTR pParameters,
            /* [defaultvalue][in][optional] */ double valuationdate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_zcspreaded )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR zcSprId,
            /* [in] */ BSTR zcInitId,
            /* [in] */ double date,
            /* [in] */ BSTR MMFreq,
            /* [in] */ BSTR SwapFreq,
            /* [in] */ BSTR ccyId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CptBaseCorrelation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOf,
            /* [in] */ BSTR name,
            /* [in] */ BSTR CalMethod,
            /* [in] */ BSTR IndexId,
            /* [in] */ VARIANT __RPC_FAR *pStrikeLow,
            /* [in] */ VARIANT __RPC_FAR *pStrikeHigh,
            /* [in] */ VARIANT __RPC_FAR *pVMktBid,
            /* [in] */ VARIANT __RPC_FAR *pVMktAsk,
            /* [in] */ VARIANT __RPC_FAR *pVUpfBid,
            /* [in] */ VARIANT __RPC_FAR *pVUpfAsk,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVInitialCorrel,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVDeltaLevrage,
            /* [defaultvalue][in][optional] */ BSTR ModelId,
            /* [defaultvalue][in][optional] */ double integrationStep,
            /* [defaultvalue][in][optional] */ double lagStartDate,
            /* [defaultvalue][in][optional] */ double creditLag,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVectorPrevIndexId,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pMatrixPrevBC,
            /* [defaultvalue][in][optional] */ double step,
            /* [defaultvalue][in][optional] */ BSTR CalMeth,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDO2 )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pEffectiveDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pPortfolio,
            /* [in] */ double pSpread,
            /* [in] */ double pSubAmount,
            /* [in] */ double pMezzAmount,
            /* [defaultvalue][in][optional] */ BSTR pFreqFeeLeg,
            /* [defaultvalue][in][optional] */ BSTR pFreqDefLeg,
            /* [defaultvalue][in][optional] */ BSTR pDayCountFrq,
            /* [defaultvalue][in][optional] */ double pFirst_period_refdate,
            /* [defaultvalue][in][optional] */ BSTR pAccruedOnDefault,
            /* [defaultvalue][in][optional] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ double pPayCreditLag,
            /* [defaultvalue][in][optional] */ BSTR pStub,
            /* [defaultvalue][in][optional] */ double pBinary,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ BSTR CrossSub,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Portfolio )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pSecuritiesID,
            /* [defaultvalue][in][optional] */ BSTR Parameters,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCreateZCSwapInt )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [in] */ BSTR pMMVsFut,
            /* [in] */ BSTR pSwapVsFut,
            /* [in] */ BSTR pRaw,
            /* [in] */ BSTR pInterp,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetInitialCurveFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in] */ BSTR pAdjOrNot,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetRate);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMTHREEMONTHFUT )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDelivery,
            /* [in] */ long pMarket,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFutPibor )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDelivery,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMIRFUT )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pIdUnderlying,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLibor )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pLiborTypeId,
            /* [in] */ BSTR pCcyId,
            /* [in] */ BSTR pResetFreqId,
            /* [in] */ BSTR pPayFreqId,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [defaultvalue][in] */ BSTR pBasis,
            /* [defaultvalue][in] */ BSTR pIntrule,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLiborSwaption )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pReceiveOrPay,
            /* [in] */ double pStrike,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pLiborType,
            /* [in] */ double pSpread,
            /* [in] */ BSTR pExerciseType,
            /* [in] */ BSTR pResetFreq,
            /* [in] */ BSTR pPayFreq,
            /* [in] */ BSTR pCcyId,
            /* [defaultvalue][in] */ BSTR pRefObj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFixedLeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pReceiveOrPay,
            /* [in] */ double pFixRate,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [defaultvalue][in][optional] */ BSTR pFreq,
            /* [defaultvalue][in][optional] */ BSTR pDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ BSTR pPayCalName,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetCorrelationMatrix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pModelMultiCurvesId,
            /* [in] */ BSTR pCorrMatrixId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DPMktDataFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR Issuer,
            /* [in] */ BSTR CurveName,
            /* [in] */ BSTR Parameter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetDPFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR Issuer,
            /* [in] */ BSTR CurveName,
            /* [in] */ BSTR ircurveId,
            /* [in] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CloneCorrMatrixBary )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrMatrixId,
            /* [in] */ double Beta,
            /* [in] */ int UpOrDown,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ConstantDefaultCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [in] */ VARIANT __RPC_FAR *pRates,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR Ccy,
            /* [in] */ BSTR Label,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ BSTR IsSummit,
            /* [defaultvalue][optional][in] */ BSTR calibrationData,
            /* [defaultvalue][optional][in] */ int lag,
            /* [optional][in] */ BSTR calibrationAlgo,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DefProbModelNew )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDefCurve,
            /* [in] */ BSTR pIRcurve,
            /* [defaultvalue][in][optional] */ BSTR VolCurve,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ZeroCouponDefaultCurveFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR bIssuer,
            /* [in] */ BSTR bCurrency,
            /* [in] */ BSTR bCvName,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR bLabel,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBsSlModel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pVolSpreadLock,
            /* [in] */ BSTR pCvCapVol,
            /* [in] */ BSTR pCvIndexVol,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFXVolFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [defaultvalue][in][optional] */ BSTR pType,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFXCorrelFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [in] */ VARIANT __RPC_FAR *pTenors,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMInfocentreConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGlobDFBS )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDomBSId,
            /* [in] */ BSTR pDomCurrId,
            /* [in] */ BSTR pFrgBSId,
            /* [in] */ BSTR pFrgCurrId,
            /* [in] */ BSTR pFxVolCrvId,
            /* [in] */ BSTR pFFxCorrId,
            /* [in] */ BSTR pRatesCorrId,
            /* [defaultvalue][in][optional] */ BSTR pFxVolModelId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetInitialCurveFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [in] */ VARIANT __RPC_FAR *pCvName,
            /* [in] */ VARIANT __RPC_FAR *pDate,
            /* [in] */ VARIANT __RPC_FAR *value,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetCorrelFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ BSTR pIndex2,
            /* [in] */ double pDate,
            /* [in] */ BSTR pCvName,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMDFFXBS )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDVolId,
            /* [in] */ BSTR pFVolId,
            /* [in] */ BSTR pDZcId,
            /* [in] */ BSTR pFZcId,
            /* [in] */ BSTR pDFxCorrId,
            /* [in] */ BSTR pFFxCorrId,
            /* [in] */ BSTR pFxVolId,
            /* [in] */ double pRatesCorr,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMTRIBSMODEL )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pModel1,
            /* [in] */ BSTR pModel2,
            /* [in] */ BSTR pDiscModel,
            /* [in] */ BSTR pFX1DiscVol,
            /* [in] */ BSTR pFX2DiscVol,
            /* [in] */ BSTR pIdx1Idx2Corr,
            /* [in] */ BSTR pIdx1DiscIdxCorr,
            /* [in] */ BSTR pIdx2DiscIdxCorr,
            /* [in] */ BSTR pIdx1FxCorr,
            /* [in] */ BSTR pIdx2FxCorr,
            /* [in] */ int pQuantoFlag,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMTRIBSDUAL )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pModel1,
            /* [in] */ BSTR pModel2,
            /* [in] */ BSTR pDiscModel,
            /* [in] */ BSTR pFX1DiscVol,
            /* [in] */ BSTR pFX2DiscVol,
            /* [in] */ BSTR pIdx1Idx2Corr,
            /* [in] */ BSTR pIdx1DiscIdxCorr,
            /* [in] */ BSTR pIdx2DiscIdxCorr,
            /* [in] */ BSTR pIdx1FxCorr,
            /* [in] */ BSTR pIdx2FxCorr,
            /* [in] */ int pQuantoFlag,
            /* [in] */ double pCorrelForAdj,
            /* [in] */ int pWithslopeflag,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcptBonibor )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pDepart,
            /* [in] */ BSTR pMatStruct,
            /* [in] */ BSTR pMatTot,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pSjUSD,
            /* [in] */ BSTR pTiming,
            /* [in] */ double pBarriere,
            /* [in] */ double pSpdPostBar,
            /* [in] */ double pMarge,
            /* [in] */ double pFunding,
            /* [in] */ BSTR pFundingFreq,
            /* [in] */ double pSpd2phase,
            /* [in] */ double pSoulte,
            /* [in] */ BSTR pYcModId,
            /* [in] */ BSTR pBsModId,
            /* [in] */ BSTR pBsModVolUSDId,
            /* [in] */ BSTR pBsModCorrPlusId,
            /* [in] */ BSTR pBsModCorrMoinsId,
            /* [in] */ BSTR pCrossModId,
            /* [in] */ BSTR pProbaMarge,
            /* [defaultvalue][in][optional] */ double pInt,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CMTranche )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pParticipationRate,
            /* [in] */ VARIANT __RPC_FAR *pMezzAmount,
            /* [in] */ VARIANT __RPC_FAR *pSubAmount,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pFreqFeeLeg,
            /* [in] */ VARIANT __RPC_FAR *pDayCount,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
            /* [in] */ VARIANT __RPC_FAR *Currency,
            /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
            /* [in] */ VARIANT __RPC_FAR *pStub,
            /* [in] */ VARIANT __RPC_FAR *pFreqDefLeg,
            /* [in] */ VARIANT __RPC_FAR *pBinary,
            /* [in] */ VARIANT __RPC_FAR *pPayCal,
            /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
            /* [defaultvalue][in][optional] */ double TradedNotional,
            /* [defaultvalue][in][optional] */ double FwdFixedDate,
            /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
            /* [defaultvalue][in] */ double pFstCpnEffDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Index )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ double YearFrac,
            /* [defaultvalue][in] */ double pSpread,
            /* [defaultvalue][in] */ BSTR Method,
            /* [defaultvalue][in] */ BSTR Basis,
            /* [defaultvalue][in] */ BSTR ResetFreq,
            /* [defaultvalue][in] */ BSTR PayFreq,
            /* [defaultvalue][in] */ BSTR ccy,
            /* [defaultvalue][in] */ BSTR DefCurve,
            /* [defaultvalue][in] */ BSTR fwdRule,
            /* [defaultvalue][in] */ BSTR resetTiming,
            /* [defaultvalue][in] */ int resetGap,
            /* [defaultvalue][in] */ BSTR payTiming,
            /* [defaultvalue][in] */ int payGap,
            /* [defaultvalue][in] */ BSTR intRule,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ int cm_resetWeekDay,
            /* [defaultvalue][in] */ int cm_resetOccur,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Parameters )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [in] */ long nbcols,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMAswPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pMaturity,
            /* [in] */ double pCpn,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pBase,
            /* [in] */ double pMargin,
            /* [defaultvalue][in][optional] */ double pRedemptionPrice,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pFixDecompfreq,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ double pMinValue,
            /* [defaultvalue][in][optional] */ double pMaxValue,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMAswMargin )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pMaturity,
            /* [in] */ double pCpn,
            /* [in] */ BSTR pFreq,
            /* [in] */ BSTR pBase,
            /* [in] */ double pPrice,
            /* [defaultvalue][in][optional] */ double pRedemptionPrice,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ BSTR pFixDecompfreq,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ double pMinValue,
            /* [defaultvalue][in][optional] */ double pMaxValue,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDSIndex )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
            /* [in] */ VARIANT __RPC_FAR *pEndDate,
            /* [in] */ VARIANT __RPC_FAR *pSpread,
            /* [in] */ VARIANT __RPC_FAR *pIndex,
            /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
            /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
            /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
            /* [in] */ VARIANT __RPC_FAR *pFixedPayerAmount,
            /* [in] */ VARIANT __RPC_FAR *pFloatingPayerAmount,
            /* [in] */ BSTR StubRule,
            /* [in] */ VARIANT __RPC_FAR *pCurrency,
            /* [defaultvalue][in] */ BSTR Adjusted,
            /* [defaultvalue][in] */ int CreditLag,
            /* [defaultvalue][in] */ BSTR IncludeMaturity,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionStartDate,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionEndDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Option )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *UnderlyingMaturity,
            /* [in] */ double OptionExpiry,
            /* [in] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ BSTR CdsAdj,
            /* [defaultvalue][in][optional] */ BSTR EndAdj,
            /* [defaultvalue][in][optional] */ double pStrike,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pOptionType,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pKoType,
            /* [defaultvalue][in][optional] */ double Notional,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FwdSpreadPricer )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ double Maturity1,
            /* [in] */ double Maturity2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ImpliedVol )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ double pMktPrice,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_VirtualCdsSpread )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pMaturity,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_BSGreeks )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ VARIANT __RPC_FAR *pGreekType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetCorrelation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pModelMultiCurvesId,
            /* [in] */ BSTR pCorrelationId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CorrelationStrike )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pVolCurves,
            /* [in] */ VARIANT __RPC_FAR *pProportions,
            /* [in] */ VARIANT __RPC_FAR *pSmileStrikeLow,
            /* [in] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
            /* [in] */ VARIANT __RPC_FAR *pIndexVector,
            /* [defaultvalue][in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Beta_Correlation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pCoefs,
            /* [in] */ double AsOf,
            /* [defaultvalue][in] */ BSTR Name,
            /* [defaultvalue][in] */ BSTR idIndex1,
            /* [defaultvalue][in] */ BSTR idIndex2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_GETINSTRUMENTFROMSUMMIT )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR SummitId,
            /* [in] */ BSTR Type,
            /* [defaultvalue][in][optional] */ double AsOf,
            /* [defaultvalue][in][optional] */ BSTR ExoticFilter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFrnPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [in] */ double pFacialMargin,
            /* [in] */ double pValoMargin,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ double pFixing,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ double pOutMode,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFrnMargin )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pDelivery,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pIndex1,
            /* [in] */ BSTR pFwdCurve1,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
            /* [in] */ double pFacialMargin,
            /* [in] */ double pPrice,
            /* [defaultvalue][in][optional] */ BSTR pCcy2,
            /* [defaultvalue][in][optional] */ BSTR pIndex2,
            /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
            /* [defaultvalue][in][optional] */ double pFixing,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ double pOutMode,
            /* [defaultvalue][in][optional] */ long pSolve,
            /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetEqStrikeDown )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetEqStrikeUp )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetCorrelStrikeDown )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrelId,
            /* [in] */ double yfmaturity,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetCorrelStrikeUp )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrelId,
            /* [in] */ double yfmaturity,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetCorrelation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR ModelId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetModelFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR IRcurve,
            /* [in] */ BSTR IDSummit,
            /* [in] */ BSTR type,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetVolatility )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR VolCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_NextCpnDate )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ double maturity,
            /* [in] */ BSTR frequency,
            /* [defaultvalue][in][optional] */ BSTR rule,
            /* [defaultvalue][in][optional] */ BSTR currency,
            /* [defaultvalue][in][optional] */ BSTR intrule,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetProportionsInfos )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR correlId,
            /* [in] */ BSTR IndexName,
            /* [in] */ double proportion,
            /* [defaultvalue][in][optional] */ double forcedstrikelow,
            /* [defaultvalue][in][optional] */ double forcedstrikehigh);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcptDigital )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodDelta,
            /* [in] */ BSTR pBsmodVega,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pPorR,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCcyIdx,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pCapOrFloor,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pStrike,
            /* [in] */ BSTR pPayOff,
            /* [in] */ BSTR pSpd,
            /* [in] */ double pResetGap,
            /* [defaultvalue][in][optional] */ double pSpreadBelow,
            /* [defaultvalue][in][optional] */ double pSpreadAbove,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMRefValue )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pdates,
            /* [in] */ VARIANT __RPC_FAR *pvalues,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pvalues2,
            /* [defaultvalue][in][optional] */ long valueType,
            /* [defaultvalue][in][optional] */ long conversion,
            /* [defaultvalue][in][optional] */ BSTR calcMethod,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCreateGenCorrelManager )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pMktTags,
            /* [in] */ VARIANT __RPC_FAR *pIntraMktTags,
            /* [in] */ VARIANT __RPC_FAR *pCorrelCurveIds,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBSConvAdjust )( 
            IARMModule __RPC_FAR * This,
            /* [defaultvalue][in][optional] */ BSTR pSUMMITFormulaeUsed,
            /* [defaultvalue][in][optional] */ BSTR pUseSABRCMS,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBsModelGen )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pYieldCurve,
            /* [in] */ BSTR pVolatility,
            /* [defaultvalue][in][optional] */ BSTR pCorrMgr,
            /* [defaultvalue][in][optional] */ BSTR pCnvxManager,
            /* [defaultvalue][in][optional] */ BSTR pCapletVol,
            /* [defaultvalue][in][optional] */ BSTR pSpreadLock,
            /* [defaultvalue][in][optional] */ BSTR pDiscCurve,
            /* [defaultvalue][in][optional] */ BSTR pCorrel,
            /* [defaultvalue][in][optional] */ BSTR pCashVol,
            /* [defaultvalue][in][optional] */ BSTR pSpreadVol,
            /* [defaultvalue][in][optional] */ BSTR pModelType,
            /* [defaultvalue][in][optional] */ BSTR pSpreadVolType,
            /* [defaultvalue][in][optional] */ BSTR pSabrMod,
            /* [defaultvalue][in][optional] */ BSTR pLnorNorVol,
            /* [defaultvalue][in][optional] */ long pNumSteps,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcptSPTQTF )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ BSTR pIndexPhase2,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ BSTR pFreqPPhase1,
            /* [in] */ BSTR pFreqPPhase2,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqPPhase3,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pResetTimingPhase2,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pCcyIdx,
            /* [in] */ BSTR pDayCount,
            /* [in] */ double pFee,
            /* [in] */ BSTR pIsRateFixedPhase2,
            /* [in] */ double pFixedRatePhase2,
            /* [in] */ BSTR pBarrier,
            /* [in] */ BSTR pSpdPhase1,
            /* [in] */ BSTR pSpdPhase1Fund,
            /* [in] */ BSTR pSpdPhase2Tf,
            /* [in] */ BSTR pSpdPhase2fund,
            /* [in] */ BSTR pSpdPhase3,
            /* [in] */ BSTR pSpdPhase3fund,
            /* [in] */ double pResetGapPhase1,
            /* [in] */ double pResetGapPhase2,
            /* [in] */ double pResetGapPhase3,
            /* [in] */ BSTR pAmort,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodDeltaCcy1,
            /* [in] */ BSTR pBsmodVegaCcy1,
            /* [defaultvalue][in][optional] */ BSTR pBsmodDeltaCcy2,
            /* [defaultvalue][in][optional] */ BSTR pBsmodVegaCcy2,
            /* [defaultvalue][in][optional] */ BSTR pBsmodFxCorrel,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMDisplaySchedule )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR legId,
            /* [in] */ BSTR dataType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMIrIndex )( 
            IARMModule __RPC_FAR * This,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ double pMaturity,
            /* [defaultvalue][in][optional] */ BSTR pCompMethod,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ double pPayGap,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIndexType,
            /* [defaultvalue][in][optional] */ double pDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwapleg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndexId,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pRecOrPay,
            /* [in] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMConstRefvalue )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pValue,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CptImplCvForCDO2 )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricerId,
            /* [in] */ BSTR Name,
            /* [in] */ BSTR Tenor,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_AddPeriod )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR Maturity,
            /* [in] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR AdjRule,
            /* [defaultvalue][in][optional] */ BSTR AdjCDS,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMShutdownETK )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetCoupons )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CdsorCdoId,
            /* [in] */ BSTR CouponsId,
            /* [in] */ BSTR TypesId,
            /* [defaultvalue][in][optional] */ BSTR PartId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_InputDefaultCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pDates,
            /* [in] */ VARIANT __RPC_FAR *pRates,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [in] */ BSTR bCurrency,
            /* [in] */ BSTR bLabel,
            /* [defaultvalue][in][optional] */ BSTR bInterpolType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetInitialVolFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [in] */ BSTR pType,
            /* [defaultvalue][in][optional] */ BSTR pMatIndex,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBond )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedempPrice,
            /* [in] */ double pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ double pSettleGap,
            /* [defaultvalue][in][optional] */ double pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFCreateOATLeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pInfIdx,
            /* [in] */ BSTR pRcvOrPay,
            /* [defaultvalue][in][optional] */ BSTR pInterpType,
            /* [defaultvalue][in][optional] */ double pLeverage,
            /* [defaultvalue][in][optional] */ double pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pResetNumGap,
            /* [defaultvalue][in][optional] */ double pResetDenomGap,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ double pPayGap,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ BSTR pFinalNotionalType,
            /* [defaultvalue][in][optional] */ double pFirstReset,
            /* [defaultvalue][in][optional] */ double pCoMultiple,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwap )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSwapleg1,
            /* [in] */ BSTR pSwapleg2,
            /* [defaultvalue][in][optional] */ double pMinPay,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMPToYield )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pPrice,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMYToPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pYield,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMYToDuration )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pBond,
            /* [in] */ double pSettleDate,
            /* [in] */ double pActuRate,
            /* [defaultvalue][in][optional] */ double pFlagCpn,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLiborleg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pLiborType,
            /* [in] */ BSTR pRecOrPay,
            /* [defaultvalue][in][optional] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFReq,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [defaultvalue][in][optional] */ BSTR pCpnDaycount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMImpliedSpread )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSwap,
            /* [in] */ BSTR pModel,
            /* [in] */ double pPrice,
            /* [in] */ double pLeg1or2,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMDiscountPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZeroCurve,
            /* [in] */ double pMatu,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFCreateCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pIndexName,
            /* [in] */ double pCPIIndexValue,
            /* [in] */ double pCPIIndexDate,
            /* [in] */ VARIANT __RPC_FAR *pMatu,
            /* [in] */ VARIANT __RPC_FAR *pRate,
            /* [defaultvalue][in][optional] */ BSTR pMonthlyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pDCFMonthly,
            /* [defaultvalue][in][optional] */ BSTR pDCFDaily,
            /* [defaultvalue][in][optional] */ BSTR pExtrapolType,
            /* [defaultvalue][in][optional] */ BSTR pResetManager,
            /* [defaultvalue][in][optional] */ BSTR pSeasonManager,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFInterpCPI )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZc,
            /* [in] */ double pCPIDate,
            /* [defaultvalue][in][optional] */ BSTR pDCFlag,
            /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
            /* [defaultvalue][in][optional] */ BSTR pCPIlag,
            /* [defaultvalue][in][optional] */ double pWeight,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFSeasonManager )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pMonthList,
            /* [in] */ VARIANT __RPC_FAR *pValues,
            /* [defaultvalue][in][optional] */ BSTR pSeasonAdjMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFResetManager )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pDatas,
            /* [in] */ double pNbIndex,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMINFYcMod )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pYieldCurve,
            /* [in] */ BSTR pInfCurve,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBaseReplicationConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetInitialFXVolFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pCvName,
            /* [defaultvalue][in][optional] */ BSTR pImpOrHist,
            /* [defaultvalue][in][optional] */ BSTR pVolType,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCreateZCFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pCvName,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pAdj,
            /* [defaultvalue][in][optional] */ BSTR pRaw,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBumpCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pZc,
            /* [in] */ double pEpsilon,
            /* [defaultvalue][in][optional] */ long pMethod,
            /* [defaultvalue][in][optional] */ BSTR pPlot,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMAccrued )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSec,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pModel,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSwitchToWSETK )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DataFromLabel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricer,
            /* [in] */ BSTR label,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeBilibor )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pDateSecondPhase,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexSF,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodFund,
            /* [in] */ BSTR pBsmodDeltaCcy1,
            /* [in] */ BSTR pBsmodDeltaFund,
            /* [in] */ BSTR pBsmodDeltaCcy2,
            /* [in] */ BSTR pBsmodFxCorrel,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqRFund,
            /* [in] */ BSTR pFreqPSF,
            /* [in] */ BSTR pFreqRSF,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pResetTimingSF,
            /* [in] */ BSTR pCcy1,
            /* [in] */ BSTR pCcy2,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pDayCountSF,
            /* [in] */ double pSpdPF,
            /* [in] */ double pSpdSF,
            /* [in] */ double pSpdfund,
            /* [in] */ double pSpdfund2,
            /* [in] */ double pResetGap,
            /* [in] */ double pResetGapSF,
            /* [in] */ BSTR pAmort,
            /* [in] */ double pRefDate,
            /* [in] */ double pFee,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenerateImpliedCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetEqStrike )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR IndexName,
            /* [in] */ BSTR UpOrLow,
            /* [out] */ VARIANT __RPC_FAR *pRetMatu,
            /* [out] */ VARIANT __RPC_FAR *pRetStrikes);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeOptilix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pDateSecondPhase,
            /* [in] */ double pEndDate,
            /* [in] */ double pNtl,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ BSTR pIndexSF,
            /* [in] */ BSTR pBsmod,
            /* [in] */ BSTR pBsmodFund,
            /* [in] */ BSTR pBsmodDeltaCcy,
            /* [in] */ BSTR pBsmodDeltaFund,
            /* [in] */ BSTR pFreqP,
            /* [in] */ BSTR pFreqR,
            /* [in] */ BSTR pFreqPFund,
            /* [in] */ BSTR pFreqRFund,
            /* [in] */ BSTR pFreqPSF,
            /* [in] */ BSTR pFreqRSF,
            /* [in] */ BSTR pResetTiming,
            /* [in] */ BSTR pResetTimingSF,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pDayCountSF,
            /* [in] */ double pSpdSF,
            /* [in] */ VARIANT pSpdfund,
            /* [in] */ double pResetGap,
            /* [in] */ double pResetGapSF,
            /* [in] */ BSTR pAmort,
            /* [in] */ double pRefDate,
            /* [in] */ double pFee,
            /* [defaultvalue][in][optional] */ BSTR pFwdRule,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
            /* [defaultvalue][in][optional] */ double pTxAmort,
            /* [defaultvalue][in][optional] */ double pAmountAmort,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DefaultIntensity )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricerId,
            /* [in] */ double Maturity,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMRiskyBond )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedemptionPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double pRepo,
            /* [defaultvalue][in][optional] */ double pSsl,
            /* [defaultvalue][in][optional] */ double pRecoveryRate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMRiskyBondWithCF )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOfDate,
            /* [in] */ double pRedemptionPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [in] */ VARIANT __RPC_FAR *pYearTerms,
            /* [in] */ VARIANT __RPC_FAR *pCashFlows,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double pRepo,
            /* [defaultvalue][in][optional] */ double pSsl,
            /* [defaultvalue][in][optional] */ double pRecoveryRate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_EmptyLeg )( 
            IARMModule __RPC_FAR * This,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMClonedAndSetNotional )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR bLegId,
            /* [defaultvalue][in][optional] */ BSTR bAmortId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_INF_GetZcFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR Index,
            /* [in] */ BSTR Ccy,
            /* [in] */ BSTR cvname,
            /* [in] */ double date,
            /* [defaultvalue][in][optional] */ BSTR seasonAdj,
            /* [defaultvalue][in][optional] */ BSTR seasonAdjMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_IRLEGTOCREDITLEG )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR SwapLegId,
            /* [in] */ BSTR LegType,
            /* [in] */ BSTR creditindexId,
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Collateral )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDSGEN )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR FeeLegId,
            /* [in] */ BSTR DefLegId,
            /* [in] */ double RcvFee,
            /* [in] */ double TradedNot,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_NTDGEN )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CdsId,
            /* [in] */ int firstnumdef,
            /* [in] */ int lastnumdef,
            /* [in] */ BSTR CollateralId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDOGEN )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CdsId,
            /* [in] */ double subamount,
            /* [in] */ BSTR CollateralId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenLeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double StartDate,
            /* [in] */ double EndDate,
            /* [in] */ double FixedRate,
            /* [in] */ double FixedNotional,
            /* [defaultvalue][in][optional] */ BSTR RefValNotional,
            /* [defaultvalue][in][optional] */ BSTR RefValRate,
            /* [defaultvalue][in][optional] */ BSTR XChangeNotional,
            /* [defaultvalue][in][optional] */ BSTR Frequency,
            /* [defaultvalue][in][optional] */ BSTR Basis,
            /* [defaultvalue][in][optional] */ BSTR payTiming,
            /* [defaultvalue][in][optional] */ BSTR intrule,
            /* [defaultvalue][in][optional] */ BSTR stubrule,
            /* [defaultvalue][in][optional] */ BSTR ccyid,
            /* [defaultvalue][in][optional] */ BSTR paycalname,
            /* [defaultvalue][in][optional] */ double refdate,
            /* [defaultvalue][in][optional] */ BSTR includematurity,
            /* [defaultvalue][in][optional] */ BSTR adjstartdate,
            /* [defaultvalue][in][optional] */ BSTR legtype,
            /* [defaultvalue][in][optional] */ BSTR indexobj,
            /* [defaultvalue][in][optional] */ int creditlag,
            /* [defaultvalue][in][optional] */ double binary,
            /* [defaultvalue][in][optional] */ BSTR name,
            /* [defaultvalue][in][optional] */ BSTR Nxchange,
            /* [defaultvalue][in][optional] */ BSTR baccruedOnDef,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CDO_SQUARE_GEN )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CdsId,
            /* [in] */ double subamount,
            /* [in] */ BSTR portfolioId,
            /* [in] */ double binary,
            /* [in] */ double rcvfee,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetLastFixingDate )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR instId,
            /* [in] */ VARIANT asofDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetPastFixing )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR instId,
            /* [in] */ VARIANT resetDate,
            /* [in] */ VARIANT fixingValue,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFixing )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR source,
            /* [in] */ BSTR index,
            /* [in] */ BSTR term,
            /* [in] */ BSTR ccy,
            /* [in] */ double date,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetRiskyProfile )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CdsorCdoId,
            /* [in] */ BSTR CouponsId,
            /* [in] */ BSTR TypesId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMHyperCube )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
            /* [in] */ VARIANT __RPC_FAR *pKeys,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetPricerForRatesComputation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR legId,
            /* [in] */ BSTR pricerId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCmsLeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double startDate,
            /* [in] */ double endDate,
            /* [in] */ BSTR cmsTypeId,
            /* [defaultvalue][in][optional] */ BSTR receiveOrPay,
            /* [defaultvalue][in][optional] */ BSTR yieldDecompFreq,
            /* [defaultvalue][in][optional] */ BSTR swapLegDayCount,
            /* [defaultvalue][in][optional] */ BSTR resetFreq,
            /* [defaultvalue][in][optional] */ BSTR payFreq,
            /* [defaultvalue][in][optional] */ long resetGap,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR ccyName,
            /* [defaultvalue][in][optional] */ BSTR resetTiming,
            /* [defaultvalue][in][optional] */ BSTR stubRule,
            /* [defaultvalue][in][optional] */ BSTR adjStartDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetMatuLabel )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurveId,
            /* [in] */ VARIANT __RPC_FAR *pMatuLabels,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputePentifix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ BSTR pFloorOrCap,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pCorrManager,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
            /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_ReplicConvAdjust_Create )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR Payoff_ReplicMode,
            /* [in] */ double Payoff_StepOrReplicPrecision,
            /* [in] */ BSTR Payoff_StopMode,
            /* [in] */ double Payoff_StopThreshold,
            /* [in] */ BSTR Sensi_ReplicMode,
            /* [in] */ double Sensi_StepOrReplicPrecision,
            /* [in] */ BSTR Sensi_StopMode,
            /* [in] */ double Sensi_StopThreshold,
            /* [in] */ BSTR UsedModelId,
            /* [in] */ double StrikeMinReplic,
            /* [in] */ double StrikeMaxReplic,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_MapConvAdjust_Create )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR LiborArrearAdj,
            /* [in] */ BSTR NaturalCMSAdj,
            /* [in] */ BSTR PaymentLagAdj,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetFees )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR securityId,
            /* [in] */ BSTR RefvalueId);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetBounds )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR securityId,
            /* [out] */ double __RPC_FAR *down,
            /* [out] */ double __RPC_FAR *up);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputePentibor )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPay,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ BSTR pIndexIndexCorrelCube,
            /* [in] */ BSTR pCorrEUR,
            /* [in] */ BSTR pInterCorr,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetRecovCoef )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurveId,
            /* [in] */ double RecovCoef,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMIndexIndexCorrelCube )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
            /* [in] */ VARIANT __RPC_FAR *pTenors1List,
            /* [in] */ VARIANT __RPC_FAR *pTenors2List,
            /* [defaultvalue][in][optional] */ BSTR pInterSurfInterp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetInterpolationType )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pVolCurveId,
            /* [in] */ VARIANT __RPC_FAR *pInterpolType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Customized_CDO )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pNotionals,
            /* [defaultvalue][in] */ BSTR Currency,
            /* [in] */ BSTR pDefaultLeg,
            /* [in] */ BSTR pPremiumLeg,
            /* [in] */ BSTR pParameters,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCreateGenCorrelatorManager )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pMktTags,
            /* [in] */ VARIANT __RPC_FAR *pHyperDiagVol,
            /* [in] */ VARIANT __RPC_FAR *pIndexIndexVol,
            /* [in] */ VARIANT __RPC_FAR *pCorrelVol,
            /* [in] */ VARIANT __RPC_FAR *pIndexVol,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CLN )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double EffectiveDate,
            /* [in] */ double EndDate,
            /* [in] */ double Spread,
            /* [defaultvalue][in][optional] */ BSTR IndexId,
            /* [defaultvalue][in][optional] */ double Refdate,
            /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
            /* [defaultvalue][in][optional] */ double Notional,
            /* [defaultvalue][in][optional] */ BSTR AccOnDef,
            /* [defaultvalue][in][optional] */ BSTR DayCount,
            /* [defaultvalue][in][optional] */ BSTR DecompFreq,
            /* [defaultvalue][in][optional] */ BSTR StubRule,
            /* [defaultvalue][in][optional] */ double resetgap,
            /* [defaultvalue][in][optional] */ BSTR Currency,
            /* [defaultvalue][in][optional] */ BSTR ResetCal,
            /* [defaultvalue][in][optional] */ BSTR PayCal,
            /* [defaultvalue][in][optional] */ BSTR Nxchange,
            /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
            /* [defaultvalue][in][optional] */ BSTR AdjustedStartDate,
            /* [defaultvalue][in][optional] */ double Binary,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputePentilix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ double pStartdatePhase1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPay,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ double pSpreadPhase1,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ double pStartDatePhase2,
            /* [in] */ double pSpreadPhase2,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ double pStartDatePhase3,
            /* [in] */ double pEndDatePhase3,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ double pSpreadPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ BSTR pIndexIndexCorrelCube,
            /* [in] */ BSTR pCorrEUR,
            /* [in] */ BSTR pInterCorr,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_RiskyPV01 )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDefCurve,
            /* [optional][in] */ VARIANT Date1,
            /* [in] */ VARIANT Date2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetResetMgrFromSummit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pSource,
            /* [defaultvalue][in] */ BSTR pCcy,
            /* [defaultvalue][in] */ BSTR pIsInflationIndex,
            /* [defaultvalue][in] */ BSTR pTerm,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetReset )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pResetMgr,
            /* [in] */ double pDate,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSetLastFixing )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSecurityId,
            /* [in] */ double pRate,
            /* [in] */ double pAsOf,
            /* [defaultvalue][in][optional] */ double pBeforeLastFixingDate,
            /* [defaultvalue][in][optional] */ double pResetDate,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Sectorial_Correlation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ DATE AsOf,
            /* [in] */ BSTR structName,
            /* [in] */ BSTR correlation_Type,
            /* [in] */ VARIANT __RPC_FAR *vLabels,
            /* [in] */ VARIANT __RPC_FAR *vector_Membership,
            /* [defaultvalue][in][optional] */ double intra_Sector_Correlation,
            /* [defaultvalue][in][optional] */ double inter_Sector_Correlation,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas_Down,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas_Down,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeReviPentix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ VARIANT pDate,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexPhase1,
            /* [in] */ VARIANT pSpread,
            /* [in] */ BSTR pDayCountPhase1,
            /* [in] */ BSTR pPayFreqPhase1,
            /* [in] */ BSTR pResetFreqPhase1,
            /* [in] */ BSTR pResetTimingPhase1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pAdjPhase1,
            /* [in] */ BSTR pStub,
            /* [in] */ BSTR pIndexPhase2DIG,
            /* [in] */ BSTR pIndexLongPhase2DIG,
            /* [in] */ BSTR pStrikePhase2DIG,
            /* [in] */ BSTR pResetTimingPhase2DIG,
            /* [in] */ BSTR pAdjPhase2DIG,
            /* [in] */ BSTR pDayCountPhase2,
            /* [in] */ BSTR pPayFreqPhase2,
            /* [in] */ BSTR pResetFreqPhase2,
            /* [in] */ BSTR pAdjPhase2,
            /* [in] */ BSTR pIndexPhase3,
            /* [in] */ BSTR pDayCountPhase3,
            /* [in] */ BSTR pPayFreqPhase3,
            /* [in] */ BSTR pResetFreqPhase3,
            /* [in] */ BSTR pResetTimingPhase3,
            /* [in] */ BSTR pAdjPhase3,
            /* [in] */ BSTR pIndexFund,
            /* [in] */ VARIANT pSpreadFund,
            /* [in] */ BSTR pDayCountFund,
            /* [in] */ BSTR pPayFreqFund,
            /* [in] */ BSTR pResetFreqFund,
            /* [in] */ BSTR pResetTimingFund,
            /* [in] */ BSTR pAdjFund,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            /* [in] */ BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ BSTR pFloorOrCap,
            /* [in] */ double pFee,
            /* [in] */ double pLevier,
            /* [in] */ double pTxFixeMax,
            /* [in] */ BSTR pIsCapped,
            /* [in] */ double pTxCap,
            /* [in] */ BSTR pVolCurvFromMatriceShift,
            /* [in] */ BSTR pVol,
            /* [in] */ BSTR pVolCub,
            /* [in] */ BSTR pCorrManager,
            /* [in] */ BSTR pConvexityManager,
            /* [in] */ BSTR pZc,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pHyperCubeCorrel,
            /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
            /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGenAmortization )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pSwaplegId,
            /* [in] */ BSTR pAmortMethod,
            /* [defaultvalue][in][optional] */ BSTR pAmortFreq,
            /* [defaultvalue][in][optional] */ double pAmortAmount,
            /* [defaultvalue][in][optional] */ BSTR pDaycount,
            /* [defaultvalue][in][optional] */ double pLegNotional,
            /* [defaultvalue][in][optional] */ double pAmortRate,
            /* [defaultvalue][in][optional] */ double pReducedMaturity,
            /* [defaultvalue][in][optional] */ BSTR pModelId,
            /* [defaultvalue][in][optional] */ double pCleanUp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCptRefvalue )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pRefValId,
            /* [in] */ double pDate,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetDPFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR irCurveId,
            /* [defaultvalue][in][optional] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *ret);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCalypsoDevConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetZCFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pTerm,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DefProbInverse )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurveId,
            /* [in] */ double dDefProba,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DPMktDataFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [in] */ VARIANT __RPC_FAR *Parameter,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ZeroCouponDefaultCurveFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pDate,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ BSTR issuer,
            /* [in] */ BSTR seniority,
            /* [in] */ BSTR ccy,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR irCurveId,
            /* [defaultvalue][in][optional] */ BSTR label,
            /* [retval][out] */ VARIANT __RPC_FAR *ret);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetInitialCurveFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pIndex,
            /* [in] */ BSTR pCurrency,
            /* [in] */ BSTR pTerm,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFile,
            /* [defaultvalue][in][optional] */ BSTR pDoAdj,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetRate);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCalypsoProdConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMCalypsoRecConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetBasketCorrelMkDataFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double date,
            /* [in] */ BSTR forceCurveName,
            /* [defaultvalue][in][optional] */ BSTR xmlFileName,
            /* [out] */ VARIANT __RPC_FAR *pRetMat,
            /* [out] */ VARIANT __RPC_FAR *pRetTenor,
            /* [out] */ VARIANT __RPC_FAR *pRetVol);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_QMatrix )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pQMatrix,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_MarketDataMng )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pstrVect,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ModelMultiCvMktDataMng )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pIRcurve,
            /* [in] */ VARIANT __RPC_FAR *pDefCurves,
            /* [in] */ VARIANT __RPC_FAR *pRecovery,
            /* [in] */ BSTR CorrelId,
            /* [in] */ BSTR MktdataMngId,
            /* [defaultvalue][in][optional] */ BSTR pVolcurve,
            /* [defaultvalue][in][optional] */ BSTR cloneorNot,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Local_ARM_ProdConnect )( 
            IARMModule __RPC_FAR * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_SetDefaultCurrency )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR isoCCy,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLivretALeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pRcvOrPay,
            /* [defaultvalue][in][optional] */ VARIANT pSpread,
            /* [defaultvalue][in][optional] */ BSTR pResetFreq,
            /* [defaultvalue][in][optional] */ BSTR pPayFreq,
            /* [defaultvalue][in][optional] */ BSTR pResetTiming,
            /* [defaultvalue][in][optional] */ BSTR pPayTiming,
            /* [defaultvalue][in][optional] */ BSTR pCcy,
            /* [defaultvalue][in][optional] */ BSTR pIntRule,
            /* [defaultvalue][in][optional] */ double pResetGap,
            /* [defaultvalue][in][optional] */ BSTR pResetCal,
            /* [defaultvalue][in][optional] */ BSTR pPayCal,
            /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
            /* [defaultvalue][in][optional] */ BSTR pNxChange,
            /* [defaultvalue][in][optional] */ BSTR pStubRule,
            /* [defaultvalue][in][optional] */ double pRefDate,
            /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
            /* [defaultvalue][in][optional] */ BSTR pDayCount,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CorrelationSmileStrike )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [in] */ VARIANT __RPC_FAR *pVolCurves,
            /* [in] */ VARIANT __RPC_FAR *pProportions,
            /* [in] */ double AsOf,
            /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeLow,
            /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
            /* [in][optional] */ VARIANT __RPC_FAR *pIndexVector,
            /* [defaultvalue][in][optional] */ BSTR Name,
            /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeLow,
            /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeUp,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMFutDelivery )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pFut,
            /* [in] */ BSTR pCcy,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMLivretACurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pAsOf,
            /* [in] */ BSTR pInfCurvId,
            /* [in] */ BSTR pEuribCurvId,
            /* [defaultvalue][in][optional] */ double pFlagRouding,
            /* [defaultvalue][in][optional] */ BSTR pInfResetMgrId,
            /* [defaultvalue][in][optional] */ BSTR pFixingLivretAId,
            /* [defaultvalue][in][optional] */ BSTR pFixingEuribId,
            /* [defaultvalue][in][optional] */ BSTR pMonthForAugust,
            /* [defaultvalue][in][optional] */ BSTR pMonthForFebruary,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeLivretA )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ double pStartDateLeg1,
            /* [in] */ double pEndDateLeg1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexLeg1,
            /* [in] */ VARIANT pSpreadLeg1,
            /* [in] */ BSTR pDayCountLeg1,
            /* [in] */ BSTR pPayFreqLeg1,
            /* [in] */ BSTR pResetFreqLeg1,
            /* [in] */ BSTR pResetTimingLeg1,
            /* [in] */ BSTR pAdjLeg1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pStub,
            /* [in] */ double pEndDateLA,
            /* [in] */ double pSpreadLeg2,
            /* [in] */ BSTR pDayCountLA,
            /* [in] */ BSTR pPayFreqLA,
            /* [in] */ BSTR pResetFreqLA,
            /* [in] */ BSTR pResetTimingLA,
            /* [in] */ BSTR pAdjLA,
            /* [in] */ BSTR pIndexLeg2,
            /* [in] */ BSTR pDayCountLeg2,
            /* [in] */ BSTR pPayFreqLeg2,
            /* [in] */ BSTR pResetFreqLeg2,
            /* [in] */ BSTR pResetTimingLeg2,
            /* [in] */ BSTR pAdjLeg2,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [in] */ BSTR pLAMod,
            /* [in] */ BSTR pLAModBump,
            /* [in] */ BSTR pLAModBumpInflation,
            /* [in] */ VARIANT __RPC_FAR *pResetMgrIds,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMGetFixingFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR source,
            /* [in] */ BSTR index,
            /* [in] */ BSTR term,
            /* [in] */ BSTR ccy,
            /* [in] */ BSTR curveName,
            /* [in] */ DATE date,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_FxConvertFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR ccy1,
            /* [in] */ BSTR ccy2,
            /* [defaultvalue][in][optional] */ BSTR pCvName,
            /* [in] */ DATE pDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FwdSpread )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR defcurveId,
            /* [in] */ double Maturity1,
            /* [in] */ double Maturity2,
            /* [defaultvalue][in][optional] */ double FwdStartDate,
            /* [defaultvalue][in][optional] */ double FwdEndDate,
            /* [defaultvalue][in][optional] */ BSTR VolId,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Flat_Correlation )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ DATE AsOf,
            /* [in] */ BSTR structName,
            /* [in] */ double correlValue,
            /* [defaultvalue][in] */ BSTR idIndex1,
            /* [defaultvalue][in] */ BSTR idIndex2,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CorridorLeg_Sche )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double Notional,
            /* [in] */ BSTR RecieveOrPay,
            /* [in] */ BSTR RefValueSpreadsInBP,
            /* [in] */ BSTR floatingIdx,
            /* [in] */ double leverageFloatIdx,
            /* [in] */ BSTR creditIdx,
            /* [in] */ BSTR refvalueKUPinBP,
            /* [in] */ BSTR refvalueKDWinBP,
            /* [in] */ BSTR ScheduleInfoId,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *accondef,
            /* [defaultvalue][in][optional] */ BSTR disc_ccy,
            /* [defaultvalue][in][optional] */ BSTR Name,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Schedule_Info )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double EffectiveDate,
            /* [in] */ double MaturityDate,
            /* [defaultvalue][in][optional] */ BSTR payFrequency,
            /* [defaultvalue][in][optional] */ BSTR ResetFreq,
            /* [defaultvalue][in][optional] */ BSTR DayCount,
            /* [defaultvalue][in][optional] */ BSTR Stubrule,
            /* [defaultvalue][in][optional] */ BSTR intRule,
            /* [defaultvalue][in][optional] */ BSTR payCalName,
            /* [defaultvalue][in][optional] */ BSTR PayTiming,
            /* [defaultvalue][in][optional] */ BSTR ResetTiming,
            /* [defaultvalue][in][optional] */ BSTR fwdRule,
            /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
            /* [defaultvalue][in][optional] */ BSTR adj,
            /* [defaultvalue][in][optional] */ BSTR intStartAdj,
            /* [defaultvalue][in][optional] */ BSTR AccDayCount,
            /* [defaultvalue][in][optional] */ double ReferenceDate,
            /* [defaultvalue][in][optional] */ double FirstCpnEffDate,
            /* [defaultvalue][in][optional] */ BSTR AdjCal,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMInfCurveSetResetMgr )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pInfCurve,
            /* [in] */ BSTR pResetMgr,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CPDO )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pRiskyLeg,
            /* [in] */ BSTR pRollLeg,
            /* [in] */ BSTR pNoRiskyLeg,
            /* [in] */ double pInitialValo,
            /* [in] */ double pTarget,
            /* [in] */ double pMaturity,
            /* [in] */ BSTR pCpnType,
            /* [in] */ double pUFFees,
            /* [in] */ double pRunningFees,
            /* [in] */ double pVExpo,
            /* [in] */ double pV0Expo,
            /* [in] */ double pAlpha,
            /* [in] */ double pBeta,
            /* [in] */ double pDesactivation,
            /* [in] */ int pNbAssets,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_PriceVector )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR pCPTTYPE,
            /* [out] */ VARIANT __RPC_FAR *pRetVectorValos);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenPrice )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pPricer,
            /* [in] */ BSTR pCPTTYPE,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeTxFixed )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pNtl,
            /* [in] */ double pStartDateLeg1,
            /* [in] */ double pEndDateLeg1,
            /* [in] */ BSTR pCcy,
            /* [in] */ BSTR pIndexLeg1,
            /* [in] */ VARIANT pSpreadLeg1,
            /* [in] */ BSTR pDayCountLeg1,
            /* [in] */ BSTR pPayFreqLeg1,
            /* [in] */ BSTR pResetFreqLeg1,
            /* [in] */ BSTR pResetTimingLeg1,
            /* [in] */ BSTR pAdjLeg1,
            /* [in] */ BSTR pRoll,
            /* [in] */ BSTR pStub,
            /* [in] */ double pEndDateFixed,
            /* [in] */ double pSpreadLeg2,
            /* [in] */ BSTR pDayCountFixed,
            /* [in] */ BSTR pPayFreqFixed,
            /* [in] */ BSTR pResetFreqFixed,
            /* [in] */ BSTR pResetTimingFixed,
            /* [in] */ BSTR pAdjFixed,
            /* [in] */ BSTR pIndexLeg2,
            /* [in] */ BSTR pDayCountLeg2,
            /* [in] */ BSTR pPayFreqLeg2,
            /* [in] */ BSTR pResetFreqLeg2,
            /* [in] */ BSTR pResetTimingLeg2,
            /* [in] */ BSTR pAdjLeg2,
            /* [in] */ double pEndDateAmort,
            /* [in] */ BSTR pDayCountAmort,
            /* [in] */ BSTR pIntRuleAmort,
            /* [in] */ double pTxAmort,
            BSTR pFreqAmort,
            /* [in] */ double pAmountAmort,
            /* [in] */ BSTR pTypeAmort,
            /* [in] */ double pFee,
            /* [in] */ BSTR pSmiledMod,
            /* [in] */ BSTR pSmiledModBump,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FixingCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pDates,
            /* [in] */ VARIANT __RPC_FAR *pValues,
            /* [defaultvalue][in][optional] */ double AsOfDate,
            /* [defaultvalue][in][optional] */ BSTR B_IndexName,
            /* [defaultvalue][in][optional] */ BSTR B_IndexID,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CptInterpolDefCurveOLD )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurve,
            /* [in] */ BSTR pTenor,
            /* [in] */ double pSlope,
            /* [in] */ double pDate,
            /* [defaultvalue][in][optional] */ double pInterpDate,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_CreateBasketCorrelMkDataFromCalypso )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricingEnv,
            /* [in] */ double date,
            /* [in] */ BSTR forceCurveName,
            /* [in] */ BSTR Ccy,
            /* [defaultvalue][in][optional] */ BSTR xmlFileName,
            /* [defaultvalue][in][optional] */ BSTR indexId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSetCalendar )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pFileName,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMInitGigaSpaces )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pUrl,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GetExpectedLoss )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pricerId,
            /* [in] */ double YearTerm,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_VariableCollateral )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT pLabels,
            /* [in] */ VARIANT pNotionals,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_IndexCompo )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR IndexName,
            /* [in] */ VARIANT __RPC_FAR *pLabels,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *YearFrac,
            /* [defaultvalue][in] */ VARIANT __RPC_FAR *pSpread,
            /* [defaultvalue][in] */ BSTR Method,
            /* [defaultvalue][in] */ BSTR Basis,
            /* [defaultvalue][in] */ BSTR ResetFreq,
            /* [defaultvalue][in] */ BSTR PayFreq,
            /* [defaultvalue][in] */ BSTR ccy,
            /* [defaultvalue][in] */ BSTR fwdRule,
            /* [defaultvalue][in] */ BSTR resetTiming,
            /* [defaultvalue][in] */ int resetGap,
            /* [defaultvalue][in] */ BSTR payTiming,
            /* [defaultvalue][in] */ int payGap,
            /* [defaultvalue][in] */ BSTR intRule,
            /* [defaultvalue][in] */ BSTR AdjCalType,
            /* [defaultvalue][in] */ int cm_resetWeekDay,
            /* [defaultvalue][in] */ int cm_resetOccur,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_createFlatCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurve,
            /* [in] */ VARIANT pTenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FunctionRegister )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ long address);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBermudanXStyle )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pxDates,
            /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pexpiryDates,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMcomputeCRA )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pFixorFloat,
            /* [in] */ double pFee,
            /* [in] */ double pAsOf,
            /* [in] */ double pStartDate,
            /* [in] */ double pEndDate,
            /* [in] */ BSTR pCcy,
            /* [in] */ double pLevelUp,
            /* [in] */ BSTR pUpSpec,
            /* [in] */ double pLevelDown,
            /* [in] */ BSTR pDownSpec,
            /* [in] */ BSTR pRefIndex,
            /* [in] */ BSTR pDayCount,
            /* [in] */ BSTR pPayFreqPayIndex,
            /* [in] */ BSTR pResetFreqRefIndex,
            /* [in] */ BSTR pPaidRstTiming,
            /* [in] */ BSTR pRefRstTiming,
            /* [in] */ BSTR pStubRule,
            /* [in] */ BSTR pPOrR,
            /* [in] */ double pStartCallDate,
            /* [in] */ BSTR pXStyle,
            /* [in] */ BSTR pFundingIndex,
            /* [in] */ BSTR pResetFreqFunding,
            /* [in] */ BSTR pPayFreqFunding,
            /* [in] */ VARIANT pSpreadFunding,
            /* [in] */ BSTR pPOrRFunding,
            /* [in] */ double pDecompPricingFlag,
            /* [in] */ double pdiscMarginFactor,
            /* [in] */ BSTR pPreInitFlag,
            /* [in] */ double pMeanReversion,
            /* [in] */ VARIANT pCalibParams,
            /* [in] */ VARIANT pCalibParamsPF,
            /* [in] */ double pKernelToGP,
            /* [in] */ VARIANT pMarkovTreeParams,
            /* [in] */ double pMarkovTreePathNumber,
            /* [in] */ BSTR pBsmodId,
            /* [in] */ BSTR pBsmodSwoptId,
            /* [in] */ BSTR pBsmodSwoptBumpId,
            /* [in] */ BSTR pzcId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Math_BivNormale )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double x,
            /* [in] */ double y,
            /* [in] */ double rho,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_GETINSTRUMENTFROMCALYPSO )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR CalypsoId,
            /* [in] */ BSTR Type,
            /* [defaultvalue][optional][in] */ double AsOf,
            /* [defaultvalue][optional][in] */ BSTR ModelType,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMSetDiscountPricingMode )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pModelId,
            /* [in] */ int pDiscountPricingMode,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Math_RandUniform )( 
            IARMModule __RPC_FAR * This,
            /* [defaultvalue][optional][in] */ double seed,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Math_Interpol )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *X,
            /* [in] */ VARIANT __RPC_FAR *Y,
            /* [in] */ double value,
            /* [defaultvalue][in][optional] */ double type,
            /* [defaultvalue][in][optional] */ double smooth,
            /* [in][optional] */ VARIANT __RPC_FAR *Weights,
            /* [defaultvalue][in][optional] */ double modeSpline,
            /* [defaultvalue][in][optional] */ double withC1condition,
            /* [defaultvalue][in][optional] */ double leftSlope,
            /* [defaultvalue][in][optional] */ double rightSlope,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Random_Generator )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR RandomType,
            /* [defaultvalue][optional][in] */ BSTR ParamId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenerateOneRandom )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR RandomId,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_GenerateRandoms )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR RandomId,
            /* [in] */ int DimVector,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_ResetRandom )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR RandomId);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_FwdSpreadAsIndex )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR DefCurveId,
            /* [in] */ double matu1,
            /* [in] */ double Matu2,
            /* [retval][out] */ double __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_createDefCurveFromBase )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pCurveCDS,
            /* [in] */ BSTR pCurveIndex,
            /* [in] */ VARIANT vBase,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_RiskyPV01AsSensitivity )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pDefCurve,
            /* [in] */ BSTR Tenor,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_SetVolCurve )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR Model,
            /* [in] */ BSTR VolCurveId,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMPF )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT __RPC_FAR *pinsts,
            /* [in] */ VARIANT __RPC_FAR *pcoeffs,
            /* [in] */ VARIANT __RPC_FAR *pmarketPrices,
            VARIANT __RPC_FAR *pprecisions,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMBondTEC )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double pIssueDate,
            /* [in] */ double pMaturityDate,
            /* [in] */ double pFirstCpnDate,
            /* [in] */ double pCpnRate,
            /* [in] */ double pRedempPrice,
            /* [in] */ long pPeriodicity,
            /* [in] */ VARIANT pDaycount,
            /* [defaultvalue][in][optional] */ long pSettleGap,
            /* [defaultvalue][in][optional] */ long pCpnDateFlag,
            /* [defaultvalue][in][optional] */ BSTR pCcyId,
            /* [defaultvalue][in][optional] */ double ptec,
            /* [defaultvalue][in][optional] */ BSTR pPFTecId,
            /* [defaultvalue][in][optional] */ BSTR pModTecId,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMPFModFit )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR pmodName,
            /* [in] */ BSTR ppf,
            /* [in] */ double psettlement,
            /* [in] */ BSTR pzc,
            /* [in] */ VARIANT __RPC_FAR *pvList,
            /* [in] */ VARIANT __RPC_FAR *pfList,
            /* [defaultvalue][in][optional] */ long nag_algo,
            /* [defaultvalue][in][optional] */ long pstep,
            /* [defaultvalue][in][optional] */ double phorizon,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_Restrikable_CDO )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double UnderlyingMatu,
            /* [in] */ double Expiry,
            /* [in] */ double Strike,
            /* [in] */ int OptionType,
            /* [in] */ BSTR pUnderlying,
            /* [in] */ double Rehauss,
            /* [in] */ double TriggerFreq,
            /* [in] */ int DiffCDO,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_PropertyList )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ VARIANT attrNames,
            /* [in] */ VARIANT attrValues,
            /* [in][optional] */ VARIANT attrTypes,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARM_Credit_DefCurveIntensityPWC )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ double AsOfDate,
            /* [in] */ VARIANT __RPC_FAR *pMatuRates,
            /* [in] */ VARIANT __RPC_FAR *pInputs,
            /* [in] */ double Type,
            /* [in] */ double Recovery,
            /* [in] */ BSTR IRCurveId,
            /* [defaultvalue][in][optional] */ BSTR bCurrency,
            /* [defaultvalue][in][optional] */ BSTR bLabel,
            /* [defaultvalue][in][optional] */ BSTR VolCurveId,
            /* [defaultvalue][in][optional] */ BSTR calibrationAlgo,
            /* [defaultvalue][in][optional] */ int lag,
            /* [retval][out] */ VARIANT __RPC_FAR *pRet);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ARMTMLeg )( 
            IARMModule __RPC_FAR * This,
            /* [in] */ BSTR ptmIxType,
            /* [in] */ double pstartDate,
            /* [in] */ double pendDate,
            /* [in] */ BSTR pPorR,
            /* [defaultvalue][in][optional] */ double pspread,
            /* [defaultvalue][in][optional] */ BSTR ppayFrequency,
            /* [defaultvalue][in][optional] */ BSTR presetFrequency,
            /* [defaultvalue][in][optional] */ BSTR pinterestRule,
            /* [defaultvalue][in][optional] */ BSTR pfwdRule,
            /* [defaultvalue][in][optional] */ BSTR pstubRule,
            /* [in] */ BSTR pccy,
            /* [retval][out] */ BSTR __RPC_FAR *pRet);
        
        END_INTERFACE
    } IARMModuleVtbl;

    interface IARMModule
    {
        CONST_VTBL struct IARMModuleVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IARMModule_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IARMModule_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IARMModule_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IARMModule_GetTypeInfoCount(This,pctinfo)	\
    (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo)

#define IARMModule_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo)

#define IARMModule_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)

#define IARMModule_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)


#define IARMModule_ARMNextBusinessDay(This,pDate,pCalendrier,pNbDays,pDate2)	\
    (This)->lpVtbl -> ARMNextBusinessDay(This,pDate,pCalendrier,pNbDays,pDate2)

#define IARMModule_ARMAdjustToBusDate(This,pDate,pCalendrier,pRule,pDate2)	\
    (This)->lpVtbl -> ARMAdjustToBusDate(This,pDate,pCalendrier,pRule,pDate2)

#define IARMModule_ARMFreeObject(This,pId,pRet)	\
    (This)->lpVtbl -> ARMFreeObject(This,pId,pRet)

#define IARMModule_ARMIsBusinessDay(This,pDate,pCalendrier,pRet)	\
    (This)->lpVtbl -> ARMIsBusinessDay(This,pDate,pCalendrier,pRet)

#define IARMModule_ARMGetZCFromSummit(This,pIndex,pCurrency,pCvName,pDate,pInterpMethod,pRet)	\
    (This)->lpVtbl -> ARMGetZCFromSummit(This,pIndex,pCurrency,pCvName,pDate,pInterpMethod,pRet)

#define IARMModule_ARMFreeAllObjects(This,pRet)	\
    (This)->lpVtbl -> ARMFreeAllObjects(This,pRet)

#define IARMModule_ARMYcMod(This,pZc,pZcDiscount,pRet)	\
    (This)->lpVtbl -> ARMYcMod(This,pZc,pZcDiscount,pRet)

#define IARMModule_ARMForwardYield(This,pZc,pMatu1,pMatu2,pMeth,pAdj,pRet)	\
    (This)->lpVtbl -> ARMForwardYield(This,pZc,pMatu1,pMatu2,pMeth,pAdj,pRet)

#define IARMModule_ARMDiscountYield(This,pZc,pMatu,pMeth,pRet)	\
    (This)->lpVtbl -> ARMDiscountYield(This,pZc,pMatu,pMeth,pRet)

#define IARMModule_ARMLiborSwap(This,pStartDate,pEndDate,pLiborType,pRecOrPay,pFixedRate,pSpread,pCcy,pDaycount,pFloatingDaycount,pRet)	\
    (This)->lpVtbl -> ARMLiborSwap(This,pStartDate,pEndDate,pLiborType,pRecOrPay,pFixedRate,pSpread,pCcy,pDaycount,pFloatingDaycount,pRet)

#define IARMModule_ARMSwapPriceToRate(This,pSwap,pDate,pPrice,pModel,pRet)	\
    (This)->lpVtbl -> ARMSwapPriceToRate(This,pSwap,pDate,pPrice,pModel,pRet)

#define IARMModule_ARMPrice(This,pSec,pModel,pRet)	\
    (This)->lpVtbl -> ARMPrice(This,pSec,pModel,pRet)

#define IARMModule_ARMBetweenDates(This,pDate1,pDate2,pDaycount,pIsYearFrac,pRet)	\
    (This)->lpVtbl -> ARMBetweenDates(This,pDate1,pDate2,pDaycount,pIsYearFrac,pRet)

#define IARMModule_ARMAddPeriod(This,pDate,pFreq,pCcy,pNbPeriods,pAdjRule,pRet)	\
    (This)->lpVtbl -> ARMAddPeriod(This,pDate,pFreq,pCcy,pNbPeriods,pAdjRule,pRet)

#define IARMModule_ARMIsoCcy(This,pCcy,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMIsoCcy(This,pCcy,pRefObj,pRet)

#define IARMModule_ARMGetSpotDays(This,pCcy,pRet)	\
    (This)->lpVtbl -> ARMGetSpotDays(This,pCcy,pRet)

#define IARMModule_ARMGetLiborIndexDaycount(This,pCcy,pRet)	\
    (This)->lpVtbl -> ARMGetLiborIndexDaycount(This,pCcy,pRet)

#define IARMModule_ARMGetLiborTerm(This,pCcy,pRet)	\
    (This)->lpVtbl -> ARMGetLiborTerm(This,pCcy,pRet)

#define IARMModule_ARMGetFixedDayCount(This,pCcy,pRet)	\
    (This)->lpVtbl -> ARMGetFixedDayCount(This,pCcy,pRet)

#define IARMModule_ARMGetFixedPayFreq(This,pCcy,pRet)	\
    (This)->lpVtbl -> ARMGetFixedPayFreq(This,pCcy,pRet)

#define IARMModule_ARMComputeVolatility(This,pVol,pmatu,pStrike,pTenor,pRet)	\
    (This)->lpVtbl -> ARMComputeVolatility(This,pVol,pmatu,pStrike,pTenor,pRet)

#define IARMModule_ARMVolCurv(This,pMatu,pStrike,pVol,pAsOf,pStrikeType,pVolType,pCcy,pIndexId,pRet)	\
    (This)->lpVtbl -> ARMVolCurv(This,pMatu,pStrike,pVol,pAsOf,pStrikeType,pVolType,pCcy,pIndexId,pRet)

#define IARMModule_ARMGetVolCubeFromSummit(This,pIndex,pCcy,pCvName,pDate,pType,pSmiles,pTypeCube,indexId,pRet)	\
    (This)->lpVtbl -> ARMGetVolCubeFromSummit(This,pIndex,pCcy,pCvName,pDate,pType,pSmiles,pTypeCube,indexId,pRet)

#define IARMModule_ARM_FxConvert(This,pccy1,pccy2,pDate,pCvName,pRet)	\
    (This)->lpVtbl -> ARM_FxConvert(This,pccy1,pccy2,pDate,pCvName,pRet)

#define IARMModule_ARM_DiscountPrice(This,pcurve,pmatu,pRet)	\
    (This)->lpVtbl -> ARM_DiscountPrice(This,pcurve,pmatu,pRet)

#define IARMModule_ARM_Credit_Delivery(This,pAsOfDate,pTenorContract,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Delivery(This,pAsOfDate,pTenorContract,pRet)

#define IARMModule_ARM_Credit_CptInterpolDefCurve(This,pCurve,pTenor,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CptInterpolDefCurve(This,pCurve,pTenor,pRet)

#define IARMModule_ARM_Credit_DefaultProba(This,pCurve,pMatu,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DefaultProba(This,pCurve,pMatu,pRet)

#define IARMModule_ARM_Credit_GetBeta(This,pPricer,pLabel,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetBeta(This,pPricer,pLabel,pRet)

#define IARMModule_ARM_Credit_Mezzanine(This,pEffectiveDate,pEndDate,pSpread,pMezzAmount,pSubAmount,pLabels,pNotionals,pFreqFeeLeg,pDayCount,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,IncludeMatu,pFstCpnEffDate,intRule,adjStartDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Mezzanine(This,pEffectiveDate,pEndDate,pSpread,pMezzAmount,pSubAmount,pLabels,pNotionals,pFreqFeeLeg,pDayCount,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,IncludeMatu,pFstCpnEffDate,intRule,adjStartDate,pRet)

#define IARMModule_ARM_Credit_ModelMultiCurves(This,pIRcurve,pDefCurves,pRecoveryRates,CorrelId,VolCurve,CpnInfCurve,CpnIRCurve,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ModelMultiCurves(This,pIRcurve,pDefCurves,pRecoveryRates,CorrelId,VolCurve,CpnInfCurve,CpnIRCurve,pRet)

#define IARMModule_ARM_Credit_FTD(This,pEffectiveDate,pEndDate,pSpread,pLabels,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pIssuerNotional,pAccruedOnDefault,pCurrency,pPayCreditLag,pStub,pFstCpnEffDate,pintRule,pstartAdj,pRet)	\
    (This)->lpVtbl -> ARM_Credit_FTD(This,pEffectiveDate,pEndDate,pSpread,pLabels,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pIssuerNotional,pAccruedOnDefault,pCurrency,pPayCreditLag,pStub,pFstCpnEffDate,pintRule,pstartAdj,pRet)

#define IARMModule_ARM_Credit_NTD(This,pEffectiveDate,pEndDate,pSpread,pFirstNumDefault,pLastNumDefault,pLabels,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pIssuerNotional,pAccruedOnDefault,pCurrency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,IncludeMatu,pFstCpnEffDate,intRule,startAdj,pRet)	\
    (This)->lpVtbl -> ARM_Credit_NTD(This,pEffectiveDate,pEndDate,pSpread,pFirstNumDefault,pLastNumDefault,pLabels,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pIssuerNotional,pAccruedOnDefault,pCurrency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,IncludeMatu,pFstCpnEffDate,intRule,startAdj,pRet)

#define IARMModule_ARM_Credit_Price(This,pPricer,pAsofDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Price(This,pPricer,pAsofDate,pRet)

#define IARMModule_ARM_Credit_RiskyDuration(This,pDefCurve,date,pRet)	\
    (This)->lpVtbl -> ARM_Credit_RiskyDuration(This,pDefCurve,date,pRet)

#define IARMModule_ARM_Credit_CDONPV(This,pPricer,pCPTTYPE,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDONPV(This,pPricer,pCPTTYPE,pRet)

#define IARMModule_ARM_Credit_CorrMatrix(This,pLabels,pCoefs,AsOf,Name,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CorrMatrix(This,pLabels,pCoefs,AsOf,Name,pRet)

#define IARMModule_ARM_Credit_ExtractCorrMatrix(This,pCorrMatrixId,pLabels,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ExtractCorrMatrix(This,pCorrMatrixId,pLabels,pRet)

#define IARMModule_ARM_Credit_Spread(This,pPricer,pMTM,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Spread(This,pPricer,pMTM,pRet)

#define IARMModule_ARM_Credit_SetLabel(This,pCurveId,pLabel,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetLabel(This,pCurveId,pLabel,pRet)

#define IARMModule_ARM_Credit_GetLabel(This,pCurveId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetLabel(This,pCurveId,pRet)

#define IARMModule_ARM_Credit_Sensitivity(This,pPricerId,pType,pPlot,pLabel,pEpsilon,epsilonGamma,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Sensitivity(This,pPricerId,pType,pPlot,pLabel,pEpsilon,epsilonGamma,pRet)

#define IARMModule_ARM_Credit_GetCleanSpreadTranche(This,pPricerId,pPlot,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetCleanSpreadTranche(This,pPricerId,pPlot,pRet)

#define IARMModule_ARM_Credit_GetDefProbTranche(This,PricerId,Yearterm,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetDefProbTranche(This,PricerId,Yearterm,pRet)

#define IARMModule_ARM_Credit_GetDuration(This,pPricerId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetDuration(This,pPricerId,pRet)

#define IARMModule_ARM_Credit_GenSchedule(This,pAccStartDate,pAccEndDate,pFixingFreq,pDayCountFrq,prefDate,pCurrency,ptypeDates,pModFol,pCreditGap,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenSchedule(This,pAccStartDate,pAccEndDate,pFixingFreq,pDayCountFrq,prefDate,pCurrency,ptypeDates,pModFol,pCreditGap,pRet)

#define IARMModule_ARM_Credit_CashFlows(This,pCoefs,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CashFlows(This,pCoefs,pRet)

#define IARMModule_ARM_View(This,pObjet,pRet)	\
    (This)->lpVtbl -> ARM_View(This,pObjet,pRet)

#define IARMModule_ARM_Credit_Version(This,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Version(This,pRet)

#define IARMModule_ARM_Credit_CDS(This,pEffectiveDate,pEndDate,pSpread,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pFixedPayerAmount,pFloatingPayerAmount,StubRule,pCurrency,Adjusted,CreditDefLag,IncludeMatu,StartProtection,EndProtection,name,binary,pFstCpnEffDate,StartAdj,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDS(This,pEffectiveDate,pEndDate,pSpread,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pFixedPayerAmount,pFloatingPayerAmount,StubRule,pCurrency,Adjusted,CreditDefLag,IncludeMatu,StartProtection,EndProtection,name,binary,pFstCpnEffDate,StartAdj,pRet)

#define IARMModule_ARMGetVolFromSummit(This,pIndex,pCcy,pCvName,pDate,pType,pMatIndex,pImpOrHist,pindexId,pRet)	\
    (This)->lpVtbl -> ARMGetVolFromSummit(This,pIndex,pCcy,pCvName,pDate,pType,pMatIndex,pImpOrHist,pindexId,pRet)

#define IARMModule_ARMParallelShift(This,pZc,pBump,pRet)	\
    (This)->lpVtbl -> ARMParallelShift(This,pZc,pBump,pRet)

#define IARMModule_ARMBumpVolatility(This,pVolCurv,pBump,pNthLine,pNthCol,pIsCumul,pRet)	\
    (This)->lpVtbl -> ARMBumpVolatility(This,pVolCurv,pBump,pNthLine,pNthCol,pIsCumul,pRet)

#define IARMModule_ARMBsSmiledModel(This,pDate,pSpot,pDividend,pDiscrate,pVolATM,pRo,pNu,pIsSABR,pBeta,pRet)	\
    (This)->lpVtbl -> ARMBsSmiledModel(This,pDate,pSpot,pDividend,pDiscrate,pVolATM,pRo,pNu,pIsSABR,pBeta,pRet)

#define IARMModule_ARMSetEtoolkit(This,pUserName,pPassWord,pDatabaseContext,pItConfigDomainDir,pItDomainName,pRet)	\
    (This)->lpVtbl -> ARMSetEtoolkit(This,pUserName,pPassWord,pDatabaseContext,pItConfigDomainDir,pItDomainName,pRet)

#define IARMModule_ARMConnectionEtoolkit(This,pRet)	\
    (This)->lpVtbl -> ARMConnectionEtoolkit(This,pRet)

#define IARMModule_ARMVolFlat(This,pVol,pDate,pCcy,pRet)	\
    (This)->lpVtbl -> ARMVolFlat(This,pVol,pDate,pCcy,pRet)

#define IARMModule_ARMVolCube(This,pATMVol,pSmileCurveIds,pTenors,pVolType,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMVolCube(This,pATMVol,pSmileCurveIds,pTenors,pVolType,pRefObj,pRet)

#define IARMModule_ARMDeconnectionEtoolkit(This,pRet)	\
    (This)->lpVtbl -> ARMDeconnectionEtoolkit(This,pRet)

#define IARMModule_ARMZcFlat(This,pZc,pDate,pCcy,pRet)	\
    (This)->lpVtbl -> ARMZcFlat(This,pZc,pDate,pCcy,pRet)

#define IARMModule_ARMBsModel(This,pDate,pSpot,pDividend,pDiscrate,pVol,pTypeStk,pRet)	\
    (This)->lpVtbl -> ARMBsModel(This,pDate,pSpot,pDividend,pDiscrate,pVol,pTypeStk,pRet)

#define IARMModule_ARMSwitchToETK(This)	\
    (This)->lpVtbl -> ARMSwitchToETK(This)

#define IARMModule_ARMSwitchToFLATFILE(This)	\
    (This)->lpVtbl -> ARMSwitchToFLATFILE(This)

#define IARMModule_ARMZCLINT(This,pMatu,pRate,pMeth,pDate,pCurrency,pInterpMethod,pRet)	\
    (This)->lpVtbl -> ARMZCLINT(This,pMatu,pRate,pMeth,pDate,pCurrency,pInterpMethod,pRet)

#define IARMModule_ARM_Credit_Pricer(This,pSecurity,pModel,pPricerType,l_nbpaths,pParameters,valuationdate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Pricer(This,pSecurity,pModel,pPricerType,l_nbpaths,pParameters,valuationdate,pRet)

#define IARMModule_ARM_zcspreaded(This,zcSprId,zcInitId,date,MMFreq,SwapFreq,ccyId,pRet)	\
    (This)->lpVtbl -> ARM_zcspreaded(This,zcSprId,zcInitId,date,MMFreq,SwapFreq,ccyId,pRet)

#define IARMModule_ARM_Credit_CptBaseCorrelation(This,AsOf,name,CalMethod,IndexId,pStrikeLow,pStrikeHigh,pVMktBid,pVMktAsk,pVUpfBid,pVUpfAsk,pVInitialCorrel,pVDeltaLevrage,ModelId,integrationStep,lagStartDate,creditLag,pVectorPrevIndexId,pMatrixPrevBC,step,CalMeth,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CptBaseCorrelation(This,AsOf,name,CalMethod,IndexId,pStrikeLow,pStrikeHigh,pVMktBid,pVMktAsk,pVUpfBid,pVUpfAsk,pVInitialCorrel,pVDeltaLevrage,ModelId,integrationStep,lagStartDate,creditLag,pVectorPrevIndexId,pMatrixPrevBC,step,CalMeth,pRet)

#define IARMModule_ARM_Credit_CDO2(This,pEffectiveDate,pEndDate,pPortfolio,pSpread,pSubAmount,pMezzAmount,pFreqFeeLeg,pFreqDefLeg,pDayCountFrq,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pBinary,pPayCal,LongOrShortRisk,TradedNotional,CrossSub,IncludeMatu,pFstCpnEffDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDO2(This,pEffectiveDate,pEndDate,pPortfolio,pSpread,pSubAmount,pMezzAmount,pFreqFeeLeg,pFreqDefLeg,pDayCountFrq,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pBinary,pPayCal,LongOrShortRisk,TradedNotional,CrossSub,IncludeMatu,pFstCpnEffDate,pRet)

#define IARMModule_ARM_Credit_Portfolio(This,pSecuritiesID,Parameters,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Portfolio(This,pSecuritiesID,Parameters,pRet)

#define IARMModule_ARMCreateZCSwapInt(This,pDate,pMatu,pRate,pMMVsFut,pSwapVsFut,pRaw,pInterp,pCcy,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMCreateZCSwapInt(This,pDate,pMatu,pRate,pMMVsFut,pSwapVsFut,pRaw,pInterp,pCcy,pRefObj,pRet)

#define IARMModule_ARMGetInitialCurveFromSummit(This,pIndex,pCurrency,pCvName,pDate,pAdjOrNot,pRetMat,pRetRate)	\
    (This)->lpVtbl -> ARMGetInitialCurveFromSummit(This,pIndex,pCurrency,pCvName,pDate,pAdjOrNot,pRetMat,pRetRate)

#define IARMModule_ARMTHREEMONTHFUT(This,pDelivery,pMarket,pCcy,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMTHREEMONTHFUT(This,pDelivery,pMarket,pCcy,pRefObj,pRet)

#define IARMModule_ARMFutPibor(This,pDelivery,pRet)	\
    (This)->lpVtbl -> ARMFutPibor(This,pDelivery,pRet)

#define IARMModule_ARMIRFUT(This,pDelivery,pIdUnderlying,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMIRFUT(This,pDelivery,pIdUnderlying,pRefObj,pRet)

#define IARMModule_ARMLibor(This,pLiborTypeId,pCcyId,pResetFreqId,pPayFreqId,pRefObj,pBasis,pIntrule,pRet)	\
    (This)->lpVtbl -> ARMLibor(This,pLiborTypeId,pCcyId,pResetFreqId,pPayFreqId,pRefObj,pBasis,pIntrule,pRet)

#define IARMModule_ARMLiborSwaption(This,pStartDate,pEndDate,pReceiveOrPay,pStrike,pMaturity,pLiborType,pSpread,pExerciseType,pResetFreq,pPayFreq,pCcyId,pRefObj,pRet)	\
    (This)->lpVtbl -> ARMLiborSwaption(This,pStartDate,pEndDate,pReceiveOrPay,pStrike,pMaturity,pLiborType,pSpread,pExerciseType,pResetFreq,pPayFreq,pCcyId,pRefObj,pRet)

#define IARMModule_ARMFixedLeg(This,pStartDate,pEndDate,pReceiveOrPay,pFixRate,pDayCount,pFreq,pDecompFreq,pPayTiming,pIntRule,pStubRule,pCcyId,pPayCalName,pNxChange,pRefDate,pAdjStartDate,pRet)	\
    (This)->lpVtbl -> ARMFixedLeg(This,pStartDate,pEndDate,pReceiveOrPay,pFixRate,pDayCount,pFreq,pDecompFreq,pPayTiming,pIntRule,pStubRule,pCcyId,pPayCalName,pNxChange,pRefDate,pAdjStartDate,pRet)

#define IARMModule_ARM_Credit_SetCorrelationMatrix(This,pModelMultiCurvesId,pCorrMatrixId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetCorrelationMatrix(This,pModelMultiCurvesId,pCorrMatrixId,pRet)

#define IARMModule_ARM_Credit_DPMktDataFromSummit(This,AsOfDate,Issuer,CurveName,Parameter,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DPMktDataFromSummit(This,AsOfDate,Issuer,CurveName,Parameter,pRet)

#define IARMModule_ARM_Credit_GetDPFromSummit(This,AsOfDate,Issuer,CurveName,ircurveId,label,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetDPFromSummit(This,AsOfDate,Issuer,CurveName,ircurveId,label,pRet)

#define IARMModule_ARM_Credit_CloneCorrMatrixBary(This,CorrMatrixId,Beta,UpOrDown,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CloneCorrMatrixBary(This,CorrMatrixId,Beta,UpOrDown,pRet)

#define IARMModule_ARM_Credit_ConstantDefaultCurve(This,AsOfDate,pTenors,pRates,Recovery,IRCurveId,Ccy,Label,AdjCalType,IsSummit,calibrationData,lag,calibrationAlgo,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ConstantDefaultCurve(This,AsOfDate,pTenors,pRates,Recovery,IRCurveId,Ccy,Label,AdjCalType,IsSummit,calibrationData,lag,calibrationAlgo,pRet)

#define IARMModule_ARM_Credit_DefProbModelNew(This,pDefCurve,pIRcurve,VolCurve,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DefProbModelNew(This,pDefCurve,pIRcurve,VolCurve,pRet)

#define IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromSummit(This,AsOfDate,bIssuer,bCurrency,bCvName,IRCurveId,bLabel,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ZeroCouponDefaultCurveFromSummit(This,AsOfDate,bIssuer,bCurrency,bCvName,IRCurveId,bLabel,pRet)

#define IARMModule_ARMBsSlModel(This,pDate,pZc,pVolSpreadLock,pCvCapVol,pCvIndexVol,pRet)	\
    (This)->lpVtbl -> ARMBsSlModel(This,pDate,pZc,pVolSpreadLock,pCvCapVol,pCvIndexVol,pRet)

#define IARMModule_ARMGetFXVolFromSummit(This,pCcy1,pCcy2,pDate,pCvName,pType,pRet)	\
    (This)->lpVtbl -> ARMGetFXVolFromSummit(This,pCcy1,pCcy2,pDate,pCvName,pType,pRet)

#define IARMModule_ARMGetFXCorrelFromSummit(This,pCcy1,pIndex,pCcy2,pDate,pCvName,pTenors,pRet)	\
    (This)->lpVtbl -> ARMGetFXCorrelFromSummit(This,pCcy1,pIndex,pCcy2,pDate,pCvName,pTenors,pRet)

#define IARMModule_ARMInfocentreConnect(This)	\
    (This)->lpVtbl -> ARMInfocentreConnect(This)

#define IARMModule_ARMGlobDFBS(This,pDomBSId,pDomCurrId,pFrgBSId,pFrgCurrId,pFxVolCrvId,pFFxCorrId,pRatesCorrId,pFxVolModelId,pRet)	\
    (This)->lpVtbl -> ARMGlobDFBS(This,pDomBSId,pDomCurrId,pFrgBSId,pFrgCurrId,pFxVolCrvId,pFFxCorrId,pRatesCorrId,pFxVolModelId,pRet)

#define IARMModule_ARM_Credit_GetInitialCurveFromSummit(This,pIndex,pCurrency,pCvName,pDate,value,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetInitialCurveFromSummit(This,pIndex,pCurrency,pCvName,pDate,value,pRet)

#define IARMModule_ARMGetCorrelFromSummit(This,pCcy1,pIndex1,pCcy2,pIndex2,pDate,pCvName,pRet)	\
    (This)->lpVtbl -> ARMGetCorrelFromSummit(This,pCcy1,pIndex1,pCcy2,pIndex2,pDate,pCvName,pRet)

#define IARMModule_ARMDFFXBS(This,pDVolId,pFVolId,pDZcId,pFZcId,pDFxCorrId,pFFxCorrId,pFxVolId,pRatesCorr,pRet)	\
    (This)->lpVtbl -> ARMDFFXBS(This,pDVolId,pFVolId,pDZcId,pFZcId,pDFxCorrId,pFFxCorrId,pFxVolId,pRatesCorr,pRet)

#define IARMModule_ARMTRIBSMODEL(This,pModel1,pModel2,pDiscModel,pFX1DiscVol,pFX2DiscVol,pIdx1Idx2Corr,pIdx1DiscIdxCorr,pIdx2DiscIdxCorr,pIdx1FxCorr,pIdx2FxCorr,pQuantoFlag,pRet)	\
    (This)->lpVtbl -> ARMTRIBSMODEL(This,pModel1,pModel2,pDiscModel,pFX1DiscVol,pFX2DiscVol,pIdx1Idx2Corr,pIdx1DiscIdxCorr,pIdx2DiscIdxCorr,pIdx1FxCorr,pIdx2FxCorr,pQuantoFlag,pRet)

#define IARMModule_ARMTRIBSDUAL(This,pModel1,pModel2,pDiscModel,pFX1DiscVol,pFX2DiscVol,pIdx1Idx2Corr,pIdx1DiscIdxCorr,pIdx2DiscIdxCorr,pIdx1FxCorr,pIdx2FxCorr,pQuantoFlag,pCorrelForAdj,pWithslopeflag,pRet)	\
    (This)->lpVtbl -> ARMTRIBSDUAL(This,pModel1,pModel2,pDiscModel,pFX1DiscVol,pFX2DiscVol,pIdx1Idx2Corr,pIdx1DiscIdxCorr,pIdx2DiscIdxCorr,pIdx1FxCorr,pIdx2FxCorr,pQuantoFlag,pCorrelForAdj,pWithslopeflag,pRet)

#define IARMModule_ARMcptBonibor(This,pDate,pDepart,pMatStruct,pMatTot,pAmort,pFreq,pSjUSD,pTiming,pBarriere,pSpdPostBar,pMarge,pFunding,pFundingFreq,pSpd2phase,pSoulte,pYcModId,pBsModId,pBsModVolUSDId,pBsModCorrPlusId,pBsModCorrMoinsId,pCrossModId,pProbaMarge,pInt,pRet)	\
    (This)->lpVtbl -> ARMcptBonibor(This,pDate,pDepart,pMatStruct,pMatTot,pAmort,pFreq,pSjUSD,pTiming,pBarriere,pSpdPostBar,pMarge,pFunding,pFundingFreq,pSpd2phase,pSoulte,pYcModId,pBsModId,pBsModVolUSDId,pBsModCorrPlusId,pBsModCorrMoinsId,pCrossModId,pProbaMarge,pInt,pRet)

#define IARMModule_ARM_Credit_CMTranche(This,pEffectiveDate,pEndDate,pParticipationRate,pMezzAmount,pSubAmount,pLabels,pNotionals,pIndex,pFreqFeeLeg,pDayCount,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,FwdFixedDate,IncludeMatu,pFstCpnEffDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CMTranche(This,pEffectiveDate,pEndDate,pParticipationRate,pMezzAmount,pSubAmount,pLabels,pNotionals,pIndex,pFreqFeeLeg,pDayCount,pFirst_period_refdate,pAccruedOnDefault,Currency,pPayCreditLag,pStub,pFreqDefLeg,pBinary,pPayCal,LongOrShortRisk,TradedNotional,FwdFixedDate,IncludeMatu,pFstCpnEffDate,pRet)

#define IARMModule_ARM_Credit_Index(This,pLabels,YearFrac,pSpread,Method,Basis,ResetFreq,PayFreq,ccy,DefCurve,fwdRule,resetTiming,resetGap,payTiming,payGap,intRule,AdjCalType,cm_resetWeekDay,cm_resetOccur,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Index(This,pLabels,YearFrac,pSpread,Method,Basis,ResetFreq,PayFreq,ccy,DefCurve,fwdRule,resetTiming,resetGap,payTiming,payGap,intRule,AdjCalType,cm_resetWeekDay,cm_resetOccur,pRet)

#define IARMModule_ARM_Credit_Parameters(This,pCoefs,nbcols,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Parameters(This,pCoefs,nbcols,pRet)

#define IARMModule_ARMAswPrice(This,pMaturity,pCpn,pFreq,pBase,pMargin,pRedemptionPrice,pAsOf,pDelivery,pFixDecompfreq,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pAmortizationId,pSolve,pMinValue,pMaxValue,pRet)	\
    (This)->lpVtbl -> ARMAswPrice(This,pMaturity,pCpn,pFreq,pBase,pMargin,pRedemptionPrice,pAsOf,pDelivery,pFixDecompfreq,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pAmortizationId,pSolve,pMinValue,pMaxValue,pRet)

#define IARMModule_ARMAswMargin(This,pMaturity,pCpn,pFreq,pBase,pPrice,pRedemptionPrice,pAsOf,pDelivery,pFixDecompfreq,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pAmortizationId,pSolve,pMinValue,pMaxValue,pRet)	\
    (This)->lpVtbl -> ARMAswMargin(This,pMaturity,pCpn,pFreq,pBase,pPrice,pRedemptionPrice,pAsOf,pDelivery,pFixDecompfreq,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pAmortizationId,pSolve,pMinValue,pMaxValue,pRet)

#define IARMModule_ARM_Credit_CDSIndex(This,pEffectiveDate,pEndDate,pSpread,pIndex,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pFixedPayerAmount,pFloatingPayerAmount,StubRule,pCurrency,Adjusted,CreditLag,IncludeMaturity,ProtectionStartDate,ProtectionEndDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDSIndex(This,pEffectiveDate,pEndDate,pSpread,pIndex,pFixingFreq,pDayCountFrq,pFirst_period_refdate,pFixedPayerAmount,pFloatingPayerAmount,StubRule,pCurrency,Adjusted,CreditLag,IncludeMaturity,ProtectionStartDate,ProtectionEndDate,pRet)

#define IARMModule_ARM_Credit_Option(This,UnderlyingMaturity,OptionExpiry,Currency,CdsAdj,EndAdj,pStrike,pOptionType,pKoType,Notional,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Option(This,UnderlyingMaturity,OptionExpiry,Currency,CdsAdj,EndAdj,pStrike,pOptionType,pKoType,Notional,pRet)

#define IARMModule_ARM_Credit_FwdSpreadPricer(This,pPricer,Maturity1,Maturity2,pRet)	\
    (This)->lpVtbl -> ARM_Credit_FwdSpreadPricer(This,pPricer,Maturity1,Maturity2,pRet)

#define IARMModule_ARM_Credit_ImpliedVol(This,pPricer,pMktPrice,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ImpliedVol(This,pPricer,pMktPrice,pRet)

#define IARMModule_ARM_Credit_VirtualCdsSpread(This,pPricer,pMaturity,pRet)	\
    (This)->lpVtbl -> ARM_Credit_VirtualCdsSpread(This,pPricer,pMaturity,pRet)

#define IARMModule_ARM_Credit_BSGreeks(This,pPricer,pGreekType,pRet)	\
    (This)->lpVtbl -> ARM_Credit_BSGreeks(This,pPricer,pGreekType,pRet)

#define IARMModule_ARM_Credit_SetCorrelation(This,pModelMultiCurvesId,pCorrelationId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetCorrelation(This,pModelMultiCurvesId,pCorrelationId,pRet)

#define IARMModule_ARM_Credit_CorrelationStrike(This,pLabels,pVolCurves,pProportions,pSmileStrikeLow,pSmileStrikeHigh,pIndexVector,AsOf,Name,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CorrelationStrike(This,pLabels,pVolCurves,pProportions,pSmileStrikeLow,pSmileStrikeHigh,pIndexVector,AsOf,Name,pRet)

#define IARMModule_ARM_Credit_Beta_Correlation(This,pLabels,pCoefs,AsOf,Name,idIndex1,idIndex2,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Beta_Correlation(This,pLabels,pCoefs,AsOf,Name,idIndex1,idIndex2,pRet)

#define IARMModule_ARM_GETINSTRUMENTFROMSUMMIT(This,SummitId,Type,AsOf,ExoticFilter,pRet)	\
    (This)->lpVtbl -> ARM_GETINSTRUMENTFROMSUMMIT(This,SummitId,Type,AsOf,ExoticFilter,pRet)

#define IARMModule_ARMFrnPrice(This,pAsOf,pDelivery,pMaturity,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pFacialMargin,pValoMargin,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pFixing,pSpread,pOutMode,pSolve,pAmortizationId,pRet)	\
    (This)->lpVtbl -> ARMFrnPrice(This,pAsOf,pDelivery,pMaturity,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pFacialMargin,pValoMargin,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pFixing,pSpread,pOutMode,pSolve,pAmortizationId,pRet)

#define IARMModule_ARMFrnMargin(This,pAsOf,pDelivery,pMaturity,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pFacialMargin,pPrice,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pFixing,pSpread,pOutMode,pSolve,pAmortizationId,pRet)	\
    (This)->lpVtbl -> ARMFrnMargin(This,pAsOf,pDelivery,pMaturity,pCcy1,pIndex1,pFwdCurve1,pDiscCurve1,pFacialMargin,pPrice,pCcy2,pIndex2,pFwdCurve2,pDiscCurve2,pFixing,pSpread,pOutMode,pSolve,pAmortizationId,pRet)

#define IARMModule_ARM_Credit_GetEqStrikeDown(This,CorrelId,IndexName,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetEqStrikeDown(This,CorrelId,IndexName,pRet)

#define IARMModule_ARM_Credit_GetEqStrikeUp(This,CorrelId,IndexName,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetEqStrikeUp(This,CorrelId,IndexName,pRet)

#define IARMModule_ARM_Credit_GetCorrelStrikeDown(This,CorrelId,yfmaturity,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetCorrelStrikeDown(This,CorrelId,yfmaturity,pRet)

#define IARMModule_ARM_Credit_GetCorrelStrikeUp(This,CorrelId,yfmaturity,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetCorrelStrikeUp(This,CorrelId,yfmaturity,pRet)

#define IARMModule_ARM_Credit_GetCorrelation(This,ModelId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetCorrelation(This,ModelId,pRet)

#define IARMModule_ARM_Credit_GetModelFromSummit(This,IRcurve,IDSummit,type,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetModelFromSummit(This,IRcurve,IDSummit,type,pRet)

#define IARMModule_ARM_Credit_SetVolatility(This,pPricer,VolCurveId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetVolatility(This,pPricer,VolCurveId,pRet)

#define IARMModule_ARM_NextCpnDate(This,AsOfDate,maturity,frequency,rule,currency,intrule,pRet)	\
    (This)->lpVtbl -> ARM_NextCpnDate(This,AsOfDate,maturity,frequency,rule,currency,intrule,pRet)

#define IARMModule_ARM_Credit_SetProportionsInfos(This,correlId,IndexName,proportion,forcedstrikelow,forcedstrikehigh)	\
    (This)->lpVtbl -> ARM_Credit_SetProportionsInfos(This,correlId,IndexName,proportion,forcedstrikelow,forcedstrikehigh)

#define IARMModule_ARMcptDigital(This,pAsOf,pStartDate,pEndDate,pNtl,pIndex,pBsmod,pBsmodDelta,pBsmodVega,pFreqP,pResetTiming,pPorR,pCcy,pCcyIdx,pDayCount,pCapOrFloor,pAmort,pStrike,pPayOff,pSpd,pResetGap,pSpreadBelow,pSpreadAbove,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRefDate,pRet)	\
    (This)->lpVtbl -> ARMcptDigital(This,pAsOf,pStartDate,pEndDate,pNtl,pIndex,pBsmod,pBsmodDelta,pBsmodVega,pFreqP,pResetTiming,pPorR,pCcy,pCcyIdx,pDayCount,pCapOrFloor,pAmort,pStrike,pPayOff,pSpd,pResetGap,pSpreadBelow,pSpreadAbove,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRefDate,pRet)

#define IARMModule_ARMRefValue(This,pdates,pvalues,pvalues2,valueType,conversion,calcMethod,pRet)	\
    (This)->lpVtbl -> ARMRefValue(This,pdates,pvalues,pvalues2,valueType,conversion,calcMethod,pRet)

#define IARMModule_ARMCreateGenCorrelManager(This,pMktTags,pIntraMktTags,pCorrelCurveIds,pRet)	\
    (This)->lpVtbl -> ARMCreateGenCorrelManager(This,pMktTags,pIntraMktTags,pCorrelCurveIds,pRet)

#define IARMModule_ARMBSConvAdjust(This,pSUMMITFormulaeUsed,pUseSABRCMS,pRet)	\
    (This)->lpVtbl -> ARMBSConvAdjust(This,pSUMMITFormulaeUsed,pUseSABRCMS,pRet)

#define IARMModule_ARMBsModelGen(This,pYieldCurve,pVolatility,pCorrMgr,pCnvxManager,pCapletVol,pSpreadLock,pDiscCurve,pCorrel,pCashVol,pSpreadVol,pModelType,pSpreadVolType,pSabrMod,pLnorNorVol,pNumSteps,pRet)	\
    (This)->lpVtbl -> ARMBsModelGen(This,pYieldCurve,pVolatility,pCorrMgr,pCnvxManager,pCapletVol,pSpreadLock,pDiscCurve,pCorrel,pCashVol,pSpreadVol,pModelType,pSpreadVolType,pSabrMod,pLnorNorVol,pNumSteps,pRet)

#define IARMModule_ARMcptSPTQTF(This,pAsOf,pStartDate,pStartDatePhase2,pStartDatePhase3,pEndDate,pNtl,pIndexPhase1,pIndexPhase2,pIndexFund,pIndexPhase3,pFreqPPhase1,pFreqPPhase2,pFreqPFund,pFreqPPhase3,pFreqR,pResetTimingPhase1,pResetTimingPhase2,pResetTimingPhase3,pCcy,pCcyIdx,pDayCount,pFee,pIsRateFixedPhase2,pFixedRatePhase2,pBarrier,pSpdPhase1,pSpdPhase1Fund,pSpdPhase2Tf,pSpdPhase2fund,pSpdPhase3,pSpdPhase3fund,pResetGapPhase1,pResetGapPhase2,pResetGapPhase3,pAmort,pBsmod,pBsmodDeltaCcy1,pBsmodVegaCcy1,pBsmodDeltaCcy2,pBsmodVegaCcy2,pBsmodFxCorrel,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRefDate,pRet)	\
    (This)->lpVtbl -> ARMcptSPTQTF(This,pAsOf,pStartDate,pStartDatePhase2,pStartDatePhase3,pEndDate,pNtl,pIndexPhase1,pIndexPhase2,pIndexFund,pIndexPhase3,pFreqPPhase1,pFreqPPhase2,pFreqPFund,pFreqPPhase3,pFreqR,pResetTimingPhase1,pResetTimingPhase2,pResetTimingPhase3,pCcy,pCcyIdx,pDayCount,pFee,pIsRateFixedPhase2,pFixedRatePhase2,pBarrier,pSpdPhase1,pSpdPhase1Fund,pSpdPhase2Tf,pSpdPhase2fund,pSpdPhase3,pSpdPhase3fund,pResetGapPhase1,pResetGapPhase2,pResetGapPhase3,pAmort,pBsmod,pBsmodDeltaCcy1,pBsmodVegaCcy1,pBsmodDeltaCcy2,pBsmodVegaCcy2,pBsmodFxCorrel,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRefDate,pRet)

#define IARMModule_ARMDisplaySchedule(This,legId,dataType,pRet)	\
    (This)->lpVtbl -> ARMDisplaySchedule(This,legId,dataType,pRet)

#define IARMModule_ARMIrIndex(This,pDaycount,pPayFreq,pMaturity,pCompMethod,pFwdRule,pResetTiming,pResetGap,pPayTiming,pPayGap,pCcy,pIndexType,pDecompFreq,pIntRule,pResetFreq,pRet)	\
    (This)->lpVtbl -> ARMIrIndex(This,pDaycount,pPayFreq,pMaturity,pCompMethod,pFwdRule,pResetTiming,pResetGap,pPayTiming,pPayGap,pCcy,pIndexType,pDecompFreq,pIntRule,pResetFreq,pRet)

#define IARMModule_ARMSwapleg(This,pIndexId,pStartDate,pEndDate,pRecOrPay,pSpread,pCcy,pDayCount,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pRet)	\
    (This)->lpVtbl -> ARMSwapleg(This,pIndexId,pStartDate,pEndDate,pRecOrPay,pSpread,pCcy,pDayCount,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pRet)

#define IARMModule_ARMConstRefvalue(This,pValue,pRet)	\
    (This)->lpVtbl -> ARMConstRefvalue(This,pValue,pRet)

#define IARMModule_ARM_Credit_CptImplCvForCDO2(This,pricerId,Name,Tenor,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CptImplCvForCDO2(This,pricerId,Name,Tenor,pRet)

#define IARMModule_ARM_Credit_AddPeriod(This,pAsOf,Maturity,pCcy,AdjRule,AdjCDS,pRet)	\
    (This)->lpVtbl -> ARM_Credit_AddPeriod(This,pAsOf,Maturity,pCcy,AdjRule,AdjCDS,pRet)

#define IARMModule_ARMShutdownETK(This)	\
    (This)->lpVtbl -> ARMShutdownETK(This)

#define IARMModule_ARM_Credit_SetCoupons(This,CdsorCdoId,CouponsId,TypesId,PartId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetCoupons(This,CdsorCdoId,CouponsId,TypesId,PartId,pRet)

#define IARMModule_ARM_Credit_InputDefaultCurve(This,AsOfDate,pDates,pRates,Recovery,IRCurveId,bCurrency,bLabel,bInterpolType,pRet)	\
    (This)->lpVtbl -> ARM_Credit_InputDefaultCurve(This,AsOfDate,pDates,pRates,Recovery,IRCurveId,bCurrency,bLabel,bInterpolType,pRet)

#define IARMModule_ARMGetInitialVolFromSummit(This,pIndex,pCurrency,pCvName,pDate,pType,pMatIndex,pRetMat,pRetTenor,pRetVol)	\
    (This)->lpVtbl -> ARMGetInitialVolFromSummit(This,pIndex,pCurrency,pCvName,pDate,pType,pMatIndex,pRetMat,pRetTenor,pRetVol)

#define IARMModule_ARMBond(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedempPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcy,pRet)	\
    (This)->lpVtbl -> ARMBond(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedempPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcy,pRet)

#define IARMModule_ARMINFCreateOATLeg(This,pStartDate,pEndDate,pInfIdx,pRcvOrPay,pInterpType,pLeverage,pSpread,pResetFreq,pDaycount,pResetCal,pFwdRule,pIntRule,pStubRule,pResetNumGap,pResetDenomGap,pPayFreq,pPayGap,pPayCal,pFinalNotionalType,pFirstReset,pCoMultiple,pRet)	\
    (This)->lpVtbl -> ARMINFCreateOATLeg(This,pStartDate,pEndDate,pInfIdx,pRcvOrPay,pInterpType,pLeverage,pSpread,pResetFreq,pDaycount,pResetCal,pFwdRule,pIntRule,pStubRule,pResetNumGap,pResetDenomGap,pPayFreq,pPayGap,pPayCal,pFinalNotionalType,pFirstReset,pCoMultiple,pRet)

#define IARMModule_ARMSwap(This,pSwapleg1,pSwapleg2,pMinPay,pRet)	\
    (This)->lpVtbl -> ARMSwap(This,pSwapleg1,pSwapleg2,pMinPay,pRet)

#define IARMModule_ARMPToYield(This,pBond,pSettleDate,pPrice,pRet)	\
    (This)->lpVtbl -> ARMPToYield(This,pBond,pSettleDate,pPrice,pRet)

#define IARMModule_ARMYToPrice(This,pBond,pSettleDate,pYield,pRet)	\
    (This)->lpVtbl -> ARMYToPrice(This,pBond,pSettleDate,pYield,pRet)

#define IARMModule_ARMYToDuration(This,pBond,pSettleDate,pActuRate,pFlagCpn,pRet)	\
    (This)->lpVtbl -> ARMYToDuration(This,pBond,pSettleDate,pActuRate,pFlagCpn,pRet)

#define IARMModule_ARMLiborleg(This,pStartDate,pEndDate,pLiborType,pRecOrPay,pSpread,pResetFReq,pPayFreq,pResetTiming,pPayTiming,pCcy,pIntRule,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pCpnDaycount,pRet)	\
    (This)->lpVtbl -> ARMLiborleg(This,pStartDate,pEndDate,pLiborType,pRecOrPay,pSpread,pResetFReq,pPayFreq,pResetTiming,pPayTiming,pCcy,pIntRule,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pCpnDaycount,pRet)

#define IARMModule_ARMImpliedSpread(This,pSwap,pModel,pPrice,pLeg1or2,pRet)	\
    (This)->lpVtbl -> ARMImpliedSpread(This,pSwap,pModel,pPrice,pLeg1or2,pRet)

#define IARMModule_ARMDiscountPrice(This,pZeroCurve,pMatu,pRet)	\
    (This)->lpVtbl -> ARMDiscountPrice(This,pZeroCurve,pMatu,pRet)

#define IARMModule_ARMINFCreateCurve(This,pAsOf,pIndexName,pCPIIndexValue,pCPIIndexDate,pMatu,pRate,pMonthlyInterpType,pDailyInterpType,pDCFMonthly,pDCFDaily,pExtrapolType,pResetManager,pSeasonManager,pRet)	\
    (This)->lpVtbl -> ARMINFCreateCurve(This,pAsOf,pIndexName,pCPIIndexValue,pCPIIndexDate,pMatu,pRate,pMonthlyInterpType,pDailyInterpType,pDCFMonthly,pDCFDaily,pExtrapolType,pResetManager,pSeasonManager,pRet)

#define IARMModule_ARMINFInterpCPI(This,pZc,pCPIDate,pDCFlag,pDailyInterpType,pCPIlag,pWeight,pRet)	\
    (This)->lpVtbl -> ARMINFInterpCPI(This,pZc,pCPIDate,pDCFlag,pDailyInterpType,pCPIlag,pWeight,pRet)

#define IARMModule_ARMINFSeasonManager(This,pMonthList,pValues,pSeasonAdjMode,pRet)	\
    (This)->lpVtbl -> ARMINFSeasonManager(This,pMonthList,pValues,pSeasonAdjMode,pRet)

#define IARMModule_ARMINFResetManager(This,pDatas,pNbIndex,pRet)	\
    (This)->lpVtbl -> ARMINFResetManager(This,pDatas,pNbIndex,pRet)

#define IARMModule_ARMINFYcMod(This,pYieldCurve,pInfCurve,pRet)	\
    (This)->lpVtbl -> ARMINFYcMod(This,pYieldCurve,pInfCurve,pRet)

#define IARMModule_ARMBaseReplicationConnect(This)	\
    (This)->lpVtbl -> ARMBaseReplicationConnect(This)

#define IARMModule_ARMGetInitialFXVolFromSummit(This,pCcy1,pCcy2,pDate,pCvName,pImpOrHist,pVolType,pRetMat,pRetTenor,pRetVol)	\
    (This)->lpVtbl -> ARMGetInitialFXVolFromSummit(This,pCcy1,pCcy2,pDate,pCvName,pImpOrHist,pVolType,pRetMat,pRetTenor,pRetVol)

#define IARMModule_ARMCreateZCFromSummit(This,pIndex,pCurrency,pCvName,pDate,pAdj,pRaw,pRet)	\
    (This)->lpVtbl -> ARMCreateZCFromSummit(This,pIndex,pCurrency,pCvName,pDate,pAdj,pRaw,pRet)

#define IARMModule_ARMBumpCurve(This,pZc,pEpsilon,pMethod,pPlot,pRet)	\
    (This)->lpVtbl -> ARMBumpCurve(This,pZc,pEpsilon,pMethod,pPlot,pRet)

#define IARMModule_ARMAccrued(This,pSec,pDate,pModel,pRet)	\
    (This)->lpVtbl -> ARMAccrued(This,pSec,pDate,pModel,pRet)

#define IARMModule_ARMSwitchToWSETK(This)	\
    (This)->lpVtbl -> ARMSwitchToWSETK(This)

#define IARMModule_ARM_Credit_DataFromLabel(This,pricer,label,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DataFromLabel(This,pricer,label,pRet)

#define IARMModule_ARMcomputeBilibor(This,pAsOf,pStartDate,pDateSecondPhase,pEndDate,pNtl,pIndex,pIndexFund,pIndexSF,pBsmod,pBsmodFund,pBsmodDeltaCcy1,pBsmodDeltaFund,pBsmodDeltaCcy2,pBsmodFxCorrel,pFreqP,pFreqR,pFreqPFund,pFreqRFund,pFreqPSF,pFreqRSF,pResetTiming,pResetTimingSF,pCcy1,pCcy2,pDayCount,pDayCountSF,pSpdPF,pSpdSF,pSpdfund,pSpdfund2,pResetGap,pResetGapSF,pAmort,pRefDate,pFee,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRet)	\
    (This)->lpVtbl -> ARMcomputeBilibor(This,pAsOf,pStartDate,pDateSecondPhase,pEndDate,pNtl,pIndex,pIndexFund,pIndexSF,pBsmod,pBsmodFund,pBsmodDeltaCcy1,pBsmodDeltaFund,pBsmodDeltaCcy2,pBsmodFxCorrel,pFreqP,pFreqR,pFreqPFund,pFreqRFund,pFreqPSF,pFreqRSF,pResetTiming,pResetTimingSF,pCcy1,pCcy2,pDayCount,pDayCountSF,pSpdPF,pSpdSF,pSpdfund,pSpdfund2,pResetGap,pResetGapSF,pAmort,pRefDate,pFee,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRet)

#define IARMModule_ARM_Credit_GenerateImpliedCurve(This,pricerId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenerateImpliedCurve(This,pricerId,pRet)

#define IARMModule_ARM_Credit_GetEqStrike(This,CorrelId,IndexName,UpOrLow,pRetMatu,pRetStrikes)	\
    (This)->lpVtbl -> ARM_Credit_GetEqStrike(This,CorrelId,IndexName,UpOrLow,pRetMatu,pRetStrikes)

#define IARMModule_ARMcomputeOptilix(This,pAsOf,pStartDate,pDateSecondPhase,pEndDate,pNtl,pIndex,pIndexFund,pIndexSF,pBsmod,pBsmodFund,pBsmodDeltaCcy,pBsmodDeltaFund,pFreqP,pFreqR,pFreqPFund,pFreqRFund,pFreqPSF,pFreqRSF,pResetTiming,pResetTimingSF,pCcy,pDayCount,pDayCountSF,pSpdSF,pSpdfund,pResetGap,pResetGapSF,pAmort,pRefDate,pFee,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRet)	\
    (This)->lpVtbl -> ARMcomputeOptilix(This,pAsOf,pStartDate,pDateSecondPhase,pEndDate,pNtl,pIndex,pIndexFund,pIndexSF,pBsmod,pBsmodFund,pBsmodDeltaCcy,pBsmodDeltaFund,pFreqP,pFreqR,pFreqPFund,pFreqRFund,pFreqPSF,pFreqRSF,pResetTiming,pResetTimingSF,pCcy,pDayCount,pDayCountSF,pSpdSF,pSpdfund,pResetGap,pResetGapSF,pAmort,pRefDate,pFee,pFwdRule,pIntRule,pStubRule,pFreqAmort,pTxAmort,pAmountAmort,pRet)

#define IARMModule_ARM_Credit_DefaultIntensity(This,pricerId,Maturity,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DefaultIntensity(This,pricerId,Maturity,pRet)

#define IARMModule_ARMRiskyBond(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedemptionPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcyId,pRepo,pSsl,pRecoveryRate,pRet)	\
    (This)->lpVtbl -> ARMRiskyBond(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedemptionPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcyId,pRepo,pSsl,pRecoveryRate,pRet)

#define IARMModule_ARMRiskyBondWithCF(This,pAsOfDate,pRedemptionPrice,pPeriodicity,pDaycount,pYearTerms,pCashFlows,pSettleGap,pCpnDateFlag,pCcyId,pRepo,pSsl,pRecoveryRate,pRet)	\
    (This)->lpVtbl -> ARMRiskyBondWithCF(This,pAsOfDate,pRedemptionPrice,pPeriodicity,pDaycount,pYearTerms,pCashFlows,pSettleGap,pCpnDateFlag,pCcyId,pRepo,pSsl,pRecoveryRate,pRet)

#define IARMModule_ARM_Credit_EmptyLeg(This,pRet)	\
    (This)->lpVtbl -> ARM_Credit_EmptyLeg(This,pRet)

#define IARMModule_ARMClonedAndSetNotional(This,bLegId,bAmortId,pRet)	\
    (This)->lpVtbl -> ARMClonedAndSetNotional(This,bLegId,bAmortId,pRet)

#define IARMModule_ARM_INF_GetZcFromSummit(This,Index,Ccy,cvname,date,seasonAdj,seasonAdjMode,pRet)	\
    (This)->lpVtbl -> ARM_INF_GetZcFromSummit(This,Index,Ccy,cvname,date,seasonAdj,seasonAdjMode,pRet)

#define IARMModule_ARM_Credit_IRLEGTOCREDITLEG(This,SwapLegId,LegType,creditindexId,pricerId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_IRLEGTOCREDITLEG(This,SwapLegId,LegType,creditindexId,pricerId,pRet)

#define IARMModule_ARM_Credit_Collateral(This,pLabels,pNotionals,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Collateral(This,pLabels,pNotionals,pRet)

#define IARMModule_ARM_Credit_CDSGEN(This,FeeLegId,DefLegId,RcvFee,TradedNot,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDSGEN(This,FeeLegId,DefLegId,RcvFee,TradedNot,pRet)

#define IARMModule_ARM_Credit_NTDGEN(This,CdsId,firstnumdef,lastnumdef,CollateralId,binary,rcvfee,pRet)	\
    (This)->lpVtbl -> ARM_Credit_NTDGEN(This,CdsId,firstnumdef,lastnumdef,CollateralId,binary,rcvfee,pRet)

#define IARMModule_ARM_Credit_CDOGEN(This,CdsId,subamount,CollateralId,binary,rcvfee,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDOGEN(This,CdsId,subamount,CollateralId,binary,rcvfee,pRet)

#define IARMModule_ARM_Credit_GenLeg(This,StartDate,EndDate,FixedRate,FixedNotional,RefValNotional,RefValRate,XChangeNotional,Frequency,Basis,payTiming,intrule,stubrule,ccyid,paycalname,refdate,includematurity,adjstartdate,legtype,indexobj,creditlag,binary,name,Nxchange,baccruedOnDef,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenLeg(This,StartDate,EndDate,FixedRate,FixedNotional,RefValNotional,RefValRate,XChangeNotional,Frequency,Basis,payTiming,intrule,stubrule,ccyid,paycalname,refdate,includematurity,adjstartdate,legtype,indexobj,creditlag,binary,name,Nxchange,baccruedOnDef,pRet)

#define IARMModule_ARM_Credit_CDO_SQUARE_GEN(This,CdsId,subamount,portfolioId,binary,rcvfee,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CDO_SQUARE_GEN(This,CdsId,subamount,portfolioId,binary,rcvfee,pRet)

#define IARMModule_ARM_Credit_GetLastFixingDate(This,instId,asofDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetLastFixingDate(This,instId,asofDate,pRet)

#define IARMModule_ARM_Credit_SetPastFixing(This,instId,resetDate,fixingValue,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetPastFixing(This,instId,resetDate,fixingValue,pRet)

#define IARMModule_ARMGetFixing(This,source,index,term,ccy,date,pRet)	\
    (This)->lpVtbl -> ARMGetFixing(This,source,index,term,ccy,date,pRet)

#define IARMModule_ARM_Credit_SetRiskyProfile(This,CdsorCdoId,CouponsId,TypesId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetRiskyProfile(This,CdsorCdoId,CouponsId,TypesId,pRet)

#define IARMModule_ARMHyperCube(This,pVolCurvId,pKeys,pRet)	\
    (This)->lpVtbl -> ARMHyperCube(This,pVolCurvId,pKeys,pRet)

#define IARMModule_ARM_Credit_SetPricerForRatesComputation(This,legId,pricerId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetPricerForRatesComputation(This,legId,pricerId,pRet)

#define IARMModule_ARMCmsLeg(This,startDate,endDate,cmsTypeId,receiveOrPay,yieldDecompFreq,swapLegDayCount,resetFreq,payFreq,resetGap,intRule,ccyName,resetTiming,stubRule,adjStartDate,pRet)	\
    (This)->lpVtbl -> ARMCmsLeg(This,startDate,endDate,cmsTypeId,receiveOrPay,yieldDecompFreq,swapLegDayCount,resetFreq,payFreq,resetGap,intRule,ccyName,resetTiming,stubRule,adjStartDate,pRet)

#define IARMModule_ARM_Credit_SetMatuLabel(This,pCurveId,pMatuLabels,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetMatuLabel(This,pCurveId,pMatuLabels,pRet)

#define IARMModule_ARMcomputePentifix(This,pNtl,pStartdatePhase1,pCcy,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFloorOrCap,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pCorrManager,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pBumpBsGenMod,pBumpVolBsGenMod,pRet)	\
    (This)->lpVtbl -> ARMcomputePentifix(This,pNtl,pStartdatePhase1,pCcy,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFloorOrCap,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pCorrManager,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pBumpBsGenMod,pBumpVolBsGenMod,pRet)

#define IARMModule_ARM_ReplicConvAdjust_Create(This,Payoff_ReplicMode,Payoff_StepOrReplicPrecision,Payoff_StopMode,Payoff_StopThreshold,Sensi_ReplicMode,Sensi_StepOrReplicPrecision,Sensi_StopMode,Sensi_StopThreshold,UsedModelId,StrikeMinReplic,StrikeMaxReplic,pRet)	\
    (This)->lpVtbl -> ARM_ReplicConvAdjust_Create(This,Payoff_ReplicMode,Payoff_StepOrReplicPrecision,Payoff_StopMode,Payoff_StopThreshold,Sensi_ReplicMode,Sensi_StepOrReplicPrecision,Sensi_StopMode,Sensi_StopThreshold,UsedModelId,StrikeMinReplic,StrikeMaxReplic,pRet)

#define IARMModule_ARM_MapConvAdjust_Create(This,LiborArrearAdj,NaturalCMSAdj,PaymentLagAdj,pRet)	\
    (This)->lpVtbl -> ARM_MapConvAdjust_Create(This,LiborArrearAdj,NaturalCMSAdj,PaymentLagAdj,pRet)

#define IARMModule_ARM_Credit_SetFees(This,securityId,RefvalueId)	\
    (This)->lpVtbl -> ARM_Credit_SetFees(This,securityId,RefvalueId)

#define IARMModule_ARM_Credit_GetBounds(This,securityId,down,up)	\
    (This)->lpVtbl -> ARM_Credit_GetBounds(This,securityId,down,up)

#define IARMModule_ARMcomputePentibor(This,pNtl,pStartdatePhase1,pCcy,pIndexPay,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pIndexIndexCorrelCube,pCorrEUR,pInterCorr,pRet)	\
    (This)->lpVtbl -> ARMcomputePentibor(This,pNtl,pStartdatePhase1,pCcy,pIndexPay,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pIndexIndexCorrelCube,pCorrEUR,pInterCorr,pRet)

#define IARMModule_ARM_Credit_SetRecovCoef(This,pCurveId,RecovCoef,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetRecovCoef(This,pCurveId,RecovCoef,pRet)

#define IARMModule_ARMIndexIndexCorrelCube(This,pVolCurvId,pTenors1List,pTenors2List,pInterSurfInterp,pRet)	\
    (This)->lpVtbl -> ARMIndexIndexCorrelCube(This,pVolCurvId,pTenors1List,pTenors2List,pInterSurfInterp,pRet)

#define IARMModule_ARM_Credit_SetInterpolationType(This,pVolCurveId,pInterpolType,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetInterpolationType(This,pVolCurveId,pInterpolType,pRet)

#define IARMModule_ARM_Credit_Customized_CDO(This,pLabels,pNotionals,Currency,pDefaultLeg,pPremiumLeg,pParameters,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Customized_CDO(This,pLabels,pNotionals,Currency,pDefaultLeg,pPremiumLeg,pParameters,pRet)

#define IARMModule_ARMCreateGenCorrelatorManager(This,pMktTags,pHyperDiagVol,pIndexIndexVol,pCorrelVol,pIndexVol,pRet)	\
    (This)->lpVtbl -> ARMCreateGenCorrelatorManager(This,pMktTags,pHyperDiagVol,pIndexIndexVol,pCorrelVol,pIndexVol,pRet)

#define IARMModule_ARM_Credit_CLN(This,EffectiveDate,EndDate,Spread,IndexId,Refdate,pFstCpnEffDate,Notional,AccOnDef,DayCount,DecompFreq,StubRule,resetgap,Currency,ResetCal,PayCal,Nxchange,IncludeMaturity,AdjustedStartDate,Binary,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CLN(This,EffectiveDate,EndDate,Spread,IndexId,Refdate,pFstCpnEffDate,Notional,AccOnDef,DayCount,DecompFreq,StubRule,resetgap,Currency,ResetCal,PayCal,Nxchange,IncludeMaturity,AdjustedStartDate,Binary,pRet)

#define IARMModule_ARMcomputePentilix(This,pNtl,pStartdatePhase1,pCcy,pIndexPay,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pIndexIndexCorrelCube,pCorrEUR,pInterCorr,pRet)	\
    (This)->lpVtbl -> ARMcomputePentilix(This,pNtl,pStartdatePhase1,pCcy,pIndexPay,pIndexPhase1,pSpreadPhase1,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pStartDatePhase2,pSpreadPhase2,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pStartDatePhase3,pEndDatePhase3,pIndexPhase3,pSpreadPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pVolCurvFromMatriceShift,pVol,pVolCub,pConvexityManager,pZc,pSmiledMod,pSmiledModBump,pHyperCubeCorrel,pIndexIndexCorrelCube,pCorrEUR,pInterCorr,pRet)

#define IARMModule_ARM_Credit_RiskyPV01(This,pDefCurve,Date1,Date2,pRet)	\
    (This)->lpVtbl -> ARM_Credit_RiskyPV01(This,pDefCurve,Date1,Date2,pRet)

#define IARMModule_ARMGetResetMgrFromSummit(This,pAsOf,pIndex,pSource,pCcy,pIsInflationIndex,pTerm,pRet)	\
    (This)->lpVtbl -> ARMGetResetMgrFromSummit(This,pAsOf,pIndex,pSource,pCcy,pIsInflationIndex,pTerm,pRet)

#define IARMModule_ARMGetReset(This,pResetMgr,pDate,pRet)	\
    (This)->lpVtbl -> ARMGetReset(This,pResetMgr,pDate,pRet)

#define IARMModule_ARMSetLastFixing(This,pSecurityId,pRate,pAsOf,pBeforeLastFixingDate,pResetDate,pRet)	\
    (This)->lpVtbl -> ARMSetLastFixing(This,pSecurityId,pRate,pAsOf,pBeforeLastFixingDate,pResetDate,pRet)

#define IARMModule_ARM_Credit_Sectorial_Correlation(This,AsOf,structName,correlation_Type,vLabels,vector_Membership,intra_Sector_Correlation,inter_Sector_Correlation,vBetas,vLambdas,vBetas_Down,vLambdas_Down,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Sectorial_Correlation(This,AsOf,structName,correlation_Type,vLabels,vector_Membership,intra_Sector_Correlation,inter_Sector_Correlation,vBetas,vLambdas,vBetas_Down,vLambdas_Down,pRet)

#define IARMModule_ARMcomputeReviPentix(This,pNtl,pDate,pCcy,pIndexPhase1,pSpread,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pIndexPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFloorOrCap,pFee,pLevier,pTxFixeMax,pIsCapped,pTxCap,pVolCurvFromMatriceShift,pVol,pVolCub,pCorrManager,pConvexityManager,pZc,pSmiledMod,pHyperCubeCorrel,pBumpBsGenMod,pBumpVolBsGenMod,pRet)	\
    (This)->lpVtbl -> ARMcomputeReviPentix(This,pNtl,pDate,pCcy,pIndexPhase1,pSpread,pDayCountPhase1,pPayFreqPhase1,pResetFreqPhase1,pResetTimingPhase1,pRoll,pAdjPhase1,pStub,pIndexPhase2DIG,pIndexLongPhase2DIG,pStrikePhase2DIG,pResetTimingPhase2DIG,pAdjPhase2DIG,pDayCountPhase2,pPayFreqPhase2,pResetFreqPhase2,pAdjPhase2,pIndexPhase3,pDayCountPhase3,pPayFreqPhase3,pResetFreqPhase3,pResetTimingPhase3,pAdjPhase3,pIndexFund,pSpreadFund,pDayCountFund,pPayFreqFund,pResetFreqFund,pResetTimingFund,pAdjFund,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFloorOrCap,pFee,pLevier,pTxFixeMax,pIsCapped,pTxCap,pVolCurvFromMatriceShift,pVol,pVolCub,pCorrManager,pConvexityManager,pZc,pSmiledMod,pHyperCubeCorrel,pBumpBsGenMod,pBumpVolBsGenMod,pRet)

#define IARMModule_ARMGenAmortization(This,pSwaplegId,pAmortMethod,pAmortFreq,pAmortAmount,pDaycount,pLegNotional,pAmortRate,pReducedMaturity,pModelId,pCleanUp,pRet)	\
    (This)->lpVtbl -> ARMGenAmortization(This,pSwaplegId,pAmortMethod,pAmortFreq,pAmortAmount,pDaycount,pLegNotional,pAmortRate,pReducedMaturity,pModelId,pCleanUp,pRet)

#define IARMModule_ARMCptRefvalue(This,pRefValId,pDate,pRet)	\
    (This)->lpVtbl -> ARMCptRefvalue(This,pRefValId,pDate,pRet)

#define IARMModule_ARM_Credit_GetDPFromCalypso(This,pDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,irCurveId,label,ret)	\
    (This)->lpVtbl -> ARM_Credit_GetDPFromCalypso(This,pDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,irCurveId,label,ret)

#define IARMModule_ARMCalypsoDevConnect(This)	\
    (This)->lpVtbl -> ARMCalypsoDevConnect(This)

#define IARMModule_ARMGetZCFromCalypso(This,pIndex,pCurrency,pTerm,pricingEnv,pDate,pInterpMethod,forceCurveName,xmlFile,pRet)	\
    (This)->lpVtbl -> ARMGetZCFromCalypso(This,pIndex,pCurrency,pTerm,pricingEnv,pDate,pInterpMethod,forceCurveName,xmlFile,pRet)

#define IARMModule_ARM_Credit_DefProbInverse(This,pCurveId,dDefProba,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DefProbInverse(This,pCurveId,dDefProba,pRet)

#define IARMModule_ARM_Credit_DPMktDataFromCalypso(This,AsOfDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,Parameter,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DPMktDataFromCalypso(This,AsOfDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,Parameter,pRet)

#define IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromCalypso(This,pDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,irCurveId,label,ret)	\
    (This)->lpVtbl -> ARM_Credit_ZeroCouponDefaultCurveFromCalypso(This,pDate,pricingEnv,issuer,seniority,ccy,forceCurveName,xmlFile,irCurveId,label,ret)

#define IARMModule_ARMGetInitialCurveFromCalypso(This,pIndex,pCurrency,pTerm,pricingEnv,pDate,forceCurveName,xmlFile,pDoAdj,pRetMat,pRetRate)	\
    (This)->lpVtbl -> ARMGetInitialCurveFromCalypso(This,pIndex,pCurrency,pTerm,pricingEnv,pDate,forceCurveName,xmlFile,pDoAdj,pRetMat,pRetRate)

#define IARMModule_ARMCalypsoProdConnect(This)	\
    (This)->lpVtbl -> ARMCalypsoProdConnect(This)

#define IARMModule_ARMCalypsoRecConnect(This)	\
    (This)->lpVtbl -> ARMCalypsoRecConnect(This)

#define IARMModule_ARM_Credit_GetBasketCorrelMkDataFromCalypso(This,pricingEnv,date,forceCurveName,xmlFileName,pRetMat,pRetTenor,pRetVol)	\
    (This)->lpVtbl -> ARM_Credit_GetBasketCorrelMkDataFromCalypso(This,pricingEnv,date,forceCurveName,xmlFileName,pRetMat,pRetTenor,pRetVol)

#define IARMModule_ARM_Credit_QMatrix(This,pQMatrix,pRet)	\
    (This)->lpVtbl -> ARM_Credit_QMatrix(This,pQMatrix,pRet)

#define IARMModule_ARM_Credit_MarketDataMng(This,pstrVect,pRet)	\
    (This)->lpVtbl -> ARM_Credit_MarketDataMng(This,pstrVect,pRet)

#define IARMModule_ARM_Credit_ModelMultiCvMktDataMng(This,pIRcurve,pDefCurves,pRecovery,CorrelId,MktdataMngId,pVolcurve,cloneorNot,pRet)	\
    (This)->lpVtbl -> ARM_Credit_ModelMultiCvMktDataMng(This,pIRcurve,pDefCurves,pRecovery,CorrelId,MktdataMngId,pVolcurve,cloneorNot,pRet)

#define IARMModule_Local_ARM_ProdConnect(This)	\
    (This)->lpVtbl -> Local_ARM_ProdConnect(This)

#define IARMModule_ARM_SetDefaultCurrency(This,isoCCy,pRet)	\
    (This)->lpVtbl -> ARM_SetDefaultCurrency(This,isoCCy,pRet)

#define IARMModule_ARMLivretALeg(This,pStartDate,pEndDate,pRcvOrPay,pSpread,pResetFreq,pPayFreq,pResetTiming,pPayTiming,pCcy,pIntRule,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pDayCount,pRet)	\
    (This)->lpVtbl -> ARMLivretALeg(This,pStartDate,pEndDate,pRcvOrPay,pSpread,pResetFreq,pPayFreq,pResetTiming,pPayTiming,pCcy,pIntRule,pResetGap,pResetCal,pPayCal,pDecompPricingFlag,pNxChange,pStubRule,pRefDate,pAdjStartDate,pDayCount,pRet)

#define IARMModule_ARM_Credit_CorrelationSmileStrike(This,pLabels,pVolCurves,pProportions,AsOf,pSmileStrikeLow,pSmileStrikeHigh,pIndexVector,Name,pFullStrikeLow,pFullStrikeUp,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CorrelationSmileStrike(This,pLabels,pVolCurves,pProportions,AsOf,pSmileStrikeLow,pSmileStrikeHigh,pIndexVector,Name,pFullStrikeLow,pFullStrikeUp,pRet)

#define IARMModule_ARMFutDelivery(This,pFut,pCcy,pRet)	\
    (This)->lpVtbl -> ARMFutDelivery(This,pFut,pCcy,pRet)

#define IARMModule_ARMLivretACurve(This,pAsOf,pInfCurvId,pEuribCurvId,pFlagRouding,pInfResetMgrId,pFixingLivretAId,pFixingEuribId,pMonthForAugust,pMonthForFebruary,pRet)	\
    (This)->lpVtbl -> ARMLivretACurve(This,pAsOf,pInfCurvId,pEuribCurvId,pFlagRouding,pInfResetMgrId,pFixingLivretAId,pFixingEuribId,pMonthForAugust,pMonthForFebruary,pRet)

#define IARMModule_ARMcomputeLivretA(This,pNtl,pStartDateLeg1,pEndDateLeg1,pCcy,pIndexLeg1,pSpreadLeg1,pDayCountLeg1,pPayFreqLeg1,pResetFreqLeg1,pResetTimingLeg1,pAdjLeg1,pRoll,pStub,pEndDateLA,pSpreadLeg2,pDayCountLA,pPayFreqLA,pResetFreqLA,pResetTimingLA,pAdjLA,pIndexLeg2,pDayCountLeg2,pPayFreqLeg2,pResetFreqLeg2,pResetTimingLeg2,pAdjLeg2,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pSmiledMod,pSmiledModBump,pLAMod,pLAModBump,pLAModBumpInflation,pResetMgrIds,pRet)	\
    (This)->lpVtbl -> ARMcomputeLivretA(This,pNtl,pStartDateLeg1,pEndDateLeg1,pCcy,pIndexLeg1,pSpreadLeg1,pDayCountLeg1,pPayFreqLeg1,pResetFreqLeg1,pResetTimingLeg1,pAdjLeg1,pRoll,pStub,pEndDateLA,pSpreadLeg2,pDayCountLA,pPayFreqLA,pResetFreqLA,pResetTimingLA,pAdjLA,pIndexLeg2,pDayCountLeg2,pPayFreqLeg2,pResetFreqLeg2,pResetTimingLeg2,pAdjLeg2,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pSmiledMod,pSmiledModBump,pLAMod,pLAModBump,pLAModBumpInflation,pResetMgrIds,pRet)

#define IARMModule_ARMGetFixingFromCalypso(This,source,index,term,ccy,curveName,date,pRet)	\
    (This)->lpVtbl -> ARMGetFixingFromCalypso(This,source,index,term,ccy,curveName,date,pRet)

#define IARMModule_ARM_FxConvertFromCalypso(This,ccy1,ccy2,pCvName,pDate,pRet)	\
    (This)->lpVtbl -> ARM_FxConvertFromCalypso(This,ccy1,ccy2,pCvName,pDate,pRet)

#define IARMModule_ARM_Credit_FwdSpread(This,defcurveId,Maturity1,Maturity2,FwdStartDate,FwdEndDate,VolId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_FwdSpread(This,defcurveId,Maturity1,Maturity2,FwdStartDate,FwdEndDate,VolId,pRet)

#define IARMModule_ARM_Credit_Flat_Correlation(This,AsOf,structName,correlValue,idIndex1,idIndex2,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Flat_Correlation(This,AsOf,structName,correlValue,idIndex1,idIndex2,pRet)

#define IARMModule_ARM_Credit_CorridorLeg_Sche(This,Notional,RecieveOrPay,RefValueSpreadsInBP,floatingIdx,leverageFloatIdx,creditIdx,refvalueKUPinBP,refvalueKDWinBP,ScheduleInfoId,accondef,disc_ccy,Name,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CorridorLeg_Sche(This,Notional,RecieveOrPay,RefValueSpreadsInBP,floatingIdx,leverageFloatIdx,creditIdx,refvalueKUPinBP,refvalueKDWinBP,ScheduleInfoId,accondef,disc_ccy,Name,pRet)

#define IARMModule_ARM_Credit_Schedule_Info(This,EffectiveDate,MaturityDate,payFrequency,ResetFreq,DayCount,Stubrule,intRule,payCalName,PayTiming,ResetTiming,fwdRule,IncludeMaturity,adj,intStartAdj,AccDayCount,ReferenceDate,FirstCpnEffDate,AdjCal,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Schedule_Info(This,EffectiveDate,MaturityDate,payFrequency,ResetFreq,DayCount,Stubrule,intRule,payCalName,PayTiming,ResetTiming,fwdRule,IncludeMaturity,adj,intStartAdj,AccDayCount,ReferenceDate,FirstCpnEffDate,AdjCal,pRet)

#define IARMModule_ARMInfCurveSetResetMgr(This,pInfCurve,pResetMgr,pRet)	\
    (This)->lpVtbl -> ARMInfCurveSetResetMgr(This,pInfCurve,pResetMgr,pRet)

#define IARMModule_ARM_Credit_CPDO(This,pRiskyLeg,pRollLeg,pNoRiskyLeg,pInitialValo,pTarget,pMaturity,pCpnType,pUFFees,pRunningFees,pVExpo,pV0Expo,pAlpha,pBeta,pDesactivation,pNbAssets,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CPDO(This,pRiskyLeg,pRollLeg,pNoRiskyLeg,pInitialValo,pTarget,pMaturity,pCpnType,pUFFees,pRunningFees,pVExpo,pV0Expo,pAlpha,pBeta,pDesactivation,pNbAssets,pRet)

#define IARMModule_ARM_Credit_PriceVector(This,pPricer,pCPTTYPE,pRetVectorValos)	\
    (This)->lpVtbl -> ARM_Credit_PriceVector(This,pPricer,pCPTTYPE,pRetVectorValos)

#define IARMModule_ARM_Credit_GenPrice(This,pPricer,pCPTTYPE,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenPrice(This,pPricer,pCPTTYPE,pRet)

#define IARMModule_ARMcomputeTxFixed(This,pNtl,pStartDateLeg1,pEndDateLeg1,pCcy,pIndexLeg1,pSpreadLeg1,pDayCountLeg1,pPayFreqLeg1,pResetFreqLeg1,pResetTimingLeg1,pAdjLeg1,pRoll,pStub,pEndDateFixed,pSpreadLeg2,pDayCountFixed,pPayFreqFixed,pResetFreqFixed,pResetTimingFixed,pAdjFixed,pIndexLeg2,pDayCountLeg2,pPayFreqLeg2,pResetFreqLeg2,pResetTimingLeg2,pAdjLeg2,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pSmiledMod,pSmiledModBump,pRet)	\
    (This)->lpVtbl -> ARMcomputeTxFixed(This,pNtl,pStartDateLeg1,pEndDateLeg1,pCcy,pIndexLeg1,pSpreadLeg1,pDayCountLeg1,pPayFreqLeg1,pResetFreqLeg1,pResetTimingLeg1,pAdjLeg1,pRoll,pStub,pEndDateFixed,pSpreadLeg2,pDayCountFixed,pPayFreqFixed,pResetFreqFixed,pResetTimingFixed,pAdjFixed,pIndexLeg2,pDayCountLeg2,pPayFreqLeg2,pResetFreqLeg2,pResetTimingLeg2,pAdjLeg2,pEndDateAmort,pDayCountAmort,pIntRuleAmort,pTxAmort,pFreqAmort,pAmountAmort,pTypeAmort,pFee,pSmiledMod,pSmiledModBump,pRet)

#define IARMModule_ARM_Credit_FixingCurve(This,pDates,pValues,AsOfDate,B_IndexName,B_IndexID,pRet)	\
    (This)->lpVtbl -> ARM_Credit_FixingCurve(This,pDates,pValues,AsOfDate,B_IndexName,B_IndexID,pRet)

#define IARMModule_ARM_Credit_CptInterpolDefCurveOLD(This,pCurve,pTenor,pSlope,pDate,pInterpDate,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CptInterpolDefCurveOLD(This,pCurve,pTenor,pSlope,pDate,pInterpDate,pRet)

#define IARMModule_ARM_Credit_CreateBasketCorrelMkDataFromCalypso(This,pricingEnv,date,forceCurveName,Ccy,xmlFileName,indexId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_CreateBasketCorrelMkDataFromCalypso(This,pricingEnv,date,forceCurveName,Ccy,xmlFileName,indexId,pRet)

#define IARMModule_ARMSetCalendar(This,pFileName,pRet)	\
    (This)->lpVtbl -> ARMSetCalendar(This,pFileName,pRet)

#define IARMModule_ARMInitGigaSpaces(This,pUrl,pRet)	\
    (This)->lpVtbl -> ARMInitGigaSpaces(This,pUrl,pRet)

#define IARMModule_ARM_Credit_GetExpectedLoss(This,pricerId,YearTerm,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GetExpectedLoss(This,pricerId,YearTerm,pRet)

#define IARMModule_ARM_Credit_VariableCollateral(This,pLabels,pNotionals,pRet)	\
    (This)->lpVtbl -> ARM_Credit_VariableCollateral(This,pLabels,pNotionals,pRet)

#define IARMModule_ARM_Credit_IndexCompo(This,IndexName,pLabels,YearFrac,pSpread,Method,Basis,ResetFreq,PayFreq,ccy,fwdRule,resetTiming,resetGap,payTiming,payGap,intRule,AdjCalType,cm_resetWeekDay,cm_resetOccur,pRet)	\
    (This)->lpVtbl -> ARM_Credit_IndexCompo(This,IndexName,pLabels,YearFrac,pSpread,Method,Basis,ResetFreq,PayFreq,ccy,fwdRule,resetTiming,resetGap,payTiming,payGap,intRule,AdjCalType,cm_resetWeekDay,cm_resetOccur,pRet)

#define IARMModule_ARM_Credit_createFlatCurve(This,pCurve,pTenor,pRet)	\
    (This)->lpVtbl -> ARM_Credit_createFlatCurve(This,pCurve,pTenor,pRet)

#define IARMModule_ARM_Credit_FunctionRegister(This,address)	\
    (This)->lpVtbl -> ARM_Credit_FunctionRegister(This,address)

#define IARMModule_ARMBermudanXStyle(This,pxDates,pexpiryDates,pRet)	\
    (This)->lpVtbl -> ARMBermudanXStyle(This,pxDates,pexpiryDates,pRet)

#define IARMModule_ARMcomputeCRA(This,pFixorFloat,pFee,pAsOf,pStartDate,pEndDate,pCcy,pLevelUp,pUpSpec,pLevelDown,pDownSpec,pRefIndex,pDayCount,pPayFreqPayIndex,pResetFreqRefIndex,pPaidRstTiming,pRefRstTiming,pStubRule,pPOrR,pStartCallDate,pXStyle,pFundingIndex,pResetFreqFunding,pPayFreqFunding,pSpreadFunding,pPOrRFunding,pDecompPricingFlag,pdiscMarginFactor,pPreInitFlag,pMeanReversion,pCalibParams,pCalibParamsPF,pKernelToGP,pMarkovTreeParams,pMarkovTreePathNumber,pBsmodId,pBsmodSwoptId,pBsmodSwoptBumpId,pzcId,pRet)	\
    (This)->lpVtbl -> ARMcomputeCRA(This,pFixorFloat,pFee,pAsOf,pStartDate,pEndDate,pCcy,pLevelUp,pUpSpec,pLevelDown,pDownSpec,pRefIndex,pDayCount,pPayFreqPayIndex,pResetFreqRefIndex,pPaidRstTiming,pRefRstTiming,pStubRule,pPOrR,pStartCallDate,pXStyle,pFundingIndex,pResetFreqFunding,pPayFreqFunding,pSpreadFunding,pPOrRFunding,pDecompPricingFlag,pdiscMarginFactor,pPreInitFlag,pMeanReversion,pCalibParams,pCalibParamsPF,pKernelToGP,pMarkovTreeParams,pMarkovTreePathNumber,pBsmodId,pBsmodSwoptId,pBsmodSwoptBumpId,pzcId,pRet)

#define IARMModule_ARM_Credit_Math_BivNormale(This,x,y,rho,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Math_BivNormale(This,x,y,rho,pRet)

#define IARMModule_ARM_GETINSTRUMENTFROMCALYPSO(This,CalypsoId,Type,AsOf,ModelType,pRet)	\
    (This)->lpVtbl -> ARM_GETINSTRUMENTFROMCALYPSO(This,CalypsoId,Type,AsOf,ModelType,pRet)

#define IARMModule_ARMSetDiscountPricingMode(This,pModelId,pDiscountPricingMode,pRet)	\
    (This)->lpVtbl -> ARMSetDiscountPricingMode(This,pModelId,pDiscountPricingMode,pRet)

#define IARMModule_ARM_Credit_Math_RandUniform(This,seed,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Math_RandUniform(This,seed,pRet)

#define IARMModule_ARM_Credit_Math_Interpol(This,X,Y,value,type,smooth,Weights,modeSpline,withC1condition,leftSlope,rightSlope,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Math_Interpol(This,X,Y,value,type,smooth,Weights,modeSpline,withC1condition,leftSlope,rightSlope,pRet)

#define IARMModule_ARM_Credit_Random_Generator(This,RandomType,ParamId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Random_Generator(This,RandomType,ParamId,pRet)

#define IARMModule_ARM_Credit_GenerateOneRandom(This,RandomId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenerateOneRandom(This,RandomId,pRet)

#define IARMModule_ARM_Credit_GenerateRandoms(This,RandomId,DimVector,pRet)	\
    (This)->lpVtbl -> ARM_Credit_GenerateRandoms(This,RandomId,DimVector,pRet)

#define IARMModule_ARM_Credit_ResetRandom(This,RandomId)	\
    (This)->lpVtbl -> ARM_Credit_ResetRandom(This,RandomId)

#define IARMModule_ARM_Credit_FwdSpreadAsIndex(This,DefCurveId,matu1,Matu2,pRet)	\
    (This)->lpVtbl -> ARM_Credit_FwdSpreadAsIndex(This,DefCurveId,matu1,Matu2,pRet)

#define IARMModule_ARM_Credit_createDefCurveFromBase(This,pCurveCDS,pCurveIndex,vBase,pRet)	\
    (This)->lpVtbl -> ARM_Credit_createDefCurveFromBase(This,pCurveCDS,pCurveIndex,vBase,pRet)

#define IARMModule_ARM_Credit_RiskyPV01AsSensitivity(This,pDefCurve,Tenor,pRet)	\
    (This)->lpVtbl -> ARM_Credit_RiskyPV01AsSensitivity(This,pDefCurve,Tenor,pRet)

#define IARMModule_ARM_Credit_SetVolCurve(This,Model,VolCurveId,pRet)	\
    (This)->lpVtbl -> ARM_Credit_SetVolCurve(This,Model,VolCurveId,pRet)

#define IARMModule_ARMPF(This,pinsts,pcoeffs,pmarketPrices,pprecisions,pRet)	\
    (This)->lpVtbl -> ARMPF(This,pinsts,pcoeffs,pmarketPrices,pprecisions,pRet)

#define IARMModule_ARMBondTEC(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedempPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcyId,ptec,pPFTecId,pModTecId,pRet)	\
    (This)->lpVtbl -> ARMBondTEC(This,pIssueDate,pMaturityDate,pFirstCpnDate,pCpnRate,pRedempPrice,pPeriodicity,pDaycount,pSettleGap,pCpnDateFlag,pCcyId,ptec,pPFTecId,pModTecId,pRet)

#define IARMModule_ARMPFModFit(This,pmodName,ppf,psettlement,pzc,pvList,pfList,nag_algo,pstep,phorizon,pRet)	\
    (This)->lpVtbl -> ARMPFModFit(This,pmodName,ppf,psettlement,pzc,pvList,pfList,nag_algo,pstep,phorizon,pRet)

#define IARMModule_ARM_Credit_Restrikable_CDO(This,UnderlyingMatu,Expiry,Strike,OptionType,pUnderlying,Rehauss,TriggerFreq,DiffCDO,pRet)	\
    (This)->lpVtbl -> ARM_Credit_Restrikable_CDO(This,UnderlyingMatu,Expiry,Strike,OptionType,pUnderlying,Rehauss,TriggerFreq,DiffCDO,pRet)

#define IARMModule_ARM_Credit_PropertyList(This,attrNames,attrValues,attrTypes,pRet)	\
    (This)->lpVtbl -> ARM_Credit_PropertyList(This,attrNames,attrValues,attrTypes,pRet)

#define IARMModule_ARM_Credit_DefCurveIntensityPWC(This,AsOfDate,pMatuRates,pInputs,Type,Recovery,IRCurveId,bCurrency,bLabel,VolCurveId,calibrationAlgo,lag,pRet)	\
    (This)->lpVtbl -> ARM_Credit_DefCurveIntensityPWC(This,AsOfDate,pMatuRates,pInputs,Type,Recovery,IRCurveId,bCurrency,bLabel,VolCurveId,calibrationAlgo,lag,pRet)

#define IARMModule_ARMTMLeg(This,ptmIxType,pstartDate,pendDate,pPorR,pspread,ppayFrequency,presetFrequency,pinterestRule,pfwdRule,pstubRule,pccy,pRet)	\
    (This)->lpVtbl -> ARMTMLeg(This,ptmIxType,pstartDate,pendDate,pPorR,pspread,ppayFrequency,presetFrequency,pinterestRule,pfwdRule,pstubRule,pccy,pRet)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMNextBusinessDay_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCalendrier,
    /* [in] */ long pNbDays,
    /* [retval][out] */ double __RPC_FAR *pDate2);


void __RPC_STUB IARMModule_ARMNextBusinessDay_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMAdjustToBusDate_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCalendrier,
    /* [in] */ BSTR pRule,
    /* [retval][out] */ double __RPC_FAR *pDate2);


void __RPC_STUB IARMModule_ARMAdjustToBusDate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFreeObject_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pId,
    /* [retval][out] */ long __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFreeObject_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMIsBusinessDay_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCalendrier,
    /* [retval][out] */ long __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMIsBusinessDay_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetZCFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pCvName,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetZCFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFreeAllObjects_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [retval][out] */ long __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFreeAllObjects_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMYcMod_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZc,
    /* [defaultvalue][in][optional] */ BSTR pZcDiscount,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMYcMod_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMForwardYield_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZc,
    /* [in] */ double pMatu1,
    /* [in] */ double pMatu2,
    /* [in] */ BSTR pMeth,
    /* [in] */ BSTR pAdj,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMForwardYield_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMDiscountYield_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pZc,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [in] */ VARIANT __RPC_FAR *pMeth,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMDiscountYield_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLiborSwap_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pStartDate,
    /* [in] */ VARIANT __RPC_FAR *pEndDate,
    /* [in] */ VARIANT __RPC_FAR *pLiborType,
    /* [in] */ VARIANT __RPC_FAR *pRecOrPay,
    /* [in] */ VARIANT __RPC_FAR *pFixedRate,
    /* [in] */ VARIANT __RPC_FAR *pSpread,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [defaultvalue][in][optional] */ BSTR pDaycount,
    /* [defaultvalue][in][optional] */ BSTR pFloatingDaycount,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLiborSwap_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwapPriceToRate_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pSwap,
    /* [in] */ VARIANT __RPC_FAR *pDate,
    /* [in] */ VARIANT __RPC_FAR *pPrice,
    /* [in] */ VARIANT __RPC_FAR *pModel,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSwapPriceToRate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pSec,
    /* [in] */ VARIANT __RPC_FAR *pModel,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBetweenDates_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pDate1,
    /* [in] */ VARIANT __RPC_FAR *pDate2,
    /* [in] */ VARIANT __RPC_FAR *pDaycount,
    /* [in] */ VARIANT __RPC_FAR *pIsYearFrac,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBetweenDates_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMAddPeriod_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pDate,
    /* [in] */ VARIANT __RPC_FAR *pFreq,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [in] */ VARIANT __RPC_FAR *pNbPeriods,
    /* [in] */ VARIANT __RPC_FAR *pAdjRule,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMAddPeriod_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMIsoCcy_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCcy,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMIsoCcy_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetSpotDays_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetSpotDays_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetLiborIndexDaycount_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetLiborIndexDaycount_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetLiborTerm_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetLiborTerm_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFixedDayCount_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFixedDayCount_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFixedPayFreq_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFixedPayFreq_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMComputeVolatility_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pVol,
    /* [in] */ VARIANT __RPC_FAR *pmatu,
    /* [in] */ VARIANT __RPC_FAR *pStrike,
    /* [in] */ VARIANT __RPC_FAR *pTenor,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMComputeVolatility_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMVolCurv_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [in] */ VARIANT __RPC_FAR *pStrike,
    /* [in] */ VARIANT __RPC_FAR *pVol,
    /* [in] */ double pAsOf,
    /* [in] */ BSTR pStrikeType,
    /* [in] */ BSTR pVolType,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR pIndexId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMVolCurv_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetVolCubeFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pCvName,
    /* [in] */ double pDate,
    /* [in] */ BSTR pType,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pSmiles,
    /* [defaultvalue][in][optional] */ BSTR pTypeCube,
    /* [defaultvalue][in][optional] */ BSTR indexId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetVolCubeFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_FxConvert_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pccy1,
    /* [in] */ VARIANT __RPC_FAR *pccy2,
    /* [in] */ VARIANT __RPC_FAR *pDate,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pCvName,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_FxConvert_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_DiscountPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pcurve,
    /* [in] */ VARIANT __RPC_FAR *pmatu,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_DiscountPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Delivery_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pAsOfDate,
    /* [in] */ VARIANT __RPC_FAR *pTenorContract,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Delivery_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CptInterpolDefCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurve,
    /* [in] */ VARIANT pTenor,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CptInterpolDefCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DefaultProba_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCurve,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DefaultProba_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetBeta_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pLabel,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetBeta_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Mezzanine_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pEffectiveDate,
    /* [in] */ double pEndDate,
    /* [in] */ double pSpread,
    /* [in] */ double pMezzAmount,
    /* [in] */ double pSubAmount,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pNotionals,
    /* [defaultvalue][in] */ BSTR pFreqFeeLeg,
    /* [defaultvalue][in] */ BSTR pDayCount,
    /* [defaultvalue][in] */ double pFirst_period_refdate,
    /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
    /* [defaultvalue][in] */ BSTR Currency,
    /* [defaultvalue][in] */ double pPayCreditLag,
    /* [defaultvalue][in] */ BSTR pStub,
    /* [defaultvalue][in] */ BSTR pFreqDefLeg,
    /* [defaultvalue][in] */ double pBinary,
    /* [defaultvalue][in] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
    /* [defaultvalue][in][optional] */ double TradedNotional,
    /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
    /* [defaultvalue][in] */ double pFstCpnEffDate,
    /* [defaultvalue][in][optional] */ BSTR intRule,
    /* [defaultvalue][in][optional] */ BSTR adjStartDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Mezzanine_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ModelMultiCurves_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pIRcurve,
    /* [in] */ VARIANT __RPC_FAR *pDefCurves,
    /* [in] */ VARIANT __RPC_FAR *pRecoveryRates,
    /* [defaultvalue][in][optional] */ BSTR CorrelId,
    /* [defaultvalue][in][optional] */ BSTR VolCurve,
    /* [defaultvalue][in][optional] */ BSTR CpnInfCurve,
    /* [defaultvalue][in][optional] */ BSTR CpnIRCurve,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ModelMultiCurves_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FTD_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
    /* [in] */ VARIANT __RPC_FAR *pEndDate,
    /* [in] */ VARIANT __RPC_FAR *pSpread,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
    /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
    /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
    /* [in] */ VARIANT __RPC_FAR *pIssuerNotional,
    /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
    /* [in] */ VARIANT __RPC_FAR *pCurrency,
    /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
    /* [in] */ VARIANT __RPC_FAR *pStub,
    /* [defaultvalue][in] */ double pFstCpnEffDate,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pintRule,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pstartAdj,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_FTD_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_NTD_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pEffectiveDate,
    /* [in] */ double pEndDate,
    /* [in] */ double pSpread,
    /* [defaultvalue][in] */ int pFirstNumDefault,
    /* [defaultvalue][in] */ int pLastNumDefault,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [defaultvalue][in] */ BSTR pFixingFreq,
    /* [defaultvalue][in] */ BSTR pDayCountFrq,
    /* [defaultvalue][in] */ double pFirst_period_refdate,
    /* [defaultvalue][in] */ double pIssuerNotional,
    /* [defaultvalue][in] */ BSTR pAccruedOnDefault,
    /* [defaultvalue][in] */ BSTR pCurrency,
    /* [defaultvalue][in] */ double pPayCreditLag,
    /* [defaultvalue][in] */ BSTR pStub,
    /* [defaultvalue][in] */ BSTR pFreqDefLeg,
    /* [defaultvalue][in] */ double pBinary,
    /* [defaultvalue][in] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
    /* [defaultvalue][in][optional] */ double TradedNotional,
    /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
    /* [defaultvalue][in] */ double pFstCpnEffDate,
    /* [defaultvalue][in][optional] */ BSTR intRule,
    /* [defaultvalue][in][optional] */ BSTR startAdj,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_NTD_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Price_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pAsofDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Price_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_RiskyDuration_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDefCurve,
    /* [in] */ VARIANT date,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_RiskyDuration_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDONPV_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pCPTTYPE,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDONPV_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CorrMatrix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pCoefs,
    /* [defaultvalue][in] */ double AsOf,
    /* [defaultvalue][in] */ BSTR Name,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CorrMatrix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ExtractCorrMatrix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCorrMatrixId,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ExtractCorrMatrix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Spread_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pMTM,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Spread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetLabel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCurveId,
    /* [in] */ VARIANT __RPC_FAR *pLabel,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetLabel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetLabel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCurveId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetLabel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Sensitivity_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricerId,
    /* [in] */ VARIANT __RPC_FAR *pType,
    /* [in] */ VARIANT __RPC_FAR *pPlot,
    /* [in] */ VARIANT __RPC_FAR *pLabel,
    /* [in] */ VARIANT __RPC_FAR *pEpsilon,
    /* [defaultvalue][in][optional] */ double epsilonGamma,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Sensitivity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetCleanSpreadTranche_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricerId,
    /* [in] */ VARIANT __RPC_FAR *pPlot,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetCleanSpreadTranche_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetDefProbTranche_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR PricerId,
    /* [in] */ double Yearterm,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetDefProbTranche_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetDuration_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricerId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetDuration_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenSchedule_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pAccStartDate,
    /* [in] */ VARIANT __RPC_FAR *pAccEndDate,
    /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
    /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
    /* [in] */ VARIANT __RPC_FAR *prefDate,
    /* [in] */ VARIANT __RPC_FAR *pCurrency,
    /* [in] */ VARIANT __RPC_FAR *ptypeDates,
    /* [in] */ VARIANT __RPC_FAR *pModFol,
    /* [in] */ VARIANT __RPC_FAR *pCreditGap,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenSchedule_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CashFlows_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCoefs,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CashFlows_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_View_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pObjet,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_View_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Version_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Version_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDS_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pEffectiveDate,
    /* [in] */ double pEndDate,
    /* [in] */ double pSpread,
    /* [defaultvalue][in] */ BSTR pFixingFreq,
    /* [defaultvalue][in] */ BSTR pDayCountFrq,
    /* [defaultvalue][in] */ double pFirst_period_refdate,
    /* [defaultvalue][in] */ double pFixedPayerAmount,
    /* [defaultvalue][in] */ double pFloatingPayerAmount,
    /* [defaultvalue][in] */ BSTR StubRule,
    /* [defaultvalue][in] */ BSTR pCurrency,
    /* [defaultvalue][in] */ BSTR Adjusted,
    /* [defaultvalue][in] */ int CreditDefLag,
    /* [defaultvalue][in] */ BSTR IncludeMatu,
    /* [defaultvalue][in] */ double StartProtection,
    /* [defaultvalue][in] */ double EndProtection,
    /* [defaultvalue][in] */ BSTR name,
    /* [defaultvalue][in] */ double binary,
    /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
    /* [defaultvalue][in][optional] */ BSTR StartAdj,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDS_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetVolFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pIndex,
    /* [in] */ VARIANT __RPC_FAR *pCcy,
    /* [in] */ VARIANT __RPC_FAR *pCvName,
    /* [in] */ VARIANT __RPC_FAR *pDate,
    /* [in] */ VARIANT __RPC_FAR *pType,
    /* [in] */ VARIANT __RPC_FAR *pMatIndex,
    /* [in] */ VARIANT __RPC_FAR *pImpOrHist,
    /* [defaultvalue][in][optional] */ BSTR pindexId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetVolFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMParallelShift_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZc,
    /* [in] */ double pBump,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMParallelShift_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBumpVolatility_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pVolCurv,
    /* [in] */ double pBump,
    /* [in] */ long pNthLine,
    /* [in] */ long pNthCol,
    /* [in] */ BSTR pIsCumul,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBumpVolatility_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBsSmiledModel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ double pSpot,
    /* [in] */ BSTR pDividend,
    /* [in] */ BSTR pDiscrate,
    /* [in] */ BSTR pVolATM,
    /* [in] */ BSTR pRo,
    /* [in] */ BSTR pNu,
    /* [in] */ BSTR pIsSABR,
    /* [defaultvalue][in][optional] */ BSTR pBeta,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBsSmiledModel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSetEtoolkit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pUserName,
    /* [in] */ BSTR pPassWord,
    /* [in] */ BSTR pDatabaseContext,
    /* [in] */ BSTR pItConfigDomainDir,
    /* [in] */ BSTR pItDomainName,
    /* [retval][out] */ long __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSetEtoolkit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMConnectionEtoolkit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMConnectionEtoolkit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMVolFlat_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pVol,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMVolFlat_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMVolCube_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pATMVol,
    /* [in] */ VARIANT __RPC_FAR *pSmileCurveIds,
    /* [in] */ VARIANT __RPC_FAR *pTenors,
    /* [defaultvalue][in] */ BSTR pVolType,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMVolCube_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMDeconnectionEtoolkit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMDeconnectionEtoolkit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMZcFlat_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pZc,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMZcFlat_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBsModel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ double pSpot,
    /* [in] */ BSTR pDividend,
    /* [in] */ BSTR pDiscrate,
    /* [in] */ BSTR pVol,
    /* [defaultvalue][in][optional] */ BSTR pTypeStk,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBsModel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwitchToETK_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMSwitchToETK_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwitchToFLATFILE_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMSwitchToFLATFILE_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMZCLINT_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [in] */ VARIANT __RPC_FAR *pRate,
    /* [in] */ BSTR pMeth,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCurrency,
    /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMZCLINT_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Pricer_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSecurity,
    /* [in] */ BSTR pModel,
    /* [in] */ BSTR pPricerType,
    /* [in] */ int l_nbpaths,
    /* [defaultvalue][in][optional] */ BSTR pParameters,
    /* [defaultvalue][in][optional] */ double valuationdate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Pricer_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_zcspreaded_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR zcSprId,
    /* [in] */ BSTR zcInitId,
    /* [in] */ double date,
    /* [in] */ BSTR MMFreq,
    /* [in] */ BSTR SwapFreq,
    /* [in] */ BSTR ccyId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_zcspreaded_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CptBaseCorrelation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOf,
    /* [in] */ BSTR name,
    /* [in] */ BSTR CalMethod,
    /* [in] */ BSTR IndexId,
    /* [in] */ VARIANT __RPC_FAR *pStrikeLow,
    /* [in] */ VARIANT __RPC_FAR *pStrikeHigh,
    /* [in] */ VARIANT __RPC_FAR *pVMktBid,
    /* [in] */ VARIANT __RPC_FAR *pVMktAsk,
    /* [in] */ VARIANT __RPC_FAR *pVUpfBid,
    /* [in] */ VARIANT __RPC_FAR *pVUpfAsk,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVInitialCorrel,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVDeltaLevrage,
    /* [defaultvalue][in][optional] */ BSTR ModelId,
    /* [defaultvalue][in][optional] */ double integrationStep,
    /* [defaultvalue][in][optional] */ double lagStartDate,
    /* [defaultvalue][in][optional] */ double creditLag,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pVectorPrevIndexId,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pMatrixPrevBC,
    /* [defaultvalue][in][optional] */ double step,
    /* [defaultvalue][in][optional] */ BSTR CalMeth,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CptBaseCorrelation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDO2_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pEffectiveDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pPortfolio,
    /* [in] */ double pSpread,
    /* [in] */ double pSubAmount,
    /* [in] */ double pMezzAmount,
    /* [defaultvalue][in][optional] */ BSTR pFreqFeeLeg,
    /* [defaultvalue][in][optional] */ BSTR pFreqDefLeg,
    /* [defaultvalue][in][optional] */ BSTR pDayCountFrq,
    /* [defaultvalue][in][optional] */ double pFirst_period_refdate,
    /* [defaultvalue][in][optional] */ BSTR pAccruedOnDefault,
    /* [defaultvalue][in][optional] */ BSTR Currency,
    /* [defaultvalue][in][optional] */ double pPayCreditLag,
    /* [defaultvalue][in][optional] */ BSTR pStub,
    /* [defaultvalue][in][optional] */ double pBinary,
    /* [defaultvalue][in][optional] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
    /* [defaultvalue][in][optional] */ double TradedNotional,
    /* [defaultvalue][in][optional] */ BSTR CrossSub,
    /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
    /* [defaultvalue][in] */ double pFstCpnEffDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDO2_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Portfolio_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pSecuritiesID,
    /* [defaultvalue][in][optional] */ BSTR Parameters,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Portfolio_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCreateZCSwapInt_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [in] */ VARIANT __RPC_FAR *pRate,
    /* [in] */ BSTR pMMVsFut,
    /* [in] */ BSTR pSwapVsFut,
    /* [in] */ BSTR pRaw,
    /* [in] */ BSTR pInterp,
    /* [in] */ BSTR pCcy,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCreateZCSwapInt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetInitialCurveFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pCvName,
    /* [in] */ double pDate,
    /* [defaultvalue][in] */ BSTR pAdjOrNot,
    /* [out] */ VARIANT __RPC_FAR *pRetMat,
    /* [out] */ VARIANT __RPC_FAR *pRetRate);


void __RPC_STUB IARMModule_ARMGetInitialCurveFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMTHREEMONTHFUT_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDelivery,
    /* [in] */ long pMarket,
    /* [in] */ BSTR pCcy,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMTHREEMONTHFUT_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFutPibor_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDelivery,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFutPibor_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMIRFUT_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDelivery,
    /* [in] */ BSTR pIdUnderlying,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMIRFUT_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLibor_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pLiborTypeId,
    /* [in] */ BSTR pCcyId,
    /* [in] */ BSTR pResetFreqId,
    /* [in] */ BSTR pPayFreqId,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [defaultvalue][in] */ BSTR pBasis,
    /* [defaultvalue][in] */ BSTR pIntrule,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLibor_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLiborSwaption_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pReceiveOrPay,
    /* [in] */ double pStrike,
    /* [in] */ double pMaturity,
    /* [in] */ BSTR pLiborType,
    /* [in] */ double pSpread,
    /* [in] */ BSTR pExerciseType,
    /* [in] */ BSTR pResetFreq,
    /* [in] */ BSTR pPayFreq,
    /* [in] */ BSTR pCcyId,
    /* [defaultvalue][in] */ BSTR pRefObj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLiborSwaption_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFixedLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pReceiveOrPay,
    /* [in] */ double pFixRate,
    /* [defaultvalue][in][optional] */ BSTR pDayCount,
    /* [defaultvalue][in][optional] */ BSTR pFreq,
    /* [defaultvalue][in][optional] */ BSTR pDecompFreq,
    /* [defaultvalue][in][optional] */ BSTR pPayTiming,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ BSTR pCcyId,
    /* [defaultvalue][in][optional] */ BSTR pPayCalName,
    /* [defaultvalue][in][optional] */ BSTR pNxChange,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFixedLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetCorrelationMatrix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pModelMultiCurvesId,
    /* [in] */ BSTR pCorrMatrixId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetCorrelationMatrix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DPMktDataFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ BSTR Issuer,
    /* [in] */ BSTR CurveName,
    /* [in] */ BSTR Parameter,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DPMktDataFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetDPFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ BSTR Issuer,
    /* [in] */ BSTR CurveName,
    /* [in] */ BSTR ircurveId,
    /* [in] */ BSTR label,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetDPFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CloneCorrMatrixBary_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrMatrixId,
    /* [in] */ double Beta,
    /* [in] */ int UpOrDown,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CloneCorrMatrixBary_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ConstantDefaultCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ VARIANT __RPC_FAR *pTenors,
    /* [in] */ VARIANT __RPC_FAR *pRates,
    /* [in] */ double Recovery,
    /* [in] */ BSTR IRCurveId,
    /* [in] */ BSTR Ccy,
    /* [in] */ BSTR Label,
    /* [defaultvalue][in] */ BSTR AdjCalType,
    /* [defaultvalue][in] */ BSTR IsSummit,
    /* [defaultvalue][optional][in] */ BSTR calibrationData,
    /* [defaultvalue][optional][in] */ int lag,
    /* [optional][in] */ BSTR calibrationAlgo,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ConstantDefaultCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DefProbModelNew_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDefCurve,
    /* [in] */ BSTR pIRcurve,
    /* [defaultvalue][in][optional] */ BSTR VolCurve,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DefProbModelNew_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ BSTR bIssuer,
    /* [in] */ BSTR bCurrency,
    /* [in] */ BSTR bCvName,
    /* [in] */ BSTR IRCurveId,
    /* [in] */ BSTR bLabel,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBsSlModel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pZc,
    /* [in] */ BSTR pVolSpreadLock,
    /* [in] */ BSTR pCvCapVol,
    /* [in] */ BSTR pCvIndexVol,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBsSlModel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFXVolFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pCcy2,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCvName,
    /* [defaultvalue][in][optional] */ BSTR pType,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFXVolFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFXCorrelFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCcy2,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCvName,
    /* [in] */ VARIANT __RPC_FAR *pTenors,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFXCorrelFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMInfocentreConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMInfocentreConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGlobDFBS_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDomBSId,
    /* [in] */ BSTR pDomCurrId,
    /* [in] */ BSTR pFrgBSId,
    /* [in] */ BSTR pFrgCurrId,
    /* [in] */ BSTR pFxVolCrvId,
    /* [in] */ BSTR pFFxCorrId,
    /* [in] */ BSTR pRatesCorrId,
    /* [defaultvalue][in][optional] */ BSTR pFxVolModelId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGlobDFBS_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetInitialCurveFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pIndex,
    /* [in] */ VARIANT __RPC_FAR *pCurrency,
    /* [in] */ VARIANT __RPC_FAR *pCvName,
    /* [in] */ VARIANT __RPC_FAR *pDate,
    /* [in] */ VARIANT __RPC_FAR *value,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetInitialCurveFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetCorrelFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex1,
    /* [in] */ BSTR pCcy2,
    /* [in] */ BSTR pIndex2,
    /* [in] */ double pDate,
    /* [in] */ BSTR pCvName,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetCorrelFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMDFFXBS_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDVolId,
    /* [in] */ BSTR pFVolId,
    /* [in] */ BSTR pDZcId,
    /* [in] */ BSTR pFZcId,
    /* [in] */ BSTR pDFxCorrId,
    /* [in] */ BSTR pFFxCorrId,
    /* [in] */ BSTR pFxVolId,
    /* [in] */ double pRatesCorr,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMDFFXBS_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMTRIBSMODEL_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pModel1,
    /* [in] */ BSTR pModel2,
    /* [in] */ BSTR pDiscModel,
    /* [in] */ BSTR pFX1DiscVol,
    /* [in] */ BSTR pFX2DiscVol,
    /* [in] */ BSTR pIdx1Idx2Corr,
    /* [in] */ BSTR pIdx1DiscIdxCorr,
    /* [in] */ BSTR pIdx2DiscIdxCorr,
    /* [in] */ BSTR pIdx1FxCorr,
    /* [in] */ BSTR pIdx2FxCorr,
    /* [in] */ int pQuantoFlag,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMTRIBSMODEL_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMTRIBSDUAL_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pModel1,
    /* [in] */ BSTR pModel2,
    /* [in] */ BSTR pDiscModel,
    /* [in] */ BSTR pFX1DiscVol,
    /* [in] */ BSTR pFX2DiscVol,
    /* [in] */ BSTR pIdx1Idx2Corr,
    /* [in] */ BSTR pIdx1DiscIdxCorr,
    /* [in] */ BSTR pIdx2DiscIdxCorr,
    /* [in] */ BSTR pIdx1FxCorr,
    /* [in] */ BSTR pIdx2FxCorr,
    /* [in] */ int pQuantoFlag,
    /* [in] */ double pCorrelForAdj,
    /* [in] */ int pWithslopeflag,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMTRIBSDUAL_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcptBonibor_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pDepart,
    /* [in] */ BSTR pMatStruct,
    /* [in] */ BSTR pMatTot,
    /* [in] */ BSTR pAmort,
    /* [in] */ BSTR pFreq,
    /* [in] */ BSTR pSjUSD,
    /* [in] */ BSTR pTiming,
    /* [in] */ double pBarriere,
    /* [in] */ double pSpdPostBar,
    /* [in] */ double pMarge,
    /* [in] */ double pFunding,
    /* [in] */ BSTR pFundingFreq,
    /* [in] */ double pSpd2phase,
    /* [in] */ double pSoulte,
    /* [in] */ BSTR pYcModId,
    /* [in] */ BSTR pBsModId,
    /* [in] */ BSTR pBsModVolUSDId,
    /* [in] */ BSTR pBsModCorrPlusId,
    /* [in] */ BSTR pBsModCorrMoinsId,
    /* [in] */ BSTR pCrossModId,
    /* [in] */ BSTR pProbaMarge,
    /* [defaultvalue][in][optional] */ double pInt,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcptBonibor_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CMTranche_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
    /* [in] */ VARIANT __RPC_FAR *pEndDate,
    /* [in] */ VARIANT __RPC_FAR *pParticipationRate,
    /* [in] */ VARIANT __RPC_FAR *pMezzAmount,
    /* [in] */ VARIANT __RPC_FAR *pSubAmount,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pNotionals,
    /* [in] */ VARIANT __RPC_FAR *pIndex,
    /* [in] */ VARIANT __RPC_FAR *pFreqFeeLeg,
    /* [in] */ VARIANT __RPC_FAR *pDayCount,
    /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
    /* [in] */ VARIANT __RPC_FAR *pAccruedOnDefault,
    /* [in] */ VARIANT __RPC_FAR *Currency,
    /* [in] */ VARIANT __RPC_FAR *pPayCreditLag,
    /* [in] */ VARIANT __RPC_FAR *pStub,
    /* [in] */ VARIANT __RPC_FAR *pFreqDefLeg,
    /* [in] */ VARIANT __RPC_FAR *pBinary,
    /* [in] */ VARIANT __RPC_FAR *pPayCal,
    /* [defaultvalue][in][optional] */ BSTR LongOrShortRisk,
    /* [defaultvalue][in][optional] */ double TradedNotional,
    /* [defaultvalue][in][optional] */ double FwdFixedDate,
    /* [defaultvalue][in][optional] */ BSTR IncludeMatu,
    /* [defaultvalue][in] */ double pFstCpnEffDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CMTranche_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Index_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [defaultvalue][in] */ double YearFrac,
    /* [defaultvalue][in] */ double pSpread,
    /* [defaultvalue][in] */ BSTR Method,
    /* [defaultvalue][in] */ BSTR Basis,
    /* [defaultvalue][in] */ BSTR ResetFreq,
    /* [defaultvalue][in] */ BSTR PayFreq,
    /* [defaultvalue][in] */ BSTR ccy,
    /* [defaultvalue][in] */ BSTR DefCurve,
    /* [defaultvalue][in] */ BSTR fwdRule,
    /* [defaultvalue][in] */ BSTR resetTiming,
    /* [defaultvalue][in] */ int resetGap,
    /* [defaultvalue][in] */ BSTR payTiming,
    /* [defaultvalue][in] */ int payGap,
    /* [defaultvalue][in] */ BSTR intRule,
    /* [defaultvalue][in] */ BSTR AdjCalType,
    /* [defaultvalue][in] */ int cm_resetWeekDay,
    /* [defaultvalue][in] */ int cm_resetOccur,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Index_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Parameters_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pCoefs,
    /* [in] */ long nbcols,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Parameters_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMAswPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pMaturity,
    /* [in] */ double pCpn,
    /* [in] */ BSTR pFreq,
    /* [in] */ BSTR pBase,
    /* [in] */ double pMargin,
    /* [defaultvalue][in][optional] */ double pRedemptionPrice,
    /* [in] */ double pAsOf,
    /* [in] */ double pDelivery,
    /* [in] */ BSTR pFixDecompfreq,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex1,
    /* [in] */ BSTR pFwdCurve1,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
    /* [defaultvalue][in][optional] */ BSTR pCcy2,
    /* [defaultvalue][in][optional] */ BSTR pIndex2,
    /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
    /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
    /* [defaultvalue][in][optional] */ long pSolve,
    /* [defaultvalue][in][optional] */ double pMinValue,
    /* [defaultvalue][in][optional] */ double pMaxValue,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMAswPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMAswMargin_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pMaturity,
    /* [in] */ double pCpn,
    /* [in] */ BSTR pFreq,
    /* [in] */ BSTR pBase,
    /* [in] */ double pPrice,
    /* [defaultvalue][in][optional] */ double pRedemptionPrice,
    /* [in] */ double pAsOf,
    /* [in] */ double pDelivery,
    /* [in] */ BSTR pFixDecompfreq,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex1,
    /* [in] */ BSTR pFwdCurve1,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
    /* [defaultvalue][in][optional] */ BSTR pCcy2,
    /* [defaultvalue][in][optional] */ BSTR pIndex2,
    /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
    /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
    /* [defaultvalue][in][optional] */ long pSolve,
    /* [defaultvalue][in][optional] */ double pMinValue,
    /* [defaultvalue][in][optional] */ double pMaxValue,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMAswMargin_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDSIndex_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pEffectiveDate,
    /* [in] */ VARIANT __RPC_FAR *pEndDate,
    /* [in] */ VARIANT __RPC_FAR *pSpread,
    /* [in] */ VARIANT __RPC_FAR *pIndex,
    /* [in] */ VARIANT __RPC_FAR *pFixingFreq,
    /* [in] */ VARIANT __RPC_FAR *pDayCountFrq,
    /* [in] */ VARIANT __RPC_FAR *pFirst_period_refdate,
    /* [in] */ VARIANT __RPC_FAR *pFixedPayerAmount,
    /* [in] */ VARIANT __RPC_FAR *pFloatingPayerAmount,
    /* [in] */ BSTR StubRule,
    /* [in] */ VARIANT __RPC_FAR *pCurrency,
    /* [defaultvalue][in] */ BSTR Adjusted,
    /* [defaultvalue][in] */ int CreditLag,
    /* [defaultvalue][in] */ BSTR IncludeMaturity,
    /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionStartDate,
    /* [defaultvalue][in] */ VARIANT __RPC_FAR *ProtectionEndDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDSIndex_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Option_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *UnderlyingMaturity,
    /* [in] */ double OptionExpiry,
    /* [in] */ BSTR Currency,
    /* [defaultvalue][in][optional] */ BSTR CdsAdj,
    /* [defaultvalue][in][optional] */ BSTR EndAdj,
    /* [defaultvalue][in][optional] */ double pStrike,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pOptionType,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pKoType,
    /* [defaultvalue][in][optional] */ double Notional,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Option_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FwdSpreadPricer_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ double Maturity1,
    /* [in] */ double Maturity2,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_FwdSpreadPricer_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ImpliedVol_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ double pMktPrice,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ImpliedVol_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_VirtualCdsSpread_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pMaturity,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_VirtualCdsSpread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_BSGreeks_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ VARIANT __RPC_FAR *pGreekType,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_BSGreeks_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetCorrelation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pModelMultiCurvesId,
    /* [in] */ BSTR pCorrelationId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetCorrelation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CorrelationStrike_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pVolCurves,
    /* [in] */ VARIANT __RPC_FAR *pProportions,
    /* [in] */ VARIANT __RPC_FAR *pSmileStrikeLow,
    /* [in] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
    /* [in] */ VARIANT __RPC_FAR *pIndexVector,
    /* [defaultvalue][in] */ double AsOf,
    /* [defaultvalue][in] */ BSTR Name,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CorrelationStrike_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Beta_Correlation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pCoefs,
    /* [in] */ double AsOf,
    /* [defaultvalue][in] */ BSTR Name,
    /* [defaultvalue][in] */ BSTR idIndex1,
    /* [defaultvalue][in] */ BSTR idIndex2,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Beta_Correlation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_GETINSTRUMENTFROMSUMMIT_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR SummitId,
    /* [in] */ BSTR Type,
    /* [defaultvalue][in][optional] */ double AsOf,
    /* [defaultvalue][in][optional] */ BSTR ExoticFilter,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_GETINSTRUMENTFROMSUMMIT_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFrnPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pDelivery,
    /* [in] */ double pMaturity,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex1,
    /* [in] */ BSTR pFwdCurve1,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
    /* [in] */ double pFacialMargin,
    /* [in] */ double pValoMargin,
    /* [defaultvalue][in][optional] */ BSTR pCcy2,
    /* [defaultvalue][in][optional] */ BSTR pIndex2,
    /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
    /* [defaultvalue][in][optional] */ double pFixing,
    /* [defaultvalue][in][optional] */ double pSpread,
    /* [defaultvalue][in][optional] */ double pOutMode,
    /* [defaultvalue][in][optional] */ long pSolve,
    /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFrnPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFrnMargin_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pDelivery,
    /* [in] */ double pMaturity,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pIndex1,
    /* [in] */ BSTR pFwdCurve1,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve1,
    /* [in] */ double pFacialMargin,
    /* [in] */ double pPrice,
    /* [defaultvalue][in][optional] */ BSTR pCcy2,
    /* [defaultvalue][in][optional] */ BSTR pIndex2,
    /* [defaultvalue][in][optional] */ BSTR pFwdCurve2,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve2,
    /* [defaultvalue][in][optional] */ double pFixing,
    /* [defaultvalue][in][optional] */ double pSpread,
    /* [defaultvalue][in][optional] */ double pOutMode,
    /* [defaultvalue][in][optional] */ long pSolve,
    /* [defaultvalue][in][optional] */ BSTR pAmortizationId,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFrnMargin_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetEqStrikeDown_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrelId,
    /* [in] */ BSTR IndexName,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetEqStrikeDown_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetEqStrikeUp_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrelId,
    /* [in] */ BSTR IndexName,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetEqStrikeUp_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetCorrelStrikeDown_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrelId,
    /* [in] */ double yfmaturity,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetCorrelStrikeDown_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetCorrelStrikeUp_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrelId,
    /* [in] */ double yfmaturity,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetCorrelStrikeUp_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetCorrelation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR ModelId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetCorrelation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetModelFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR IRcurve,
    /* [in] */ BSTR IDSummit,
    /* [in] */ BSTR type,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetModelFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetVolatility_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ BSTR VolCurveId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetVolatility_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_NextCpnDate_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ double maturity,
    /* [in] */ BSTR frequency,
    /* [defaultvalue][in][optional] */ BSTR rule,
    /* [defaultvalue][in][optional] */ BSTR currency,
    /* [defaultvalue][in][optional] */ BSTR intrule,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_NextCpnDate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetProportionsInfos_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR correlId,
    /* [in] */ BSTR IndexName,
    /* [in] */ double proportion,
    /* [defaultvalue][in][optional] */ double forcedstrikelow,
    /* [defaultvalue][in][optional] */ double forcedstrikehigh);


void __RPC_STUB IARMModule_ARM_Credit_SetProportionsInfos_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcptDigital_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ double pNtl,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pBsmod,
    /* [in] */ BSTR pBsmodDelta,
    /* [in] */ BSTR pBsmodVega,
    /* [in] */ BSTR pFreqP,
    /* [in] */ BSTR pResetTiming,
    /* [in] */ BSTR pPorR,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pCcyIdx,
    /* [in] */ BSTR pDayCount,
    /* [in] */ BSTR pCapOrFloor,
    /* [in] */ BSTR pAmort,
    /* [in] */ BSTR pStrike,
    /* [in] */ BSTR pPayOff,
    /* [in] */ BSTR pSpd,
    /* [in] */ double pResetGap,
    /* [defaultvalue][in][optional] */ double pSpreadBelow,
    /* [defaultvalue][in][optional] */ double pSpreadAbove,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
    /* [defaultvalue][in][optional] */ double pTxAmort,
    /* [defaultvalue][in][optional] */ double pAmountAmort,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcptDigital_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMRefValue_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pdates,
    /* [in] */ VARIANT __RPC_FAR *pvalues,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pvalues2,
    /* [defaultvalue][in][optional] */ long valueType,
    /* [defaultvalue][in][optional] */ long conversion,
    /* [defaultvalue][in][optional] */ BSTR calcMethod,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMRefValue_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCreateGenCorrelManager_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pMktTags,
    /* [in] */ VARIANT __RPC_FAR *pIntraMktTags,
    /* [in] */ VARIANT __RPC_FAR *pCorrelCurveIds,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCreateGenCorrelManager_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBSConvAdjust_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [defaultvalue][in][optional] */ BSTR pSUMMITFormulaeUsed,
    /* [defaultvalue][in][optional] */ BSTR pUseSABRCMS,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBSConvAdjust_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBsModelGen_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pYieldCurve,
    /* [in] */ BSTR pVolatility,
    /* [defaultvalue][in][optional] */ BSTR pCorrMgr,
    /* [defaultvalue][in][optional] */ BSTR pCnvxManager,
    /* [defaultvalue][in][optional] */ BSTR pCapletVol,
    /* [defaultvalue][in][optional] */ BSTR pSpreadLock,
    /* [defaultvalue][in][optional] */ BSTR pDiscCurve,
    /* [defaultvalue][in][optional] */ BSTR pCorrel,
    /* [defaultvalue][in][optional] */ BSTR pCashVol,
    /* [defaultvalue][in][optional] */ BSTR pSpreadVol,
    /* [defaultvalue][in][optional] */ BSTR pModelType,
    /* [defaultvalue][in][optional] */ BSTR pSpreadVolType,
    /* [defaultvalue][in][optional] */ BSTR pSabrMod,
    /* [defaultvalue][in][optional] */ BSTR pLnorNorVol,
    /* [defaultvalue][in][optional] */ long pNumSteps,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBsModelGen_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcptSPTQTF_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pStartDate,
    /* [in] */ double pStartDatePhase2,
    /* [in] */ double pStartDatePhase3,
    /* [in] */ double pEndDate,
    /* [in] */ double pNtl,
    /* [in] */ BSTR pIndexPhase1,
    /* [in] */ BSTR pIndexPhase2,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ BSTR pIndexPhase3,
    /* [in] */ BSTR pFreqPPhase1,
    /* [in] */ BSTR pFreqPPhase2,
    /* [in] */ BSTR pFreqPFund,
    /* [in] */ BSTR pFreqPPhase3,
    /* [in] */ BSTR pFreqR,
    /* [in] */ BSTR pResetTimingPhase1,
    /* [in] */ BSTR pResetTimingPhase2,
    /* [in] */ BSTR pResetTimingPhase3,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pCcyIdx,
    /* [in] */ BSTR pDayCount,
    /* [in] */ double pFee,
    /* [in] */ BSTR pIsRateFixedPhase2,
    /* [in] */ double pFixedRatePhase2,
    /* [in] */ BSTR pBarrier,
    /* [in] */ BSTR pSpdPhase1,
    /* [in] */ BSTR pSpdPhase1Fund,
    /* [in] */ BSTR pSpdPhase2Tf,
    /* [in] */ BSTR pSpdPhase2fund,
    /* [in] */ BSTR pSpdPhase3,
    /* [in] */ BSTR pSpdPhase3fund,
    /* [in] */ double pResetGapPhase1,
    /* [in] */ double pResetGapPhase2,
    /* [in] */ double pResetGapPhase3,
    /* [in] */ BSTR pAmort,
    /* [in] */ BSTR pBsmod,
    /* [in] */ BSTR pBsmodDeltaCcy1,
    /* [in] */ BSTR pBsmodVegaCcy1,
    /* [defaultvalue][in][optional] */ BSTR pBsmodDeltaCcy2,
    /* [defaultvalue][in][optional] */ BSTR pBsmodVegaCcy2,
    /* [defaultvalue][in][optional] */ BSTR pBsmodFxCorrel,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
    /* [defaultvalue][in][optional] */ double pTxAmort,
    /* [defaultvalue][in][optional] */ double pAmountAmort,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcptSPTQTF_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMDisplaySchedule_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR legId,
    /* [in] */ BSTR dataType,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMDisplaySchedule_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMIrIndex_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [defaultvalue][in][optional] */ BSTR pDaycount,
    /* [defaultvalue][in][optional] */ BSTR pPayFreq,
    /* [defaultvalue][in][optional] */ double pMaturity,
    /* [defaultvalue][in][optional] */ BSTR pCompMethod,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pResetTiming,
    /* [defaultvalue][in][optional] */ double pResetGap,
    /* [defaultvalue][in][optional] */ BSTR pPayTiming,
    /* [defaultvalue][in][optional] */ double pPayGap,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR pIndexType,
    /* [defaultvalue][in][optional] */ double pDecompFreq,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pResetFreq,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMIrIndex_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwapleg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndexId,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pRecOrPay,
    /* [in] */ VARIANT pSpread,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR pDayCount,
    /* [defaultvalue][in][optional] */ double pResetGap,
    /* [defaultvalue][in][optional] */ BSTR pResetCal,
    /* [defaultvalue][in][optional] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
    /* [defaultvalue][in][optional] */ BSTR pNxChange,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSwapleg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMConstRefvalue_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pValue,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMConstRefvalue_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CptImplCvForCDO2_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricerId,
    /* [in] */ BSTR Name,
    /* [in] */ BSTR Tenor,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CptImplCvForCDO2_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_AddPeriod_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ BSTR Maturity,
    /* [in] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR AdjRule,
    /* [defaultvalue][in][optional] */ BSTR AdjCDS,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_AddPeriod_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMShutdownETK_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMShutdownETK_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetCoupons_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CdsorCdoId,
    /* [in] */ BSTR CouponsId,
    /* [in] */ BSTR TypesId,
    /* [defaultvalue][in][optional] */ BSTR PartId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetCoupons_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_InputDefaultCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ VARIANT __RPC_FAR *pDates,
    /* [in] */ VARIANT __RPC_FAR *pRates,
    /* [in] */ double Recovery,
    /* [in] */ BSTR IRCurveId,
    /* [in] */ BSTR bCurrency,
    /* [in] */ BSTR bLabel,
    /* [defaultvalue][in][optional] */ BSTR bInterpolType,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_InputDefaultCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetInitialVolFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pCvName,
    /* [in] */ double pDate,
    /* [in] */ BSTR pType,
    /* [defaultvalue][in][optional] */ BSTR pMatIndex,
    /* [out] */ VARIANT __RPC_FAR *pRetMat,
    /* [out] */ VARIANT __RPC_FAR *pRetTenor,
    /* [out] */ VARIANT __RPC_FAR *pRetVol);


void __RPC_STUB IARMModule_ARMGetInitialVolFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBond_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pIssueDate,
    /* [in] */ double pMaturityDate,
    /* [in] */ double pFirstCpnDate,
    /* [in] */ double pCpnRate,
    /* [in] */ double pRedempPrice,
    /* [in] */ double pPeriodicity,
    /* [in] */ VARIANT pDaycount,
    /* [defaultvalue][in][optional] */ double pSettleGap,
    /* [defaultvalue][in][optional] */ double pCpnDateFlag,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBond_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFCreateOATLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pInfIdx,
    /* [in] */ BSTR pRcvOrPay,
    /* [defaultvalue][in][optional] */ BSTR pInterpType,
    /* [defaultvalue][in][optional] */ double pLeverage,
    /* [defaultvalue][in][optional] */ double pSpread,
    /* [defaultvalue][in][optional] */ BSTR pResetFreq,
    /* [defaultvalue][in][optional] */ BSTR pDaycount,
    /* [defaultvalue][in][optional] */ BSTR pResetCal,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ double pResetNumGap,
    /* [defaultvalue][in][optional] */ double pResetDenomGap,
    /* [defaultvalue][in][optional] */ BSTR pPayFreq,
    /* [defaultvalue][in][optional] */ double pPayGap,
    /* [defaultvalue][in][optional] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ BSTR pFinalNotionalType,
    /* [defaultvalue][in][optional] */ double pFirstReset,
    /* [defaultvalue][in][optional] */ double pCoMultiple,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFCreateOATLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwap_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSwapleg1,
    /* [in] */ BSTR pSwapleg2,
    /* [defaultvalue][in][optional] */ double pMinPay,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSwap_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMPToYield_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pBond,
    /* [in] */ double pSettleDate,
    /* [in] */ double pPrice,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMPToYield_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMYToPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pBond,
    /* [in] */ double pSettleDate,
    /* [in] */ double pYield,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMYToPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMYToDuration_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pBond,
    /* [in] */ double pSettleDate,
    /* [in] */ double pActuRate,
    /* [defaultvalue][in][optional] */ double pFlagCpn,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMYToDuration_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLiborleg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pLiborType,
    /* [in] */ BSTR pRecOrPay,
    /* [defaultvalue][in][optional] */ VARIANT pSpread,
    /* [defaultvalue][in][optional] */ BSTR pResetFReq,
    /* [defaultvalue][in][optional] */ BSTR pPayFreq,
    /* [defaultvalue][in][optional] */ BSTR pResetTiming,
    /* [defaultvalue][in][optional] */ BSTR pPayTiming,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ double pResetGap,
    /* [defaultvalue][in][optional] */ BSTR pResetCal,
    /* [defaultvalue][in][optional] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
    /* [defaultvalue][in][optional] */ BSTR pNxChange,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
    /* [defaultvalue][in][optional] */ BSTR pCpnDaycount,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLiborleg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMImpliedSpread_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSwap,
    /* [in] */ BSTR pModel,
    /* [in] */ double pPrice,
    /* [in] */ double pLeg1or2,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMImpliedSpread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMDiscountPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZeroCurve,
    /* [in] */ double pMatu,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMDiscountPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFCreateCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ BSTR pIndexName,
    /* [in] */ double pCPIIndexValue,
    /* [in] */ double pCPIIndexDate,
    /* [in] */ VARIANT __RPC_FAR *pMatu,
    /* [in] */ VARIANT __RPC_FAR *pRate,
    /* [defaultvalue][in][optional] */ BSTR pMonthlyInterpType,
    /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
    /* [defaultvalue][in][optional] */ BSTR pDCFMonthly,
    /* [defaultvalue][in][optional] */ BSTR pDCFDaily,
    /* [defaultvalue][in][optional] */ BSTR pExtrapolType,
    /* [defaultvalue][in][optional] */ BSTR pResetManager,
    /* [defaultvalue][in][optional] */ BSTR pSeasonManager,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFCreateCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFInterpCPI_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZc,
    /* [in] */ double pCPIDate,
    /* [defaultvalue][in][optional] */ BSTR pDCFlag,
    /* [defaultvalue][in][optional] */ BSTR pDailyInterpType,
    /* [defaultvalue][in][optional] */ BSTR pCPIlag,
    /* [defaultvalue][in][optional] */ double pWeight,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFInterpCPI_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFSeasonManager_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pMonthList,
    /* [in] */ VARIANT __RPC_FAR *pValues,
    /* [defaultvalue][in][optional] */ BSTR pSeasonAdjMode,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFSeasonManager_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFResetManager_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pDatas,
    /* [in] */ double pNbIndex,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFResetManager_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMINFYcMod_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pYieldCurve,
    /* [in] */ BSTR pInfCurve,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMINFYcMod_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBaseReplicationConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMBaseReplicationConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetInitialFXVolFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pCcy2,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pCvName,
    /* [defaultvalue][in][optional] */ BSTR pImpOrHist,
    /* [defaultvalue][in][optional] */ BSTR pVolType,
    /* [out] */ VARIANT __RPC_FAR *pRetMat,
    /* [out] */ VARIANT __RPC_FAR *pRetTenor,
    /* [out] */ VARIANT __RPC_FAR *pRetVol);


void __RPC_STUB IARMModule_ARMGetInitialFXVolFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCreateZCFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pCvName,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pAdj,
    /* [defaultvalue][in][optional] */ BSTR pRaw,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCreateZCFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBumpCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pZc,
    /* [in] */ double pEpsilon,
    /* [defaultvalue][in][optional] */ long pMethod,
    /* [defaultvalue][in][optional] */ BSTR pPlot,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBumpCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMAccrued_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSec,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pModel,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMAccrued_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSwitchToWSETK_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMSwitchToWSETK_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DataFromLabel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricer,
    /* [in] */ BSTR label,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DataFromLabel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeBilibor_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pStartDate,
    /* [in] */ double pDateSecondPhase,
    /* [in] */ double pEndDate,
    /* [in] */ double pNtl,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ BSTR pIndexSF,
    /* [in] */ BSTR pBsmod,
    /* [in] */ BSTR pBsmodFund,
    /* [in] */ BSTR pBsmodDeltaCcy1,
    /* [in] */ BSTR pBsmodDeltaFund,
    /* [in] */ BSTR pBsmodDeltaCcy2,
    /* [in] */ BSTR pBsmodFxCorrel,
    /* [in] */ BSTR pFreqP,
    /* [in] */ BSTR pFreqR,
    /* [in] */ BSTR pFreqPFund,
    /* [in] */ BSTR pFreqRFund,
    /* [in] */ BSTR pFreqPSF,
    /* [in] */ BSTR pFreqRSF,
    /* [in] */ BSTR pResetTiming,
    /* [in] */ BSTR pResetTimingSF,
    /* [in] */ BSTR pCcy1,
    /* [in] */ BSTR pCcy2,
    /* [in] */ BSTR pDayCount,
    /* [in] */ BSTR pDayCountSF,
    /* [in] */ double pSpdPF,
    /* [in] */ double pSpdSF,
    /* [in] */ double pSpdfund,
    /* [in] */ double pSpdfund2,
    /* [in] */ double pResetGap,
    /* [in] */ double pResetGapSF,
    /* [in] */ BSTR pAmort,
    /* [in] */ double pRefDate,
    /* [in] */ double pFee,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
    /* [defaultvalue][in][optional] */ double pTxAmort,
    /* [defaultvalue][in][optional] */ double pAmountAmort,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeBilibor_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenerateImpliedCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricerId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenerateImpliedCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetEqStrike_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CorrelId,
    /* [in] */ BSTR IndexName,
    /* [in] */ BSTR UpOrLow,
    /* [out] */ VARIANT __RPC_FAR *pRetMatu,
    /* [out] */ VARIANT __RPC_FAR *pRetStrikes);


void __RPC_STUB IARMModule_ARM_Credit_GetEqStrike_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeOptilix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ double pStartDate,
    /* [in] */ double pDateSecondPhase,
    /* [in] */ double pEndDate,
    /* [in] */ double pNtl,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ BSTR pIndexSF,
    /* [in] */ BSTR pBsmod,
    /* [in] */ BSTR pBsmodFund,
    /* [in] */ BSTR pBsmodDeltaCcy,
    /* [in] */ BSTR pBsmodDeltaFund,
    /* [in] */ BSTR pFreqP,
    /* [in] */ BSTR pFreqR,
    /* [in] */ BSTR pFreqPFund,
    /* [in] */ BSTR pFreqRFund,
    /* [in] */ BSTR pFreqPSF,
    /* [in] */ BSTR pFreqRSF,
    /* [in] */ BSTR pResetTiming,
    /* [in] */ BSTR pResetTimingSF,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pDayCount,
    /* [in] */ BSTR pDayCountSF,
    /* [in] */ double pSpdSF,
    /* [in] */ VARIANT pSpdfund,
    /* [in] */ double pResetGap,
    /* [in] */ double pResetGapSF,
    /* [in] */ BSTR pAmort,
    /* [in] */ double pRefDate,
    /* [in] */ double pFee,
    /* [defaultvalue][in][optional] */ BSTR pFwdRule,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ BSTR pFreqAmort,
    /* [defaultvalue][in][optional] */ double pTxAmort,
    /* [defaultvalue][in][optional] */ double pAmountAmort,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeOptilix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DefaultIntensity_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricerId,
    /* [in] */ double Maturity,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DefaultIntensity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMRiskyBond_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pIssueDate,
    /* [in] */ double pMaturityDate,
    /* [in] */ double pFirstCpnDate,
    /* [in] */ double pCpnRate,
    /* [in] */ double pRedemptionPrice,
    /* [in] */ long pPeriodicity,
    /* [in] */ VARIANT pDaycount,
    /* [defaultvalue][in][optional] */ long pSettleGap,
    /* [defaultvalue][in][optional] */ long pCpnDateFlag,
    /* [defaultvalue][in][optional] */ BSTR pCcyId,
    /* [defaultvalue][in][optional] */ double pRepo,
    /* [defaultvalue][in][optional] */ double pSsl,
    /* [defaultvalue][in][optional] */ double pRecoveryRate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMRiskyBond_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMRiskyBondWithCF_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOfDate,
    /* [in] */ double pRedemptionPrice,
    /* [in] */ long pPeriodicity,
    /* [in] */ VARIANT pDaycount,
    /* [in] */ VARIANT __RPC_FAR *pYearTerms,
    /* [in] */ VARIANT __RPC_FAR *pCashFlows,
    /* [defaultvalue][in][optional] */ long pSettleGap,
    /* [defaultvalue][in][optional] */ long pCpnDateFlag,
    /* [defaultvalue][in][optional] */ BSTR pCcyId,
    /* [defaultvalue][in][optional] */ double pRepo,
    /* [defaultvalue][in][optional] */ double pSsl,
    /* [defaultvalue][in][optional] */ double pRecoveryRate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMRiskyBondWithCF_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_EmptyLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_EmptyLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMClonedAndSetNotional_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR bLegId,
    /* [defaultvalue][in][optional] */ BSTR bAmortId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMClonedAndSetNotional_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_INF_GetZcFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR Index,
    /* [in] */ BSTR Ccy,
    /* [in] */ BSTR cvname,
    /* [in] */ double date,
    /* [defaultvalue][in][optional] */ BSTR seasonAdj,
    /* [defaultvalue][in][optional] */ BSTR seasonAdjMode,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_INF_GetZcFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_IRLEGTOCREDITLEG_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR SwapLegId,
    /* [in] */ BSTR LegType,
    /* [in] */ BSTR creditindexId,
    /* [in] */ BSTR pricerId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_IRLEGTOCREDITLEG_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Collateral_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pNotionals,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Collateral_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDSGEN_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR FeeLegId,
    /* [in] */ BSTR DefLegId,
    /* [in] */ double RcvFee,
    /* [in] */ double TradedNot,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDSGEN_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_NTDGEN_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CdsId,
    /* [in] */ int firstnumdef,
    /* [in] */ int lastnumdef,
    /* [in] */ BSTR CollateralId,
    /* [in] */ double binary,
    /* [in] */ double rcvfee,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_NTDGEN_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDOGEN_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CdsId,
    /* [in] */ double subamount,
    /* [in] */ BSTR CollateralId,
    /* [in] */ double binary,
    /* [in] */ double rcvfee,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDOGEN_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double StartDate,
    /* [in] */ double EndDate,
    /* [in] */ double FixedRate,
    /* [in] */ double FixedNotional,
    /* [defaultvalue][in][optional] */ BSTR RefValNotional,
    /* [defaultvalue][in][optional] */ BSTR RefValRate,
    /* [defaultvalue][in][optional] */ BSTR XChangeNotional,
    /* [defaultvalue][in][optional] */ BSTR Frequency,
    /* [defaultvalue][in][optional] */ BSTR Basis,
    /* [defaultvalue][in][optional] */ BSTR payTiming,
    /* [defaultvalue][in][optional] */ BSTR intrule,
    /* [defaultvalue][in][optional] */ BSTR stubrule,
    /* [defaultvalue][in][optional] */ BSTR ccyid,
    /* [defaultvalue][in][optional] */ BSTR paycalname,
    /* [defaultvalue][in][optional] */ double refdate,
    /* [defaultvalue][in][optional] */ BSTR includematurity,
    /* [defaultvalue][in][optional] */ BSTR adjstartdate,
    /* [defaultvalue][in][optional] */ BSTR legtype,
    /* [defaultvalue][in][optional] */ BSTR indexobj,
    /* [defaultvalue][in][optional] */ int creditlag,
    /* [defaultvalue][in][optional] */ double binary,
    /* [defaultvalue][in][optional] */ BSTR name,
    /* [defaultvalue][in][optional] */ BSTR Nxchange,
    /* [defaultvalue][in][optional] */ BSTR baccruedOnDef,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CDO_SQUARE_GEN_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CdsId,
    /* [in] */ double subamount,
    /* [in] */ BSTR portfolioId,
    /* [in] */ double binary,
    /* [in] */ double rcvfee,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CDO_SQUARE_GEN_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetLastFixingDate_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR instId,
    /* [in] */ VARIANT asofDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetLastFixingDate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetPastFixing_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR instId,
    /* [in] */ VARIANT resetDate,
    /* [in] */ VARIANT fixingValue,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetPastFixing_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFixing_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR source,
    /* [in] */ BSTR index,
    /* [in] */ BSTR term,
    /* [in] */ BSTR ccy,
    /* [in] */ double date,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFixing_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetRiskyProfile_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CdsorCdoId,
    /* [in] */ BSTR CouponsId,
    /* [in] */ BSTR TypesId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetRiskyProfile_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMHyperCube_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
    /* [in] */ VARIANT __RPC_FAR *pKeys,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMHyperCube_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetPricerForRatesComputation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR legId,
    /* [in] */ BSTR pricerId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetPricerForRatesComputation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCmsLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double startDate,
    /* [in] */ double endDate,
    /* [in] */ BSTR cmsTypeId,
    /* [defaultvalue][in][optional] */ BSTR receiveOrPay,
    /* [defaultvalue][in][optional] */ BSTR yieldDecompFreq,
    /* [defaultvalue][in][optional] */ BSTR swapLegDayCount,
    /* [defaultvalue][in][optional] */ BSTR resetFreq,
    /* [defaultvalue][in][optional] */ BSTR payFreq,
    /* [defaultvalue][in][optional] */ long resetGap,
    /* [defaultvalue][in][optional] */ BSTR intRule,
    /* [defaultvalue][in][optional] */ BSTR ccyName,
    /* [defaultvalue][in][optional] */ BSTR resetTiming,
    /* [defaultvalue][in][optional] */ BSTR stubRule,
    /* [defaultvalue][in][optional] */ BSTR adjStartDate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCmsLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetMatuLabel_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurveId,
    /* [in] */ VARIANT __RPC_FAR *pMatuLabels,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetMatuLabel_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputePentifix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ double pStartdatePhase1,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexPhase1,
    /* [in] */ double pSpreadPhase1,
    /* [in] */ BSTR pDayCountPhase1,
    /* [in] */ BSTR pPayFreqPhase1,
    /* [in] */ BSTR pResetFreqPhase1,
    /* [in] */ BSTR pResetTimingPhase1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pAdjPhase1,
    /* [in] */ BSTR pStub,
    /* [in] */ BSTR pIndexPhase2DIG,
    /* [in] */ BSTR pIndexLongPhase2DIG,
    /* [in] */ BSTR pStrikePhase2DIG,
    /* [in] */ BSTR pResetTimingPhase2DIG,
    /* [in] */ BSTR pAdjPhase2DIG,
    /* [in] */ double pStartDatePhase2,
    /* [in] */ double pSpreadPhase2,
    /* [in] */ BSTR pDayCountPhase2,
    /* [in] */ BSTR pPayFreqPhase2,
    /* [in] */ BSTR pResetFreqPhase2,
    /* [in] */ BSTR pAdjPhase2,
    /* [in] */ double pStartDatePhase3,
    /* [in] */ double pEndDatePhase3,
    /* [in] */ BSTR pIndexPhase3,
    /* [in] */ double pSpreadPhase3,
    /* [in] */ BSTR pDayCountPhase3,
    /* [in] */ BSTR pPayFreqPhase3,
    /* [in] */ BSTR pResetFreqPhase3,
    /* [in] */ BSTR pResetTimingPhase3,
    /* [in] */ BSTR pAdjPhase3,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ VARIANT pSpreadFund,
    /* [in] */ BSTR pDayCountFund,
    /* [in] */ BSTR pPayFreqFund,
    /* [in] */ BSTR pResetFreqFund,
    /* [in] */ BSTR pResetTimingFund,
    /* [in] */ BSTR pAdjFund,
    /* [in] */ double pEndDateAmort,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    /* [in] */ BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ BSTR pFloorOrCap,
    /* [in] */ double pFee,
    /* [in] */ BSTR pVolCurvFromMatriceShift,
    /* [in] */ BSTR pVol,
    /* [in] */ BSTR pVolCub,
    /* [in] */ BSTR pCorrManager,
    /* [in] */ BSTR pConvexityManager,
    /* [in] */ BSTR pZc,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pSmiledModBump,
    /* [in] */ BSTR pHyperCubeCorrel,
    /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
    /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputePentifix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_ReplicConvAdjust_Create_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR Payoff_ReplicMode,
    /* [in] */ double Payoff_StepOrReplicPrecision,
    /* [in] */ BSTR Payoff_StopMode,
    /* [in] */ double Payoff_StopThreshold,
    /* [in] */ BSTR Sensi_ReplicMode,
    /* [in] */ double Sensi_StepOrReplicPrecision,
    /* [in] */ BSTR Sensi_StopMode,
    /* [in] */ double Sensi_StopThreshold,
    /* [in] */ BSTR UsedModelId,
    /* [in] */ double StrikeMinReplic,
    /* [in] */ double StrikeMaxReplic,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_ReplicConvAdjust_Create_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_MapConvAdjust_Create_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR LiborArrearAdj,
    /* [in] */ BSTR NaturalCMSAdj,
    /* [in] */ BSTR PaymentLagAdj,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_MapConvAdjust_Create_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetFees_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR securityId,
    /* [in] */ BSTR RefvalueId);


void __RPC_STUB IARMModule_ARM_Credit_SetFees_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetBounds_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR securityId,
    /* [out] */ double __RPC_FAR *down,
    /* [out] */ double __RPC_FAR *up);


void __RPC_STUB IARMModule_ARM_Credit_GetBounds_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputePentibor_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ double pStartdatePhase1,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexPay,
    /* [in] */ BSTR pIndexPhase1,
    /* [in] */ double pSpreadPhase1,
    /* [in] */ BSTR pDayCountPhase1,
    /* [in] */ BSTR pPayFreqPhase1,
    /* [in] */ BSTR pResetFreqPhase1,
    /* [in] */ BSTR pResetTimingPhase1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pAdjPhase1,
    /* [in] */ BSTR pStub,
    /* [in] */ BSTR pIndexPhase2DIG,
    /* [in] */ BSTR pIndexLongPhase2DIG,
    /* [in] */ BSTR pStrikePhase2DIG,
    /* [in] */ BSTR pResetTimingPhase2DIG,
    /* [in] */ BSTR pAdjPhase2DIG,
    /* [in] */ double pStartDatePhase2,
    /* [in] */ double pSpreadPhase2,
    /* [in] */ BSTR pDayCountPhase2,
    /* [in] */ BSTR pPayFreqPhase2,
    /* [in] */ BSTR pResetFreqPhase2,
    /* [in] */ BSTR pAdjPhase2,
    /* [in] */ double pStartDatePhase3,
    /* [in] */ double pEndDatePhase3,
    /* [in] */ BSTR pIndexPhase3,
    /* [in] */ double pSpreadPhase3,
    /* [in] */ BSTR pDayCountPhase3,
    /* [in] */ BSTR pPayFreqPhase3,
    /* [in] */ BSTR pResetFreqPhase3,
    /* [in] */ BSTR pResetTimingPhase3,
    /* [in] */ BSTR pAdjPhase3,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ VARIANT pSpreadFund,
    /* [in] */ BSTR pDayCountFund,
    /* [in] */ BSTR pPayFreqFund,
    /* [in] */ BSTR pResetFreqFund,
    /* [in] */ BSTR pResetTimingFund,
    /* [in] */ BSTR pAdjFund,
    /* [in] */ double pEndDateAmort,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    /* [in] */ BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ double pFee,
    /* [in] */ BSTR pVolCurvFromMatriceShift,
    /* [in] */ BSTR pVol,
    /* [in] */ BSTR pVolCub,
    /* [in] */ BSTR pConvexityManager,
    /* [in] */ BSTR pZc,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pSmiledModBump,
    /* [in] */ BSTR pHyperCubeCorrel,
    /* [in] */ BSTR pIndexIndexCorrelCube,
    /* [in] */ BSTR pCorrEUR,
    /* [in] */ BSTR pInterCorr,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputePentibor_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetRecovCoef_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurveId,
    /* [in] */ double RecovCoef,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetRecovCoef_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMIndexIndexCorrelCube_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pVolCurvId,
    /* [in] */ VARIANT __RPC_FAR *pTenors1List,
    /* [in] */ VARIANT __RPC_FAR *pTenors2List,
    /* [defaultvalue][in][optional] */ BSTR pInterSurfInterp,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMIndexIndexCorrelCube_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetInterpolationType_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pVolCurveId,
    /* [in] */ VARIANT __RPC_FAR *pInterpolType,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetInterpolationType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Customized_CDO_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pNotionals,
    /* [defaultvalue][in] */ BSTR Currency,
    /* [in] */ BSTR pDefaultLeg,
    /* [in] */ BSTR pPremiumLeg,
    /* [in] */ BSTR pParameters,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Customized_CDO_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCreateGenCorrelatorManager_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pMktTags,
    /* [in] */ VARIANT __RPC_FAR *pHyperDiagVol,
    /* [in] */ VARIANT __RPC_FAR *pIndexIndexVol,
    /* [in] */ VARIANT __RPC_FAR *pCorrelVol,
    /* [in] */ VARIANT __RPC_FAR *pIndexVol,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCreateGenCorrelatorManager_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CLN_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double EffectiveDate,
    /* [in] */ double EndDate,
    /* [in] */ double Spread,
    /* [defaultvalue][in][optional] */ BSTR IndexId,
    /* [defaultvalue][in][optional] */ double Refdate,
    /* [defaultvalue][in][optional] */ double pFstCpnEffDate,
    /* [defaultvalue][in][optional] */ double Notional,
    /* [defaultvalue][in][optional] */ BSTR AccOnDef,
    /* [defaultvalue][in][optional] */ BSTR DayCount,
    /* [defaultvalue][in][optional] */ BSTR DecompFreq,
    /* [defaultvalue][in][optional] */ BSTR StubRule,
    /* [defaultvalue][in][optional] */ double resetgap,
    /* [defaultvalue][in][optional] */ BSTR Currency,
    /* [defaultvalue][in][optional] */ BSTR ResetCal,
    /* [defaultvalue][in][optional] */ BSTR PayCal,
    /* [defaultvalue][in][optional] */ BSTR Nxchange,
    /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
    /* [defaultvalue][in][optional] */ BSTR AdjustedStartDate,
    /* [defaultvalue][in][optional] */ double Binary,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CLN_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputePentilix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ double pStartdatePhase1,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexPay,
    /* [in] */ BSTR pIndexPhase1,
    /* [in] */ double pSpreadPhase1,
    /* [in] */ BSTR pDayCountPhase1,
    /* [in] */ BSTR pPayFreqPhase1,
    /* [in] */ BSTR pResetFreqPhase1,
    /* [in] */ BSTR pResetTimingPhase1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pAdjPhase1,
    /* [in] */ BSTR pStub,
    /* [in] */ BSTR pIndexPhase2DIG,
    /* [in] */ BSTR pIndexLongPhase2DIG,
    /* [in] */ BSTR pStrikePhase2DIG,
    /* [in] */ BSTR pResetTimingPhase2DIG,
    /* [in] */ BSTR pAdjPhase2DIG,
    /* [in] */ double pStartDatePhase2,
    /* [in] */ double pSpreadPhase2,
    /* [in] */ BSTR pDayCountPhase2,
    /* [in] */ BSTR pPayFreqPhase2,
    /* [in] */ BSTR pResetFreqPhase2,
    /* [in] */ BSTR pAdjPhase2,
    /* [in] */ double pStartDatePhase3,
    /* [in] */ double pEndDatePhase3,
    /* [in] */ BSTR pIndexPhase3,
    /* [in] */ double pSpreadPhase3,
    /* [in] */ BSTR pDayCountPhase3,
    /* [in] */ BSTR pPayFreqPhase3,
    /* [in] */ BSTR pResetFreqPhase3,
    /* [in] */ BSTR pResetTimingPhase3,
    /* [in] */ BSTR pAdjPhase3,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ VARIANT pSpreadFund,
    /* [in] */ BSTR pDayCountFund,
    /* [in] */ BSTR pPayFreqFund,
    /* [in] */ BSTR pResetFreqFund,
    /* [in] */ BSTR pResetTimingFund,
    /* [in] */ BSTR pAdjFund,
    /* [in] */ double pEndDateAmort,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    /* [in] */ BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ double pFee,
    /* [in] */ BSTR pVolCurvFromMatriceShift,
    /* [in] */ BSTR pVol,
    /* [in] */ BSTR pVolCub,
    /* [in] */ BSTR pConvexityManager,
    /* [in] */ BSTR pZc,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pSmiledModBump,
    /* [in] */ BSTR pHyperCubeCorrel,
    /* [in] */ BSTR pIndexIndexCorrelCube,
    /* [in] */ BSTR pCorrEUR,
    /* [in] */ BSTR pInterCorr,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputePentilix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_RiskyPV01_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDefCurve,
    /* [optional][in] */ VARIANT Date1,
    /* [in] */ VARIANT Date2,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_RiskyPV01_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetResetMgrFromSummit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pSource,
    /* [defaultvalue][in] */ BSTR pCcy,
    /* [defaultvalue][in] */ BSTR pIsInflationIndex,
    /* [defaultvalue][in] */ BSTR pTerm,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetResetMgrFromSummit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetReset_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pResetMgr,
    /* [in] */ double pDate,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetReset_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSetLastFixing_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSecurityId,
    /* [in] */ double pRate,
    /* [in] */ double pAsOf,
    /* [defaultvalue][in][optional] */ double pBeforeLastFixingDate,
    /* [defaultvalue][in][optional] */ double pResetDate,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSetLastFixing_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Sectorial_Correlation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ DATE AsOf,
    /* [in] */ BSTR structName,
    /* [in] */ BSTR correlation_Type,
    /* [in] */ VARIANT __RPC_FAR *vLabels,
    /* [in] */ VARIANT __RPC_FAR *vector_Membership,
    /* [defaultvalue][in][optional] */ double intra_Sector_Correlation,
    /* [defaultvalue][in][optional] */ double inter_Sector_Correlation,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vBetas_Down,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *vLambdas_Down,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Sectorial_Correlation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeReviPentix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ VARIANT pDate,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexPhase1,
    /* [in] */ VARIANT pSpread,
    /* [in] */ BSTR pDayCountPhase1,
    /* [in] */ BSTR pPayFreqPhase1,
    /* [in] */ BSTR pResetFreqPhase1,
    /* [in] */ BSTR pResetTimingPhase1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pAdjPhase1,
    /* [in] */ BSTR pStub,
    /* [in] */ BSTR pIndexPhase2DIG,
    /* [in] */ BSTR pIndexLongPhase2DIG,
    /* [in] */ BSTR pStrikePhase2DIG,
    /* [in] */ BSTR pResetTimingPhase2DIG,
    /* [in] */ BSTR pAdjPhase2DIG,
    /* [in] */ BSTR pDayCountPhase2,
    /* [in] */ BSTR pPayFreqPhase2,
    /* [in] */ BSTR pResetFreqPhase2,
    /* [in] */ BSTR pAdjPhase2,
    /* [in] */ BSTR pIndexPhase3,
    /* [in] */ BSTR pDayCountPhase3,
    /* [in] */ BSTR pPayFreqPhase3,
    /* [in] */ BSTR pResetFreqPhase3,
    /* [in] */ BSTR pResetTimingPhase3,
    /* [in] */ BSTR pAdjPhase3,
    /* [in] */ BSTR pIndexFund,
    /* [in] */ VARIANT pSpreadFund,
    /* [in] */ BSTR pDayCountFund,
    /* [in] */ BSTR pPayFreqFund,
    /* [in] */ BSTR pResetFreqFund,
    /* [in] */ BSTR pResetTimingFund,
    /* [in] */ BSTR pAdjFund,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    /* [in] */ BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ BSTR pFloorOrCap,
    /* [in] */ double pFee,
    /* [in] */ double pLevier,
    /* [in] */ double pTxFixeMax,
    /* [in] */ BSTR pIsCapped,
    /* [in] */ double pTxCap,
    /* [in] */ BSTR pVolCurvFromMatriceShift,
    /* [in] */ BSTR pVol,
    /* [in] */ BSTR pVolCub,
    /* [in] */ BSTR pCorrManager,
    /* [in] */ BSTR pConvexityManager,
    /* [in] */ BSTR pZc,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pHyperCubeCorrel,
    /* [in] */ VARIANT __RPC_FAR *pBumpBsGenMod,
    /* [in] */ VARIANT __RPC_FAR *pBumpVolBsGenMod,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeReviPentix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGenAmortization_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pSwaplegId,
    /* [in] */ BSTR pAmortMethod,
    /* [defaultvalue][in][optional] */ BSTR pAmortFreq,
    /* [defaultvalue][in][optional] */ double pAmortAmount,
    /* [defaultvalue][in][optional] */ BSTR pDaycount,
    /* [defaultvalue][in][optional] */ double pLegNotional,
    /* [defaultvalue][in][optional] */ double pAmortRate,
    /* [defaultvalue][in][optional] */ double pReducedMaturity,
    /* [defaultvalue][in][optional] */ BSTR pModelId,
    /* [defaultvalue][in][optional] */ double pCleanUp,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGenAmortization_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCptRefvalue_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pRefValId,
    /* [in] */ double pDate,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMCptRefvalue_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetDPFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ BSTR issuer,
    /* [in] */ BSTR seniority,
    /* [in] */ BSTR ccy,
    /* [defaultvalue][in][optional] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFile,
    /* [defaultvalue][in][optional] */ BSTR irCurveId,
    /* [defaultvalue][in][optional] */ BSTR label,
    /* [retval][out] */ VARIANT __RPC_FAR *ret);


void __RPC_STUB IARMModule_ARM_Credit_GetDPFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCalypsoDevConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMCalypsoDevConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetZCFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pTerm,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR pInterpMethod,
    /* [defaultvalue][in][optional] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFile,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetZCFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DefProbInverse_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurveId,
    /* [in] */ double dDefProba,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DefProbInverse_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DPMktDataFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ BSTR issuer,
    /* [in] */ BSTR seniority,
    /* [in] */ BSTR ccy,
    /* [defaultvalue][in][optional] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFile,
    /* [in] */ VARIANT __RPC_FAR *Parameter,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DPMktDataFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pDate,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ BSTR issuer,
    /* [in] */ BSTR seniority,
    /* [in] */ BSTR ccy,
    /* [defaultvalue][in][optional] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFile,
    /* [defaultvalue][in][optional] */ BSTR irCurveId,
    /* [defaultvalue][in][optional] */ BSTR label,
    /* [retval][out] */ VARIANT __RPC_FAR *ret);


void __RPC_STUB IARMModule_ARM_Credit_ZeroCouponDefaultCurveFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetInitialCurveFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pIndex,
    /* [in] */ BSTR pCurrency,
    /* [in] */ BSTR pTerm,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFile,
    /* [defaultvalue][in][optional] */ BSTR pDoAdj,
    /* [out] */ VARIANT __RPC_FAR *pRetMat,
    /* [out] */ VARIANT __RPC_FAR *pRetRate);


void __RPC_STUB IARMModule_ARMGetInitialCurveFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCalypsoProdConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMCalypsoProdConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMCalypsoRecConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_ARMCalypsoRecConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetBasketCorrelMkDataFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ double date,
    /* [in] */ BSTR forceCurveName,
    /* [defaultvalue][in][optional] */ BSTR xmlFileName,
    /* [out] */ VARIANT __RPC_FAR *pRetMat,
    /* [out] */ VARIANT __RPC_FAR *pRetTenor,
    /* [out] */ VARIANT __RPC_FAR *pRetVol);


void __RPC_STUB IARMModule_ARM_Credit_GetBasketCorrelMkDataFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_QMatrix_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pQMatrix,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_QMatrix_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_MarketDataMng_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pstrVect,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_MarketDataMng_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ModelMultiCvMktDataMng_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pIRcurve,
    /* [in] */ VARIANT __RPC_FAR *pDefCurves,
    /* [in] */ VARIANT __RPC_FAR *pRecovery,
    /* [in] */ BSTR CorrelId,
    /* [in] */ BSTR MktdataMngId,
    /* [defaultvalue][in][optional] */ BSTR pVolcurve,
    /* [defaultvalue][in][optional] */ BSTR cloneorNot,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_ModelMultiCvMktDataMng_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_Local_ARM_ProdConnect_Proxy( 
    IARMModule __RPC_FAR * This);


void __RPC_STUB IARMModule_Local_ARM_ProdConnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_SetDefaultCurrency_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR isoCCy,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_SetDefaultCurrency_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLivretALeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pRcvOrPay,
    /* [defaultvalue][in][optional] */ VARIANT pSpread,
    /* [defaultvalue][in][optional] */ BSTR pResetFreq,
    /* [defaultvalue][in][optional] */ BSTR pPayFreq,
    /* [defaultvalue][in][optional] */ BSTR pResetTiming,
    /* [defaultvalue][in][optional] */ BSTR pPayTiming,
    /* [defaultvalue][in][optional] */ BSTR pCcy,
    /* [defaultvalue][in][optional] */ BSTR pIntRule,
    /* [defaultvalue][in][optional] */ double pResetGap,
    /* [defaultvalue][in][optional] */ BSTR pResetCal,
    /* [defaultvalue][in][optional] */ BSTR pPayCal,
    /* [defaultvalue][in][optional] */ double pDecompPricingFlag,
    /* [defaultvalue][in][optional] */ BSTR pNxChange,
    /* [defaultvalue][in][optional] */ BSTR pStubRule,
    /* [defaultvalue][in][optional] */ double pRefDate,
    /* [defaultvalue][in][optional] */ BSTR pAdjStartDate,
    /* [defaultvalue][in][optional] */ BSTR pDayCount,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLivretALeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CorrelationSmileStrike_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [in] */ VARIANT __RPC_FAR *pVolCurves,
    /* [in] */ VARIANT __RPC_FAR *pProportions,
    /* [in] */ double AsOf,
    /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeLow,
    /* [in][optional] */ VARIANT __RPC_FAR *pSmileStrikeHigh,
    /* [in][optional] */ VARIANT __RPC_FAR *pIndexVector,
    /* [defaultvalue][in][optional] */ BSTR Name,
    /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeLow,
    /* [in][optional] */ VARIANT __RPC_FAR *pFullStrikeUp,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CorrelationSmileStrike_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMFutDelivery_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pFut,
    /* [in] */ BSTR pCcy,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMFutDelivery_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMLivretACurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pAsOf,
    /* [in] */ BSTR pInfCurvId,
    /* [in] */ BSTR pEuribCurvId,
    /* [defaultvalue][in][optional] */ double pFlagRouding,
    /* [defaultvalue][in][optional] */ BSTR pInfResetMgrId,
    /* [defaultvalue][in][optional] */ BSTR pFixingLivretAId,
    /* [defaultvalue][in][optional] */ BSTR pFixingEuribId,
    /* [defaultvalue][in][optional] */ BSTR pMonthForAugust,
    /* [defaultvalue][in][optional] */ BSTR pMonthForFebruary,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMLivretACurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeLivretA_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ double pStartDateLeg1,
    /* [in] */ double pEndDateLeg1,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexLeg1,
    /* [in] */ VARIANT pSpreadLeg1,
    /* [in] */ BSTR pDayCountLeg1,
    /* [in] */ BSTR pPayFreqLeg1,
    /* [in] */ BSTR pResetFreqLeg1,
    /* [in] */ BSTR pResetTimingLeg1,
    /* [in] */ BSTR pAdjLeg1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pStub,
    /* [in] */ double pEndDateLA,
    /* [in] */ double pSpreadLeg2,
    /* [in] */ BSTR pDayCountLA,
    /* [in] */ BSTR pPayFreqLA,
    /* [in] */ BSTR pResetFreqLA,
    /* [in] */ BSTR pResetTimingLA,
    /* [in] */ BSTR pAdjLA,
    /* [in] */ BSTR pIndexLeg2,
    /* [in] */ BSTR pDayCountLeg2,
    /* [in] */ BSTR pPayFreqLeg2,
    /* [in] */ BSTR pResetFreqLeg2,
    /* [in] */ BSTR pResetTimingLeg2,
    /* [in] */ BSTR pAdjLeg2,
    /* [in] */ double pEndDateAmort,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ double pFee,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pSmiledModBump,
    /* [in] */ BSTR pLAMod,
    /* [in] */ BSTR pLAModBump,
    /* [in] */ BSTR pLAModBumpInflation,
    /* [in] */ VARIANT __RPC_FAR *pResetMgrIds,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeLivretA_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMGetFixingFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR source,
    /* [in] */ BSTR index,
    /* [in] */ BSTR term,
    /* [in] */ BSTR ccy,
    /* [in] */ BSTR curveName,
    /* [in] */ DATE date,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMGetFixingFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_FxConvertFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR ccy1,
    /* [in] */ BSTR ccy2,
    /* [defaultvalue][in][optional] */ BSTR pCvName,
    /* [in] */ DATE pDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_FxConvertFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FwdSpread_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR defcurveId,
    /* [in] */ double Maturity1,
    /* [in] */ double Maturity2,
    /* [defaultvalue][in][optional] */ double FwdStartDate,
    /* [defaultvalue][in][optional] */ double FwdEndDate,
    /* [defaultvalue][in][optional] */ BSTR VolId,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_FwdSpread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Flat_Correlation_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ DATE AsOf,
    /* [in] */ BSTR structName,
    /* [in] */ double correlValue,
    /* [defaultvalue][in] */ BSTR idIndex1,
    /* [defaultvalue][in] */ BSTR idIndex2,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Flat_Correlation_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CorridorLeg_Sche_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double Notional,
    /* [in] */ BSTR RecieveOrPay,
    /* [in] */ BSTR RefValueSpreadsInBP,
    /* [in] */ BSTR floatingIdx,
    /* [in] */ double leverageFloatIdx,
    /* [in] */ BSTR creditIdx,
    /* [in] */ BSTR refvalueKUPinBP,
    /* [in] */ BSTR refvalueKDWinBP,
    /* [in] */ BSTR ScheduleInfoId,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *accondef,
    /* [defaultvalue][in][optional] */ BSTR disc_ccy,
    /* [defaultvalue][in][optional] */ BSTR Name,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CorridorLeg_Sche_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Schedule_Info_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double EffectiveDate,
    /* [in] */ double MaturityDate,
    /* [defaultvalue][in][optional] */ BSTR payFrequency,
    /* [defaultvalue][in][optional] */ BSTR ResetFreq,
    /* [defaultvalue][in][optional] */ BSTR DayCount,
    /* [defaultvalue][in][optional] */ BSTR Stubrule,
    /* [defaultvalue][in][optional] */ BSTR intRule,
    /* [defaultvalue][in][optional] */ BSTR payCalName,
    /* [defaultvalue][in][optional] */ BSTR PayTiming,
    /* [defaultvalue][in][optional] */ BSTR ResetTiming,
    /* [defaultvalue][in][optional] */ BSTR fwdRule,
    /* [defaultvalue][in][optional] */ BSTR IncludeMaturity,
    /* [defaultvalue][in][optional] */ BSTR adj,
    /* [defaultvalue][in][optional] */ BSTR intStartAdj,
    /* [defaultvalue][in][optional] */ BSTR AccDayCount,
    /* [defaultvalue][in][optional] */ double ReferenceDate,
    /* [defaultvalue][in][optional] */ double FirstCpnEffDate,
    /* [defaultvalue][in][optional] */ BSTR AdjCal,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Schedule_Info_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMInfCurveSetResetMgr_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pInfCurve,
    /* [in] */ BSTR pResetMgr,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMInfCurveSetResetMgr_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CPDO_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pRiskyLeg,
    /* [in] */ BSTR pRollLeg,
    /* [in] */ BSTR pNoRiskyLeg,
    /* [in] */ double pInitialValo,
    /* [in] */ double pTarget,
    /* [in] */ double pMaturity,
    /* [in] */ BSTR pCpnType,
    /* [in] */ double pUFFees,
    /* [in] */ double pRunningFees,
    /* [in] */ double pVExpo,
    /* [in] */ double pV0Expo,
    /* [in] */ double pAlpha,
    /* [in] */ double pBeta,
    /* [in] */ double pDesactivation,
    /* [in] */ int pNbAssets,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CPDO_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_PriceVector_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ BSTR pCPTTYPE,
    /* [out] */ VARIANT __RPC_FAR *pRetVectorValos);


void __RPC_STUB IARMModule_ARM_Credit_PriceVector_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenPrice_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pPricer,
    /* [in] */ BSTR pCPTTYPE,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenPrice_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeTxFixed_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pNtl,
    /* [in] */ double pStartDateLeg1,
    /* [in] */ double pEndDateLeg1,
    /* [in] */ BSTR pCcy,
    /* [in] */ BSTR pIndexLeg1,
    /* [in] */ VARIANT pSpreadLeg1,
    /* [in] */ BSTR pDayCountLeg1,
    /* [in] */ BSTR pPayFreqLeg1,
    /* [in] */ BSTR pResetFreqLeg1,
    /* [in] */ BSTR pResetTimingLeg1,
    /* [in] */ BSTR pAdjLeg1,
    /* [in] */ BSTR pRoll,
    /* [in] */ BSTR pStub,
    /* [in] */ double pEndDateFixed,
    /* [in] */ double pSpreadLeg2,
    /* [in] */ BSTR pDayCountFixed,
    /* [in] */ BSTR pPayFreqFixed,
    /* [in] */ BSTR pResetFreqFixed,
    /* [in] */ BSTR pResetTimingFixed,
    /* [in] */ BSTR pAdjFixed,
    /* [in] */ BSTR pIndexLeg2,
    /* [in] */ BSTR pDayCountLeg2,
    /* [in] */ BSTR pPayFreqLeg2,
    /* [in] */ BSTR pResetFreqLeg2,
    /* [in] */ BSTR pResetTimingLeg2,
    /* [in] */ BSTR pAdjLeg2,
    /* [in] */ double pEndDateAmort,
    /* [in] */ BSTR pDayCountAmort,
    /* [in] */ BSTR pIntRuleAmort,
    /* [in] */ double pTxAmort,
    BSTR pFreqAmort,
    /* [in] */ double pAmountAmort,
    /* [in] */ BSTR pTypeAmort,
    /* [in] */ double pFee,
    /* [in] */ BSTR pSmiledMod,
    /* [in] */ BSTR pSmiledModBump,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeTxFixed_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FixingCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pDates,
    /* [in] */ VARIANT __RPC_FAR *pValues,
    /* [defaultvalue][in][optional] */ double AsOfDate,
    /* [defaultvalue][in][optional] */ BSTR B_IndexName,
    /* [defaultvalue][in][optional] */ BSTR B_IndexID,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_FixingCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CptInterpolDefCurveOLD_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurve,
    /* [in] */ BSTR pTenor,
    /* [in] */ double pSlope,
    /* [in] */ double pDate,
    /* [defaultvalue][in][optional] */ double pInterpDate,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CptInterpolDefCurveOLD_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_CreateBasketCorrelMkDataFromCalypso_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricingEnv,
    /* [in] */ double date,
    /* [in] */ BSTR forceCurveName,
    /* [in] */ BSTR Ccy,
    /* [defaultvalue][in][optional] */ BSTR xmlFileName,
    /* [defaultvalue][in][optional] */ BSTR indexId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_CreateBasketCorrelMkDataFromCalypso_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSetCalendar_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pFileName,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSetCalendar_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMInitGigaSpaces_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pUrl,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMInitGigaSpaces_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GetExpectedLoss_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pricerId,
    /* [in] */ double YearTerm,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GetExpectedLoss_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_VariableCollateral_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT pLabels,
    /* [in] */ VARIANT pNotionals,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_VariableCollateral_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_IndexCompo_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR IndexName,
    /* [in] */ VARIANT __RPC_FAR *pLabels,
    /* [defaultvalue][in] */ VARIANT __RPC_FAR *YearFrac,
    /* [defaultvalue][in] */ VARIANT __RPC_FAR *pSpread,
    /* [defaultvalue][in] */ BSTR Method,
    /* [defaultvalue][in] */ BSTR Basis,
    /* [defaultvalue][in] */ BSTR ResetFreq,
    /* [defaultvalue][in] */ BSTR PayFreq,
    /* [defaultvalue][in] */ BSTR ccy,
    /* [defaultvalue][in] */ BSTR fwdRule,
    /* [defaultvalue][in] */ BSTR resetTiming,
    /* [defaultvalue][in] */ int resetGap,
    /* [defaultvalue][in] */ BSTR payTiming,
    /* [defaultvalue][in] */ int payGap,
    /* [defaultvalue][in] */ BSTR intRule,
    /* [defaultvalue][in] */ BSTR AdjCalType,
    /* [defaultvalue][in] */ int cm_resetWeekDay,
    /* [defaultvalue][in] */ int cm_resetOccur,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_IndexCompo_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_createFlatCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurve,
    /* [in] */ VARIANT pTenor,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_createFlatCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FunctionRegister_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ long address);


void __RPC_STUB IARMModule_ARM_Credit_FunctionRegister_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBermudanXStyle_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pxDates,
    /* [defaultvalue][in][optional] */ VARIANT __RPC_FAR *pexpiryDates,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBermudanXStyle_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMcomputeCRA_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pFixorFloat,
    /* [in] */ double pFee,
    /* [in] */ double pAsOf,
    /* [in] */ double pStartDate,
    /* [in] */ double pEndDate,
    /* [in] */ BSTR pCcy,
    /* [in] */ double pLevelUp,
    /* [in] */ BSTR pUpSpec,
    /* [in] */ double pLevelDown,
    /* [in] */ BSTR pDownSpec,
    /* [in] */ BSTR pRefIndex,
    /* [in] */ BSTR pDayCount,
    /* [in] */ BSTR pPayFreqPayIndex,
    /* [in] */ BSTR pResetFreqRefIndex,
    /* [in] */ BSTR pPaidRstTiming,
    /* [in] */ BSTR pRefRstTiming,
    /* [in] */ BSTR pStubRule,
    /* [in] */ BSTR pPOrR,
    /* [in] */ double pStartCallDate,
    /* [in] */ BSTR pXStyle,
    /* [in] */ BSTR pFundingIndex,
    /* [in] */ BSTR pResetFreqFunding,
    /* [in] */ BSTR pPayFreqFunding,
    /* [in] */ VARIANT pSpreadFunding,
    /* [in] */ BSTR pPOrRFunding,
    /* [in] */ double pDecompPricingFlag,
    /* [in] */ double pdiscMarginFactor,
    /* [in] */ BSTR pPreInitFlag,
    /* [in] */ double pMeanReversion,
    /* [in] */ VARIANT pCalibParams,
    /* [in] */ VARIANT pCalibParamsPF,
    /* [in] */ double pKernelToGP,
    /* [in] */ VARIANT pMarkovTreeParams,
    /* [in] */ double pMarkovTreePathNumber,
    /* [in] */ BSTR pBsmodId,
    /* [in] */ BSTR pBsmodSwoptId,
    /* [in] */ BSTR pBsmodSwoptBumpId,
    /* [in] */ BSTR pzcId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMcomputeCRA_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Math_BivNormale_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double x,
    /* [in] */ double y,
    /* [in] */ double rho,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Math_BivNormale_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_GETINSTRUMENTFROMCALYPSO_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR CalypsoId,
    /* [in] */ BSTR Type,
    /* [defaultvalue][optional][in] */ double AsOf,
    /* [defaultvalue][optional][in] */ BSTR ModelType,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_GETINSTRUMENTFROMCALYPSO_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMSetDiscountPricingMode_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pModelId,
    /* [in] */ int pDiscountPricingMode,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMSetDiscountPricingMode_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Math_RandUniform_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [defaultvalue][optional][in] */ double seed,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Math_RandUniform_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Math_Interpol_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *X,
    /* [in] */ VARIANT __RPC_FAR *Y,
    /* [in] */ double value,
    /* [defaultvalue][in][optional] */ double type,
    /* [defaultvalue][in][optional] */ double smooth,
    /* [in][optional] */ VARIANT __RPC_FAR *Weights,
    /* [defaultvalue][in][optional] */ double modeSpline,
    /* [defaultvalue][in][optional] */ double withC1condition,
    /* [defaultvalue][in][optional] */ double leftSlope,
    /* [defaultvalue][in][optional] */ double rightSlope,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Math_Interpol_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Random_Generator_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR RandomType,
    /* [defaultvalue][optional][in] */ BSTR ParamId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Random_Generator_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenerateOneRandom_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR RandomId,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenerateOneRandom_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_GenerateRandoms_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR RandomId,
    /* [in] */ int DimVector,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_GenerateRandoms_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_ResetRandom_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR RandomId);


void __RPC_STUB IARMModule_ARM_Credit_ResetRandom_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_FwdSpreadAsIndex_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR DefCurveId,
    /* [in] */ double matu1,
    /* [in] */ double Matu2,
    /* [retval][out] */ double __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_FwdSpreadAsIndex_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_createDefCurveFromBase_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pCurveCDS,
    /* [in] */ BSTR pCurveIndex,
    /* [in] */ VARIANT vBase,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_createDefCurveFromBase_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_RiskyPV01AsSensitivity_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pDefCurve,
    /* [in] */ BSTR Tenor,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_RiskyPV01AsSensitivity_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_SetVolCurve_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR Model,
    /* [in] */ BSTR VolCurveId,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_SetVolCurve_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMPF_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT __RPC_FAR *pinsts,
    /* [in] */ VARIANT __RPC_FAR *pcoeffs,
    /* [in] */ VARIANT __RPC_FAR *pmarketPrices,
    VARIANT __RPC_FAR *pprecisions,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMPF_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMBondTEC_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double pIssueDate,
    /* [in] */ double pMaturityDate,
    /* [in] */ double pFirstCpnDate,
    /* [in] */ double pCpnRate,
    /* [in] */ double pRedempPrice,
    /* [in] */ long pPeriodicity,
    /* [in] */ VARIANT pDaycount,
    /* [defaultvalue][in][optional] */ long pSettleGap,
    /* [defaultvalue][in][optional] */ long pCpnDateFlag,
    /* [defaultvalue][in][optional] */ BSTR pCcyId,
    /* [defaultvalue][in][optional] */ double ptec,
    /* [defaultvalue][in][optional] */ BSTR pPFTecId,
    /* [defaultvalue][in][optional] */ BSTR pModTecId,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMBondTEC_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMPFModFit_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR pmodName,
    /* [in] */ BSTR ppf,
    /* [in] */ double psettlement,
    /* [in] */ BSTR pzc,
    /* [in] */ VARIANT __RPC_FAR *pvList,
    /* [in] */ VARIANT __RPC_FAR *pfList,
    /* [defaultvalue][in][optional] */ long nag_algo,
    /* [defaultvalue][in][optional] */ long pstep,
    /* [defaultvalue][in][optional] */ double phorizon,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMPFModFit_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_Restrikable_CDO_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double UnderlyingMatu,
    /* [in] */ double Expiry,
    /* [in] */ double Strike,
    /* [in] */ int OptionType,
    /* [in] */ BSTR pUnderlying,
    /* [in] */ double Rehauss,
    /* [in] */ double TriggerFreq,
    /* [in] */ int DiffCDO,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_Restrikable_CDO_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_PropertyList_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ VARIANT attrNames,
    /* [in] */ VARIANT attrValues,
    /* [in][optional] */ VARIANT attrTypes,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_PropertyList_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARM_Credit_DefCurveIntensityPWC_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ double AsOfDate,
    /* [in] */ VARIANT __RPC_FAR *pMatuRates,
    /* [in] */ VARIANT __RPC_FAR *pInputs,
    /* [in] */ double Type,
    /* [in] */ double Recovery,
    /* [in] */ BSTR IRCurveId,
    /* [defaultvalue][in][optional] */ BSTR bCurrency,
    /* [defaultvalue][in][optional] */ BSTR bLabel,
    /* [defaultvalue][in][optional] */ BSTR VolCurveId,
    /* [defaultvalue][in][optional] */ BSTR calibrationAlgo,
    /* [defaultvalue][in][optional] */ int lag,
    /* [retval][out] */ VARIANT __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARM_Credit_DefCurveIntensityPWC_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IARMModule_ARMTMLeg_Proxy( 
    IARMModule __RPC_FAR * This,
    /* [in] */ BSTR ptmIxType,
    /* [in] */ double pstartDate,
    /* [in] */ double pendDate,
    /* [in] */ BSTR pPorR,
    /* [defaultvalue][in][optional] */ double pspread,
    /* [defaultvalue][in][optional] */ BSTR ppayFrequency,
    /* [defaultvalue][in][optional] */ BSTR presetFrequency,
    /* [defaultvalue][in][optional] */ BSTR pinterestRule,
    /* [defaultvalue][in][optional] */ BSTR pfwdRule,
    /* [defaultvalue][in][optional] */ BSTR pstubRule,
    /* [in] */ BSTR pccy,
    /* [retval][out] */ BSTR __RPC_FAR *pRet);


void __RPC_STUB IARMModule_ARMTMLeg_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IARMModule_INTERFACE_DEFINED__ */



#ifndef __LOCAL_DLLARMLib_LIBRARY_DEFINED__
#define __LOCAL_DLLARMLib_LIBRARY_DEFINED__

/* library LOCAL_DLLARMLib */
/* [helpstring][version][uuid] */ 


EXTERN_C const IID LIBID_LOCAL_DLLARMLib;

EXTERN_C const CLSID CLSID_ARMModule;

#ifdef __cplusplus

class DECLSPEC_UUID("47CC0BC7-4EDF-4853-838C-89B47E579AD4")
ARMModule;
#endif
#endif /* __LOCAL_DLLARMLib_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

unsigned long             __RPC_USER  BSTR_UserSize(     unsigned long __RPC_FAR *, unsigned long            , BSTR __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  BSTR_UserMarshal(  unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, BSTR __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  BSTR_UserUnmarshal(unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, BSTR __RPC_FAR * ); 
void                      __RPC_USER  BSTR_UserFree(     unsigned long __RPC_FAR *, BSTR __RPC_FAR * ); 

unsigned long             __RPC_USER  VARIANT_UserSize(     unsigned long __RPC_FAR *, unsigned long            , VARIANT __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  VARIANT_UserMarshal(  unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, VARIANT __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  VARIANT_UserUnmarshal(unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, VARIANT __RPC_FAR * ); 
void                      __RPC_USER  VARIANT_UserFree(     unsigned long __RPC_FAR *, VARIANT __RPC_FAR * ); 

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif
