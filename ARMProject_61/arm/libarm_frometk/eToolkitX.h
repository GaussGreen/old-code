

/* this ALWAYS GENERATED file contains the definitions for the interfaces */


 /* File created by MIDL compiler version 6.00.0361 */
/* at Wed Apr 20 16:24:40 2005
 */
/* Compiler settings for .\eToolkitX.idl:
    Oicf, W1, Zp8, env=Win32 (32b run)
    protocol : dce , ms_ext, c_ext, robust
    error checks: allocation ref bounds_check enum stub_data 
    VC __declspec() decoration level: 
         __declspec(uuid()), __declspec(selectany), __declspec(novtable)
         DECLSPEC_UUID(), MIDL_INTERFACE()
*/
//@@MIDL_FILE_HEADING(  )

#pragma warning( disable: 4049 )  /* more than 64k source lines */


/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 475
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

#ifndef __eToolkitX_h__
#define __eToolkitX_h__

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

/* Forward Declarations */ 

#ifndef __IeToolkit_FWD_DEFINED__
#define __IeToolkit_FWD_DEFINED__
typedef interface IeToolkit IeToolkit;
#endif 	/* __IeToolkit_FWD_DEFINED__ */


#ifndef __eToolkit_FWD_DEFINED__
#define __eToolkit_FWD_DEFINED__

#ifdef __cplusplus
typedef class eToolkit eToolkit;
#else
typedef struct eToolkit eToolkit;
#endif /* __cplusplus */

#endif 	/* __eToolkit_FWD_DEFINED__ */


/* header files for imported files */
#include "oaidl.h"
#include "ocidl.h"

#ifdef __cplusplus
extern "C"{
#endif 

void * __RPC_USER MIDL_user_allocate(size_t);
void __RPC_USER MIDL_user_free( void * ); 

#ifndef __IeToolkit_INTERFACE_DEFINED__
#define __IeToolkit_INTERFACE_DEFINED__

/* interface IeToolkit */
/* [unique][helpstring][dual][uuid][object] */ 


EXTERN_C const IID IID_IeToolkit;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("457080BA-F587-496D-AA87-E66B5C8852F8")
    IeToolkit : public IDispatch
    {
    public:
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE Connect( 
            VARIANT *userName,
            VARIANT *password,
            VARIANT *applicationName,
            VARIANT *programMode,
            VARIANT *databaseContext) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE Disconnect( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE Execute( 
            VARIANT *command,
            VARIANT *XMLRequest,
            VARIANT *XMLResponse,
            VARIANT *messageList) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetContextList( 
            /* [retval][out] */ VARIANT *result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetIsin( 
            VARIANT *SecId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Isin) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetBondDef( 
            VARIANT *SecId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *XMLResponse) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetCodeTiersFM( 
            VARIANT *Cust,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *XMLResponse) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetCommsetName( 
            VARIANT *Name,
            VARIANT *Name2,
            VARIANT *AsOfDate,
            VARIANT *Type,
            VARIANT *Id,
            VARIANT *XMLResponse,
            VARIANT *messageList) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetSecId( 
            VARIANT *Isin,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *SecId) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetNbConnections( 
            /* [retval][out] */ VARIANT *nb) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE ShutDown( void) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getSecIds( 
            VARIANT *Sec,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getCustomers( 
            VARIANT *custType,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE GetCurrentBase( 
            /* [retval][out] */ VARIANT *pDatabaseContext) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getCustType( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getTradingCustomers( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE tradeIdExiste( 
            VARIANT *tradeId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFixings( 
            VARIANT *startDate,
            VARIANT *endDate,
            VARIANT *Desk,
            VARIANT *Ccy,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getContracts( 
            VARIANT *Exchange,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getWriteDesk( 
            VARIANT *UserId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getBrokerCom( 
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getBookFils( 
            VARIANT *parentDesk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getWriteBook( 
            VARIANT *UserId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFrancoClearBroker( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getListExchange( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getDRule( 
            VARIANT *Contrat,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFavoriteContract( 
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFutDesk( 
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getPrefBrokerCom( 
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getBadLeadTrade( 
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getOptDesk( 
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getListExchangeOpt( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getPrefBrokerComOpt( 
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getContractsOpt( 
            VARIANT *Exchange,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getERule( 
            VARIANT *Contrat,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getBrokerComOpt( 
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFrancoClearBrokerOpt( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getFavoriteContractOpt( 
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getCustFromMMO( 
            VARIANT *CodeMMO,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getDefProbIssuerList( 
            VARIANT *CurveId,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getAppliqueClearBrokerOpt( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getAppliqueClearBroker( 
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
        virtual /* [helpstring][id] */ HRESULT STDMETHODCALLTYPE getREFRATE( 
            VARIANT *Source,
            VARIANT *Ccy,
            VARIANT *Index,
            VARIANT *Term,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IeToolkitVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE *QueryInterface )( 
            IeToolkit * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void **ppvObject);
        
        ULONG ( STDMETHODCALLTYPE *AddRef )( 
            IeToolkit * This);
        
        ULONG ( STDMETHODCALLTYPE *Release )( 
            IeToolkit * This);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfoCount )( 
            IeToolkit * This,
            /* [out] */ UINT *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetTypeInfo )( 
            IeToolkit * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo **ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE *GetIDsOfNames )( 
            IeToolkit * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR *rgszNames,
            /* [in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE *Invoke )( 
            IeToolkit * This,
            /* [in] */ DISPID dispIdMember,
            /* [in] */ REFIID riid,
            /* [in] */ LCID lcid,
            /* [in] */ WORD wFlags,
            /* [out][in] */ DISPPARAMS *pDispParams,
            /* [out] */ VARIANT *pVarResult,
            /* [out] */ EXCEPINFO *pExcepInfo,
            /* [out] */ UINT *puArgErr);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *Connect )( 
            IeToolkit * This,
            VARIANT *userName,
            VARIANT *password,
            VARIANT *applicationName,
            VARIANT *programMode,
            VARIANT *databaseContext);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *Disconnect )( 
            IeToolkit * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *Execute )( 
            IeToolkit * This,
            VARIANT *command,
            VARIANT *XMLRequest,
            VARIANT *XMLResponse,
            VARIANT *messageList);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetContextList )( 
            IeToolkit * This,
            /* [retval][out] */ VARIANT *result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetIsin )( 
            IeToolkit * This,
            VARIANT *SecId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Isin);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetBondDef )( 
            IeToolkit * This,
            VARIANT *SecId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *XMLResponse);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetCodeTiersFM )( 
            IeToolkit * This,
            VARIANT *Cust,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *XMLResponse);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetCommsetName )( 
            IeToolkit * This,
            VARIANT *Name,
            VARIANT *Name2,
            VARIANT *AsOfDate,
            VARIANT *Type,
            VARIANT *Id,
            VARIANT *XMLResponse,
            VARIANT *messageList);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetSecId )( 
            IeToolkit * This,
            VARIANT *Isin,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *SecId);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetNbConnections )( 
            IeToolkit * This,
            /* [retval][out] */ VARIANT *nb);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *ShutDown )( 
            IeToolkit * This);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getSecIds )( 
            IeToolkit * This,
            VARIANT *Sec,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getCustomers )( 
            IeToolkit * This,
            VARIANT *custType,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *GetCurrentBase )( 
            IeToolkit * This,
            /* [retval][out] */ VARIANT *pDatabaseContext);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getCustType )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getTradingCustomers )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *tradeIdExiste )( 
            IeToolkit * This,
            VARIANT *tradeId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFixings )( 
            IeToolkit * This,
            VARIANT *startDate,
            VARIANT *endDate,
            VARIANT *Desk,
            VARIANT *Ccy,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getContracts )( 
            IeToolkit * This,
            VARIANT *Exchange,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getWriteDesk )( 
            IeToolkit * This,
            VARIANT *UserId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getBrokerCom )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getBookFils )( 
            IeToolkit * This,
            VARIANT *parentDesk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getWriteBook )( 
            IeToolkit * This,
            VARIANT *UserId,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFrancoClearBroker )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getListExchange )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getDRule )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFavoriteContract )( 
            IeToolkit * This,
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFutDesk )( 
            IeToolkit * This,
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getPrefBrokerCom )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getBadLeadTrade )( 
            IeToolkit * This,
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getOptDesk )( 
            IeToolkit * This,
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *Desk,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getListExchangeOpt )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getPrefBrokerComOpt )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getContractsOpt )( 
            IeToolkit * This,
            VARIANT *Exchange,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getERule )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getBrokerComOpt )( 
            IeToolkit * This,
            VARIANT *Contrat,
            VARIANT *BrokerType,
            VARIANT *PreferedBrokers,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFrancoClearBrokerOpt )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getFavoriteContractOpt )( 
            IeToolkit * This,
            VARIANT *LoginUser,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getCustFromMMO )( 
            IeToolkit * This,
            VARIANT *CodeMMO,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getDefProbIssuerList )( 
            IeToolkit * This,
            VARIANT *CurveId,
            VARIANT *AsOfDate,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getAppliqueClearBrokerOpt )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getAppliqueClearBroker )( 
            IeToolkit * This,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        /* [helpstring][id] */ HRESULT ( STDMETHODCALLTYPE *getREFRATE )( 
            IeToolkit * This,
            VARIANT *Source,
            VARIANT *Ccy,
            VARIANT *Index,
            VARIANT *Term,
            VARIANT *messageList,
            /* [retval][out] */ VARIANT *Result);
        
        END_INTERFACE
    } IeToolkitVtbl;

    interface IeToolkit
    {
        CONST_VTBL struct IeToolkitVtbl *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IeToolkit_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IeToolkit_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IeToolkit_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IeToolkit_GetTypeInfoCount(This,pctinfo)	\
    (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo)

#define IeToolkit_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo)

#define IeToolkit_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)

#define IeToolkit_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)


#define IeToolkit_Connect(This,userName,password,applicationName,programMode,databaseContext)	\
    (This)->lpVtbl -> Connect(This,userName,password,applicationName,programMode,databaseContext)

#define IeToolkit_Disconnect(This)	\
    (This)->lpVtbl -> Disconnect(This)

#define IeToolkit_Execute(This,command,XMLRequest,XMLResponse,messageList)	\
    (This)->lpVtbl -> Execute(This,command,XMLRequest,XMLResponse,messageList)

#define IeToolkit_GetContextList(This,result)	\
    (This)->lpVtbl -> GetContextList(This,result)

#define IeToolkit_GetIsin(This,SecId,messageList,Isin)	\
    (This)->lpVtbl -> GetIsin(This,SecId,messageList,Isin)

#define IeToolkit_GetBondDef(This,SecId,messageList,XMLResponse)	\
    (This)->lpVtbl -> GetBondDef(This,SecId,messageList,XMLResponse)

#define IeToolkit_GetCodeTiersFM(This,Cust,messageList,XMLResponse)	\
    (This)->lpVtbl -> GetCodeTiersFM(This,Cust,messageList,XMLResponse)

#define IeToolkit_GetCommsetName(This,Name,Name2,AsOfDate,Type,Id,XMLResponse,messageList)	\
    (This)->lpVtbl -> GetCommsetName(This,Name,Name2,AsOfDate,Type,Id,XMLResponse,messageList)

#define IeToolkit_GetSecId(This,Isin,messageList,SecId)	\
    (This)->lpVtbl -> GetSecId(This,Isin,messageList,SecId)

#define IeToolkit_GetNbConnections(This,nb)	\
    (This)->lpVtbl -> GetNbConnections(This,nb)

#define IeToolkit_ShutDown(This)	\
    (This)->lpVtbl -> ShutDown(This)

#define IeToolkit_getSecIds(This,Sec,messageList,Result)	\
    (This)->lpVtbl -> getSecIds(This,Sec,messageList,Result)

#define IeToolkit_getCustomers(This,custType,messageList,Result)	\
    (This)->lpVtbl -> getCustomers(This,custType,messageList,Result)

#define IeToolkit_GetCurrentBase(This,pDatabaseContext)	\
    (This)->lpVtbl -> GetCurrentBase(This,pDatabaseContext)

#define IeToolkit_getCustType(This,messageList,Result)	\
    (This)->lpVtbl -> getCustType(This,messageList,Result)

#define IeToolkit_getTradingCustomers(This,messageList,Result)	\
    (This)->lpVtbl -> getTradingCustomers(This,messageList,Result)

#define IeToolkit_tradeIdExiste(This,tradeId,messageList,Result)	\
    (This)->lpVtbl -> tradeIdExiste(This,tradeId,messageList,Result)

#define IeToolkit_getFixings(This,startDate,endDate,Desk,Ccy,messageList,Result)	\
    (This)->lpVtbl -> getFixings(This,startDate,endDate,Desk,Ccy,messageList,Result)

#define IeToolkit_getContracts(This,Exchange,messageList,Result)	\
    (This)->lpVtbl -> getContracts(This,Exchange,messageList,Result)

#define IeToolkit_getWriteDesk(This,UserId,messageList,Result)	\
    (This)->lpVtbl -> getWriteDesk(This,UserId,messageList,Result)

#define IeToolkit_getBrokerCom(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)	\
    (This)->lpVtbl -> getBrokerCom(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)

#define IeToolkit_getBookFils(This,parentDesk,messageList,Result)	\
    (This)->lpVtbl -> getBookFils(This,parentDesk,messageList,Result)

#define IeToolkit_getWriteBook(This,UserId,messageList,Result)	\
    (This)->lpVtbl -> getWriteBook(This,UserId,messageList,Result)

#define IeToolkit_getFrancoClearBroker(This,messageList,Result)	\
    (This)->lpVtbl -> getFrancoClearBroker(This,messageList,Result)

#define IeToolkit_getListExchange(This,messageList,Result)	\
    (This)->lpVtbl -> getListExchange(This,messageList,Result)

#define IeToolkit_getDRule(This,Contrat,messageList,Result)	\
    (This)->lpVtbl -> getDRule(This,Contrat,messageList,Result)

#define IeToolkit_getFavoriteContract(This,LoginUser,AsOfDate,messageList,Result)	\
    (This)->lpVtbl -> getFavoriteContract(This,LoginUser,AsOfDate,messageList,Result)

#define IeToolkit_getFutDesk(This,LoginUser,AsOfDate,Desk,messageList,Result)	\
    (This)->lpVtbl -> getFutDesk(This,LoginUser,AsOfDate,Desk,messageList,Result)

#define IeToolkit_getPrefBrokerCom(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)	\
    (This)->lpVtbl -> getPrefBrokerCom(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)

#define IeToolkit_getBadLeadTrade(This,Desk,messageList,Result)	\
    (This)->lpVtbl -> getBadLeadTrade(This,Desk,messageList,Result)

#define IeToolkit_getOptDesk(This,LoginUser,AsOfDate,Desk,messageList,Result)	\
    (This)->lpVtbl -> getOptDesk(This,LoginUser,AsOfDate,Desk,messageList,Result)

#define IeToolkit_getListExchangeOpt(This,messageList,Result)	\
    (This)->lpVtbl -> getListExchangeOpt(This,messageList,Result)

#define IeToolkit_getPrefBrokerComOpt(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)	\
    (This)->lpVtbl -> getPrefBrokerComOpt(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)

#define IeToolkit_getContractsOpt(This,Exchange,messageList,Result)	\
    (This)->lpVtbl -> getContractsOpt(This,Exchange,messageList,Result)

#define IeToolkit_getERule(This,Contrat,messageList,Result)	\
    (This)->lpVtbl -> getERule(This,Contrat,messageList,Result)

#define IeToolkit_getBrokerComOpt(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)	\
    (This)->lpVtbl -> getBrokerComOpt(This,Contrat,BrokerType,PreferedBrokers,messageList,Result)

#define IeToolkit_getFrancoClearBrokerOpt(This,messageList,Result)	\
    (This)->lpVtbl -> getFrancoClearBrokerOpt(This,messageList,Result)

#define IeToolkit_getFavoriteContractOpt(This,LoginUser,AsOfDate,messageList,Result)	\
    (This)->lpVtbl -> getFavoriteContractOpt(This,LoginUser,AsOfDate,messageList,Result)

#define IeToolkit_getCustFromMMO(This,CodeMMO,messageList,Result)	\
    (This)->lpVtbl -> getCustFromMMO(This,CodeMMO,messageList,Result)

#define IeToolkit_getDefProbIssuerList(This,CurveId,AsOfDate,messageList,Result)	\
    (This)->lpVtbl -> getDefProbIssuerList(This,CurveId,AsOfDate,messageList,Result)

#define IeToolkit_getAppliqueClearBrokerOpt(This,messageList,Result)	\
    (This)->lpVtbl -> getAppliqueClearBrokerOpt(This,messageList,Result)

#define IeToolkit_getAppliqueClearBroker(This,messageList,Result)	\
    (This)->lpVtbl -> getAppliqueClearBroker(This,messageList,Result)

#define IeToolkit_getREFRATE(This,Source,Ccy,Index,Term,messageList,Result)	\
    (This)->lpVtbl -> getREFRATE(This,Source,Ccy,Index,Term,messageList,Result)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_Connect_Proxy( 
    IeToolkit * This,
    VARIANT *userName,
    VARIANT *password,
    VARIANT *applicationName,
    VARIANT *programMode,
    VARIANT *databaseContext);


void __RPC_STUB IeToolkit_Connect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_Disconnect_Proxy( 
    IeToolkit * This);


void __RPC_STUB IeToolkit_Disconnect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_Execute_Proxy( 
    IeToolkit * This,
    VARIANT *command,
    VARIANT *XMLRequest,
    VARIANT *XMLResponse,
    VARIANT *messageList);


void __RPC_STUB IeToolkit_Execute_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetContextList_Proxy( 
    IeToolkit * This,
    /* [retval][out] */ VARIANT *result);


void __RPC_STUB IeToolkit_GetContextList_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetIsin_Proxy( 
    IeToolkit * This,
    VARIANT *SecId,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Isin);


void __RPC_STUB IeToolkit_GetIsin_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetBondDef_Proxy( 
    IeToolkit * This,
    VARIANT *SecId,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *XMLResponse);


void __RPC_STUB IeToolkit_GetBondDef_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetCodeTiersFM_Proxy( 
    IeToolkit * This,
    VARIANT *Cust,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *XMLResponse);


void __RPC_STUB IeToolkit_GetCodeTiersFM_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetCommsetName_Proxy( 
    IeToolkit * This,
    VARIANT *Name,
    VARIANT *Name2,
    VARIANT *AsOfDate,
    VARIANT *Type,
    VARIANT *Id,
    VARIANT *XMLResponse,
    VARIANT *messageList);


void __RPC_STUB IeToolkit_GetCommsetName_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetSecId_Proxy( 
    IeToolkit * This,
    VARIANT *Isin,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *SecId);


void __RPC_STUB IeToolkit_GetSecId_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetNbConnections_Proxy( 
    IeToolkit * This,
    /* [retval][out] */ VARIANT *nb);


void __RPC_STUB IeToolkit_GetNbConnections_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_ShutDown_Proxy( 
    IeToolkit * This);


void __RPC_STUB IeToolkit_ShutDown_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getSecIds_Proxy( 
    IeToolkit * This,
    VARIANT *Sec,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getSecIds_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getCustomers_Proxy( 
    IeToolkit * This,
    VARIANT *custType,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getCustomers_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_GetCurrentBase_Proxy( 
    IeToolkit * This,
    /* [retval][out] */ VARIANT *pDatabaseContext);


void __RPC_STUB IeToolkit_GetCurrentBase_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getCustType_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getCustType_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getTradingCustomers_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getTradingCustomers_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_tradeIdExiste_Proxy( 
    IeToolkit * This,
    VARIANT *tradeId,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_tradeIdExiste_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFixings_Proxy( 
    IeToolkit * This,
    VARIANT *startDate,
    VARIANT *endDate,
    VARIANT *Desk,
    VARIANT *Ccy,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFixings_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getContracts_Proxy( 
    IeToolkit * This,
    VARIANT *Exchange,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getContracts_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getWriteDesk_Proxy( 
    IeToolkit * This,
    VARIANT *UserId,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getWriteDesk_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getBrokerCom_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *BrokerType,
    VARIANT *PreferedBrokers,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getBrokerCom_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getBookFils_Proxy( 
    IeToolkit * This,
    VARIANT *parentDesk,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getBookFils_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getWriteBook_Proxy( 
    IeToolkit * This,
    VARIANT *UserId,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getWriteBook_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFrancoClearBroker_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFrancoClearBroker_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getListExchange_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getListExchange_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getDRule_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getDRule_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFavoriteContract_Proxy( 
    IeToolkit * This,
    VARIANT *LoginUser,
    VARIANT *AsOfDate,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFavoriteContract_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFutDesk_Proxy( 
    IeToolkit * This,
    VARIANT *LoginUser,
    VARIANT *AsOfDate,
    VARIANT *Desk,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFutDesk_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getPrefBrokerCom_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *BrokerType,
    VARIANT *PreferedBrokers,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getPrefBrokerCom_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getBadLeadTrade_Proxy( 
    IeToolkit * This,
    VARIANT *Desk,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getBadLeadTrade_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getOptDesk_Proxy( 
    IeToolkit * This,
    VARIANT *LoginUser,
    VARIANT *AsOfDate,
    VARIANT *Desk,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getOptDesk_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getListExchangeOpt_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getListExchangeOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getPrefBrokerComOpt_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *BrokerType,
    VARIANT *PreferedBrokers,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getPrefBrokerComOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getContractsOpt_Proxy( 
    IeToolkit * This,
    VARIANT *Exchange,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getContractsOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getERule_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getERule_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getBrokerComOpt_Proxy( 
    IeToolkit * This,
    VARIANT *Contrat,
    VARIANT *BrokerType,
    VARIANT *PreferedBrokers,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getBrokerComOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFrancoClearBrokerOpt_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFrancoClearBrokerOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getFavoriteContractOpt_Proxy( 
    IeToolkit * This,
    VARIANT *LoginUser,
    VARIANT *AsOfDate,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getFavoriteContractOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getCustFromMMO_Proxy( 
    IeToolkit * This,
    VARIANT *CodeMMO,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getCustFromMMO_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getDefProbIssuerList_Proxy( 
    IeToolkit * This,
    VARIANT *CurveId,
    VARIANT *AsOfDate,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getDefProbIssuerList_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getAppliqueClearBrokerOpt_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getAppliqueClearBrokerOpt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getAppliqueClearBroker_Proxy( 
    IeToolkit * This,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getAppliqueClearBroker_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [helpstring][id] */ HRESULT STDMETHODCALLTYPE IeToolkit_getREFRATE_Proxy( 
    IeToolkit * This,
    VARIANT *Source,
    VARIANT *Ccy,
    VARIANT *Index,
    VARIANT *Term,
    VARIANT *messageList,
    /* [retval][out] */ VARIANT *Result);


void __RPC_STUB IeToolkit_getREFRATE_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IeToolkit_INTERFACE_DEFINED__ */



#ifndef __ETOOLKITXLib_LIBRARY_DEFINED__
#define __ETOOLKITXLib_LIBRARY_DEFINED__

/* library ETOOLKITXLib */
/* [helpstring][version][uuid] */ 


EXTERN_C const IID LIBID_ETOOLKITXLib;

EXTERN_C const CLSID CLSID_eToolkit;

#ifdef __cplusplus

class DECLSPEC_UUID("2273F4EF-F17C-4540-8B45-F3F34C23C8FF")
eToolkit;
#endif
#endif /* __ETOOLKITXLib_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

unsigned long             __RPC_USER  VARIANT_UserSize(     unsigned long *, unsigned long            , VARIANT * ); 
unsigned char * __RPC_USER  VARIANT_UserMarshal(  unsigned long *, unsigned char *, VARIANT * ); 
unsigned char * __RPC_USER  VARIANT_UserUnmarshal(unsigned long *, unsigned char *, VARIANT * ); 
void                      __RPC_USER  VARIANT_UserFree(     unsigned long *, VARIANT * ); 

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif


