/*-----------------------------------------------------------------------------*

     Global E-ToolKit header

 *-----------------------------------------------------------------------------*/
#ifndef ARM_LOCAL_ETOOLKIT_H
#define ARM_LOCAL_ETOOLKIT_H


// because of windows.h
#ifdef GetObject
#undef GetObject
#endif

#include <libCCTools++\CCVector.h>
#include <libCCTools++\CCString.h>

#include <string>
using namespace std;

class ARM_Date;
class ARM_result;


#define FFRETRIEVER		0
#define ETKRETRIEVER	1
#define WSETKRETRIEVER	2

// Flat File Fallback ETK
#define FFETKRETRIEVER	11


extern CCString RESCOMMSET; 


extern int DATARETRIEVER;



extern void checkETK();

extern long deconnection_etoolkit();

extern long shutdown_etoolkit();

extern long kill_etoolkit();

extern void set_etoolkit(const CCString& username,
				         const CCString& password,
				         const CCString& databasecontext,
				         const CCString& itConfigDomainsDir,
				         const CCString& itDomainName);

extern void connection_etoolkit();

extern void connection_etoolkit(const CCString& username,
						        const CCString& password,
						        const CCString& databasecontext,
						        const CCString& itConfigDomainsDir,
						        const CCString& itDomainName);

extern void connection_wsetoolkit(const CCString& pProd,
								  const CCString& pBase,
								  const CCString& pUserName);

extern long etoolkit_execute(CCString command,
					         CCString xmlRequest,
					         CCString & xmlResponse,
					         CCString & messageList);

extern long etoolkit_getcommsetname(const CCString& index,
							        const CCString& ccy,
							        const CCString& cvname,
							        const CCString& type,
 							        const CCString& typeData,
							        ARM_Date asof,
							        CCString & xmlResponse,
							        CCString & messageList,
							        const char* matuIndex = NULL,
									const char* CallorPut = "C");

extern long etoolkit_newgetcommsetname(const CCString& index,
									   const CCString& ccy,
									   const CCString& cvname,
									   const CCString& type,
									   const CCString& typeData,
									   ARM_Date asof,
									   CCString & xmlResponse,
									   CCString & messageList,
									   const char* matuIndex = NULL,
									   const char* CallorPut = "C");

extern long etoolkit_getrefrate(const CCString& source,
								const CCString& ccy,
								const CCString& index,
								const CCString& tenor,
								CCString& xmlResponse,
								CCString& messageList);

extern long etoolkit_connecte();

extern long etoolkit_getasof();

extern long etoolkit_setCurveId(const CCString& cvname);

extern long etoolkit_getlastdatewarm(CCString& xmlResponse);

extern long etoolkit_getVolCurveByTenor(const CCString& sortRequete,
	   							        const CCString& cvname,
								        ARM_Date asof,
								        const CCString& typeForGVC,
										const CCString& impOrHist,
								        CCString& xmlResponse,
								        CCString& messageList);

extern long etoolkit_getSmileCurveByStrike(const CCString& strike,
									       const CCString& request,
									       const CCString& cvname,
									       ARM_Date asof,
									       const CCString& type,
									       CCString& xmlResponse,
									       CCString& messageList);

extern long etoolkit_getStringVolCurveByTenor(const CCString& sortRequete,
									          const CCString& cvname,
									          ARM_Date asof,
									          const CCString& typeForGVC,
									          CCString& xmlResponse,
									          CCString& messageList);

extern long etoolkit_getCorrelCurveByMatu(const CCString& sortRequete,
								          const CCString& ccy2,
								          const CCString& index2,
								          const CCString& cvname,
								          ARM_Date asof,
								          CCString& xmlResponse,
								          CCString& messageList);

extern long etoolkit_getfxvol(const CCString& ccy1,
					   const CCString& ccy2,
					   const CCString& cvName,
					   ARM_Date date,
					   const CCString& impOrHist,
					   CCString& xmlResponse,
					   CCString& messageList);

CCString etoolkit_getXMLZCFromSummit(const CCString& index,
											const CCString& currency,
											const CCString& cvName,
											ARM_Date aSdate);

void etoolkit_getXMLZCFromCalypso(const std::string& index,
								  const std::string& ccy,
								  const std::string& term,
								  const std::string& forceCurveName,
									const std::string& pricingEnv,
									const ARM_Date& date,
									const std::string& xmlFileName,
									std::string& ret); 
										

extern CCString etoolkit_getXMLMYAndZCFromSummit(const CCString& index,
												 const CCString& currency,
												 const CCString& cvName,
												 ARM_Date aSdate);

extern CCString etoolkit_getXMLObjectFromSummit(const CCString& idSummit,
												const CCString& typeId);

extern long etoolkit_getfxcorrel(const CCString& ccy1,
						  const CCString& index,
						  const CCString& term,
						  const CCString& ccy2,
						  const CCString& cvName,
						  ARM_Date date,
						  CCString& xmlResponse,
						  CCString& messageList);

extern long etoolkit_getcorrel(const CCString& ccy1,
						const CCString& index1,
						const CCString& ccy2,
						const CCString& index2,
						const CCString& cvName,
						ARM_Date date,
		 				CCString& xmlResponse,
						CCString& messageList);

extern long etoolkit_getMeanRev(const CCString& ccy,
						 const CCString& index,
						 const CCString& cvName,
						 const CCString& C_NumFactor,
						 ARM_Date date,
						 CCString& xmlResponse,
						 CCString& messageList);

extern long etoolkit_getCutOff(const CCString& ccy,
						const CCString& index,
						const CCString& cvName,
						const CCString& NumFactor,
						ARM_Date date,
						CCString& xmlResponse,
						CCString& messageList);

extern long etoolkit_getFilter(const CCString& filter,
							   CCString& xmlResponse,
							   CCString& messageList);

extern long etoolkit_getFilterByStructId(const CCString& structId,
										 const CCString& desk,
										 CCString& xmlResponse,
										 CCString& messageList);

extern long etoolkit_iterateTradeList(const CCString& handle,
							   const CCString& action,
							   CCString& xmlResponse,
							   CCString& messageList);

extern long etoolkit_releaseIterateTradeList(const CCString& handle,
									  CCString& xmlResponse,
									  CCString& messageList);

extern long etoolkit_getQmodParam(const CCString& ccy,
						   const CCString& index,
						   const CCString& cvName,
						   ARM_Date date,
						   CCString& xmlResponse,
						   CCString& messageList);

extern long etoolkit_getQFXParam(const CCString& ccy,
						   const CCString& index,
						   const CCString& cvName,
						   ARM_Date date,
						   CCString& xmlResponse,
						   CCString& messageList);

extern long etoolkit_getFxSmileByStrike(const CCString& strike,	// ex "1000.000"
								const CCString& request,	//"SMILE/" & Dev1 & "/" & Dev2 & "/FXOPT/C/"
								const CCString& cvname,		//"MO"
								ARM_Date asof,
								const CCString& type,		//"SMILE"
								CCString& xmlResponse,
								CCString& messageList);

extern long etoolkit_getFxSmile(const CCString& request,	//"SMILE/" & Dev1 & "/" & Dev2 & "/FXOPT/C/"
								const CCString& cvname,
								ARM_Date asof,
								CCString& xmlResponse,
								CCString& messageList);

extern bool IsETKVersion();
extern int GetFallBackDataRetrieverVersion();
extern int GetDataRetrieverVersion();

extern void switchToETK(int withFallBack = 0);

extern void switchToWSETK();

extern void switchToFLATFILE();

extern CCString etoolkit_getXMLCapCashFlowFromSummit(const CCString& idSummit,
													const CCString& typeId);

extern long etoolkit_getFixing(const CCString& source,
							   const CCString& index,
							   const CCString& term,
							   const CCString& ccy,
							   ARM_Date asof,
							   CCString& xmlResponse,
							   CCString& messageList);

extern long etoolkit_getFactor(const CCString& ccy,
							   const CCString& index,
							   const CCString& cvName,
							   const CCString& Type,
							   const CCString& NumFactor,
							   ARM_Date date,
							   CCString& xmlResponse,
							   CCString& messageList);

extern CCString etoolkit_getXMLSeasonMgrFromSummit(const CCString& index,
												   const CCString& ccy,
												   const CCString& cvname,
												   ARM_Date date);
							   
extern VECTOR<string> etoolkit_getTradeListByStructId(const CCString& structureId,
													  const CCString& desk);

extern VECTOR<string> etoolkit_getTradeListByBook(const CCString& book);

extern CCString etoolkit_getXMLResetMgrFromSummit(const CCString& index,
												  const CCString& source,
												  const CCString& ccy,
												  const CCString& term);

extern CCString etoolkit_getListTenors(const CCString& index,
									   const CCString& ccy,
									   const CCString& cvname,
									   ARM_Date date,
									   const CCString& cvtype);

extern long etoolkit_getModelParam (ARM_Date date,
									const CCString& model,
									const CCString& type,
									const CCString& factorName,
									const CCString& ccy,
									const CCString& index,
									const CCString& cvName,
									CCString& xmlResponse,
									CCString& messageList);

#endif
/*--------------------------------------------------------------------------------*/
/*---- End Of File ----*/
