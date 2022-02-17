#ifndef ARM_LOCAL_PARSEXML_NT_SECURITY_H
#define ARM_LOCAL_PARSEXML_NT_SECURITY_H

#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
#include <ICMKernel\inst\icm_mez.h>
#include <ICMKernel\inst\icm_nthtd.h>
#include <ICMKernel\inst\icm_option.h>
#include <ICMKernel\util\icm_pgcd.h>
#include <ICMKernel\inst\icm_cdo2.h>

void ARMLOCAL_XML_STDINFOCREDIT(const char* chaineXML, 
								std::string& bookName, 
								std::string& structureId, 
								std::string& custId, 
								std::string& dealId,
								std::string& summitProdType);

extern ICM_Nthtd* ARMLOCAL_XML_NTD(const char* chaineXML, 
								   CCString& bookName, 
								   CCString& custId, 
								   CCString& dealId,
								   ARM_ReferenceValue*fees);

extern ARM_ReferenceValue* ARMLOCAL_XML_EVENT(const char* chaineXML, CCString EventType);

void ARMLOCAL_XML_EVENT_FIXING(const char* chaineXML, CCString EventType,vector<double>& FixDates,vector<double>& Fixings);

/**
extern ICM_Mez* ARMLOCAL_XML_CDO(const char* chaineXML, 
								 CCString& bookId, 
								 CCString& custId, 
								 CCString& dealId,
								 ARM_ReferenceValue*fees);
								 **/ 

extern ICM_Option* ARMLOCAL_XML_CDS_OPTION(const char* chaineXML, 
									   CCString& bookName, 
									   CCString& custId, 
									   CCString& dealId,
									   ICM_Pricer* UnderlyingPricer);

extern ICM_Cds* ARMLOCAL_XML_CDS_SIMPLIFIED(const char* chaineXML,
										CCString& bookId, 
										CCString& custId, 
										CCString& dealId);

extern ICM_Cds* ARMLOCAL_XML_CDS(const char* chaineXML, 
								 CCString& bookId, 
								 CCString& custId, 
								 CCString& dealId,
								 ARM_ReferenceValue*fees); // might be null

void ConversionMnts(int nbnames, double* not_init,double sub_init, double end_init,
					double*& not_out, double& sub_out, double& tranche_out);


/** 
extern ICM_Credit_Index* ARMLOCAL_XML_CREDIT_INDEX(const char* chaineXML, 
												CCString& bookId, 
												CCString& custId, 
												CCString& dealId);
												**/ 

/**
ICM_Cdo2* ARMLOCAL_XML_CDO2_EMPTY(const char* chaineXML, 
								CCString& bookId, 
								CCString& custId, 
								CCString& dealId,
								ICM_Portfolio* Pf);
								**/ 

void ARMLOCAL_XML_COLLATERAL(const char* chaineXML, 
								  vector<string>& IssuersNames);


/**
extern ICM_Cdo2* ARMLOCAL_XML_CDO2(const char* chaineXML, 
								CCString& bookId, 
								CCString& custId, 
								CCString& dealId,
								ARM_ReferenceValue*fees);
								**/ 

extern void SearchTradeIdFromList(const vector<string>& TradeList,
						   vector<string>& TradeIdList,
						   vector<string>& TradeType);

#endif 
