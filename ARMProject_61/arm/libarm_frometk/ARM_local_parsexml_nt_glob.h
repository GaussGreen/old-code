#ifndef ARM_LOCAL_PARSEXML_NT_GLOB_H
#define ARM_LOCAL_PARSEXML_NT_GLOB_H

#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
#include <ICMKernel\glob\icm_corrmatrix.h>
#include <ICMKernel\glob\icm_smile_correlation.h>
#include <ICMKernel\util\icm_matrix.h>

extern ICM_CorrMatrix* ARMLOCAL_XML_CORRMATRIX(const char* chaineXML);
extern ICM_Smile_Correlation* ARMLOCAL_XML_CORRELATION_SMILE(const char* chaineXML);
extern ICM_Matrix<ARM_Vector>* ARMLOCAL_XML_CASHFLOWS(const char* chaineXML, 
										  const CCString& Path, 
										  CCString& bookName, 
										  CCString& structureId, 
										  CCString& custId, 
										  CCString& dealId);

/**
extern ICM_Matrix<ARM_Vector>* ARMLOCAL_XML_PRICING_CMCDS(const char* chaineXML, 
											  const CCString& Path, 
											  CCString& bookName, 
											  CCString& structureId, 
											  CCString& custId, 
											  CCString& dealId);
											  **/ 

extern ARM_Vector* ARMLOCAL_XML_VOLCURVE_DEFPROB_BY_TENOR(const char* chaineXML,
											  const string& DefProbName,
											  const string& TENOR,
											  string RatesOrDate /*R ou D*/,
											  ARM_Date AsOf);

// FIXMEFRED: mig.vc8 (30/05/2007 17:39:51):missing return type
extern void ARMLOCAL_XML_VOLCURVE_DEFPROB_MATURITIES(const char* chaineXML,
												const string& DefProbName,
											    vector<string>& Maturities,
												ARM_Vector*& YearFrac,
												ARM_Date& AsOf);

extern ARM_VolCurve* ARMLOCAL_XML_VOLCURVE_DEFPROB(const char* chaineXML,const string& DefProbName);

#endif 
