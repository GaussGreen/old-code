#ifndef ARM_LOCAL_PARSEXML_FOR_ICM_H
#define ARM_LOCAL_PARSEXML_FOR_ICM_H

#include <ARM\libarm_local\ARM_local_glob.h>
#include <ICMKernel\crv\icm_defaultcurve.h>
#include <ICMKernel\mod\modelmulticurves.h>

#include "ARM_local_etoolkit.h"
//#include "arm_local_parsexml_common.h"

class ICM_Cdo2;
double ConvertSecToRecovery(char* secured);

ARM_Vector* ARMLOCAL_RecoveryCurve (const char* chaineXML_Recovery);

ICM_DefaultCurve* ARMLOCAL_ParseDefProbCurveCalypso (const std::string& chaineXML_defcurve,
											const ARM_Date &CurveDate, 
											ARM_ZeroCurve *zcpy,
											const std::string& label,
											const std::string& Seniority,
											const std::string& Currency,
											const std::string& ForceCurveName,
											map<string,ARM_ZeroCurve*>* MercureMap=NULL);
extern ICM_DefaultCurve* ARMLOCAL_ParseDefProbCurve (const char* chaineXML_defcurve,
											ARM_Date CurveDate, 
											ARM_ZeroCurve *zcpy,
											const char* label,
											map<string,ARM_ZeroCurve*>* MercureMap = NULL);

extern long  ARMLOCAL_ParseDefProbCurveMktData (const char* chaineXML,
										VECTOR<CCString>& matu,
										VECTOR<double>& spread,
										VECTOR<double>& recovery,
										CCString& currency,
										CCString& indexname,
										long& AIMMADJ);

void ARMLOCAL_ParseDefProbCurveMktDataCalypso (const std::string & chaineXML,
											   std::vector<std::string >& matu,
											   std::vector<double>& spread,
											   double & recovery,
											   std::string & currency,
											   std::string & indexname,
											   qCDS_ADJ & AIMMADJ);

extern double ICMLOCAL_ParseCorrCredit(const char* chaineXML,const CCString& bookName,const CCString& structureId);


extern ARM_Model* ICMLOCAL_ParseModelCredit(ARM_ZeroCurve* DiscountCurve, 
									 const char* chaineXML, 
									 const CCString& SummitId, 
									 const CCString& Type,
									 const CCString& CurveId,
									 const CCString& CorrCurveId);

extern long  ARMLOCAL_GetSecidForRecovery (const char* chaineXML,CCString& SecId);

ARM_Model* ICMLOCAL_ParseModelForSpreadOptions(ARM_ZeroCurve* DiscountCurve, 
												const char* chaineXML);


void 
ICMLOCAL_ParseBasketCorrelMkDataFromCalypso(const std::string& xml,
										  std::vector<std::string>& matus,
										  std::vector<double>&attach,
										  ICM_QMatrix<double>& correls) ;
#endif 
