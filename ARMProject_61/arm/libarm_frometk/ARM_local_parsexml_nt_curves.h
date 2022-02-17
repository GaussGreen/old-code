#ifndef ARM_LOCAL_PARSEXML_NT_CURVES_H
#define ARM_LOCAL_PARSEXML_NT_CURVES_H

#include <ARM\libarm_frometk\arm_local_parsexml_common.h>
#include <ARM\libarm_frometk\ARM_local_etoolkit_for_ICM.h>
#include <ARM\libarm_frometk\ARM_local_parsexml_for_ICM.h>
#include <ICMKernel\crv\icm_constant_piecewise.h>

extern ARM_ZeroCurve* ARMLOCAL_XML_ZCPY_with_summit_stripper(const char* chaineXML,
																 const char* chaineCCY);

extern ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_summit_stripper(const char* chaineXML,
																 const char* chaineDEFPROB,
																 ARM_ZeroCurve* ircurve_input = NULL);

ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_summit_stripper_Etk (const char* chaineXML,
																 ARM_ZeroCurve* ircurve_input,
																 CCString OverLabel);
ICM_DefaultCurve* ARMLOCAL_XML_DEFPROB_with_calypso_stripper(const std::string& chaineXML,
															 const ARM_ZeroCurve& ircurve,
															 const std::string& OverLabel);
#endif 
