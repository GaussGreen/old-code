#ifndef CCATL_H
#define CCATL_H


#ifdef STL_WIN32
#include <vector>
#define VECTOR	std::vector
#else
#include <vector.h>
#define VECTOR	vector
#endif	// STL_WIN32

#include <CCString.h>



HRESULT VARIANT2CCString (VARIANT var, CCString& str);
HRESULT VARIANT2Double (VARIANT var, double* outVar);
HRESULT VARIANT2Long (VARIANT var, long* outVar);
HRESULT Long2VARIANT (long inVar, VARIANT* var);
HRESULT Double2VARIANT (double inVar, VARIANT* var);
HRESULT CCString2VARIANT (const CCString& inVar, VARIANT* var);
HRESULT VARIANT2VECTORDOUBLE (VARIANT var, VECTOR<double>& outVar, long& size);
HRESULT VARIANT2VECTORLONG (VARIANT var, VECTOR<long>& outVar, long& size);
HRESULT VARIANT2VECTORCCSTRING (VARIANT var, VECTOR<CCString>& outVar, long& size);
HRESULT MSG_printf_variant (VARIANT var);
HRESULT VECTORDOUBLE2VARIANT (VECTOR<double> inVar, VARIANT* var);
HRESULT VECTORCCSTRING2VARIANT (VECTOR<CCString> inVar, VARIANT* var);
HRESULT VECTORDOUBLE2MATRIXVARIANT (VECTOR<double> inVar,int sizemat ,VARIANT* var);
HRESULT VECTORCCSTRING2MATRIXVARIANT (VECTOR<CCString> inVar,int sizemat,VARIANT* var);



#define V_readNumVariant(V_var,C_var,error_msg,result)			if((error = VARIANT2Double (V_var, &C_var)) != S_OK)\
																{\
	  																	CCString local_msg (error_msg);\
																		result.setMsg (local_msg);\
																		ARM_ERR();\
																		return &V_result;\
		  														}
															
#define V_readStrVariant(V_Str,C_Str,error_msg,result)			if((error = VARIANT2CCString (V_Str, C_Str)) != S_OK)\
																{\
																	CCString local_msg (error_msg);\
																	result.setMsg (local_msg);\
																	ARM_ERR();\
																	return &V_result;\
		  														}



#endif	// CCATL_H
