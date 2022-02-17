#include <comdef.h>
#include <ARM\libarm_local\undef_va_vars.h>
#include <CCString.h>
#include <CCmessage.h>


//Copie de ComUtil.h
//CLEAR commenté car PLANTE quand Variant pas initialisé !!!???
void SetString(VARIANT *varparam,const char* pSrc) throw(_com_error)
{
	//
	// Free up previous VARIANT
	//
	//Clear(); 

	V_VT(varparam) = VT_BSTR;
	V_BSTR(varparam) = _com_util::ConvertStringToBSTR(pSrc);

	if (V_BSTR(varparam) == NULL && pSrc != NULL) {
		_com_issue_error(E_OUTOFMEMORY);
	}
}


HRESULT VARIANT2CCString (VARIANT var, CCString& str)
{
	//MSG_printf_message (MSG_INFO, "VARIANT2CCString (): entree dans VARIANT2CCString");

	_variant_t val (var);
	try
	{
		val.ChangeType(VT_BSTR);
	}
	catch(...)
	{
		throw("Pb in converting Variant To CCString");
		//MSG_printf_message (MSG_INFO, "VARIANT2CCString (): sortie de VARIANT2CCString: FALSE");
		return S_FALSE;
	}
		
	str.Set ((const char*)(_bstr_t)val);

	//MSG_printf_message (MSG_INFO, "VARIANT2CCString (): str = %s", (const char*)str);

	//MSG_printf_message (MSG_INFO, "VARIANT2CCString (): sortie de VARIANT2CCString: TRUE");

	return S_OK;
}

HRESULT VARIANT2Double (VARIANT var, double* outVar)
{
	_variant_t val (var);
	try
	{
		val.ChangeType(VT_R8);
	}
	catch(...)
	{
		return S_FALSE;
	}
		
	*outVar = (double)val;

	return S_OK;
}

HRESULT VARIANT2Long (VARIANT var, long* outVar)
{
	_variant_t val (var);
	try
	{
		val.ChangeType(VT_I4);
	}
	catch(...)
	{
		return S_FALSE;
	}
		
	*outVar = (long)val;

	return S_OK;
}

HRESULT Long2VARIANT (long inVar, VARIANT* var)
{
	VariantClear (var);
	var->vt = VT_I4;
	var->lVal = inVar;

	return S_OK;
}

HRESULT Double2VARIANT (double inVar, VARIANT* var)
{
	VariantClear (var);
	var->vt = VT_R8;
	var->dblVal = inVar;

	return S_OK;
}

HRESULT CCString2VARIANT (const CCString& inVar, VARIANT* var)
{
	_variant_t   tmp("rr");
	//VARIANT test;
	tmp.Attach(*var);
//	tmp.ChangeType(VT_BSTR);
//	SetString(&tmp,inVar.GetStr());
//	*var=tmp.Detach();
	SetString(&tmp,(const char*) inVar);
	// tmp.SetString((const char*)inVar); 
	*var=tmp.Detach();

	return S_OK;
}

HRESULT VARIANT2VECTORDOUBLE (VARIANT var, VECTOR<double>& outVar, long& size)
{
#ifdef _MD
	MSG_printf_message (MSG_INFO, "VARIANT2VECTORDOUBLE (): vt = %d", var.vt);
#endif
	SAFEARRAY* array = NULL;

	if((var.vt & VT_BYREF) != 0)
	{
		array = *(var.pparray);
	}
	else
	{
		array = var.parray;
	}

	if(array != NULL)
	{
#ifdef _MD
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORDOUBLE (): cDims = %d", array->cDims);
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORDOUBLE (): fFeatures = %d", array->fFeatures);
#endif
						
		if((array->cDims != 1) || (array->fFeatures & FADF_VARIANT == 0))
		{
			return S_FALSE;
		}
		int indexMax = array->rgsabound[0].lLbound + array->rgsabound[0].cElements;
		double tmp;

		size = array->rgsabound[0].cElements;
//		outVar = new double[size];

		for(long cur = array->rgsabound[0].lLbound; cur < indexMax; cur++)
		{
			VARIANT2Double (((VARIANT*)(array->pvData))[cur], &tmp);
			outVar.push_back(tmp);
		}
	}
    
	return S_OK;
}

HRESULT VARIANT2VECTORLONG (VARIANT var, VECTOR<long>& outVar, long& size)
{
#ifdef _MD
	MSG_printf_message (MSG_INFO, "VARIANT2VECTORLONG (): vt = %d", var.vt);
#endif
	SAFEARRAY* array = NULL;

	if((var.vt & VT_BYREF) != 0)
	{
		array = *(var.pparray);
	}
	else
	{
		array = var.parray;
	}

	if(array != NULL)
	{
#ifdef _MD
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORLONG (): cDims = %d", array->cDims);
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORLONG (): fFeatures = %d", array->fFeatures);
#endif
						
		if((array->cDims != 1) || (array->fFeatures & FADF_VARIANT == 0))
		{
			return S_FALSE;
		}
		int indexMax = array->rgsabound[0].lLbound + array->rgsabound[0].cElements;
		long tmp;

		size = array->rgsabound[0].cElements;
//		outVar = new double[size];

		for(long cur = array->rgsabound[0].lLbound; cur < indexMax; cur++)
		{
			VARIANT2Long (((VARIANT*)(array->pvData))[cur], &tmp);
			outVar.push_back(tmp);
		}
	}
    
	return S_OK;
}

HRESULT VARIANT2VECTORCCSTRING (VARIANT var, VECTOR<CCString>& outVar, long& size)
{
	#ifdef _MD
	MSG_printf_message (MSG_INFO, "VARIANT2VECTORCCSTRING (): vt = %d", var.vt);
	#endif

	SAFEARRAY* array = NULL;

	if((var.vt & VT_BYREF) != 0)
	{
		array = *(var.pparray);
	}
	else
	{
		array = var.parray;
	}

	if(array != NULL)
	{
		#ifdef _MD
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORCCSTRING (): cDims = %d", array->cDims);
		MSG_printf_message (MSG_INFO, "VARIANT2VECTORCCSTRING (): fFeatures = %d", array->fFeatures);
		#endif
						
		if((array->cDims != 1) || (array->fFeatures & FADF_VARIANT == 0))
		{
			return S_FALSE;
		}
		int indexMax = array->rgsabound[0].lLbound + array->rgsabound[0].cElements;
		CCString tmp;

		size = array->rgsabound[0].cElements;
//		outVar = new double[size];

		for(long cur = array->rgsabound[0].lLbound; cur < indexMax; cur++)
		{
			VARIANT2CCString (((VARIANT*)(array->pvData))[cur], tmp);
			outVar.push_back(tmp);
		}
	}
    
	return S_OK;
}

HRESULT MSG_printf_variant (VARIANT var)
{
	#ifdef _MD
	MSG_printf_message (MSG_INFO, "entree dans MSG_printf_variant ()");
	MSG_printf_message (MSG_INFO, "vt = %ld", var.vt);

	switch(var.vt)
	{
		case VT_BOOL:
			MSG_printf_message (MSG_INFO, "boolVal = %ld", var.boolVal);
			break;
		case VT_I2:
			MSG_printf_message (MSG_INFO, "iVal = %d", var.iVal);
			break;
		case VT_I4:
			MSG_printf_message (MSG_INFO, "lVal = %ld", var.lVal);
			break;
		case VT_R8:
			MSG_printf_message (MSG_INFO, "dblVal = %lf", var.dblVal);
			break;
		case VT_BYREF|VT_R8:
			MSG_printf_message (MSG_INFO, "pdblVal = %lf", *(var.pdblVal));
			break;
		case VT_BSTR:
			MSG_printf_message (MSG_INFO, "bstrVal = %s", (const char*)var.bstrVal);
			break;
		case VT_DISPATCH:
			MSG_printf_message (MSG_INFO, "pdispVal = %ld", var.pdispVal);
			break;
		case VT_BYREF|VT_BSTR:
			MSG_printf_message (MSG_INFO, "pbstrVal = %s", (const char*)*(var.pbstrVal));
			break;
		case VT_BYREF|VT_VARIANT:
			MSG_printf_message (MSG_INFO, "pvarVal = %ld", var.pvarVal);
			MSG_printf_variant (*(var.pvarVal));
			break;
		default:
			break;
	}

	MSG_printf_message (MSG_INFO, "sortie de MSG_printf_variant ()");
	#endif

	return S_OK;
}

HRESULT VECTORDOUBLE2VARIANT (VECTOR<double> inVar, VARIANT* var)
{
   VARIANT tmpVariant ;
   SAFEARRAYBOUND bound[1];    
   int size = inVar.size();

  // initialisation du safe array
  bound[0].lLbound = 0;
  bound[0].cElements = size;

  var->parray = SafeArrayCreate( VT_VARIANT , 1 , (SAFEARRAYBOUND*) &bound ) ;
  var->vt = VT_ARRAY | VT_VARIANT ;
 
  for (int j = 0;j<size; j++)
  {
      tmpVariant.dblVal = inVar[j]; 
      tmpVariant.vt = VT_R8;
      SafeArrayPutElement( var->parray , (long*)&j , &tmpVariant ) ;
   }

  return S_OK;

}

HRESULT VECTORCCSTRING2VARIANT (VECTOR<CCString> inVar, VARIANT* var)
{
   VARIANT tmpVariant ;
   SAFEARRAYBOUND bound[1];    
   int size = inVar.size();

  // initialisation du safe array
  bound[0].lLbound = 0;
  bound[0].cElements = size;

  var->parray = SafeArrayCreate( VT_VARIANT , 1 , (SAFEARRAYBOUND*) &bound ) ;
  var->vt = VT_ARRAY | VT_VARIANT ;
 
  for (int j = 0;j<size; j++)
  {
      tmpVariant.bstrVal = (_bstr_t)inVar[j]; 
      tmpVariant.vt =  VT_BSTR ;
      SafeArrayPutElement( var->parray , (long*)&j , &tmpVariant ) ;
   }

  return S_OK;
}


HRESULT VECTORDOUBLE2MATRIXVARIANT (VECTOR<double> inVar,int sizemat ,VARIANT* var)
{
   int size = inVar.size();
   int il = 0;

  if (size != sizemat*sizemat)
	  return S_FALSE;

   SAFEARRAYBOUND* bound = new SAFEARRAYBOUND[sizemat];    

  // initialisation du safe array
  for (il = 0; il<sizemat; il++)
  {
	(bound + il)->lLbound = (long) (sizemat - 1);
	(bound + il)->cElements = sizemat;
  }

  var->parray = SafeArrayCreate( VT_VARIANT , sizemat , (SAFEARRAYBOUND*) bound ) ;
  var->vt = VT_ARRAY | VT_VARIANT ;

  for (il = 0; il< sizemat; il++)
  {	
   VARIANT tmpVariant ;
	for (int j = 0;j<sizemat; j++)
	{
      tmpVariant.dblVal = inVar[j+il*sizemat]; 
      tmpVariant.vt = VT_R8;
	}

    SafeArrayPutElement( var->parray , (long*)&il , &tmpVariant ) ;
  }	

  if (bound)
	  delete[] bound;

  return S_OK;

}

HRESULT VECTORCCSTRING2MATRIXVARIANT (VECTOR<CCString> inVar,int sizemat,VARIANT* var)
{
   int size = inVar.size();
   int il = 0;

  if (size != sizemat*sizemat)
	  return S_FALSE;

   SAFEARRAYBOUND* bound = new SAFEARRAYBOUND[sizemat];    

  // initialisation du safe array
  for (il = 0; il<sizemat; il++)
  {
	(bound + il)->lLbound = (long) (sizemat - 1);
	(bound + il)->cElements = sizemat;
  }

  var->parray = SafeArrayCreate( VT_VARIANT , sizemat , (SAFEARRAYBOUND*) bound ) ;
  var->vt = VT_ARRAY | VT_VARIANT ;

  for (il = 0; il< sizemat; il++)
  {	
   VARIANT tmpVariant ;
	for (int j = 0;j<sizemat; j++)
	{
      tmpVariant.bstrVal = (_bstr_t)inVar[j]; 
      tmpVariant.vt =  VT_BSTR ;
	}

      SafeArrayPutElement( var->parray , (long*)&il , &tmpVariant ) ;
  }	

  if (bound)
	  delete[] bound;

  return S_OK;
}
