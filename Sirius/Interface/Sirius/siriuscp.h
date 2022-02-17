#ifndef _SIRIUSCP_H_
#define _SIRIUSCP_H_

//#import "c:\progra~1\sirius\Debug\Sirius.dll" raw_interfaces_only, raw_native_types, no_namespace, named_guids	//"Import typelib"
extern "C" const GUID __declspec(selectany) __DIID__ISiriusApplicationEvents =
    {0xd22b96b4,0x1e36,0x4bbb,{0xa7,0xb4,0xb4,0xc0,0x15,0x75,0xd4,0xe8}};
// #include "Sirius_i.c"

//const IID _DIID__ISiriusApplicationEvents = {0xD22B96B4,0x1E36,0x4bbb,{0xA7,0xB4,0xB4,0xC0,0x15,0x75,0xD4,0xE8}};


template <class T>
class CProxy_ISiriusApplicationEvents : public IConnectionPointImpl<T, &__DIID__ISiriusApplicationEvents , CComDynamicUnkArray>
{
	//Warning this class may be recreated by the wizard.
public:
	HRESULT Fire_OnError(BSTR bstrError)
	{
		CComVariant varResult;
		T* pT = static_cast<T*>(this);
		int nConnectionIndex;
		CComVariant* pvars = new CComVariant[1];
		int nConnections = m_vec.GetSize();
		
		for (nConnectionIndex = 0; nConnectionIndex < nConnections; nConnectionIndex++)
		{
			pT->Lock();
			CComPtr<IUnknown> sp = m_vec.GetAt(nConnectionIndex);
			pT->Unlock();
			IDispatch* pDispatch = reinterpret_cast<IDispatch*>(sp.p);
			if (pDispatch != NULL)
			{
				VariantClear(&varResult);
				pvars[0] = bstrError;
				DISPPARAMS disp = { pvars, NULL, 1, 0 };
				pDispatch->Invoke(0x1, IID_NULL, LOCALE_USER_DEFAULT, DISPATCH_METHOD, &disp, &varResult, NULL, NULL);
			}
		}
		delete[] pvars;
		return varResult.scode;
	
	}
};
#endif