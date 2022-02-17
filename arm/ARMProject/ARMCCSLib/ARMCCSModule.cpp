// ARMCCSModule.cpp : Implementation of ARMCCSModule
#include "stdafx.h"
#include "ARMCCSLib.h"
#include "ARMCCSModule.h"

/////////////////////////////////////////////////////////////////////////////
// ARMCCSModule

STDMETHODIMP ARMCommonModule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IARMCCSModule
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}
