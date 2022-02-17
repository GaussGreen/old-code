// ARMFIModule.cpp : Implementation of ARMFIModule
#include "stdafx.h"
#include "ARMFILib.h"
#include "ARMFIModule.h"

/////////////////////////////////////////////////////////////////////////////
// ARMFIModule

STDMETHODIMP ARMCommonModule::InterfaceSupportsErrorInfo(REFIID riid)
{
	static const IID* arr[] = 
	{
		&IID_IARMFIModule
	};
	for (int i=0; i < sizeof(arr) / sizeof(arr[0]); i++)
	{
		if (InlineIsEqualGUID(*arr[i],riid))
			return S_OK;
	}
	return S_FALSE;
}
