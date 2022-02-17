//	siriusbase.h: Include this in all projects dependent on sirius
//
//	author:		  David Cuin
//
////////////////////////////////////////////////////////////////////////

#ifndef _SIRIUSBASE_H
#define _SIRIUSBASE_H

#include <limits>

const IID	IID_ExcelRange				= {0x00020846,0x0000,0x0000,{0xC0,0x00,0x00,0x00,0x00,0x00,0x00,0x46}};			// excel range object
const IID	IID_ExcelApplication		= {0x000208D5,0x0000,0x0000,{0xC0,0x00,0x00,0x00,0x00,0x00,0x00,0x46}};			// excel application object
const IID	_IID_ISiriusApplication		= {0x47230C12,0x55D8,0x4BC4,{0xBA,0x43,0xA8,0x39,0x44,0xCC,0x67,0x44}};
const CLSID _CLSID_SiriusApplication	= {0x82941A4A,0xE9B8,0x466F,{0x98,0x1D,0x5A,0xAF,0x92,0x3C,0xD4,0x08}};

#define VARIANT_MISSING										CComVariant(DISP_E_PARAMNOTFOUND, VT_ERROR)					// you can therefore take the address of this!
#define for if(false) {} else for		// Emulates the correct C++ grammar so variable instances declared within the for(...) statement are only within the for loop scope.

#include <MlEqHandles.h>
#include <ComDispatchDriverEx.h>
#include <ParameterMap.h>

#endif
