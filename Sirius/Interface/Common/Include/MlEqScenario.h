//	MlEqScenario.h :		 Scenario handle class
//							 THIS SHOULD NOT BE IN THE ANALYTICS LIBRARY
//                           SINCE IT CONTAINS VARIANTS ETC.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQSCENARIO_H
#define _MLEQSCENARIO_H

#pragma once

#include "smart.h"

class MlEqScenario : public RCObject
{
public:
	void								GetValue(CComVariant* pVal) const;
	void								PutValue(const CComVariant& newVal);

protected:
	CParameterMap						m_pm;
};

#endif