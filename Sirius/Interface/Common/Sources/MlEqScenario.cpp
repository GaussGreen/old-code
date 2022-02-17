//	MlEqScenario.cpp :		 Implementation of MlEqScenario
//							 THIS SHOULD NOT BE IN THE ANALYTICS LIBRARY
//                           SINCE IT CONTAINS VARIANTS ETC.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqScenario.h"

void MlEqScenario::GetValue(CComVariant* pVal) const
{
	m_pm.GetValue(pVal);
}
void MlEqScenario::PutValue(const CComVariant& newVal)
{
	m_pm.SetValue(newVal);
}
