//	MlEqSpotSchedule.cpp : Implementation of the Spot Schedule
//
//	Author :			   David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mleqspotschedule.h"


DataSourceEnum MlEqSpotSchedule::GetDataSource(void) const
{
	return m_ds;
}

estring MlEqSpotSchedule::GetName(void) const
{
	return m_szName;
}

// specialisation to deal with the USD.USD case
double MlEqSpotSchedule::GetValueAt(long date) const
{
	try {
		return MlEqSchedule<long, double>::GetValueAt(date);
	} catch (const std::string& szError) {
		if (estring::CompareNoCase("USD.USD", m_szName)) throw szError;
		return 1.0;
	}
}

void MlEqSpotSchedule::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

void MlEqSpotSchedule::PutName(const std::string& szName)
{
	m_szName.assign(szName);
}
