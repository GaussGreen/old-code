//	MlEqPuff.cpp: implementation of the MlEqPuff class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqPuff.h"

MlEqPuff::MlEqPuff()
{
	m_fMultiplier = 2.3;
	m_fStrike = -1.9;
	m_nDate = 32784;
}

double MlEqPuff::GetStrike() const
{
	return m_fStrike;
}

void MlEqPuff::PutStrike(double f)
{
	m_fStrike = f;
}


long MlEqPuff::GetMaturity() const
{
	return m_nDate;
}

void MlEqPuff::PutMaturity(long n)
{
	m_nDate = n;
}

double MlEqPuff::GetMultiplier() const
{
	return m_fMultiplier;
}

void MlEqPuff::PutMultiplier(double f)
{
	m_fMultiplier = f;
}

MlEqAssetHandle MlEqPuff::GetUnderlying() const
{
	return m_hUnderlying;
}

void MlEqPuff::PutUnderlying(MlEqAssetHandle h)
{ 
	m_hUnderlying = h;
}




double MlEqPuff::GetPrice() const
{
	
	if (!m_hUnderlying) throw "No underlying defined";


	double nToday = m_hUnderlying->GetDateHandle()->GetDate();

	double fForward = m_hUnderlying->GetQuantoForward(nToday, m_nDate, false);



	return fForward * m_fMultiplier;

}