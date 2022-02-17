//	MlEqDeal.cpp :			Implementation of MlEqDeal
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqDeal.h"
#include "MlEqDate.h"

MlEqDeal::MlEqDeal(void) : m_fNotional(1.0), m_nDate(0), m_ds(NoDataSource)
{
}

CComPtr<IResult> MlEqDeal::Evaluate(const BSTR& Calculate) const
{
	CComPtr<IResult>					spResult;

	if (m_spPositions->Evaluate(Calculate, &spResult)) propagate_error;
	if (spResult->Multiply(m_fNotional)) propagate_error;
	return spResult;
}

std::string MlEqDeal::GetBook(void) const
{
	return m_szBook;
}

DataSourceEnum MlEqDeal::GetDataSource(void) const
{
	return m_ds;
}

long MlEqDeal::GetDate(void) const
{
	return m_nDate;
}

std::string MlEqDeal::GetName(void) const
{
	return m_szName;
}

double MlEqDeal::GetNotional(void) const
{
	return m_fNotional;
}

CComPtr<IPositions> MlEqDeal::GetPositions(void) const
{
	return m_spPositions;
}

void MlEqDeal::PutBook(const std::string& szBook)
{
	m_szBook = szBook;
}

void MlEqDeal::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

void MlEqDeal::PutDate(long nDate)
{	
	m_nDate = nDate;
}

void MlEqDeal::PutName(const std::string& szName)
{
	m_szName = szName;
}

void MlEqDeal::PutNotional(double fNotional)
{
	m_fNotional = fNotional;
}

void MlEqDeal::PutPositions(CComPtr<IPositions> spPositions)
{
	m_spPositions = spPositions;
}