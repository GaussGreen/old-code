//	MlEqPosition.cpp :		 Implementation of MlEqPosition
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqPosition.h"
#include "MlEqDictionary.h"
#include "MlEqDate.h"
#include "comobjectcollectionserialisablekey.h"

MlEqPosition::MlEqPosition(void) : m_fNotional(1.0), m_nDate(0), m_ds(NoDataSource)
{
}

CComPtr<IResult> MlEqPosition::Evaluate(const BSTR& Calculate)
{
	CComPtr<IResult>					spResult;
		
	if (m_spProducts->Evaluate(Calculate, &spResult)) propagate_error;
	if (spResult->Multiply(m_fNotional)) propagate_error;
	return m_spResult = spResult;
}

std::string MlEqPosition::GetBook(void) const
{
	return m_szBook;
}

std::string MlEqPosition::GetComment(void) const
{
	return m_szComment;
}

std::string MlEqPosition::GetCounterparty(void) const
{
	return m_szCounterparty;
}

DataSourceEnum MlEqPosition::GetDataSource(void) const
{
	return m_ds;
}

long MlEqPosition::GetDate(void) const
{
	return m_nDate;
}

std::string MlEqPosition::GetDescription(void) const
{
	return m_szDescription;
}

std::string MlEqPosition::GetName(void) const
{
	return m_szName;
}

double MlEqPosition::GetNotional(void) const
{
	return m_fNotional;
}

CComPtr<IProducts> MlEqPosition::GetProducts(void) const
{
	return m_spProducts;
}

CComPtr<IResult> MlEqPosition::GetResult(void) const
{
	return m_spResult;
}

std::string MlEqPosition::GetStrategy(void) const
{
	return m_szStrategy;
}

//	returns a collection of assets on which the product depends
CComPtr<IAssets> MlEqPosition::GetUnderlyings(CComPtr<IDispatch> spObject) const
{
	CComPtr<IAssets>					spAssets;

	if (m_spProducts->GetUnderlyings(&spAssets)) propagate_error;
	return spAssets;
}

void MlEqPosition::PutBook(const std::string& szBook)
{
	m_szBook = szBook;
}

void MlEqPosition::PutComment(const std::string& szComment)
{
	m_szComment = szComment;
}

void MlEqPosition::PutCounterparty(const std::string& szCounterparty)
{
	m_szCounterparty = szCounterparty;
}

void MlEqPosition::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

void MlEqPosition::PutDate(long nDate)
{	
	m_nDate = nDate;
}

void MlEqPosition::PutDescription(const std::string& szDescription)
{
	m_szDescription = szDescription;
}

void MlEqPosition::PutName(const std::string& szName)
{
	m_szName = szName;
}

void MlEqPosition::PutNotional(double fNotional)
{
	m_fNotional = fNotional;
}

void MlEqPosition::PutProducts(CComPtr<IProducts> spProducts)
{
	m_spProducts = spProducts;
}

void MlEqPosition::PutResult(CComPtr<IResult> spResult)
{
	m_spResult = spResult;
}

void MlEqPosition::PutStrategy(const std::string& szStrategy)
{
	m_szStrategy = szStrategy;
	estring::trim(&m_szStrategy);
}