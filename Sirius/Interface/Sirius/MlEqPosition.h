//	MlEqPosition.h :		 Position handle base class.
//							 This is an interface class.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQPOSITION_H
#define _MLEQPOSITION_H

#pragma once

class MlEqPosition;
typedef RCPtr<MlEqPosition>						MlEqPositionHandle;

#include "smart.h"

class MlEqPosition : public RCObject
{
public:	
	MlEqPosition(void);	
	CComPtr<IResult>							Evaluate(const BSTR& Calculate);
	std::string									GetBook(void) const;
	std::string									GetComment(void) const;
	std::string									GetCounterparty(void) const;
	DataSourceEnum								GetDataSource(void) const;
	long										GetDate(void) const;
	std::string									GetDescription(void) const;
	std::string									GetName(void) const;
	double										GetNotional(void) const;
	CComPtr<IProducts>							GetProducts(void) const;
	CComPtr<IResult>							GetResult(void) const;
	std::string									GetStrategy(void) const;
	CComPtr<IAssets>							GetUnderlyings(CComPtr<IDispatch> spObject) const;
	void										PutBook(const std::string& szBook);
	void										PutComment(const std::string& szComment);
	void										PutCounterparty(const std::string& szCounterparty);
	void										PutDataSource(DataSourceEnum ds);
	void										PutDate(long nDate);
	void										PutDescription(const std::string& szDescription);
	void										PutName(const std::string& szName);
	void										PutNotional(double fNotional);
	void										PutProducts(CComPtr<IProducts> spProducts);
	void										PutResult(CComPtr<IResult> spResult);
	void										PutStrategy(const std::string& szStrategy);

protected:	
	CComPtr<IProducts>							m_spProducts;
	std::string									m_szName;
	std::string									m_szDescription;
	std::string									m_szCounterparty;
	std::string									m_szBook;
	std::string									m_szStrategy;
	double										m_fNotional;
	std::string									m_szComment;
	CComPtr<IResult>							m_spResult;
	DataSourceEnum								m_ds;
	long										m_nDate;
};

#endif
