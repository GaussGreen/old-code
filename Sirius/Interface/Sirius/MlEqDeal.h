//	MlEqDeal.h :			 Deal handle base class.
//							 This is an interface class.
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQDEAL_H
#define _MLEQDEAL_H

#pragma once

class MlEqDeal;
typedef RCPtr<MlEqDeal>							MlEqDealHandle;

#include "smart.h"

class MlEqDeal : public RCObject
{
public:	
	MlEqDeal(void);
	CComPtr<IResult>							Evaluate(const BSTR& Calculate) const;
	std::string									GetBook(void) const;
	DataSourceEnum								GetDataSource(void) const;
	long										GetDate(void) const;
	std::string									GetName(void) const;
	double										GetNotional(void) const;
	CComPtr<IPositions>							GetPositions(void) const;
	void										PutBook(const std::string& szBook);
	void										PutDataSource(DataSourceEnum ds);
	void										PutDate(long nDate);
	void										PutName(const std::string& szName);
	void										PutNotional(double fNotional);
	void										PutPositions(CComPtr<IPositions> spPositions);

protected:	
	DataSourceEnum								m_ds;
	double										m_fNotional;
	long										m_nDate;
	std::string									m_szBook;
	std::string									m_szName;
	CComPtr<IPositions>							m_spPositions;
};

#endif
