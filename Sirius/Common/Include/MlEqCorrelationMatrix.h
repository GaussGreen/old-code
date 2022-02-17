//	MlEqCorrelationMatrix.h: correlation matrix class
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQCORRELATIONMATRIX_H_
#define _MLEQCORRELATIONMATRIX_H_

#include "smart.h"

class MlEqCorrelationMatrix : public RCObject
{
public:
	class helper	// contans a correlation element (two strings such that the first one is not lexographically larger than the second)
	{		
	public:
		helper(const helper& h);		
		helper(const std::string& sz1, const std::string& sz2);		
		operator std::string(void) const;		
		bool operator<(const helper& h) const;		
		bool operator>(const helper& h) const;
		const helper& operator=(const helper& h);		
		const std::string&				H(void) const;
		const std::string&				L(void) const;
		
	protected:
		std::string						szH;							// not lexographically smaller than szL
		std::string						szL;							// not lexographically larger than szH
		
		void Set(const std::string& sz1, const std::string& sz2);		
	};

public:
	MlEqCorrelationMatrix(void);
	void								Add(MlEqCorrelationMatrixHandle h, bool bReplace);
	void								Clear();
	double								GetCorrelation(const std::string& sz1, const std::string& sz2) const;
	DataSourceEnum						GetDataSource(void) const;
	long								GetDate(void) const;
	std::string							GetName(void) const;
	const std::map<helper, double>&		GetMap(void) const;
	double								GetVolVolCorrelation(const std::string& sz1, const std::string& sz2) const;
	bool								IsCorrelationDefined(const std::string& sz1, const std::string& sz2) const;
	void								PutDataSource(DataSourceEnum ds);
	void								PutDate(long nDate);
	void								PutName(const std::string& szName);
	void								SetCorrelation(const std::string& sz1, const std::string& sz2, double fCorrelation);	

protected:
	std::map<helper, double>			m_map;							// each entry in the map is a correlation value
	DataSourceEnum						m_ds;
	long								m_nDate;
	std::string							m_szName;
};

#endif