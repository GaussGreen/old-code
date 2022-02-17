//	MlEqCorrelationMatrix.cpp : Implementation of the Correlation Matrix
//
//	Author :				    David Cuin
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MlEqCorrelationMatrix.h"
#include "utility.h"
#include "MlEqDate.h"


/////////////////////////////////////////////////////////////////////////////
//	MlEqCorrelationMatrix::helper implementation
//
MlEqCorrelationMatrix::helper::helper(const helper& h)
{
	if (this != &h){
		szL.assign(h.szL);
		szH.assign(h.szH);
	}
}

MlEqCorrelationMatrix::helper::helper(const std::string& sz1, const std::string& sz2)
{
	Set(sz1, sz2);
}		

MlEqCorrelationMatrix::helper::operator std::string(void) const
{
	return szL + "/" + szH;			
}

bool MlEqCorrelationMatrix::helper::operator<(const helper& h) const
{
	if (szL < h.szL) return true;
	if (szL > h.szL) return false;
	if (szH < h.szH) return true;
	return false;
}

bool MlEqCorrelationMatrix::helper::operator>(const helper& h) const
{
	if (szL > h.szL) return true;
	if (szL < h.szL) return false;
	if (szH > h.szH) return true;
	return false;
}	

const MlEqCorrelationMatrix::helper& MlEqCorrelationMatrix::helper::operator=(const helper& h)
{
	if (this != &h){
		szL.assign(h.szL);
		szH.assign(h.szH);
	}
	return *this;	
}

const std::string& MlEqCorrelationMatrix::helper::H(void) const
{
	return szH;
}

const std::string& MlEqCorrelationMatrix::helper::L(void) const
{
	return szL;
}

void MlEqCorrelationMatrix::helper::Set(const std::string& sz1, const std::string& sz2)
{						
	if (sz1 < sz2){
		szL.assign(sz1);
		szH.assign(sz2);
	} else {
		szL.assign(sz2);
		szH.assign(sz1);
	}
	trim(&szH);
	trim(&szL);	
}

/////////////////////////////////////////////////////////////////////////////
//	MlEqCorrelationMatrix::helper implementation
//
MlEqCorrelationMatrix::MlEqCorrelationMatrix(void) : m_ds(NoDataSource), m_nDate(0)
{
}

//	Merges an input correlation matrix with this matrix.
void MlEqCorrelationMatrix::Add(MlEqCorrelationMatrixHandle h, bool bReplace)
{
	if (bReplace){
		m_map = h->m_map;
		return;
	}
	for (std::map<helper, double>::const_iterator it = h->m_map.begin(); it != h->m_map.end(); it++){
		m_map[it->first] = it->second;
	}
}

void MlEqCorrelationMatrix::Clear(void)
{
	m_map.clear();
}

//	Retrieves a value from the correlation matrix.
double MlEqCorrelationMatrix::GetCorrelation(const std::string& sz1, const std::string& sz2) const
{
	if (sz1 == sz2) return 1.0;	
	if (sz1 == "USD.USD" || sz2 == "USD.USD") return 0.0;
	helper											h(sz1, sz2);			
	std::map<helper, double>::const_iterator		it = m_map.find(h);
	if (it == m_map.end()){
		throw "No correlation has been defined between '" + sz1 + "' and '" + sz2 + "'.";
	}
	double f = it->second;
	return f;
}

DataSourceEnum MlEqCorrelationMatrix::GetDataSource(void) const
{
	return m_ds;
}

long MlEqCorrelationMatrix::GetDate(void) const
{
	return m_nDate;
}

//	Returns a (strictly; not to be changed by any programmer!) constant reference to the map
const std::map<MlEqCorrelationMatrix::helper, double>& MlEqCorrelationMatrix::GetMap(void) const
{
	return m_map;
}

std::string MlEqCorrelationMatrix::GetName(void) const
{
	return m_szName;
}

//	Retrieves a volatility-volatiliy value from the correlation matrix.
double MlEqCorrelationMatrix::GetVolVolCorrelation(const std::string& sz1, const std::string& sz2) const
{
	if (sz1 == sz2) return 1.0;
	if (sz1 == "USD.USD" || sz2 == "USD.USD") return 0.0;
	helper											h("~" + sz1, "~" + sz2);
	std::map<helper, double>::const_iterator		it = m_map.find(h);
	if (it == m_map.end()){
		throw "No correlation has been defined between '~" + sz1 + "' and '~" + sz2 + "'.";
	}
	double f = it->second;
	return f;	
}
//	Returns TRUE if there are correlation data for sz1 and sz2
bool MlEqCorrelationMatrix::IsCorrelationDefined(const std::string& sz1, const std::string& sz2) const
{
	helper h(sz1, sz2);
	std::map<helper, double>::const_iterator it = m_map.find(h);	
	return it != m_map.end();
}

void MlEqCorrelationMatrix::PutDataSource(DataSourceEnum ds)
{
	m_ds = ds;
}

void MlEqCorrelationMatrix::PutDate(long nDate)
{
	m_nDate = nDate;
}

void MlEqCorrelationMatrix::PutName(const std::string& szName)
{
	m_szName.assign(szName);
}

//	Sets a value in the correlation matrix.
void MlEqCorrelationMatrix::SetCorrelation(const std::string& sz1, const std::string& sz2, double fCorrelation)
{
	if (sz1 == sz2){
		// must set fCorrelation to 1.0
		fCorrelation = 1.0;
	}
	helper h(sz1, sz2);
	m_map[h] = fCorrelation;
}