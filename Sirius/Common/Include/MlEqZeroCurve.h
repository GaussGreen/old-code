//	mleqzerocurve.h : Rate curve class. (Essentially a GDA wrapper but with a
//					  crude caching mechanism.)
//
//	author :          David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef __MLEQZEROCURVE_H_
#define __MLEQZEROCURVE_H_

#include "smart.h"

typedef std::map<int, double>					MlEqZeroCurveDiscountFactorCache;
												
class MlEqZeroCurve : public RCObject			
{												
public:											
	MlEqZeroCurve(void);
	explicit									MlEqZeroCurve(GDA::Functor_const_ref hYieldCurve);
	virtual										~MlEqZeroCurve(void);
												
	static std::string							Default(const std::string& szCurrency, const std::string& szRequest);
	double										GetCashRate(const std::string& szTermOrDate) const;
	DataSourceEnum								GetDataSource(void) const;
	double										GetDiscountFactor(long nMaturity) const;
	double										GetDiscountFactor(long nFrom, long nTo) const;
	double										GetFxSpot(void) const;
	std::string									GetInterpolator(void) const;
	std::string									GetName(void) const;
	long										GetReferenceDate(void) const;
	GDA::HDElement								GetYieldCurve(void) const;	
												
	void										PutDataSource(DataSourceEnum ds);
	void										PutFxSpot(double fFxSpot);
	void										PutInterpolator(const std::string& szInterpolator);
	void										PutReferenceDate(long nDate);
	void										PutShift(double fShift);
	void										PutShift(const MlEqInterpolatorHandle& hInterpolator);
	void										PutYieldCurve(GDA::Functor_const_ref hYieldCurve);
	void										Reset(void);
	void										Stick(void);

protected:		
	DataSourceEnum								m_ds;
	mutable GDA::Functor_const_ref				m_hYieldCurve;
	mutable MlEqZeroCurveDiscountFactorCache	m_mapDiscountFactors;			// map of precomputed discount factors
	GDA::Functor_const_ref						m_hYieldCurve_original;			// this is set only by PutYieldCurve and is therefore the mechanism for Reset()
	mutable	long								m_nCachedReferenceDate;			// cached from GDA	
	#ifdef _DEBUG
		std::string								m_szName;						// solely to facilitate debugging
	#endif

private:
	void										ClearCache(void) const;
};											
											
#endif
