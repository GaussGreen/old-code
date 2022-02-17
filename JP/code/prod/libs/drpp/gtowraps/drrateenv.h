// drrateenv.h: interface for the drrateenv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRRATEENV_H__5DD73672_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRRATEENV_H__5DD73672_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drstruct.h"
#include "drswapvols.h"
#include "drbasevols.h"
#include "drzerocurves.h"

class DRRateEnv  
{
public:
	DRRateEnv (MbsIO&);
	DRRateEnv (DRDate& envDate, DRZeroCurves&, DRBaseVols&, DRSwapVols&, 
		double meanReversion);

	DRDate& envDate ();
	DRSwapVols& swapVols();
	DRBaseVols& baseVols();
	DRZeroCurves& zeroCurves();
	DArray& calibParams();

	friend ostream& operator<<(ostream& s, DRRateEnv& a);
	bool operator==(const DRRateEnv&) const;

protected:
	DRDate m_envDate;

	DRSwapVols m_swapVols;
	DRBaseVols m_baseVols;
	DRZeroCurves m_zeroCurves;

	DArray m_calibParams;

};

inline DRDate& DRRateEnv::envDate () {return m_envDate;}
inline DRSwapVols& DRRateEnv::swapVols() {return m_swapVols;}
inline DRBaseVols& DRRateEnv::baseVols() {return m_baseVols;}
inline DRZeroCurves& DRRateEnv::zeroCurves() {return m_zeroCurves;}
inline DArray& DRRateEnv::calibParams() {return m_calibParams;}

typedef DRPtr<DRRateEnv> DRRateEnvPtr;

#endif // !defined(AFX_DRRATEENV_H__5DD73672_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
