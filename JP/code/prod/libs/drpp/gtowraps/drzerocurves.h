// drzerocurves.h: interface for the drzerocurves class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRZEROCURVES_H__5DD73673_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRZEROCURVES_H__5DD73673_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drgtowrapenum.h"
#include "drcurve.h"
#include "drparrates.h"
#include "drcurvearray.h"
#include "drrate.h"

class DRZeroCurves  
{
public:
	DRZeroCurves() {}
	DRZeroCurves (TCurve* liborCurve, TCurve* cmtCurve);

	DRZeroCurves(FileType, DRDate& envDate, DRString cmtFile = "cmtspds.usd", 
		DRString yldFile = "yldcrvlc.usd");

	DRZeroCurves(DRDate& envDate, MbsIO&);

	DRCurve& discountCurve();
	DRCurveArray& indexCurves();

	void TweakParRate(Tweak tweak, double bpTweak);
	double GetForwardRate (DRRate&, DRDate&);
	bool operator==(const DRZeroCurves&) const;

	friend ostream& operator<< (ostream&, DRZeroCurves&);

protected:
	DRCurve	m_discountCurve;
	DRCurveArray m_indexCurves;
	void SetZerosFromParRates(DRDate& envDate);

	bool m_usingParRates;
	DRParRatesPtr m_parRatesPtr;
};

inline DRCurve& DRZeroCurves::discountCurve() {return m_discountCurve;}
inline DRCurveArray& DRZeroCurves::indexCurves() {return m_indexCurves;}


#endif // !defined(AFX_DRZEROCURVES_H__5DD73673_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
