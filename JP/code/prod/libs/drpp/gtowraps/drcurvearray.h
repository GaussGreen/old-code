// drcurvearray.h: interface for the drcurvearray class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRCURVEARRAY_H__5DD73674_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRCURVEARRAY_H__5DD73674_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000


#include "drcurve.h"

typedef vector < DRCurve, MYALLOC (DRCurve) > DRCurveVector;

class DRCurveArray  : public DRCurveVector
{
public:
	DRCurveArray(int n = 0) : DRCurveVector (n) {}
	operator TCurve**();
	friend ostream& operator<< (ostream&, const DRCurveArray&);

protected:
	DRValarray<TCurve*> m_tcurves;
};



#endif // !defined(AFX_DRCURVEARRAY_H__5DD73674_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
