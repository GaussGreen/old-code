// drcurve.h: interface for the drcurve class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRCURVE_H__5DD73654_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRCURVE_H__5DD73654_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drdatelist.h"
#include "drstruct.h"

extern "C" {
#include "tcurve.h"
}


class DRCurve
{
public:
	DRCurve();
	DRCurve(TCurve*);

	DRCurve(DRDateList& dates, DArray& rates, DRDate& baseDate, double basis, 
		long dayCountConv);

	DRCurve(DRDate& envDate, DArray& parRates, DRDateList& dates);
	// creates a DRCurve from a mixed combination of par/swap rates 
	// assumes the long convention of 25 points

	DRCurve(DRString filename);	// assumes kapital file!!
	
	virtual ~DRCurve();
	
	DRCurve(const DRCurve&);
	DRCurve& operator=(const DRCurve&);
	
	DRDateList GetDates () const;
	DArray GetRates () const;
	
	operator TCurve* ();
	
	int size() const;
	DRDate baseDate() const;
	double basis() const;
	long dayCountConv() const;
	bool operator==(const DRCurve&) const;
	
	friend ostream& operator<<(ostream&, const DRCurve&);
	
protected:
	virtual void Print(ostream&) const ;
	TCurve* m_tcurve;
};

inline int DRCurve::size() const {return m_tcurve->fNumItems;}
inline DRDate DRCurve::baseDate() const {return DRDate(m_tcurve->fBaseDate);}
inline double DRCurve::basis() const {return m_tcurve->fBasis;}
inline long DRCurve::dayCountConv() const {return m_tcurve->fDayCountConv;}
inline DRCurve::operator TCurve* () {return m_tcurve;}



#endif // !defined(AFX_DRCURVE_H__5DD73654_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)

