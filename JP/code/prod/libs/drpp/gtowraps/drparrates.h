// drparrates.h: interface for the drparrates class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRPARRATES_H__5DD73671_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)
#define AFX_DRPARRATES_H__5DD73671_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drstruct.h"
#include "drdatelist.h"
#include "tweaks.h"

class DRParRates  
{
public:
	DRParRates(DRDate& envDate, DRString cmtspdlFile = "cmtspd.usd", 
		DRString yldcrvlcFile = "yldcrvlc.usd");

	const DRDateList& dates() const;
	const DArray& parRates() const;
	const DArray& cmt10Yhumps() const;
	const DArray& swapSpreads () const;

	void TweakParRate(Tweak, double);
	friend ostream& operator<<(ostream& s, const DRParRates& a);

protected:
	DRDateList	m_rateDates;
	DArray		m_swapParRates;
	DArray		m_humpRates;
	DArray		m_liborSpreads;

	void LoadCmtspds (DRString filename);
	void LoadYldcrvlc (DRString filename, TDate startDate);
};

inline
const DRDateList& DRParRates::dates() const {return m_rateDates;}

inline
const DArray& DRParRates::parRates() const {return m_swapParRates;}

inline
const DArray& DRParRates::cmt10Yhumps() const {return m_humpRates;}

inline
const DArray& DRParRates::swapSpreads () const {return m_liborSpreads;}

typedef DRPtr<DRParRates> DRParRatesPtr;

#endif // !defined(AFX_DRPARRATES_H__5DD73671_A7A6_11D1_80A4_00C04FB91C08__INCLUDED_)

