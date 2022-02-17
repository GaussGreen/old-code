#ifndef __drrate__h
#define __drrate__h

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drstring.h"
#include <map>

extern "C" {
#include "ldate.h"
//#include "mbsutils.h"
#include "fltratea.h"
}

// Wrapped description of a rate
// There's some weird Gto representation called TFloatRateArray, which is too
// complicated for me.  Thus, I keep my rates around as LIBOR1M, etc.

// In theory, I have created as extern global variables all the rates you need.
// You should NOT ever create your own rates.  Simply, use the global variables or
// add to the list as needed

class DRRate {
public:
	DRRate(string name, long cpnsPerYear, long matInMonths, long dayCountConv, 
		long curveIndex, long numSettleDays);

	~DRRate();

	operator TFloatRateArray*() const;	// cast into the gto interp

	bool operator==(const DRRate&) const;	// comparisons are done by looking at the
	bool operator!=(const DRRate&) const;	// address only -- remember I create all 
											// of the rates -- do not make your own copies

	long cpnsPerYear() const {return m_cpnsPerYear;}
	long matInMonths() const {return m_matInMonths;}
	long dayCountConv() const {return m_dayCountConv;}
	long curveIndex() const {return m_curveIndex;}
	long numSettleDays() const {return m_numSettleDays;}
	const upper_string& name() const {return m_name;}

	friend ostream& operator<< (ostream&, const DRRate&);
//	friend istream& operator>> (istream&, DRRate*);
private:
	DRRate& operator=(const DRRate&);
	DRRate(const DRRate&);

	TFloatRateArray* m_floatRateArray;
	void SetFloatArray (long cpnsPerYear, long matInMonths, long dayCountConv, 
		long curveIndex, long numSettleDays);

	long m_cpnsPerYear;
	long m_matInMonths;
	long m_dayCountConv;
	long m_curveIndex;
	long m_numSettleDays;
	upper_string m_name;
	static map <upper_string, DRRate*, less <upper_string>, MYALLOC(DRRate*) > m_allRates;
};

// Note libor is always given as the gto_curve_discount and cmt is the
// gto_curve_index1 -- I can compute zeros off either, but in the virtual
// tree, I will diffuse libor
extern DRRate LIBOR1M;
extern DRRate LIBOR3M;
extern DRRate LIBOR1Y;
extern DRRate LIBOR2Y;
extern DRRate LIBOR10Y;

extern DRRate CMT1M;
extern DRRate CMT3M;
extern DRRate CMT1Y;
extern DRRate CMT2Y;
extern DRRate CMT10Y;

extern DRRate REFIRATE;

const DRRate* TranslateRate (DRString);
DRString TranslateRate (const DRRate&);

inline bool DRRate::operator==(const DRRate &rhs) const {return (this == &rhs);}

inline bool DRRate::operator!=(const DRRate &rhs) const {return !(*this == rhs);}

inline DRRate::operator TFloatRateArray*() const {return m_floatRateArray;}

inline DRRate::~DRRate () {GtoFloatRateArrayFree (m_floatRateArray);}

#endif





