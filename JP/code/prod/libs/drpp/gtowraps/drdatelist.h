// drdatelist.h: interface for the drdatelist class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRDATELIST_H__CC81DA9A_9984_11D1_80A2_00C04FB91C08__INCLUDED_)
#define AFX_DRDATELIST_H__CC81DA9A_9984_11D1_80A2_00C04FB91C08__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drdatein.h"
#include "drdate.h"
#include <vector>

extern "C" {
#include "convert.h"
#include "datelist.h"
}

enum DateIndexType {EXACT, ROUND_EARLIER, ROUND_LATER};

// This class wraps Gto TDateList (sorta of)
// My interpretation of datelist is that is a mainly a read-only structure.
// Note that element access does not return
// dates by reference, so you CANNOT change dates in the datelist via operator[]
//
// I do cast into a TDateList*, but it is for read purposes only.
// Any changes are not reflected in the original date list.


typedef KVector(DRDate) DRDateVector;

class DRDateList {
public:
	DRDateList ();				// constructors
	DRDateList (const DRDate&);
	DRDateList (TDateList*);
	DRDateList (const DRDateList&);
	DRDateList (const TDate*, int size);
	~DRDateList ();

	DRDateList& operator=(const DRDateList&);

	bool IsDate (DRDate&);	// checks whether this date is on the list

	operator TDateList*();	// casts into TDatelist* or TDate*
	operator TDate*();

	DRDate operator[](int) const;	// element access
	DRDate operator()(int) const;	// note that I did return a const DRDate&, cause
									// I'm too sloppy with my const's to make it work

	bool operator==(const DRDateList&) const;
	
	int size() const;

	DRDate First() const;
	DRDate Last() const;

	DRDateList& AddDate (const DRDate&);	// append single date
	DRDateList& AddDates (const DRDateList&);	// append a datelist

	// appends dates starting at start dates and apply dateIn, numDates times
	DRDateList& AddDates (DRDate startDate, int numDates, const DRDateIn&);

	// appends dates starting at startdate, ending before at at endDate, and
	// moving forward via dateIn
	DRDateList& AddDates (DRDate startDate, const DRDate& endDate, const DRDateIn&);
	
	DRDateList& Delete (const DRDate&);

	// Get the date index for a date.  The date index corresponds to the
	// access of the datelist via operator[] or operator()
	// The various search types are EXACT -- find only exact matches, error otherwise
	//								ROUND_EARLIER
	//								ROUND_LATER
	int GetDateIndex(DRDate, DateIndexType searchType = EXACT);

	// Performs date search ignoring day of month
	int GetDateIndexIgnoreDOM (DRDate);

	friend ostream& operator<<(ostream&, const DRDateList&);
	void print_short (ostream&) const;

private:
	DRDateVector m_dateVec;
	
	bool m_dateListGood;
	TDateList* m_dateList;

	void BuildTDateList();
	pair<bool, int> PerformBinarySearch (DRDate& date, DateIndexType searchType = EXACT);
	DRDateVector Merge (const DRDateVector& a1, const DRDateVector& a2);
};

inline DRDate DRDateList::operator() (int loc) const {return m_dateVec[loc];}

inline DRDate DRDateList::operator[] (int loc) const {return m_dateVec[loc];}

inline DRDate DRDateList::First() const {return m_dateVec.front();}

inline DRDate DRDateList::Last() const {return m_dateVec.back();}

inline int DRDateList::size() const {return m_dateVec.size();}

inline DRDateList::operator TDateList* ()
{
	BuildTDateList();
	return m_dateList;
}

inline DRDateList::operator TDate*()
{
	BuildTDateList();
	return m_dateList->fArray;
}

#endif // !defined(AFX_DRDATELIST_H__CC81DA9A_9984_11D1_80A2_00C04FB91C08__INCLUDED_)
