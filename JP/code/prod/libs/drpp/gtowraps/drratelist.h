#ifndef __drratelist__h
#define __drratelist__h

#include "drrate.h"
#include <vector>

// list of rates
// point to note:  the REFIRATE is always the last rate in the rate list

class DRRateList {
public:
	const DRRate& operator[](int) const;
	int size() const;
	void AddRate (const DRRate&);
	void AddRates (const DRRateList&);
	bool IsRate (const DRRate&);

	const DRRate& Last();
	const DRRate& First();

	int GetIndex (const DRRate&);

	bool operator==(const DRRateList&) const;

	friend ostream& operator<<(ostream&, const DRRateList&);

private:
	typedef vector<const DRRate*, MYALLOC(const DRRate*)> DRRateVector;
	DRRateVector m_rateVector;
};

inline const DRRate& DRRateList::operator[](int loc) const {return *(m_rateVector[loc]);}

inline const DRRate& DRRateList::Last() {return *(m_rateVector[m_rateVector.size()-1]);}

inline const DRRate& DRRateList::First() {return *(m_rateVector[0]);}

inline int DRRateList::size() const {return m_rateVector.size();}

#ifdef _TEST_CODE
bool TestDRRateList();
#endif

#endif

