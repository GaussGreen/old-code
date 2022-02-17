#include "drratelist.h"
#include "drexception.h"

int DRRateList::GetIndex (const DRRate& rate)
{	
	for (int i = 0; i < size(); i++) {
		if (&rate == m_rateVector[i]) return i;
	}
	throw DRException ("Failed to find rate in ratelist");
	return 0;
}

bool DRRateList::IsRate(const DRRate& rate)
{
	bool found = false;

	DRRateVector::const_iterator iter = m_rateVector.begin();
	
	while (iter != m_rateVector.end()) {
		if (*iter == &rate) {
			found = true;
			break;
		}
		iter++;
	}

	return found;
}

void DRRateList::AddRate(const DRRate& rate)
{
	bool found = false;

	DRRateVector::const_iterator iter = m_rateVector.begin();
	
	while (iter != m_rateVector.end()) {
		if (*iter == &rate) {
			found = true;
			break;
		}
		iter++;
	}

	// make sure the REFIRATE is last
	if (!found) {
		if (m_rateVector.size() > 0 && m_rateVector.back() == &REFIRATE) {							

			m_rateVector.pop_back();
			m_rateVector.push_back(&rate);
			m_rateVector.push_back(&REFIRATE);
		}
		else {
			m_rateVector.push_back(&rate);
		}
	}
}

bool DRRateList::operator==(const DRRateList& a) const
{
	return (m_rateVector == a.m_rateVector);
}


void DRRateList::AddRates (const DRRateList& rates)
{
	DRRateVector::const_iterator iter;
	for (iter = rates.m_rateVector.begin(); iter != rates.m_rateVector.end(); iter++) {
		AddRate(**iter);
	}
	// Wow, why is this **?  Like who can understand that?  This can't happen in C++!?
	// Well, AddRate takes a rate.  iter is an iterator.	
	// The first * gives an elemenet in the vector.
	// The vector contains DRRate*, so the second * gives me the element
	// Ugly? yes
}

ostream& operator<<(ostream& s, const DRRateList& a)
{
	DRRateList::DRRateVector::const_iterator theIter;

	for (theIter = a.m_rateVector.begin(); theIter !=a.m_rateVector.end(); theIter++) {
		s << TranslateRate(**theIter) << " ";
	}
	s << endl;
	return s;
}

