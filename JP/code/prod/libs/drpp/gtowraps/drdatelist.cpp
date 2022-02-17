// drdatelist.cpp: implementation of the drdatelist class.
//
//////////////////////////////////////////////////////////////////////

#include "drdatelist.h"
#include "drexception.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
DRDateList::DRDateList() 
: m_dateListGood(false), m_dateList(NULL), m_dateVec()
{
}

DRDateList::DRDateList (const DRDate& date)
: m_dateListGood(false), m_dateList(NULL), m_dateVec()
{
	m_dateVec.push_back(date);
}

DRDateList::DRDateList (TDateList* datelist)
: m_dateListGood(false), m_dateList(NULL), m_dateVec()
{
	for (int i = 0; i < datelist->fNumItems; i++) {
		m_dateVec.push_back (datelist->fArray[i]);
	}
	
}

DRDateList::DRDateList (const TDate* p, int size)
: m_dateListGood(false), m_dateList(NULL), m_dateVec(size)
{
	for (int i = 0; i < size; i++) 
		m_dateVec[i] = p[i];
}


DRDateList::DRDateList (const DRDateList& dates)
{
	m_dateVec = dates.m_dateVec;
	m_dateListGood = false;
	m_dateList = NULL;
}

DRDateList::~DRDateList()
{
	if (m_dateList)
		GtoFreeDateList (m_dateList);
}

DRDateList& DRDateList::operator= (const DRDateList& dates)
{
	if (this == &dates) return *this;
	
	m_dateVec = dates.m_dateVec;
	m_dateListGood = false;
	m_dateList = NULL;
	return *this;
}

bool DRDateList::operator==(const DRDateList& rhs) const
{
	if (m_dateVec.size() != rhs.m_dateVec.size()) return false;
	
	for (int i = 0; i < m_dateVec.size(); i++) {
		if (m_dateVec[i] != rhs.m_dateVec[i]) return false;
	}
	
	return true;
}

DRDateList& DRDateList::AddDate (const DRDate& date)
{
	DRDateVector temp;
	temp.push_back(date);
	
	m_dateVec = Merge (m_dateVec, temp);
	m_dateListGood = false;
	return *this;
}

DRDateList& DRDateList::AddDates (const DRDateList& dates)
{
	m_dateVec = Merge (m_dateVec, dates.m_dateVec);
	m_dateListGood = false;
	return *this;
}

DRDateList& DRDateList::AddDates (DRDate startDate, int numDates, const DRDateIn& dateInt)
{
	DRDateVector temp;
	for (int i = 0; i < numDates; i++) {
		temp.push_back(startDate);
		startDate += dateInt;
	}
	
	m_dateVec = Merge (m_dateVec, temp);
	m_dateListGood = false;
	return *this;
}

DRDateList& DRDateList::AddDates (DRDate startDate, const DRDate& endDate, const DRDateIn& dateInt)
{
	DRDateVector temp;
	while (startDate <= endDate) {
		temp.push_back(startDate);
		startDate += dateInt;
	}
	
	m_dateVec = Merge (m_dateVec, temp);
	m_dateListGood = false;
	return *this;
}

DRDateList& DRDateList::Delete (const DRDate& date)
{
	
	DRDateVector res = m_dateVec;
	int count = 0;
	for (int i = 0; i < m_dateVec.size(); i++) {
		if (m_dateVec[i] != date) {
			res[count++] = m_dateVec[i];
		}
	}
	m_dateVec = res;
	m_dateVec.resize(count);
	
	m_dateListGood = false;
	return *this;
}

bool DRDateList::IsDate (DRDate& date)
{
	return PerformBinarySearch(date).first;
}

pair<bool, int> DRDateList::PerformBinarySearch (DRDate& date, DateIndexType searchType)
{
	if (m_dateVec.size() == 0) return pair<bool, int> (false, 0);
	
	int begin = 0;
	int end = m_dateVec.size() - 1;
	
	while (end - begin > 1) {
		int mid = (begin + end) / 2;
		if (m_dateVec[mid] < date) begin = mid;
		else if (m_dateVec[mid] > date) end = mid;
		else {
			begin = end = mid;
			break;
		}
	}
	
	if (searchType == EXACT) {
		if (m_dateVec[begin] == date) return pair<bool,int> (true, begin);
		else if (m_dateVec[end] == date) return pair <bool,int> (true, end);
		else return pair<bool,int>(false, 0);
	}
	else if (searchType == ROUND_EARLIER) {
		if (m_dateVec[end] == date) return pair<bool,int> (true, end);
		else return pair<bool, int> (true, begin);
	}
	else {
		if (m_dateVec[begin] == date) return pair<bool,int> (true, begin);
		else return pair<bool, int> (true, end);
	}
}		


int DRDateList::GetDateIndex(DRDate date, DateIndexType searchType)
{
	pair<bool, int> pairVal = PerformBinarySearch(date, searchType);
	
	if (!pairVal.first) 
		throw DRException("Failure in GetDateIndex ") << date << " in " << *this;
	
	return pairVal.second;
}

int DRDateList::GetDateIndexIgnoreDOM (DRDate date)
{
	pair<bool, int> pairVal = PerformBinarySearch(date, ROUND_EARLIER);
	
	if (!pairVal.first) 
		throw DRException("Failure in GetDateIndexIgnoreDOM ") << date << " in " << *this;
	
	int index = pairVal.second;
	
	if (date.GetMonthYear() == m_dateVec[index].GetMonthYear()) return index;
	else if (date.GetMonthYear() == m_dateVec[index+1].GetMonthYear()) return index + 1;
	
	throw DRException("Failure to find date with good month year");
	
	return 0;
}

void DRDateList::BuildTDateList()
{
	if (m_dateListGood) return;
	
	if (m_dateList) {
		GtoFreeDateList(m_dateList);
	}
	
	m_dateList = GtoNewEmptyDateList (m_dateVec.size());
	
	for (int i = 0; i < m_dateVec.size(); i++) {
		m_dateList->fArray[i] = (TDate) m_dateVec[i];
	}
	m_dateListGood = true;
}

DRDateVector DRDateList::Merge (const DRDateVector& a1, const DRDateVector& a2)
{
	DRDateVector ans (a1.size() + a2.size());
	
	int count1 = 0;
	int count2 = 0;
	int count  = 0;
	DRDate newDate;
	
	while (count1 < a1.size() && count2 < a2.size()) {
		if (a1[count1] < a2[count2]) {
			newDate = a1[count1++];
		}
		else if (a1[count1] == a2 [count2]) {
			newDate = a1[count1++];
			count2++;
		}
		else {
			newDate = a2[count2++];
		}
		ans[count++] = newDate;
	}
	
	while (count1 < a1.size()) ans[count++] = a1[count1++];
	
	while (count2 < a2.size()) ans[count++] = a2[count2++];
	
	ans.resize (count);
	return ans;
}

void DRDateList::print_short (ostream& s) const
{
	if (size() < 3) {
		for (int i = 0; i < size(); i++) s << operator[](i) << "\t";
		s << endl;
	}
	else {
		for (int i = 0; i < 2; i++) s << operator[](i) << "\t";
		s << "...\t" << operator[](size() -1) << endl;
	}
}


ostream& operator<<(ostream& stream, const DRDateList& a)
{
	stream << "<" << a.size() << ">: ";
	for (int i = 0; i < a.size(); i++) {
		stream << GtoFormatDate((TDate) a.m_dateVec[i]) << " " ;
	}
	stream << endl;
	
	return stream;
}


