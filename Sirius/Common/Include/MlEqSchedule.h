//	MlEqSchedule.h:			Schedule of dates and values class.
//
//	author:					David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _MLEQSCHEDULE_H
#define _MLEQSCHEDULE_H

#pragma once
#include "smart.h"
#include "MlEqDate.h"
template <class TDate, class TValue>				// typically, TDate is set to long or DATE
class MlEqSchedule : public std::map<TDate, TValue>, public RCObject
{	
public:	
	MlEqSchedule(bool (*func_extrapolate)(void) = NULL) : m_func_extrapolate(func_extrapolate), m_date(NULL) {}
	
	explicit MlEqSchedule( const std::map<TDate, TValue>& map, TDate date, bool (*func_extrapolate)(void) = NULL) : 
		std::map<TDate, TValue>(map), 
		m_date(date),
		m_func_extrapolate(func_extrapolate)
	{
	}
		
	TDate GetDate(void) const
	{
		return m_date;
	}
	
	// extract the dates from the schedule
	void GetDates(std::vector<TDate>& vector) const
	{
		vector.clear();
		for (std::map<TDate, TValue>::const_iterator itCount = begin(); itCount != end(); ++itCount){
			vector.push_back(itCount->first);
		}
		if (Extrapolate()){
			vector.push_back(m_date);
		}		
	}
	GVector<long> GetDates(void) const
	{
		GVector<long>		gv(size() + (Extrapolate() ? 1 : 0));
		long				n = 0;		
		for (std::map<TDate, TValue>::const_iterator itCount = begin(); itCount != end(); ++itCount){
			gv[n++] = itCount->first;
		}
		if (Extrapolate()){
			gv[n++] = m_date;
		}
		return gv;
	}
	void GetDates(GVector<TDate>& vector) const
	{
		vector.resize(size() + (Extrapolate() ? 1 : 0));
		long				n = 0;
		for (std::map<TDate, TValue>::const_iterator itCount = begin(); itCount != end(); ++itCount){
			vector[n++] = itCount->first;
		}
		if (Extrapolate()){		
			vector[n++] = m_date;
		}
	}
			
	// extract the values from the schedule
	void GetValues(std::vector<TValue>& vector) const
	{
		vector.clear();
		for (std::map<TDate, TValue>::const_iterator itCount = begin(); itCount != end(); ++itCount){
			vector.push_back(itCount->second);
		}
		if (Extrapolate()){
			vector.push_back(GetLastValue());
		}
	}

	// extract the dates and the values from the schedule
	void GetDatesAndValues(std::vector<TDate>& vectorDates, std::vector<TValue>& vectorValues) const
	{
		vectorDates.clear();
		vectorValues.clear();
		for (std::map<TDate, TValue>::const_iterator itCount = begin(); itCount != end(); ++itCount){
			vectorDates.push_back(itCount->first);
			vectorValues.push_back(itCount->second);
		}

		if (Extrapolate()){
			vectorDates.push_back(m_date);
			vectorValues.push_back(GetLastValue());
		}
	}

	TDate GetFirstDate(void) const
	{
		// No need to consider extrapolation
		if (!size()) throw "The schedule is empty.";
		return begin()->first;
	}

	TValue GetFirstValue(void) const
	{
		// No need to consider extrapolation
		if (!size()) throw "The schedule is empty.";
		return begin()->second;
	}

	TDate GetLastDate(void) const
	{
		if (!size()) throw "The schedule is empty.";
		if (Extrapolate()){
			return m_date;
		}
		std::map<TDate, TValue>::const_iterator it = end();
		return (--it)->first;		
	}

	TDate GetPenultimateDate(void) const
	{
		if (!size()) throw "The schedule is empty.";
		if (size() == 1) throw "The schedule contains only one element.";
		
		std::map<TDate, TValue>::const_iterator it = --end();		

		if (Extrapolate()){
			// Return the last date in the actual schedule.
			return it->first;			
		} else {
			return (--it)->first;
		}				
	}

	TValue GetLastValue(void) const
	{
		// There is no need to consider the extrapolation feature here since all it
		// would to is reappend the last value in the schedule; therefore the result
		// of this function would not change.
		if (!size()) throw "The schedule is empty.";				
		std::map<TDate, TValue>::const_iterator it = end();
		return (--it)->second;
	}

	virtual TValue GetValueAt(TDate date) const
	{
		if (!size()) throw "The schedule is empty.";
		std::map<TDate, TValue>::const_iterator it = find(date);
		if (it == end()){			
			if (Extrapolate() && (date == m_date)) return GetLastValue();			
			throw "No value is defined at " + MlEqDate(date).GetString() + ".";
		}
		return it->second;		
	}

	void GetValuesBeforeDate(long nDate, const std::vector<TDate>& vectorDates, std::vector<TValue>* pOut, const std::string& szName) const
	// szName - some name that we can associate with the schedule for error handling
	{
		if (!size()) throw "The schedule is empty.";
		pOut->clear();
		for (std::vector<TDate>::const_iterator it = vectorDates.begin(); it != vectorDates.end(); ++it){		
			if (*it >= nDate) break;
			try {
				pOut->push_back(GetValueAt(*it));
			} catch (const std::string&){
				if (szName.size()){
					throw "No value is defined at " + MlEqDate(*it).GetString() + " for schedule '" + szName + "'";
				} else {
					throw "No value is defined at " + MlEqDate(*it).GetString();
				}
			}
		}
	}

	bool HasValueAt(TDate date)
	{
		if (find(date) != end()) return true;
		if (Extrapolate() && date == m_date) return true;
		return false;
	}

	void Merge(const MlEqSchedule<TDate, TValue>& s)
	// values in the input schedule (s) take precedence
	{
		for (std::map<TDate, TValue>::const_iterator it = s.begin(); it != s.end(); ++it){
			(*this)[it->first] = it->second;
		}
	}
	
	void Multiply(const TValue& Amount)
	{		
		for (std::map<TDate, TValue>::iterator it = begin(); it != end(); ++it){
			it->second *= Amount;
		}
	}

	void PutDate(const TDate& date)
	{
		m_date = date;
	}

	void PutFirstValue(const TValue& v)
	{
		if (!size()) throw "The schedule is empty.";
		begin()->second = v;
	}
	
	void PutLastValue(const TValue& v)
	{
		if (!size()) throw "The schedule is empty.";
		std::map<TDate, TValue>::iterator it = end();
		(--it)->second = v;
	}

	void PutSchedule(const std::map<TDate, TValue>& map)
	{
		static_cast<std::map<TDate, TValue>& >(*this) = map;
	}

	void RemoveAfter(TDate date)
	{
		for (std::map<TDate, TValue>::iterator it = begin(); it != end(); ++it){
			if (it->first > date) erase(it);
		}
	}

	void RemoveAt(TDate date, bool bThrow)
	{
		if (!size()){
			if (bThrow) throw "The schedule is empty.":
			return;
		}

		std::map<TDate, TValue>::iterator it = find(date);
		if (it == end()){
			if (bThrow) throw "No value is defined at " + MlEqDate(date).GetString() + ".";
			return;
		}
		erase(it);
	}

	void RemoveOnOrAfter(TDate date)
	{
		for (std::map<TDate, TValue>::iterator it = begin(); it != end(); ++it){
			if (it->first >= date) erase(it);
		}
	}

protected:
	TDate								m_date;									// This is used by Extrapolate() and is to be understood as the current date.
	bool								(*m_func_extrapolate)(void);			// this function returns true if we are extrapolating today's value or not 
	bool								Extrapolate(void) const
	{
		if (!m_date) return false;		// can't do anything if the date has not been set up
		if (!size()) return false;
		if (!m_func_extrapolate || !m_func_extrapolate()) return false;
		TDate date = m_date;
		std::map<TDate, TValue>::const_iterator it = end();
		return (--it)->first < date;
	}
};

#endif