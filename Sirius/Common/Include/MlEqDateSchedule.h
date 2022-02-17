//	MlEqDateSchedule.h :	Date schedule class. Holds a vector of doubles
//							against each date.
//
//	Author :				David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQDATESCHEDULE_H_
#define _MLEQDATESCHEDULE_H_

#include "mleqschedule.h"

#undef max

class MlEqDateSchedule : public MlEqSchedule<long, std::vector<double> >
{
public:
	void								Create(long nStartDate, long nEndDate, const std::string& szFrequency, const std::string& szCalendar, BusinessDayConventionEnum bdc);
	long								GetWidth(long nIncludeBefore = std::numeric_limits<long>::max(), long* pnElementsConsidered = NULL) const;
	void								GetColumnBeforeToday(long nToday, long nColumn, std::vector<double>& af) const;
	void								GetColumnBeforeToday(long nToday, long nColumn, CMatrix& m) const;
	void								GetColumnsBeforeToday(long nToday, CMatrix& m) const;
};

#endif