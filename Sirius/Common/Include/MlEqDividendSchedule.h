//	MlEqDividendSchedule.h : dividend schedule class
//
//	Authors :				 David Cuin / Romain Barc
//
//
//	The dividend schedule maintains two sets of date / value schedules; one 
//  holds discrete dividends, the other holds continuous dividends.
//
//  The discrete case is simple. It is of the form
//		t_1			v_1
//		t_2			v_2
//		t_3			v_3
//		...			...
//		t_n			v_n
//  where v_i is the value of the dividend paid at t_1. The units of v_i are
//  understood to be the same as any asset that becomes associated with the
//  dividend schedule object.
//
//	For the continuous case, we still maintain the above form. However, the
//  value v_i defines the dividend yield from t_i-1 up to and including t_i.
//  Specifically v_1 is the yield from time 0 to t_1. The yield is zero for
//  all times after t_n.
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQDIVIDENDSCHEDULE_H_
#define _MLEQDIVIDENDSCHEDULE_H_

#include "mleqschedule.h" 
#include "smart.h"

class MlEqDividendSchedule : public RCObject
{
public:	
	explicit MlEqDividendSchedule(void);
	explicit MlEqDividendSchedule(const std::string& szName, DataSourceEnum ds, long nDate, DayCountConventionEnum dcc, MlEqZeroCurveHandle hForeignCurve, const std::map<long, double> mapDiscrete, const std::map<long, double> mapContinuous);
		
	void								ApplyGrowth(double fGrowth, const std::string& szContinuousYieldAfter);
	void								DiscreteShift(double fGrowth, const std::string& szTenor);	
	DataSourceEnum						GetDataSource(void) const;
	long								GetDate(void) const;
	DayCountConventionEnum				GetDayCountConvention(void) const;
	estring								GetDayCountConventionStr(void) const;	
	std::string							GetName(void) const;	
	static double						GetForward(long nStart,long nMaturity, double forwardStartSpot, MlEqZeroCurveHandle hzc, MlEqDividendScheduleHandle hDivs, MlEqZeroCurveHandle hzcFXYield);
	double								GetPresentValue(MlEqZeroCurveHandle hzc, long nFrom, long nTo, long nDiscountFromDate = 0);
	double								GetYield(long nFrom, long nTo);
	MlEqZeroCurveHandle					GetForeignCurve(void) const;
	void								CalibrateYield(double fSpot, MlEqZeroCurveHandle hzc, const std::string& szTenor, const std::string& szAddAt);
	void								PutContinuous(const std::map<long, double>& mapContinuous);
	void								PutContinuousValueAt(long nDate, double fValue);
	void								PutDataSource(DataSourceEnum ds);
	void								PutDate(long nDate);
	void								PutDayCountConvention(DayCountConventionEnum dcc);
	void								PutDiscrete(const std::map<long, double>& mapDiscrete);
	void								PutDiscreteValueAt(long nDate, double fValue);
	void								PutName(const std::string& szName);	
	void								PutForeignCurve(MlEqZeroCurveHandle hzc);
	void								Reset(void);
	void								Shift(double fYield);
	void								Stick(void);
	double								GetRate(long nDate);

public:
	// These constant references make data reading less obscure but still respect encapsulation.
	const MlEqSchedule<long, double>&	Discrete;
	const MlEqSchedule<long, double>&	Continuous;
	
protected:
	static const long					s_nMaxDate;
	MlEqSchedule<long, double>			m_discrete;						// Discrete dividends.
	MlEqSchedule<long, double>			m_continuous;					// Continuous dividends.
	DataSourceEnum						m_ds;
	DayCountConventionEnum				m_dcc;
	MlEqZeroCurveHandle					m_hForeignCurve;				// This allows us to associated a 'foreign' zero curve with a dividend schedule to act like a yield term.
	std::string							m_szName;
	
	// Scenario members
	std::map<long, double>				m_discrete_original;
	std::map<long, double>				m_continuous_original;
};

#endif

