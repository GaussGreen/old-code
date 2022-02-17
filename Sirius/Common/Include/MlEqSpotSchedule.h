//	MlEqSpotSchedule.h :	spot schedule class
//
//	Author :				David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQSPOTSCHEDULE_H_
#define _MLEQSPOTSCHEDULE_H_

#include "mleqschedule.h"
#include "smart.h"

class MlEqSpotSchedule : public MlEqSchedule<long, double>						 
{
public:
	MlEqSpotSchedule(bool (*func_extrapolate)(void) = NULL) : MlEqSchedule<long, double>(func_extrapolate) , m_ds(NoDataSource) {}
	DataSourceEnum						GetDataSource(void) const;	
	estring								GetName(void) const;
	double								GetValueAt(long date) const;
	void								PutDataSource(DataSourceEnum ds);	
	void								PutName(const std::string& szName);

protected:
	estring								m_szName;
	DataSourceEnum						m_ds;
};

#endif