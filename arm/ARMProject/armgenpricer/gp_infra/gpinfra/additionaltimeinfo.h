/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file additionaltimeinfo.h
 *
 *  \brief 
 *	\author  A.Schauly
 *	\version 1.0
 *	\date June 2005
 */

#ifndef _INGPINFRA_ADDITIONALTIMEINFO_H
#define _INGPINFRA_ADDITIONALTIMEINFO_H

#include "gpbase/port.h"
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_AdditionalTimeInfo
{
	enum AdditionalTimeInfoType
	{
		ADTNONE,
		ADTCPI,
	};

	virtual AdditionalTimeInfoType getTimeInfoType() { return ADTNONE; }

};

struct CPITimeInfo : public ARM_AdditionalTimeInfo
{
private:
	double itsPrevCPITime;
	double itsNextCPITime;
	double itsPrevPublishTime;
	double itsNextPublishTime;

	double itsCPITime;
	double itsPublishTime;

	string itsInfcurveName;
	string itsDCFLag;
	long itsDailyInterp;
	string itsResetLag;

public:
	inline void setCPITime( double CPITime ) { itsCPITime = CPITime; }
	inline void setPublishTime( double publishTime ) { itsPublishTime = publishTime; }
	inline void setPrevCPITime( double prevCPITime ) { itsPrevCPITime = prevCPITime; }
	inline void setNextCPITime( double nextCPITime ) { itsNextCPITime = nextCPITime; }
	inline void setPrevPublishTime( double prevPublishTime ) { itsPrevPublishTime = prevPublishTime; }
	inline void setNextPublishTime( double nextPublishTime ) { itsNextPublishTime = nextPublishTime; }

	inline double getCPITime() { return itsCPITime; }	
	inline double getPublishTime() { return itsPublishTime; }
	inline double getPrevCPITime() { return itsPrevCPITime; }
	inline double getNextCPITime() { return itsNextCPITime; }
	inline double getPrevPublishTime() { return itsPrevPublishTime; }
	inline double getNextPublishTime() { return itsNextPublishTime; }

	inline string getInfcurveName() { return itsInfcurveName; }
	inline string getDCFLag() { return itsDCFLag; }
	inline long getDailyInterp() { return itsDailyInterp; }
	inline string getResetLag() { return itsResetLag; }

	inline void setResetLag( const string& ResetLag ) { itsResetLag = itsResetLag; }
	inline void setDCFLag( const string& DCFLag ) { itsDCFLag = DCFLag; }
	inline void setInfcurvename( const string& infcurveName ) { itsInfcurveName = infcurveName; }
	inline void setDailyInterp( long dailyInterp ) { itsDailyInterp = dailyInterp; }

	virtual AdditionalTimeInfoType getTimeInfoType() { return ADTCPI; }

	CPITimeInfo( const string& InfcurveName, const string& DCFLag, const string& resetLag, long dailyInterp, 
		double CPITime, double publishTime,	double prevCPITime, double nextCPITime, double prevPublishTime, 
		double nextPublishTime) : 
	itsInfcurveName(InfcurveName), itsDCFLag(DCFLag), itsResetLag(resetLag), itsDailyInterp(dailyInterp), itsCPITime(CPITime), itsPublishTime(publishTime),
	itsPrevCPITime(prevCPITime), itsNextCPITime(nextCPITime),itsPrevPublishTime(prevPublishTime),itsNextPublishTime(nextPublishTime) {}

	CPITimeInfo( ) : 
	itsInfcurveName(""), itsDCFLag(""), itsResetLag(""), itsDailyInterp(0), itsCPITime(0), itsPublishTime(0),
	itsPrevCPITime(0), itsNextCPITime(0),itsPrevPublishTime(0),itsNextPublishTime(0) {}

	~CPITimeInfo() {}

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

