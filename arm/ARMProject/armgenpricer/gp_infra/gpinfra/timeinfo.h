/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: timeinfo.h,v $
 * Revision 1.3  2003/10/24 07:38:43  jmprie
 * Initial revision
 */

/*! \file timeinfo.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_TIMEINFO_H
#define _INGPINFRA_TIMEINFO_H

#include "gpbase/port.h"
#include "typedef.h"
#include "additionaltimeinfo.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
/// \class ARM_TimeInfo
/// \brief
/// TimeInfo class is a helper class for time information
/// from the parser. The parser is responsible for
/// converting date to time!
/////////////////////////////////////////////////////
class ARM_TimeInfo
{
private:
	double									itsEventTime;
	ARM_VectorPtr							itsPayTimes;
	ARM_AdditionalTimeInfoPtrVector*		itsAdditionalTimeInfos;

public:
	inline ARM_TimeInfo() : itsAdditionalTimeInfos( new ARM_AdditionalTimeInfoPtrVector(0) ) { }
	~ARM_TimeInfo() { delete itsAdditionalTimeInfos;}

	inline ARM_TimeInfo( double eventTime, const ARM_VectorPtr& payTimes )
	:	itsEventTime( eventTime ), itsPayTimes( payTimes ), 
	itsAdditionalTimeInfos( new ARM_AdditionalTimeInfoPtrVector(0) ) {}

	inline ARM_TimeInfo( double eventTime, const ARM_VectorPtr& payTimes, const ARM_AdditionalTimeInfoPtrVector& AdditionalTimeInfoPtrVector )
	:	itsEventTime( eventTime ), itsPayTimes( payTimes ), itsAdditionalTimeInfos( new ARM_AdditionalTimeInfoPtrVector(AdditionalTimeInfoPtrVector) ){}

    /// Accessors
	inline void SetEventTime(double time)  { itsEventTime = time;}
    inline double GetEventTime() const {return itsEventTime;}
	inline const std::vector<double> GetPayTimes() const{ return vector<double>(0);};// {return (*itsPayTimes).GetValues();}
	//inline std::vector<double>& GetPayTimes() const {return itsPayTimes->GetValues();}

	inline void AddAdditionalTimeInfo( const ARM_AdditionalTimeInfoPtr& adt ) { itsAdditionalTimeInfos->push_back(adt);}
	inline ARM_AdditionalTimeInfoPtrVector& getAdditionalTimeInfos() { return *itsAdditionalTimeInfos; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

