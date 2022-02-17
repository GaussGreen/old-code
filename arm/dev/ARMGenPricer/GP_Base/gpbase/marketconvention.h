/*!
 *
 * Copyright (c) IXIS CI January 2005 Paris
 *
 *	\file marketconvention.h
 *
 *  \brief file for the market convention object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPBASE_MARKETCONVENTION_H
#define _INGPBASE_MARKETCONVENTION_H
 
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include <string>
CC_USING_NS(std,string)
#include <map>
CC_USING_NS(std,map)
CC_USING_NS(std,less)

CC_BEGIN_NAMESPACE( ARM )

struct ARM_MarketConventionData
{
	/// currency name
	char* itsName;

	/// adjust payment date to the previous (-1) 
	/// or next (+1) business day
    int itsFwdRule;

	/// Term of Libor rate
    int itsLiborTerm;

	/// day count basis for libor
	int itsLiborIndexDayCount;
	
	//// Nb of days after (>0) or prior (<0) 
	// payment will be made
    int itsSpotDays;
	
	/// day count basis for MM
	int itsMMDayCount;

	/// Payment frequency of the fixed leg
	/// of plain vanilla swap
    int itsFixedPayFreq;

	/// day count basis for fixed leg 
    int itsFixedDayCount;   

	/// Payment frequency for bond
    int itsCashPayFreq;

	/// day count basis for bond
    int itsCashDayCount;    
};

/// ARM_MarketConvention class
class ARM_MarketConvention : public ARM_RootObject
{
private:
	typedef map< string , ARM_MarketConventionData, less< string > >  MapStringToCurrencyData;
	static MapStringToCurrencyData* CreateCcyDataTable();
	static void ReleaseTheCcyDataTable();
	static MapStringToCurrencyData* TheCcyDataTable;
	ARM_MarketConventionData itsData;
	
public:
	ARM_MarketConvention( const string& name );

	/// accessors
	string GetName() const { return itsData.itsName; }
	inline int GetFwdRule() const { return itsData.itsFwdRule; }
    inline int GetLiborTerm() const { return itsData.itsLiborTerm; }
	inline int GetLiborIndexDayCount() const { return itsData.itsLiborIndexDayCount; }
    inline int GetSpotDays() const { return itsData.itsSpotDays; }
	inline int GetMMDayCount() const { return itsData.itsMMDayCount; }
    inline int GetFixedPayFreq() const { return itsData.itsFixedPayFreq; }
    inline int GetFixedDayCount() const { return itsData.itsFixedDayCount; }
    inline int GetCashPayFreq() const { return itsData.itsCashPayFreq; }
    inline int GetCashDayCount() const { return itsData.itsCashDayCount; }

	/// standard ARM_Object support
    virtual string toString(const string& indent="", const string& nextIndent="") const;
    virtual ARM_Object* Clone() const { return new ARM_MarketConvention(*this); }
};

/// equality test by the name of the market convention!
bool operator==(const ARM_MarketConvention& lhs, const ARM_MarketConvention& rhs )
{	return lhs.GetName() == rhs.GetName(); }

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

