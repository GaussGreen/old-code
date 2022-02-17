/*! \file ARM_local_gp_inflationCapFloor.h
 *
 *
 *	\author  F. Poitou
 *
 */


#ifndef _INFCAPFLOOR_INTERFACE_H
#define _INFCAPFLOOR_INTERFACE_H


#include "firstToBeIncluded.h"

#include "glob\dates.h"
#include "ccy\currency.h"
#include "GP_Base\gpbase\countedptr.h"				//ARM_CountedPtr
#include "GP_Inflation\gpinflation\infcapfloor_.h"
#include <libCCtools++\CCstring.h>
#include "ARM_result.h"


/*
class InputDateStructure
{
		virtual std::string type() const { return "DateType"; }
};

class StringDate : public InputDateStructure
{
	public:
		StringDate(const std::string& s) : sDate_(s) {}
		virtual ~StringDate() {}
		const std::string& sDate() const { return sDate_;}
	protected:
		std::string sDate_;
};

class DoubleDate : public InputDateStructure
{
	public:
		DoubleDate(double afDate) : fDate_(afDate) {}
		virtual ~DoubleDate() {}
		double dDate() const { return fDate_;}
	protected:
		double fDate_;

};

struct DateFactory
{
	static ARM::ARM_Date create(	const ARM::ARM_CountedPtr<InputDateStructure>& ids);
	static ARM::ARM_Date create(	const ARM::ARM_CountedPtr<InputDateStructure>& ids,
												ARM::ARM_Date refDate, char* calendar,
												int settlementDays = 0,
												int rollingConvention = K_MOD_FOLLOWING);
	static ARM::ARM_Date create(const ARM::ARM_CountedPtr<InputDateStructure>& ids, ARM::ARM_Date refDate);
	static ARM::ARM_Date create(const ARM::ARM_CountedPtr<InputDateStructure>& ids, ARM::ARM_Date refDate, int settlementDays);
	static ARM::ARM_Date today() ;
};

class GenericDate
{
	public :
		GenericDate(const LPXLOPER& vData);
 		virtual ~GenericDate(){}

		ARM_CountedPtr<InputDateStructure> GetInputDateStructure() const { return itsIds ;}

	protected :
		ARM_CountedPtr<InputDateStructure> itsIds ;
};
*/

extern long ARMLOCAL_InfCapFloor_Load(
	double startDate,
	double endDate,
	const CCString& indexName,
	int capOrFloor,
	double strike,
	int interpType,
	int resetFreq,
	int dayCount,
	const CCString& resetCalendar,
	int fwdRule,
	int intRule,
	int stubRule,
	int resetGap,
	int payFreq,
	int payGap,
	const CCString& payCalendar,
	double firstReset,
	ARM_result&	result,
	long objId 	);




#endif


