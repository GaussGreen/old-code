#ifndef _INGPCALCULATORS_CSOQUANTOCALCULATOR_H
#define _INGPCALCULATORS_CSOQUANTOCALCULATOR_H

#include "callablequantocalculator.h"

#include "gpcalib/typedef.h"


CC_BEGIN_NAMESPACE( ARM )


struct ScheduleInformation {
	int frequency; 
	int isAdjusted; 
	int dayCounter; 
	int resetType; 
	int resetGap; 
	string resetCalendar;
	string paymentCalendar;
}; 


//! CallableQuantoSpreadOptionCreator allows to create a GenSecurity corresponding to a callable QuantoSpreadOption
/*!
	Its constructor needs all the financial parameters defining the product:
	- the start date
	- the end date
	- schedule information (frequency, daycounters, calendars)
	- the underlying information 
	- the payment currency
	- payoff information (strike, level and notional)
	- callability information
 */
class CallableQuantoSpreadOptionCreator 
{
public:
	CallableQuantoSpreadOptionCreator(const ARM_Date& startDate, const ARM_Date& endDate); 

	void setScheduleInformation(const ScheduleInformation& scheduleInfo); 
	
	void setUnderlyings(const ARM_Currency& cpnCcy, 
						ARM_INDEX_TYPE shortIndex,
						ARM_INDEX_TYPE longIndex,
						const ARM_Currency& fundCcy,
						ARM_INDEX_TYPE fundIndex); 

	void setPayoffInformation(const ARM_Curve&	itsNominal,
							const ARM_Curve&	itsCpnStrikes,
							const ARM_Curve&	itsCpnLeverageShort,
							const ARM_Curve&	itsCpnLeverageLong,
							const ARM_Curve&	itsCpnMin,
							const ARM_Curve&	itsCpnMax,
							const ARM_Curve&	itsFundMargin); 

	void setCallabilityInformation(int				exerciseFreq,
								int					noticeGap,
								const ARM_Curve&	fees); 
								
	
	ARM_GenSecurityPtr create() const;
	
protected:
	
	static const string CQSOColNamesTable [] ;


	// Schedule info
	ARM_Date			itsStartDate;
	ARM_Date			itsEndDate; 
	ScheduleInformation itsScheduleInformation; 

	// Underlying info
	ARM_Currency		itsCpnCcy;
	ARM_INDEX_TYPE		itsShortIndex;
	ARM_INDEX_TYPE		itsLongIndex;
	ARM_Currency		itsFundCcy;
	ARM_INDEX_TYPE		itsFundIndex;

	// Payoff info
	ARM_Curve			itsNominal;
	ARM_Curve			itsCpnStrikes;
	ARM_Curve			itsCpnLeverageShort;
	ARM_Curve			itsCpnLeverageLong;
	ARM_Curve			itsCpnMin;
	ARM_Curve			itsCpnMax;
	ARM_Curve			itsFundMargin;
	
	// Callability info
	ARM_Curve			itsFees;

}; 


CC_END_NAMESPACE()

#endif