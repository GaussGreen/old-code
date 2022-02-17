/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file holidays.cpp
 *
 *  \brief file for the holidays object
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#if defined(__USE_BASE_LIBRARY)

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/env.h"
#include "gpbase/holidays.h"
#include "gpbase/stringmanip.h"

/// ARM Kernel
#include <glob/dates.h>
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////
///	Class  : ARM_Holidays
///	Routine: Constructor
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_Holidays::ARM_Holidays(const string& vacationsFile)
{
	FILE* vacFile = fopen(vacationsFile.c_str(), "r");
	if( !vacFile )
		throw Exception(__LINE__, __FILE__, ERR_PB_OPEN_FILE, 
			ARM_USERNAME + ": Vacations File " + vacationsFile + " Not found");

	char isoName[4];
	char ccyGroup[4];
	int Day, Month, Year;
	ARM_GP_Vector vacations;
	vacations.reserve(1000);

	int rc = fscanf(vacFile, "%3s%*[ ]%4d%*c%2d%*c%2d%*[^\n]%[\n]", isoName, &Year, &Month, &Day );
	while ( rc != EOF )
	{
		vacations.resize(0);
        strcpy(ccyGroup, isoName);
		
		while( strcmp(ccyGroup,isoName)==0 && (rc!=EOF) )
		{
			ARM_Date vacDay(Day, Month, Year);
			double julDate = vacDay.GetJulian();

			vacations.push_back(julDate);
			rc = fscanf(vacFile, "%3s%*[ ]%4d%*c%2d%*c%2d%*[^\n]%[\n]", isoName, &Year, &Month, &Day );
		}
		
		string ccyGroupString(ccyGroup);
		ARM_VacationData data( vacations, false );
		if( !itsData.insert( CC_NS(std,pair)< const string, ARM_VacationData >( ccyGroupString, data ) ).second )
			throw Exception(__LINE__, __FILE__, ERR_PB_COUNTRY_VACATIONS, 
				" Problem inserting calendar " + ccyGroupString );
	}
	
	fclose(vacFile);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_Holidays
///	Routine: Is a business day?
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

bool ARM_Holidays::IsBusinessDay(const string& isoName, double julianDate) const
{
	DataMap::const_iterator iter = itsData.find( isoName );
	
	/// no calendar?
	/// just return whether a date is a Week End Day
	if( iter == itsData.end() )
	{
	    ARM_Date tmpDate(julianDate);
		return !tmpDate.IsWeekEndDay();
	}
	else
	{
	    if( ((*iter).second).itsIncludeWeekEnd )
		{
			ARM_Date tmpDate(julianDate);
			if( tmpDate.IsWeekEndDay() )
				return false;
		}
		
		ARM_GP_Vector::const_iterator dateIter = CC_NS(std,lower_bound)( (*iter).second.itsVacationsDates.begin(), 
			(*iter).second.itsVacationsDates.end(), julianDate );
		return dateIter == (*iter).second.itsVacationsDates.end() || *dateIter != julianDate;
	}
}


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/