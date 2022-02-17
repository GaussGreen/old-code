/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 * $Log: calend.cpp,v $
 * Revision 1.8  2003/10/16 16:59:41  ebenhamou
 * make the table const for safetyness
 *
 * Revision 1.7  2003/10/16 16:10:19  ebenhamou
 * remove dotNet bug
 *
 * Revision 1.6  2003/08/22 12:43:31  ebenhamou
 * just to make purify happy
 *
 * Revision 1.5  2003/08/15 16:25:25  ebenhamou
 * case insensitive for weekend calendar
 *
 * Revision 1.4  2003/08/06 15:49:47  ebenhamou
 * added log
 *
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : calend.cpp                                                   */
/*                                                                            */
/*                                                                            */
/* DESCRIPTION : Vacation days Calendar encapsulation                         */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#include <string.h>
#include "expt.h"
#include "dates.h"
#include "calend.h"


/// add here the noWeekEndCalendars
/// case insensitive!
const char* noWeekEndCalendar[] = 
{
        "INF",
};

static int NOWEEKENDCALNB = sizeof( noWeekEndCalendar ) / sizeof( noWeekEndCalendar[ 0 ] );


ARM_Calendar::ARM_Calendar(char* vacationsFile,bool isString)
{
	Init();

	if (isString == false)
	{
		FILE* vacFile = NULL;
		CountryVacations* countryVac;
		char isoName[ARM_ISO_NAME_LENGTH];
		char ccyGroup[ARM_ISO_NAME_LENGTH];
		double vacations[6000];
		int  rc; 
   
		vacFile = fopen(vacationsFile, "r");

		if ( vacFile == NULL )
		{
		   throw Exception(__LINE__, __FILE__, ERR_PB_OPEN_FILE,
						   "==> Vacations File Not found");
		}    

		int Day, Month, Year;

		char ignoreChar[1];

		rc = fscanf(vacFile, "%3s", isoName );
		rc = fscanf(vacFile, "%*[ ]", ignoreChar );
		rc = fscanf(vacFile, "%4d", &Year );
		rc = fscanf(vacFile, "%*c", ignoreChar );
		rc = fscanf(vacFile, "%2d", &Month );
		rc = fscanf(vacFile, "%*c", ignoreChar );
		rc = fscanf(vacFile, "%2d", &Day);
		rc = fscanf(vacFile, "%*[^\n]%[\n]", ignoreChar );
                     

		while ( rc != EOF )
		{
			int sz = 0;

			strcpy(ccyGroup, isoName);

			while (( strcmp(ccyGroup, isoName) == 0 ) && ( rc != EOF ))
			{
				ARM_Date vacDay(Day, Month, Year);
				double julDate = vacDay.GetJulian();
 
				vacations[sz++] = julDate;

				rc = fscanf(vacFile, "%3s", isoName );
				rc = fscanf(vacFile, "%*[ ]", ignoreChar );
				rc = fscanf(vacFile, "%4d", &Year );
				rc = fscanf(vacFile, "%*c", ignoreChar );
				rc = fscanf(vacFile, "%2d", &Month );
				rc = fscanf(vacFile, "%*c", ignoreChar );
				rc = fscanf(vacFile, "%2d", &Day);
				rc = fscanf(vacFile, "%*[^\n]%[\n]", ignoreChar );
			}

			countryVac = new CountryVacations(ccyGroup, sz, vacations);
			SetCountryVacationDays(countryVac);
		}

		fclose(vacFile);
	}
	else
	{
		CountryVacations* countryVac;
		char isoName[ARM_ISO_NAME_LENGTH];
		char ccyGroup[ARM_ISO_NAME_LENGTH];
		double vacations[6000];
		int  rc; 

		if ( vacationsFile == NULL )
		{
		   throw Exception(__LINE__, __FILE__, ERR_PB_OPEN_FILE,
						   "==> Vacations Not found");
		}    

		int Day, Month, Year;

		char ignoreChar[1];
		char* buf;

		rc = sscanf(vacationsFile, "%3s", isoName );
		buf = vacationsFile+3;
		rc = sscanf(buf, "%*[ ]", ignoreChar );
		buf++;
		rc = sscanf(buf, "%4d", &Year );
		buf = buf+4;
		rc = sscanf(buf, "%*c", ignoreChar );
		buf++;
		rc = sscanf(buf, "%2d", &Month );
		buf = buf+2;
		rc = sscanf(buf, "%*c", ignoreChar );
		buf++;
		rc = sscanf(buf, "%2d", &Day);
		buf = buf+2;
		rc = sscanf(buf, "%*[^\n]%[\n]", ignoreChar );
		buf++;
                     

		while ( rc != EOF )
		{
			int sz = 0;

			strcpy(ccyGroup, isoName);

			while (( strcmp(ccyGroup, isoName) == 0 ) && ( rc != EOF ))
			{
				ARM_Date vacDay(Day, Month, Year);
				double julDate = vacDay.GetJulian();
 
				vacations[sz++] = julDate;

				rc = sscanf(buf, "%3s", isoName );
				if (rc != EOF)
				{
					buf = buf+3;
					rc = sscanf(buf, "%*[ ]", ignoreChar );
					buf++;
					rc = sscanf(buf, "%4d", &Year );
					buf = buf+4;
					rc = sscanf(buf, "%*c", ignoreChar );
					buf++;
					rc = sscanf(buf, "%2d", &Month );
					buf = buf+2;
					rc = sscanf(buf, "%*c", ignoreChar );
					buf++;
					rc = sscanf(buf, "%2d", &Day);
					buf = buf+2;
					rc = sscanf(buf, "%*[^\n]%[\n]", ignoreChar );
					buf++;
				}
			}

			countryVac = new CountryVacations(ccyGroup, sz, vacations);
			SetCountryVacationDays(countryVac);
		}
	}
}



CountryVacations* ARM_Calendar::GetVacationDaysByCountry(char* isoName)
{
    if ( vacationDays == NULL )
    {
       return((CountryVacations *) NULL);
    }
            
    CountryVacations* curCountryVac = NULL;

    curCountryVac = (CountryVacations *) vacationDays->Start();

    int found = 0;

    while ( curCountryVac && !(found) )
    {
        if ( strcmp(curCountryVac->ISOName, isoName) == 0 )
        { 
           found = 1;
        }
        else
        {
           curCountryVac = (CountryVacations *) vacationDays->Next();
        }
    }

    if (found)
    {
       return(curCountryVac);
    }
    else
    {
       return((CountryVacations *) NULL);
    }
}



void ARM_Calendar::SetCountryVacationDays(CountryVacations* countryVacDays)
{
    if ( countryVacDays == NULL )
    {
       return;
    }

    if ( vacationDays == NULL )
    {
       vacationDays = new ARM_Container();
    }

    vacationDays->Append(countryVacDays);
}



int TestIsBusinessDay(CountryVacations* countryVac, char* isoName, 
                      double julianDate)
{
    int i;
    int sz = countryVac->size;
 
    
    i = 0;

    int found = 0;
 
    ARM_Val value = (ARM_Val) julianDate;

    if (countryVac->julianVacationDays->SearchValue(&value))
    {
       found = 1;
    }

    return(!(found));
}



////////////////////////////////////////////////////////
//// checks whether weekends are holidays for a calendar
//// if weekends are holidays return true
//// otherwise return false
//// the list of the special calendars is defined
//// in noWeekEndCalendar
////////////////////////////////////////////////////////

bool ARM_Calendar::IncludeWeekEndCal(char* isoName )
{
    int i;
    char bufname[50];
    strcpy(bufname,isoName);

    (void) ARM_UPPERCASE(bufname);

    for (i = 0; i < NOWEEKENDCALNB; ++ i)
        if ( strcmp( bufname, noWeekEndCalendar[i] ) == 0 )
           return false;

    /// otherCase
    return true;
}



int ARM_Calendar::IsBusinessDay(char* isoName, double julianDate)
{
    ARM_Date tmpDate(julianDate);

    CountryVacations* countryVac = NULL;

    //// checks for specific type calendar
    if ( IncludeWeekEndCal( isoName )
        && tmpDate.IsWeekEndDay() )
       return(0);

    countryVac = GetVacationDaysByCountry(isoName);

    if ( countryVac == NULL )
    {
       return(1);
    } 

    int yesOrNo = TestIsBusinessDay(countryVac, isoName, julianDate);

    return(yesOrNo);
}
 



#undef NOWEEKENDCALNB

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
