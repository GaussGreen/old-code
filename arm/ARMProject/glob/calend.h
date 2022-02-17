/*
 * $Log: calend.h,v $
 * Revision 1.3  2003/08/06 15:45:20  ebenhamou
 * added inflation calendar.. no week end
 *
 * Revision 1.2  2003/07/01 19:42:28  jpriaudel
 * suppression of include armglob
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : calend.h                                                     */
/*                                                                            */
/*                                                                            */
/* DESCRIPTION : Vacation days Calendar encapsulation                         */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _CALEND_H
#define _CALEND_H
 


#include "containr.h"
#include "tree.h"



#define ARM_ISO_NAME_LENGTH		4



struct CountryVacations : public ARM_Object
{
    char    ISOName[ARM_ISO_NAME_LENGTH];

    int     size;

    ARM_ValueTable* julianVacationDays;




    CountryVacations(void)
    {
        MEMSET(ISOName, 0, sizeof(ISOName)); 

        size = 0;

        julianVacationDays = (ARM_ValueTable *) NULL;
    }

    CountryVacations(char* isoName, int sz, double* julVacDays)
    {
        strcpy(ISOName, isoName);

        size = sz;

        julianVacationDays = new ARM_ValueTable(new ARM_Val(julVacDays[0]));


        for (int i = 1; i < size; i++)
        {
            ARM_Val* newVal = new ARM_Val(julVacDays[i]);

            julianVacationDays = julianVacationDays->InsertValue(newVal);
        }
    }

   ~CountryVacations(void)
    {
        if (julianVacationDays)
        {
           delete julianVacationDays;
        }
 
        size = 0;

        julianVacationDays = (ARM_ValueTable *) NULL;

        MEMSET(ISOName, 0, sizeof(ISOName)); 
     }
};



class ARM_Date;

class ARM_Calendar: public ARM_Object
{
    private:
 
        ARM_Container* vacationDays;


        void SetCountryVacationDays(CountryVacations* countryVacDays);

        CountryVacations* GetVacationDaysByCountry(char* isoName);

        void Init(void)
        {
            vacationDays = NULL;
        }

		bool IncludeWeekEndCal(char* isoName );

    public:
 
        ARM_Calendar(void)
        {
            Init();
        }

       ~ARM_Calendar(void)
        {
            if (vacationDays)
            {
               vacationDays->Destroy();
            }

            vacationDays = (ARM_Container *) NULL;
        }

        ARM_Calendar(ARM_Container* vacDays)
        {
            vacationDays = vacDays;
        }

        ARM_Calendar(char* vacationsFile,bool isString = false);

        int IsBusinessDay(char* isoName, double julianDate);
};











#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
