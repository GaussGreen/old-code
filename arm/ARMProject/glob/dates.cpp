/*
 * $Log: dates.cpp,v $
 * Revision 1.33  2004/03/23 13:58:12  ebenhamou
 * rremove ugly this->!
 *
 * Revision 1.32  2003/11/27 14:16:47  ebenhamou
 * make function returns real throw instead of buggy 3 month!
 *
 * Revision 1.30  2003/11/13 11:14:13  emezzine
 * correct a bug in AddPeriod
 *
 * Revision 1.29  2003/11/13 10:35:13  ebenhamou
 * bug correction for negative period
 *
 * Revision 1.27  2003/11/07 11:39:31  ebenhamou
 * added no base
 *
 * Revision 1.26  2003/11/06 13:11:05  ebenhamou
 * change the toString
 *
 * Revision 1.25  2003/10/13 10:19:52  emezzine
 * rset dates.cpp
 * Added "CptEndDates" and "CptInterestTerms" to calculate swaprate in model
 *
 * Revision 1.24  2003/09/26 07:39:04  ebenhamou
 * move iomanip header also to the non unix part since not used
 *
 * Revision 1.22  2003/09/02 12:41:59  jpriaudel
 * modif du format de la date en string
 *
 * Revision 1.21  2003/09/01 12:55:47  jpriaudel
 * ajout de JulianToStrDateDay
 * + modif dans Previous et NextRegularDate
 *
 * Revision 1.20  2003/08/15 07:01:29  ebenhamou
 * added CountYearsWithoutException
 *
 * Revision 1.19  2003/07/31 06:45:08  ebenhamou
 * put if for clarity reason instead of ternaire operator
 *
 * Revision 1.18  2003/07/29 08:45:43  ebenhamou
 * cas ZEROCOUPON
 *
 * Revision 1.17  2003/07/22 13:59:07  mab
 * const supressed in operator -
 *
 * Revision 1.16  2003/07/22 13:35:55  mab
 * const suppressed in :
 * double operator - ( const ARM_Date& d)
 *
 * Revision 1.15  2003/07/17 06:53:42  ebenhamou
 * added DCGetName for view method
 *
 * Revision 1.14  2003/07/16 06:55:55  ebenhamou
 * addPeriod support period with string
 *
 * Revision 1.13  2003/07/01 20:14:23  jpriaudel
 * correction
 *
 * Revision 1.12  2003/04/01 17:19:12  jpriaudel
 * Ajout des fonctions DaysBetweenDatesWithoutException
 * et CountYearsWithoutException qui permettent d'avoir
 * des resultats negatifs
 *
 * Revision 1.11  2003/03/11 16:26:58  mab
 * Added : parameter adjFirstdate
 *  in CptStartDates
 *
 * Revision 1.10  2002/11/15 10:12:15  mab
 * SysToday(void) added
 *
 * Revision 1.9  2002/03/01 14:40:23  mab
 * ADDED : ARM_Date& ARM_Date::AddPeriod(int period, int cday, char *ccy)
 *
 * Revision 1.8  2002/01/17 13:59:23  mab
 * Rajout de : if (( frequency != K_DAILY ) && ( frequency != K_WEEKLY )
 * ds CptStartDates
 *
 * Revision 1.7  2001/12/28 14:33:13  arm
 * Generalisation de : IsBusinessDay
 * pour les CROSS Calendars
 * Rajout de : AddPeriodMult
 *
 */


/*-----------------------------------------------------------------------*/
/*                                                                       */
/* FILE        : dates.cpp                                               */
/*                                                                       */
/* DESCRIPTION : Methods of date Object                                  */
/*                                                                       */
/* DATE        : Wed Apr 17 1996                                         */
/*                                                                       */
/*-----------------------------------------------------------------------*/




#include "expt.h"

#include <string>
using std::string;

#define SIZEOFSHORTWEEKDAYS 3
                                                   

#include "dates.h" 

#include "calend.h"



/*-----------------------------------------------------------------------*/


ARM_Calendar* ARM_Date::vacationsCalendar = NULL;



int MonthDays[MAXMONTH] = { 31, 28, 31, 30, 31, 30, 31, 31, 
                                    30, 31, 30, 31 };

char ShortMonths[MAXMONTH][4] = {"Jan", "Feb", "Mar", "Apr", "May", 
                                 "Jun", "Jul", "Aug", "Sep", "Oct", 
                                 "Nov", "Dec"};

char ShortMonthsCapitalLetter[MAXMONTH][4] = 
                                {"JAN", "FEB", "MAR", "APR", "MAY", 
                                 "JUN", "JUL", "AUG", "SEP", "OCT", 
                                 "NOV", "DEC"};

static char LongMonths[MAXMONTH][10] = {"January", "February", "March", 
                                        "April", "May", "June", 
                                        "July", "August", "September", 
                                        "October", "November", "December"};

static char ShortWeekDays[MAXWEEKDAY+1][4] = {"Sun", "Mon", "Tue", "Wed", 
                                            "Thu", "Fri", "Sat"};

static char LongWeekDays[MAXWEEKDAY+1][10] = {"Sunday", "Monday", "Tuesday", 
                                              "Wednesday", "Thursday", 
                                                "Friday", "Saturday"};





/*-----------------------------------------------------------------------*/
/* Utlities functions                                                    */
/*-----------------------------------------------------------------------*/


void ARM_Date::InitCalendar(char* calFile, bool isString)
{
    ARM_Calendar* cal = new ARM_Calendar(calFile,isString);

    SetCalendar(cal);
}

void ARM_Date::SetCalendar(ARM_Calendar* cal)
{
    ARM_Date::vacationsCalendar = cal;
}

ARM_Calendar* ARM_Date::GetCalendar(void)
{
    return(ARM_Date::vacationsCalendar);
}



/*----------------------------------------------------------------------------*/
/* Methods                                                                    */
/*----------------------------------------------------------------------------*/

char* ARM_Date::GetLongDayName(int d)
{
    return(LongWeekDays[d]);
}



char* ARM_Date::GetLongMonthName(int m)
{
    return(LongMonths[m-1]);
}



char* ARM_Date::GetShortDayName(int d)
{
    return(ShortWeekDays[d]);
}



char* ARM_Date::GetShortMonthName(int m)
{
     return(ShortMonths[m-1]);
}



/*----------------------------------------------------------------------------*/
/* Days Number in a given month                                               */
/*----------------------------------------------------------------------------*/

int ARM_Date::GetDaysInMonth(int m, int y)
{
    return(DaysInMonth(m, y));
}




/*----------------------------------------------------------------------------*/
/*   Friend Function Operators                                                */
/*----------------------------------------------------------------------------*/


ARM_Date& operator + (ARM_Date& d, long x)
{
    ARM_Date sum;


    sum = d;
    sum.Julian += x;
 
    if (( sum.Julian < MINDATE ) || ( sum.Julian > MAXDATE ))
    {
       sum.Julian = BADDATE;

       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
             "Invalid Date");

       throw x;

       d = sum;

       return(d);
    }
    else
    {
       JulianToDMY(sum.Julian, &(sum.Day), &(sum.Month), &(sum.Year));

       sum.DayOfWeek = JulianToDayOfWeek(sum.Julian);
    }

    d = sum;

    return(d);
}



ARM_Date& operator + (long x, ARM_Date &d)
{
    return(d + x);
}






/*----------------------------------------------------------------------------*/
/* Operators methods                                                          */
/*----------------------------------------------------------------------------*/




double ARM_Date::operator - (const ARM_Date& d) const 
{
    return(Julian - d.Julian);
}


ARM_Date& ARM_Date::operator - (long x)
{ 
    return(*this + (-x));
}


ARM_Date& ARM_Date::operator ++ (void)
{
    ARM_Date d;


    d =  (*this) + 1L;

    (*this) = d;

    return(*this);
}



ARM_Date& ARM_Date::operator -- (void)
{
    ARM_Date d;
 
 
    d =  (*this) + (-1L);
 
    (*this) = d;
 
    return(*this);
}



ARM_Date& ARM_Date::operator += (long x)
{
    ARM_Date d;
 
 
    d =  (*this) + x;
 
    (*this) = d;
 
    return(*this);
}



ARM_Date& ARM_Date::operator -= (long x)
{
    ARM_Date d;
 
 
    d =  (*this) + (-x);
 
    (*this) = d;

    return(*this);
}



/*----------------------------------------------------------------------------*/
/* Other utilities                                                            */
/*----------------------------------------------------------------------------*/


ARM_Date& ARM_Date::Paque(int cYear)
{
    ARM_Date ItIsPaque;
 

    //****   Gestion des erreurs de la procedure entre 1900 et 2099

    if ( cYear == 1954 )
    {
       ItIsPaque.ChangeDate(18,04,1954);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }
 
    if ( cYear == 1981 )
    {
       ItIsPaque.ChangeDate(19,04,1981);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }
 
    if ( cYear == 2049 )
    {
       ItIsPaque.ChangeDate(18,04,2049);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }
 
    if ( cYear == 2076 )
    {
       ItIsPaque.ChangeDate(19,04,2076);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }

    // Fin de la  gestion des erreurs de la procedure entre 1900 et 2099
 
    int a1,a2,a3,a4,a5;
 
    a1 = cYear % 19;
 
    a2 = cYear % 4 ;
 
    a3 = cYear % 7 ;
 
    a4 = a1*19+24 ;
 
    a4 = a4 % 30 ;
 
    a5 = ((2*a2) + (4*a3) + (6*a4) + 5);
 
    a5 = a5 % 7;
 
    if ( (a4+a5) <= 9 )
    {
       ItIsPaque.ChangeDate(a4+a5+22,03,cYear);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }
    else
    {
       ItIsPaque.ChangeDate(a4+a5-9,04,cYear);
 
       (*this) = ItIsPaque;
 
       return(*this);
    }
}



int ARM_Date::IsBusinessDay(char* Country)
{
    char curCal[4];
    int sz = strlen(Country);


    curCal[3] = '\0';

    if ( sz > 3 )
    {
       if ( sz % (int(3)) != 0 )
       {
          throw Exception(__LINE__, __FILE__, ERR_PB_COUNTRY_VACATIONS,
           "ARM_Date::IsBusinessDay: Cross Calendar not valid");
       }

       int i = 0;

       int boolRes = 1;

       while ( i < sz ) 
       {
           strncpy(curCal, &Country[i], 3); 


           boolRes *= IsBusinessDay(curCal);


           if (!(boolRes))
           {
              return(0);
           }

           i += 3;
       }

       return(boolRes);
    }

    /*------ Check Calendar first -----*/

    if (vacationsCalendar)
    {
       int predic;

       predic = vacationsCalendar->IsBusinessDay(Country,
                                                 this->GetJulian());

       return(predic);
    }
 
    if (( DayOfWeek == 0 ) || ( DayOfWeek == 6 )) 
       return(0); // WeekEnd

    if (( Day == 1 ) && ( Month == 1 )) 
       return(0); // 1er Janvier
 
 
    if ( !strcmp("FRF",Country) )
    {
       if (( Day == 1 ) && ( Month == 11 ))
          return(0); // 1er Nov. Toussaint

       if (( Day == 1 ) && ( Month == 5 )) 
          return(0); // 1er Mai
 
       if (( Day == 8 ) && ( Month == 5 )) 
          return(0); // 8 Mai
 
       if (( Day == 14 ) && ( Month == 7 )) 
          return(0); // 14 Juillet
 
       if (( Day == 15 ) && ( Month == 8 )) 
          return(0); // 15 Aout
 
       if (( Day == 11 ) && ( Month == 11 )) 
          return(0); // Armistice
 
       if (( Day == 25 ) && ( Month == 12 )) 
          return(0); // Noel
 
       if (( Month > 2 ) || ( Month < 7 ))
       {
          ARM_Date Mobile;

          Mobile.Paque(Year);
          Mobile.AddDays(-2);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(1);  // Vendredi saint
 
          Mobile.AddDays(+3);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Lundi de Paque
 
          Mobile.AddDays(38);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Ascension
 
          Mobile.AddDays(11);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month))
             return(0);  // Lundi de Pentecote
       }

       return(1);
    }
    else if  ( !strcmp("USD",Country) )
    {
       // Jours feries en janvier

       if ( Month == 1 )
       {
          if (( Day == 2 )  && ( DayOfWeek == 1 )) 
             return (0);

          if (( DayOfWeek == 1) && (( Day > 14 ) && ( Day < 22 ))) 
             return (0);

          return (1);
       }

       // Jours feries en fevrier
       if ( Month == 2 )
       {
          if ( (DayOfWeek == 1) && ( (Day>14) && (Day<22) ) ) 
             return (0);

          return (1);
       }

       //  Jours feries en Mai

       if ( Month == 5 )
       {
          if ((DayOfWeek == 1) && ( Day > 24)) 
             return (0);

          return (1);
       }
 
       //  Jours feries en Juillet

       if ( Month == 7 )
       {
          if ( Day == 4 ) 
             return (0);

          if (( Day == 3 ) && ( DayOfWeek == 5 ) ) 
             return (0);

          if (( Day == 5 ) && ( DayOfWeek == 1 )) 
             return (0);

          return(1);
       }
 
       //  Jours feries en septembre

       if ( Month == 9 )
       {
          if (( DayOfWeek == 1 ) && ( Day < 8 )) 
             return (0);

          return (1);
       }
 
       // Jours feries en Octobre

       if ( Month == 10 )
       {
          if (( DayOfWeek == 1) && (( Day > 7 ) && ( Day < 15 )))
             return(0);

          return (1);
       }

       //  Jours feries en Novembre

       if ( Month == 11 )
       {
          if ( Day == 11 ) 
             return (0);

          if (( Day == 10 ) && ( DayOfWeek == 5 ) ) 
             return(0);

          if (( Day == 12 ) && ( DayOfWeek == 1 )) 
             return (0);

          return (1);
       }
 
       if ( Month == 12 )
       {
          if ( Day == 25 ) 
             return (0);

          if (( Day == 26 ) && ( DayOfWeek == 1 )) 
             return (0);

          return (1);
       }
    } // Fin USA
    else if   ( !strcmp("DEM",Country) )
    {  
       ARM_Date Mobile;
 
       if (( Day == 1 ) && ( Month == 5 )) 
          return(0); // 1er Mai
 
       if (( Day == 5 ) && ( Month == 10 )) 
          return(0); // 3 Octobre
 
       if (( Day == 1 ) && ( Month == 11 )) 
          return(0); // 1er Novembre
 
       if (( Day == 24 ) && ( Month == 12 )) 
          return(0); // 24 Decembre

       if (( Day == 25 ) && ( Month == 12 )) 
          return(0); // 25 Decembre

       if (( Day == 26 ) && ( Month == 12 )) 
          return(0); // 26 Decembre
 
       if (( Month > 2 ) || ( Month < 7 ))
       {
          Mobile.Paque(Year);
          Mobile.AddDays(-2);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Vendredi saint
 
          Mobile.AddDays(+3);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Lundi de Paque
 
          Mobile.AddDays(38);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Ascension
  
          Mobile.AddDays(11);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month))
             return(0);  // Lundi de Pentecote
 
          Mobile.AddDays(10);

          if ( ( Day == Mobile.Day ) && ( Month == Mobile.Month) )
             return(0);  // ??????????????????
       }

       return(1);
    }
    else if ( !strcmp("CHF",Country) )
    {
       ARM_Date Mobile;

       if (( Day == 2 ) && ( Month == 1 )) 
          return(0); // 2 Janvier

       if (( Day == 1 ) && ( Month == 8 )) 
          return(0); // 1er Aout

       if (( Day == 25 ) && ( Month == 12 )) 
          return(0); // 25 Decembre

       if (( Day == 26 ) && ( Month == 12 )) 
          return(0); // 26 Decembre
 
       if (( Month > 2 ) || ( Month < 7 ))
       {
          Mobile.Paque(Year);

          Mobile.AddDays(-2);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Vendredi saint
 
          Mobile.AddDays(+3);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Lundi de Paque
 
          Mobile.AddDays(38);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Ascension
 
          Mobile.AddDays(11);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month))
             return(0);  // Lundi de Pentecote
       }

       return (1);
    }
    else if ( !strcmp("XEU",Country) || !(strcmp("ECU",Country)))
    {
       ARM_Date Mobile;
 

       if (( Day == 1 ) && ( Month == 5 )) 
          return(0); // 1er Mai

       if (( Day == 15 ) && ( Month == 8 )) 
          return(0); // 15 Aout

       if (( Day == 1 ) && ( Month == 11 )) 
          return(0); // 1er Novembre

       if (( Day == 25 ) && ( Month == 12 )) 
          return(0); // 24 Decembre

       if (( Day == 26 ) && ( Month == 12 )) 
          return(0); // 24 Decembre
 
       if (( Month > 2 ) || ( Month < 7 ))
       {
          Mobile.Paque(Year);
          Mobile.AddDays(-2);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Vendredi saint
 
          Mobile.AddDays(+3);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Lundi de Paque
 
          Mobile.AddDays(38);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month ))
             return(0);  // Ascension
 
          Mobile.AddDays(11);

          if (( Day == Mobile.Day ) && ( Month == Mobile.Month))
             return(0);  // Lundi de Pentecote
       }

       return (1);
    }
    else if ( !strcmp("JPY",Country) )
    {
       ARM_Date Mobile;


       if (( Day == 2 ) && ( Month == 1 )) 
          return(0); // 2 Janvier

       if (( Day == 3 ) && ( Month == 1 )) 
          return(0); // 3 Janvier
 
       if (( Day == 15 ) && ( Month == 1 )) 
          return(0); //  15 Fevrier

       if (( Day == 16 ) && ( Month == 1 ) && (DayOfWeek==1)) 
          return(0); // 15 Janvier
 
       if (( Day == 11 ) && ( Month == 2 )) 
          return(0); //  11 Fevrier

       if (( Day == 12 ) && ( Month == 2 ) && (DayOfWeek==1)) 
          return(0); // 11 Fevrier
 
       if (( Day == 21 ) && ( Month == 3 )) 
          return(0); //  21 Mars

       if (( Day == 22 ) && ( Month == 3 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 29 ) && ( Month == 4 )) 
          return(0); //  29 Avril

       if (( Day == 30 ) && ( Month == 4 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 3 ) && ( Month == 5 )) 
          return(0); //  3 Mai

       if (( Day == 4 ) && ( Month == 5 )) 
          return(0); //  4 Mai

       if (( Day == 5 ) && ( Month == 5 )) 
          return(0); //  5 Mai

       if (( Day == 6 ) && ( Month == 5 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 15 ) && ( Month == 9 )) 
          return(0); //  15 Septembre

       if (( Day == 16 ) && ( Month == 9 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 23 ) && ( Month == 9 )) 
          return(0); //  23 Septembre

       if (( Day == 24 ) && ( Month == 9 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 10 ) && ( Month == 10 )) 
          return(0); //  10 Octobre

       if (( Day == 11 ) && ( Month == 10 ) && ( DayOfWeek == 1 )) 
          return(0); //

       if (( Day == 3 ) && ( Month == 11 )) 
          return(0); //  03 Novembre

       if (( Day == 4 ) && ( Month == 11 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
 
       if (( Day == 23 ) && ( Month == 11 )) 
          return(0); //  23 Novembre

       if (( Day == 24 ) && ( Month == 11 ) && ( DayOfWeek == 1 )) 
          return(0); //
 
       if (( Day == 23 ) && ( Month == 12 )) 
          return(0); //  23 Decembre
 
       return (1);
    }
    else if   ( !strcmp("GBP",Country) )
    {
       ARM_Date Mobile;
 

       // Un jour Ferie debut Mai
       // ???????????????????????
 
       if ( Month == 5 )
       {   // Dernier Lundi en Mai

          if ( (DayOfWeek == 1) && (Day>24) ) 
             return (0);

          return (1);
       }

       if ( Month == 8 )
       {   // Dernier Lundi d'Aout

          if ( (DayOfWeek == 1) && (Day>24) ) 
             return (0);

          return (1);
       }
 
       if (( Day == 25 ) && ( Month == 12 )) 
          return(0); // 25 Decembre

       if (( Day == 26 ) && ( Month == 12 )) 
          return(0); // 26 Decembre

       if (( Day == 27 ) && (( Month == 12 ) && (DayOfWeek==1) )) 
          return(0);

       if (( Day == 27 ) && (( Month == 12 ) && (DayOfWeek==2) )) 
          return(0);

       if (( Day == 28 ) && (( Month == 12 ) && (DayOfWeek==2) )) 
          return(0);

 
       if ( (Month > 2) || (Month <7) )
       {
          Mobile.Paque(Year);
          Mobile.AddDays(-2);

          if ( ( Day == Mobile.Day ) && ( Month == Mobile.Month ) )
             return(0);  // Vendredi saint
 
          Mobile.AddDays(+3);

          if ( ( Day == Mobile.Day ) && ( Month == Mobile.Month ) )
             return(0);  // Lundi de Paque
       }

       return(1);
    }

    return(1);  // Ouvert par defaut
}



ARM_Date& ARM_Date::SubstractYears(int n)
{
    Year -= n;



    if ( ( Month == 2 ) // February
         &&
         ( Day == 29 )
       )
    {
       if (!(IsLeapYear()))
       {
          Day = 28;
       }
    }
    else
    {
       if ( Day > MonthDays[Month-1] )
       {
          Day = MonthDays[Month-1];
       }
    } 

    Julian = DMYToJulian();

    return(*this);
}



ARM_Date& ARM_Date::AddHours(int h) 
{ 
    Julian += ((double) h) / 24.0;

    JulianDateToDMY();

    DayOfWeek = JulianToDayOfWeek(Julian);

    return(*this);
}



ARM_Date& ARM_Date::AddDays(int d) 
{
    Julian += (double) d; 

    JulianDateToDMY();

    DayOfWeek = JulianToDayOfWeek(Julian);

    return(*this);
}





ARM_Date& ARM_Date::AddYears(int nYears)
{
    if (IsLeapYear())
    {
        if (Month == 2 && Day == 29)
            Day--;
    }

    Year += nYears;
 
    Julian = DMYToJulian();

    DayOfWeek = JulianToDayOfWeek(Julian);

    return(*this);
}



/*----
    Add one period length of time (K_DAILY = 1 day,
    K_MONTHLY = 1 month, etc..) to the current date 
----*/

ARM_Date& ARM_Date::AddPeriod(int period, char *ccy, int AdjustDaily)
{   
    int freq = abs(period);
    
	if (freq == K_ZEROCOUPON)
		return (*this);

    int signedFreq = 12/period; 
   
    switch (freq) 
    {
        case K_ANNUAL :
        case K_SEMIANNUAL :
        case K_QUARTERLY :
        case K_BIMONTHLY :
        case K_MONTHLY :
        {
            AddMonths(signedFreq);
        };
        break;

        case K_WEEKLY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            AddDays(7*pSign);
        };
        break;

        case K_DAILY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            AddDays(pSign);
            if (AdjustDaily == 1)
            {
                while (!(IsBusinessDay(ccy))) //gere les week-ends
                    AddDays(pSign);
            }
        };
        break;

        default :
        {
            /// tries to find out if there is a way to get back a real date!
            if( K_WEEKLY % freq == 0 )
            {
                int pSign = ( (period >= 0) ? 1 : -1 );
                AddDays(7*pSign*K_WEEKLY/freq);
            }
            else
            {
                if( K_MONTHLY % freq == 0 )
                {
                    int pSign = ( (period >= 0) ? 1 : -1 );
                    AddMonths(K_MONTHLY/freq * pSign);
                }
                else
                {
                    if( K_DAILY % freq == 0 )
                    {
                        int i = 0;
                        int pSign    = ( (period >= 0) ? 1 : -1 );
                        int mult    = K_DAILY/freq;
                        while (i < mult)
                        {
                            AddDays(pSign);
                            while (!(IsBusinessDay(ccy))) //to handle business days!
                                AddDays(pSign);
                            i++;
                        }
                    }
                    else
                    {
                        char msg[250];
                        sprintf( msg, "invalid %d freq", freq );
                        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
                    }
                }
            }
        }
        break;
    }

    return(*this);
}


ARM_Date& ARM_Date::AddPeriodMult(int period, int mult, char *ccy, int GOTO_END_OF_MONTH)
{
    int freq = abs(period);
    
	if (freq == K_ZEROCOUPON)
		return (*this);

    int signedFreq = 12/period; 
   
    switch (freq) 
    {
        case K_ANNUAL :
        case K_SEMIANNUAL :
        case K_QUARTERLY :
        case K_BIMONTHLY :
        case K_MONTHLY :
        {
            AddMonths(signedFreq*mult, GOTO_END_OF_MONTH);
        };
        break;

        case K_WEEKLY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            AddDays(7*pSign*mult);
        };
        break;

        case K_DAILY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            int i = 0;

            while (i < mult)
            {
                this->AddDays(pSign);

                while (!(IsBusinessDay(ccy))) //gere les week-ends
                    AddDays(pSign);

                i++;
            }
        };
        break;

        default :
        {
            /// tries to find out if there is a way to get back a real date!
            if( K_WEEKLY % freq == 0 )
            {
                int pSign = ( (period >= 0) ? 1 : -1 );
                AddDays(7*pSign*K_WEEKLY/freq * mult);
            }
            else
            {
                if( K_MONTHLY % freq == 0 )
                {
                    int pSign = ( (period >= 0) ? 1 : -1 );
                    AddMonths(K_MONTHLY/freq * pSign * mult, GOTO_END_OF_MONTH);
                }
                else
                {
                    if( K_DAILY % freq == 0 )
                    {
                        int i = 0;
                        int pSign    = ( (period >= 0) ? 1 : -1 ) * (mult>0? 1: -1);
                        int total    = K_DAILY/freq*abs(mult);
                        while (i < total)
                        {
                            AddDays(pSign);
                            while (!(IsBusinessDay(ccy))) //to handle business days!
                                AddDays(pSign);
                            i++;
                        }
                    }
                    else
                    {
                        char msg[250];
                        sprintf( msg, "invalid %d freq", freq );
                        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
                    }
                }
            }
        }

        break;
    }
    return(*this);
}

ARM_Date& ARM_Date::AddPeriodMult(int period, int mult, const std::string& ccy, int GOTO_END_OF_MONTH)
{
	return AddPeriodMult(period,mult,(char*)ccy.c_str(), GOTO_END_OF_MONTH); 
}
/*
 * This addins gives the ability to addperiod
 * using string...example are 1y, 20y, 3m, 1w, 1d etc
 */
ARM_Date& ARM_Date::AddPeriod( const string& periodStr, const char* calendar) 
{
    /// parses the periodStr as a number and
    /// a period ... could be yYmMwWdD
    string maturity( "yYMmwWdD" );
    string::size_type pos = periodStr.find_first_of( maturity );
    if( pos == string::npos )
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, string ( "period ") + periodStr 
            + string( " invalid: should be a nb followed by either d/w/m/y case insensitive") );

    string nbString    = periodStr.substr( 0, pos );
    int nb            = atoi( nbString.c_str() );

    // test null period
    if( nb == 0 )
        return *this;

    /// test that is only nb + something of type maturity
    string matuString    = periodStr.substr( pos, string::npos );
    if( matuString.size() != 1 )
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, string( "period " ) + periodStr 
            + string( " invalid: period can be either d/w/m/y while it founds " ) + matuString );
    int period;
    
    /// handle the various maturity type
    switch( matuString[ 0 ] )
    {
        case 'D': case 'd':
            period = K_DAILY;
            break;
        case 'W': case 'w':
            period = K_WEEKLY;
            break;
        case 'M': case 'm':
            period = K_MONTHLY;
            break;
        case 'Y': case 'y':
            period = K_ANNUAL;
            break;
        default:
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                     "Unknown maturity");
    }
    /// because of the non const correctness
    /// we were forced to const_cast

    /// in the case of negative number,
    /// we change the period to be the inverse of the period
    if( nb < 0 )
      {
        period *= -1;
        nb *=-1;
      }
    return AddPeriodMult( period, nb, const_cast<char *>(calendar) );
}



ARM_Date& ARM_Date::AddPeriod(int period, int cday, char *ccy)
{
    int freq = abs(period);

    if (freq == K_ZEROCOUPON)
		return (*this);

    int signedFreq = 12/period; 
   
    switch (freq) 
    {
        case K_ANNUAL :
        case K_SEMIANNUAL :
        case K_QUARTERLY :
        case K_BIMONTHLY :
        case K_MONTHLY :
        {
            this->AddMonths(signedFreq);
            int tday = DaysInMonth(Month,Year);
            if (cday>0 && cday != Day && tday != Day )
            {
                if(tday >= cday)
                    AddDays(cday-Day);
                else
                    AddDays(tday-Day);
            }
        };
        break;

        case K_WEEKLY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );

            this->AddDays(7*pSign);
        };
        break;

        case K_DAILY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            this->AddDays(pSign);
            while (!(this->IsBusinessDay(ccy))) //gere les week-ends
                this->AddDays(pSign);
        };
        break;

        default :
        {
            char msg[250];
            sprintf( msg, "invalid %d freq", freq );
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
        }
        break;
    }
    return(*this);
}

/* FP : new addperiod method */
ARM_Date& ARM_Date::AddPeriodNoAdj(int unit, int length) 
{
	ARM_Date date(*this); 
    switch (unit) 
    {
        case K_ANNUAL :
        case K_SEMIANNUAL :
        case K_QUARTERLY :
        case K_BIMONTHLY :
        case K_MONTHLY :
        {
            return date.AddMonths(12/unit*length, 0);
        };
        break;

        case K_WEEKLY :
        {
            return date.AddDays(7*length);
        };
        break;

        case K_DAILY :
        {
			return date.AddDays(length);
        };
        break;

        default :
        	ARMTHROW(ERR_INVALID_ARGUMENT,"wrong unit in AddPeriodNoAdj()");
        break;
    }
}

ARM_Date& ARM_Date::AddPeriodAdj(int unit, int length, char* ccy, int rollingConvention )
{
	ARM_Date date(*this) ; 
    date = date.AddPeriodNoAdj(unit, length) ;
	return date.AdjustToBusDate(ccy, rollingConvention) ;
}

ARM_Date& ARM_Date::AddPeriodNoAdj(const Period& period)
{
	return AddPeriodNoAdj(period.GetUnit(), period.GetLength()) ;
}
ARM_Date& ARM_Date::AddPeriodAdj(const Period& period, char* ccy, int rollingConvention)
{
	return AddPeriodAdj(period.GetUnit(), period.GetLength(), ccy,rollingConvention) ;
}

/* FP : new addperiod method */

int ARM_Date::IsOrdinaryYear(void)
{
    return(!this->IsLeapYear());
}



int ARM_Date::IsWeekEndDay(void)
{
    return(( DayOfWeek == 0 ) || ( DayOfWeek == 6 ));
}



// Ajoute d jours calendaires puis renvoie le jour ouvres suivant

ARM_Date& ARM_Date::NextNonWeekEndDay(int d)
{
    ARM_Date busDate;


    busDate = (*this);

    busDate.AddDays(d);
 
    while (busDate.IsWeekEndDay())
    {
        (void) busDate.AddDays(1);
    }

    (*this) = busDate;
 
    return(*this);
}

// Roll to d non WeekEnd days Previous(d<0) or Later(d>0) 
// if date is week end day 
// Note : if d = 0 the method returns the same date

// Ajout (ou soustrait) d jours ouvres

ARM_Date& ARM_Date::GapNonWeekEndDay(int d)
{
    if (d == 0)
    {
        return(*this);
    }
    else
    {
        ARM_Date busDate;

        int dayRoll, nbBusDays = 0;

        dayRoll = (d > 0 ? 1 : -1);

        busDate = (*this);

        while (nbBusDays < abs(d))
        {
            (void) busDate.AddDays(dayRoll);
            if (!busDate.IsWeekEndDay()) nbBusDays++;
        }

        (*this) = busDate;
 
        return(*this);
    }
}

// Ajout (ou soustrait) d jours ouvres

ARM_Date& ARM_Date::GapBusinessDay(int d, char *ccy, int fwdRule)
{
    ARM_Date busDate;

    busDate = (*this);
 
    if (d == 0)
    {
        busDate.GoodBusinessDay(fwdRule, ccy);

        (*this) = busDate;
 
        return(*this);
    }
    else
    {
        int dayRoll, nbBusDays = 0;

        dayRoll = (d > 0 ? 1 : -1);

        while (nbBusDays < abs(d))
        {
            (void) busDate.AddDays(dayRoll);
            if (busDate.IsBusinessDay(ccy)) nbBusDays++;
        }

        (*this) = busDate;
 
        return(*this);
    }
}

// Compute Previous(d<0) or Next(d>0) non week end day
// if date is week end day 
// Note : if d = 0 the method returns the same date


ARM_Date& ARM_Date::NonWeekEndDay(int d)
{
    if (d == 0)
    {
        return(*this);
    }
    else
    {
        ARM_Date busDate;

        int dayRoll;

        dayRoll = (d > 0 ? 1 : -1);

        busDate = (*this);
 
        while (busDate.IsWeekEndDay())
        {
            (void) busDate.AddDays(dayRoll);
        }

        (*this) = busDate;
 
        return(*this);
    }
}




ARM_Date& ARM_Date::GoodBusinessDay(int rule, char *ccy)
{

    ARM_Date busDate;

    busDate = (*this);

    if ( rule == 0 )
    {
       return((*this));
    }

    if (this->IsBusinessDay(ccy))
    {
        return(*this);
    }
    else if (abs(rule) == 1)  // Following or Previous
    {
        if ( rule > 0 )
           busDate.NextBusinessDay(ccy);
        else
           busDate.PreviousBusinessDay(ccy);
 
        (*this) = busDate;
 
        return(*this);
    }
    else if (abs(rule) == 2)  // Modify (Following or Previous)
    {
        if (rule > 0)
            busDate.NextBusinessDay(ccy);
        else
            busDate.PreviousBusinessDay(ccy);

        if (busDate.GetMonth() != (*this).GetMonth())
        {
            if (rule > 0)
                busDate.PreviousBusinessDay(ccy);
            else
                busDate.NextBusinessDay(ccy);
        }

        (*this) = busDate;
 
        return(*this);
    }
    else // Default
    {
        busDate.NonWeekEndDay(rule);

        (*this) = busDate;

        return(*this);
    }
}
 


ARM_Date& ARM_Date::GetLivretAFixingDate()
{
    ARM_Date tmpDate;

	if (GetMonth() < 2)
	{
		tmpDate = ARM_Date(1,8,GetYear()-1);
	}
	else if (GetMonth() < 8)
	{
		tmpDate = ARM_Date(1,2,GetYear());
	}
	else
	{
		tmpDate = ARM_Date(1,8,GetYear());
	}

	(*this) = tmpDate;

	return (*this);
}
 
 
 
void ARM_Date::GetTime(int* h, int* mn, int* s)
{
    double t; 
 
 
    t  = Julian - floor(Julian);
 
    t *= 24.0;
 
    *h = (int) floor(t);
 
    t -= (double) (*h);
    t *= 60.0;
    *mn = (int) floor(t);
 
    t -= (double) (*mn);
    t *= 60.0;
    *s = (int) floor(t);
}


void ARM_Date::SysToday(void)
{
    struct tm *t = NULL;
    double julianOfToDay = (double) -1;

    time_t ltime;
    int h, mn, sec;


    if ( t == NULL )
    {
       time(&ltime);

       t = localtime(&ltime);
    }

    Year = t->tm_year + 1900 ;
    Day = t->tm_mday;
    Month = t->tm_mon + 1;

    h = t->tm_hour;
    mn = t->tm_min;
    sec = t->tm_sec;

    DayOfWeek = t->tm_wday;


    if ( julianOfToDay == (double) -1 )
    {
       julianOfToDay = DMYToJulian();

       Julian = julianOfToDay;
    }
    else
    {
       Julian = julianOfToDay;
    }
}


void ARM_Date::Today(void)
{
	static int initialized = 0;
	static struct tm t;
    static double julianOfToDay = (double) -1;

    time_t ltime;
    int h, mn, sec;
 
 
    if ( initialized == 0)
    {
       time(&ltime);

       t = *localtime(&ltime);

	   initialized = 1;
    }
 
    Year = t.tm_year + 1900 ;
    Day = t.tm_mday;
    Month = t.tm_mon + 1;
 
    h = t.tm_hour;
    mn = t.tm_min;
    sec = t.tm_sec;
 
    DayOfWeek = t.tm_wday;
 

    if ( julianOfToDay == (double) -1 )
    {
       julianOfToDay = DMYToJulian();

       Julian = julianOfToDay;
    }
    else
    {
       Julian = julianOfToDay;
    }
}



 
/*----------------------------------------------------------------------------*/
/* Format date                                                                */
/*----------------------------------------------------------------------------*/

int CountLetters(char** s)
{
    int n, c;


    for (n = 0, c = **s; **s == c; n++, (*s)++);

    return(n);
}




void GetToday(char* Today) // Format : DD/MM/YYYY hh:mm:ss
{
     ARM_Date today;
     int h, mn, sec;
     char s[20];


     today.Today();

     today.GetTime(&h, &mn, &sec);
   

     sprintf(s, "%2d%2d%4d %2d:%2d:%2d", today.GetDay(), today.GetMonth(), 
                 today.GetYear(), h, mn, sec);

     strcpy(Today, s);
}



/*----------------------------------------------------------------------------*/
/* Internal function for computing the number of days between date1           */
/*        and date2 using 30/360 day count convention.                        */
/*----------------------------------------------------------------------------*/

static double DaysBetweenDates30_360(ARM_Date& date1, ARM_Date& date2, int meth)
{
    double  nDays;
    int     d1, d2, y1, y2;


    // get dates with int[] format
    y1 = date1.GetYear();
    y2 = date2.GetYear();
 
    //   date1.Dayifies to account for 30/360 rule and compute number of days

    d1 = date1.GetDay();

    d2 = date2.GetDay();

    if ( d1 == 31 ) 
    {
       d1 = 30;
    }

	if (meth == 0) // US convention by default
	{
		if (( d1 == 30 ) && ( d2 == 31 )) 
		{
			d2 = 30;
		}
	}
	else // European convention also called 30_360E
	{
        if (( d1 <= 30 ) && ( d2 == 31 )) 
		{
			d2 = 30;
		}
	}

    nDays = (double) (y2 - y1) * (double) 360
                + (double) ( date2.GetMonth() - date1.GetMonth() ) 
                * (double) 30 + (double) d2 - (double) d1;
 
    // add day fractions

    nDays += date2.GetFracDay() - date1.GetFracDay();
 
    return(nDays);
}



double DaysBetweenDates(int DayCount, ARM_Date& date1, ARM_Date& date2)
{
    double  nDays=0.0;       
 

    if ( date2 < date1 )
    {
       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
             "Invalid Date");
 
       throw x;
    }

    switch (DayCount) 
    {
        case KACTUAL_360 :
        case KACTUAL_365 :
        case KACTUAL_ACTUAL :
        case KACTUAL_FEB29 :
        case KACTUAL_REAL :
        {
            nDays = date2.GetJulian() - date1.GetJulian();
        };
        break;

        case K30_360 :
        {
            nDays = DaysBetweenDates30_360(date1, date2, 0);
        };
        break;

		case K30_360E :
        {
            nDays = DaysBetweenDates30_360(date1, date2, 1);
        };
        break;

        default :
            nDays = -1;
        break;
    }
 
    return(nDays);
}

const char* DCGetName( long DCF )
{
    switch( DCF )
    {
        case    KACTUAL_ACTUAL:        return "ACTUAL_ACTUAL";
        case    KACTUAL_365:        return "ACTUAL_365"; 
        case    KACTUAL_360:        return "ACTUAL_360";
        case    K30_360:            return "30_360"; 
        case    KACTUAL_REAL:       return "ACTUAL_REAL";
        case    KACTUAL_FEB29:        return "ACTUAL_FEB29";
        case    KACTUAL_ISMA:        return "ACTUAL_ISMA";
        case    K30_360E:            return "30_360E"; 
        default:
            throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID, "Invalid Day Count");
    }
}


double DaysBetweenDates(int DayCount, double JulianDate1, double JulianDate2)
{
    double  nDays=0.0;       
 

    if ( JulianDate2 < JulianDate1 )
    {
       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
             "Invalid Dates");
 
       throw x;
    }

    switch (DayCount) 
    {
        case KACTUAL_360 :
        case KACTUAL_365 :
        case KACTUAL_ACTUAL :
        case KACTUAL_FEB29 :
        case KACTUAL_REAL :
        {
            nDays = JulianDate2 - JulianDate1;
        };
        break;

        case K30_360 :
        {
            nDays = DaysBetweenDates30_360(ARM_Date(JulianDate1),
                                           ARM_Date(JulianDate2), 0);
        };
        break;

		case K30_360E :
        {
            nDays = DaysBetweenDates30_360(ARM_Date(JulianDate1),
                                           ARM_Date(JulianDate2), 1);
        };
        break;

        default :
            nDays = -1;
        break;
    }
 
    return(nDays);
}

double DaysBetweenDatesWithoutException(int DayCount, double JulianDate1, double JulianDate2)
{
    double  nDays=0.0;

    switch (DayCount)
    {
        case KACTUAL_360 :
        case KACTUAL_365 :
        case KACTUAL_ACTUAL :
        case KACTUAL_FEB29 :
        case KACTUAL_REAL :
        {
            nDays = JulianDate2 - JulianDate1;
        };
        break;

        case K30_360 :
        {
            nDays = DaysBetweenDates30_360(ARM_Date(JulianDate1),
                                           ARM_Date(JulianDate2), 0);
        };
        break;

		case K30_360E :
        {
            nDays = DaysBetweenDates30_360(ARM_Date(JulianDate1),
                                           ARM_Date(JulianDate2), 1);
        };
        break;

        default :
            nDays = -1;
        break;
    }

    return(nDays);
}

static double baseDates(ARM_Date& d1, ARM_Date& d2)
{
    double base = 365.0;

    int year = 1992;
    
    if (d2.IsLeapYear())
    {
        year = d2.GetYear();
    }
 
    if (d1.IsLeapYear())
    {
        year = d1.GetYear();
    }
    
    ARM_Date feb29(29,2,year);

    if (feb29>=d1 && feb29<=d2)
       base = 366;

    return(base);
}

// - no need to have fromDate < toDate
// - fromDate and toDate have not to be business day
int CountBusinessDays(const ARM_Date& fromDate, const ARM_Date& toDate, char* calendar)
{
    ARM_Date minDate, maxDate, tmpDate;
    int sign, nbDays = 0;

    if (fromDate == toDate)
	{
		return 0;
	}
    else if (fromDate < toDate)
    {
        minDate = fromDate;
		tmpDate = fromDate;
        maxDate = toDate;
        sign = 1;
    }
    else
    {
        minDate = toDate;
		tmpDate = toDate;
        maxDate = fromDate;
        sign = -1;
    }

    while (minDate < maxDate)
    {
        minDate.NextBusinessDay(calendar);
        nbDays++;
    }

	if (tmpDate.IsBusinessDay(calendar) == 0)
		nbDays--;

    return sign*nbDays;
}

double CountYears(int DayCount, ARM_Date& date1, ARM_Date& date2)
{
    double  yt, yt1, yt2;
 
    switch (DayCount) 
    {
        case KNOBASE:
        {
            if( date1 > date2 )
                throw Exception(__LINE__, __FILE__, ERR_PB_COUNTRY_VACATIONS,
                string( "date 1: " ) + date1.toString() + string( "after date 2: ") +
                date2.toString() );
            return 1.0;
        }
        break;

        case K30_360 :
        case K30_360E :
        case KACTUAL_360 :
        {
            yt = DaysBetweenDates(DayCount, date1, date2);

            return(yt / 360.0);
        };
        break;
 
        case KACTUAL_365 :
        {
            yt = DaysBetweenDates(DayCount, date1, date2);
            return(yt / 365.0);
        };
        break;
 
        case KACTUAL_FEB29 :
        {
            yt = DaysBetweenDates(DayCount, date1, date2);

            double base = baseDates(date1, date2);

            return(yt / base);
        };
        break;
 
        case KACTUAL_ACTUAL :
        {
            ARM_Date beginYear1(1,1, date1.Year);
            ARM_Date beginYear2(1,1, date2.Year);

            yt1 = - DaysBetweenDates(DayCount, beginYear1, date1);

            yt2 = DaysBetweenDates(DayCount, beginYear2, date2);
 
            yt1 /=  (date1.IsOrdinaryYear() ? 365.0 : 366.0);
            yt2 /=  (date2.IsOrdinaryYear() ? 365.0 : 366.0);
 
            yt = yt1+yt2 + (double) (date2.Year - date1.Year);
 
            return(yt);
       };
       break;
 
       case KACTUAL_REAL : 
       {
            ARM_Date tmpDate;
            ARM_Date date2Sup(date2.Day, date2.Month, date1.Year);
            ARM_Date date2Inf;

           
            date2Inf = date2Sup;

            if ( date2Sup > date1 )
            {
               (void) date2Inf.SubstractYears(1); 
            }
            else
            {
               (void) date2Sup.AddYears(1);
            }
 
            yt = date2.Year - date2Sup.Year;

            yt += DaysBetweenDates(KACTUAL_ACTUAL, date1, 
                    date2Sup)/DaysBetweenDates(KACTUAL_ACTUAL, 
                                        date2Inf, date2Sup); 
 
            return(yt);
       };
       break;

       default :
           return (1000000);
           //// FIX FIX
           //// should rather be an exception!
           //// throw Exception(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
           ////        "Invalid dayConut");
    }
}

double CountYears(int DayCount, double JulianDate1, double JulianDate2)
{
    double  yt, yt1, yt2;
 
 
    switch (DayCount) 
    {
        case KNOBASE:
        {
            if( JulianDate1 > JulianDate2 )
                throw Exception(__LINE__, __FILE__, ERR_PB_COUNTRY_VACATIONS,
                string( "date 1: " ) + ARM_Date( JulianDate1 ).toString() + string( "after date 2: ") +
                ARM_Date( JulianDate2 ).toString() );
            return 1.0;
        }
        break;

        case KCOUPON:
        {
            if( JulianDate1 > JulianDate2 )
                throw Exception(__LINE__, __FILE__, ERR_PB_COUNTRY_VACATIONS,
                string( "date 1: " ) + ARM_Date( JulianDate1 ).toString() + string( "after date 2: ") +
                ARM_Date( JulianDate2 ).toString() );
            yt = (JulianDate2 - JulianDate1) / 365.;
            if (yt > 0.9) return roundsimple(yt,0);
            else return 1./roundsimple(1./yt,0);
        };
        break;

        case K30_360 :
        case K30_360E :
        case KACTUAL_360 :
        {
            yt = DaysBetweenDates(DayCount, JulianDate1, JulianDate2);

            return(yt / 360.0);
        };
        break;
 
        case KACTUAL_365 :
        {
            yt = DaysBetweenDates(DayCount, JulianDate1, JulianDate2);
            return(yt / 365.0);
        };
        break;
 
        case KACTUAL_FEB29 :
        {
            yt = DaysBetweenDates(DayCount, JulianDate1, JulianDate2);

            double base = baseDates(ARM_Date(JulianDate1),
                                    ARM_Date(JulianDate2));

            return(yt / base);
        };
        break;
 
        case KACTUAL_ACTUAL :
        {
            ARM_Date beginYear1(1,1, ARM_Date(JulianDate1).Year);
            ARM_Date beginYear2(1,1, ARM_Date(JulianDate2).Year);

            yt1 = - DaysBetweenDates(DayCount, beginYear1.GetJulian(), JulianDate1);

            yt2 = DaysBetweenDates(DayCount, beginYear2.GetJulian(), JulianDate2);
 
            yt1 /=  (beginYear1.IsOrdinaryYear() ? 365.0 : 366.0);
            yt2 /=  (beginYear2.IsOrdinaryYear() ? 365.0 : 366.0);
 
            yt = yt1+yt2 + (double) (beginYear2.Year - beginYear1.Year);
 
            return(yt);
       };
       break;
 
       case KACTUAL_REAL : 
       {
            ARM_Date date1 = ARM_Date(JulianDate1);
            ARM_Date date2 = ARM_Date(JulianDate2);
            ARM_Date date2Sup(date2.Day, date2.Month, date1.Year);
            ARM_Date date2Inf;

           
            date2Inf = date2Sup;

            if ( date2Sup > date1 )
            {
               (void) date2Inf.SubstractYears(1); 
            }
            else
            {
               (void) date2Sup.AddYears(1);
            }
 
            yt = date2.Year - date2Sup.Year;

            yt += DaysBetweenDates(KACTUAL_ACTUAL, JulianDate1,
                                   date2Sup.GetJulian())
                 / DaysBetweenDates(KACTUAL_ACTUAL, date2Inf.GetJulian(),
                                    date2Sup.GetJulian()); 

            return(yt);
       };
       break;

       default :
            return(10000000);
    }
}


double CountYearsWithoutException(int DayCount, double JulianDate1, double JulianDate2)
{
    double  yt, yt1, yt2;
 
    switch (DayCount) 
    {
        case KNOBASE:
            return JulianDate1 <= JulianDate2;

        case KCOUPON:
        {
            yt = (JulianDate2 - JulianDate1) / 365.;
            if (yt > 0.9) return roundsimple(yt,0);
            else return 1./roundsimple(1./yt,0);
        };
        break;

        case K30_360 :
        case K30_360E :
        case KACTUAL_360 :
        {
            yt = DaysBetweenDatesWithoutException(DayCount, JulianDate1, JulianDate2);

            return(yt / 360.0);
        };
        break;
 
        case KACTUAL_365 :
        {
            yt = DaysBetweenDatesWithoutException(DayCount, JulianDate1, JulianDate2);
            return(yt / 365.0);
        };
        break;
 
        case KACTUAL_FEB29 :
        {
            yt = DaysBetweenDatesWithoutException(DayCount, JulianDate1, JulianDate2);

            double base = baseDates(ARM_Date(JulianDate1),
                                    ARM_Date(JulianDate2));

            return(yt / base);
        };
        break;
 
        case KACTUAL_ACTUAL :
        {
            ARM_Date beginYear1(1,1, ARM_Date(JulianDate1).Year);
            ARM_Date beginYear2(1,1, ARM_Date(JulianDate2).Year);

            yt1 = - DaysBetweenDatesWithoutException(DayCount, beginYear1.GetJulian(), JulianDate1);

            yt2 = DaysBetweenDatesWithoutException(DayCount, beginYear2.GetJulian(), JulianDate2);
 
            yt1 /=  (beginYear1.IsOrdinaryYear() ? 365.0 : 366.0);
            yt2 /=  (beginYear2.IsOrdinaryYear() ? 365.0 : 366.0);
 
            yt = yt1+yt2 + (double) (beginYear2.Year - beginYear1.Year);
 
            return(yt);
       };
       break;
 
       case KACTUAL_REAL : 
       {
            ARM_Date date1 = ARM_Date(JulianDate1);
            ARM_Date date2 = ARM_Date(JulianDate2);

            /// handle leap year
            int day = date2.Day;
            if( date2.Month == 2 && day > 28 && !date1.IsLeapYear() )
                day = 28;
            ARM_Date date2Sup(day, date2.Month, date1.Year);
            ARM_Date date2Inf;

           
            date2Inf = date2Sup;

            if ( date2Sup > date1 )
            {
               (void) date2Inf.SubstractYears(1); 
            }
            else
            {
               (void) date2Sup.AddYears(1);
            }
 
            yt = date2.Year - date2Sup.Year;

            yt += DaysBetweenDatesWithoutException(KACTUAL_ACTUAL, JulianDate1,
                                   date2Sup.GetJulian())
                 / DaysBetweenDatesWithoutException(KACTUAL_ACTUAL, date2Inf.GetJulian(),
                                    date2Sup.GetJulian()); 

            return(yt);
       };
       break;

       default :
            return(10000000);
    }
}


double CountYearsWithoutException(int DayCount, const ARM_Date& date1, const ARM_Date& date2)
{
    return CountYearsWithoutException( DayCount, date1.GetJulian(), date2.GetJulian() );
}


/*----------------------------------------------------------------------------*/
/*  Compute and returns in nextDate the smallest regular date STRICLY after   */
/*   settlement. Periodicity is given by frequency regular dates per year     */
/*        and refDate is a regular date.                                      */ 
/*----------------------------------------------------------------------------*/

void NextRegularDate(int freq, ARM_Date& refDate, ARM_Date& settlement, 
            ARM_Date& nextDate)
{
    nextDate = refDate;
    nextDate.Year = settlement.Year;

    nextDate.AddMonths(-12);
 
    int i = 1;
    ARM_Date tmpDate(nextDate);

    while ( nextDate <= settlement ) 
    {
        nextDate = tmpDate;
        nextDate.AddPeriodMult(freq,i);
        i++;
    //    nextDate.AddMonths(12/freq);
    }
}



/*----------------------------------------------------------------------------*/
/*       Compute and returns in prevDate the largest regular date before      */
/*    settlement. Periodicity is given by frequency regular dates per year    */
/*        and refDate is a regular date.                                      */
/*----------------------------------------------------------------------------*/

void PreviousRegularDate(int freq, const ARM_Date& refDate, const ARM_Date& settlement, 
                ARM_Date& prevDate)
{
    prevDate = refDate;
    prevDate.Year = settlement.Year;
  
    prevDate.AddMonths(12);

    int i = 1;
    ARM_Date tmpDate(prevDate);
 
    while ( prevDate > settlement ) 
    {
        prevDate = tmpDate;
        prevDate.AddPeriodMult(freq,-i);
        i++;
   //     prevDate.AddMonths(-(12/freq));
    }
}



double FromDayCountToBasis(int dayCount)
{
    double basis =365.0;
 
    switch (dayCount)
    {
        case KACTUAL_360 :
        case K30_360E :
        case K30_360 :
        {
            basis = 360.0;
            break;
        }
 
        case KACTUAL_ACTUAL :
        case KACTUAL_365 :
        {
            basis = 365.0;
            break;
        }
 
        default :
            basis = 365.0;
    }
    return(basis);
}




/*----------------------------------------------------------------------------*/
/*       Returns the year term between date1 and date2 using itsDayCount      */
/*        day count convention and regular flow dates with flow frequency     */
/*     int frequency and reference flow date refDate.                         */
/*----------------------------------------------------------------------------*/

double PeriodicYearTerm(int DayCount, int freq, ARM_Date& refDate, 
                        ARM_Date& date1, ARM_Date& date2)
{
    ARM_Date prevDate1, nextDate1, prevDate2, nextDate2;
    int nMonths;
    double  term, tt1, tt2;
 
 
    nMonths = 12 / freq;
 
    // compute term between date1 and first periodic date after date1
    // refDate == firstCouponDate
    // date1 == dateCalc
    // date2 == flowDate
 
    if ( refDate < date1 ) 
    {
       nextDate1 = refDate;

       while ( nextDate1 < date1 ) 
       {
           nextDate1.AddMonths(nMonths);
       }

       prevDate1 = nextDate1;
       prevDate1.AddMonths(-nMonths); 

       // here : prevDate1 <= date1 <= nextDate1
    }
    else 
    {
       prevDate1 = refDate;

       while ( prevDate1 > date1 ) 
       {
           prevDate1.AddMonths(-nMonths);
       }

       nextDate1 = prevDate1;
       nextDate1.AddMonths(nMonths);
    }

    tt1 = DaysBetweenDates(DayCount, date1, nextDate1);
    tt1 /= ((double) freq) * DaysBetweenDates(DayCount, prevDate1, nextDate1);
  
    //tt1 /= ((double) freq) * FromDayCountToBasis(DayCount);
    // Compute term between first periodic date before date2 and date2
    if ( refDate < date2 ) 
    {
       nextDate2 = refDate;

       while ( nextDate2 < date2 ) 
       {
           nextDate2.AddMonths(nMonths);
       }

       prevDate2 = nextDate2;
       prevDate2.AddMonths(-nMonths);
    }
    else 
    {
       prevDate2 = refDate;

       while ( prevDate2 > date2 ) 
       {
           prevDate2.AddMonths(-nMonths);
       }

       nextDate2 = prevDate2;
       nextDate2.AddMonths(nMonths);
    }


    tt2 = DaysBetweenDates(DayCount, prevDate2, date2);
    tt2 /= ((double) freq)*DaysBetweenDates(DayCount, prevDate2, nextDate2);

    // TMP tt2 /= ((double) freq)*FromDayCountToBasis(DayCount);
    term = tt1 + tt2 + (double) (prevDate2.Year - nextDate1.Year)  
                    + ((double) (prevDate2.Month - nextDate1.Month)) / 12.0;


    return(term);
}



ARM_Vector* ComputeCashFlowTerms(ARM_Date refDate, ARM_Vector* flowJulDates,
                            int dayCount, int frequency, ARM_Date couponDate)
{
    ARM_Vector *cfYearTerm;
    int size,
        i;
 
    size = flowJulDates->GetSize();
    cfYearTerm = new ARM_Vector(size);
 
    for (i = 0; i < size; i++)
    {
        (*cfYearTerm)[i] = PeriodicYearTerm(dayCount, frequency,
                           couponDate,
                           refDate,
                           (ARM_Date) (*flowJulDates)[i]);
    }
 
    return(cfYearTerm);
}
 
 



/******************************************************************
 *
 ******************************************************************/
ARM_Vector* CptStartDates(ARM_Date& inStartDate, 
                          ARM_Date& inEndDate, int frequency,
                          int fwdRule, int TypeStub, 
                          int intRule, char* ccy,
                          int adjFirstdate,
						  long rollDay)
{
    ARM_Date StartDate(inStartDate);
    ARM_Date EndDate(inEndDate);

    /// part to handle correctly the case of zero coupon
    ARM_Date AdjStartDate;

    int theRule;

    if ( intRule == 0 )
       theRule = 0;
    else
       theRule = fwdRule;

    if ((adjFirstdate) && ( theRule != 0 ))
    {
       AdjStartDate = StartDate;
       AdjStartDate.GoodBusinessDay(theRule, ccy);
    }
    else
       AdjStartDate = StartDate;

    if ( frequency == K_ZEROCOUPON )
       return new ARM_Vector(1, StartDate.GetJulian());

    /// other cases
    int fin, nbFlow = 0;
    int IsStub = 0;
    ARM_Date flowDate, AdjFlowDate, AdjEndDate;
    ARM_Date tmpDate;
    int i, SizeMax;

    AdjEndDate = EndDate;
    AdjEndDate.GoodBusinessDay(theRule, ccy);

    // Case of zero coupon : put in result the right start date
    SizeMax = int(((EndDate - StartDate)/365+1)*frequency);

    double *result = new double[SizeMax];
    memset(result, '\0', SizeMax * sizeof(double));

   // if the stub must be at the begining, we decrease <period> months
   // from EndDate until StartDate
    int step = 2;

    if ( TypeStub == K_SHORTSTART || TypeStub == K_LONGSTART )
    {
        flowDate = EndDate;

        flowDate.AddPeriod(-1*frequency, ccy, 0);

        AdjFlowDate = flowDate;

        AdjFlowDate.GoodBusinessDay(theRule, ccy);

        while ( AdjStartDate.Julian < AdjFlowDate.Julian )
        {
            result[nbFlow] = flowDate.Julian;

            nbFlow++;

            if (( frequency != K_DAILY )
                &&
                ( frequency != K_WEEKLY )
               )
            {
               flowDate = EndDate;
             
               flowDate.AddPeriodMult(-1*frequency, step, ccy);
            }
            else
            {
               flowDate.AddPeriod(-1*frequency, ccy, 0);
            }

            AdjFlowDate = flowDate;

            AdjFlowDate.GoodBusinessDay(theRule, ccy);

            step++;
        }

        if (AdjFlowDate.Julian == StartDate.Julian)
        {
            IsStub = 0;

            result[nbFlow] = flowDate.Julian;

            nbFlow++;
        }
        else 
        {
            IsStub = 1; 

            result[nbFlow] = StartDate.Julian;

            nbFlow++;
        }

    }
   // if the stub must be at the end, we increase <period> months
   // from StartDate until EndDate
    else
    {
        flowDate = StartDate;

        result[nbFlow] = flowDate.Julian;

        flowDate.AddPeriod(frequency, ccy, 0);

        AdjFlowDate = flowDate;

        AdjFlowDate.GoodBusinessDay(theRule, ccy);

        nbFlow++;

        while ( AdjFlowDate.Julian < AdjEndDate.Julian)
        {
            result[nbFlow] = flowDate.Julian;

            nbFlow++;

            flowDate = StartDate;
            
            flowDate.AddPeriodMult(frequency, step, ccy);

            AdjFlowDate = flowDate;

            AdjFlowDate.GoodBusinessDay(theRule, ccy);

            step++;
        }

        if (AdjFlowDate.Julian == AdjEndDate.Julian)
        {
            IsStub = 0;
        }
        else
        {
            IsStub = 1;
        }
    }

    fin = nbFlow;

	if ( (rollDay > 0) && (TypeStub == K_LONGSTART) )
	{
		IsStub	=	0;
	}

    // if TypeStub is LONG (start or end) we merge the last buckets
    if (IsStub && TypeStub == K_LONGSTART && nbFlow >= 2)
    {
        nbFlow--;

        result[nbFlow-1] = result[nbFlow];
    }

    if (IsStub && TypeStub == K_LONGEND && nbFlow >= 2)
    {
        nbFlow--;
    }

    //  Adjusting the dates
    double* resultAdj = new double[nbFlow];

    if ( TypeStub == K_SHORTSTART || TypeStub == K_LONGSTART )
    {
        for (i = 0; i < nbFlow; i++)
        {
            tmpDate = ARM_Date(result[i]);
            resultAdj[nbFlow-(i+1)] =
                tmpDate.GoodBusinessDay(theRule, ccy).Julian;
        }
    }
    else if (TypeStub == K_SHORTEND || TypeStub == K_LONGEND)
    {
        for (i = 0; i < nbFlow; i++)
        {
            tmpDate = ARM_Date(result[i]);
            resultAdj[i] = tmpDate.GoodBusinessDay(theRule, ccy).Julian;
        }
    }

	if ( rollDay > 0 )
	{
		long nbDays, month, year;
		long origin(0), coeff(1);

		if ( TypeStub == K_SHORTSTART || TypeStub == K_LONGSTART )
		{
			origin	=	nbFlow - 1;
			coeff	=	-1;
		}

		for ( i = 0; i < nbFlow; i++ )
		{
			tmpDate	=	result[i];

			year	=	tmpDate.GetYear();
			month	=	tmpDate.GetMonth();

			nbDays	=	tmpDate.GetDaysInMonth(month, year);

			if ( rollDay <= nbDays )
			{
			   tmpDate.ChangeDate(rollDay, month, year);
			}
			else
			{
			   tmpDate.ChangeDate(nbDays, month, year);
			}

			resultAdj[(origin+coeff*i)]	=	tmpDate.GoodBusinessDay(theRule, ccy).Julian;
		}
	}

    //  Delete doublons
    int size = 0;
    double prevDate = 0.0;

    double* lastResult = new double[nbFlow];

    for (i = 0; i < nbFlow; i++)
    {
        if ( prevDate != resultAdj[i] )
        {
            prevDate = resultAdj[i];
            lastResult[size] = prevDate;
            size++;
        }
    }

    // Create StartDates vector

    ARM_Vector* StartDates = new ARM_Vector(size, lastResult);

    delete [] result;
    delete [] resultAdj;
    delete [] lastResult;

    return(StartDates);
}



ARM_Vector* CptEndDates(ARM_Vector* startDates,ARM_Date& EndDate,
                        int fwdRule,int intRule,char* PayCalendar)
{     
    int size = startDates->GetSize();
    ARM_Vector* endDates = new ARM_Vector(size);
	int i = 0;
    for (; i < size-1; i++)
        endDates->Elt(i) = startDates->Elt(i+1);

    EndDate.GoodBusinessDay(fwdRule*intRule, PayCalendar);
    endDates->Elt(i) = EndDate.GetJulian();

    return endDates;
}



ARM_Vector* CptInterestTerms(ARM_Vector* startDates,ARM_Vector* endDates, int dayCount)
{
    int nFlows = startDates->GetSize();

    ARM_Vector* InterestTerms = new ARM_Vector(nFlows);

    for (int i = 0; i < nFlows-1; i++) 
    {
        // Perf
        switch (dayCount)
        {
            case KACTUAL_360 :
                InterestTerms->Elt(i) = (endDates->Elt(i)
                                          -startDates->Elt(i))/360.0;
            break;

            case KACTUAL_365 :
                InterestTerms->Elt(i) = (endDates->Elt(i)
                                          -startDates->Elt(i))/365.0;
            break;

            default : 
                InterestTerms->Elt(i) = CountYears(dayCount, 
                                             startDates->Elt(i), 
                                             endDates->Elt(i));
        }
    }
    
    switch (dayCount)
    {
        case KACTUAL_360 :
             InterestTerms->Elt(nFlows-1) = (endDates->Elt(nFlows-1)
                                     -startDates->Elt(nFlows-1))/360.0;
        break;

        case KACTUAL_365 :
             InterestTerms->Elt(nFlows-1) = (endDates->Elt(nFlows-1)
                                     -startDates->Elt(nFlows-1))/365.0;
        break;

        default : 
             InterestTerms->Elt(nFlows-1) = CountYears(dayCount,
                                           startDates->Elt(nFlows-1), 
                                           endDates->Elt(nFlows-1));
    }
    
    return InterestTerms;
}



void ARM_Date::JulianToStrDateDay(char* strDate) // Format : CCC DD/MM/YYYY
{
    sprintf(strDate, "%s %02d.%02d.%04d", ShortWeekDays[DayOfWeek],
                                          Day, Month, Year);
}

void ARM_Date::JulianToSpeDateDay(char* strDate) // Format : CCC DD/MM/YYYY
{
    sprintf(strDate, "%02d/%02d/%04d", Day, Month, Year);
}

string ARM_Date::toString() const
{
    char msg[255];

    sprintf(msg, "%s %02d.%02d.%04d", ShortWeekDays[DayOfWeek],
            Day, Month, Year);

    return string( msg );
}

string ARM_Date::toString(char aSeparator) const
{
	char vDate[11];

	sprintf(vDate, "%02d%c%02d%c%04d", Day, aSeparator, Month, aSeparator, Year);

	return	string(vDate);
}



/// function to get the corresponding month
/// very crude search! step by step
int GetNumMonthFromStr(char* m)
{
    int i;
    char month[4];


    sscanf(m, "%3s", month);

    (void) ARM_UPPERCASE(month);

    i = 0;

    while ( strcmp(month, ShortMonthsCapitalLetter[i]) && i < 11 )
    {
    ++i;
    }

    if ( strcmp( month, ShortMonthsCapitalLetter[i]) == 0 )
    {
       return(i+1);
    }
    else
    {
       return(-1);
    }
}



// function to get the number of days in the date's month
int GetMonthDays(ARM_Date& date)
{
    int month = date.GetMonth();
    int monthdays = MonthDays[month-1];
    
    if ( month == 2 && date.IsLeapYear())
    {
        monthdays = 29;
    }

    return monthdays;
}

/*
		pour pouvoir instancier ARM_GP_T_Vector<ARM_Date>
		FIXME : TODO : implement this method 
					
*/
std::ostream& operator<<(std::ostream&o,const ARM_Date&d)
{
	o<<d.toString('/'); 
	return o ;
}
/*----------------------------------------------------------------------*/
/*---- End of file ----*/
