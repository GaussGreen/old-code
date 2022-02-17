/*
 * $Log: dates.h,v $
 * Revision 1.37  2004/04/08 09:05:23  jpriaudel
 * modif pour  inline int isMinDate(void) const
 *
 * Revision 1.36  2003/12/12 12:53:18  jpriaudel
 * ajout de IsMinDate
 *
 * Revision 1.35  2003/11/25 10:11:39  rguillemot
 * Inflation Season Manager
 *
 * Revision 1.34  2003/10/16 16:10:55  ebenhamou
 * remove dotnet but on uppercase on static char*
 *
 * Revision 1.32  2003/10/07 15:14:00  ebenhamou
 * make toString virtual
 *
 * Revision 1.31  2003/09/26 07:30:02  ebenhamou
 * added conversion to string
 *
 * Revision 1.30  2003/09/01 12:55:27  jpriaudel
 * ajout de JulianToStrDateDay
 *
 * Revision 1.29  2003/08/15 06:59:37  ebenhamou
 * added CountYearsWithoutException for ARM_Date
 *
 * Revision 1.28  2003/07/22 13:58:29  mab
 * const supressed in operator -
 *
 * Revision 1.27  2003/07/22 13:35:17  mab
 * const suppressed in :
 * double operator - ( const ARM_Date& d)
 *
 * Revision 1.26  2003/07/21 10:24:22  jpriaudel
 * ajout d'un format dans la creation d'une date
 *
 * Revision 1.25  2003/07/17 06:52:47  ebenhamou
 * remove uncommented code
 *
 * Revision 1.24  2003/07/16 06:53:16  ebenhamou
 * version with AddPeriod with string
 *
 * Revision 1.23  2003/07/01 20:34:38  jpriaudel
 * include linalg
 *
 * Revision 1.22  2003/07/01 20:02:05  jpriaudel
 * suppression of include armglob
 *
 * Revision 1.21  2003/06/30 16:21:43  ebenhamou
 * for const correctness say that fction GetJulian is const
 *
 * Revision 1.20  2003/06/26 17:54:34  ebenhamou
 * create operator with bool
 *
 * Revision 1.19  2003/06/20 11:52:25  ebenhamou
 * proper rewritting of the code
 *
 * Revision 1.18  2003/04/01 17:17:33  jpriaudel
 * Ajout des fonctions DaysBetweenDatesWithoutException
 * et CountYearsWithoutException qui permettent d'avoir
 * des resultats negatifs
 *
 * Revision 1.17  2003/03/11 16:26:30  mab
 * Added : parameter adjFirstdate
 * in CptStartDates
 *
 * Revision 1.16  2003/02/10 17:31:01  mab
 * Added : US date Format : if (strcmp(ARM_DEFAULT_COUNTRY,"USD") == 0)
 *
 * Revision 1.15  2002/11/15 10:10:34  mab
 * Update of date format when USD
 *
 * Revision 1.13  2002/03/01 14:38:38  mab
 * Added : GetStrDate PrtStrDate
 * ARM_Date& AddPeriod(int, int cday, char* ccy = ARM_DEFAULT_COUNTRY);
 * void JulianToStrDateN(char* strDate)
 *
 * Revision 1.12  2002/01/29 15:46:01  mab
 * Annulation de la correction de :  ARM_Date::PiborDelivery(void)
 * La delivery reste : Le Mercredi
 *
 * Revision 1.11  2001/12/28 14:32:19  arm
 * Rajout de AddPeriodMult
 *
 * Revision 1.10  2001/03/12 09:57:30  abizid
 * Optimisation Code
 *
 * Revision 1.9  2001/01/30 09:54:01  smysona
 * Today supprime en inline, passe en .cpp
 * Pb compilo SUN
 *
 * Revision 1.8  2000/11/14 08:35:20  mab
 * Fonction AddPeriod predn en compte la Currency
 *
 * Revision 1.7  2000/10/25 10:01:03  smysona
 * Les operateurs de test sont maintenant declares avec des parametres const.
 *
 * Revision 1.6  2000/04/25 15:09:35  mab
 * extern "C" enleve de la section includes
 *
 * Revision 1.5  2000/02/02 08:43:07  mab
 * Rajout : AdjustToBusDate()
 *
 * Revision 1.4  1999/03/02 10:36:05  ypilchen
 * Rajout test Ds inline ARM_Date::ARM_Date(char* stringDate, char* format)
 *
 * Revision 1.3  1999/03/02 10:19:08  ypilchen
 * Rajout de test Ds : inline ARM_Date::ARM_Date(char* stringDate, char* format)
 *
 * Revision 1.2  1999/02/18 18:46:47  nicolasm
 * Ajout du $Log
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : dates.h                                                      */
/*                                                                            */
/* DESCRIPTION : Object dealing with dates                                    */
/*                                                                            */
/* DATE        : Wed Apr 17 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _DATES_H
#define _DATES_H

#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

#include "linalg.h"
#include "expt.h"

#include "period.h"

//#include <iostream>

class ARM_Calendar;



/*----------------------------------------------------------------------------*/
#define IsLeap(y) ((y % 4 == 0) && (y % 4000 != 0) && ((y % 100 != 0) \
                           || (y % 400 == 0)))
 
#define MODULO(x,y)  (( x >= 0 ) ? ( x % y ) : ((x+y*(1+(int)((-x)/y))) % y))
 

#define MINMONTH      1
#define MAXMONTH      12
#define MINDAY        1
#define MAXDAY        31
#define MINWEEKDAY    0
#define MAXWEEKDAY    6
#define WEEKDAYS      7
#define YEARMONTHS    12
#define BADDATE       -1L
#define MINYEAR       -4700L
#define MAXYEAR       25000L
#define MINDATE       4749L      // 1 janvier 4700 av JC
#define MAXDATE       10852487L  // 31 decembre 25000


extern int MonthDays[MAXMONTH];

extern char ShortMonths[MAXMONTH][4];
extern char ShortMonthsCapitalLetter[MAXMONTH][4];


/*----------------------------------------------------------------------------*/

extern double FromDayCountToBasis(int dayCount);


class ARM_Date : public ARM_Object
{
    private:

        double Julian;

        int    Day;        // 1 a 31
        int    Month;      // 1 a 12
        int    Year;       // < 0 si avant JC

        int    DayOfWeek;  // 0.. 6


        static ARM_Calendar* vacationsCalendar;

    public:
	
		
		ARM_Date(void);

        // The Next Constructor add a period specified by p parameter
        // formatted as <Number>< D | W | M | Y>
        // Exp : 1D, 20D, 3W, 2M, 3Y 

        // ARM_Date(char* p);


        ARM_Date(int d, int m, int y);         
        ARM_Date(int d, int m, int y, int h, int mn, int sec);         


        // Format : DD/MM/YYYY is also possible
        // the separator '/' may be replaced by '-' or '.'
        // Other Formats should be provided in further versions

        ARM_Date(char* stringDate, char* format = "DD.MM.YYYY");
		ARM_Date(const std::string& date,const std::string& format); 
        ARM_Date(double julianDate);

        ARM_Date(const ARM_Date& d);

		virtual ARM_Object*	Clone( ) { return new ARM_Date(*this);	}

		virtual ARM_Object*	Clone( ) const { return (const_cast<ARM_Date*>(this))->Clone();	}

        ~ARM_Date(void) {}

        static void InitCalendar(char* calFile,bool isString = false);

        static void SetCalendar(ARM_Calendar* cal);

        static ARM_Calendar* GetCalendar(void);

//        operator double() { return(Julian); }

        inline int GetDay(void) const { return Day; }
        inline int GetMonth(void) const { return Month; }
        inline int GetYear(void) const { return Year;}

        inline void SetDay(int d) { Day = d; }
        inline void SetMonth(int m) { Month = m; }
        inline void SetYear(int y) { Year = y;}

        inline int GetDayOfWeek(void) const { return(DayOfWeek); }
 
		/// need to be const to use it properly
		/// when using const
        inline double GetJulian(void) const { return(Julian); }
 
		inline void SetJulian(double d) { 
			Julian = d; 
			ARM_Date tmp = ARM_Date(d);
			Day		= tmp.GetDay();        // 1 a 31
			Month	= tmp.GetMonth();      // 1 a 12
			Year	= tmp.GetYear();   
		}

        inline double GetFracDay(void) const { return(Julian-floor(Julian)); }

        char* GetStrDate(void) const 
        {
            char* stringDate = new char[50];
            this->JulianToStrDateN(stringDate);
            return stringDate;
        }

        char* PrtStrDate(void) const 
        {
            char* stringDate = new char[50];
            this->JulianToStrDate(stringDate);
            return stringDate;
        }

        virtual string toString() const;

		virtual string toString(char aSeparator) const;


        // Operators

        ARM_Date& operator = (const ARM_Date& d);

        double operator - (const ARM_Date& d) const ;

        ARM_Date& operator - (long x);

        friend ARM_Date& operator + (long x, ARM_Date&);
        friend ARM_Date& operator + (ARM_Date&, long x);

        ARM_Date& operator ++ (void);
        ARM_Date& operator -- (void);

        ARM_Date& operator += (long x);
        ARM_Date& operator -= (long x);

        bool operator != (const ARM_Date& d) const { return(Julian != d.Julian); } 
        bool operator == (const ARM_Date& d) const { return !operator!=(d); }

        bool operator >  (const ARM_Date& d) const { return(Julian >  d.Julian); }
        bool operator <= (const ARM_Date& d) const { return !operator>(d); }
        bool operator <  (const ARM_Date& d) const { return(Julian <  d.Julian); }
        bool operator >= (const ARM_Date& d) const { return !operator<(d); }
        
        ARM_Date& SubstractYears(int n);

        ARM_Date& AddHours(int n);
        ARM_Date& AddDays(int n);
        ARM_Date& AddMonths(int months, int GOTO_END_OF_MONTH = 0);
        ARM_Date& AddYears(int);
        ARM_Date& AddPeriod(int, char* ccy = ARM_DEFAULT_COUNTRY, int AdjustDaily = 1);
        ARM_Date& AddPeriod(int, int cday, char* ccy = ARM_DEFAULT_COUNTRY);
        /// version with string
        ARM_Date& AddPeriod(const string& s, const char* ccy = ARM_DEFAULT_COUNTRY);

        ARM_Date& AddPeriodMult(int, int mult = 1, char* ccy = ARM_DEFAULT_COUNTRY, int GOTO_END_OF_MONTH = 0);
		ARM_Date& AddPeriodMult(int, int mult , const std::string& , int GOTO_END_OF_MONTH = 0);

		/*****************************************
			FP : addperiod new methods : 
			
			AddPeriodNoAdj(int unit, int length) 
			add a positive or negative period without adjusting to the next ou previous business day

			AddPeriodAdj(int unit, int length, char* ccy, int rollingConvention, int GOTO_END_OF_MONTH = 0)
			add a positive or negative period and adjust to the next ou previous business day
		
			unit : "K_ANNUAL, K_SEMIANNUAL, K_QUARTERLY, K_BIMONTHLY, K_MONTHLY, K_WEEKLY, K_DAILY"
			rollingConvention : "K_PREVIOUS, K_MOD_PREVIOUS, K_FOLLOWING, K_MOD_FOLLOWING"

			THOSE METHODS DO NOT MODIFY THE CURRENT OBJECT (*this)

		*****************************************/
		ARM_Date& AddPeriodNoAdj(int unit, int length) ;
		ARM_Date& AddPeriodAdj(int unit, int length, char* ccy, int rollingConvention);
		
		ARM_Date& AddPeriodNoAdj(const Period& period) ;
		ARM_Date& AddPeriodAdj(const Period& period, char* ccy, int rollingConvention);
		
		/* end of FP */
        
		
		ARM_Date& NextNonWeekEndDay(int d);

        ARM_Date& AdjustToBusDate(char* ccy, long rule);

        ARM_Date& NextBusinessDay(char* ccy = ARM_DEFAULT_COUNTRY);// Not Holidays, 
                                                                   // Not week end 

        ARM_Date& NextBusinessDay(int nb, 
                                  char* ccy = ARM_DEFAULT_COUNTRY);// Not Holidays, 
                                                                   // Not week end

        ARM_Date& NextBusinessDay(int nb, const std::string&ccy) ;

        ARM_Date& UnadjustDate(int nb);

        // Not Holidays, Not week end
        ARM_Date& PreviousBusinessDay(char* ccy = ARM_DEFAULT_COUNTRY);
 
        // Not Holidays, Not week end
        ARM_Date& PreviousBusinessDay(int nb, 
                                  char* ccy = ARM_DEFAULT_COUNTRY);

        ARM_Date& GoodBusinessDay(int rule=0, char *ccy = ARM_DEFAULT_COUNTRY);

        ARM_Date& GapBusinessDay(int gap, char* ccy = ARM_DEFAULT_COUNTRY, int fwdRule = K_FOLLOWING);

        ARM_Date& NonWeekEndDay(int d);
        ARM_Date& GapNonWeekEndDay(int d);

        ARM_Date& GetLivretAFixingDate();

        void GetTime(int* h, int* mn, int* s);

        double DMYToJulian(void);
        void JulianDateToDMY(void); 

        int IsOrdinaryYear(void);
        int IsWeekEndDay(void);


        int isValid(void) { return( Julian != BADDATE ); }

        inline int isMinDate(void) const { return( Julian == MINDATE ); }

        void ChangeDate(int d, int m, int y);        
		void ChangeDate(double julianDate);

        void SysToday(void);

        void Today(void);

        void JulianToStrDate(char* strDate) const // Format : DD/MM/YYYY
        {
            sprintf(strDate, "%02d.%02d.%04d", Day, Month, Year);
        } 
        
        void JulianToStrDateN(char* strDate) const // Format : YYYY/MM/DD
        {
            sprintf(strDate, "%4d%02d%02d", Year, Month, Day);
        } 

        void JulianToStrDateDay(char* strDate); // Format : CCC DD/MM/YYYY
 
		void JulianToSpeDateDay(char* strDate); // Format : CCC DD/MM/YYYY

        void GetCompleteStrDate(char* d)
        {
            int h, mn, sec;
                     

            GetTime(&h, &mn, &sec);

            sprintf(d, "%2d.%2d.%4d %2d:%2d:%2d", GetDay(), GetMonth(), 
                                                  GetYear(), h, mn, sec);
        }

        int IsLeapYear(void)
        {
            return(IsLeap(this->Year));
        }
 
        int   GetDaysInMonth(int m, int y);
        static char* GetLongDayName(int d);
        static char* GetLongMonthName(int m);
        static char* GetShortDayName(int d);
        static char* GetShortMonthName(int m);

        // Friend Methods

        friend void GetToday(char* Today);

        friend double DaysBetweenDates(int DayCount, ARM_Date& d1, 
                                       ARM_Date& d2);
        friend double DaysBetweenDates(int DayCount, double Juliand1, 
                                       double Juliand2);
        friend double DaysBetweenDatesWithoutException(int DayCount, double Juliand1, 
                                       double Juliand2);

		friend int CountBusinessDays(const ARM_Date& fromDate, const ARM_Date& toDate, char* calendar);

        friend double CountYears(int DayCount, ARM_Date& d1, ARM_Date& d2); 
        friend double CountYears(int DayCount, double JulianDate1,
                                 double JulianDate2);
        friend double CountYearsWithoutException(int DayCount, double JulianDate1,
								double JulianDate2);
		friend double CountYearsWithoutException(int DayCount, const ARM_Date& date1, 
								const ARM_Date& date2);
        
        friend void NextRegularDate(int freq, ARM_Date& refDate, 
                                    ARM_Date& settlement, 
                                    ARM_Date& nextDate);

        friend void PreviousRegularDate(int freq, const ARM_Date& refDate, 
                                        const ARM_Date& settlement, 
                                        ARM_Date& prevDate);

        friend double PeriodicYearTerm(int DayCount, int freq, 
                                       ARM_Date& refDate, 
                                       ARM_Date& date1, ARM_Date& date2);

        friend ARM_Vector* ComputeCashFlowTerms(ARM_Date refDate, 
                            ARM_Vector* flowJulDates,
                            int dayCount, int frequency,
                            ARM_Date couponDate);
 
        friend ARM_Vector* CptStartDates(ARM_Date& StartDate,
                            ARM_Date& EndDate, int frequency, int FwdRule,
                            int TypeStub = K_LONGSTART, 
                            int intRule = K_ADJUSTED,
                            char* ccy = ARM_DEFAULT_COUNTRY,
                            int adjFirstdate = 1,
							long rollDay = -1);

         friend ARM_Vector* CptEndDates(ARM_Vector* startDates,
                                        ARM_Date& EndDate,
                                        int fwdRule,
                                        int intRule = K_ADJUSTED,
                                        char* PayCalendar = ARM_DEFAULT_COUNTRY);

        friend ARM_Vector* CptInterestTerms(ARM_Vector* startDates,
                                            ARM_Vector* endDates, 
                                            int dayCount=KACTUAL_365);

        ARM_Date& Paque(int Annee);
        int IsBusinessDay(char* Country);

        ARM_Date& PiborDelivery(const char* ccy = NULL);

        void Calc_MATIF_BDFutureDelivery(int m, int y)
        {
            this->ChangeDate(1, m, y);

            *this = this->PiborDelivery();
        }   

        void Calc_MATIF_IRFutureDelivery(int m, int y)
        {
            this->ChangeDate(1, m, y);

            *this = this->PiborDelivery();
        } 

        void Calc_BUND_Delivery(int m, int y)
        {
            this->ChangeDate(10, m, y);
 
            while (!(this->IsBusinessDay("DEM")))
            {
                this->AddDays(1);
            }
        }

        void Calc_GILT_Delivery(int m, int y)
        {
            this->ChangeDate(1, m, y);

            this->ChangeDate(GetDaysInMonth(m, y), m, y);

 
            while (!(this->IsBusinessDay("GBP")))
            {
                this->AddDays(-1);
            }
        }
};

/// in order to print in the view method Day Count long arg
const char* DCGetName( long DCF );

int GetNumMonthFromStr(char* m);

inline void GetMonthYearFromExpiryDate(char* str, int* m, int* y)
{
    char month[20];
    char year[10];



    sscanf(str, "%[^0123456789]%[0123456789]", month, year);

    *m = GetNumMonthFromStr(month);
            
    if ( (*m) < 0 )
    {
	   char msg [200];
	   sprintf(msg, "Invalid Month in %s: Must be like : DEC, JAN, FEB, ...", str);
       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID, msg);
       throw x;
    }

    *y = atoi(year);
 
    if ( strlen(year) == 2 )
    {
       if ( *y < 89 )
       {
          *y = 2000+(*y);
       }
       else
       {
          *y = 1900+(*y); 
       }
    }
}



// IR Future : Third Wednesday of 3,6,9,12th Month 

inline ARM_Date& ARM_Date::PiborDelivery(const char* ccy)
{
    ARM_Date NextEch;
    int m;
 
 
 
    NextEch.ChangeDate(Day, Month, Year);
 
    m = Month;
    m = m % 3;
 
    if ( m == 0 )
    {
       NextEch.AddDays(-NextEch.GetDay()+1);
 
       if ( NextEch.GetDayOfWeek() <= 3 )
       {
          NextEch.AddDays(17-NextEch.GetDayOfWeek());
       }
       else
       {
          NextEch.AddDays(24-NextEch.GetDayOfWeek());
       }
 
	   if (ccy)
	   {
			if (strcmp(ccy,"CHF") == 0)
				NextEch.AddDays(2);
	   }

       // Critere pour determiner qd on passe au contrat suivant
 
       if ( Day < NextEch.GetDay() )
       {
          (*this) = NextEch;
 
          return(*this);
       }
    }
 
    m = 3-m;

    NextEch.AddMonths(m);

    NextEch.AddDays(-NextEch.GetDay()+1);
 
    if ( NextEch.GetDayOfWeek() <= 3 )
    {
       NextEch.AddDays(17 - NextEch.GetDayOfWeek());
    }
    else
    {
       NextEch.AddDays(24 - NextEch.GetDayOfWeek() );
    }
 
    // Correction

	if (ccy)
	{
		if (strcmp(ccy,"CHF") == 0)
			NextEch.AddDays(2);
	}

    (*this) = NextEch;

    return (*this);
}



inline ARM_Date& ARM_Date::NextBusinessDay(char* ccy)
{
    *this = this->AddDays(1);

    while (!(this->IsBusinessDay(ccy)))
    {
        *this = this->AddDays(1);
    }
 
    return(*this);
}



inline ARM_Date& ARM_Date::AdjustToBusDate(char* ccy, long rule)
{
    ARM_Date busDate;
	
    busDate = *this;
	
    if (this->IsBusinessDay(ccy))
		return *this;
	
	switch( rule )
	{
	/// no adjustment
	case 0:
		break;
	/// Following
	case 1: 
		busDate = this->NextBusinessDay(ccy);
		break;
	/// previous
	case -1 :
		// ----- Get the previous Day
		busDate = this->PreviousBusinessDay(ccy);
		break;
	/// modified following
	case 2:
       busDate = busDate.NextBusinessDay(ccy);
       if ( busDate.GetMonth() != this->GetMonth() )
             busDate = busDate.PreviousBusinessDay(ccy);
	   break;
	/// modified previous
	case -2:
		// ----- Get the previous Day
		busDate = busDate.PreviousBusinessDay(ccy);
        if ( busDate.GetMonth() != this->GetMonth() )
             busDate = busDate.NextBusinessDay(ccy);
 	   break;
	default:
          throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Unknown rule");
	} 
     
	*this = busDate;
    return *this;
}


inline ARM_Date& ARM_Date::UnadjustDate(int nb)
{
    long nbDays = nb - GetDay();
	
	if (nbDays > 15)
	{
		int m = GetMonth();
		
		while (GetMonth() == m)
			this->AddDays(-1);
		
		while (GetDay() > nb) 
			this->AddDays(-1);
	}
	else if (nbDays < -15)
	{
		int m = GetMonth();
		
		while (GetMonth() == m)
			this->AddDays(1);
		
		while (GetDay() < nb) 
			this->AddDays(1);
	}
	else if (nbDays > 0)
	{
		int m = GetMonth();
		
		while ((GetMonth() == m) && (GetDay() < nb))
			this->AddDays(1);

		if (GetMonth() != m)
			this->AddDays(-1);
	}
	else if (nbDays < 0)
	{
		int m = GetMonth();
		
		while ((GetMonth() == m) && (GetDay() > nb))
			this->AddDays(-1);

		if (GetMonth() != m)
			this->AddDays(1);
	}

    return(*this);
}


inline ARM_Date& ARM_Date::NextBusinessDay(int nb, char* ccy)
{
    for (int i = 0; i < nb; i++)
    {
        this->NextBusinessDay(ccy);
    }

    return(*this);
}


inline ARM_Date& ARM_Date::NextBusinessDay(int nb, const std::string& ccy)
{
	return NextBusinessDay(nb,(char*)ccy.c_str()); 
}


inline ARM_Date& ARM_Date::PreviousBusinessDay(char* ccy)
{
    *this = this->AddDays(-1);
 
    while (!(this->IsBusinessDay(ccy)))
    {
        *this = this->AddDays(-1);
    }
 
    return(*this);
}
 
 
 
inline ARM_Date& ARM_Date::PreviousBusinessDay(int nb, char* ccy)
{
    for (int i = 0; i < nb; i++)
    {
        this->PreviousBusinessDay(ccy);
    }
 
    return(*this);
}



/*----------------------------------------------------------------------------*/
/* Conversion Julian <-> DMY                                                  */
/* Algorithm swipped directly from "Numerical Recipes in C"                   */
/* by Press, Flannery, Teukolsky and vettering, p10.                          */
/*----------------------------------------------------------------------------*/
 
 
 
 
#define GREGOR  2299161L
 
inline void JulianToDMY(double Jul, int* d, int* m, int* y)
{
    long Ja, JAlpha, Jb, Jc, Jd, Je;
 
    if (( Jul != BADDATE ) && ( Jul >= MINDATE ) && ( Jul <= MAXDATE ))
    {
       if ( Jul >= GREGOR )
       {
          JAlpha = (long) (((double)(Jul - 1867216L) - 0.25) / 36524.25);
          Ja = (long) (Jul + 1 + JAlpha - (long)(0.25*JAlpha));
       }
       else
       {
          Ja = (long) Jul;
       }
 
       Jb = Ja + 1524;
       Jc = (long) (6680.0 + ((double)(Jb - 2439870L) - 122.1) / 365.25);
       Jd = (long) (365*Jc + (0.25*Jc));
       Je = (long)((Jb - Jd)/30.6001);
       *d = (int)(Jb - Jd - (int)(30.6001*Je));
       *m = (int)Je - 1;
 
       if (*m > 12)
          *m -= 12;
 
       *y = (int)(Jc - 4715);
 
       if (*m > 2)
          --(*y);
 
       if ( *y <= 0)
          --(*y);
    }
    else
    {
	   char msg [500];
	   sprintf(msg, "Invalid Date : Converting Julian Date (%.0f) to dd/mm/yyyy failed", Jul);
       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID, msg);
 
       throw x;
    }
}



inline int DaysInMonth(int m, int y)
{
    if (( m == 2 ) && IsLeap(y))
    {
      return(29);
    }
    else
    {
       return(MonthDays[m - 1]);
    }
}



inline int CheckForValidDate(int d, int m, int y)
{
    if ( y / 1960.0 < 1.0 )
    {
       return(0);
    }

    if (( y < MINYEAR ) || ( y > MAXYEAR ) || ( y == 0 )
        || 
        ( m < MINMONTH ) || ( m > MAXMONTH ))
    {
       return(0);
    }
 
    return (d <= DaysInMonth(m, y));
}



inline int JulianToDayOfWeek(double Jul)
{
    int dow;
    int vTmp;
 
 
    vTmp = int(floor(Jul + 1));
 
    dow = vTmp % 7L;
 
    return(dow);
}



/*----------------------------------------------------------------------------*/
/* Conversion Julian <-> DMY                                                  */
/* Algorithm swipped directly from "Numerical Recipes in C"                   */
/* by Press, Flannery, Teukolsky and vettering, p10.                          */
/*----------------------------------------------------------------------------*/
 
#define IGREG (15+31L*(10+12L*1582))
 
inline double DMYtoJulian(int d, int m, int y)
{
    long Ja, Jm, Jy;
    double Jul;
    
    char buf[100];
 
 
    if (!CheckForValidDate(d, m, y))
    {
       sprintf(buf, "Invalid Date : %d-%d-%d", d, m, y);

       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID, buf);
 
 
       throw x;
 
       return(BADDATE);
    }
 
    if ( y < 0 )
    {
       y++;
    }
 
    if ( m > 2 )
    {
       Jy = y;
       Jm = m + 1;
    }
    else
    {
       Jy = y - 1;
       Jm = m + 13;
    }
 
    Jul = (floor(365.25*(double)Jy)+floor(30.6001*(double)Jm)+d+1720995L);
 
    if ( d+31L*(m+12L*y) >= IGREG )
    {
       Ja = (long) (0.01*Jy);
       Jul += (2-Ja+(long)(0.25*Ja));
    }

    return(Jul);
}



/*----------------------------------------------------------------------------*/
/*  Change date                                                               */
/*----------------------------------------------------------------------------*/
 
inline void ARM_Date::ChangeDate(int d, int m, int y)
{
    double julianDate;
 
 
    julianDate = DMYtoJulian(d, m, y);
 
    if ( julianDate != BADDATE )
    {
       Julian = julianDate;
       Day = d;
       Month = m;
       Year = y;
 
       DayOfWeek = JulianToDayOfWeek(Julian);
    }
    else
    {
       Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
                       "Invalid Date");
       throw x;
    }
}



inline void ARM_Date::ChangeDate(double julianDate)
{
    Julian = julianDate;
 
    JulianDateToDMY();
}



inline ARM_Date& ARM_Date::operator = (const ARM_Date& d)
{
    Day       = d.Day;
    Month     = d.Month;
    Year      = d.Year;
    Julian    = d.Julian;
    DayOfWeek = d.DayOfWeek;
  
 
    return(*this);
}



inline void ARM_Date::JulianDateToDMY(void)
{
    int d;
    int m;
    int y;
 
 
    JulianToDMY(Julian, &d, &m, &y);
 
    Day = d;
    Month = m;
    Year = y;
 
    DayOfWeek = JulianToDayOfWeek(Julian);
}



inline double ARM_Date::DMYToJulian(void)
{
    double julianDate;
 
 
    julianDate = DMYtoJulian(Day, Month, Year);
 
    return(julianDate);
}
 



inline ARM_Date::ARM_Date(void)
{
    SetName(ARM_DATE);
 
    Today();
}



inline ARM_Date::ARM_Date(double julianDate)
{
    SetName(ARM_DATE);
 
    Julian = julianDate;
 
    JulianDateToDMY();
}
 
 
 
inline ARM_Date::ARM_Date(int d, int m, int y)
{
    SetName(ARM_DATE);
 
    ChangeDate(d, m, y);
}
 
 
 
inline ARM_Date::ARM_Date(int d, int m, int y,
                          int h, int mn, int sec)
{
 
    SetName(ARM_DATE);
 
    ChangeDate(d, m, y);
 
    // Add hours, minutes and seconds
 
    double dayFrac;
 
 
    dayFrac = ((double) h ) + ((double) mn) / 60.0 + ((double) sec) / 3600.0;
    dayFrac /= 24.0;
 
    Julian += dayFrac;
}


inline  ARM_Date::ARM_Date(const std::string& date,const std::string& format) 
{
	new(this)ARM_Date((char*)date.c_str(),(char*)format.c_str()); 
}
inline ARM_Date::ARM_Date(char* stringDate, char* format)
{
    int d;
    int m;
    int h = 0, mn = 0, s = 0;
    int y;
    char c1, c2;
 
 
 
    SetName(ARM_DATE);
 
    if (( strcmp(format, "DD.MM.YYYY hh:mm:ss") == 0 )
        ||
        ( strcmp(format, "DD/MM/YYYY hh:mm:ss") == 0 )
        ||
        ( strcmp(format, "DD-MM-YYYY hh:mm:ss") == 0 )
       )
    {
       sscanf(stringDate, "%2d/%2d/%4d %2d:%2d:%2d", &d, &m, &y, &h, &mn, &s);
 
       ARM_Date tmpDate(d, m, y, h, mn, s);
 
       *this = tmpDate;
 
    }
	//	Calypso format (for the time beeing...) 
	else if ( strcmp(format, "YYYYMMDD hh:mm") == 0 ) 
	{
		sscanf(stringDate, "%4d/%2d/%2d %2d:%2d", &d, &m, &y, &h, &mn);
		s=0; 
		new(this)ARM_Date(d, m, y, h, mn, s);
	}
    else // Other possible format
    {
       if ( ( strcmp(format, "DD/MM/YYYY") == 0 )
            || ( strcmp(format, "DD-MM-YYYY") == 0 )
            || ( strcmp(format, "DD.MM.YYYY") == 0 )
          )
       {
          if ( strlen(stringDate) < 10 )
          {
             throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date");
          }
 
          sscanf(stringDate, "%d%c%d%c%4d", &d, &c1, &m, &c2, &y);
 
          ARM_Date tmpDate(d, m, y, h, mn, s);
 
          *this = tmpDate;
       }
       else if ( ( strcmp(format, "YYYY.MM.DD") == 0 )
                 ||
                 ( strcmp(format, "YYYY-MM-DD") == 0 )
                 || 
                 ( strcmp(format, "YYYY/MM/DD") == 0 )
               )
       {
          if ( strlen(stringDate) < 10 )
          {
             throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date");
          }

          sscanf(stringDate, "%4d%c%d%c%d", &y, &c1, &m, &c2, &d);
 
          ARM_Date tmpDate(d, m, y, h, mn, s);
 
          *this = tmpDate;
       }
       else if ( ( strcmp(format, "MM.DD.YYYY") == 0 )
                 ||
                 ( strcmp(format, "MM-DD-YYYY") == 0 )
                 || 
                 ( strcmp(format, "MM/DD/YYYY") == 0 )
               )
       {
          if ( strlen(stringDate) < 10 )
          {
             throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date");
          }

          sscanf(stringDate, "%d%c%d%c%4d", &m, &c1, &d, &c2, &y);
 
          ARM_Date tmpDate(d, m, y, h, mn, s);
 
          *this = tmpDate;
       }
       else if ( strcmp(format, "YYYYMMDD") == 0 )
       {
          if ( strlen(stringDate) < 8 )
          {
             throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date");
          }

          sscanf(stringDate, "%4d%2d%2d", &y, &m, &d);
 
          ARM_Date tmpDate(d, m, y, h, mn, s);
 
          *this = tmpDate;
       }
       else
       {
          if ( strlen(stringDate) < 10 )
          {
             throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date");
          }

          throw Exception(__LINE__, __FILE__,
                     ERR_DATE_NOT_VALID, "Invalid Date FORMAT");
       }
    }
}


inline ARM_Date::ARM_Date(const ARM_Date& d):ARM_Object(d) 
{
	SetName(ARM_DATE);
 
    Julian = d.Julian;
 
    Day    = d.Day;
 
    Month  = d.Month;
 
    Year   = d.Year;
 
    DayOfWeek = d.DayOfWeek;
}



inline ARM_Date& ARM_Date::AddMonths(int nMonths, int GOTO_END_OF_MONTH)
{
    int k, m1, m2;
    int prevMonth;
    int prevYear;
 
 
 
    prevMonth = Month;
    prevYear  = Year;
 
    m1 = Month + nMonths-1;
    m2 = MODULO(m1, 12);
    k = (m1 - m2) / 12;
 
    Year += k;
 
    Month = m2+1;
 
 
    if ( Month == 2 ) // February
    {
       if ( Day >= 29 )
       {
          if (IsLeap(Year))
          {
             Day = 29;
          }
          else
          {
             Day = 28;
          }
       }
    }
    else
    {
       if ( MonthDays[Month-1] < Day )
       {
          Day = MonthDays[Month-1];
       }
    }

    /* If the option GOTO_END_OF_MONTH is set, which means */
    /* if we are at the end of the origin month, we must   */
    /* stay at the end of the new calculated month         */
    /* 28 fev 1995 + 1 Month = 31 Mars 1995                */
 
    if (GOTO_END_OF_MONTH)
    {
       if ( Day == DaysInMonth(prevMonth, prevYear) )  // At the End of a month
       {
          Day = DaysInMonth(Month, Year);
       }
    }
 
    Julian = DMYToJulian();
 
    DayOfWeek = JulianToDayOfWeek(Julian);
 
    return(*this);
}

/*
		pour pouvoir instancier ARM_GP_T_Vector<ARM_Date>
		FIXME : TODO : implement this method 
					
*/
std::ostream& operator<<(std::ostream&o,const ARM_Date&d);



#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
