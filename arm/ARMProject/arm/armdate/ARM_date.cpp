#include "ARM_Date.h"

int LocalMonthDays[MAXMONTH] = { 31, 28, 31, 30, 31, 30, 31, 31, 
                            30, 31, 30, 31 };

ARM_Local_Date& ARM_Local_Date::AddPeriod(int period)
{
    int freq = abs(period);
    
    int signedFreq = 12 / period; 
   
    switch (freq) 
    {
        case K_ANNUAL:
        case K_SEMIANNUAL:
        case K_QUARTERLY:
        case K_BIMONTHLY:
        case K_MONTHLY:
        {
            this->AddMonths (signedFreq);
        };
        break;

        case K_WEEKLY :
		{
			int pSign = ((period >= 0) ? 1 : -1);
			this->AddDays (7*pSign);
		};
        break;

        case K_DAILY :
		{
			int pSign = ((period >= 0) ? 1 : -1);
			this->AddDays (pSign);
		};
		break;

        default :
			this->AddMonths (3);
			break;
    }
	
    return (*this);
}


ARM_Local_Date& ARM_Local_Date::AddPeriodMult(int period, int mult)
{
    int freq = abs(period);
    
    int signedFreq = 12/period; 
   
    switch (freq) 
    {
        case K_ANNUAL :
        case K_SEMIANNUAL :
        case K_QUARTERLY :
        case K_BIMONTHLY :
        case K_MONTHLY :
        {
            this->AddMonths(signedFreq*mult);
        };
        break;

        case K_WEEKLY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );

            this->AddDays(7*pSign*mult);
        };
        break;

        case K_DAILY :
        {
            int pSign = ( (period >= 0) ? 1 : -1 );
            this->AddDays(pSign);
			int i = 0;
			while (i < mult)
            {
				this->AddDays(pSign);
				i++;
			}
        };
        break;

        default :
           this->AddMonths(3);
        break;
    }
    return(*this);
}


ARM_Local_Date& ARM_Local_Date::AddDays (int d)
{
    Julian += (double)d; 

    JulianDateToDMY ();

    DayOfWeek = JulianToDayOfWeek (Julian);

    return (*this);
}

ARM_Local_Date& ARM_Local_Date::AddMonths(int nMonths, int GOTO_END_OF_MONTH)
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
  
    if(Month == 2)
    {
       if (Day >= 29)
       {
          if(IsLeap (Year))
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
       if(LocalMonthDays[Month-1] < Day)
       {
          Day = LocalMonthDays[Month-1];
       }
    }

    /* If the option GOTO_END_OF_MONTH is set, which means */
    /* if we are at the end of the origin month, we must   */
    /* stay at the end of the new calculated month         */
    /* 28 fev 1995 + 1 Month = 31 Mars 1995                */
 
    if(GOTO_END_OF_MONTH)
    {
       if(Day == DaysInMonth (prevMonth, prevYear))  // At the End of a month
       {
          Day = DaysInMonth (Month, Year);
       }
    }
 
    Julian = DMYToJulian ();
 
    DayOfWeek = JulianToDayOfWeek (Julian);
 
    return (*this);
}

int ARM_Local_Date::JulianToDayOfWeek (double Jul)
{
    int dow;
    int vTmp;
  
	vTmp = (int)(floor (Jul + 1));
	dow = vTmp % 7L;
 
    return (dow);
}

void ARM_Local_Date::JulianDateToDMY (void)
{
    int d;
    int m;
    int y;
  
    JulianToDMY (Julian, &d, &m, &y);
 
    Day = d;
    Month = m;
    Year = y;
 
    DayOfWeek = JulianToDayOfWeek (Julian);
}

double ARM_Local_Date::DMYToJulian(void)
{
    double julianDate;
 
    julianDate = DMYtoJulian (Day, Month, Year);
 
    return (julianDate);
}

double ARM_Local_Date::DMYtoJulian (int d, int m, int y)
{
    long Ja, Jm, Jy;
    double Jul;
    
    char buf[100];
 
    if(!CheckForValidDate (d, m, y))
    {
       sprintf(buf, "Invalid Date : %d-%d-%d", d, m, y);

	   //Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID, buf);
 
       //throw x;
 
       return (BADDATE);
    }
 
    if(y < 0)
    {
       y++;
    }
 
    if(m > 2)
    {
       Jy = y;
       Jm = m + 1;
    }
    else
    {
       Jy = y - 1;
       Jm = m + 13;
    }
 
    Jul = (floor (365.25 * (double)Jy) + floor (30.6001 * (double)Jm) + d + 1720995L);
 
    if(d + 31L * (m + 12L * y) >= IGREG)
    {
       Ja = (long) (0.01*Jy);
       Jul += (2 - Ja + (long)(0.25 * Ja));
    }

    return (Jul);
}

void ARM_Local_Date::JulianToDMY (double Jul, int* d, int* m, int* y)
{
    long Ja, JAlpha, Jb, Jc, Jd, Je;
  
    if((Jul != BADDATE) && (Jul >= MINDATE) && (Jul <= MAXDATE))
    {
		if (Jul >= GREGOR)
		{
			JAlpha = (long)(((double)(Jul - 1867216L) - 0.25) / 36524.25);
			Ja = (long)(Jul + 1 + JAlpha - (long)(0.25 * JAlpha));
		}
		else
		{
			Ja = (long)Jul;
		}
 
		Jb = Ja + 1524;
		Jc = (long)(6680.0 + ((double)(Jb - 2439870L) - 122.1) / 365.25);
		Jd = (long)(365 * Jc + (0.25 * Jc));
		Je = (long)((Jb - Jd) / 30.6001);
		*d = (int)(Jb - Jd - (int)(30.6001 * Je));
		*m = (int)Je - 1;
 
		if(*m > 12)
		{
			*m -= 12;
		}
 
		*y = (int)(Jc - 4715);
 
		if(*m > 2)
		{
			--(*y);
		}
 
		if(*y <= 0)
		{
			--(*y);
		}
	}
	else
	{
		//MSG_printf_message (MSG_ERROR, "Invalid Date : Converting Julian Date to dd/mm/yyyy failed");

		/*Exception x(__LINE__, __FILE__, ERR_DATE_NOT_VALID,
			"Invalid Date : Converting Julian Date to dd/mm/yyyy failed");
 
		throw x;*/
	}
}

int ARM_Local_Date::DaysInMonth (int m, int y)
{
    if((m == 2) && IsLeap (y))
    {
		return (29);
    }
    else
    {
       return (LocalMonthDays[m - 1]);
    }
}

int ARM_Local_Date::CheckForValidDate (int d, int m, int y)
{
    if(y / 1960.0 < 1.0)
    {
       return (0);
    }

    if((y < MINYEAR) || (y > MAXYEAR) || (y == 0)
		|| (m < MINMONTH) || (m > MAXMONTH))
    {
       return(0);
    }
 
    return (d <= DaysInMonth (m, y));
}

ARM_Local_Date::ARM_Local_Date (double julianDate)
{
    //SetName(ARM_Local_Date);
 
    Julian = julianDate;
 
    JulianDateToDMY();
}

ARM_Local_Date::ARM_Local_Date () : Julian (0.0), Day (0), Month (0), Year (0), DayOfWeek (0)
{
    //SetName(ARM_Local_Date);
}

ARM_Local_Date::ARM_Local_Date (int d, int m, int y)
{
    //SetName(ARM_Local_Date);
 
    ChangeDate (d, m, y);
}

void ARM_Local_Date::ChangeDate (int d, int m, int y)
{
    double julianDate;
 
     julianDate = DMYtoJulian (d, m, y);
 
    if (julianDate != BADDATE)
    {
       Julian = julianDate;
       Day = d;
       Month = m;
       Year = y;
 
       DayOfWeek = JulianToDayOfWeek (Julian);
    }
    else
    {
    }
}
