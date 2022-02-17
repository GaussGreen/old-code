#ifndef ARM_Date_H
#define ARM_Date_H

#include <math.h>
#include <stdio.h>

//#include "CCmessage.h"
#include "armdef.h"

#define MODULO(x,y)		(( x >= 0 ) ? ( x % y ) : ((x+y*(1+(int)((-x)/y))) % y))
#define IsLeap(y)		((y % 4 == 0) && (y % 4000 != 0) && ((y % 100 != 0) \
							|| (y % 400 == 0)))
#define IGREG			(15+31L*(10+12L*1582))

#define MINMONTH		1
#define MAXMONTH		12
#define MINDAY			1
#define MAXDAY			31
#define MINWEEKDAY		0
#define MAXWEEKDAY		6
#define WEEKDAYS		7
#define YEARMONTHS		12
#define BADDATE			-1L
#define MINYEAR			-4700L
#define MAXYEAR			25000L
#define MINDATE			4749L      // 1 janvier 4700 av JC
#define MAXDATE			10852487L  // 31 decembre 25000
#define GREGOR			2299161L

#define ERR_DATE_NOT_VALID	115

class ARM_Local_Date
{
private:
	double Julian;
	int	Day;        
    int Month;      
    int Year;
	int	DayOfWeek;

public:
	ARM_Local_Date ();
	ARM_Local_Date (double julianDate);
	ARM_Local_Date::ARM_Local_Date (int d, int m, int y);

	ARM_Local_Date& AddPeriod (int period);
    ARM_Local_Date& AddPeriodMult(int period, int mult = 1);
	ARM_Local_Date& AddDays (int d);
	ARM_Local_Date& AddMonths(int nMonths, int GOTO_END_OF_MONTH = 0);
	int JulianToDayOfWeek (double Jul);
	void JulianDateToDMY (void);
	double DMYToJulian(void);
	double DMYtoJulian (int d, int m, int y);
	void JulianToDMY (double Jul, int* d, int* m, int* y);
	int DaysInMonth (int m, int y);
	int CheckForValidDate (int d, int m, int y);
	void ChangeDate (int d, int m, int y);

	inline int GetDay (void) 
	{ 
		return Day;
	}
    
    inline int GetMonth (void)
	{ 
		return Month; 
	}
        
	int GetYear (void)
	{ 
		return Year;
	}

	inline operator double()
	{ 
		return (Julian); 
	}
};

#endif	// ARM_Date_H