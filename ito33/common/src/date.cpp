/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/date.cpp
// Purpose:     implementation of Date class
// Author:      Vadim Zeitlin
// Created:     03.09.02
// RCS-ID:      $Id: date.cpp,v 1.25 2006/06/15 17:15:40 wang Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/error.h"

#include "ito33/date.h"

#ifdef _WIN32
  #include "ito33/win32/winwrap.h"
  #include <ole2.h>
#endif

#include <ctype.h>
#include <time.h>
#include <math.h>

extern const ito33::Error ITO33_BAD_PARAM,
                          ITO33_BAD_DATE,
                          ITO33_INVALID_DAYCOUNTCONVENTION;

namespace ito33
{

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

// this array contains the cumulated number of days in all previous months for
// normal and leap years
static const Date::Day_t gs_cumulatedDays[2][Date::MONTHS_IN_YEAR] =
{
  { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 },
  { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 }
};

// the max JDN we support: it could be greater if the algorithm in GetTm() were
// slightly modified to avoid overflow when calculating 'n' there but this is
// enough for our purposes
static const unsigned long JDN_MAX = ULONG_MAX / 4 - 68569;

// ============================================================================
// helper functions
// ============================================================================

// ============================================================================
// Date implementation
// ============================================================================

// ----------------------------------------------------------------------------
// utility functions
// ----------------------------------------------------------------------------

/* static */
void Date::ThrowBadDate()
{
  throw EXCEPTION(ITO33_BAD_DATE);
}

void ThrowBadDCC()
{
  throw EXCEPTION(ITO33_INVALID_DAYCOUNTCONVENTION);
}

void Date::Check() const
{
  if ( !IsValid() )
    ThrowBadDate();
}

// ----------------------------------------------------------------------------
// calendar functions
// ----------------------------------------------------------------------------

// get the number of days in the given month of the given year
/* static */
Date::Day_t Date::GetNumOfDaysInMonth(Date::Year_t year, Date::Month month)
{
  // the number of days in month in Julian/Gregorian calendar: the first line
  // is for normal years, the second one is for the leap ones
  static Day_t daysInMonth[2][MONTHS_IN_YEAR] =
  {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 },
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }
  };

  return daysInMonth[IsLeapYear(year)][month - 1];
}

Date& Date::SetTM(const struct tm& tm)
{
  return Set(1900 + tm.tm_year, static_cast<Month>(Jan + tm.tm_mon), tm.tm_mday);
}

Date& Date::SetTicks(time_t t)
{
#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  struct tm local;
  localtime_s(&local, &t);

  return SetTM(local);
#else
  struct tm *tm = localtime(&t);

  return SetTM(*tm);
#endif
}

/* static */
Date Date::Today()
{
  Date date;

  return date.SetTicks(time(NULL));
}

// ----------------------------------------------------------------------------
// JDN <-> (year, month, day)
// ----------------------------------------------------------------------------

/*
   the algorithm used here to convert between JDN and the Gregorian calender
   is due to Henry F. Fliegel and Thomas C. Van Flandern, published in
   Communications of the ACM (CACM, volume 11, number 10, October 1968, p.657))
 */

Date::Tm Date::GetTm() const
{
  unsigned long l = GetJDN();

  if ( l > JDN_MAX ) {
    // JDN out of range, it would result in an overflow below
    throw EXCEPTION(ITO33_BAD_PARAM);
  }

  l += 68569l;
  unsigned long n = ( 4 * l ) / 146097l;
  l -= ( 146097l * n + 3 ) / 4;
  unsigned long i = ( 4000 * ( l + 1 ) ) / 1461001l;
  l -= ( 1461 * i ) / 4 - 31;
  unsigned long j = ( 80 * l ) / 2447;

  Tm tm;
  tm.mday = (Day_t)(l - ( 2447 * j ) / 80);
  l = j / 11;
  tm.mon = (Month)(j + 2 - ( 12 * l ));
  tm.year = (Year_t)(100 * ( n - 49 ) + i + l);

  return tm;
}

Date& Date::Set(Date::Year_t year, Date::Month month, Date::Day_t day)
{
  if ( (day <= 0) || (day > GetNumOfDaysInMonth(year, month)) ) {
    throw EXCEPTION(ITO33_BAD_PARAM);
  }

  const long m = (month - 14) / 12;

  return SetJDN(( 1461 * ( year + 4800 + m) ) / 4 +
         ( 367 * ( month - 2 - 12 * m ) ) / 12 -
         ( 3 * ( ( year + 4900 + m ) / 100 ) ) / 4 +
         day - 32075);
}

Date& Date::SetYearDay(Year_t year, Day_t day)
{
  *this = Date(year, Jan, 1);

  return AddDays(day - 1);
}

Date::Day_t Date::GetDayOfYear() const
{
  const Tm tm = GetTm();

  return (Day_t)(gs_cumulatedDays[IsLeapYear(tm.year)][tm.mon - 1] + tm.mday);
}

Date::Day_t Date::GetWeekOfYear(Date::WeekFlags flags) const
{
  if ( flags == Default_First ) {
    // TODO: don't always assume USA, query the system
    flags = Sunday_First;
  }

  // int, not Day_t as it may become negative below
  int yday = GetDayOfYear();

  // Jan 4 of the same year
  Date dateJan4;
  dateJan4.SetJDN(GetJDN() - yday + 4);

  int wd = dateJan4.GetDayOfWeek();
  if ( flags == Monday_First ) {
    if ( wd == Sun )
      wd = 6;
    else
      wd--;
  }

  // the first ISO week always contains Jan 4, so offset by it to be sure
  // that we're always in the first week mod 7
  yday -= 4;
  Day_t week = (Day_t)(yday >= 0 ? yday / DAYS_PER_WEEK + 1 : 0);
  if ( yday % 7 + wd >= DAYS_PER_WEEK )
    week++;

  return week;
}

time_t Date::GetTicks() const
{
    // ticks is the number of seconds since the Unix epoch
    return DaysDiff(Date(1970, Jan, 1), *this) * 60 * 60 * 24;
}

// ----------------------------------------------------------------------------
// to/from Excel serial value
// ----------------------------------------------------------------------------

/*
   Excel uses the number of days since Jan 1, 1900 + 1 to represent dates, i.e.
   Excel day 1 == Jan 1, 1900.

   There is an added twist, however: Excel date code has an intentional bug
   and considers the year 1900 as a leap one (which it wasn't). We must account
   for this here to get correct results.
 */

// JDN of Dec 31, 1899
static const unsigned long JDN_EXCEL_EPOCH = 2415020ul;

// number of days between Jan 1 and Mar 1, 1900 (date from which Excel starts
// to count correctly)
static const unsigned long EXCEL_CUTOFF = 59;

Date& Date::SetExcel(unsigned long serial)
{
  if ( serial > EXCEL_CUTOFF ) {
    // Excel has counted an extra (nonexisting) day
    serial--;
  }

  return SetJDN(JDN_EXCEL_EPOCH + serial);
}

unsigned long Date::GetExcel() const
{
  unsigned long jdn = GetJDN();
  if ( jdn < JDN_EXCEL_EPOCH ) {
    throw EXCEPTION(ITO33_BAD_PARAM);
  }

  jdn -= JDN_EXCEL_EPOCH;

  if ( jdn > EXCEL_CUTOFF ) {
    // do as Excel, count Feb 29, 1900
    jdn++;
  }

  return jdn;
}

// ----------------------------------------------------------------------------
// Calendar calculations
// ----------------------------------------------------------------------------

Date::Tm Date::DoAddMonths(int months) const
{
  Tm tm = GetTm();

  // we want month to be 0-based here, hence subtract Jan
  int month = tm.mon - Jan + months;

  tm.year += month / MONTHS_IN_YEAR;
  month %= MONTHS_IN_YEAR;

  // special adjustment for negative arguments
  if ( month < 0 ) {
    tm.year--;
    month += MONTHS_IN_YEAR;
  }

  // now make it 1-based as usual
  tm.mon = (Month)(month + Jan);

  return tm;
}

Date& Date::EndOfMonth(int months)
{
  // this functions works in the same way as Excel EOMONTH()
  const Tm tm = DoAddMonths(months);

  return Set(tm.year, tm.mon, GetNumOfDaysInMonth(tm.year, tm.mon));
}

Date& Date::AddMonths(int months)
{
  // this function works in the same way as Excel EDATE(), i.e. it returns
  // the closest date to this one in the month which is monthsDiff from now
  Tm tm = DoAddMonths(months);

  const Day_t dayLast = GetNumOfDaysInMonth(tm.year, tm.mon);
  if ( tm.mday >= dayLast )
    tm.mday = dayLast;

  return SetTm(tm);
}

/* static */
long Date::DaysDiffWithDayCount(const Date& date1,
                                const Date& date2,
                                DayCountConvention dcc)
{
  date1.Check();
  date2.Check();

  Date
    firstDate,
    secondDate;

  // Indicates if date1 <= date2
  bool bPositive;

  if (date1 == date2)
  {
    return 0;
  }
  else if (date1 < date2)
  {
    bPositive = true;
    firstDate = date1;
    secondDate = date2;
  }
  else
  {
    bPositive = false; 
    firstDate = date2;
    secondDate = date1;
  }

  long iNbDays;

  switch ( dcc ) 
  {
    case DayCountConvention_30360:
    case DayCountConvention_30360_NO_EOM:
      {
        const Tm 
          tm1 = firstDate.GetTm(),
          tm2 = secondDate.GetTm();

        Day_t 
          day1 = tm1.mday,
          day2 = tm2.mday;

        if ( day1 == 31 ) 
          day1 = 30;

        if ( day2 == 31 && day1 == 30 ) 
          day2 = 30;

        iNbDays = (day2 - day1) + 30 * (tm2.mon - tm1.mon) 
                + 360 * (tm2.year - tm1.year);
      }
      break;

    case DayCountConvention_30E360:
    case DayCountConvention_30E360_NO_EOM:
      {
        const Tm 
          tm1 = firstDate.GetTm(),
          tm2 = secondDate.GetTm();

        Day_t 
          day1 = tm1.mday,
          day2 = tm2.mday;

        if ( day1 == 31 )
          day1 = 30;

        if ( day2 == 31 )
          day2 = 30;

        iNbDays = (day2 - day1) + 30 * (tm2.mon - tm1.mon) 
                + 360 * (tm2.year - tm1.year);
      }
      break;

    case DayCountConvention_30U360:
    case DayCountConvention_30U360_NO_EOM:
      {
        const Tm 
          tm1 = firstDate.GetTm(),
          tm2 = secondDate.GetTm();

        Day_t 
          day1 = tm1.mday,
          day2 = tm2.mday;

        if ( day1 == 31 ) 
          day1 = 30;

        if ( day2 == 31 && day1 == 30 ) 
          day2 = 30;

        // Treatment for february (this is the difference with the classical
        // case of the DayCountConvention_30360)
        Date dateTmp = firstDate.EndOfMonth(0);
        if ( tm1.mon == Feb && day1 == dateTmp.GetDay() )
        {
          day1 = 30;
          
          dateTmp = secondDate.EndOfMonth(0);
          if ( tm2.mon == Feb && day2 == dateTmp.GetDay() )
            day2 = 30;
        }

        iNbDays = (day2 - day1) + 30 * (tm2.mon - tm1.mon) 
                + 360 * (tm2.year - tm1.year);
      }
      break;
      
    case DayCountConvention_ActAct:
    //case DayCountConvention_ActAct_ISDA:
    case DayCountConvention_Act360:
    case DayCountConvention_Act365:
    case DayCountConvention_Act365L:      
    case DayCountConvention_ActAct_NO_EOM:
    //case DayCountConvention_ActAct_ISDA_NO_EOM:
    case DayCountConvention_Act360_NO_EOM:
    case DayCountConvention_Act365_NO_EOM:
    case DayCountConvention_Act365L_NO_EOM:
      {
        iNbDays = DaysDiff(firstDate, secondDate);
      }
      break;

    default:
      throw EXCEPTION(ITO33_BAD_PARAM);
  }

  return ( bPositive ? iNbDays : -iNbDays );
}

/* static */
double Date::YearsDiff(const Date& date1,
                       const Date& date2,
                       DayCountConvention dcc)
{
  double fraction;

  switch ( dcc ) 
  {
    case DayCountConvention_30360:
    case DayCountConvention_30E360:
    case DayCountConvention_30U360:
    case DayCountConvention_Act360:
    case DayCountConvention_30360_NO_EOM:
    case DayCountConvention_30E360_NO_EOM:
    case DayCountConvention_30U360_NO_EOM:
    case DayCountConvention_Act360_NO_EOM:
      {
        fraction = DaysDiffWithDayCount(date1, date2, dcc) / 360.;
      }
      break;

    case DayCountConvention_Act365:
    case DayCountConvention_Act365_NO_EOM:
      {
        fraction = DaysDiffWithDayCount(date1, date2, dcc) / 365.;
      }
      break;

    case DayCountConvention_ActAct:
    //case DayCountConvention_ActAct_ISDA:
    case DayCountConvention_Act365L:
    case DayCountConvention_ActAct_NO_EOM:
    //case DayCountConvention_ActAct_ISDA_NO_EOM:
    case DayCountConvention_Act365L_NO_EOM:
      {
        const Year_t 
          y1 = date1.GetYear(),
          y2 = date2.GetYear();

        // micro optimization for the case of the same year (the "else"
        // branch would still give the same result but with slightly
        // more work)
        if ( y1 == y2 ) 
        {
          fraction = DaysDiff(date1, date2);
          fraction /= GetNumOfDaysInYear(y1);
        }
        else 
        { // different years
          // add the fraction of the days left in the first year and
          // the fraction of the days passed in the second one and
          // the number of full years in between
          const double 
            daysTotal1 = GetNumOfDaysInYear(y1),
            daysTotal2 = GetNumOfDaysInYear(y2);

          fraction = y2 - y1 - 1 
                   + (daysTotal1 - date1.GetDayOfYear()) / daysTotal1 
                   + date2.GetDayOfYear() / daysTotal2;
        }
      }
      break;

    default:
      throw EXCEPTION(ITO33_BAD_PARAM);
  }

  return fraction;
}

// ----------------------------------------------------------------------------
// formatting dates
// ----------------------------------------------------------------------------

// for GetWeekDayName, GetMonthName and GetWeekDayFromName, GetMonthFromName
enum NameFlags
{
  Name_Full = 0x01,       // return full name
  Name_Abbr = 0x02        // return abbreviated name
};

// this function is a wrapper around strftime(3)
static std::string CallStrftime(const char *format, const tm* tm)
{
  char buf[4096];
  if ( !strftime(buf, SIZEOF(buf), format, tm) ) {
    FAIL( "strftime() failed" );
  }

  return buf;
}

// fll the struct tm with default values
static void InitTm(struct tm& tm)
{
  // struct tm may have etxra fields (undocumented and with unportable
  // names) which, nevertheless, must be set to 0
  memset(&tm, 0, sizeof(struct tm));

  tm.tm_mday = 1;   // mday 0 is invalid
  tm.tm_year = 76;  // any valid year
  tm.tm_isdst = -1; // auto determine
}

static std::string GetMonthName(Date::Month month, NameFlags flags)
{
  CHECK( month != Date::Inv_Month, "", "invalid month" );

  // notice that we must set all the fields to avoid confusing libc (GNU one
  // gets confused to a crash if we don't do this)
  tm tm;
  InitTm(tm);
  tm.tm_mon = month - Date::Jan;      // tm_mon is 0-based

  return CallStrftime(flags == Name_Abbr ? "%b" : "%B", &tm);
}

static std::string GetWeekDayName(Date::WeekDay wday, NameFlags flags)
{
  CHECK( wday != Date::Inv_WeekDay, "", "invalid weekday" );

  // take some arbitrary Sunday
  tm tm;
  InitTm(tm);
  tm.tm_mday = 28;
  tm.tm_mon = Date::Nov - Date::Jan;  // tm_mon is 0-based
  tm.tm_year = 99;

  // and offset it by the number of days needed to get the correct wday
  tm.tm_mday += wday - Date::Sun;

  // call mktime() to normalize it...
  (void)mktime(&tm);

  // ... and call strftime()
  return CallStrftime(flags == Name_Abbr ? "%a" : "%A", &tm);
}

// returns the string corresponding to "%x" format specificator
static std::string FormatDefault(const Date& date)
{
  // TODO: use real current format, not hardcoded US one
  return date.Format("%Y/%m/%d");
}

std::string Date::Format(const char *format) const
{
  std::string res;
  CHECK( format, res, "NULL format in Date::Format()" );

  const Tm tm = GetTm();

  for ( const char *p = format; *p; p++ ) {
    if ( *p != '%' ) {
      // copy as is
      res += *p;

      continue;
    }

    // set the default format
    std::string fmt;
    switch ( *++p ) {
      case 'Y':               // year has 4 digits
        fmt = "%04d";
        break;

      case 'j':               // day of year has 3 digits
        fmt = "%03d";
        break;

      case 'w':               // week day as number has only one
        fmt = "%d";
        break;

      default:
        // it's either another valid format specifier in which case
        // the format is "%02d" (for all the rest) or we have the
        // field width preceding the format in which case it will
        // override the default format anyhow
        fmt = "%02d";
    }

    bool restart = true;
    while ( restart ) {
      restart = false;

      // start of the format specification
      switch ( *p ) {
        case 'a':       // a weekday name
        case 'A':
          // second parameter should be true for abbreviated names
          res += GetWeekDayName(GetDayOfWeek(),
                     *p == 'a' ? Name_Abbr : Name_Full);
          break;

        case 'b':       // a month name
        case 'B':
          res += GetMonthName(tm.mon,
                    *p == 'b' ? Name_Abbr : Name_Full);
          break;

        case 'x':       // locale default date representation
          res += FormatDefault(*this);
          break;

        case 'd':       // day of a month (01-31)
          res += String::Printf(fmt.c_str(), tm.mday);
          break;

        case 'j':       // day of the year
          res += String::Printf(fmt.c_str(), GetDayOfYear());
          break;

        case 'm':       // month as a number (01-12)
          res += String::Printf(fmt.c_str(), tm.mon);
          break;

        case 'U':       // week number in the year (Sunday 1st week day)
          res += String::Printf(fmt.c_str(), GetWeekOfYear(Sunday_First));
          break;

        case 'W':       // week number in the year (Monday 1st week day)
          res += String::Printf(fmt.c_str(), GetWeekOfYear(Monday_First));
          break;

        case 'w':       // weekday as a number (0-6), Sunday = 0
          res += String::Printf(fmt.c_str(), GetDayOfWeek());
          break;

        case 'y':       // year without century (00-99)
          res += String::Printf(fmt.c_str(), tm.year % 100);
          break;

        case 'Y':       // year with century
          res += String::Printf(fmt.c_str(), tm.year);
          break;

        default:
          // is it the format width?
          for ( fmt = "";
             *p == '-' || *p == '+' || *p == '#'
              || isdigit(*p); p++ ) {
            fmt += *p;
          }

          if ( !fmt.empty() ) {
            // we've only got the flags and width so far in fmt
            fmt = std::string("%") + fmt + 'd';

            restart = true;

            break;
          }

          // no, it wasn't the width
          FAIL("unknown format specificator");

          // fall through and just copy it nevertheless

        case '%':       // a percent sign
          res += *p;
          break;

        case 0:             // the end of string
          FAIL("missing format at the end of string");

          // just put the '%' which was the last char in format
          res += '%';
          break;
      }
    }
  }

  return res;
}

// ----------------------------------------------------------------------------
// parsing dates
// ----------------------------------------------------------------------------

// return the month if the string is a month name or Inv_Month otherwise
static Date::Month GetMonthFromName(const std::string& name, NameFlags flags)
{
  Date::Month mon;
  for ( mon = Date::Jan; mon < Date::Inv_Month; Date::NextMonth(mon) ) {
    if ( String::Stricmp(name, GetMonthName(mon, flags)) == 0 )
      break;
  }

  return mon;
}

// return the weekday if the string is a weekday name or Inv_WeekDay otherwise
static Date::WeekDay GetWeekDayFromName(const std::string& name, NameFlags flags)
{
  Date::WeekDay wd;
  for ( wd = Date::Sun; wd < Date::Inv_WeekDay; Date::NextWeekDay(wd) ) {
    if ( String::Stricmp(name, GetWeekDayName(wd, flags)) == 0 )
      break;
  }

  return wd;
}

// scans all digits (but no more than len) and returns the resulting number
static bool GetNumericToken(size_t len, const char*& p, unsigned long *number)
{
  size_t n = 1;
  std::string s;
  while ( isdigit(*p) ) {
    s += *p++;

    if ( len && ++n > len )
      break;
  }

  return !s.empty() && String::ToULong(s, number);
}

// scans all alphabetic characters and returns the resulting string
static std::string GetAlphaToken(const char*& p)
{
  std::string s;
  while ( isalpha(*p) ) {
    s += *p++;
  }

  return s;
}

const char *Date::Parse(const char *date, const char *format)
{
  CHECK( date && format, NULL, "NULL pointer in Date::Parse()" );

  Clear();

  std::string str;
  unsigned long num;

  // is the date already set?
  bool done = false;

  // what fields have we found?
  bool haveWDay = false,
    haveYDay = false,
    haveDay = false,
    haveMon = false,
    haveYear = false;

  // and the value of the items we have (init them to get rid of warnings)
  WeekDay wday = Inv_WeekDay;
  Day_t yday = 0,
     mday = 0;
  Date::Month mon = Inv_Month;
  Year_t year = 0;

  const char *input = date;
  for ( const char *fmt = format; *fmt; fmt++ ) {
    if ( *fmt != '%' ) {
      if ( isspace(*fmt) ) {
        // a white space in the format string matches 0 or more white
        // spaces in the input
        while ( isspace(*input) ) {
          input++;
        }
      }
      else  { // !space
        // any other character (not whitespace, not '%') must be
        // matched by itself in the input
        if ( *input++ != *fmt ) {
          // no match
          return NULL;
        }
      }

      // done with this format char
      continue;
    }

    // start of a format specification

    // parse the optional width
    size_t width = 0;
    while ( isdigit(*++fmt) ) {
      width *= 10;
      width += *fmt - '0';
    }

    // the default widths for the various fields
    if ( !width )
    {
      switch ( *fmt )
      {
        case ('Y'):               // year has 4 digits
            width = 4;
            break;

        case ('j'):               // day of year has 3 digits
        case ('l'):               // milliseconds have 3 digits
            width = 3;
            break;

        case ('w'):               // week day as number has only one
            width = 1;
            break;

        default:
            // default for all other fields
            width = 2;
      }
    }

    // then the format itself
    switch ( *fmt ) {
      case 'a':       // a weekday name
      case 'A':
        {
          NameFlags flag = *fmt == 'a' ? Name_Abbr : Name_Full;
          wday = GetWeekDayFromName(GetAlphaToken(input), flag);
          if ( wday == Inv_WeekDay ) {
            // no match
            return NULL;
          }
        }
        haveWDay = true;
        break;

      case 'b':       // a month name
      case 'B':
        {
          NameFlags flag = *fmt == 'b' ? Name_Abbr : Name_Full;
          mon = GetMonthFromName(GetAlphaToken(input), flag);
          if ( mon == Inv_Month ) {
            // no match
            return NULL;
          }
        }
        haveMon = true;
        break;

      case 'd':       // day of a month (01-31)
        if ( !GetNumericToken(width, input, &num) ||
            (num > 31) || (num < 1) ) {
          // no match
          return NULL;
        }

        // we can't check whether the day range is correct yet, will
        // do it later - assume ok for now
        haveDay = true;
        mday = (Day_t)num;
        break;

      case 'j':       // day of the year
        if ( !GetNumericToken(width, input, &num)
            || !num || (num > 366) ) {
          // no match
          return NULL;
        }

        haveYDay = true;
        yday = (Day_t)num;
        break;

      case 'm':       // month as a number (01-12)
        if ( !GetNumericToken(width, input, &num)
            || !num || (num > MONTHS_IN_YEAR) ) {
          // no match
          return NULL;
        }

        haveMon = true;
        mon = (Month)num;
        break;

      case 'w':       // weekday as a number (0-6), Sunday = 0
        if ( !GetNumericToken(width, input, &num) || (wday > 6) ) {
          // no match
          return (char *)NULL;
        }

        haveWDay = true;
        wday = (WeekDay)num;
        break;

      case 'x':       // locale default date representation
        // TODO: use real current format, not hardcoded US one
        input = Parse(input, "%Y/%m/%d");
        if ( !input )
          return NULL;

        // we're already set, no need to change anything
        done = true;
        break;

      case 'y':       // year without century (00-99)
        if ( !GetNumericToken(width, input, &num) || (num > 99) ) {
          // no match
          return NULL;
        }

        haveYear = true;

        // TODO should have an option for roll over date instead of
        //      hard coding it here
        year = (num > 30 ? 1900 : 2000) + (Year_t)num;
        break;

      case 'Y':       // year with century
        if ( !GetNumericToken(width, input, &num) ) {
          // no match
          return NULL;
        }

        haveYear = true;
        year = (Year_t)num;
        break;

      case '%':       // a percent sign
        if ( *input++ != '%' ) {
          // no match
          return NULL;
        }
        break;

      case 0:             // the end of string
        FAIL("unexpected format end");

        // fall through

      default:            // not a known format spec
        return NULL;
    }
  }

  if ( done ) {
    // we have already set this date to the correct value because of the
    // recursive call above
    return input;
  }

  // format matched, try to construct a date from what we have now
  if ( !haveYear ) {
    FAIL( "invalid format -- should include the year" );

    return NULL;
  }

  // TODO we don't check here that the values are consistent, if both year
  //      day and month/day were found, we just ignore the year day then
  if ( haveMon && haveDay ) {
    if ( mday > GetNumOfDaysInMonth(year, mon) ) {
      return NULL;
    }

    Set(year, mon, mday);
  }
  else if ( haveYDay ) {
    if ( yday > GetNumOfDaysInYear(year) ) {
      return (char *)NULL;
    }

    SetYearDay(year, yday);
  }
  else {
    FAIL( "invalid format -- should include either month or year day" );

    return NULL;
  }

  // check that the week day is consistent if we have it
  if ( haveWDay ) {
    if ( GetDayOfWeek() != wday ) {
      Clear();

      return NULL;
    }
  }

  return input;
}

// ----------------------------------------------------------------------------
// COM date format
// ----------------------------------------------------------------------------

#ifdef _WIN32

double Date::GetOleDate() const
{
  SYSTEMTIME st;
  ::ZeroMemory(&st, sizeof(st));
  st.wYear = static_cast<WORD>(GetYear());
  st.wMonth = static_cast<WORD>(GetMonth());
  st.wDayOfWeek = static_cast<WORD>(GetDayOfWeek());
  st.wDay = static_cast<WORD>(GetDay());

  double d;
  if ( !::SystemTimeToVariantTime(&st, &d) )
    throw EXCEPTION(ITO33_BAD_DATE);

  return d;
}

Date& Date::SetOleDate(double d)
{
  SYSTEMTIME st;
  if ( !::VariantTimeToSystemTime(d, &st) )
    throw EXCEPTION(ITO33_BAD_DATE);

  return Set(st.wYear, static_cast<Date::Month>(st.wMonth), st.wDay);
}

#endif // _WIN32

} // namesapce ito33
