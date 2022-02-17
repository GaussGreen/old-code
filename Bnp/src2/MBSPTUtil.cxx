#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "MBSPTUtil.h"

#include "utdates.h"

// there is a local static array of ints flags
// if readOnly>0        , then just store flags[FlagNames] = Value and return
// value else read it off
int MBS_Flags(MBS_FlagNames FlagNames, int value, int readOnly) {
  static int mbs_flag_values[1];
  static int initial;
  if (initial == 0) {
    initial = 1;
    mbs_flag_values[_USESRTDATES] = 1;
  }
  if (!readOnly)
    mbs_flag_values[FlagNames] = value;
  return (mbs_flag_values[FlagNames]);
} // MBS_Flags

//***mbs_err_handler***
// for an error text        , if write==1        , try to open and append the
// text        , in any case        , return the error text exitmode = 1 , exit
//

char *mbs_err_handler(char *error_text, int write, int exit_mode) {
  FILE *outError = 0;
  char string[MAXSTRBUFSZ];

  if (write) {
    string[0] = '\0';
    strcat(string, TBAOUTPUTDIR);
    strcat(string, "ErrorMessages.txt");
    outError = fopen(string, "at");
  } // if

  if (outError) {
    fprintf(outError, "%s\n", error_text);
  } // if

  if (outError)
    fclose(outError);
  if (exit_mode)
    exit(EXIT_FAILURE);
  return (error_text);
} // mbs_err_handler

//*****  nrerror  ***********************************************************

char *nrerror(char *error_text) {
  return (mbs_err_handler(error_text, 0, 0));
} // nrerror

//***mbspt_lvector***
// allocate memory for an array of longs
// starting with index start        , ending with end        , inclusive

long *mbspt_lvector(int start, int end)

{
  long *v;

#ifdef CHECK_MEMALLOC_LIM
  if (end < start)
    return 0;
#endif // CHECK_MEMALLOC_LIM

  v = (long *)calloc((unsigned)(end - start + 1), sizeof(long));
  if (!v) {
    mbs_err_handler("Can't calloc lvector", MBS_WRITE_FILE, 1);
    return (0);
  } // if

  return (v - start);

} // mbspt_lvector */

//***mbspt_ivector***
// allocate memory for an array of integers
// starting with index start        , ending with end        , inclusive

int *mbspt_ivector(int start, int end)

{
  int *v;

#ifdef CHECK_MEMALLOC_LIM
  if (end < start)
    return 0;
#endif // CHECK_MEMALLOC_LIM

  v = (int *)calloc((unsigned)(end - start + 1), sizeof(int));
  if (!v) {
    mbs_err_handler("Can't calloc ivector", MBS_WRITE_FILE, 1);
    return (0);
  } // if

  return (v - start);

} // mbspt_ivector */

//***mbspt_dvector***
// allocate memory for an array of doubles
// starting with index start        , ending with end        , inclusive

double *mbspt_dvector(int start, int end) {

  double *v;

#ifdef CHECK_MEMALLOC_LIM
  if (end < start)
    return 0;
#endif // CHECK_MEMALLOC_LIM

  v = (double *)calloc((unsigned)(end - start + 1), sizeof(double));
  if (!v) {
    mbs_err_handler("Can't calloc dvector", MBS_WRITE_FILE, EXIT_FAILURE);
    return (0);
  } // if

  return (v - start);

} // mbspt_dvector

//****mbspt_dmatrix***
// Allocation of memory for a matrix of doubles [rowStart        , colStart] x
// [colStart        , colEnd]

double **mbspt_dmatrix(int rowStart, int rowEnd, int colStart, int colEnd) {
  int i;
  double **m;

#ifdef CHECK_MEMALLOC_LIM
  if ((colStart < rowStart) || (colEnd < colEnd))
    return 0;
#endif // CHECK_MEMALLOC_LIM

  m = (double **)calloc((unsigned)(rowEnd - rowStart + 1), sizeof(double *));
  if (!m) {
    mbs_err_handler("Can't calloc (double **) for dmatrix", MBS_WRITE_FILE,
                    EXIT_FAILURE);
    return (0);
  } // if
  m -= rowStart;

  for (i = rowStart; i <= rowEnd; i++) {
    m[i] = (double *)calloc((unsigned)(colEnd - colStart + 1), sizeof(double));
    if (!m[i]) {
      mbs_err_handler("Can't calloc *(double *) in dmatrix", MBS_WRITE_FILE,
                      EXIT_FAILURE);
      return (0);
    } // if
    else
      m[i] -= colStart;
  } // for i

  return (m);

} // mbspt_dmatrix

//****mbspt_free_ivector***
// Free array of int v [start x end]

void mbspt_free_lvector(long *v, int start, int end) {
  free((long *)(v + start));

} // mbspt_free_lvector

//****mbspt_free_ivector***
// Free array of int v [start x end]

void mbspt_free_ivector(int *v, int start, int end) {
  free((char *)(v + start));

} // mbspt_free_ivector

//****mbspt_free_dvector***
// Free array of double v [start x end]

void mbspt_free_dvector(double *v, int start, int end) {
  free((char *)(v + start));
} // mbspt_dvector

/***mbsptfree_dmatrix***
//Free matrix of doubles m [rowStart        , colStart] x [colStart        ,
colEnd]
*/
void mbspt_free_dmatrix(double **m, int rowStart, int rowEnd, int colStart,
                        int colEnd) {
  int i;

  for (i = rowEnd; i >= rowStart; i--)
    free((char *)(m[i] + colStart));

  free((char *)(m + rowStart));

} // mbspt_free_dmatrix

DATE_UNIT *mbspt_DateUnitVector(int start, int end) {
  DATE_UNIT *v;

#ifdef CHECK_MEMALLOC_LIM
  if (end < start)
    return 0;
#endif // CHECK_MEMALLOC_LIM

  v = (DATE_UNIT *)calloc((unsigned)(end - start + 1), sizeof(DATE_UNIT));
  if (!v) {
    mbs_err_handler("Can't calloc DateUnitVector", MBS_WRITE_FILE, 1);
    return (0);
  } // if

  return (v - start);
} // mbspt_DateUnitVector

void mbspt_free_DateUnitVector(DATE_UNIT *v, int start, int end) {
  free((char *)v);
}

// code must have '/' like "ACT/ACT"
char *readBasis(YBASIS *ybasis, DIFFBASIS *diffbasis, char *code) {
  char *ptr;
  char *err = 0;

  ptr = strchr(code, '/');
  if (!ptr) {
    *ybasis = _BADYBASIS;
    *diffbasis = _BADDIFFBASIS;
    return ("Invalid basis");
  }
  //
  *ptr = '\0';
  *diffbasis = readDIFFBASIS(code);
  *ybasis = readYBASIS(ptr + 1);
  *ptr = '/';
  if (*ybasis == _BADYBASIS || *diffbasis == _BADDIFFBASIS)
    return ("Bad basis");
  //
  return (err);

} // readBasis

YBASIS readYBASIS(char *code) {

  if (strcmp("ACT", code) == 0)
    return (_YACT);
  if (strcmp("360", code) == 0)
    return (_360);
  if (strcmp("365", code) == 0)
    return (_365);
  return (_BADYBASIS);
} // readYBASIS

DIFFBASIS readDIFFBASIS(char *code) {
  if (strcmp("ACT", code) == 0)
    return (_DACT);
  if (strcmp("30", code) == 0)
    return (_30);
  return (_BADDIFFBASIS);
} // readDIFFBASIS

void Dsplit(long d, long *yy_o, long *mm_o, long *dd_o) {
  if (MBS_Flags(_USESRTDATES, 0, 1)) {
    *yy_o = year(d);
    *mm_o = month(d);
    *dd_o = day(d);
  } else {
    *yy_o = d / 10000;
    *mm_o = (d - *yy_o * 10000) / 100;
    *dd_o = d - *yy_o * 10000 - *mm_o * 100;
  }
} // Dsplit

//***month_diff
// get num of months between the first of 2 two dates' months

int month_diff(long d1, long d2) {
  long Y1, Y2, M1, M2, dummy;
  Dsplit(d1, &Y1, &M1, &dummy);
  Dsplit(d2, &Y2, &M2, &dummy);
  return (12 * (Y2 - Y1) + (M2 - M1));
} // month_diff

//***Isleap***
// Returns zero if 4 digit argument is not a leap year
// otherwise returns 1.
//
int Isleap(long year) {
  if (MBS_Flags(_USESRTDATES, 0, 1))
    return (leap_year(year));
  else {
    if (year % 4 != 0)
      return (0); // not divisible by 4
    if (year % 100 != 0)
      return (1); // divisible by 4 but not 100
    if (year % 400 != 0)
      return (0); // divisible by 100 but not 400
    return (1);   // divisible by 400        , so is a leap year
  }
} // Isleap

// Dayofwk
// Calculates day of week (0-6) of YYYYMMDD formatted date: Sunday <-> 0 ,
// Saturday 6
long Dayofwk(long d) {
  static long noleap[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  static long leap[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  long mm, dd, yy;
  long *daysin;

  if (MBS_Flags(_USESRTDATES, 0, 1))
    return week_day(d);
  else {
    // FACT: 1/1/2000 is Sunday (Excel varifies it)        , between 20000101
    // and 00000101 there are 2000x365 days        , plus (2000/4 - 20 + 5) leap
    // days        , mod 7 gives 0        , i.e.        , 00000101 is Sunday

    // last days on consecutive years        , if the first year is non-leap ,
    // are off by one week day; if leap        , by 2 days e.g.        ,
    // 12/31/2001 is Monday while 12/31/2002 is Tuesday

    Dsplit(d, &yy, &mm, &dd);
    daysin = Isleap(yy) ? leap : noleap;

    while (mm > 1)
      dd += daysin[--mm];
    // we have YYYYMMDD - last day of previous year
    if (yy > 0) // whether we count from from 19000101 or 00000101 makes no
                // difference
    {
      --yy;
      dd += yy;                           // number of years
      dd += yy / 4 - yy / 100 + yy / 400; // number of leap years
    }
    return (dd % 7);
  }
} // Dayofwk

//  Converts a 2 digit year (84) to a 4 digit year (1984);
//  assumes that 2 digit years less than 50 are after year 2000
long YY2YYYY(long year_i) {

  long y4_o, y2;

  y2 = year_i;
  if (y2 >= 100)
    y4_o = y2;
  else {
    if (y2 <= 50)
      y2 = y2 + 100;
    y4_o = y2 + 1900;
  }
  return (y4_o);
} // YY2YYYYY

//  Packs mm_i        ,dd_i        ,yy_i into YYYYMMDD format
//  calls YY2YYYY on yy_i before packing
long DateYMMDD(long Y, long MM, long DD) {
  if (MBS_Flags(_USESRTDATES, 0, 1))
    return date(DD, MM, YY2YYYY(Y));
  else
    return (10000 * YY2YYYY(Y) + 100 * MM + DD);
} // DateYMMDD

//*  Returns a date in YYYYMMDD format based on moving
//*  forward or backwards a number of calendar days
long Nxtday(long date, long days) {
  static long noleap[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  static long leap[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  long nxtm, nxtd, nxty, edays;
  long *daysin;
  if (MBS_Flags(_USESRTDATES, 0, 1))
    return (date + days);
  else {
    Dsplit(date, &nxty, &nxtm, &nxtd);
    edays = nxtd + days; // positive or negative
    daysin = Isleap(nxty) ? leap : noleap;
    if (days >= 0) // move forward
    {
      while (edays > daysin[nxtm]) {
        edays -= daysin[nxtm];
        nxtm++;
        // see if into next year
        if (nxtm > 12) {
          nxtm = 1;
          nxty++;
          daysin = Isleap(nxty) ? leap : noleap;
        }
      }    // nxtm        , edays        , nxty contain the proper date
    } else // move backwards
    {
      while (edays < 1) {
        nxtm--;
        if (nxtm < 1) {
          nxty--;
          nxtm = 12;
          daysin = Isleap(nxty) ? leap : noleap;
        }
        edays += daysin[nxtm];
      } // nxtm        , edays        , nxty contain the proper date
    }
    // return result
    return (DateYMMDD(nxty, nxtm, edays));
  }
} // Nxtday

//***Nxtwkday***
// This routine returns a new date after advancing or going
// backward a number of weekdays (Monday - Friday).
// The format for date is YYYYMMDD.
//
long Nxtwkday(long date, long advance) {
  long dow, work_date;

  if (MBS_Flags(_USESRTDATES, 0, 1))
    return add_unit(date, advance, SRT_BDAY /*no holiday!!*/,
                    NO_BUSDAY_CONVENTION);
  else {
    work_date = date;
    if (advance == 0)
      return (date);
    else if (advance < 0) // from date        , go backward
    {
      while (advance != 0) {
        dow = Dayofwk(work_date);
        advance++;
        if (dow == 1) // Monday
          work_date = Nxtday(work_date, (long)-3);
        else if (dow == 0) // Sunday
          work_date = Nxtday(work_date, (long)-2);
        else /* Tuesday - Saturday */
          work_date = Nxtday(work_date, (long)-1);
      }
    } else // from date        , go forward
    {
      while (advance != 0) {
        dow = Dayofwk(work_date);
        advance--;
        if (dow == 5) // Friday
          work_date = Nxtday(work_date, (long)3);
        else if (dow == 6) // Saturday
          work_date = Nxtday(work_date, (long)2);
        else // Sunday - Thursday
          work_date = Nxtday(work_date, (long)1);
      }
    }
    return (work_date);
  }
} // Nxtwkday

//***Nxtmth***
//*  Returns a date in YYYYMMDD format based on moving
//*  forward or backwards (mths) from (date).
//*  if eom > 0 and date is the last date in a month
//*  then the date returned will be the last day of the month
//*  (mths) before or after (date)
//*  a positive entry for mths => move forward
//*  a negative entry for mths => move backward.
//**************/

long Nxtmth(long date_, long mths, long eom) {
  static long noleap[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  static long leap[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  long *daysin, *daysout;
  long mm, dd, yy;
  long sign, nxtmm, nxtdd, nxtyy;
  long nxtdte;

  Dsplit(date_, &yy, &mm, &dd);

  if (mths > 0)
    sign = 1;
  else
    sign = -1;
  /* move forward by mths */
  nxtmm = mm + mths;
  nxtyy = yy;
  while (nxtmm > 12 || nxtmm < 1) {
    nxtmm -= 12 * sign;
    nxtyy += sign;
  }
  // now nxtyy and nxtmm have appropriate values
  // see if have to adjust days for end of month

  daysin = Isleap(yy) ? leap : noleap;
  daysout = Isleap(nxtyy) ? leap : noleap;
  nxtdd = dd;

  if (eom) {
    if ((dd == daysin[mm]) || (nxtdd > daysout[nxtmm]))
      nxtdd = daysout[nxtmm];
  } else {
    if (nxtdd > daysout[nxtmm]) {
      nxtdd = nxtdd - daysout[nxtmm];
      ++nxtmm;
      //  Year will never be incremented because
      //  December has 31 days.
    }
  }

  if (MBS_Flags(_USESRTDATES, 0, 1))
    return (date(nxtdd, nxtmm, nxtyy));
  else {
    nxtdte = nxtyy * 10000 + nxtmm * 100 + nxtdd;
    return (nxtdte);
  }
} // Nxtmth

char *HolidaysConstruct(Holidays *holidays, char *dir) {
  static char *fileNames = "Holiday.dat";

  FILE *stream = 0;
  char string[MAXSTRBUFSZ];

  string[0] = '\0';
  strcat(string, dir);
  strcat(string, "Holiday.dat");
  stream = fopen(string, "r");

  if (!stream)
    return (mbs_err_handler("Can't open HOLIDAY.dat", 0, 0));

  holidays->numHolidays = 0;
  while (fgets(string, 80, stream)) {
    if (holidays->numHolidays == MAXHOLIDAYS) {
      fclose(stream);
      return ("Too Many Holidays");
    }
    string[8] = '\0';
    holidays->Holidays[holidays->numHolidays++] = YYYYMMDD2XlDate(atol(string));
  }
  fclose(stream);
  return (0);
} // HolidaysConstruct

//***isHoliday***
//*This routine checks whether the date passed in is a holiday
//*or not according to the HOLIDAYS file. It returns (0) if
//*it is not and (1) if it is. The format is YYYYMMDD.
//
int isHoliday(long date_i) {
  static Holidays holiday;
  static int initialized;
  int i;

  if (!initialized) {
    initialized = 1;
    HolidaysConstruct(&holiday, TBADATADIR);
  }

  if (MBS_Flags(_USESRTDATES, 0, 1)) {
    for (i = 0; i < holiday.numHolidays; ++i)
      if (holiday.Holidays[i] == date_i)
        return (TRUE);
  } else {
    for (i = 0; i < holiday.numHolidays; ++i)
      if (holiday.Holidays[i] == YYYYMMDD2XlDate(date_i))
        return (TRUE);
  }

  return (FALSE);
} // isHoliday

int isBusDay(long date_i) {
  int dd = Dayofwk(date_i);

  if (MBS_Flags(_USESRTDATES, 0, 1)) {
    if (dd == SAT || dd == SUN)
      return (0); // if Sunday or Sat.
  } else {
    if (dd == 6 || dd == 0)
      return (0); // if Sunday or Sat.
  }
  if (isHoliday(date_i))
    return (0); // if Holiday
  return (1);
} // isBusDay

//*  This routine returns a new date after advancing or going
//*  backward a number of business days (Monday - Friday &
//*  holidays). It checks for holidays in HOLIDAY file.
//*  The format for date is YYYYMMDD.
//*/
long MBSNxtbusday(long date, long advance) {
  long work_date;

  work_date = date;
  if (!MBS_Flags(_USESRTDATES, 0, 1)) {
    if (advance == 0)
      return (date);
    else if (advance < 0) // go backward
    {
      while (advance != 0) {
        advance++;
        work_date = Nxtwkday(work_date, (long)-1);
        while (isHoliday(work_date))
          work_date = Nxtwkday(work_date, (long)-1);
      }
    } else // from date        , go forward
    {
      while (advance != 0) {
        advance--;
        work_date = Nxtwkday(work_date, (long)1);
        while (isHoliday(work_date))
          work_date = Nxtwkday(work_date, (long)1);
      }
    }
    return (work_date);
  } else {
    work_date = Nxtwkday(work_date, (long)0);
    while (isHoliday(work_date))
      work_date = Nxtwkday(work_date, (long)1);
    if (advance == 0)
      return work_date;
    else if (advance < 0) // go backward
    {
      while (advance != 0) {
        advance++;
        work_date = Nxtwkday(work_date, (long)-1);
        while (isHoliday(work_date))
          work_date = Nxtwkday(work_date, (long)-1);
      }
    } else // from date        , go forward
    {
      while (advance > 0) {
        advance--;
        work_date = Nxtwkday(work_date, (long)1);
        while (isHoliday(work_date))
          work_date = Nxtwkday(work_date, (long)1);
      }
    }
    return (work_date);
  }
} // MBSNxtbusday

// Date_Add
// if unit = day        , bus is needed: if bus        , advance by len bus-days
// , else by len calendar days if unit = month or year        , BOTH bus and eom
// needed: first advance by suitale number of months        ,
//             and eom affects how to handle end of month
//			then        , if bus        , then move to the next bus day
//if it is not a bus day
long Date_Add(long Ref_Date, DATE_UNIT unit, int len, int bus, int eom) {
  long ret;
  long mo;

  if (unit == DAY) {
    if (bus)
      return (MBSNxtbusday(Ref_Date, len));
    else
      return (Nxtday(Ref_Date, len));
  } else {
    mo = len;
    if (unit == YEAR)
      mo *= 12;
    ret = Nxtmth(Ref_Date, mo, eom);
    if (bus && !isBusDay(ret))
      ret = MBSNxtbusday(ret, 1);
    return (ret);
  }
} // Date_Add

long YYYYMMDD2XlDate(long d) {
  long yy, mm, dd;
  yy = d / 10000;
  mm = d - yy * 10000;
  mm = mm / 100;
  dd = d - mm * 100 - yy * 10000;
  return date(dd, mm, yy);
} // YYYYMMDD2XlDate

long fromYYYYMMDD(long d) {
  if (MBS_Flags(_USESRTDATES, 0, 1)) {
    return YYYYMMDD2XlDate(d);
  } else {
    return d;
  }

} // fromYYYYMMDD

long YYYYMM2FirstOfMon(long d) { return (fromYYYYMMDD(d * 100 + 1)); }

long DateDiff(long date1_i, long date2_i, DIFFBASIS base) {
  if (!MBS_Flags(_USESRTDATES, 0, 1)) {

    date1_i = YYYYMMDD2XlDate(date1_i);
    date2_i = YYYYMMDD2XlDate(date2_i);
  }
  if (base == _DACT)
    return (day_count_date(date1_i, date2_i, BASIS_ACT_USD));
  else
    return (day_count_date(date1_i, date2_i, BASIS_30_360));
}

//***Dateok***
//  Checks for validity of YYYYMMDD formatted dates
// ifndef USESRTDATES
// returns 0 if all is well
//  returns 1 if bad year
//  returns 2 if bad month
// returns 3 if bad day
// else return 0 is valid 1 if not
//
int Dateok(long date) {
  static long noleap[13] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  static long leap[13] = {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  long *daysin;
  long mm, dd, yy;
  if (MBS_Flags(_USESRTDATES, 0, 1)) {
    return 1 - valid_date(day(date), month(date), year(date));
  } else {

    Dsplit(date, &yy, &mm, &dd);
    if (yy < 0 || yy > 3000)
      return (1); // bad year
    if (mm < 1 || mm > 12)
      return (2); // bad month
    // see if leap year
    daysin = Isleap(yy) ? leap : noleap;
    if (dd < 1 || dd > daysin[mm])
      return (3); // bad day
    // all is well if we are here
    return (0);
  }
} // Dateok

void linterp(double x, double *y, double x0, double x1, double y0, double y1) {
  double a, b, l0, l1;

  a = x - x1;
  b = x0 - x1;
  l0 = a / b;

  a = x - x0;
  b = x1 - x0;
  l1 = a / b;

  *y = l0 * y0 + l1 * y1;

  return;
} /*linterp*/

double linear_interp(int n, double *xs, double *ys, double x,
                     int linear_extrapolate) {
  int j;
  j = 0;
  if (n == 0)
    return ys[0];
  while ((j < n) && (xs[j] < x))
    j++;
  if (j == 0) {
    if (linear_extrapolate) {
      j = 1;
    } else
      return ys[j];
  } else if (j == n) {
    if (linear_extrapolate) {
      j = n - 1;
    } else
      return ys[j - 1];
  } else {
    // do nothing
  }
  return ys[j - 1] +
         (ys[j] - ys[j - 1]) / (xs[j] - xs[j - 1]) * (x - xs[j - 1]);
} // linear_interp

double linterpFlat(int n, double *xs, double *ys, double x) {
  return (linear_interp(n, xs, ys, x, 0));
} // linterpFlat

// if x is in [x(i)        , x(i+1))        , return i+1
// if x < X(0)        , return 1
// if x >= x(n-1)        , return n-1
// if n = 0        , return 0
int interp_search(int n, double *xs, double x) {
  int j;
  j = 0;
  if (n == 0)
    return 0;
  while ((j < n) && (xs[j] < x))
    j++;
  if (j == 0) {
    return 1;
  } else if (j == n) {
    return (n - 1);
  } else {
    return (j);
  }
} // interp_search

double loglinear_interp(int n, double *xs, double *ys, double x,
                        int loglinear_extrapolate) {
  int j;
  j = 0;
  if (n == 0)
    return ys[0];
  while ((j < n) && (xs[j] < x))
    j++;
  if (j == 0) {
    if (loglinear_extrapolate) {
      j = 1;
    } else
      return ys[j];
  } else if (j == n) {
    if (loglinear_extrapolate) {
      j = n - 1;
    } else
      return ys[j - 1];
  } else {
    // do nothing
  }
  return ys[j - 1] *
         pow(ys[j] / ys[j - 1], (x - xs[j - 1]) / (xs[j] - xs[j - 1]));
} // loglinear_interp

double loglinterpFlat(int n, double *xs, double *ys, double x) {
  return (loglinear_interp(n, xs, ys, x, 0));
} // loglinterpFlat

//
// INPUTS to QINTERP
//
// price_n is a positive integer at least 4        ,
// price_x is an array of price_n variable values ("x values")        ,
// price is an array of price_n prices at x values given by price_x        ,
// dens_nb is a positive integer at least 2        ,
// dens_x is an array of dens_nb (generally different) x values        ,
// dens is an array of dens_nb densities for the x values given by dens_x.
//
// WARNING: price_x and dens_x must each be in numerical order (from lowest to
// highest).
//          Further        , the interpolation is meant to be applied only when
//          price_x[1] <= dens_x[1] and price_x[price_n-2] >= dens_x[dens_nb-1]
//                , i.e.        , there are at least two price_x values on each
//                side of the
//          entire distribution of dens_x values.
// WARNING: Further        , the dens arrays must always be at least as fine as
// the price arrays:
//          Between two adjacent dens_x values there should be at most one
//          price_x value.
//
// OUTPUT of QINTERPI
//
// if price_n = N        , price_x is x[0..(N-1)] price y[0..(N-1)]
// dens_nb = M        , dens_x z[0..(M-1)] dens d[0..(M-1)]
// for each z_i        , associates with x_j        , where z_i in (x_(j-1) ,
// x_j] or j=0 or j=(N-2) get the corresponding result of interpolating z_i
// using qinterps with x[(j-1)        ,(j+1)] (and let x_(-1)=x_0 x_N=x_(N-1))
// say the result is t_i
//       , then return sum_i z_i * t_i use: say we have density function r(x)
//       with
// values r[1..m] on [z[1..m]        , with m >> 1
// and f(x) with f[1..n] on x[1..n]        , with n moderate
// we seed int_(z_1)^(z_m) f(x)r(x) dx as follows:
// get F(x) with F[1..m] at z[1..m]        , with F[i] obtained via qintreps
// then we use sum r_i F_i as approximation.

double qinterpp(int price_n, double *price_x, double *price, int dens_nb,
                double *dens_x, double *dens) {
  int dens2pri[1000]; // Pointer from density array to price array.
  double tmpx[3], tmpp[3];
  double totval = 0.;
  int i, j, jcurr = 0;
  //
  // Fill the pointer from density array to price array.
  // i-th value gives index of price array whose value is just below that of
  // i-th in the density array.

  for (i = 0; i < dens_nb; i++) {
    while (price_x[jcurr] < dens_x[i])
      jcurr++;
    dens2pri[i] = jcurr;
    if (jcurr == 0)
      dens2pri[i] =
          1; // if price_x is less than smallest dens_x        , use the first
    if (jcurr >= (price_n - 1))
      dens2pri[i] =
          price_n - 2; // p_x is too large        , use the last 3 points
  }
  //
  // For each segment of density array        , add its contribution to the
  // average price.

  for (i = 0; i < dens_nb; i++) {
    if (dens2pri[i] == 1)
      totval +=
          dens[i] *
          (price[0] +
           (dens_x[i] - price_x[0]) * (price[1] - price[0]) /
               (price_x[1] -
                price_x[0])); // linear interpolation for the first interval
    else {
      jcurr = dens2pri[i];
      for (j = 0; j < 3; j++) {
        tmpx[j] = price_x[jcurr - 1 + j];
        tmpp[j] = price[jcurr - 1 + j];
      }
      totval += qinterpi(tmpx, tmpp, &dens_x[i], &dens[i]);
    }
  } // quadratic interpolation for other intervals

  return totval;
} // qinterpp

//
// INPUTS to QINTERPI
//
// price_x is an array of 4 variable values ("x values")        , only the first
// 3 are used price is an array of 4 prices at x values given by price_x , only
// the first 3 are used dens_x is an array of 2 (generally different) x values ,
// only the first is used dens is an array of 2 densities for the x values given
// by dens_x        , only the first is used
//
// WARNING: price_x and dens_x must each be in numerical order (from lowest to
// highest).
//          Further        , the interpolation is meant to be applied only when
//          the dens_x values are both between price_x[1] and price_x[2].
//
// OUTPUT of QINTERPI
// get the result y of interpolating dens_x[0] using qinterps on price_x[0..2]
// and price[0..2]
// and return y * dens[0]

double qinterpi(double *price_x, double *price, double *dens_x, double *dens) {
  double avgpri = 0.;
  double coeff[6];
  double xjunk;

  //
  // Quadratic interpolations: the first three values        , then the last
  // three values (out of four).

  xjunk = qinterps(price_x, price, *dens_x, coeff);
  avgpri = *dens * xjunk;

  return avgpri;
} // qinterpi

// INPUTS:
// x is an array of 3 variable values        ,
// y is an array of 3 prices at values given by x        ,
// x1 is a particular variable value.
//
// Returns the price at x1 based on quadratic interpolation:
// if x[1]==x[0]        , or x[2]==x[1]: linear intrep using x[0] and x[2]
// else        , use the Lagrange polynomial interpolation formula to get the
// quadratic passing through points (x{i]        , y[i])        , i=0        ,1
// ,2 amd read off the value at x

double qinterps(double *x, double *y, double x1, double *coeff) {
  double result;
  double intmed1 = (y[2] - y[0]) / (x[2] - x[0]);
  double intmed2;

  result = y[0] + intmed1 * (x1 - x[0]); // Linear between endpoints.
  coeff[0] = y[0] - x[0] * intmed1;
  coeff[1] = intmed1;
  coeff[2] = 0.;
  //
  // If not degenerate        , make the adjustment for full quadratic
  // interpolation.

  if ((x[0] != x[1]) && (x[1] != x[2])) {
    intmed2 = (y[1] - (y[0] + intmed1 * (x[1] - x[0]))) /
              ((x[1] - x[0]) * (x[1] - x[2]));
    result += intmed2 * (x1 - x[0]) * (x1 - x[2]);
    coeff[0] += intmed2 * x[0] * x[2];
    coeff[1] -= (x[0] + x[2]) * intmed2;
    coeff[2] = intmed2;
  }
  return result;
} // qinterps

FILE *_OpenMBSFile(char *data_dir, TBA_DATFile datfile) {
  static FILE *apFiles[11];
  static char *fileNames[] = {
      "Prepay.dat",   "Inhist.dat",     "Deal.dat",       "Term.dat",
      "Inweight.dat", "Holiday.dat",    "PrepayFunc.dat", "PPSpeedDens.dat",
      "BKTree.dat",   "TermStruct.dat", "PTDeals.dat"};

  FILE *pFile = 0;
  char string[MAXSTRBUFSZ];

  if (datfile < 0 || datfile >= sizeof(apFiles) / sizeof(apFiles[0]))
    return 0;
  pFile = apFiles[datfile];

  // close open file?
  if (pFile)
    fclose(pFile);
  pFile = 0;

  string[0] = '\0';
  strcat(string, data_dir);
  strcat(string, fileNames[datfile]);
  pFile = fopen(string, "r");
  apFiles[datfile] = pFile;

  return pFile;
}

int round(double x) {
  if (x - floor(x) < 0.5)
    return ((int)x);
  else
    return ((int)ceil(x));
} // round

// set vec[begin..end] = x
void set_vec_to_val(double *vec, double x, int begin, int end) {
  int i = begin;
  while (i++ <= end)
    *(vec++) = x;
} // set_vec_to_val