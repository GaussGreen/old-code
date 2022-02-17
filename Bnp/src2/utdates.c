/*****************************************************************************
        FILE_NAME:	utdates.c
******************************************************************************/

#include "utallhdr.h"

/*****************************************************************************
        Check that  number input can be considered as a date
*****************************************************************************/
Err srt_f_gen_test_date(long theDate) {
  if (theDate < 30317) {
    /* theDate is before 01-jan-1983 */
    return (serror("date = %d too far into the past", theDate));
  } else if (theDate > 73051) {
    /* the date is after 01-jan-2100 */
    /* changed to allow 50Y deals */
    return (serror("date = %d too far into the future (01/01/2100)", theDate));
  } else {
    /* theDate is in-range */
    return (NULL);
  }
}

/*****************************************************************************
        Adds number of days  , months or years to date1
*****************************************************************************/
Date add_unit(Date date1, int nb_unit, SrtUnit unit, SrtBusDayConv conv) {
  int d, m, y;
  int dirn;
  int i, q, r;

  d = day(date1);
  m = month(date1);
  y = year(date1);

  if (unit == SRT_WEEK) {
    nb_unit *= 7;
  }

  switch (unit) {
  case SRT_DAY:
  case SRT_WEEK:
    date1 += nb_unit;

    if (conv == SUCCEEDING) {
      if (week_day(date1) == SATURDAY)
        return date1 + 2;
      if (week_day(date1) == SUNDAY)
        return date1 + 1;
    } else if (conv == MODIFIED_SUCCEEDING) {
      d = day(date1);
      m = month(date1);
      y = year(date1);

      if (week_day(date1) == SATURDAY) {
        d = d + 2;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      } else if (week_day(date1) == SUNDAY) {
        d = d + 1;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      }

      date1 = date(d, m, y);
    }

    break;

  case SRT_BDAY:
    dirn = (nb_unit >= 0 ? 1 : -1);

    q = (int)(abs(nb_unit) / 5);
    r = (int)(abs(nb_unit) % 5);
    date1 += 7 * q * dirn;

    /* we have to go to the next business day first */
    /* unless we just want the next bus day in modified succeeding (nb_unit ==
     * 0) */
    if ((nb_unit == 0) && (conv == MODIFIED_SUCCEEDING)) {
      d = day(date1);
      m = month(date1);
      y = year(date1);

      if (week_day(date1) == SATURDAY) {
        d = d + 2;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      } else if (week_day(date1) == SUNDAY) {
        d = d + 1;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      }

      date1 = date(d, m, y);
    }
    /* other 'succeeding' */
    else if (week_day(date1) == SATURDAY)
      date1 = date1 + 2;
    else if (week_day(date1) == SUNDAY)
      date1 = date1 + 1;

    i = 0;
    while (i < r) {
      date1 += 1 * dirn;
      if ((week_day(date1) != SATURDAY) && (week_day(date1) != SUNDAY))
        i++;
    }

    break;

  case SRT_MONTH:
    m += nb_unit;
    nb_unit = 0;
    /* the number of months went from 1 to 12 */
    while (m > 12) {
      m -= 12;
      y++;
    }

    while (m < 1) {
      m += 12;
      y--;
    }

  case SRT_YEAR:
    y += nb_unit;

    while (valid_date(d, m, y) == 0)
      d--;

    date1 = date(d, m, y);

    if (conv == SUCCEEDING) {
      if (week_day(date1) == SATURDAY)
        return date1 + 2;
      if (week_day(date1) == SUNDAY)
        return date1 + 1;
    } else if (conv == MODIFIED_SUCCEEDING) {
      if (week_day(date1) == SATURDAY) {
        d = d + 2;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      } else if (week_day(date1) == SUNDAY) {
        d = d + 1;
        if (valid_date(d, m, y) == 0)
          d = d - 3;
      }
    }

    date1 = date(d, m, y);

    break;

  default:
    return 0;
    break;
  }

  return date1;
}

/* ------------------------------------------------------------- */
/* Add a tenor to a start date */
/* ------------------------------------------------------------- */
Err add_tenor(Date start, String tenor, SrtBusDayConv conv, Date *end) {
  int ny = 0, nm = 0, nd = 0;
  Err err;

  strupper(tenor);
  strip_white_space(tenor);

  if (err = interp_tenor_string(tenor, &ny, &nm, &nd)) {
    return err;
  }

  start = add_unit(start, ny, SRT_YEAR, NO_BUSDAY_CONVENTION);
  start = add_unit(start, nm, SRT_MONTH, NO_BUSDAY_CONVENTION);
  *end = add_unit(start, nd, SRT_DAY, conv);

  return NULL;
}

/* ------------------------------------------------------------- */

/*
        interpret a tenor string into
        a number of years  , of months  , and of days.
        if days  , there cannot be months or years.
*/

#define ERRR serror("bad tenor %s", tenor)

Err interp_tenor_string(String tenor, int *ny, int *nm, int *nd) {
  String endp;
  int num1, num2, num3;
  char let1, let2, let3;
  int num_part, num_years = 0, num_months = 0, num_days = 0;

  /* first number and letter */

  num1 = (int)strtol(tenor, &endp, 10);
  if (!(*endp) || (num1 < 0))
    return ERRR;

  let1 = *endp;
  if ((let1 != 'M') && (let1 != 'Y') && (let1 != 'D') && (let1 != 'W'))
    return ERRR;

  endp++;
  if (!(*endp))
    num_part = 1;

  /* second number and letter */

  else {
    num2 = (int)strtol(endp, &endp, 10);
    if (!(*endp) || (num2 < 0))
      return ERRR;

    let2 = *endp;
    if (let1 == let2)
      return ERRR;

    endp++;
    if (!(*endp))
      num_part = 2;

    /* third number and letter */

    else {
      num_part = 3;

      num3 = (int)strtol(endp, &endp, 10);
      if (!(*endp) || (num3 < 0))
        return ERRR;

      if (num3 < 0)
        return ERRR;

      let3 = *endp;
      if ((let3 == let2) || (let3 == let1))
        return ERRR;

      endp++;
      if (*endp != '\0')
        return ERRR;
    }
  }

  switch (num_part) {
  case 1:
    if (let1 == 'Y')
      num_years = num1;
    if (let1 == 'M')
      num_months = num1;
    if (let1 == 'D')
      num_days = num1;
    if (let1 == 'W')
      num_days = num1 * 7;
    break;

  case 2:
    if (let1 == 'Y')
      num_years = num1;
    if (let1 == 'M')
      num_months = num1;
    if (let1 == 'D')
      num_days = num1;

    if (let2 == 'Y')
      num_years = num2;
    if (let2 == 'M')
      num_months = num2;
    if (let2 == 'D')
      num_days = num2;
    break;

  case 3:
    if (let1 == 'Y')
      num_years = num1;
    if (let1 == 'M')
      num_months = num1;
    if (let1 == 'D')
      num_days = num1;

    if (let2 == 'Y')
      num_years = num2;
    if (let2 == 'M')
      num_months = num2;
    if (let2 == 'D')
      num_days = num2;

    if (let3 == 'Y')
      num_years = num3;
    if (let3 == 'M')
      num_months = num3;
    if (let3 == 'D')
      num_days = num3;
    break;

  default:
    return ERRR;
  }

  *ny = num_years;
  *nm = num_months;
  *nd = num_days;

  return NULL;
}

#undef ERRR

/*****************************************************************************
        Returns number of days between date1 and date2  ,
        under the conventions implied by basis_code and currency

        date1 should be before date2
*****************************************************************************/
int day_count_date(Date date1, Date date2, SrtBasisCode basis_code) {
  int d1 = day(date1), d2 = day(date2);
  int m1 = month(date1), m2 = month(date2);
  int y1 = year(date1), y2 = year(date2);
  int num_days = 0;

  switch (basis_code) {
  case BASIS_30_360:
    if (d1 == 31)
      d1 = 30;
    if ((d1 == 30) && (d2 == 31))
      d2 = 30;
    num_days = d2 - d1 + 30 * (m2 - m1) + 360 * (y2 - y1);
    break;

  case BASIS_30_360E:
    if (d1 == 31)
      d1 = 30;
    if (d2 == 31)
      d2 = 30;
    num_days = d2 - d1 + 30 * (m2 - m1) + 360 * (y2 - y1);
    break;

  case BASIS_ACT_360:
  case BASIS_ACT_365:
  case BASIS_ACT_ACT:
  case BASIS_ACT_USD:
  default:
    num_days = date(d2, m2, y2) - date(d1, m1, y1);
    break;
  }

  return (num_days);
}

/*****************************************************************************
    Returns the coverage  , i.e num_days/ basis

    Assumes that date2 >= date1 and year(date2)<= year(date1)+1
        This function is wrong for BASIS_ACT_USD but there is no way
        we can do otherwise
*****************************************************************************/

double coverage(Date date1, Date date2, SrtBasisCode basis_code) {
  int d1 = day(date1), d2 = day(date2);
  int m1 = month(date1), m2 = month(date2);
  int y1 = year(date1), y2 = year(date2);
  double res;

  switch (basis_code) {
  case BASIS_ACT_ACT:
    res = coverage_ACT_ACT(date1, date2, basis_code);
    break;
  case BASIS_ACT_360:
  case BASIS_30_360:
  case BASIS_30_360E:
    res = day_count_date(date1, date2, basis_code) / 360.0;
    break;
  case BASIS_ACT_365:
  case BASIS_ACT_USD:
  default:
    res = day_count_date(date1, date2, basis_code) / 365.0;
    break;
  }

  return (res);
}

/*****************************************************************************
    Returns the coverage for the case of ACT_ACT  , which is more complex than
        the other basis types.

    Assumes that year(date1) <= year(date2) <=  year(date1)+1

   If(year(date1) = year(date2)
             cover = (date2 - date1) / num_days_year(date1)
   Else  if year(date2) < year(date1)+1
       cover = num_of_days_in_date1_year / num_days_year(date1)
                       +  num_of_days_in_date1_year / num_days_year(date1)

*****************************************************************************/
double coverage_ACT_ACT(Date date1, Date date2, SrtBasisCode basis_code) {
  int d1 = day(date1), d2 = day(date2);
  int m1 = month(date1), m2 = month(date2);
  int y1 = year(date1), y2 = year(date2);
  int d_jan1 = 1, m_jan1 = 1;
  Date begin_year;
  double res = 0;

  if (y1 != y2) {
    begin_year = date(d_jan1, m_jan1, y1 + 1);
    res = ((double)day_count_date(date1, begin_year, basis_code)) /
          num_days_year(date1);

    y1++;
    while (y1 != y2) {
      res += 1;
      y1++;
    }

    begin_year =
        date(d_jan1, m_jan1, y2); // Stanley Mrose (11.11.2002): have to update
                                  // start date when calculating final stub
    res += ((double)day_count_date(begin_year, date2, basis_code)) /
           num_days_year(date2);
  } else
    res = ((double)day_count_date(date1, date2, basis_code)) /
          num_days_year(date1);

  return (res);
}

/**************************************************************/
/* GET THE 'CURRENT BUSINESS DAY' DEPENDING ON THE CONVENTION */
/**************************************************************/
Date bus_date_method(Date date1, SrtBusDayConv conv) {
  date1 = add_unit(date1, 0, SRT_DAY, conv);

  return (date1);
}

/***********/
/* WEEKDAY */
/***********/
int weekday(int day, int month, int year) {
  Date date1;

  date1 = date(day, month, year);

  return week_day(date1);
}

/*************************/
/* LEAP YEAR CALCULATION */
/*************************/
int leap_year(int year) {
  year += 1900;

  if (((year % 4 == 0) && (year % 100 != 0)) || (year % 400 == 0)) {
    return 1;
  } else
    return 0;
}

/***************************************/
/* RETURN THE NUMBER OF DAYS IN A YEAR */
/***************************************/
double num_days_year(Date date) {
  if (leap_year(year(date)) == 1)
    return (366.0);
  else
    return (365.0);
}

/**************/
/* VALID DATE */
/**************/
int valid_date(int day, int month, int year) {
  static int days[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if (month > 12 || month < 1 || day > 31 || day < 1)
    return 0;

  if (day > days[month]) {
    if (month == 2 && day == 29) {
      if (leap_year(year)) {
        return 1;
      }
    }

    return 0;
  }

  return 1;
}

/*********************************/
/* DATE FUNCTIONS AS LOTUS DATES */
/*********************************/

/* Transform a date in string into a number */
Date str_to_date(char *dt) {
  int y, m, d;

  d = atoi(dt);
  m = atoi(dt + 3);
  y = atoi(dt + 6);

  return date(d, m, y);
}

/* Transform a date in number into a string dd/mm/yyyy */
void date_to_str(Date dt, char *str) {
  int d, m, y;

  d = day(dt);
  m = month(dt);
  y = year(dt) + 1900;

  sprintf(str, "%02d/%02d/%04d", d, m, y);
}

/*********************************/
/* DATE FUNCTIONS AS LOTUS DATES */
/*********************************/

/* TRANSFORM IN DATE */
Date date(int d, int m, int y) {
  static int cummonth[] = {0,   0,   31,  59,  90,  120, 151,
                           181, 212, 243, 273, 304, 334, 365};
  static int mlength[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  if (y > 1899)
    y -= 1900;

  if (d < 1 || d > 31 || m < 1 || m > 12 || y < 0 || y > 199)
    return -1;

  if ((m != 2 && d > mlength[m]) || (m == 2 && d > 28 + (y % 4 ? 0 : 1)))
    return -1;

  return ((long)y / 4) * 1461 + (y % 4) * 365 + cummonth[m] + d +
         (y % 4 == 0 && m < 3 ? 0 : 1);
}

/* YEAR */
int year(Date n) {
  int y;

  if (n < 1 || n > 73050)
    return -1;

  y = (int)(((n - 1) / 1461) * 4);
  y = y + (int)((n - (long)y / 4 * 1461 - 2) / 365);

  if (y < 0 || y > 199)
    return -1;

  return y;
}

/* MONTH */
int month(Date n) {
  static int cummonth[] = {0,   0,   31,  59,  90,  120, 151,
                           181, 212, 243, 273, 304, 334, 365};
  static int mlength[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  int m = 1;

  if (n < 1 || n > 73050)
    return -1;

  n = (n - 1) % 1461;

  if (n < 31)
    return 1;
  if (n < 60)
    return 2;

  n = (n - 1) % 365 + 1;

  while (n > cummonth[m + 1])
    m++;

  return m;
}

/* DAY */
int day(Date n) {
  if (n < 0 || n > 73050)
    return -1;

  return (int)(n - date(1, month(n), year(n)) + 1);
}

/***************************************
        return third wednesday of a given month and year
***************************************/
Date third_wednesday(int m, int y) {
  int num_wed;
  Date d1;

  /* first day of month */
  d1 = date(1, m, y);

  num_wed = 2 + (week_day(d1) == WEDNESDAY ? 1 : 0);
  d1 += 14;

  if (num_wed == 2)
    while (week_day(d1) != WEDNESDAY)
      d1++;

  return d1;
}

/***************************************
        return second_friday of a given month and year
***************************************/
Date second_friday(int m, int y) {
  Date d1;

  /* first day of month */
  d1 = date(1, m, y);

  while (week_day(d1) != FRIDAY)
    d1++;

  d1 += 7;

  return d1;
}

/***************************************
        return eleventh thursday of a given quarter
        (1st mo of quarter is given) and year
***************************************/
Date eleventh_thursday(int m, int y) {
  int num_thu;
  Date d1;

  /* first day of month */
  d1 = date(1, m, y);

  num_thu = 10 + (week_day(d1) == THURSDAY ? 1 : 0);
  d1 += 70;

  if (num_thu == 10)
    while (week_day(d1) != THURSDAY)
      d1++;

  return d1;
}

/*********************************************************************
     Amended by julia matsumoto
*/

Date bus_hld_prv_fct(Date d, double hlds[30][2], int nrow, int ncol, int sun) {
  int j, flag;

  flag = 1;

  while (flag != 0) {
    flag = 0;

    for (j = 0; j < nrow; j++) {
      if ((week_day(d) != SATURDAY) && (week_day(d) != SUNDAY) &&
          (month(d) == (int)hlds[j][0]) && (day(d) == (int)hlds[j][1])) {
        d--;
        flag++;
      }

      if ((week_day(d - 1) == SUNDAY) && (month(d - 1) == (int)hlds[j][0]) &&
          (day(d - 1) == (int)hlds[j][1]) && (sun == 1)) {
        d -= 3;
        flag++;
      }
    }

    if (week_day(d) == SATURDAY) {
      d--;
      flag++;
    }

    if (week_day(d) == SUNDAY) {
      d -= 2;
      flag++;
    }
  }

  return d;
}

Date bus_hld_fct(Date d, double hlds[30][2], int nrow, int ncol, int sun,
                 int end_m) {
  int j, flag;
  Date dold;

  dold = d;
  flag = 1;

  while (flag != 0) {
    flag = 0;

    for (j = 0; j < nrow; j++) {
      if ((week_day(d) != SATURDAY) && (week_day(d) != SUNDAY) &&
          (month(d) == (int)hlds[j][0]) && (day(d) == (int)hlds[j][1])) {
        d++;
        flag++;
      }

      if ((week_day(d) == SUNDAY) && (month(d) == (int)hlds[j][0]) &&
          (day(d) == (int)hlds[j][1]) && (sun == 1)) {
        d += 2;
        flag++;
      }

      if ((week_day(d) == SATURDAY) && (month(d + 1) == (int)hlds[j][0]) &&
          (day(d + 1) == (int)hlds[j][1]) && (sun == 1)) {
        d += 3;
        flag++;
      }
    }

    if (week_day(d) == SATURDAY) {
      d += 2;
      flag = 1;
    }
    if (week_day(d) == SUNDAY) {
      d++;
      flag = 1;
    }
  }

  if ((month(d) != month(dold)) && (end_m)) {
    d = bus_hld_prv_fct(dold, hlds, nrow, ncol, sun);
  }

  return d;
}

/*****************************************************************************
        Adds number of business days (excluding holidays) to date1
*****************************************************************************/
Date add_days_hlds(Date d, int nbd, double hlds[30][2], int nrow, int ncol,
                   int sun) {
  int i;
  Date dn;

  dn = bus_hld_fct(d, hlds, nrow, ncol, sun, 0);
  if (dn != d) {
    d = dn;
    nbd--;
  }

  for (i = 0; i < nbd; i++) {
    d++;
    d = bus_hld_fct(d, hlds, nrow, ncol, sun, 0);
  }

  return (d);
}

/*****************************************************************************
        Adds number of business days (excluding holidays) to date1
*****************************************************************************/
Date less_days_hlds(Date d, int nbd, double hlds[30][2], int nrow, int ncol,
                    int sun) {
  int i;
  Date dn;

  dn = bus_hld_prv_fct(d, hlds, nrow, ncol, sun);

  if (dn != d) {
    d = dn;
    nbd--;
  }

  for (i = 0; i < nbd; i++) {
    d--;
    d = bus_hld_prv_fct(d, hlds, nrow, ncol, sun);
  }

  return (d);
}
/*
     end of amended by julia matsumoto
*/
