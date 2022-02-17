/* ===============================================================
   FILENAME : utdates.h

   PURPOSE:   all the necessary date functions
   =============================================================== */

#ifndef UTDATES_H
#define UTDATES_H

/* Include files */
#include "utError.h"
#include "utTypes.h"
#include "utString.h"

/* ----------------------------------------------------------
   The DATE type is a LONG									   
   ---------------------------------------------------------- */

typedef    long      Date;
typedef    double    Ddate;

/* ---------------------------------------------------------- */

Err		srt_f_gen_test_date(long theDate);

Date	add_unit(Date date1, int nb_unit, SrtUnit unit, SrtBusDayConv conv);

int     day_count_date(Date date1, Date date2, SrtBasisCode basis_code) ;

double  coverage(Date date1, Date date2, SrtBasisCode basis_code) ;  
double  coverage_ACT_ACT(Date date1, Date date2, SrtBasisCode basis_code);
double  num_days_year(Date date);
        
int     leap_year(int year);

int     weekday(int day, int month, int year);

int     valid_date(int day, int month, int year);

Date    bus_date_method(Date date1, SrtBusDayConv conv);

Date	str_to_date(char *dt);

void	date_to_str(Date dt, char *str);

Date    date(int d,int m,int y);

int     year(Date  n);
int     month(Date  n);
int     day(Date  n);

Date    third_wednesday(int m, int y);
Date    second_friday(int m, int y);
Date    eleventh_thursday(int m, int y);

Date    add_days_hlds(Date d, int nbd, double hlds[30][2], int nrow, int ncol, 
                      int sun);

Date    less_days_hlds(Date d, int nbd, double hlds[30][2], int nrow, int ncol, 
                       int sun);
Date    bus_hld_prv_fct(Date d, double hlds[30][2], int nrow, int ncol, 
                        int sun);
Date    bus_hld_fct(Date d, double hlds[30][2], int nrow, int ncol, int sun, 
                    int end_m);

Err		interp_tenor_string(String tenor,int *ny, int *nm, int *nd);

Err		add_tenor(Date start, String tenor, SrtBusDayConv conv, Date *end);

#define week_day(x)  ((x)%7)

#endif
