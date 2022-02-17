#ifndef MBSPTUTIL_H
#define MBSPTUTIL_H

#pragma warning(1 : 4100 4101) // NB: force warnings for unused arguments and
                               // local parameters
#pragma warning(1 : 4700 4701) // NB: force warnings for uninitialised variables

#include "math.h"
#include "memory.h"
#include "stdio.h"

///************************************************
/// MACROS
///************************************************

//**compilation controls***

//#define TIEOUT

#ifdef TIEOUT

#define TREEREFDATEISSETTLEDATE 1

#else

#define MBSPTCORRECTDAYCT 1
#define TREEREFDATEISSETTLEDATE 1

#endif // TIEOUT

// if we want to check that memory allocation/de-alloc's input limits are
// sensible

#ifndef CHECK_MEMALLOC_LIM
#define CHECK_MEMALLOC_LIM 1
#endif // CHECK_MEMALLOC_LIM

//**MAX SIZES**

#define MAXSTRBUFSZ 128
#define MAXHOLIDAYS 5000

//***Constants***
#define TRUE 1
#define FALSE 0
#define TBAOUTPUTDIR "O:/Options/Tba/outfiles/"
#define TBADATADIR "O:/Options/Tba/datfiles/"
#define MBS_SMALL (1E-7)

//*MACROS******

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MBSUNUSED
#define MBSUNUSED(x) x
#endif

#define MBS_INIT_ZERO(__var__) memset(&__var__, 0, sizeof(__var__))
#define _MBS_ALL_ZEROS                                                         \
  { 0 }
///************************************************
/// Program running mode
///************************************************

#define MBS_WRITE_FILE 1

// Control Flags

typedef enum MBS_FlagNames { _USESRTDATES = 0 } MBS_FlagNames;

int MBS_Flags(MBS_FlagNames FlagNames, int value, int readOnly);

///************************************************
/// Error Handling
///************************************************

char *mbs_err_handler(char *error_text, int write, int exit);

char *nrerror(char *error_text);
///************************************************
/// Resources management
//*************************************************

//***basic memory***

long *mbspt_lvector(int nl, int nh);

int *mbspt_ivector(int nl, int nh);

double *mbspt_dvector(int nl, int nh);

double **mbspt_dmatrix(int nrl, int nrh, int ncl, int nch);

void mbspt_free_lvector(long *v, int nl, int nh);

void mbspt_free_ivector(int *v, int nl, int nh);

void mbspt_free_dvector(double *v, int nl, int nh);

void mbspt_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

enum DATE_UNIT *mbspt_DateUnitVector(int start, int end);
void mbspt_free_DateUnitVector(enum DATE_UNIT *v, int start, int end);
///*********************************************
/// Date Stuffs
///*********************************************

// BASIS and UNITS

typedef enum YBASIS {
  _BADYBASIS = -1,
  _YACT = 0,  //"ACT"
  _360 = 360, //"360"
  _365 = 365, //"365"
} YBASIS;

typedef enum DIFFBASIS {
  _BADDIFFBASIS = -1,
  _DACT = 0, //"ACT"
  _30,       //"30"
} DIFFBASIS;

typedef enum DATE_UNIT { DAY = 0, MONTH, YEAR } DATE_UNIT;

// code must have '/' like "ACT/ACT"
char *readBasis(YBASIS *ybasis, DIFFBASIS *diffbasis, char *code);

YBASIS readYBASIS(char *code);

DIFFBASIS readDIFFBASIS(char *code);

// Date Functions
void Dsplit(long date_i, long *yy_o, long *dmm, long *dd_o);
int month_diff(long YMMDD1, long YMMDD2);

int Isleap(long year);
long Dayofwk(long date);

long YY2YYYY(long Y);

long DateYMMDD(long Y, long MM, long DD);

long Nxtday(long date, long days);
long Nxtwkday(long date, long advance);
long Nxtmth(long date, long mths, long eom);

typedef struct {
  int numHolidays;
  long Holidays[MAXHOLIDAYS];
} Holidays;

char *HolidaysConstruct(Holidays *holidays, char *dir);
int isHoliday(long date_i);

int isBusDay(long date_i);

long MBSNxtbusday(long date, long advance);

long Date_Add(long Ref_Date, DATE_UNIT unit, int len, int bus, int eom);

long YYYYMMDD2XlDate(long d);
long fromYYYYMMDD(long d);

long YYYYMM2FirstOfMon(long d);

long DateDiff(long date1_i, long date2_i, DIFFBASIS base);
int Dateok(long date);

// interpolations
void linterp(double x, double *y, double x0, double x1, double y0, double y1);

// to be used in conjunction with the reurn value of interp_search  , if the
// return value is j  , x0=xs[j-1]
#define LogInterp(x, x0, x1, y0, y1)                                           \
  (y0) * pow(((y1) / (y0)), (((x) - (x0)) / ((x1) - (x0))));

double
linear_interp(int n, double *xs, double *ys, double x,
              int linear_extrapolate); // interpolate or do flat extrapolation

double linterpFlat(int n, double *xs, double *ys,
                   double x); // interpolate or do flat extrapolation

int interp_search(int n, double *xs, double x);

double loglinear_interp(
    int n, double *xs, double *ys, double x,
    int loglinear_extrapolate); // interpolate or do flat extrapolation

double loglinterpFlat(int n, double *xs, double *ys,
                      double x); // interpolate or do flat extrapolation

double qinterpp(int price_n, double *price_x, double *price, int dens_nb,
                double *dens_x, double *dens);

double qinterpi(double *price_x, double *price, double *dens_x, double *dens);

double qinterps(double *x, double *y, double x1, double *coeff);

// FILE Handling

typedef enum TBA_DATFILE {
  FilePrepay = 0,
  FileInhist, // stay
  FileDeal,
  FileTerm,
  FileInweight,
  FileHolidays, // stay
  FilePPFUNCS,
  FilePPSpeedDens,
  FileBKTree,
  FileTS,
  FilePTDeals
} TBA_DATFile;

FILE *_OpenMBSFile(char *data_dir, TBA_DATFile datfile);
// FILE * MBS_Sequential_File_Manager( char * data_dir  , char * fileName );

// misc
int round(double x);
void set_vec_to_val(double *vec, double x, int begin, int end);

#endif // MBSPTUTIL_H