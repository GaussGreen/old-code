#ifndef MBSPTPPFUNCS_H
#define MBSPTPPFUNCS_H

//**MAX****
#define NBPOINTS 9 /* Number of points used to define the prepayment function  \
                    */
#define MAXSFSNB 200 /* Maximum number of SFS structures that can be priced */
#define MAXFIXEDPREPAY 180 /* Maximum number of fixed prepayment inputs */
#define MAXNUMP 6
#define NUMRAMPS 3
#define MAXOTERM 360
#define MAXNBWEIGHTS 5000
//************************

//**CONSTANTS**
#define REFI_THRESHOLD 0.9
//***************

//**MACROS***
#define transformIncentive(incent, gwac, refi, CMS10Y)                         \
  (incent * gwac - refi + CMS10Y)
//**********

// densities and histories
typedef struct {

  long dates[1000];    // historical date
  double hrates[1000]; // rates at thoese dates
  int Number;          // number of dates or rates in the file
} HISTORY;

char *history_reader(char *data_dir, HISTORY *history);

void history_copy(HISTORY *in, HISTORY *out);

char *HISTORY_Check(HISTORY *history);

typedef struct {
  double irates[20]; // irates at which density is known initially
  double IncentKnots[20];
  int nirates;            // number of input rates at which the density is known
  double iweights[20];    // input weights
  double fastweights[20]; // fast density at irates
  double sloweights[20];  // slow density at irates
  int Nbweights;          // number of weights used in prepayment function
  double urates[MAXNBWEIGHTS];
  double ufweights[MAXNBWEIGHTS];
  double usweights[MAXNBWEIGHTS];
} MBSDENSITY;

char *density_reader(char *data_dir, MBSDENSITY *density);

void density_copy(MBSDENSITY *in, MBSDENSITY *out);

char *MBSDENSITY_Check(MBSDENSITY *density);

void adjust_density(MBSDENSITY *density, double M, double gwac, double refiRate,
                    double Rate10Y,
                    long Nbweights); // potentially also refiRate and Rate10Y
void adjust_history(HISTORY *history, int term, double refiRate,
                    double Rate10Y); // potentially also refiRate and Rate10Y

void spdclass(double *fastweights, double M, int NUMWTS);

void neweights(double *rates, double *weights, double *inrates,
               double *inweights, int N);

void ldensity(double *irates,   /* input rates  */
              double *iweights, /* input weights at the corresponding rates */
              double *rate,     /* rate of interest */
              double *weight);  /*  weight at rate of interest */

// 1st payment <-> payment_nb=1; wam = oterm  , calc the scheduled principal
// paydown during the payment_nb-th payment  , per $1 principal remaining before
// the payment
double SCHEDULED_AMORT_RATE(int wam, double gwac);

double TOTAL_AMORT_RATE(int wam, double gwac, double smm);

#endif // MBSPTPPFUNCS_H