#include "math.h"
#include "stdio.h"
#include "string.h"

#include "MBSPPFuncs.h"
#include "MBSPTUtil.h"

char *history_reader(char *data_dir, HISTORY *history) {
  int i, rcode;
  double l[20] = {0.0}, sum = 0.0, r;
  long d;
  char string[128];
  FILE *stream;

  stream = _OpenMBSFile(data_dir, FileInhist);
  if (!stream)
    return (nrerror("Can't open Inhist.dat"));
  i = 0;
  while ((rcode = fscanf(stream, "%ld %lf\n", &d, &r)) != EOF) {
    if (rcode != 2) {
      sprintf(string, "Can't read history.dates[%d] or hrates[%d]", i, i);
      return (nrerror(string));
    }

    history->dates[i] = d;
    history->hrates[i] = r;
    i++;
  }
  history->Number = i;

  return (0);
} // History_Reader

void history_copy(HISTORY *in, HISTORY *out) {
  int i;
  if (in && out) {
    out->Number = in->Number;
    for (i = 0; i < in->Number; ++i) {
      out->dates[i] = in->dates[i];
      out->hrates[i] = in->hrates[i];
    }
  }
} // history_copy

char *HISTORY_Check(HISTORY *history) {
  int i;
  long prev = 0, curr;
  char string[MAXSTRBUFSZ];
  if (!history)
    return (0);
  if (history->Number > 0) {
    prev = history->dates[0] * 100 + 1;
    prev = fromYYYYMMDD(prev);

    if (Dateok(prev) != 0)
      return (nrerror("History initial date wrong format"));
  }
  prev = Nxtmth(prev, -1, 1);
  for (i = 0; i < history->Number; ++i) {
    if (history->hrates[i] < 0.0)
      return (nrerror("Some historical rates < 0"));
    curr = history->dates[i] * 100 + 1;
    curr = fromYYYYMMDD(curr);

    if (curr != Nxtmth(prev, 1, 1)) {
      string[0] = '\0';
      sprintf(string, "%ld", prev);
      return (nrerror(
          strcat("History dates are not consecutive at        , ", string)));
    }
    prev = curr;
  }
  return (0);
} // HISTORY_Check

char *density_reader(char *data_dir, MBSDENSITY *density) {
  int i, rcode;
  double l[20] = {0.0}, tmp, sum = 0.0;
  char string[128];
  FILE *stream;

  stream = _OpenMBSFile(data_dir, FileInweight);
  if (!stream)
    return (nrerror("Can't open Inweight.dat"));

  fgets(string, 80, stream);
  i = 0;
  while ((rcode = fscanf(stream, "%lf %lf\n", &(density->IncentKnots[i]),
                         &(density->iweights[i]))) != EOF) {
    if (rcode != 2) {
      sprintf(string, "Can't read density.IncentKnots[%d] or iweights[%d]", i,
              i);
      return (nrerror(string));
    }
    sum += density->iweights[i];
    i++;
  }
  density->nirates = i;
  for (i = 0; i < density->nirates / 2; ++i) {
    tmp = density->iweights[density->nirates - 1 - i];
    density->iweights[density->nirates - 1 - i] = density->iweights[i];
    density->iweights[i] = tmp;
    tmp = density->IncentKnots[density->nirates - 1 - i];
    density->IncentKnots[density->nirates - 1 - i] = density->IncentKnots[i];
    density->IncentKnots[i] = tmp;
  }
  for (i = 0; i < density->nirates; i++) {
    density->irates[i] = density->IncentKnots[i];
    density->fastweights[i] = density->iweights[i] / sum;
    density->sloweights[i] = density->fastweights[i];
  }
  return (0);
} // density_reader
void density_copy(MBSDENSITY *in, MBSDENSITY *out) {
  int i;
  if (in && out) {
    out->nirates = in->nirates;
    for (i = 0; i < in->nirates; ++i) {
      out->IncentKnots[i] = in->IncentKnots[i];
      out->irates[i] = in->irates[i];
      out->iweights[i] = in->iweights[i];
      out->fastweights[i] = in->fastweights[i];
      out->sloweights[i] = in->sloweights[i];
    }
    out->Nbweights = in->Nbweights;
    for (i = 0; i < in->Nbweights; ++i) {
      out->urates[i] = in->urates[i];
      out->ufweights[i] = in->ufweights[i];
      out->usweights[i] = in->usweights[i];
    }
  }
} // density_copy

char *MBSDENSITY_Check(MBSDENSITY *density) {
  double prev = 0.0;
  int i;

  if (!density)
    return (0);
  if (density->nirates < 1)
    return (nrerror("Invalid density.nirates"));
  if (density->Nbweights < 1)
    return (nrerror("Invalid Nbweights"));
  for (i = 0; i < density->nirates; ++i) {
    if (i > 0 && density->irates[i] <= prev)
      return (nrerror("density.irates must be increasing"));
    prev = density->irates[i];
    if (density->fastweights[i] < 0.0 || density->sloweights[i] < 0.0)
      return (nrerror("Fast/Slow weights >=0!"));
  }
  return (0);
} // MBSDENSITY_Check

void adjust_density(MBSDENSITY *density, double M, double gwac, double refiRate,
                    double Rate10Y, long Nbweights) {
  int i;
  spdclass(density->fastweights, M, density->nirates);
  for (i = 0; i < density->nirates; i++) {
    density->sloweights[i] = density->fastweights[i];
    density->irates[i] =
        transformIncentive(density->IncentKnots[i], gwac, refiRate,
                           Rate10Y); // convert ratio into absolute 10Y rate
  }
  density->Nbweights = Nbweights;
  neweights(density->urates, density->ufweights, density->IncentKnots,
            density->fastweights,
            density->Nbweights); /* linear interpolation of density */
  neweights(density->urates, density->usweights, density->IncentKnots,
            density->sloweights, density->Nbweights);
  for (i = 0; i < density->Nbweights; ++i)
    density->urates[i] = transformIncentive(
        density->urates[i], gwac, refiRate,
        Rate10Y); // transform urates from mtgRate space to 10YRate space
} // adjust_density

void adjust_history(HISTORY *history, int term, double refiRate,
                    double Rate10Y) {
  int i;
  double agency_adj = 0.0;
  // FIX: this has to be reconciled with speedup adjustment in
  // prepay_engine_setparams
  if (term == 180)
    agency_adj = 0.5;

  for (i = 0; i < history->Number; ++i) {
    history->hrates[i] = transformIncentive(history->hrates[i], 1.0, refiRate,
                                            Rate10Y); // convert to 10Y rate
    history->hrates[i] -= agency_adj;
  }
} // adjust_history

//***spdclass***
//*Tweak the density function of refinancing incentives
//
void spdclass(double *fastweights, double M, int NUMWTS) {
  int i;

  double l[2000], fract[2000];

  if (M != 1.0) {
    if (M < 1E-14)
      M = 1E-14;
    l[0] = fastweights[0];
    for (i = 1; i < NUMWTS; i++)
      l[i] =
          l[i - 1] + fastweights[i]; /* Fractions "left" (not triggered yet). */
    fract[0] = 1.0;
    for (i = 1; i < NUMWTS; i++)
      fract[i] =
          fastweights[i] / l[i]; /* Fraction of remaining at this trigger. */
    for (i = 0; i < NUMWTS; i++)
      fract[i] = 1. - (pow(1. - fract[i],
                           M)); /* Speed up these triggers by a factor of M. */

    l[NUMWTS - 1] = 1.0; /* We will now keep track of updated fractions "left"
                            after the speed-up. Same value of 1.0 here. */
    for (i = NUMWTS - 1; i >= 0; i--) {
      fastweights[i] = fract[i] * l[i]; /* New weight at this trigger. */
      if (i > 0)
        l[i - 1] =
            l[i] - fastweights[i]; /* New fraction "left" for next trigger. */
    }
  }
} // spdclass

void neweights(double *rates, double *weights, double *inrates,
               double *inweights, int N) {
  int i;
  double sum = 0.;
  for (i = 0; i < N; i++) {
    rates[i] =
        0.92 * (double)i /
        ((double)N - 1.); /* generate uniformly distributed trigger rates */

    ldensity(inrates, inweights, &rates[i],
             &weights[i]); /* linear interpolation */
    sum += weights[i];
  }
  for (i = 0; i < N; i++) {
    weights[i] /= sum; /* normalization */
  }
} // neweight

void ldensity(double *irates,   /* input rates  */
              double *iweights, /* input weights at the corresponding rates */
              double *rate,     /* rate of interest */
              double *weight)   /*  weight at rate of interest */

{
  int i;

  i = 0;

  if (*rate >= irates[16])
    *weight = 0.;
  else {
    while (*rate >= irates[i])
      i += 1; /* find the right segment */
    *weight = iweights[i - 1] +
              (*rate - irates[i - 1]) * (iweights[i - 1] - iweights[i]) /
                  (irates[i - 1] - irates[i]); /* linear interpolation */
  }
} // ldensity

// 5% gwac (annualized)<-> gwac = 5
// wam is in month
// wam = o_term - payment_nb + 1
// where
// 1st payment <-> payment_nb=1; calc the scheduled principal paydown during the
// payment_nb-th payment        , per $1 principal remaining before the payment
double SCHEDULED_AMORT_RATE(int wam, double gwac) {
  double c = (1.0 + gwac / 1200.0);
  return ((c - 1.0) / (pow(c, (double)(wam)) - 1.0));
} // SCHEDULED_AMORT_RATE

// 1% smm <-> smm=1
double TOTAL_AMORT_RATE(int wam, double gwac, double smm) {
  return (1.0 - (1.0 - SCHEDULED_AMORT_RATE(wam, gwac)) * (1.0 - smm / 100.0));
} // TOTAL_AMORT_RATE