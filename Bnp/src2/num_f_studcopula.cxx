
// C. Hunter: 22 January 2001
// Created Source File

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_studcopula.h"
#define NRANSI

/**********************************************************************************************************************************/
//
// Stud_Deg_Fit
// calculates the best-fit studnet degree and linear correlation parameters
//
/**********************************************************************************************************************************/
void Stud_Deg_Fit(
    double *data1,       // input historical data series
    double *data2,       // input historical data series
    unsigned long n,     // input length of data series
    unsigned long &m,    // output best-fit Student degree
    double &linear_corr, // output linear correlation obtained by best fit
    unsigned long
        Max_Stud_Degree, // The maximum allowed value of the student degree
    unsigned long RankToLin_acc, // Accuracy parameter for calculating the
                                 // linear correlation from rank
    unsigned long Max_Steps) // The maximum number of iterations allowed in the
                             // minimization loop
{

  // Calculate the vector of changes.
  double *changes_1, *changes_2;
  changes_1 = dvector(1, n);
  changes_2 = dvector(1, n);
  for (unsigned long i = 1; i < n; i++) {
    changes_1[i] = data1[i] - data1[i - 1];
    changes_2[i] = data2[i] - data2[i - 1];
  }

  // Compute the historical linear correlation (set the last two arguments to 2
  // to calculate the correlation of the data series)        , the historical
  // rank correlation and the the maximum likelihood for the first 10 degrees.
  linear_corr = srt_f_corr(changes_1, changes_2, n - 1, 2, 2);
  double rank_corr = spear_rho(changes_1, changes_2, n - 1);
  unsigned long MLE = Student_degree_MLE(linear_corr, changes_1, changes_2,
                                         n - 1, Max_Stud_Degree);

  // Calculate the linear correlation compatible with the MLE value of m and the
  // historical rank correlation. Use it to re-estimate the MLE.  Continue
  // iterating until the MLE value for m stops changing or we have looped too
  // many times
  unsigned long new_MLE = MLE;
  unsigned long steps = 0;
  do {
    MLE = new_MLE;
    linear_corr = RankToLinCorr(rank_corr, MLE, RankToLin_acc);
    new_MLE = Student_degree_MLE(linear_corr, changes_1, changes_2, n - 1);
  } while (new_MLE != MLE && ++steps < Max_Steps);

  // Free the memory and return the MLE
  free_dvector(changes_1, 1, n);
  free_dvector(changes_2, 1, n);

  m = new_MLE;
  //	return err;
}

// For a given value of the `linear correlation' parameter        , calculate
// the value of m between 0 (Gaussian) and Max_Stud_Degree that maximizes the
// log-likelihood function
unsigned long Student_degree_MLE(double linear_corr, double *changes_1,
                                 double *changes_2, unsigned long len,
                                 unsigned long Max_Stud_Degree) {
  double max_likelihood =
      log_Stud_ML(0, linear_corr, changes_1, changes_2, len);
  unsigned long MLE = 0;

  for (unsigned long j = 1; j <= Max_Stud_Degree; j++)
    if (log_Stud_ML(j, linear_corr, changes_1, changes_2, len) > max_likelihood)
      MLE = j;

  return MLE;
}

// Calculates the log-likelihood function for a bivariate Student copula
double log_Stud_ML(unsigned long degree, double linear_corr, double *x1,
                   double *x2, unsigned long NumData) {
  double sum = 0.0;
  for (unsigned long i = 0; i < NumData; i++)
    sum += log_Student_copula_dist(x1[i], x2[i], degree, linear_corr) +
           log_Student_dist(x1[i], degree) + log_Student_dist(x2[i], degree);
  return sum;
}

/**************************************************************************************************************************************/
//
// log_Student_copula_dist
//
// Returns the log of Student's t bivariate copula
//
// Inputs: x1		value at which the probability density is to be
// evaluated (can be any real number)
//		   x2		value at which the probability density is to be
// evaluated (can be any real number) 		   m		degree of the
// distribution (0 is Gaussian        , otherwise can be any positive integer)
//         rho		correlation parameter
//
/**************************************************************************************************************************************/
double log_Student_copula_dist(double x1, double x2, unsigned long m,
                               double rho) {
  if (m == 0)
    return rho * x1 * x2 / (1 - rho * rho) - 0.5 * log(1 - rho * rho);
  else
    return gammln(0.5 * (m + 2.0)) + gammln(0.5 * m) -
           2.0 * gammln(0.5 * (m + 1.0)) - 0.5 * log(1 - rho * rho) -
           0.5 * (m + 1.0) *
               log(1 + (x1 * x1 - 2.0 * rho * x1 * x2 + x2 * x2) /
                           (1 - rho * rho) / m) +
           0.5 * (m + 1.0) * log((1 + x1 * x1 / m) * (1 + x1 * x1 / m));
}

/**************************************************************************************************************************************/
//
// log_Student_dist
//
// Returns the log of Student's t probability density
//
// Inputs: x		value at which the probability density is to be
// evaluated (can be any real number)
//		   m		degree of the distribution (0 is Gaussian ,
// otherwise ca be any positive integer)
//
/**************************************************************************************************************************************/
double log_Student_dist(double x, unsigned long m) {
  if (m == 0)
    return -0.5 * (x * x + log(2.0 * SRT_PI));
  else
    return gammln(0.5 * (m + 1.0)) - gammln(0.5 * m) - 0.5 * log(SRT_PI * m) -
           0.5 * (m + 1.0) * log(1 + x * x / m);
}

#undef NRANSI
