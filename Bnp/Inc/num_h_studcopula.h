
#ifndef num_h_studcopula_h
#define num_h_studcopula_h

/**********************************************************************************************************************************/
//
// Stud_Deg_Fit
// calculates the best-fit studnet degree and linear correlation parameters
//
/**********************************************************************************************************************************/
void Stud_Deg_Fit(
    double*        data1,            // input historical data series
    double*        data2,            // input historical data series
    unsigned long  n,                // input length of data series
    unsigned long* m,                // output best-fit Student degree
    double*        linear_corr,      // output linear correlation obtained by best fit
    unsigned long  Max_Stud_Degree,  // The maximum allowed value of the student degree
    unsigned long
        RankToLin_acc,        // Accuracy parameter for calculating the linear correlation from rank
    unsigned long Max_Steps)  // The maximum number of iterations allowed in the minimization loop
;
// For a given value of the `linear correlation' parameter, calculate the value of m between 0
// (Gaussian) and Max_Stud_Degree that maximizes the log-likelihood function
unsigned long Student_degree_MLE(
    double        linear_corr,
    double*       changes_1,
    double*       changes_2,
    unsigned long len,
    unsigned long Max_Stud_Degree)
;
// Calculates the log-likelihood function for a bivariate Student copula
double log_Stud_ML(
    unsigned long degree, double linear_corr, double* x1, double* x2, unsigned long NumData);

/**************************************************************************************************************************************/
//
// log_Student_copula_dist
//
// Returns the log of Student's t bivariate copula
//
// Inputs: x1		value at which the probability density is to be evaluated (can be any real
// number)
//		   x2		value at which the probability density is to be evaluated (can be
//any real number) 		   m		degree of the distribution (0 is Gaussian, otherwise can be any
// positive integer)
//         rho		correlation parameter
//
/**************************************************************************************************************************************/
double log_Student_copula_dist(double x1, double x2, unsigned long m, double rho);

/**************************************************************************************************************************************/
//
// log_Student_dist
//
// Returns the log of Student's t probability density
//
// Inputs: x		value at which the probability density is to be evaluated (can be any real
// number)
//		   m		degree of the distribution (0 is Gaussian, otherwise ca be any
// positive integer)
//
/**************************************************************************************************************************************/
double log_Student_dist(double x, unsigned long m);

#endif

