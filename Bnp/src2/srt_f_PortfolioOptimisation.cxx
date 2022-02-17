/**********************************************************************************************
Portfolio Optimisation for a certain amount of Contracts given their Average
Daily Return and Daily St Dev and Correlation Matrix. We optimise the Average
Return of the Portfolio for a given maximum loss and probaility of losing more
than that maximum loss (15%        , 2.5% and 0.5%). Author : Ezra Nahum (May
2001)
*********************************************************************************************/

#include "num_h_allhdr.h"
#include "srt_h_all.h"
#include "utallhdr.H>

double ComputeFunction(double lambda, long size, double *mu, double *sigma,
                       double **Correl, long ndays, double DrawDown, int level,
                       double *PortfComp) {
  int i, j;
  double **MatrixSystem = NULL;
  double **MatrixSystemInverted = NULL;
  double **righthandvector = NULL;
  double **systemsolution = NULL;
  double function;
  double average;
  double variance;

  /*Memory allocation*/
  MatrixSystem = dmatrix(1, size, 1, size);
  righthandvector = dmatrix(1, size, 1, 1);

  for (i = 1; i <= size; i++) {
    righthandvector[i][1] = (lambda - ndays) * mu[i] * DrawDown;
  }

  for (i = 1; i <= size; i++) {
    for (j = 1; j <= size; j++) {
      if (i == j) {
        MatrixSystem[i][j] = (lambda - ndays) * mu[i] * mu[i] * ndays +
                             level * level * ndays * sigma[i] * sigma[i];
      }

      else {
        MatrixSystem[i][j] =
            (lambda - ndays) * mu[i] * mu[j] * ndays +
            level * level * ndays * sigma[i] * sigma[j] * Correl[i][j];
      }
    }
  }

  /*Memory allocated within the function        , make sure to free it!*/
  MatrixSystemInverted = inverse_matrix(MatrixSystem, 1, size);

  /*Memory allocated within the function        , make sure to free it!*/
  systemsolution = product_matrix(MatrixSystemInverted, 1, size, 1, size,
                                  righthandvector, 1, size, 1, 1);

  average = 0.0;

  for (i = 1; i <= size; i++) {
    average += systemsolution[i][1] * mu[i];
  }

  variance = 0.0;
  for (i = 1; i <= size; i++) {
    for (j = 1; j <= size; j++) {
      variance += systemsolution[i][1] * systemsolution[j][1] * sigma[i] *
                  sigma[j] * Correl[i][j];
    }
  }

  variance = max(variance, 0);

  function = ndays * average - level * sqrt(ndays * variance) - DrawDown;

  for (i = 1; i <= size; i++) {
    PortfComp[i] = systemsolution[i][1];
  }

  /*free memory*/
  free_dmatrix(MatrixSystem, 1, size, 1, size);
  free_dmatrix(righthandvector, 1, size, 1, 1);
  free_dmatrix(MatrixSystemInverted, 1, size, 1, size);
  free_dmatrix(systemsolution, 1, size, 1, 1);

  return function;
}

Err srt_f_PortfolioOpt(long size, double *mu, double *sigma, double **Correl,
                       long ndays, double DrawDown, int level,
                       double *PortfComp) {
  Err err = NULL;
  int i, j;
  double fx, dfx;
  double precision = 0.00000001;
  /*double sigmamax;
  double mumin;*/
  double lambda, lambdamin, lambdamax;
  double step;
  double *lambdas = NULL;
  double *funcs = NULL;
  double *sgn = NULL;
  int chosgn;
  int iter = 0;
  int maxiter = 100;
  double bump = 0.0000001;
  double totalaverage1, totalaverage2;
  double lambda1, lambda2;

  /*
  sigmamax = sigma[1];
  mumin = fabs(mu[1]);
  for (i=2;i<=size;i++)
  {
          if (sigma[i] >= sigmamax) sigmamax = sigma[i];
          if (fabs(mu[i]) <= mumin) mumin = fabs(mu[i]);
  */
  /*Memory Allocation*/
  lambdas = dvector(1, 1000);
  funcs = dvector(1, 1000);
  sgn = dvector(1, 1000);

  lambdamin = -100;
  lambdamax = 100;
  step = (lambdamax - lambdamin) / 1000;

  for (j = 1; j <= 1000; j++) {
    lambdas[j] = lambdamin + (j - 1) * step;
    funcs[j] = ComputeFunction(lambdas[j], size, mu, sigma, Correl, ndays,
                               DrawDown, level, PortfComp);
    if (funcs[j] > 0.0) {
      sgn[j] = 1.0;
    } else {
      sgn[j] = -1.0;
    }
  }

  chosgn = 0;
  for (j = 2; j <= 1000; j++) {

    if (sgn[j] != sgn[j - 1]) {
      if (chosgn == 0) {
        lambdamin = lambdas[j];
        chosgn++;
      }

      else {
        lambdamax = lambdas[j];
      }
    }
  }
  /*Get first guess for lambda*/
  /*
  average = 0.0;

  for (i=1;i<=size;i++)
  {
          average+=mu[i];
  }

  variance = 0.0;
  for (i=1;i<=size;i++)
  {
          for(j=1;j<=size;j++)
          {
                  variance+=sigma[i]*sigma[j]*Correl[i][j];
          }
  }

  n0 = DrawDown/(ndays*average-level*sqrt(ndays*variance));
  */

  /*
  lambdamax = ndays;
  lambdamin = ndays-sigmamax/mumin;

   */
  /*
 lambda = (lambdamin+lambdamax)/2;
 */

  lambda = lambdamin;

  fx = ComputeFunction(lambda, size, mu, sigma, Correl, ndays, DrawDown, level,
                       PortfComp);

  dfx = (ComputeFunction(lambda + bump, size, mu, sigma, Correl, ndays,
                         DrawDown, level, PortfComp) -
         fx) /
        bump;

  while ((fabs(fx) > precision) && (iter < maxiter)) {
    lambda = lambda - fx / dfx;

    fx = ComputeFunction(lambda, size, mu, sigma, Correl, ndays, DrawDown,
                         level, PortfComp);

    dfx = (ComputeFunction(lambda + bump, size, mu, sigma, Correl, ndays,
                           DrawDown, level, PortfComp) -
           fx) /
          bump;

    iter++;
  }

  lambda1 = lambda;
  totalaverage1 = 0.0;
  for (i = 1; i <= size; i++) {
    totalaverage1 += PortfComp[i] * mu[i];
  }

  iter = 0;

  lambda = lambdamax;

  fx = ComputeFunction(lambda, size, mu, sigma, Correl, ndays, DrawDown, level,
                       PortfComp);

  dfx = (ComputeFunction(lambda + bump, size, mu, sigma, Correl, ndays,
                         DrawDown, level, PortfComp) -
         fx) /
        bump;

  while ((fabs(fx) > precision) && (iter < maxiter)) {
    lambda = lambda - fx / dfx;

    fx = ComputeFunction(lambda, size, mu, sigma, Correl, ndays, DrawDown,
                         level, PortfComp);

    dfx = (ComputeFunction(lambda + bump, size, mu, sigma, Correl, ndays,
                           DrawDown, level, PortfComp) -
           fx) /
          bump;

    iter++;
  }

  lambda2 = lambda;
  totalaverage2 = 0.0;
  for (i = 1; i <= size; i++) {
    totalaverage2 += PortfComp[i] * mu[i];
  }

  if (totalaverage1 > totalaverage2) {

    fx = ComputeFunction(lambda1, size, mu, sigma, Correl, ndays, DrawDown,
                         level, PortfComp);
  }

  if ((totalaverage1 < 0) && (totalaverage2 < 0)) {
    err = "no good";
  }

  /*Memory Allocation*/
  free_dvector(lambdas, 1, 1000);
  free_dvector(funcs, 1, 1000);
  free_dvector(sgn, 1, 1000);

  return err;
}

/*
Err srt_f_PortfolioOpt(long size        ,double *mu        , double *sigma ,
double
**Correl        ,long ndays        ,double DrawDown        ,int level        ,
double *PortfComp)
{
Err err = NULL;
int i;
double funcleft        , funcright        , funcmid;
double lambdaleft        , lambdaright        , lambdamid;
double lambdamidleft        ,lambdamidright;
double funcmidleft        ,funcmidright;
double precision;
double sigmamax;
double mumin;
double lambdamin        ,lambdamax;
int iter = 0;
int maxiter = 200;

sigmamax = sigma[1];
mumin = fabs(mu[1]);
for (i=2;i<=size;i++)
{
        if (sigma[i] >= sigmamax) sigmamax = sigma[i];
        if (fabs(mu[i]) <= mumin) mumin = fabs(mu[i]);
}




lambdamax = ndays;
lambdamin = ndays-sigmamax/mumin;


lambdaleft = lambdamin;
lambdaright = lambdamax;
lambdamid = (lambdamin+lambdamax)/2;
precision = 1.0;
funcleft =  ComputeFunction(lambdaleft        ,size        ,mu        ,sigma
,Correl        ,ndays ,DrawDown        ,level        , PortfComp); funcright =
ComputeFunction(lambdaright        ,size ,mu        ,sigma        ,Correl ,ndays
,DrawDown        ,level        , PortfComp);



        while ((iter < maxiter)&&(precision>=0.0000001))
        {

                funcmid = ComputeFunction(lambdamid        ,size        ,mu
,sigma        ,Correl ,ndays        ,DrawDown        ,level        , PortfComp);

                lambdamidleft = (lambdaleft+lambdamid)/2;
                lambdamidright = (lambdaright+lambdamid)/2;

                funcmidleft = ComputeFunction(lambdamidleft        ,size ,mu
,sigma ,Correl        ,ndays        ,DrawDown        ,level        , PortfComp);
funcmidright = ComputeFunction(lambdamidright        ,size        ,mu ,sigma
,Correl        ,ndays        ,DrawDown ,level        , PortfComp);

                if ((funcmid<funcmidleft)&&(funcmid<funcmidright))
                {
                        lambdaleft = lambdamidleft;
                        lambdaright = lambdamidright;
                }
                else if ((funcmidleft<funcmid)&&(funcmidleft<funcmidright))
                {
                        lambdaright = lambdamid;
                }

                else
                {
                        lambdaleft = lambdamid;
                }


                lambdamid = (lambdaleft+lambdaright)/2;
                precision = funcmidright;
                iter++;


        }

return err;

}
*/
