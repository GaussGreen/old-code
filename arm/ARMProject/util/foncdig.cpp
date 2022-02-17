/*
 * $Log: foncdig.cpp,v $
 * Revision 1.4  2002/11/25 14:38:20  mab
 * Replaced & by &&
 *
 * Revision 1.3  2002/05/30 13:55:06  mab
 * comments : double Probamax(double A1, double A2, double B1, double B2, double B)
 *
 */


/*--------------------------------------------------------------------*/

#include "foncstd.h"
#include "foncdig.h"
#include "math.h"


#define EPSILON 0.000001                      
                   




/*--------------------------------------------------------------------*/
/*  calcul du prix d'une option digitale                              */  
/*--------------------------------------------------------------------*/

double Spread_dig(double S1, 
                  double S2, 
                  double vol1, 
                  double vol2,
                  double d1,
                  double d2, 
                  double rho, 
                  double r, 
                  double K, 
                  double t)                                                    

{ 
    double A1, A2, B1, B2;
    double borne = 0.99;
 


    // on price avec une correlation
    // de -1, les spreads ayant une correlation comprise entre -0.99 et -1
    
    if (( rho >= -1.0 ) && ( rho <= -borne ))
    {
       B1 = vol2 * sqrt(t);
       B2 = -vol1 * sqrt(t);
       A1 = S2 * exp((r-d2)*t - vol2 * vol2 * 0.5 * t);
       A2 = S1 * exp((r-d1)*t - vol1 * vol1 * 0.5 * t);

       return((Probamin(A1, A2, B1, B2, K)*exp(-r*t))); 
    }
    else
    {
       // on price avec une corr‰lation de 1 et de 0.99 et on interpole

       if (( rho <= 1.0 ) && ( rho >= borne ))
       {
          double limit_proba, up_proba, proba;

          B1 = vol2 * sqrt(t);
          B2 = vol1 * sqrt(t);
          A1 = S2 * exp((r-d2)*t - vol2 * vol2 * 0.5 * t);
          A2 = S1 * exp((r-d1)*t - vol1 * vol1 * 0.5 * t);

          limit_proba=Proba0(S1, S2, vol1, vol2, d1, d2, borne, r, K, t);
          up_proba=Probamax(A1, A2, B1, B2, K);
          proba = limit_proba * (up_proba - limit_proba)/(1.0-borne)*(rho-borne);
  

          return(Probamax(A1, A2, B1, B2, K)*exp(-r*t));
       }
       else
       {
          // correlation differente de 1 en valeur absolue

          return(Proba0(S1, S2, vol1, vol2, d1, d2, rho, r, K, t)*exp(-r*t));
       }
    }
}


/*--------------------------------------------------------------------------*/
/*---- End Of File ----*/
