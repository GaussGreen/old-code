/*
 * $Log: foncdig.h,v $
 * Revision 1.2  2002/11/25 14:38:04  mab
 * Formatting
 *
 */

/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
#ifndef FONCDIG_H_
#define FONCDIG_H_





#define EPSILON 0.000001                      




/*-----------------------------------------------------------------------*/
/*  calcul du prix d'une option digitale                                 */
/*-----------------------------------------------------------------------*/

double Spread_dig(double S1, double S2, double vol1, double vol2,
                  double d1, double d2, double rho, 
                  double r, double K, double t);                                                    


        
#endif     
/*-----------------------------------------------------------------------*/
/*---- End Of File ----*/
