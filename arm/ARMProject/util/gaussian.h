/*
 * $Log: gaussian.h,v $
 * Revision 1.7  2004/02/02 16:20:15  ebenhamou
 * formatting
 *
 * Revision 1.6  2004/02/02 16:19:14  ebenhamou
 * remove duplicate line!
 *
 * Revision 1.5  2003/09/16 16:19:29  mab
 * ajout de CDFInvNormal
 *
 * Revision 1.4  2003/06/24 14:21:14  jmprie
 * ajout fct de calcul des moments ou des cumul des puissance de X qd
 *  X est gaussienne
 *
 * Revision 1.3  2002/11/25 16:31:54  mab
 * Formatting
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : gaussian.h                                                   */
/*                                                                            */
/* DESCRIPTION : gaussian utilities                                           */
/*                                                                            */
/* DATE        : Wed Apr 17 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _GAUSSIAN_H
#define _GAUSSIAN_H





extern  double  dNormal(double);
extern  double  cdfNormal(double);
extern  double  qCDFNormal(const double& z__);
extern  double  d2Normal(double, double , double);
extern  double  cdf2Normal(double, double, double);
extern  double  f2Normal(double, double, double, double, double);
extern  double  m2Normal(double, double, double);
extern  double  INV_PART_FUNC_NOR(double);
extern  double  maxLogNormalProba(double, double, double, double); 
extern  double cumulNormalMoment(double h,int rank,double mean,double sigma,
                                 double *moment);
extern  double normalMoment(int rank,double mean, double sigma,
                            double *moment);


extern double CDFInvNormal(const double& u);


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
