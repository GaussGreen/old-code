/*
 * $Log: bondmath.h,v $
 * Revision 1.2  2003/05/06 09:48:17  mab
 * RCS Comments
 *
 */

/*----------------------------------------------------------------------------*
     
     bondmath.h
 
    Header for use of bond math analytics in "bondmath.cc"
 
*----------------------------------------------------------------------------*/
#ifndef _BONDMATH_H
#define _BONDMATH_H


/*---- local constants ----*/
 
#define kMaxYield       200.0
#define kMinYield       -99.0
#define kYieldErr       1.0e-7
#define kIterMax         30
 
 
#define kTooBig 1e+30
 



extern  double accruedInterest(double, int, double);

extern  double regularPresentValueToYield(double, double, 
                     double, int, double, double);

extern  double* approxRegularPvToYield(double, double, 
                    double, int, double, double, double *);

extern  double  regularYieldToPresentValue(double, double, 
                     double, int, double, double);

extern  double    regularDpDy(double, double, double, 
                             int, double, double);

extern  double    regularD2pDy2(double, double, double, int, double, double);
extern  double    regularForwardValue(double, double, 
                           double, int, double, double, int mmFlag=1);

extern    double    regularDfDr(double, double, double, 
                          int, double, double);

extern    double    regularReinvestRate(double, double, 
                        double, int, double, double, int mmFlag=1);

extern    double    dummyRegularYieldToPresentValue(double &, 
                         double &, double, void **);

extern    double    dummyRegularForwardValue(double &, double &, double, void **);


extern    double    oddFirstPresentValueToYield(double, double, 
                                double, int, double, double, double);

extern    double    dummyOddFirstYieldToPresentValue(double &, 
                                   double &, double, void **);

extern    double    oddFirstYieldToPresentValue(double, double, 
                                double, int, double, double, double);

extern    double    oddFirstDpDy(double, double, double, int, 
                                    double, double, double);

extern    double    oddFirstD2pDy2(double, double, double, int, 
                                 double, double, double);





#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
