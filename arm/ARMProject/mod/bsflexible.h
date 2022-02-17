/*
 * $Log: bsflexible.h,v $
 * Revision 1.2  2003/09/09 09:06:47  mab
 * RCS Comments
 *
 */

/*----------------------------------------------------------------------------*
 
    bsflexible.h
 
    Header for use of options analytic functions
 
    Copyright (c) 2003 
 
*----------------------------------------------------------------------------*/
#ifndef _BSFLEXIBLE_H
#define _BSFLEXIBLE_H



extern double bsOption2(double spot,
                       double strike,
                       double volatility,
                       double dividend,
                       double discountRate,
                       double maturity,
                       double CallPut);

extern double bsflexible(double forward,
                       double totalvolatility,
                       double bondprice,
                       double strike,
                       double CallPut);

#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
