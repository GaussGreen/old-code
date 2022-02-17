/*
 * $Log: mcrnpricer.h,v $
 * Revision 1.4  2003/06/18 15:11:29  ebenhamou
 * remove include
 *
 * Revision 1.3  2003/06/18 11:31:00  ebenhamou
 * remove STL warning
 *
 */


/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : mcrnpricer.h                                                 */
/*                                                                            */
/* DESCRIPTION : classe pricers pour montecarlo forward neutre                */
/*                                                                            */
/* DATE        :     Sep 15 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _MCRNPRICER_H
#define _MCRNPRICER_H


#include "mcpricer.h"

#define SUMMATION_INIT_VALUE 1000
#define SUMMATION_GROWTH_FACTOR 10


class ARM_RNMonteCarloPricer : public ARM_MonteCarloPricer
{
    public :
        ARM_RNMonteCarloPricer(void)
        {
        };

        ARM_RNMonteCarloPricer(ARM_Security *sec, ARM_Model *mod) : 
                                                 ARM_MonteCarloPricer(sec, mod) 
        {
        };

        double Price(void);
        

};

#endif
