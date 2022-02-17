/*
 * $Log: mcfnpricer.h,v $
 * Revision 1.3  2003/06/18 15:09:34  ebenhamou
 * remove include
 *
 * Revision 1.2  2003/06/18 11:29:11  ebenhamou
 * remove STL warning
 *
 * Revision 1.1  1998/11/20 11:08:04  nicolasm
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : mcfnpricer.h                                                 */
/*                                                                            */
/* DESCRIPTION : classe pricers pour montecarlo forward neutre                */
/*                                                                            */
/* DATE        :     Sep 15 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _MCFNPRICER_H
#define _MCFNPRICER_H

#include "mcpricer.h"

#define SUMMATION_INIT_VALUE 1000
#define SUMMATION_GROWTH_FACTOR 10


class ARM_FNMonteCarloPricer : public ARM_MonteCarloPricer
{
    public :
        ARM_FNMonteCarloPricer(ARM_Security *sec, ARM_Model *mod) : 
                                                     ARM_MonteCarloPricer(sec, mod) 
        {
        };

        double Price(void);
        

};

#endif
