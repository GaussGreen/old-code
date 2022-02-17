/*
 * $Log: mcpricer.h,v $
 * Revision 1.6  2003/06/18 15:10:40  ebenhamou
 * remove include
 *
 * Revision 1.5  2003/06/18 11:30:00  ebenhamou
 * remove STL warning
 *
 * Revision 1.4  2001/07/30 09:12:26  smysona
 * Ajout de methodes pour les grecques
 *
 * Revision 1.3  2000/11/13 17:52:27  sgasquet
 * Retrait fonction virtuelle pure
 *
 * Revision 1.2  1999/11/05 15:41:19  nicolasm
 * Ajout constructeur par defaut
 *
 * Revision 1.1  1998/11/20 11:07:44  nicolasm
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : mcpricer.h                                                   */
/*                                                                            */
/* DESCRIPTION : classe pricers pour montecarlo                               */
/*                                                                            */
/* DATE        : Fri Sep 11 1998                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _MCPRICER_H
#define _MCPRICER_H

#include "pricer.h"
#include "armglob.h"

class ARM_MonteCarloPricer : public ARM_Pricer
{
    public :
        ARM_MonteCarloPricer(void)
        {
        };

        ARM_MonteCarloPricer(ARM_Security *sec, ARM_Model *mod) : 
                                                     ARM_Pricer(sec, mod) 
        {
            itsIsValidSdev = K_FALSE;

            itsSdev = 0.0;
        };

        double Price(void)
  {
            throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                  "Unimplemented <Price> method");
  }


        double GetEstimatedSdev(void);

        void SetEstimatedSdev(double new_sdev);

        inline double GetDelta(void)
        {
            return itsDelta;
        }

        inline void SetDelta(double x)
        {
            itsDelta = x;
        }

        inline double GetVega(void)
        {
            return itsVega;
        }

        inline void SetVega(double x)
        {
            itsVega = x;
        }

        inline double GetGamma(void)
        {
            return itsGamma;
        }

        inline void SetGamma(double x)
        {
            itsGamma = x;
        }


 




        

    private :

        ARM_Booleen itsIsValidSdev;
        
        double itsSdev;

        double itsDelta; 

        double itsVega;

        double itsGamma;

        
};

#endif
