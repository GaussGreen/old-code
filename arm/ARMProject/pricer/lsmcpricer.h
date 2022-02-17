/*
 * $Log: lsmcpricer.h,v $
 * Revision 1.8  2003/06/18 15:01:18  ebenhamou
 * remove included
 *
 * Revision 1.7  2003/06/18 11:05:38  ebenhamou
 * remove STL warning
 *
 * Revision 1.6  2001/07/30 09:12:10  smysona
 * ajout AssignHedge
 *
 * Revision 1.5  2001/05/23 17:43:16  smysona
 * inversion dans le cpp
 *
 * Revision 1.4  2001/04/03 12:02:31  nicolasm
 * cleanit des regressions
 *
 * Revision 1.3  2000/11/13 13:30:37  sgasquet
 * Ajout Constr par defaut et include exercise.h
 *
 * Revision 1.2  2000/10/24 12:26:16  sgasquet
 * First revision
 *
 * Revision 1.1  2000/10/24 12:20:40  sgasquet
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : LONGSTAFF SCHWARTZ mcpricer.h                                */
/*                                                                            */
/* DESCRIPTION : classe pricers pour montecarlo                               */
/*                                                                            */
/* DATE        : 2000														  */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef _LS_MCPRICER_H
#define _LS_MCPRICER_H

#include "mcpricer.h"
#include "armglob.h"
#include "linalg.h"
#include "exercise.h"


/***************************************************************************

   We propagate a multi dim pay off corresponding to a multi ex security
   (ex: flex caps)
   This pay off is propagated forward and backward

***************************************************************************/

class ARM_LSMonteCarloPricer : public ARM_MonteCarloPricer
{
private :

    int itsMode;

    public :
        ARM_LSMonteCarloPricer(ARM_Security *sec, ARM_Model *mod) 
			:ARM_MonteCarloPricer(sec, mod)
        {
        Init();
	};

    ARM_LSMonteCarloPricer(void)
    {
        Init();
    }

        double Price(void);

    long ComputeRegression(	ARM_Vector* OutParams, ARM_Vector* FittedVector, 
		ARM_Matrix* Value, long nbRelevantInputs, long nbRelevantFactors);

			
    inline int GetMode(void) 
			{
        return itsMode;
			}

    inline void SetMode(int m)
		{
        itsMode = m;
		}

    void Init(void);
			
    void EndPricing(void);


    void AssignHedge(ARM_Security* itsSecurity, ARM_Portfolio* controleVariatesPort,
                     double* coeffPortfolio, bool* vanillaOK);    
};

#endif
