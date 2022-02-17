/*
*****************************************************************************
** Product header for bond price option pricer.
*****************************************************************************
*/

#ifndef CX_BONDPRCOPT_H
#define CX_BONDPRCOPT_H

#include <crxflow/include/crxdata.h>

/*f
***************************************************************************
** Bond price option calculator
***************************************************************************
*/
CrxTBondPriceOptionCalc* CrxBondPriceOptionCalc(
CrxTBondPriceOption*   deal,             /* (I) */
long                   distType,         /* (I) */
TDate                  today,            /* (I) */
CrxTBondPriceVolCurve* volCurve,         /* (I) */
CrxTBondPrice*         bondPrice,        /* (I) */
CrxTBondRepoCurve*     repoCurve,        /* (I) */
TCurve*                discCurve         /* (I) */
);

void CrxBondPriceOptionCalcLogInputs (
    double                  notional,
    CrxTBondPriceOption*    deal,             
    long                    distType,         
    TDate                   today,            
    CrxTBondPriceVolCurve*  volCurve,
    CrxTBondPrice*          bondPrice,
    CrxTBondRepoCurve*      repoCurve,        
    TCurve*                 discCurve);

#endif

