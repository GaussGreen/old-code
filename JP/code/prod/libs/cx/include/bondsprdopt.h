/*
*****************************************************************************
** Product header for bond spread option pricer.
*****************************************************************************
*/

#ifndef CX_BONDSPRDOPT_F_H
#define CX_BONDSPRDOPT_F_H

#include <crxflow/include/crxdata.h>

/*f
***************************************************************************
** Bond spread option calculator
***************************************************************************
*/
CrxTBondSpreadOptionCalc* CrxBondSpreadOptionCalc(
CrxTBondSpreadOption*   deal,             /* (I) */
long                    distType,         /* (I) */
TDate                   today,            /* (I) */
CrxTBondSpreadVolCurve* volCurve,         /* (I) */
CrxTBondPrice*          bondPrice,        /* (I) */
CrxTBondPrice*          refBondPrice,     /* (I) */
CrxTBondRepoCurve*      repoCurve,        /* (I) */
CrxTBondRepoCurve*      refRepoCurve,     /* (I) */
TCurve*                 discCurve         /* (I) */
);

void CrxBondSpreadOptionCalcLogInputs (
    double                   notional,
    CrxTBondSpreadOption*    deal,             
    long                     distType,         
    TDate                    today,            
    CrxTBondSpreadVolCurve*  volCurve,
    CrxTBondPrice*           bondPrice,
    CrxTBondPrice*           refBondPrice,
    CrxTBondRepoCurve*       repoCurve,        
    CrxTBondRepoCurve*       refRepoCurve,     
    TCurve*                  discCurve);

#endif

