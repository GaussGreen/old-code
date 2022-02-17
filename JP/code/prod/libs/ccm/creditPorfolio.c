#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "proba_utils.h"
#include "error2.h"
#include "creditPortfolio.h"

/* ---------------------------------------------------------------------------
// CreditPortfolioCreate
//
*/
CREDIT_PORTFOLIO CreditPortfolioCreate(     double *notional,
                                            double *recovery,
                                            double *nameMaturity,
                                            long nbNames,
                                            double strike1,
                                            double strike2)
{
    int status = FAILURE;
    static char routine[] = "CreditPortfolioCreate";
    CREDIT_PORTFOLIO cp;
    memset(&cp,0,sizeof(CREDIT_PORTFOLIO));

    // malloc
    cp.notional = malloc(nbNames*sizeof(double));
    if(cp.notional == NULL) goto RETURN;
    cp.recovery = malloc(nbNames*sizeof(double));
    if(cp.recovery == NULL) goto RETURN;
    cp.nameMaturity = malloc(nbNames*sizeof(double));
    if(cp.nameMaturity == NULL) goto RETURN;

    // memcpy
    memcpy(cp.notional, notional, nbNames*sizeof(double));
    memcpy(cp.nameMaturity, nameMaturity, nbNames*sizeof(double));
    memcpy(cp.recovery, recovery, nbNames*sizeof(double));
    cp.strike1 = strike1;
    cp.strike2 = strike2;
    cp.nbNames = nbNames;
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CreditPortfolioFree(&cp);
    }
    return cp;
}

/* ---------------------------------------------------------------------------
// CreditPortfolioFree
//
*/
void CreditPortfolioFree(CREDIT_PORTFOLIO *cp)
{
    if(cp->nameMaturity) free(cp->nameMaturity);
    if(cp->notional) free (cp->notional);
    if(cp->recovery) free(cp->recovery);
}

