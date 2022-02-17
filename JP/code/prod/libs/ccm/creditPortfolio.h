#ifndef _CREDIT_PORTFOLIO_
#define _CREDIT_PORTFOLIO_
/* $Header$ */
/* Credit Portfolio structure */
typedef struct
{
    long    nbNames;
    double *notional;               /* [nbNames] */
    double *recovery;               /* [nbNames] */
    double *nameMaturity;           /* [nbNames] */
    double  strike1;
    double  strike2;
} CREDIT_PORTFOLIO;

CREDIT_PORTFOLIO CreditPortfolioCreate(
        double *notional,           /* (I) [nbNames] */
        double *recovery,           /* (I) [nbNames] */
        double *nameMaturity,       /* (I) [nbNames] */
        long nbNames,               
        double strike1,             
        double strike2);

void CreditPortfolioFree(CREDIT_PORTFOLIO *cp);
#endif
