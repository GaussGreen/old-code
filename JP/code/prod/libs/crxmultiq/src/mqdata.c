/*
***************************************************************************
** SOURCE FILE: mqdata.c
**
** Data manipulation functions for the MQDATA structure.
***************************************************************************
*/

#include "mqdata.h"

#include <ctype.h>
#include <crxflow/include/crxmacros.h>

/**
***************************************************************************
** Constructor for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataMake(
double          optExpy,             /* (I) */
double          fwdRate,             /* (I) */
double          sigATM,              /* (I) */
double          optATM,              /* (I) */
long            optType,             /* (I) */
long            nbQL,                /* (I) */
double*         kL,                  /* (I) [nbQL] */
double*         dL,                  /* (I) [nbQL] */
double*         xL,                  /* (I) [nbQL+1] */
double*         bL,                  /* (I) [nbQL] */
double*         qL,                  /* (I) [nbQL] */
long            nbQR,                /* (I) */
double*         kR,                  /* (I) [nbQR] */
double*         dR,                  /* (I) [nbQR] */
double*         xR,                  /* (I) [nbQR+1] */
double*         bR,                  /* (I) [nbQR] */
double*         qR,                  /* (I) [nbQR] */
double          sigMQ,               /* (I) */
double          optMQ,               /* (I) */
double          C,                   /* (I) */
double          K,                   /* (I) */
long            calibType,           /* (I) */
long            sSteps,              /* (I) */
double          sDelta,              /* (I) */
double          fwdTol,              /* (I) */
double          atmTol,              /* (I) */
long            calcFwd,             /* (I) */
double          muMQ                 /* (I) */
)
{
    static char routine[] = "CrxMqdataMake";
    int status = FAILURE;

    MQDATA* p = NULL;

    REQUIRE(nbQL > 0);
    REQUIRE(kL != NULL);
    REQUIRE(dL != NULL);
    REQUIRE(xL != NULL);
    REQUIRE(bL != NULL);
    REQUIRE(qL != NULL);
    REQUIRE(nbQR > 0);
    REQUIRE(kR != NULL);
    REQUIRE(dR != NULL);
    REQUIRE(xR != NULL);
    REQUIRE(bR != NULL);
    REQUIRE(qR != NULL);

    p = NEW(MQDATA);
    if (p==NULL) goto done;

    p->optExpy         = optExpy;
    p->fwdRate         = fwdRate;
    p->sigATM          = sigATM;
    p->optATM          = optATM;
    p->optType         = optType;
    p->nbQL            = nbQL;
    COPY_ARRAY (&(p->kL[0]), kL, double, p->nbQL);
    COPY_ARRAY (&(p->dL[0]), dL, double, p->nbQL);
    COPY_ARRAY (&(p->xL[0]), xL, double, p->nbQL+1);
    COPY_ARRAY (&(p->bL[0]), bL, double, p->nbQL);
    COPY_ARRAY (&(p->qL[0]), qL, double, p->nbQL);
    p->nbQR            = nbQR;
    COPY_ARRAY (&(p->kR[0]), kR, double, p->nbQR);
    COPY_ARRAY (&(p->dR[0]), dR, double, p->nbQR);
    COPY_ARRAY (&(p->xR[0]), xR, double, p->nbQR+1);
    COPY_ARRAY (&(p->bR[0]), bR, double, p->nbQR);
    COPY_ARRAY (&(p->qR[0]), qR, double, p->nbQR);
    p->sigMQ           = sigMQ;
    p->optMQ           = optMQ;
    p->C               = C;
    p->K               = K;
    p->calibType       = calibType;
    p->sSteps          = sSteps;
    p->sDelta          = sDelta;
    p->fwdTol          = fwdTol;
    p->atmTol          = atmTol;
    p->calcFwd         = calcFwd;
    p->muMQ            = muMQ;

    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxMqdataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Memory allocator for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataMakeEmpty(
long            nbQL,                /* (I) */
long            nbQR                 /* (I) */
)
{
    static char routine[] = "CrxMqdataMakeEmpty";
    int status = FAILURE;

    MQDATA* p = NULL;

    REQUIRE(nbQL > 0);
    REQUIRE(nbQR > 0);

    p = NEW(MQDATA);
    if (p==NULL) goto done;

    p->nbQL            = nbQL;

    p->nbQR            = nbQR;


    status = SUCCESS;

done:
    if (status != SUCCESS)
    {
        GtoErrMsg("%s: Failed!\n", routine);
        CrxMqdataFree(p);
        p=NULL;
    }

    return p;
}

/**
***************************************************************************
** Copy constructor for MQDATA
***************************************************************************
*/
MQDATA* CrxMqdataCopy(MQDATA* src)
{
    MQDATA* dst = NULL;
    if (src==NULL) return NULL;

    dst = CrxMqdataMake(src->optExpy,
                        src->fwdRate,
                        src->sigATM,
                        src->optATM,
                        src->optType,
                        src->nbQL,
                        src->kL,
                        src->dL,
                        src->xL,
                        src->bL,
                        src->qL,
                        src->nbQR,
                        src->kR,
                        src->dR,
                        src->xR,
                        src->bR,
                        src->qR,
                        src->sigMQ,
                        src->optMQ,
                        src->C,
                        src->K,
                        src->calibType,
                        src->sSteps,
                        src->sDelta,
                        src->fwdTol,
                        src->atmTol,
                        src->calcFwd,
                        src->muMQ);

    return dst;
}

/**
***************************************************************************
** Destructor for MQDATA
***************************************************************************
*/
void CrxMqdataFree(MQDATA *p)
{
    if (p != NULL)
    {
        FREE(p);
    }
}

