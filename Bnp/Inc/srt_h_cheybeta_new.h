/*--------------------------------------------------------------
        FILE: srt_h_cheybeta_new.h
        PURPOSE: Cheyette beta model interface (new)
        AUTHOR: Dimitri Mayevski
        DATE: 20/09/2002
  --------------------------------------------------------------*/

#ifndef __SRT_H_CHEYBETA_NEW_H__
#define __SRT_H_CHEYBETA_NEW_H__

typedef struct _SCheyBeta
{
    double  lambda;
    int     nsig;
    double *sig, *sigtms, *beta;
    char*   ycname;
    long    today;
} SCheyBeta;

double CheyBeta_Vol(
    double x,
    double phi,
    double sig,
    double beta,
    double gamT,
    double fr,
    double xcut_u,
    double xcut_d);

Err CheyBeta_Init(SCheyBeta* pmdl, char* und_name);
Err CheyBeta_Free(SCheyBeta* pmdl);

#endif /* ifndef __SRT_H_CHEYBETA_NEW_H__ */