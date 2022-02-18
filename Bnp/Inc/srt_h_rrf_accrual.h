#ifndef SRT_H_RRF_ACCRUAL_H
#define SRT_H_RRF_ACCRUAL_H

#include "utallhdr.h"
#include "uterror.h"

Err srt_f_rrf_accrual(
    String  ResetType,
    long    Today,
    long    FirstFixingDate,
    long    LastFixingDate,
    double  Margin,
    double  BandWidth,
    double  DRS0T1,
    double  DRS0Tn,
    double  Correlation,
    String  VolType,
    double  CapBsVolT1,
    double  CapBsVolTn,
    double  FwdVolTn,
    double  DfPaymentDate,
    int     NumHermitePoints,
    double* Rrf_Accrual_Premium);

#endif
